"""TRAPI client for querying NCATS Translator APIs.

Handles:
- Gene normalization using TCT name_resolver
- Batch TRAPI queries with neighborhood discovery pattern
- Response caching and rate limiting
- Parallel API queries with graceful degradation
- Progress callbacks for UI integration


"""

import json
from pathlib import Path
from typing import List, Dict, Any, Optional, Callable
from datetime import datetime
from pydantic import BaseModel, Field
import logging

try:
    from TCT import TCT, name_resolver, translator_metakg, translator_kpinfo, translator_query
    TCT_AVAILABLE = True
except ImportError:
    TCT_AVAILABLE = False

logger = logging.getLogger(__name__)


class TRAPIResponse(BaseModel):
    """Structured TRAPI query response."""

    query_id: str = Field(description="Unique query identifier")
    input_genes: List[str] = Field(description="Input gene CURIEs")
    target_disease: Optional[str] = Field(default=None, description="Target disease CURIE")
    edges: List[Dict[str, Any]] = Field(default_factory=list, description="Retrieved edges")
    metadata: Dict[str, Any] = Field(default_factory=dict, description="Query metadata")
    timestamp: datetime = Field(default_factory=datetime.now)
    apis_queried: int = Field(default=0, description="Number of APIs queried")
    apis_succeeded: int = Field(default=0, description="Number of successful APIs")


class TRAPIClient:
    """Client for querying NCATS Translator APIs using TCT library.

    Example:
        >>> client = TRAPIClient(cache_dir=Path("data/cache"))
        >>> genes = ["CD6", "IFITM3", "DPP4"]
        >>> response = client.query_gene_neighborhood(genes, disease_curie="MONDO:0100096")
        >>> print(f"Found {len(response.edges)} edges")
    """

    def __init__(
        self,
        cache_dir: Path = Path("data/cache"),
        timeout: int = 30,
        max_workers: int = 4,
    ):
        """Initialize TRAPI client.

        Args:
            cache_dir: Directory for caching responses
            timeout: API timeout in seconds
            max_workers: Number of parallel workers for batch queries
        """
        if not TCT_AVAILABLE:
            raise ImportError("TCT library not installed. Run: pip install TCT")

        self.cache_dir = cache_dir
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.timeout = timeout
        self.max_workers = max_workers

        # Translator resources (loaded lazily)
        self.APInames = None
        self.metaKG = None
        self.Translator_KP_info = None
        self._resources_loaded = False

    def _load_translator_resources(self) -> None:
        """Load Translator MetaKG and API information.

        Extracted from notebook cell 8. Uses fallback strategy if Plover APIs timeout.
        """
        if self._resources_loaded:
            return

        logger.info("Loading Translator resources...")

        try:
            # Try full resource loading
            self.APInames, self.metaKG, self.Translator_KP_info = (
                translator_metakg.load_translator_resources()
            )
            logger.info(f"Loaded resources: {len(self.APInames)} APIs, {len(self.metaKG)} MetaKG edges")

        except Exception as e:
            logger.warning(f"Error loading resources: {e}")
            logger.info("Using fallback: loading without Plover APIs...")

            # Fallback: load components separately
            self.Translator_KP_info, self.APInames = translator_kpinfo.get_translator_kp_info()
            self.metaKG = translator_metakg.get_KP_metadata(self.APInames)
            logger.info(f"Loaded resources (fallback): {len(self.APInames)} APIs, {len(self.metaKG)} MetaKG edges")

        self._resources_loaded = True

    def normalize_genes(self, gene_symbols: List[str]) -> Dict[str, str]:
        """Normalize gene symbols to CURIEs using TCT name_resolver.

        Extracted from notebook cell 4.

        Args:
            gene_symbols: List of HUGO gene symbols (e.g., ['CD6', 'IFITM3'])

        Returns:
            Dictionary mapping gene symbols to CURIEs
        """
        logger.info(f"Normalizing {len(gene_symbols)} genes using TCT name_resolver...")

        # Use TCT's batch_lookup for efficient normalization
        gene_info_dict = name_resolver.batch_lookup(gene_symbols)

        results = {}
        for gene_symbol in gene_symbols:
            if gene_symbol in gene_info_dict:
                info = gene_info_dict[gene_symbol]
                if hasattr(info, "curie") and info.curie:
                    results[gene_symbol] = info.curie
                    logger.debug(f"✓ {gene_symbol} → {info.curie}")
                else:
                    logger.warning(f"✗ {gene_symbol} - no CURIE found")
            else:
                logger.warning(f"✗ {gene_symbol} - not found in name resolver")

        logger.info(f"Successfully normalized {len(results)}/{len(gene_symbols)} genes")
        return results

    def query_gene_neighborhood(
        self,
        gene_symbols: List[str],
        disease_curie: Optional[str] = None,
        intermediate_categories: Optional[List[str]] = None,
        predicates: Optional[List[str]] = None,
        progress_callback: Optional[Callable[[str], None]] = None,
    ) -> TRAPIResponse:
        """Query Translator APIs for gene-disease paths through intermediate entities.

        Supports two query patterns:
        1. If intermediate_categories=None: Gene → [Anything] (neighborhood discovery, 1-hop)
        2. If intermediate_categories specified AND disease_curie: Gene → [Intermediate] → Disease (2-hop)

        Args:
            gene_symbols: List of gene symbols to query
            disease_curie: Disease CURIE for 2-hop queries (required if intermediate_categories provided)
            intermediate_categories: Optional list of biolink categories for intermediate nodes
                                   (e.g., ["biolink:Protein", "biolink:ChemicalEntity"])
            predicates: Optional list of predicates to filter (default: auto-discover)
            progress_callback: Optional callback for progress updates

        Returns:
            TRAPIResponse with edges and metadata
        """
        # Load Translator resources if not already loaded
        self._load_translator_resources()

        # Step 1: Normalize genes
        if progress_callback:
            progress_callback(f"Normalizing {len(gene_symbols)} genes...")

        gene_curie_map = self.normalize_genes(gene_symbols)
        if not gene_curie_map:
            raise ValueError("No genes were successfully normalized")

        input_gene_curies = list(gene_curie_map.values())

        # Step 2: Setup query parameters
        if progress_callback:
            progress_callback("Setting up TRAPI query...")

        input_gene_categories = ["biolink:Gene"]

        # Determine query pattern: 1-hop vs 2-hop
        use_two_hop = intermediate_categories and disease_curie

        if use_two_hop:
            # 2-HOP QUERY: Gene → [Intermediate] → Disease
            logger.info(f"Using 2-hop query: Gene → {intermediate_categories} → Disease")

            disease_categories = ["biolink:Disease"]

            # Auto-discover predicates for each hop
            hop1_predicates = list(
                set(
                    TCT.select_concept(
                        sub_list=input_gene_categories,
                        obj_list=intermediate_categories,
                        metaKG=self.metaKG,
                    )
                )
            )

            hop2_predicates = list(
                set(
                    TCT.select_concept(
                        sub_list=intermediate_categories,
                        obj_list=disease_categories,
                        metaKG=self.metaKG,
                    )
                )
            )

            logger.info(f"Hop 1 predicates (gene→intermediate): {len(hop1_predicates)}")
            logger.info(f"Hop 2 predicates (intermediate→disease): {len(hop2_predicates)}")

            # Create custom 2-hop TRAPI query_graph
            # Structure: n00 (genes) → e00 → n01 (intermediate) → e01 → n02 (disease)
            query_json = {
                'message': {
                    'query_graph': {
                        'nodes': {
                            'n00': {
                                'ids': input_gene_curies,
                                'categories': input_gene_categories
                            },
                            'n01': {
                                # No IDs = creative query (find any matching nodes)
                                'categories': intermediate_categories
                            },
                            'n02': {
                                'ids': [disease_curie],
                                'categories': disease_categories
                            }
                        },
                        'edges': {
                            'e00': {
                                'subject': 'n00',
                                'object': 'n01',
                                'predicates': hop1_predicates
                            },
                            'e01': {
                                'subject': 'n01',
                                'object': 'n02',
                                'predicates': hop2_predicates
                            }
                        }
                    }
                }
            }

            # Select APIs based on hop1 (most restrictive)
            selected_APIs = TCT.select_API(
                sub_list=input_gene_categories,
                obj_list=intermediate_categories,
                metaKG=self.metaKG,
            )

        else:
            # 1-HOP QUERY: Gene → [Anything] (neighborhood discovery)
            logger.info("Using 1-hop neighborhood discovery")

            target_categories = ["biolink:Disease"] if not intermediate_categories else intermediate_categories

            # Auto-discover predicates if not provided
            if predicates is None:
                predicates = list(
                    set(
                        TCT.select_concept(
                            sub_list=input_gene_categories,
                            obj_list=target_categories,
                            metaKG=self.metaKG,
                        )
                    )
                )
                logger.info(f"Auto-discovered {len(predicates)} predicates")

            # Select APIs capable of answering queries
            selected_APIs = TCT.select_API(
                sub_list=input_gene_categories,
                obj_list=target_categories,
                metaKG=self.metaKG,
            )

            # Format 1-hop query using TCT helper
            query_json = TCT.format_query_json(
                input_gene_curies,
                [],  # Empty target = neighborhood discovery
                input_gene_categories,
                target_categories,
                predicates,
            )

        logger.info(f"Selected {len(selected_APIs)} APIs for querying")

        # Step 4: Build API predicates dictionary
        API_predicates = {}
        API_withMetaKG = list(set(self.metaKG["API"]))
        for api in API_withMetaKG:
            API_predicates[api] = list(set(self.metaKG[self.metaKG["API"] == api]["Predicate"]))

        # Step 5: Execute parallel queries
        if progress_callback:
            progress_callback(f"Querying {len(selected_APIs)} Translator APIs...")

        logger.info(f"Querying {len(selected_APIs)} APIs in parallel (max_workers={self.max_workers})...")

        query_results = translator_query.parallel_api_query(
            query_json=query_json,
            select_APIs=selected_APIs,
            APInames=self.APInames,
            API_predicates=API_predicates,
            max_workers=self.max_workers
        )

        # Step 6: Convert results to edges list with provenance
        edges = []
        edge_sources = {}  # Track which API provided each edge
        query_result_id = 0  # Assign unique ID to each query result
        successful_apis = set()  # Track unique APIs that returned edges

        for k, v in query_results.items():
            if isinstance(v, dict):
                # Assign sequential query_result_id to this result
                v['query_result_id'] = query_result_id
                query_result_id += 1

                # Extract edge key (subject-predicate-object)
                edge_key = f"{v.get('subject', '')}_{v.get('predicate', '')}_{v.get('object', '')}"

                # Track source API for this edge (k might be API name or edge ID)
                # Store sources for deduplication
                if edge_key not in edge_sources:
                    edge_sources[edge_key] = []
                    edges.append(v)

                # Add knowledge_source metadata if not already present
                if 'knowledge_source' not in v:
                    v['knowledge_source'] = []

                # Track successful APIs from sources field (only count aggregators we queried)
                for source in v.get('sources', []):
                    if isinstance(source, dict) and 'resource_id' in source:
                        if source.get('resource_role') == 'aggregator_knowledge_source':
                            successful_apis.add(source['resource_id'])

        # Create CURIE to symbol reverse mapping for labels
        curie_to_symbol = {curie: symbol for symbol, curie in gene_curie_map.items()}

        # Create response
        response = TRAPIResponse(
            query_id=datetime.now().strftime("%Y%m%d_%H%M%S"),
            input_genes=input_gene_curies,
            target_disease=disease_curie,
            edges=edges,
            metadata={
                "gene_symbols": gene_symbols,
                "normalized_genes": gene_curie_map,
                "curie_to_symbol": curie_to_symbol,  # For label lookup
                "predicates_used": len(hop1_predicates) + len(hop2_predicates) if use_two_hop else len(predicates) if predicates else 0,
                "intermediate_categories": intermediate_categories,  # Track filter
                "query_pattern": "2-hop" if use_two_hop else "1-hop",
            },
            apis_queried=len(selected_APIs),
            apis_succeeded=len(successful_apis),
        )

        logger.info(f"Query complete: {len(edges)} edges from {response.apis_succeeded}/{response.apis_queried} APIs")

        # Cache response
        cache_file = self._cache_response(response)
        logger.info(f"Cached response to {cache_file}")

        if progress_callback:
            progress_callback(f"Query complete: {len(edges)} edges found")

        return response

    def _cache_response(self, response: TRAPIResponse) -> Path:
        """Cache TRAPI response to disk.

        Args:
            response: TRAPIResponse to cache

        Returns:
            Path to cached file
        """
        cache_file = self.cache_dir / f"tct_results_{response.timestamp.strftime('%Y%m%d_%H%M%S')}.json"
        with open(cache_file, "w") as f:
            json.dump(response.model_dump(mode="json"), f, indent=2)
        return cache_file

    def _load_cached_response(self, cache_file: Path) -> Optional[TRAPIResponse]:
        """Load cached TRAPI response.

        Args:
            cache_file: Path to cached response file

        Returns:
            TRAPIResponse if valid, None otherwise
        """
        if not cache_file.exists():
            return None

        try:
            with open(cache_file, "r") as f:
                data = json.load(f)
            return TRAPIResponse(**data)
        except Exception:
            return None

    def list_cached_queries(self) -> List[Dict[str, Any]]:
        """List all cached query files with metadata.

        Returns:
            List of dicts with cache file info (path, size, timestamp, etc.)
        """
        if not self.cache_dir.exists():
            return []

        cache_files = sorted(
            self.cache_dir.glob("tct_results_*.json"),
            key=lambda x: x.stat().st_mtime,
            reverse=True
        )

        results = []
        for f in cache_files:
            size_mb = f.stat().st_size / (1024 * 1024)
            timestamp_str = f.stem.replace("tct_results_", "")

            try:
                timestamp = datetime.strptime(timestamp_str, "%Y%m%d_%H%M%S")
                date_str = timestamp.strftime("%b %d, %Y %I:%M %p")
            except ValueError:
                date_str = timestamp_str
                timestamp = None

            results.append({
                "path": f,
                "size_mb": size_mb,
                "timestamp": timestamp,
                "timestamp_str": timestamp_str,
                "label": f"{date_str} ({size_mb:.1f} MB)"
            })

        return results
