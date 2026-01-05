"""TRAPI client for querying NCATS Translator APIs.

Handles:
- Gene normalization using TCT name_resolver
- Batch TRAPI queries with neighborhood discovery pattern
- Response caching and rate limiting
- Parallel API queries with graceful degradation
- Progress callbacks for UI integration


"""

import json
import requests
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from copy import deepcopy
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import List, Dict, Any, Optional, Callable, Tuple
from datetime import datetime
from pydantic import BaseModel, Field
import logging

# Node Normalizer API for CURIE-to-name resolution
NODE_NORMALIZER_URL = "https://nodenormalization-sri.renci.org/1.4/get_normalized_nodes"

try:
    from TCT import TCT, name_resolver, translator_metakg, translator_kpinfo
    TCT_AVAILABLE = True
except ImportError:
    TCT_AVAILABLE = False

from biograph_explorer.utils.biolink_predicates import (
    filter_predicates_by_granularity,
    filter_disease_bp_predicates,
    DISEASE_BP_INFORMATIVE_PREDICATES,
    DEFAULT_BP_INTERMEDIATE_CATEGORIES,
)

logger = logging.getLogger(__name__)

# CURIE prefix to biolink category mapping (matches graph_builder.py classification)
CURIE_PREFIX_TO_CATEGORY = {
    "NCBIGene:": "biolink:Gene",
    "HGNC:": "biolink:Gene",
    "UniProtKB:": "biolink:Protein",
    "PR:": "biolink:Protein",
    "CHEBI:": "biolink:ChemicalEntity",
    "CHEMBL:": "biolink:ChemicalEntity",
    "MONDO:": "biolink:Disease",
    "DOID:": "biolink:Disease",
    "GO:": "biolink:BiologicalProcess",
    "REACT:": "biolink:Pathway",
    "HP:": "biolink:PhenotypicFeature",
    "UBERON:": "biolink:AnatomicalEntity",
    "CL:": "biolink:Cell",
}


def _classify_curie_category(curie: str) -> str:
    """Classify a CURIE to its biolink category based on prefix."""
    for prefix, category in CURIE_PREFIX_TO_CATEGORY.items():
        if curie.startswith(prefix):
            return category
    return "biolink:NamedThing"  # Generic fallback


def _analyze_category_mismatches(
    edges: List[Dict[str, Any]],
    input_gene_curies: List[str],
    disease_curie: Optional[str],
    requested_categories: List[str],
) -> Dict[str, Any]:
    """Analyze intermediate nodes for category mismatches.

    Uses CURIE prefix to classify actual categories and compare against
    requested categories from the query.

    Args:
        edges: List of TRAPI edges
        input_gene_curies: List of input gene CURIEs (to exclude)
        disease_curie: Target disease CURIE (to exclude)
        requested_categories: List of requested intermediate categories

    Returns:
        Dict with category stats including mismatched counts
    """
    if not requested_categories:
        return {}

    input_set = set(input_gene_curies)

    # Identify intermediate nodes (not query genes, not disease)
    intermediate_nodes = set()
    for edge in edges:
        subj = edge.get('subject', '')
        obj = edge.get('object', '')
        if subj and subj not in input_set and subj != disease_curie:
            intermediate_nodes.add(subj)
        if obj and obj not in input_set and obj != disease_curie:
            intermediate_nodes.add(obj)

    # Classify and count
    actual_category_counts: Dict[str, int] = {}
    mismatched_count = 0

    for node in intermediate_nodes:
        actual_category = _classify_curie_category(node)
        actual_category_counts[actual_category] = actual_category_counts.get(actual_category, 0) + 1

        if actual_category not in requested_categories:
            mismatched_count += 1

    return {
        "requested_categories": requested_categories,
        "actual_category_counts": actual_category_counts,
        "mismatched_count": mismatched_count,
        "total_intermediates": len(intermediate_nodes),
    }


@dataclass
class APITiming:
    """Timing information for a single API query."""
    api_name: str
    duration_seconds: float
    success: bool
    edge_count: int
    error: Optional[str] = None


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

    # Node annotations from Node Annotator API (added after initial query)
    node_annotations: Optional[Dict[str, Dict[str, Any]]] = Field(
        default=None,
        description="Node annotation features by CURIE, for filtering and display"
    )
    annotation_metadata: Optional[Dict[str, Any]] = Field(
        default=None,
        description="Annotation metadata (filterable_attributes, searchable_attributes, etc.)"
    )


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
        timeout: int = 600,
    ):
        """Initialize TRAPI client.

        Args:
            cache_dir: Directory for caching responses
            timeout: API timeout in seconds (default: 240s)
        """
        if not TCT_AVAILABLE:
            raise ImportError("TCT library not installed. Run: pip install TCT")

        self.cache_dir = cache_dir
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.timeout = timeout

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

    def _optimize_query_json(
        self, query_json: Dict[str, Any], api_name: str, api_predicates: Dict[str, List[str]]
    ) -> Dict[str, Any]:
        """Optimize query JSON by filtering to predicates supported by the API.

        Matches TCT's optimize_query_json behavior - only optimizes edge 'e00'.

        Args:
            query_json: TRAPI query JSON
            api_name: Name of the API to query
            api_predicates: Dict mapping API names to their supported predicates

        Returns:
            Optimized query JSON with filtered predicates
        """
        # Use shallow copy like TCT does
        query_json_cur = query_json.copy()

        # Only optimize e00 predicates (matching TCT behavior)
        if api_name in api_predicates:
            try:
                e00_predicates = query_json_cur['message']['query_graph']['edges']['e00']['predicates']
                shared_predicates = list(
                    set(api_predicates[api_name]).intersection(e00_predicates)
                )
                if shared_predicates:
                    query_json_cur['message']['query_graph']['edges']['e00']['predicates'] = shared_predicates
            except KeyError:
                # No e00 edge or no predicates - keep original
                pass

        return query_json_cur

    def _query_api_with_timing(
        self,
        api_name: str,
        query_json: Dict[str, Any],
        api_predicates: Dict[str, List[str]],
    ) -> Tuple[Optional[Dict[str, Any]], APITiming]:
        """Query a single API and record timing information.

        Matches TCT's query_KP behavior but adds timing instrumentation.

        Args:
            api_name: Name of the API to query
            query_json: TRAPI query JSON
            api_predicates: Dict mapping API names to their supported predicates

        Returns:
            Tuple of (result dict or None, APITiming with timing info)
        """
        api_url = self.APInames[api_name]
        # Deep copy per-thread to avoid race conditions (matching TCT's query_KP behavior)
        query_copy = deepcopy(query_json)
        query_json_cur = self._optimize_query_json(query_copy, api_name, api_predicates)

        start_time = time.time()
        error_msg = None
        result = None
        edge_count = 0

        try:
            response = requests.post(api_url, json=query_json_cur, timeout=self.timeout)
            duration = time.time() - start_time

            if response.status_code == 200:
                data = response.json().get("message", {})
                kg = data.get("knowledge_graph", {})
                edges = kg.get("edges", {})

                if edges:
                    edge_count = len(edges)
                    result = data
                    logger.info(f"API '{api_name}' completed in {duration:.2f}s ({edge_count} edges)")
                elif "knowledge_graph" in data:
                    # API returned but with no edges (matches TCT behavior)
                    logger.info(f"API '{api_name}' completed in {duration:.2f}s (0 edges)")
                else:
                    logger.info(f"API '{api_name}' completed in {duration:.2f}s (no KG)")
            else:
                error_msg = f"HTTP {response.status_code}"
                logger.info(f"API '{api_name}' failed in {duration:.2f}s ({error_msg})")

        except requests.Timeout:
            duration = time.time() - start_time
            error_msg = "Timeout"
            logger.info(f"API '{api_name}' timed out after {duration:.2f}s")
        except Exception as e:
            duration = time.time() - start_time
            error_msg = str(e)
            logger.info(f"API '{api_name}' error in {duration:.2f}s: {error_msg}")

        timing = APITiming(
            api_name=api_name,
            duration_seconds=round(duration, 3),
            success=result is not None,
            edge_count=edge_count,
            error=error_msg,
        )

        return result, timing

    def _parallel_api_query_with_timing(
        self,
        query_json: Dict[str, Any],
        selected_apis: List[str],
        api_predicates: Dict[str, List[str]],
    ) -> Tuple[Dict[str, Any], List[APITiming]]:
        """Execute parallel API queries with timing instrumentation.

        Matches TCT's parallel_api_query behavior but adds timing.

        Args:
            query_json: TRAPI query JSON
            selected_apis: List of API names to query
            api_predicates: Dict mapping API names to their supported predicates

        Returns:
            Tuple of (merged results dict, list of APITiming for each API)
        """
        results = []
        timings = []

        with ThreadPoolExecutor(max_workers=len(selected_apis)) as executor:
            # Each thread does its own deepcopy in _query_api_with_timing (matching TCT behavior)
            future_to_api = {
                executor.submit(
                    self._query_api_with_timing, api_name, query_json, api_predicates
                ): api_name
                for api_name in selected_apis
            }

            for future in as_completed(future_to_api):
                api_name = future_to_api[future]
                try:
                    result, timing = future.result()
                    timings.append(timing)
                    if result and 'knowledge_graph' in result:
                        results.append(result)
                except Exception as e:
                    logger.warning(f"API '{api_name}' generated exception: {e}")
                    timings.append(APITiming(
                        api_name=api_name,
                        duration_seconds=0.0,
                        success=False,
                        edge_count=0,
                        error=str(e),
                    ))

        # Filter results to only include those with edges (matching TCT's included_KP_ID logic)
        valid_results = []
        for result in results:
            kg = result.get('knowledge_graph')
            if kg is not None and 'edges' in kg and len(kg['edges']) > 0:
                valid_results.append(result)

        # Merge all knowledge graph edges
        merged_edges = {}
        for result in valid_results:
            merged_edges.update(result['knowledge_graph']['edges'])

        # Sort timings by duration (slowest first) for logging
        timings.sort(key=lambda t: t.duration_seconds, reverse=True)

        return merged_edges, timings

    def normalize_genes(self, gene_symbols: List[str]) -> Dict[str, str]:
        """Normalize gene symbols to CURIEs using TCT name_resolver.

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

    def lookup_disease(self, query: str, max_results: int = 10) -> List[Dict[str, Any]]:
        """Look up disease by name using TCT name resolver.

        Uses autocomplete for partial matching (e.g., "Alzheimer" matches "Alzheimer disease").

        Args:
            query: Disease name or partial name to search
            max_results: Maximum number of results to return

        Returns:
            List of dicts with: curie, label, types, synonyms
        """
        if not query or len(query.strip()) < 2:
            return []

        logger.info(f"Looking up disease: '{query}'")

        try:
            # Use name_resolver with autocomplete for partial matching
            results = name_resolver.lookup(
                query.strip(),
                return_top_response=False,  # Get all matches
                return_synonyms=True,
                autocomplete=True,
            )

            # Handle single result (returns TranslatorNode, not list)
            if not isinstance(results, list):
                results = [results] if results else []

            # Convert to list of dicts, filtering for disease-like entities
            disease_results = []
            seen_curies = set()

            for node in results[:max_results * 2]:  # Get extra to allow filtering
                if not hasattr(node, 'curie') or not node.curie:
                    continue

                # Skip duplicates
                if node.curie in seen_curies:
                    continue
                seen_curies.add(node.curie)

                # Get types/categories
                types = getattr(node, 'types', []) or []

                # Filter for disease-related entities (prioritize MONDO, DOID, HP, OMIM)
                curie_prefix = node.curie.split(':')[0] if ':' in node.curie else ''
                is_disease_curie = curie_prefix in ['MONDO', 'DOID', 'HP', 'OMIM', 'ORPHANET', 'MESH']
                is_disease_type = any('Disease' in t or 'Phenotypic' in t for t in types)

                if is_disease_curie or is_disease_type:
                    disease_results.append({
                        'curie': node.curie,
                        'label': getattr(node, 'label', node.curie),
                        'types': types,
                        'synonyms': getattr(node, 'synonyms', []) or [],
                    })

                if len(disease_results) >= max_results:
                    break

            logger.info(f"Found {len(disease_results)} disease matches for '{query}'")
            return disease_results

        except Exception as e:
            logger.warning(f"Disease lookup failed for '{query}': {e}")
            return []

    def _resolve_curie_names(self, curies: List[str], timeout: int = 60) -> Dict[str, str]:
        """Resolve CURIEs to human-readable names using Node Normalizer API.

        Args:
            curies: List of CURIE strings (e.g., ['NCBIGene:7124', 'CHEBI:15377'])
            timeout: Request timeout in seconds

        Returns:
            Dict mapping CURIE to name. Falls back to CURIE itself if not found.
        """
        curie_to_name = {}

        if not curies:
            return curie_to_name

        try:
            response = requests.post(
                NODE_NORMALIZER_URL,
                json={"curies": curies},
                timeout=timeout
            )
            response.raise_for_status()
            data = response.json()

            for curie in curies:
                result = data.get(curie)
                if result and result.get('id', {}).get('label'):
                    curie_to_name[curie] = result['id']['label']
                else:
                    curie_to_name[curie] = curie  # Fallback to CURIE
            logger.info(f"Resolved names for {len(curie_to_name)} nodes ({sum(1 for c, n in curie_to_name.items() if c != n)} with labels)")

        except Exception as e:
            logger.warning(f"Name resolution failed: {e}, using CURIEs as fallback")
            curie_to_name = {curie: curie for curie in curies}

        return curie_to_name

    def query_gene_neighborhood(
        self,
        gene_symbols: List[str],
        disease_curie: Optional[str] = None,
        intermediate_categories: Optional[List[str]] = None,
        predicates: Optional[List[str]] = None,
        predicate_min_depth: int = 2,
        exclude_literature: bool = True,
        exclude_coexpression: bool = True,
        exclude_homology: bool = True,
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
            predicate_min_depth: Minimum depth in biolink hierarchy (0=all, 1=exclude root,
                               2=standard, 3=specific only). Default: 2
            exclude_literature: If True, exclude 'occurs_together_in_literature_with'. Default: True
            exclude_coexpression: If True, exclude 'coexpressed_with'. Default: True
            exclude_homology: If True, exclude 'homologous_to' and related predicates. Default: True
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

        # Step 5: Execute parallel queries with timing
        if progress_callback:
            progress_callback(f"Querying {len(selected_APIs)} Translator APIs...")

        query_start = time.time()
        logger.info(f"Starting parallel query: {len(selected_APIs)} APIs")

        query_results, api_timings = self._parallel_api_query_with_timing(
            query_json=query_json,
            selected_apis=selected_APIs,
            api_predicates=API_predicates,
        )

        query_duration = time.time() - query_start
        successful_count = sum(1 for t in api_timings if t.success)
        logger.info(f"Parallel query completed in {query_duration:.1f}s: {len(query_results)} edges from {successful_count}/{len(api_timings)} APIs")

        # Step 6: Convert results to edges list with provenance
        edges = []
        edge_sources = {}  # Track which API provided each edge
        query_result_id = 0  # Assign unique ID to each query result

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

        # Step 7: Post-filter edges by predicate granularity
        # APIs may return edges with predicates we didn't request, so filter them here
        edges_before_filter = len(edges)

        if progress_callback:
            progress_callback(f"Filtering {edges_before_filter} edges by predicate granularity...")

        filtered_edges = []
        for edge in edges:
            edge_predicate = edge.get('predicate', '')
            # Check if this predicate passes our filter
            filtered = filter_predicates_by_granularity(
                [edge_predicate], predicate_min_depth, exclude_literature,
                exclude_coexpression, exclude_homology
            )
            if filtered:
                filtered_edges.append(edge)

        edges = filtered_edges
        edges_removed = edges_before_filter - len(edges)
        logger.info(f"Post-filtered edges: {len(edges)} (removed {edges_removed} by predicate filter)")

        # Step 8: For 2-hop queries, remove orphan intermediate nodes
        # An intermediate node is orphan if it's not connected to both a query gene AND the disease
        if use_two_hop and disease_curie:
            # Build connectivity sets
            nodes_connected_to_genes = set()  # intermediates with edge from query gene
            nodes_connected_to_disease = set()  # intermediates with edge to disease

            for edge in edges:
                subj = edge.get('subject', '')
                obj = edge.get('object', '')

                # Edge from query gene to intermediate
                if subj in input_gene_curies:
                    nodes_connected_to_genes.add(obj)
                # Edge from intermediate to disease
                if obj == disease_curie:
                    nodes_connected_to_disease.add(subj)

            # Valid intermediates are connected to BOTH query genes AND disease
            valid_intermediates = nodes_connected_to_genes & nodes_connected_to_disease
            logger.info(f"Valid intermediates (connected to both gene and disease): {len(valid_intermediates)}")

            # Filter edges to only include those involving valid intermediates (or direct gene-disease)
            edges_before_orphan_filter = len(edges)
            valid_edges = []
            for edge in edges:
                subj = edge.get('subject', '')
                obj = edge.get('object', '')

                # Keep edge if:
                # 1. It's gene → valid_intermediate
                # 2. It's valid_intermediate → disease
                # 3. It's a direct gene → disease edge (rare but possible)
                if subj in input_gene_curies and obj in valid_intermediates:
                    valid_edges.append(edge)
                elif subj in valid_intermediates and obj == disease_curie:
                    valid_edges.append(edge)
                elif subj in input_gene_curies and obj == disease_curie:
                    valid_edges.append(edge)

            orphans_removed = edges_before_orphan_filter - len(valid_edges)
            edges = valid_edges
            edges_removed += orphans_removed
            logger.info(f"Removed {orphans_removed} edges to orphan intermediates, {len(edges)} edges remaining")

        # Step 8.5: Analyze category mismatches for 2-hop queries
        category_mismatch_stats = {}
        if use_two_hop and intermediate_categories:
            category_mismatch_stats = _analyze_category_mismatches(
                edges=edges,
                input_gene_curies=input_gene_curies,
                disease_curie=disease_curie,
                requested_categories=intermediate_categories,
            )
            if category_mismatch_stats.get("mismatched_count", 0) > 0:
                logger.warning(
                    f"Category mismatch: {category_mismatch_stats['mismatched_count']}/{category_mismatch_stats['total_intermediates']} "
                    f"intermediates don't match requested categories. Actual: {category_mismatch_stats['actual_category_counts']}"
                )

        # Step 9: Resolve common names for all nodes in filtered edges
        if progress_callback:
            progress_callback("Resolving common names for nodes...")

        unique_nodes = set()
        for edge in edges:
            unique_nodes.add(edge.get('subject', ''))
            unique_nodes.add(edge.get('object', ''))
        unique_nodes.discard('')

        curie_to_name = self._resolve_curie_names(list(unique_nodes))

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
                "curie_to_symbol": curie_to_symbol,  # For label lookup (query genes only)
                "curie_to_name": curie_to_name,  # Human-readable names for all nodes
                "predicates_used": len(hop1_predicates) + len(hop2_predicates) if use_two_hop else len(predicates) if predicates else 0,
                "intermediate_categories": intermediate_categories,  # Track filter
                "query_pattern": "2-hop" if use_two_hop else "1-hop",
                "predicate_min_depth": predicate_min_depth,
                "exclude_literature": exclude_literature,
                "exclude_coexpression": exclude_coexpression,
                "exclude_homology": exclude_homology,
                "edges_before_filter": edges_before_filter,
                "edges_removed": edges_removed,
                "api_timings": [asdict(t) for t in api_timings],
                "total_query_duration": round(query_duration, 3),
                "query_json": query_json,  # TRAPI query structure for download
                "category_mismatch_stats": category_mismatch_stats,  # Category mismatch tracking
            },
            apis_queried=len(selected_APIs),
            apis_succeeded=sum(1 for t in api_timings if t.success),
        )

        logger.info(f"Query complete: {len(edges)} edges from {response.apis_succeeded}/{response.apis_queried} APIs")

        # Cache response
        cache_file = self._cache_response(response)
        logger.info(f"Cached response to {cache_file}")

        if progress_callback:
            progress_callback(f"Query complete: {len(edges)} edges found")

        return response

    def query_disease_bioprocesses(
        self,
        disease_curie: str,
        filter_predicates: bool = True,
        timeout_override: Optional[int] = None,
        progress_callback: Optional[Callable[[str], None]] = None,
    ) -> Tuple[List[str], Dict[str, Any]]:
        """Query TRAPI for BiologicalProcesses associated with a disease.

        Stage 1 of the Gene → [Intermediate] → Disease-associated BiologicalProcess pattern.
        Discovers BiologicalProcesses linked to a disease and filters to informative predicates.

        Args:
            disease_curie: Disease CURIE (e.g., 'MONDO:0004975' for Alzheimer's)
            filter_predicates: If True, filter to informative predicates (default: True)
            timeout_override: Optional timeout override in seconds (default: use self.timeout)
            progress_callback: Optional callback for progress updates

        Returns:
            Tuple of:
                - List of unique BiologicalProcess CURIEs
                - Metadata dict with edge counts, predicate distribution, etc.
        """
        # Load Translator resources if not already loaded
        self._load_translator_resources()

        if progress_callback:
            progress_callback(f"Discovering BiologicalProcesses for {disease_curie}...")

        # Auto-discover predicates for Disease → BiologicalProcess
        disease_categories = ["biolink:Disease"]
        bp_categories = ["biolink:BiologicalProcess"]

        predicates = list(
            set(
                TCT.select_concept(
                    sub_list=disease_categories,
                    obj_list=bp_categories,
                    metaKG=self.metaKG,
                )
            )
        )
        logger.info(f"Disease→BP: {len(predicates)} predicates discovered")

        # Select APIs capable of Disease → BiologicalProcess queries
        selected_APIs = TCT.select_API(
            sub_list=disease_categories,
            obj_list=bp_categories,
            metaKG=self.metaKG,
        )
        logger.info(f"Disease→BP: {len(selected_APIs)} APIs selected")

        if progress_callback:
            progress_callback(f"Querying {len(selected_APIs)} APIs for disease BiologicalProcesses...")

        # Build 1-hop TRAPI query: Disease → BiologicalProcess
        query_json = {
            'message': {
                'query_graph': {
                    'nodes': {
                        'n00': {
                            'ids': [disease_curie],
                            'categories': disease_categories
                        },
                        'n01': {
                            'categories': bp_categories
                        }
                    },
                    'edges': {
                        'e00': {
                            'subject': 'n00',
                            'object': 'n01',
                            'predicates': predicates
                        }
                    }
                }
            }
        }

        # Build API predicates dictionary
        API_predicates = {}
        API_withMetaKG = list(set(self.metaKG["API"]))
        for api in API_withMetaKG:
            API_predicates[api] = list(set(self.metaKG[self.metaKG["API"] == api]["Predicate"]))

        # Use timeout override if provided
        original_timeout = self.timeout
        if timeout_override:
            self.timeout = timeout_override

        # Execute parallel queries
        query_start = time.time()
        try:
            query_results, api_timings = self._parallel_api_query_with_timing(
                query_json=query_json,
                selected_apis=selected_APIs,
                api_predicates=API_predicates,
            )
        finally:
            # Restore original timeout
            self.timeout = original_timeout

        query_duration = time.time() - query_start
        successful_count = sum(1 for t in api_timings if t.success)
        logger.info(f"Disease→BP query completed in {query_duration:.1f}s: {len(query_results)} edges from {successful_count}/{len(api_timings)} APIs")

        # Convert results to edges list
        edges = []
        for k, v in query_results.items():
            if isinstance(v, dict) and 'subject' in v and 'object' in v:
                edges.append(v)

        edges_before_filter = len(edges)

        # Filter to informative predicates if requested
        if filter_predicates:
            if progress_callback:
                progress_callback(f"Filtering {len(edges)} edges to informative predicates...")
            edges = filter_disease_bp_predicates(edges, use_informative_only=True)

        edges_after_filter = len(edges)
        logger.info(f"Disease→BP predicate filter: {edges_before_filter} → {edges_after_filter} edges")

        # Extract unique BiologicalProcess CURIEs
        bp_curies = list(set(edge.get('object', '') for edge in edges if edge.get('object', '')))

        # Extract Disease→BP edges with full details for tracking relationships
        disease_bp_edges = []
        bp_to_predicates: Dict[str, List[str]] = {}
        for edge in edges:
            bp_curie = edge.get('object', '')
            if bp_curie:
                disease_bp_edges.append({
                    'bp_curie': bp_curie,
                    'predicate': edge.get('predicate', ''),
                    'disease_curie': edge.get('subject', ''),
                    'sources': edge.get('sources', []),
                })
                # Group predicates by BP CURIE for easy lookup
                if bp_curie not in bp_to_predicates:
                    bp_to_predicates[bp_curie] = []
                pred = edge.get('predicate', '')
                if pred and pred not in bp_to_predicates[bp_curie]:
                    bp_to_predicates[bp_curie].append(pred)

        # Calculate predicate distribution
        predicate_counts = {}
        for edge in edges:
            pred = edge.get('predicate', '').replace('biolink:', '')
            predicate_counts[pred] = predicate_counts.get(pred, 0) + 1

        # Resolve BP names
        if progress_callback:
            progress_callback(f"Resolving names for {len(bp_curies)} BiologicalProcesses...")

        bp_names = self._resolve_curie_names(bp_curies)

        metadata = {
            'disease_curie': disease_curie,
            'total_edges_before_filter': edges_before_filter,
            'total_edges_after_filter': edges_after_filter,
            'unique_bioprocesses': len(bp_curies),
            'predicate_distribution': predicate_counts,
            'filter_applied': filter_predicates,
            'apis_queried': len(selected_APIs),
            'apis_succeeded': successful_count,
            'query_duration': round(query_duration, 3),
            'api_timings': [asdict(t) for t in api_timings],
            'bp_names': bp_names,  # CURIE → name mapping
            'bp_to_predicates': bp_to_predicates,  # CURIE → [predicates] for tracking relationships
            'disease_bp_edges': disease_bp_edges,  # Full edge list for optional rendering
        }

        if progress_callback:
            progress_callback(f"Found {len(bp_curies)} BiologicalProcesses (filtered from {edges_before_filter} edges)")

        return bp_curies, metadata

    def _cache_disease_bp_results(
        self,
        disease_curie: str,
        bp_curies: List[str],
        metadata: Dict[str, Any],
    ) -> Path:
        """Cache Disease → BiologicalProcess discovery results.

        Args:
            disease_curie: Disease CURIE
            bp_curies: List of discovered BP CURIEs
            metadata: Query metadata

        Returns:
            Path to cached file
        """
        # Sanitize disease CURIE for filename
        safe_curie = disease_curie.replace(":", "_").replace("/", "_")
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        cache_file = self.cache_dir / f"disease_bp_{safe_curie}_{timestamp}.json"

        cache_data = {
            'disease_curie': disease_curie,
            'bioprocess_curies': bp_curies,
            'metadata': metadata,
            'timestamp': timestamp,
        }

        with open(cache_file, "w") as f:
            json.dump(cache_data, f, indent=2)

        logger.info(f"Cached Disease→BP results to {cache_file}")
        return cache_file

    def list_cached_disease_bp_results(self, disease_curie: Optional[str] = None) -> List[Dict[str, Any]]:
        """List cached Disease → BiologicalProcess discovery results.

        Args:
            disease_curie: Optional filter by disease CURIE

        Returns:
            List of cache file info dicts
        """
        if not self.cache_dir.exists():
            return []

        pattern = "disease_bp_*.json"
        cache_files = sorted(
            self.cache_dir.glob(pattern),
            key=lambda x: x.stat().st_mtime,
            reverse=True
        )

        results = []
        for f in cache_files:
            try:
                with open(f, "r") as fp:
                    data = json.load(fp)

                cached_disease = data.get('disease_curie', '')

                # Filter by disease if specified
                if disease_curie and cached_disease != disease_curie:
                    continue

                bp_count = len(data.get('bioprocess_curies', []))
                timestamp_str = data.get('timestamp', '')

                try:
                    timestamp = datetime.strptime(timestamp_str, "%Y%m%d_%H%M%S")
                    date_str = timestamp.strftime("%b %d, %Y %I:%M %p")
                except ValueError:
                    date_str = timestamp_str
                    timestamp = None

                results.append({
                    'path': f,
                    'disease_curie': cached_disease,
                    'bp_count': bp_count,
                    'timestamp': timestamp,
                    'timestamp_str': timestamp_str,
                    'label': f"{cached_disease} - {bp_count} BPs ({date_str})",
                })

            except Exception as e:
                logger.warning(f"Error reading cache file {f}: {e}")
                continue

        return results

    def load_cached_disease_bp_results(self, cache_file: Path) -> Tuple[List[str], Dict[str, Any]]:
        """Load cached Disease → BiologicalProcess results.

        Args:
            cache_file: Path to cache file

        Returns:
            Tuple of (bp_curies list, metadata dict)

        Raises:
            FileNotFoundError: If cache file doesn't exist
            ValueError: If cache file is invalid
        """
        if not cache_file.exists():
            raise FileNotFoundError(f"Cache file not found: {cache_file}")

        with open(cache_file, "r") as f:
            data = json.load(f)

        bp_curies = data.get('bioprocess_curies', [])
        metadata = data.get('metadata', {})

        logger.info(f"Loaded {len(bp_curies)} cached BPs from {cache_file}")
        return bp_curies, metadata

    def query_gene_to_bioprocesses(
        self,
        gene_symbols: List[str],
        disease_curie: str,
        intermediate_categories: Optional[List[str]] = None,
        bp_curies: Optional[List[str]] = None,
        bp_metadata: Optional[Dict[str, Any]] = None,
        filter_disease_bp: bool = True,
        predicate_min_depth: int = 2,
        exclude_literature: bool = True,
        exclude_coexpression: bool = True,
        exclude_homology: bool = True,
        timeout_override: Optional[int] = None,
        progress_callback: Optional[Callable[[str], None]] = None,
    ) -> TRAPIResponse:
        """Execute two-stage Gene → [Intermediate] → Disease-associated BiologicalProcess query.

        Stage 1: Disease → BiologicalProcess (discover disease-associated BPs)
        Stage 2: Gene → [Intermediate] → BiologicalProcess (find pathways to those BPs)

        Args:
            gene_symbols: List of gene symbols to query
            disease_curie: Disease CURIE for Stage 1 BP discovery
            intermediate_categories: Intermediate entity types (default: DEFAULT_BP_INTERMEDIATE_CATEGORIES)
            bp_curies: Optional pre-computed BP CURIEs (skips Stage 1 if provided)
            bp_metadata: Optional metadata from cached Stage 1 results
            filter_disease_bp: Filter Stage 1 to informative predicates (default: True)
            predicate_min_depth: Granularity filter for Stage 2 (default: 2)
            exclude_literature: Exclude literature co-occurrence predicates (default: True)
            exclude_coexpression: Exclude 'coexpressed_with' predicate (default: True)
            exclude_homology: Exclude 'homologous_to' and related predicates (default: True)
            timeout_override: Optional timeout override in seconds
            progress_callback: Progress callback

        Returns:
            TRAPIResponse with merged Stage 2 edges and combined metadata
        """
        # Load Translator resources if not already loaded
        self._load_translator_resources()

        # Use default intermediate categories if not specified
        if intermediate_categories is None:
            intermediate_categories = DEFAULT_BP_INTERMEDIATE_CATEGORIES

        # Stage 1: Get disease-associated BiologicalProcesses (or use provided)
        if bp_curies is None:
            if progress_callback:
                progress_callback("Stage 1: Discovering disease BiologicalProcesses...")

            bp_curies, bp_metadata = self.query_disease_bioprocesses(
                disease_curie=disease_curie,
                filter_predicates=filter_disease_bp,
                timeout_override=timeout_override,
                progress_callback=progress_callback,
            )

            # Cache Stage 1 results
            self._cache_disease_bp_results(disease_curie, bp_curies, bp_metadata)

        if not bp_curies:
            raise ValueError(f"No BiologicalProcesses found for disease {disease_curie}")

        logger.info(f"Stage 1 complete: {len(bp_curies)} BiologicalProcesses")

        # Stage 2: Normalize genes
        if progress_callback:
            progress_callback(f"Stage 2: Normalizing {len(gene_symbols)} genes...")

        gene_curie_map = self.normalize_genes(gene_symbols)
        if not gene_curie_map:
            raise ValueError("No genes were successfully normalized")

        input_gene_curies = list(gene_curie_map.values())
        input_gene_categories = ["biolink:Gene"]
        bp_categories = ["biolink:BiologicalProcess"]

        # Stage 2: Discover predicates for each hop
        if progress_callback:
            progress_callback("Stage 2: Discovering predicates for 2-hop query...")

        # Hop 1: Gene → Intermediate
        hop1_predicates = list(
            set(
                TCT.select_concept(
                    sub_list=input_gene_categories,
                    obj_list=intermediate_categories,
                    metaKG=self.metaKG,
                )
            )
        )

        # Hop 2: Intermediate → BiologicalProcess
        hop2_predicates = list(
            set(
                TCT.select_concept(
                    sub_list=intermediate_categories,
                    obj_list=bp_categories,
                    metaKG=self.metaKG,
                )
            )
        )

        logger.info(f"Stage 2 predicates - Hop1 (Gene→Int): {len(hop1_predicates)}, Hop2 (Int→BP): {len(hop2_predicates)}")

        # Select APIs based on hop 1 (most restrictive)
        selected_APIs = TCT.select_API(
            sub_list=input_gene_categories,
            obj_list=intermediate_categories,
            metaKG=self.metaKG,
        )
        logger.info(f"Stage 2: {len(selected_APIs)} APIs selected")

        # Build 2-hop TRAPI query: Gene → [Intermediate] → BiologicalProcess
        # Use all BP CURIEs as endpoints
        if progress_callback:
            progress_callback(f"Stage 2: Querying {len(selected_APIs)} APIs for {len(input_gene_curies)} genes → {len(bp_curies)} BPs...")

        query_json = {
            'message': {
                'query_graph': {
                    'nodes': {
                        'n00': {
                            'ids': input_gene_curies,
                            'categories': input_gene_categories
                        },
                        'n01': {
                            # No IDs = creative query (find any matching intermediates)
                            'categories': intermediate_categories
                        },
                        'n02': {
                            'ids': bp_curies,
                            'categories': bp_categories
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

        # Build API predicates dictionary
        API_predicates = {}
        API_withMetaKG = list(set(self.metaKG["API"]))
        for api in API_withMetaKG:
            API_predicates[api] = list(set(self.metaKG[self.metaKG["API"] == api]["Predicate"]))

        # Use timeout override if provided
        original_timeout = self.timeout
        if timeout_override:
            self.timeout = timeout_override

        # Execute parallel queries
        query_start = time.time()
        try:
            query_results, api_timings = self._parallel_api_query_with_timing(
                query_json=query_json,
                selected_apis=selected_APIs,
                api_predicates=API_predicates,
            )
        finally:
            self.timeout = original_timeout

        query_duration = time.time() - query_start
        successful_count = sum(1 for t in api_timings if t.success)
        logger.info(f"Stage 2 query completed in {query_duration:.1f}s: {len(query_results)} edges from {successful_count}/{len(api_timings)} APIs")

        # Convert results to edges list
        edges = []
        for k, v in query_results.items():
            if isinstance(v, dict):
                edges.append(v)

        edges_before_filter = len(edges)

        # Post-filter edges by predicate granularity
        if progress_callback:
            progress_callback(f"Filtering {edges_before_filter} edges by predicate granularity...")

        filtered_edges = []
        for edge in edges:
            edge_predicate = edge.get('predicate', '')
            filtered = filter_predicates_by_granularity(
                [edge_predicate], predicate_min_depth, exclude_literature,
                exclude_coexpression, exclude_homology
            )
            if filtered:
                filtered_edges.append(edge)

        edges = filtered_edges
        edges_removed = edges_before_filter - len(edges)
        logger.info(f"Stage 2 post-filter: {len(edges)} edges (removed {edges_removed})")

        # Filter orphan intermediates (not connected to both gene AND BP)
        nodes_connected_to_genes = set()
        nodes_connected_to_bps = set()
        bp_curies_set = set(bp_curies)

        for edge in edges:
            subj = edge.get('subject', '')
            obj = edge.get('object', '')

            # Edge from query gene to intermediate
            if subj in input_gene_curies:
                nodes_connected_to_genes.add(obj)
            # Edge from intermediate to BP
            if obj in bp_curies_set:
                nodes_connected_to_bps.add(subj)

        # Valid intermediates are connected to BOTH query genes AND BPs
        valid_intermediates = nodes_connected_to_genes & nodes_connected_to_bps
        logger.info(f"Valid intermediates (connected to both gene and BP): {len(valid_intermediates)}")

        # Filter edges to only include those involving valid intermediates
        edges_before_orphan = len(edges)
        valid_edges = []
        for edge in edges:
            subj = edge.get('subject', '')
            obj = edge.get('object', '')

            # Keep edge if:
            # 1. Gene → valid_intermediate
            # 2. Valid_intermediate → BP
            # 3. Direct Gene → BP (rare but possible)
            if subj in input_gene_curies and obj in valid_intermediates:
                valid_edges.append(edge)
            elif subj in valid_intermediates and obj in bp_curies_set:
                valid_edges.append(edge)
            elif subj in input_gene_curies and obj in bp_curies_set:
                valid_edges.append(edge)

        orphans_removed = edges_before_orphan - len(valid_edges)
        edges = valid_edges
        edges_removed += orphans_removed
        logger.info(f"Removed {orphans_removed} orphan edges, {len(edges)} edges remaining")

        # Analyze category mismatches
        category_mismatch_stats = {}
        if intermediate_categories:
            category_mismatch_stats = _analyze_category_mismatches(
                edges=edges,
                input_gene_curies=input_gene_curies,
                disease_curie=None,  # BP query doesn't have disease as endpoint
                requested_categories=intermediate_categories,
            )
            if category_mismatch_stats.get("mismatched_count", 0) > 0:
                logger.warning(
                    f"Category mismatch: {category_mismatch_stats['mismatched_count']}/{category_mismatch_stats['total_intermediates']} "
                    f"intermediates don't match requested categories. Actual: {category_mismatch_stats['actual_category_counts']}"
                )

        # Resolve common names for all nodes in filtered edges
        if progress_callback:
            progress_callback("Resolving common names for nodes...")

        unique_nodes = set()
        for edge in edges:
            unique_nodes.add(edge.get('subject', ''))
            unique_nodes.add(edge.get('object', ''))
        unique_nodes.discard('')

        curie_to_name = self._resolve_curie_names(list(unique_nodes))

        # Merge with BP names from Stage 1 if available
        if bp_metadata and 'bp_names' in bp_metadata:
            curie_to_name.update(bp_metadata['bp_names'])

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
                "curie_to_symbol": curie_to_symbol,
                "curie_to_name": curie_to_name,
                "query_pattern": "gene-intermediate-disease-bp",
                "intermediate_categories": intermediate_categories,
                "bioprocess_endpoints": bp_curies,
                "bioprocess_count": len(bp_curies),
                "predicate_min_depth": predicate_min_depth,
                "exclude_literature": exclude_literature,
                "exclude_coexpression": exclude_coexpression,
                "exclude_homology": exclude_homology,
                "edges_before_filter": edges_before_filter,
                "edges_removed": edges_removed,
                "valid_intermediates": len(valid_intermediates),
                "api_timings": [asdict(t) for t in api_timings],
                "total_query_duration": round(query_duration, 3),
                "stage1_metadata": bp_metadata,
                "query_json": query_json,  # TRAPI query structure for download
                "category_mismatch_stats": category_mismatch_stats,  # Category mismatch tracking
            },
            apis_queried=len(selected_APIs),
            apis_succeeded=successful_count,
        )

        logger.info(f"Stage 2 complete: {len(edges)} edges from {response.apis_succeeded}/{response.apis_queried} APIs")

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
