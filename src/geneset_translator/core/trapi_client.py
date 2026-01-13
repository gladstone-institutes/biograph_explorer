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

from geneset_translator.utils.biolink_predicates import filter_predicates_by_granularity

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
    error_type: Optional[str] = None  # "client_timeout", "server_timeout", "server_error", "http_error", "exception"


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
    hpa_metadata: Optional[Dict[str, Any]] = Field(
        default=None,
        description="HPA annotation metadata (specificity categories, cell types, etc.)"
    )


class TRAPIClient:
    """Client for querying NCATS Translator APIs using TCT library.

    Example:
        >>> client = TRAPIClient(cache_dir=Path("data/cache"))
        >>> genes = ["CD6", "IFITM3", "DPP4"]
        >>> response = client.query_gene_neighborhood(genes, disease_curie="MONDO:0100096")
        >>> print(f"Found {len(response.edges)} edges")
    """
    # TODO: Identify features that we can port over to TCT
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

    def _filter_apis_by_exact_categories(
        self,
        selected_apis: List[str],
        subject_categories: List[str],
        object_categories: List[str],
        predicates: Optional[List[str]] = None,
        is_multihop: bool = False,
    ) -> Tuple[List[str], Dict[str, str]]:
        """Filter APIs to only those supporting exact requested categories and predicates.

        Filters out:
        - APIs that only support broader categories like biolink:NamedThing
        - APIs that have no matching predicates for the subject→object transition
        - Single-edge-only APIs for multi-hop queries

        Args:
            selected_apis: APIs from TCT.select_API()
            subject_categories: Subject categories (e.g., ["biolink:Gene"])
            object_categories: Object categories (e.g., ["biolink:Protein"])
            predicates: Optional list of predicates to check support for
            is_multihop: Whether this is a multi-hop query (filters single-edge APIs)

        Returns:
            Tuple of:
            - filtered_apis: List of APIs that support exact categories
            - filtered_out: Dict mapping API name to reason it was filtered
        """
        if self.metaKG is None:
            return selected_apis, {}

        filtered_apis = []
        filtered_out = {}

        # Known single-edge-only APIs (don't support multi-hop queries)
        single_edge_apis = {
            "Service Provider TRAPI",  # "smartAPI/team-specific endpoints only support single-edge queries"
        }

        # Broader categories to exclude when specific ones are requested
        broad_categories = {"biolink:NamedThing", "biolink:Entity", "biolink:ThingWithTaxon"}

        for api in selected_apis:
            # Check for single-edge-only APIs in multi-hop queries
            if is_multihop and api in single_edge_apis:
                filtered_out[api] = "Single-edge queries only (no multi-hop support)"
                continue

            # Get MetaKG entries for this API
            api_entries = self.metaKG[self.metaKG["API"] == api]

            if len(api_entries) == 0:
                filtered_out[api] = "No MetaKG entries found"
                continue

            # Get entries matching the subject→object transition
            transition_entries = api_entries[
                (api_entries["Subject"].isin(subject_categories)) &
                (api_entries["Object"].isin(object_categories))
            ]

            if len(transition_entries) == 0:
                # No direct subject→object entries - check what categories it supports
                api_objects = set(api_entries["Object"].unique())
                exact_matches = api_objects & set(object_categories)

                if not exact_matches:
                    broad_only = api_objects & broad_categories
                    if broad_only:
                        filtered_out[api] = f"Only broad categories: {broad_only}"
                    else:
                        supported = {c.replace('biolink:', '') for c in api_objects if len(c) < 50}
                        filtered_out[api] = f"Objects: {supported if len(supported) <= 5 else f'{len(api_objects)} types'}"
                else:
                    # Show all requested object categories in the message
                    subj_str = subject_categories[0].replace('biolink:', '')
                    obj_strs = [c.replace('biolink:', '') for c in object_categories]
                    if len(obj_strs) == 1:
                        obj_str = obj_strs[0]
                    else:
                        obj_str = f"[{'/'.join(obj_strs)}]"
                    filtered_out[api] = f"No {subj_str}→{obj_str} entries"
                continue

            # Check if API supports any of the requested predicates
            if predicates:
                api_predicates = set(transition_entries["Predicate"].unique())
                matching_predicates = api_predicates & set(predicates)

                if not matching_predicates:
                    # API has no matching predicates for this transition
                    sample_preds = list(api_predicates)[:3]
                    sample_str = ", ".join(p.replace("biolink:", "") for p in sample_preds)
                    filtered_out[api] = f"No matching predicates (has: {sample_str}...)"
                    continue

            filtered_apis.append(api)

        return filtered_apis, filtered_out

    def _extract_server_error(self, response_json: Dict[str, Any]) -> Optional[Dict[str, str]]:
        """Extract server-side errors/timeouts from TRAPI response.

        Translator APIs may return HTTP 200 but include error information
        in the response body when internal components fail or timeout.

        Args:
            response_json: Full JSON response from TRAPI API

        Returns:
            Dict with "message" and "type" if error found, None otherwise.
            type is one of: "server_timeout", "server_error"
        """
        # Common patterns for server-side timeout messages
        timeout_patterns = ["timed out", "timeout", "exceeded time limit"]

        # Check logs/status fields commonly used by Translator APIs
        error_sources = []

        # ARAX-style: check "logs" array for error/warning messages
        logs = response_json.get("logs", [])
        if isinstance(logs, list):
            for log in logs:
                if isinstance(log, dict):
                    level = log.get("level", "").lower()
                    message = log.get("message", "")
                    if level in ("error", "warning") and message:
                        error_sources.append(message)

        # Check "status" field (some APIs use this)
        status = response_json.get("status")
        if isinstance(status, str) and status.lower() not in ("success", "ok", "done"):
            error_sources.append(status)

        # Check message-level status
        message_status = response_json.get("message", {}).get("status")
        if isinstance(message_status, str) and message_status.lower() not in ("success", "ok", "done", ""):
            error_sources.append(message_status)

        # Check for description field with errors
        description = response_json.get("description", "")
        if isinstance(description, str) and any(word in description.lower() for word in ["error", "fail", "timeout"]):
            error_sources.append(description)

        if not error_sources:
            return None

        # Combine all error messages
        combined_message = "; ".join(error_sources[:3])  # Limit to first 3

        # Determine if it's a timeout or general error
        is_timeout = any(
            pattern in combined_message.lower()
            for pattern in timeout_patterns
        )

        return {
            "message": combined_message,
            "type": "server_timeout" if is_timeout else "server_error"
        }

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
        error_type = None
        result = None
        edge_count = 0

        try:
            response = requests.post(api_url, json=query_json_cur, timeout=self.timeout)
            duration = time.time() - start_time

            if response.status_code == 200:
                response_json = response.json()
                data = response_json.get("message", {})
                kg = data.get("knowledge_graph", {})
                edges = kg.get("edges", {})

                # Check for server-side errors/timeouts in the response
                server_error = self._extract_server_error(response_json)
                if server_error:
                    error_msg = server_error["message"]
                    error_type = server_error["type"]
                    if error_type == "server_timeout":
                        logger.info(f"API '{api_name}' reported internal timeout in {duration:.2f}s: {error_msg}")
                    else:
                        logger.info(f"API '{api_name}' reported server error in {duration:.2f}s: {error_msg}")

                if edges:
                    edge_count = len(edges)
                    result = data
                    logger.info(f"API '{api_name}' completed in {duration:.2f}s ({edge_count} edges)")
                elif "knowledge_graph" in data:
                    # API returned but with no edges (matches TCT behavior)
                    logger.info(f"API '{api_name}' completed in {duration:.2f}s (0 edges)")
                else:
                    if not error_msg:  # Only log if we haven't already logged a server error
                        logger.info(f"API '{api_name}' completed in {duration:.2f}s (no KG)")
            else:
                error_msg = f"HTTP {response.status_code}"
                error_type = "http_error"
                logger.info(f"API '{api_name}' failed in {duration:.2f}s ({error_msg})")

        except requests.Timeout:
            duration = time.time() - start_time
            error_msg = f"Client timeout after {self.timeout}s"
            error_type = "client_timeout"
            logger.info(f"API '{api_name}' client timeout after {duration:.2f}s (limit: {self.timeout}s)")
        except Exception as e:
            duration = time.time() - start_time
            error_msg = str(e)
            error_type = "exception"
            logger.info(f"API '{api_name}' error in {duration:.2f}s: {error_msg}")

        timing = APITiming(
            api_name=api_name,
            duration_seconds=round(duration, 3),
            success=result is not None,
            edge_count=edge_count,
            error=error_msg,
            error_type=error_type,
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
                        error_type="exception",
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

    def _log_api_query_summary(
        self,
        query_id: str,
        api_timings: List[APITiming],
        filtered_out_apis: Dict[str, str],
        intermediate_categories: Optional[List[str]],
        log_file: Optional[Path] = None,
    ) -> None:
        """Log API query results summary to file.

        Args:
            query_id: Unique identifier for this query
            api_timings: List of APITiming results
            filtered_out_apis: Dict of APIs filtered out and reasons
            intermediate_categories: Categories requested
            log_file: Path to log file (default: data/logs/api_queries_{date}.log)
        """
        from geneset_translator.config.settings import get_settings

        if log_file is None:
            settings = get_settings()
            log_file = settings.logs_dir / f"api_queries_{datetime.now().strftime('%Y%m%d')}.log"

        log_file.parent.mkdir(parents=True, exist_ok=True)

        with open(log_file, "a") as f:
            f.write(f"\n{'='*70}\n")
            f.write(f"Query ID: {query_id}\n")
            f.write(f"Timestamp: {datetime.now().isoformat()}\n")
            f.write(f"Intermediate Categories: {intermediate_categories}\n")
            f.write(f"{'='*70}\n\n")

            # APIs that were filtered out before querying
            if filtered_out_apis:
                f.write(f"FILTERED OUT ({len(filtered_out_apis)} APIs):\n")
                f.write("-" * 50 + "\n")
                for api, reason in sorted(filtered_out_apis.items()):
                    f.write(f"  {api}: {reason}\n")
                f.write("\n")

            # API query results
            f.write(f"QUERIED ({len(api_timings)} APIs):\n")
            f.write("-" * 50 + "\n")
            f.write(f"{'API Name':<45} {'Status':<12} {'Edges':>8} {'Duration':>10}\n")
            f.write("-" * 50 + "\n")

            successful = 0
            total_edges = 0
            for timing in sorted(api_timings, key=lambda t: t.edge_count, reverse=True):
                status = "SUCCESS" if timing.success else (timing.error_type or "FAILED").upper()[:12]
                f.write(f"{timing.api_name[:45]:<45} {status:<12} {timing.edge_count:>8} {timing.duration_seconds:>9.2f}s\n")

                if timing.success:
                    successful += 1
                    total_edges += timing.edge_count

            f.write("-" * 50 + "\n")
            f.write(f"Summary: {successful}/{len(api_timings)} APIs succeeded, {total_edges} total edges\n")
            f.write(f"{'='*70}\n")

        logger.info(f"API query summary logged to {log_file}")

    def normalize_genes(self, gene_symbols: List[str]) -> Dict[str, str]:
        """Normalize gene symbols to CURIEs using TCT name_resolver.

        Args:
            gene_symbols: List of HUGO gene symbols (e.g., ['CD6', 'IFITM3'])

        Returns:
            Dictionary mapping gene symbols to CURIEs
        """
        logger.info(f"Normalizing {len(gene_symbols)} genes using TCT name_resolver...")

        results = {}
        batch_result = name_resolver.batch_lookup(gene_symbols)

        for gene_symbol in gene_symbols:
            if gene_symbol in batch_result:
                info = batch_result[gene_symbol]
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
        """Execute pathfinder query: Gene → [Intermediate] → Disease.

        Discovers connections between genes of interest and a disease through
        biological intermediates (proteins, chemicals, other genes, etc.).

        Args:
            gene_symbols: List of gene symbols to query
            disease_curie: Disease CURIE (e.g., 'MONDO:0004975' for Alzheimer's)
            intermediate_categories: List of biolink categories for intermediate nodes
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

            # Filter to APIs that support exact requested categories and predicates
            logger.info(f"TCT.select_API returned {len(selected_APIs)} APIs")
            filtered_APIs, filtered_out = self._filter_apis_by_exact_categories(
                selected_apis=selected_APIs,
                subject_categories=input_gene_categories,
                object_categories=intermediate_categories,
                predicates=hop1_predicates,
                is_multihop=True,
            )

            # Fallback: if filtering removes all APIs, use original selection
            if not filtered_APIs:
                logger.warning(
                    f"API filtering removed all {len(selected_APIs)} APIs. "
                    f"Falling back to original selection."
                )
                filtered_APIs = selected_APIs
                filtered_out = {}
            else:
                logger.info(
                    f"Filtered to {len(filtered_APIs)} APIs "
                    f"(removed {len(filtered_out)} without matching categories/predicates)"
                )

            selected_APIs = filtered_APIs

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

            # Filter to APIs that support exact requested categories and predicates
            logger.info(f"TCT.select_API returned {len(selected_APIs)} APIs")
            filtered_APIs, filtered_out = self._filter_apis_by_exact_categories(
                selected_apis=selected_APIs,
                subject_categories=input_gene_categories,
                object_categories=target_categories,
                predicates=predicates,
                is_multihop=False,
            )

            # Fallback: if filtering removes all APIs, use original selection
            if not filtered_APIs:
                logger.warning(
                    f"API filtering removed all {len(selected_APIs)} APIs. "
                    f"Falling back to original selection."
                )
                filtered_APIs = selected_APIs
                filtered_out = {}
            else:
                logger.info(
                    f"Filtered to {len(filtered_APIs)} APIs "
                    f"(removed {len(filtered_out)} without matching categories/predicates)"
                )

            selected_APIs = filtered_APIs

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

        # Log API query summary to file
        query_id = datetime.now().strftime("%Y%m%d_%H%M%S")
        self._log_api_query_summary(
            query_id=query_id,
            api_timings=api_timings,
            filtered_out_apis=filtered_out,
            intermediate_categories=intermediate_categories,
        )

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
                "filtered_out_apis": filtered_out,  # APIs filtered out by exact category matching
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
