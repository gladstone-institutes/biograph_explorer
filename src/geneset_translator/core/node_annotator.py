"""Node Annotator API client for enriching graph nodes with metadata.

This module provides a client for the NCATS Translator Node Annotator API
(https://annotator.transltr.io/) to fetch annotations for graph nodes.

Follows patterns from TRAPIClient:
- Parallel batch queries
- Response caching
- Graceful degradation on API failures
- Progress callbacks for UI integration
"""

import json
import logging
import time
import urllib.parse
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

import networkx as nx
import requests

from .translator_node import TranslatorAttribute, TranslatorNode

logger = logging.getLogger(__name__)

# Node Annotator API base URL
NODE_ANNOTATOR_URL = "https://annotator.transltr.io/"


@dataclass
class AnnotationTiming:
    """Timing information for a batch annotation request."""

    batch_id: int
    curie_count: int
    duration_seconds: float
    success: bool
    results_count: int = 0
    error: Optional[str] = None


class NodeAnnotator:
    """Client for NCATS Translator Node Annotator API.

    Enriches graph nodes with metadata from the Node Annotator API,
    including gene types, chromosome locations, organism information,
    and other annotations.

    Example:
        >>> annotator = NodeAnnotator(cache_dir=Path("data/cache/annotations"))
        >>> graph, metadata = annotator.annotate_graph(graph)
        >>> print(f"Annotated {metadata['annotated_count']} nodes")
    """

    def __init__(
        self,
        cache_dir: Path = Path("data/cache/annotations"),
        timeout: int = 120,
        batch_size: int = 500,
        max_workers: int = 4,
    ):
        """Initialize Node Annotator client.

        Args:
            cache_dir: Directory for caching annotation responses.
            timeout: API timeout in seconds per request.
            batch_size: Maximum CURIEs per batch (API limit ~1000).
            max_workers: Maximum parallel batch requests.
        """
        self.cache_dir = cache_dir
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.timeout = timeout
        self.batch_size = batch_size
        self.max_workers = max_workers
        self.base_url = NODE_ANNOTATOR_URL

    def status(self) -> dict:
        """Check Node Annotator API status.

        Returns:
            Status response from the API.
        """
        response = requests.get(f"{self.base_url}status", timeout=10)
        response.raise_for_status()
        return response.json()

    def annotate_nodes(
        self,
        node_curies: List[str],
        fields: str = "all",
        progress_callback: Optional[Callable[[str], None]] = None,
    ) -> Dict[str, TranslatorNode]:
        """Fetch annotations for multiple node CURIEs.

        Args:
            node_curies: List of CURIEs to annotate.
            fields: "all" or comma-separated field names.
            progress_callback: Optional callback for progress updates.

        Returns:
            Dict mapping CURIEs to TranslatorNode objects.
        """
        if not node_curies:
            return {}

        logger.info(f"Annotating {len(node_curies)} nodes...")
        if progress_callback:
            progress_callback(f"Starting annotation of {len(node_curies)} nodes")

        # Execute parallel batch requests
        results, timings = self._parallel_batch_annotate(
            node_curies, fields, progress_callback
        )

        # Log summary
        successful = sum(1 for t in timings if t.success)
        total_time = sum(t.duration_seconds for t in timings)
        logger.info(
            f"Annotation complete: {len(results)} nodes annotated "
            f"({successful}/{len(timings)} batches succeeded, {total_time:.2f}s total)"
        )

        return results

    def _extract_annotation_features(
        self,
        translator_node: TranslatorNode,
    ) -> Dict[str, Any]:
        """Extract flat, clustering-ready features from a TranslatorNode.

        Returns a dict with only scalar values or lists of scalars (strings).
        All nested structures are flattened to simple key-value pairs.
        This ensures the features are JSON-serializable and suitable for
        filtering, clustering, and Cytoscape display.

        Args:
            translator_node: TranslatorNode with raw API attributes.

        Returns:
            Dict with flattened features (scalars and lists of strings only).
        """
        features: Dict[str, Any] = {}

        if not translator_node.attributes:
            return features

        for attr in translator_node.attributes:
            key = attr.attribute_type_id
            value = attr.value

            # Handle scalar values directly
            if isinstance(value, (str, int, float, bool)):
                features[key] = value

            # Handle GO terms (nested structure → flat lists of term names)
            elif key == 'go' and isinstance(value, dict):
                for ontology in ['BP', 'MF', 'CC']:
                    terms = value.get(ontology, [])
                    if terms and isinstance(terms, list):
                        term_names = []
                        for term_entry in terms:
                            if isinstance(term_entry, dict):
                                term_name = term_entry.get('term')
                                if term_name:
                                    term_names.append(term_name)
                        if term_names:
                            features[f'go_{ontology.lower()}'] = term_names

            # Handle lists of scalars (e.g., aliases)
            elif isinstance(value, list) and value:
                # Check if all items are scalars
                if all(isinstance(v, (str, int, float)) for v in value):
                    features[key] = [str(v) for v in value]
                # Skip lists of complex objects (not suitable for features)

            # Skip complex objects (dicts, nested structures) - not suitable for features

        return features

    def annotate_graph(
        self,
        graph: nx.DiGraph,
        fields: str = "all",
        progress_callback: Optional[Callable[[str], None]] = None,
    ) -> Tuple[nx.DiGraph, Dict[str, Any]]:
        """Annotate all nodes in a NetworkX graph.

        Adds annotation attributes to each node:
        - translator_node: TranslatorNode object with full annotations
        - annotation_features: Dict of flattened features (for filtering/clustering)

        The annotation_features dict contains only scalars and lists of strings,
        making it safe for JSON serialization and Cytoscape display.

        Args:
            graph: NetworkX DiGraph to annotate.
            fields: "all" or comma-separated field names.
            progress_callback: Optional callback for progress updates.

        Returns:
            Tuple of (annotated graph, metadata dict).
        """
        node_curies = list(graph.nodes())
        logger.info(f"Annotating graph with {len(node_curies)} nodes")

        # Fetch annotations
        annotations = self.annotate_nodes(node_curies, fields, progress_callback)

        # Add annotations to graph nodes and collect all attribute keys
        annotated_count = 0
        all_attribute_keys: Dict[str, set] = {}  # key -> set of unique values

        for curie in node_curies:
            if curie in annotations:
                translator_node = annotations[curie]

                # Store raw TranslatorNode (excluded from Cytoscape serialization)
                graph.nodes[curie]['translator_node'] = translator_node

                # Store flattened features (for clustering + filtering)
                features = self._extract_annotation_features(translator_node)
                graph.nodes[curie]['annotation_features'] = features

                # Collect unique values for metadata
                for key, value in features.items():
                    if key not in all_attribute_keys:
                        all_attribute_keys[key] = set()
                    if isinstance(value, list):
                        all_attribute_keys[key].update(value)
                    else:
                        all_attribute_keys[key].add(value)

                annotated_count += 1

        # Debug: log all discovered attribute keys
        logger.debug(f"All discovered attribute keys: {list(all_attribute_keys.keys())}")
        for key in ['go_bp', 'go_mf', 'go_cc']:
            if key in all_attribute_keys:
                logger.debug(f"  {key}: {len(all_attribute_keys[key])} unique values")

        # Build metadata with discovered attributes
        # List-valued attributes (GO terms) stored separately for different UI handling
        list_valued_keys = {'go_bp', 'go_mf', 'go_cc', 'alias'}

        # Scalar attributes with 2-100 unique values (good for dropdown filtering)
        # Increased upper bound from 50 to 100 to capture more useful attributes
        filterable_attributes = {}
        for key, values in all_attribute_keys.items():
            if key not in list_valued_keys and 2 <= len(values) <= 100:
                filterable_attributes[key] = sorted(str(v) for v in values)

        # Always include type_of_gene if present (useful for gene type filtering)
        if 'type_of_gene' in all_attribute_keys and 'type_of_gene' not in filterable_attributes:
            filterable_attributes['type_of_gene'] = sorted(str(v) for v in all_attribute_keys['type_of_gene'])

        # List-valued attributes (GO terms) - include all for searchable filtering
        searchable_attributes = {}
        for key in list_valued_keys:
            if key in all_attribute_keys and len(all_attribute_keys[key]) >= 2:
                searchable_attributes[key] = sorted(str(v) for v in all_attribute_keys[key])

        metadata = {
            'total_nodes': len(node_curies),
            'annotated_count': annotated_count,
            'annotation_rate': annotated_count / len(node_curies) if node_curies else 0,
            'all_attribute_keys': sorted(all_attribute_keys.keys()),
            'filterable_attributes': filterable_attributes,
            'searchable_attributes': searchable_attributes,
            'timestamp': datetime.now().isoformat(),
        }

        logger.info(
            f"Graph annotation complete: {annotated_count}/{len(node_curies)} nodes "
            f"({metadata['annotation_rate']:.1%}), "
            f"{len(filterable_attributes)} filterable attributes, "
            f"{len(searchable_attributes)} GO term attributes (searchable)"
        )

        return graph, metadata

    def annotate_with_hpa(
        self,
        graph: nx.DiGraph,
        progress_callback: Optional[Callable[[str], None]] = None,
    ) -> Tuple[nx.DiGraph, Dict[str, Any]]:
        """Add Human Protein Atlas cell type specificity annotations to gene nodes.

        This is an optional enrichment step that adds HPA single-cell RNA data
        to gene nodes. Should be called after annotate_graph().

        Args:
            graph: NetworkX DiGraph with gene nodes (already annotated).
            progress_callback: Optional callback for progress updates.

        Returns:
            Tuple of (annotated graph, HPA metadata dict).
        """
        from .hpa_client import HPAClient

        hpa_client = HPAClient(
            cache_dir=self.cache_dir.parent / "hpa",
            timeout=30,
            max_workers=4,
        )

        return hpa_client.annotate_graph(graph, progress_callback)

    def _batch_annotate(
        self,
        batch_id: int,
        curies: List[str],
        fields: str,
    ) -> Tuple[Dict[str, TranslatorNode], AnnotationTiming]:
        """Annotate a single batch of CURIEs.

        Args:
            batch_id: Batch identifier for logging.
            curies: List of CURIEs to annotate (max batch_size).
            fields: Fields parameter for API.

        Returns:
            Tuple of (results dict, timing info).
        """
        start_time = time.time()
        results = {}
        error_msg = None

        try:
            path = urllib.parse.urljoin(self.base_url, 'curie')
            response = requests.post(
                path,
                json={'ids': curies},
                params={'fields': fields} if fields != "all" else {},
                timeout=self.timeout,
            )
            response.raise_for_status()

            raw_results = response.json()
            duration = time.time() - start_time

            # Parse results to TranslatorNode objects
            for curie, data in raw_results.items():
                translator_node = self._parse_annotation_response(curie, data)
                if translator_node:
                    results[curie] = translator_node

            logger.debug(
                f"Batch {batch_id}: {len(curies)} CURIEs → {len(results)} results "
                f"in {duration:.2f}s"
            )

        except requests.Timeout:
            duration = time.time() - start_time
            error_msg = "Timeout"
            logger.warning(f"Batch {batch_id}: Timeout after {duration:.2f}s")

        except requests.HTTPError as e:
            duration = time.time() - start_time
            error_msg = f"HTTP {e.response.status_code}"
            logger.warning(f"Batch {batch_id}: {error_msg}")

        except Exception as e:
            duration = time.time() - start_time
            error_msg = str(e)
            logger.warning(f"Batch {batch_id}: Error - {error_msg}")

        timing = AnnotationTiming(
            batch_id=batch_id,
            curie_count=len(curies),
            duration_seconds=round(duration, 3),
            success=error_msg is None,
            results_count=len(results),
            error=error_msg,
        )

        return results, timing

    def _parallel_batch_annotate(
        self,
        curies: List[str],
        fields: str,
        progress_callback: Optional[Callable[[str], None]] = None,
    ) -> Tuple[Dict[str, TranslatorNode], List[AnnotationTiming]]:
        """Execute parallel batch annotation requests.

        Args:
            curies: List of all CURIEs to annotate.
            fields: Fields parameter for API.
            progress_callback: Optional progress callback.

        Returns:
            Tuple of (merged results dict, list of timings).
        """
        # Split into batches
        batches = [
            curies[i:i + self.batch_size]
            for i in range(0, len(curies), self.batch_size)
        ]
        total_batches = len(batches)

        logger.info(f"Processing {len(curies)} CURIEs in {total_batches} batches")

        all_results = {}
        all_timings = []
        completed = 0

        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            future_to_batch = {
                executor.submit(
                    self._batch_annotate, i, batch, fields
                ): i
                for i, batch in enumerate(batches)
            }

            for future in as_completed(future_to_batch):
                batch_id = future_to_batch[future]
                try:
                    results, timing = future.result()
                    all_results.update(results)
                    all_timings.append(timing)
                except Exception as e:
                    logger.error(f"Batch {batch_id} raised exception: {e}")
                    all_timings.append(AnnotationTiming(
                        batch_id=batch_id,
                        curie_count=len(batches[batch_id]),
                        duration_seconds=0.0,
                        success=False,
                        error=str(e),
                    ))

                completed += 1
                if progress_callback:
                    progress_callback(
                        f"Annotated {completed}/{total_batches} batches "
                        f"({len(all_results)} nodes)"
                    )

        # Sort timings by batch_id
        all_timings.sort(key=lambda t: t.batch_id)

        return all_results, all_timings

    def _parse_annotation_response(
        self,
        curie: str,
        node_data: Any,
    ) -> Optional[TranslatorNode]:
        """Parse Node Annotator API response to TranslatorNode.

        Args:
            curie: The CURIE being annotated.
            node_data: Raw response data for this CURIE.

        Returns:
            TranslatorNode object, or None if parsing fails.
        """
        if not node_data:
            return None

        # Unwrap list if needed (API sometimes returns [data] instead of data)
        if isinstance(node_data, list):
            if not node_data:
                return None
            node_data = node_data[0]

        if not isinstance(node_data, dict):
            return None

        # Check for "not found" marker
        if node_data.get('notfound'):
            return None

        # Extract label
        label = (
            node_data.get('name') or
            node_data.get('symbol') or
            node_data.get('label') or
            curie
        )

        # Extract types
        types = []
        if 'type_of_gene' in node_data:
            types.append(f"gene_type:{node_data['type_of_gene']}")

        # Extract synonyms
        synonyms = node_data.get('alias', [])
        if isinstance(synonyms, str):
            synonyms = [synonyms]

        # Extract taxa
        taxa = []
        if 'taxid' in node_data:
            taxa.append(f"NCBITaxon:{node_data['taxid']}")

        # Convert all fields to TranslatorAttributes
        attributes = []
        skip_keys = {'_id', '_score', 'query', 'notfound'}
        for key, value in node_data.items():
            if key not in skip_keys:
                attributes.append(TranslatorAttribute(
                    attribute_type_id=key,
                    value=value,
                ))

        return TranslatorNode(
            curie=curie,
            label=label,
            types=types if types else None,
            synonyms=synonyms if synonyms else None,
            attributes=attributes if attributes else None,
            taxa=taxa if taxa else None,
        )

    def save_annotations_cache(
        self,
        annotations: Dict[str, TranslatorNode],
        filename: Optional[str] = None,
    ) -> Path:
        """Save annotations to cache file.

        Args:
            annotations: Dict of CURIE -> TranslatorNode.
            filename: Optional filename (auto-generated if not provided).

        Returns:
            Path to saved cache file.
        """
        if filename is None:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            filename = f"node_annotations_{timestamp}.json"

        cache_path = self.cache_dir / filename

        # Convert to serializable format
        cache_data = {}
        for curie, node in annotations.items():
            cache_data[curie] = {
                'curie': node.curie,
                'label': node.label,
                'types': node.types,
                'synonyms': node.synonyms,
                'taxa': node.taxa,
                'attributes': [
                    {
                        'attribute_type_id': attr.attribute_type_id,
                        'value': attr.value,
                    }
                    for attr in (node.attributes or [])
                ],
            }

        with open(cache_path, 'w') as f:
            json.dump(cache_data, f, indent=2, default=str)

        logger.info(f"Saved annotations cache to {cache_path}")
        return cache_path

    def load_annotations_cache(
        self,
        filename: str,
    ) -> Dict[str, TranslatorNode]:
        """Load annotations from cache file.

        Args:
            filename: Cache filename to load.

        Returns:
            Dict of CURIE -> TranslatorNode.
        """
        cache_path = self.cache_dir / filename

        with open(cache_path) as f:
            cache_data = json.load(f)

        annotations = {}
        for curie, data in cache_data.items():
            attributes = [
                TranslatorAttribute(
                    attribute_type_id=attr['attribute_type_id'],
                    value=attr['value'],
                )
                for attr in data.get('attributes', [])
            ]
            annotations[curie] = TranslatorNode(
                curie=data['curie'],
                label=data.get('label'),
                types=data.get('types'),
                synonyms=data.get('synonyms'),
                taxa=data.get('taxa'),
                attributes=attributes if attributes else None,
            )

        logger.info(f"Loaded {len(annotations)} annotations from {cache_path}")
        return annotations


def lookup_curie(curie: str, **kwargs) -> Dict[str, Any]:
    """Look up a single CURIE annotation.

    Convenience function for quick annotation lookups.

    Args:
        curie: CURIE to annotate.
        **kwargs: Additional arguments for the API.

    Returns:
        Annotation data dict for the CURIE.
    """
    return lookup_curies([curie], **kwargs).get(curie, {})


def lookup_curies(curies: List[str], **kwargs) -> Dict[str, Any]:
    """Look up annotations for multiple CURIEs.

    Convenience function that returns raw API response.

    Args:
        curies: List of CURIEs to annotate.
        **kwargs: Additional arguments (raw, fields, include_extra).

    Returns:
        Dict mapping CURIEs to annotation data.
    """
    path = urllib.parse.urljoin(NODE_ANNOTATOR_URL, 'curie')
    response = requests.post(path, json={'ids': curies, **kwargs})
    response.raise_for_status()

    results = response.json()
    if not results:
        raise LookupError(f'No matching CURIE found for: {curies}')

    # Unwrap single-item lists
    for curie in results:
        if isinstance(results[curie], list) and len(results[curie]) == 1:
            results[curie] = results[curie][0]

    return results
