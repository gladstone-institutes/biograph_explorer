"""Publication analysis utilities for GeneSet Translator."""
import re
from typing import Dict, List, Any, Optional
import networkx as nx

# Regex patterns for publication IDs
PMID_PATTERN = re.compile(r'PMID:?\s*(\d+)', re.IGNORECASE)
PMC_PATTERN = re.compile(r'(?:PMC|PMCID):?\s*(PMC)?(\d+)', re.IGNORECASE)


def get_publication_frequency(graph: nx.DiGraph | nx.MultiDiGraph) -> Dict[str, int]:
    """Count publication occurrences across all edges.

    Args:
        graph: NetworkX graph (DiGraph or MultiDiGraph) with edge 'publications' attributes

    Returns:
        Dictionary mapping normalized publication ID to edge count
    """
    pub_counts: Dict[str, int] = {}
    for u, v, data in graph.edges(data=True):
        pubs = data.get('publications', [])
        if pubs:
            for pub in pubs:
                normalized = normalize_publication_id(pub)
                if normalized:
                    pub_counts[normalized] = pub_counts.get(normalized, 0) + 1
    return pub_counts


def get_publication_frequency_by_category(
    graph: nx.DiGraph | nx.MultiDiGraph
) -> Dict[str, Dict[str, int]]:
    """Count publication occurrences by intermediate node category.

    For each publication, counts how many edges cite it broken down by the
    category of the intermediate node (non-query-gene endpoint).

    Args:
        graph: NetworkX graph with edge 'publications' and node 'category' attributes

    Returns:
        Nested dict: {pub_id: {category: count}}
        Example: {"PMID:12345": {"Protein": 3, "Gene": 2}}
    """
    pub_category_counts: Dict[str, Dict[str, int]] = {}

    # Handle both DiGraph and MultiDiGraph
    if isinstance(graph, nx.MultiDiGraph):
        edges_iter = graph.edges(keys=True, data=True)
    else:
        edges_iter = ((u, v, None, d) for u, v, d in graph.edges(data=True))

    for u, v, _key, data in edges_iter:
        pubs = data.get('publications', [])
        if not pubs:
            continue

        # Determine intermediate category (non-query-gene node)
        # Check target first, then source
        target_is_query = graph.nodes.get(v, {}).get('is_query_gene', False)
        source_is_query = graph.nodes.get(u, {}).get('is_query_gene', False)

        if not target_is_query:
            category = graph.nodes.get(v, {}).get('category', 'Other')
        elif not source_is_query:
            category = graph.nodes.get(u, {}).get('category', 'Other')
        else:
            # Both are query genes - use 'Gene'
            category = 'Gene'

        for pub in pubs:
            normalized = normalize_publication_id(pub)
            if normalized:
                if normalized not in pub_category_counts:
                    pub_category_counts[normalized] = {}
                pub_category_counts[normalized][category] = (
                    pub_category_counts[normalized].get(category, 0) + 1
                )

    return pub_category_counts


def normalize_publication_id(pub_id: str) -> Optional[str]:
    """Normalize publication ID format for consistent grouping.

    Args:
        pub_id: Raw publication ID string

    Returns:
        Normalized publication ID or None if invalid
    """
    if not pub_id or not isinstance(pub_id, str):
        return None
    pub_id = pub_id.strip()

    if not pub_id:
        return None

    # Already normalized PMID format
    if pub_id.upper().startswith('PMID:'):
        return pub_id.upper()

    # PMC/PMCID formats - normalize to PMC:XXXXXXX
    if pub_id.upper().startswith(('PMC:', 'PMCID:')):
        match = PMC_PATTERN.match(pub_id)
        if match:
            return f"PMC:{match.group(2)}"
        return pub_id.upper()

    # URLs - keep as-is
    if pub_id.startswith(('http://', 'https://')):
        return pub_id

    return pub_id


def format_publication_display(pub_id: str) -> str:
    """Format publication ID for UI display.

    Args:
        pub_id: Normalized publication ID

    Returns:
        Human-readable display string
    """
    if not pub_id:
        return "Unknown"

    if pub_id.startswith('PMID:'):
        return f"PMID {pub_id[5:]}"
    if pub_id.startswith('PMC:'):
        return f"PMC {pub_id[4:]}"
    if pub_id.startswith(('http://', 'https://')):
        # Extract filename or last path segment
        parts = pub_id.rstrip('/').split('/')
        return parts[-1][:30] if parts else pub_id[:30]

    return pub_id[:40]


def _extract_publications_from_attributes(attributes: List[Dict[str, Any]]) -> List[str]:
    """Extract publications using the same logic as GraphBuilder._extract_publications_robust().

    This mirrors the extraction logic to validate what SHOULD be extracted.

    Args:
        attributes: List of TRAPI attribute dictionaries

    Returns:
        List of publication IDs that should be extracted
    """
    pubs = []

    for attr in attributes:
        # Pattern 1 & 2: Top-level publications
        if attr.get('attribute_type_id') == 'biolink:publications':
            value = attr.get('value', [])
            if isinstance(value, list):
                pubs.extend(value)
            elif value:
                pubs.append(value)

        # Pattern 3: Nested in has_supporting_study_result
        if attr.get('attribute_type_id') == 'biolink:has_supporting_study_result':
            for nested in attr.get('attributes', []):
                if nested.get('attribute_type_id') == 'biolink:publications':
                    value = nested.get('value', [])
                    if isinstance(value, list):
                        pubs.extend(value)
                    elif value:
                        pubs.append(value)

    return list(set(filter(None, pubs)))


def validate_publication_extraction(edges: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Detect duplicate edges in TRAPI data.

    With MultiDiGraph, edges with different predicates between the same nodes are
    preserved (they have different keys). This function validates that no edges
    with identical (subject, object, predicate) tuples exist, which would indicate
    duplicate data from TRAPI.

    Note: With the current edge key strategy (predicate + index), even true
    duplicates are preserved. This validation helps identify upstream data issues.

    Args:
        edges: List of raw TRAPI edge dictionaries

    Returns:
        Validation report dictionary with:
        - total_edges: Total raw TRAPI edges
        - edges_with_publications: Edges that have extractable publications
        - unique_node_pairs: Number of unique (subject, object, predicate) tuples
        - edges_lost_to_collisions: Count of duplicate edges (same S,O,P tuple)
        - publications_at_risk: Publications on duplicate edges
        - sample_collisions: Sample duplicate edge groups
        - has_issues: Boolean indicating if duplicates were found
    """
    total_edges = len(edges)
    edges_with_pubs = 0

    # Track edges by (subject, object, predicate) to detect true duplicates
    # With MultiDiGraph, edges with different predicates are preserved (different keys)
    # Only edges with identical (S, O, P) would be duplicates
    edge_tuples: Dict[tuple, List[Dict[str, Any]]] = {}

    for edge in edges:
        attributes = edge.get('attributes', [])
        extractable_pubs = _extract_publications_from_attributes(attributes)

        if extractable_pubs:
            edges_with_pubs += 1

        # Group by (subject, object, predicate) - matches MultiDiGraph key behavior
        edge_tuple = (edge.get('subject'), edge.get('object'), edge.get('predicate'))
        if edge_tuple not in edge_tuples:
            edge_tuples[edge_tuple] = []
        edge_tuples[edge_tuple].append({
            'predicate': edge.get('predicate'),
            'publications': extractable_pubs,
            'subject': edge.get('subject'),
            'object': edge.get('object'),
        })

    # Analyze duplicates (same subject, object, predicate tuple appearing multiple times)
    unique_tuples = len(edge_tuples)
    duplicate_groups = {k: v for k, v in edge_tuples.items() if len(v) > 1}
    duplicate_count = sum(len(v) - 1 for v in duplicate_groups.values())

    # Find publications on duplicate edges
    # Note: With predicate+index edge keys, all edges are preserved in the graph,
    # but duplicates may indicate upstream data issues worth flagging
    publications_on_duplicates = 0
    sample_duplicates = []

    for edge_tuple, edge_list in duplicate_groups.items():
        # Count publications on duplicate edges (all but first occurrence)
        for edge_info in edge_list[1:]:
            publications_on_duplicates += len(edge_info['publications'])

        # Collect samples
        if len(sample_duplicates) < 5:
            pubs_by_edge = [
                f"{e['predicate']}: {len(e['publications'])} pubs"
                for e in edge_list if e['publications']
            ]
            if pubs_by_edge:  # Only include if there are publications involved
                sample_duplicates.append({
                    'subject': edge_list[0]['subject'],
                    'object': edge_list[0]['object'],
                    'predicate': edge_list[0]['predicate'],
                    'duplicate_count': len(edge_list),
                    'pubs_by_edge': pubs_by_edge,
                })

    return {
        'total_edges': total_edges,
        'edges_with_publications': edges_with_pubs,
        'unique_node_pairs': unique_tuples,
        'edges_lost_to_collisions': duplicate_count,  # Kept for backward compatibility
        'collision_pairs_count': len(duplicate_groups),
        'publications_at_risk': publications_on_duplicates,  # Kept for backward compatibility
        'sample_collisions': sample_duplicates,  # Kept for backward compatibility
        'has_issues': duplicate_count > 0,
    }
