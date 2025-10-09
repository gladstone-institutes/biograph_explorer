"""Network visualization using st-link-analysis component.

Features:
- Interactive Cytoscape.js graph rendering
- Multiple layout algorithms (cose, fcose, circle, grid, etc.)
- Node sizing by centrality or gene frequency
- Node coloring and icons by category
- Edge styling with predicates
- Material Icons for node types
- Built-in fullscreen and JSON export

Phase 2 Status: Implemented with st-link-analysis
"""

from typing import Optional, List, Dict, Any, Tuple
import networkx as nx
from networkx.readwrite import json_graph
from pathlib import Path
import logging
import base64

try:
    from st_link_analysis import NodeStyle, EdgeStyle
    ST_LINK_ANALYSIS_AVAILABLE = True
except ImportError:
    ST_LINK_ANALYSIS_AVAILABLE = False

logger = logging.getLogger(__name__)

# Color scheme for node categories
CATEGORY_COLORS = {
    "Gene": "#E74C3C",  # Red
    "Disease": "#9B59B6",  # Purple
    "Protein": "#4ECDC4",  # Cyan
    "ChemicalEntity": "#F39C12",  # Yellow/Orange
    "BiologicalProcess": "#2ECC71",  # Green
    "Cluster": "#16A085",  # Teal (for cluster meta-nodes)
    "Other": "#3498DB",  # Blue
}

# Get the absolute path to the assets directory
_ASSETS_DIR = Path(__file__).parent.parent / "assets"


def _svg_to_data_uri(svg_path: Path) -> str:
    """Convert SVG file to data URI for embedding in browser.

    Args:
        svg_path: Path to SVG file

    Returns:
        Data URI string (e.g., "url('data:image/svg+xml;base64,...')")
    """
    try:
        with open(svg_path, 'rb') as f:
            svg_data = f.read()
        encoded = base64.b64encode(svg_data).decode('utf-8')
        return f"url('data:image/svg+xml;base64,{encoded}')"
    except Exception as e:
        logger.warning(f"Failed to load SVG icon {svg_path}: {e}")
        return None


# SVG icons for node categories (converted to data URIs)
CATEGORY_ICONS = {
    "Gene": _svg_to_data_uri(_ASSETS_DIR / "genetics.svg"),
    "Disease": _svg_to_data_uri(_ASSETS_DIR / "sick.svg"),
    "ChemicalEntity": _svg_to_data_uri(_ASSETS_DIR / "experiment.svg"),
    "Protein": None,
    "BiologicalProcess": None,
    "Cluster": None,
    "Other": None,
}


def prepare_cytoscape_elements(
    graph: nx.DiGraph,
    query_genes: List[str],
    sizing_metric: str = "gene_frequency",
) -> Dict[str, List[Dict]]:
    """Convert NetworkX graph to Cytoscape.js elements format.

    Uses nx.cytoscape_data() as base, then enriches with custom attributes.

    Args:
        graph: NetworkX DiGraph to convert
        query_genes: List of query gene node IDs (for highlighting)
        sizing_metric: Metric for node sizing (gene_frequency, pagerank, betweenness, degree)

    Returns:
        Dictionary with "nodes" and "edges" lists in Cytoscape.js format
    """
    # Get base Cytoscape format from NetworkX
    cyto_data = json_graph.cytoscape_data(graph)
    elements = cyto_data["elements"]

    # Calculate metric values for sizing
    metric_values = {}
    if sizing_metric in ["pagerank", "betweenness", "degree"]:
        for node in graph.nodes():
            metric_values[node] = graph.nodes[node].get(sizing_metric, 0)
    elif sizing_metric == "gene_frequency":
        for node in graph.nodes():
            metric_values[node] = graph.nodes[node].get("gene_frequency", 0)
    else:
        # Default to degree
        metric_values = dict(graph.degree())

    # Normalize to 0-1 range
    if metric_values:
        max_val = max(metric_values.values()) if max(metric_values.values()) > 0 else 1
        normalized_metrics = {k: v / max_val for k, v in metric_values.items()}
    else:
        normalized_metrics = {node: 0.5 for node in graph.nodes()}

    # Enrich nodes with custom attributes
    for node_element in elements["nodes"]:
        node_id = node_element["data"]["id"]
        node_attrs = graph.nodes[node_id]

        # Add category as label (for NodeStyle matching)
        category = node_attrs.get("category", "Other")
        node_element["data"]["label"] = category  # Used by NodeStyle for matching

        # Add display name (for caption)
        original_symbol = node_attrs.get("original_symbol", "")
        display_label = node_attrs.get("label", node_id)
        node_element["data"]["name"] = original_symbol if original_symbol else display_label

        # Add query gene flag
        is_query_gene = node_attrs.get("is_query_gene", False)
        node_element["data"]["is_query_gene"] = is_query_gene

        # Calculate node size (15-45px range from PROJECT_PLAN.md)
        base_size = 15 + (normalized_metrics.get(node_id, 0.5) * 30)

        # Query genes get larger minimum size
        if is_query_gene:
            base_size = max(base_size, 25)

        node_element["data"]["size"] = base_size

        # Add metrics for tooltip/inspection
        node_element["data"]["gene_frequency"] = node_attrs.get("gene_frequency", 0)
        node_element["data"]["pagerank"] = node_attrs.get("pagerank", 0)
        node_element["data"]["betweenness"] = node_attrs.get("betweenness", 0)
        node_element["data"]["degree"] = graph.degree(node_id)

    # Enrich edges with custom attributes
    for idx, edge_element in enumerate(elements["edges"]):
        source = edge_element["data"]["source"]
        target = edge_element["data"]["target"]

        # Add unique ID for st-link-analysis (REQUIRED)
        edge_element["data"]["id"] = f"e{idx}"

        # Get edge data from graph
        if graph.has_edge(source, target):
            edge_attrs = graph[source][target]

            # Add predicate (cleaned)
            predicate = edge_attrs.get("predicate", "")
            edge_element["data"]["label"] = predicate.replace("biolink:", "")

            # Add knowledge sources
            sources = edge_attrs.get("knowledge_source", [])
            edge_element["data"]["knowledge_sources"] = sources

    return {"nodes": elements["nodes"], "edges": elements["edges"]}


def create_node_styles(graph: nx.DiGraph) -> List["NodeStyle"]:
    """Create NodeStyle objects for each category in the graph.

    Args:
        graph: NetworkX graph containing nodes with "category" attribute

    Returns:
        List of NodeStyle objects, one per unique category
    """
    if not ST_LINK_ANALYSIS_AVAILABLE:
        raise ImportError("st-link-analysis not installed. Run: pip install st-link-analysis")

    # Extract unique categories
    categories = set()
    for node in graph.nodes():
        category = graph.nodes[node].get("category", "Other")
        categories.add(category)

    # Create NodeStyle for each category
    node_styles = []
    for category in sorted(categories):
        color = CATEGORY_COLORS.get(category, CATEGORY_COLORS["Other"])
        icon = CATEGORY_ICONS.get(category)

        # Use "name" as caption (displays gene symbol or node label)
        # Only include icon parameter if icon is not None
        if icon:
            node_styles.append(
                NodeStyle(
                    label=category,
                    color=color,
                    caption="name",
                    icon=icon
                )
            )
        else:
            node_styles.append(
                NodeStyle(
                    label=category,
                    color=color,
                    caption="name"
                )
            )

    return node_styles


def create_edge_styles(graph: nx.DiGraph) -> List["EdgeStyle"]:
    """Create EdgeStyle objects for each predicate in the graph.

    Args:
        graph: NetworkX graph containing edges with "predicate" attribute

    Returns:
        List of EdgeStyle objects, one per unique predicate
    """
    if not ST_LINK_ANALYSIS_AVAILABLE:
        raise ImportError("st-link-analysis not installed. Run: pip install st-link-analysis")

    # Extract unique predicates
    predicates = set()
    for _, _, edge_attrs in graph.edges(data=True):
        predicate = edge_attrs.get("predicate", "").replace("biolink:", "")
        if predicate:
            predicates.add(predicate)

    # Create EdgeStyle for each predicate
    edge_styles = []
    for predicate in sorted(predicates):
        edge_styles.append(
            EdgeStyle(
                label=predicate,     # Matches edge data["label"]
                caption="label",     # Display the label attribute
                directed=True        # Show arrows
            )
        )

    # If no predicates found, add a default style
    if not edge_styles:
        edge_styles.append(
            EdgeStyle(
                label="default",     # Default label
                caption="label",     # Display the label attribute
                directed=True        # Show arrows
            )
        )

    return edge_styles


def get_layout_config(layout_name: str = "dagre") -> Dict[str, Any]:
    """Get layout configuration for Cytoscape.js.

    Args:
        layout_name: Name of layout algorithm (dagre, fcose, cola, cose-bilkent, breadthfirst, cose, circle, grid, concentric)

    Returns:
        Layout configuration dictionary
    """
    # Default configuration based on demos
    layout = {
        "name": layout_name,
        "animate": "end",
        "nodeDimensionsIncludeLabels": False
    }

    # Add layout-specific options
    if layout_name == "dagre":
        layout.update({
            "rankDir": "TB",  # Top to bottom (genes → intermediates → disease)
            "ranker": "network-simplex",  # Optimal ranking algorithm
            "nodeSep": 50,  # Separation between nodes on same rank
            "edgeSep": 10,  # Separation between edges
            "rankSep": 75,  # Separation between ranks (levels)
            "fit": True,
            "padding": 30
        })
    elif layout_name == "cose":
        layout.update({
            "nodeRepulsion": 400000,
            "idealEdgeLength": 100,
            "edgeElasticity": 100,
            "nestingFactor": 5,
            "gravity": 80,
            "numIter": 1000,
            "initialTemp": 200,
            "coolingFactor": 0.95,
            "minTemp": 1.0
        })
    elif layout_name == "fcose":
        layout.update({
            "quality": "default",
            "randomize": True,
            "animate": "end",
            "fit": True,
            "padding": 30,
            "nodeSeparation": 75,
            "idealEdgeLength": 50,
            "edgeElasticity": 0.45,
            "nestingFactor": 0.1,
            "gravity": 0.25,
            "numIter": 2500,
            "tile": True,
            "tilingPaddingVertical": 10,
            "tilingPaddingHorizontal": 10,
            "gravityRangeCompound": 1.5,
            "gravityCompound": 1.0,
            "gravityRange": 3.8
        })
    elif layout_name == "cola":
        layout.update({
            "animate": True,
            "refresh": 1,
            "maxSimulationTime": 4000,
            "ungrabifyWhileSimulating": False,
            "fit": True,
            "padding": 30,
            "nodeDimensionsIncludeLabels": False,
            "randomize": False,
            "avoidOverlap": True,
            "handleDisconnected": True,
            "convergenceThreshold": 0.01,
            "nodeSpacing": 10,
            "flow": None,
            "alignment": None,
            "gapInequalities": None
        })
    elif layout_name == "cose-bilkent":
        layout.update({
            "quality": "default",
            "nodeDimensionsIncludeLabels": False,
            "refresh": 30,
            "fit": True,
            "padding": 30,
            "randomize": True,
            "nodeRepulsion": 4500,
            "idealEdgeLength": 50,
            "edgeElasticity": 0.45,
            "nestingFactor": 0.1,
            "gravity": 0.25,
            "numIter": 2500,
            "tile": True,
            "tilingPaddingVertical": 10,
            "tilingPaddingHorizontal": 10,
            "gravityRangeCompound": 1.5,
            "gravityCompound": 1.0,
            "gravityRange": 3.8,
            "initialEnergyOnIncremental": 0.5
        })

    return layout


def render_network_visualization(
    graph: nx.DiGraph,
    query_genes: List[str],
    sizing_metric: str = "gene_frequency",
    layout: str = "dagre",
    max_intermediates: int = 200,
    highlight_nodes: Optional[List[str]] = None,
    disease_curie: Optional[str] = None,
) -> Dict[str, Any]:
    """Prepare network visualization data for st-link-analysis component.

    Args:
        graph: NetworkX DiGraph to visualize
        query_genes: List of query gene node IDs (highlighted)
        sizing_metric: Node sizing metric (gene_frequency, pagerank, betweenness, degree)
        layout: Layout algorithm name (cose, fcose, circle, grid, breadthfirst, concentric)
        max_intermediates: Maximum intermediate nodes to display (default: 200)
        highlight_nodes: Optional list of node IDs to highlight (for citations)
        disease_curie: Optional disease CURIE (always included if present)

    Returns:
        Dictionary with elements, node_styles, edge_styles, and layout config

    Example:
        >>> viz_data = render_network_visualization(graph, ["NCBIGene:1803"], max_intermediates=200)
        >>> st_link_analysis(viz_data["elements"], viz_data["layout"],
        ...                  viz_data["node_styles"], viz_data["edge_styles"])
    """
    if not ST_LINK_ANALYSIS_AVAILABLE:
        logger.error("st-link-analysis not installed")
        return None

    if graph.number_of_nodes() == 0:
        logger.warning("Empty graph - nothing to visualize")
        return None

    try:
        # Sample graph by top intermediates
        graph = sample_graph_for_visualization(
            graph,
            query_genes,
            max_intermediates=max_intermediates,
            disease_curie=disease_curie
        )

        # Prepare Cytoscape elements
        elements = prepare_cytoscape_elements(graph, query_genes, sizing_metric)

        # Create node and edge styles
        node_styles = create_node_styles(graph)
        edge_styles = create_edge_styles(graph)

        # Get layout configuration
        layout_config = get_layout_config(layout)

        return {
            "elements": elements,
            "node_styles": node_styles,
            "edge_styles": edge_styles,
            "layout": layout_config
        }

    except Exception as e:
        logger.error(f"Error creating visualization: {e}")
        import traceback
        traceback.print_exc()
        return None


def sample_graph_for_visualization(
    graph: nx.DiGraph,
    query_genes: List[str],
    max_intermediates: int = 200,
    disease_curie: Optional[str] = None,
) -> nx.DiGraph:
    """Sample graph by filtering top intermediate nodes ranked by query gene connectivity.

    Sampling strategy:
    1. Always include ALL query genes
    2. Identify intermediate nodes (non-query genes)
    3. Score intermediates by number of edges connecting to query genes
    4. Keep top-N intermediates based on max_intermediates parameter
    5. Include disease node if present
    6. Build subgraph with selected nodes + their edges

    Args:
        graph: Full NetworkX graph
        query_genes: Query gene node IDs (always included)
        max_intermediates: Maximum number of intermediate nodes to include (default: 200)
        disease_curie: Optional disease CURIE (always included if present)

    Returns:
        Sampled subgraph with query genes + top intermediates by connectivity
    """
    # Identify query genes present in graph
    genes_in_graph = set(g for g in query_genes if g in graph.nodes())

    if not genes_in_graph:
        logger.warning("No query genes found in graph")
        return graph

    logger.info(f"Found {len(genes_in_graph)} query genes in graph")

    # Identify intermediate nodes (exclude query genes and disease)
    intermediate_nodes = []
    for node in graph.nodes():
        if node not in genes_in_graph and node != disease_curie:
            intermediate_nodes.append(node)

    total_intermediates = len(intermediate_nodes)
    logger.info(f"Found {total_intermediates} intermediate nodes")

    # If we have fewer intermediates than the limit, return full graph
    if total_intermediates <= max_intermediates:
        logger.info(f"Total intermediates ({total_intermediates}) <= max ({max_intermediates}) - returning full graph")
        return graph

    logger.info(f"Sampling intermediates: {total_intermediates} → {max_intermediates}")

    # Score each intermediate by connections to query genes
    intermediate_scores = {}
    for node in intermediate_nodes:
        # Count edges to/from query genes
        score = 0
        for gene in genes_in_graph:
            if graph.has_edge(gene, node):
                score += 1
            if graph.has_edge(node, gene):
                score += 1

        intermediate_scores[node] = score

    # Sort intermediates by score (descending) and take top-N
    top_intermediates = sorted(
        intermediate_scores.items(),
        key=lambda x: x[1],
        reverse=True
    )[:max_intermediates]

    logger.info(
        f"Top intermediate scores - max: {top_intermediates[0][1] if top_intermediates else 0}, "
        f"min: {top_intermediates[-1][1] if top_intermediates else 0}"
    )

    # Build set of nodes to include
    nodes_to_include = set(genes_in_graph)  # Always include query genes
    nodes_to_include.update(node for node, _ in top_intermediates)  # Add top intermediates

    # Always include disease if present
    if disease_curie and disease_curie in graph.nodes():
        nodes_to_include.add(disease_curie)
        logger.info(f"Including disease node: {disease_curie}")

    # Create subgraph with selected nodes
    sampled_graph = graph.subgraph(nodes_to_include).copy()

    logger.info(
        f"Sampled graph: {sampled_graph.number_of_nodes()} nodes, "
        f"{sampled_graph.number_of_edges()} edges "
        f"(query genes: {len(genes_in_graph)}, intermediates: {len(top_intermediates)})"
    )

    return sampled_graph


def get_node_details(node_id: str, graph: nx.DiGraph) -> Dict[str, Any]:
    """Extract detailed information about a node.

    Args:
        node_id: Node identifier (CURIE)
        graph: NetworkX graph containing the node

    Returns:
        Dictionary with node details including edges and sources
    """
    if node_id not in graph.nodes():
        return {"error": f"Node {node_id} not found in graph"}

    attrs = graph.nodes[node_id]

    # Get edges
    in_edges = list(graph.in_edges(node_id, data=True))
    out_edges = list(graph.out_edges(node_id, data=True))

    # Group edges by predicate
    edges_by_predicate = {}
    all_sources = set()

    for src, tgt, data in in_edges + out_edges:
        predicate = data.get("predicate", "unknown").replace("biolink:", "")
        sources = data.get("knowledge_source", [])

        if predicate not in edges_by_predicate:
            edges_by_predicate[predicate] = []

        edge_info = {
            "direction": "incoming" if tgt == node_id else "outgoing",
            "source" if tgt == node_id else "target": src if tgt == node_id else tgt,
            "source_label" if tgt == node_id else "target_label": graph.nodes[src if tgt == node_id else tgt].get(
                "label", src if tgt == node_id else tgt
            ),
            "sources": sources,
        }
        edges_by_predicate[predicate].append(edge_info)

        # Collect all sources
        all_sources.update(sources)

    return {
        "node_id": node_id,
        "label": attrs.get("label", node_id),
        "original_symbol": attrs.get("original_symbol", ""),
        "category": attrs.get("category", "Unknown"),
        "is_query_gene": attrs.get("is_query_gene", False),
        "metrics": {
            "gene_frequency": attrs.get("gene_frequency", 0),
            "pagerank": attrs.get("pagerank", 0),
            "betweenness": attrs.get("betweenness", 0),
            "degree": graph.degree(node_id),
            "in_degree": graph.in_degree(node_id),
            "out_degree": graph.out_degree(node_id),
        },
        "edges_by_predicate": edges_by_predicate,
        "total_edges": len(in_edges) + len(out_edges),
        "knowledge_sources": list(all_sources),
    }


def create_clustered_graph(
    graph: nx.DiGraph, clustering_results, query_genes: List[str]
) -> nx.DiGraph:
    """Create a meta-graph where nodes represent clusters.

    Args:
        graph: Original NetworkX graph
        clustering_results: ClusteringResults from clustering_engine
        query_genes: List of query gene node IDs

    Returns:
        Meta-graph with cluster nodes
    """
    meta_graph = nx.DiGraph()

    # Create a node for each community
    for community in clustering_results.communities:
        cluster_id = f"cluster_{community.community_id}"

        # Get cluster statistics
        cluster_nodes = community.nodes
        subgraph = graph.subgraph(cluster_nodes)

        # Count query genes in cluster
        query_genes_in_cluster = [n for n in cluster_nodes if n in query_genes]

        # Get top nodes by PageRank
        top_nodes = community.top_nodes[:5] if community.top_nodes else []
        top_labels = [node_info.get("label", "") for node_info in top_nodes]

        # Add meta-node
        meta_graph.add_node(
            cluster_id,
            label=f"Cluster {community.community_id}\n({len(cluster_nodes)} nodes)",
            original_symbol="",  # No original symbol for clusters
            category="Cluster",  # Special category for clusters
            is_query_gene=False,
            size=len(cluster_nodes),
            density=community.density,
            top_nodes=top_labels,
            query_gene_count=len(query_genes_in_cluster),
            node_ids=cluster_nodes,  # Store original node IDs
            # Metrics for visualization
            gene_frequency=len(query_genes_in_cluster),
            pagerank=0,
            betweenness=0,
        )

    # Add edges between clusters
    community_map = {}  # node_id → community_id
    for community in clustering_results.communities:
        for node in community.nodes:
            community_map[node] = community.community_id

    # Count inter-cluster edges
    inter_cluster_edges = {}  # (comm1, comm2) → count

    for src, tgt in graph.edges():
        src_comm = community_map.get(src)
        tgt_comm = community_map.get(tgt)

        if src_comm is not None and tgt_comm is not None and src_comm != tgt_comm:
            edge_key = tuple(sorted([src_comm, tgt_comm]))
            inter_cluster_edges[edge_key] = inter_cluster_edges.get(edge_key, 0) + 1

    # Add inter-cluster edges to meta-graph
    for (comm1, comm2), count in inter_cluster_edges.items():
        meta_graph.add_edge(
            f"cluster_{comm1}",
            f"cluster_{comm2}",
            predicate=f"{count} edges",
            weight=count,
            label=f"{count} edges",
        )

    return meta_graph
