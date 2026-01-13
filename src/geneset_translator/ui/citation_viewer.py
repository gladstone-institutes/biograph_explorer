"""Citation graph viewer for visualizing cited nodes and edges.

This module provides interactive visualization of citation subgraphs,
highlighting nodes and edges that support specific LLM-generated claims.
"""

import logging
from typing import List, Optional

import networkx as nx
import streamlit as st
from streamlit_cytoscape import streamlit_cytoscape

from geneset_translator.core.llm_summarizer import CitationGraph
from geneset_translator.ui.network_viz import render_network_visualization

logger = logging.getLogger(__name__)


def render_citation_viewer(
    citation: CitationGraph,
    graph: nx.MultiDiGraph,
    query_genes: Optional[List[str]] = None,
    category: Optional[str] = None
):
    """Render interactive citation graph with highlighted cited nodes/edges.

    Args:
        citation: Citation data with node/edge IDs
        graph: Full knowledge graph
        query_genes: Query gene CURIEs for highlighting
        category: Category name (unused, kept for API compatibility)
    """
    st.caption(f"**Claim:** {citation.claim}")

    # Extract subgraph with citation markings
    try:
        subgraph = extract_citation_subgraph(citation, graph)

        if subgraph.number_of_nodes() == 0:
            st.warning("No nodes found for this citation")
            return

        # Layout options
        layout_col1, layout_col2 = st.columns([1, 3])
        with layout_col1:
            layout_name = st.selectbox(
                "Layout",
                options=['dagre', 'cose', 'fcose', 'circle', 'grid'],
                index=0,
                key=f"citation_layout_{citation.citation_id}"
            )

        # Prepare visualization using the same method as main network view
        viz_data = render_network_visualization(
            subgraph,
            query_genes=[],  # Don't highlight query genes in citation view
            sizing_metric="gene_frequency",
            layout=layout_name,
            max_intermediates=1000,  # Include all nodes (citation subgraphs are small)
            use_metric_sizing=True,  # Use metric sizing like main view
            base_node_size=30,  # Match main view
            edge_width=2,  # Match main view
        )

        if not viz_data:
            st.error("Failed to prepare citation graph visualization")
            return

        # Render with same approach as main network view (no custom styling)
        streamlit_cytoscape(
            viz_data["elements"],
            layout=viz_data["layout"],
            node_styles=viz_data["node_styles"],
            edge_styles=viz_data["edge_styles"],
            key=f"citation_graph_{citation.citation_id}",
            hide_underscore_attrs=True,
        )

        st.caption("""
        **:material/lightbulb: How to explore:**
        - **Drag** to pan • **Scroll** to zoom • **Click** node or edge to view information
        - **Fullscreen** in top-right
        """)

        # Stats
        col1, col2 = st.columns(2)
        with col1:
            st.metric("Nodes", subgraph.number_of_nodes(), border=True)
        with col2:
            st.metric("Edges", subgraph.number_of_edges(), border=True)

    except Exception as e:
        logger.error(f"Failed to render citation viewer: {e}")
        st.error(f"Failed to render citation graph: {e}")


def extract_citation_subgraph(citation: CitationGraph, graph: nx.MultiDiGraph) -> nx.MultiDiGraph:
    """Extract subgraph containing ONLY cited nodes and cited edges.

    Args:
        citation: Citation data
        graph: Full knowledge graph

    Returns:
        Subgraph with only cited nodes and cited edges, with _cited attributes marked
    """
    subgraph_nodes = set()

    # Add cited nodes
    for node_id in citation.node_ids:
        if node_id in graph.nodes():
            subgraph_nodes.add(node_id)

    # Parse cited edge IDs and add nodes connected by cited edges
    # Edge ID format: subject→predicate→object (e.g., "NCBIGene:6776→has_part→UMLS:C1999216")
    cited_edge_set = set()
    for edge_id in citation.edge_ids:
        try:
            parts = edge_id.split('→')
            if len(parts) == 3:
                subject, predicate, obj = parts
                # Find matching edges in graph by subject, object, and predicate prefix
                for u, v, key, data in graph.edges(keys=True, data=True):
                    if u == subject and v == obj:
                        # Check if predicate matches (key format is "predicate_N")
                        predicate_from_key = key.rsplit('_', 1)[0] if '_' in key else key
                        if predicate_from_key == predicate:
                            cited_edge_set.add((u, v, key))
                            break
                # Add nodes from cited edges
                if subject in graph.nodes():
                    subgraph_nodes.add(subject)
                if obj in graph.nodes():
                    subgraph_nodes.add(obj)
        except Exception as e:
            logger.warning(f"Failed to parse cited edge {edge_id}: {e}")

    # Extract subgraph
    subgraph = nx.MultiDiGraph()

    # Add nodes with their attributes
    for node_id in subgraph_nodes:
        if node_id in graph.nodes():
            subgraph.add_node(node_id, **graph.nodes[node_id])
            # Mark as cited if in citation node list
            subgraph.nodes[node_id]['_cited'] = node_id in citation.node_ids

    # Add ONLY cited edges
    for u, v, key in cited_edge_set:
        if graph.has_edge(u, v, key):
            subgraph.add_edge(u, v, key=key, **graph.edges[u, v, key])
            subgraph.edges[u, v, key]['_cited'] = True

    return subgraph
