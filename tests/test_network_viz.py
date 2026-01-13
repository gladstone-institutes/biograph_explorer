"""Tests for network visualization stability and collapse state preservation."""
import pytest
import networkx as nx
from geneset_translator.ui.network_viz import (
    prepare_cytoscape_elements,
    create_edge_styles,
)


def test_edge_elements_stable_across_width_changes():
    """Verify edge elements don't change when edge_width changes.

    This test ensures that edge elements remain stable when the edge_width
    slider changes, which is critical for preserving collapse state.
    """
    # Create a simple MultiDiGraph with parallel edges
    G = nx.MultiDiGraph()
    G.add_node("NCBIGene:348", category="Gene", label="APOE", gene_frequency=5)
    G.add_node("NCBIGene:351", category="Gene", label="APP", gene_frequency=4)

    # Add parallel edges (multiple predicates between same nodes)
    # Sources should be dicts with resource_id, as expected by the real code
    G.add_edge("NCBIGene:348", "NCBIGene:351", key=0,
               predicate="biolink:interacts_with",
               sources=[{"resource_id": "infores:source1", "resource_role": "primary_knowledge_source"}])
    G.add_edge("NCBIGene:348", "NCBIGene:351", key=1,
               predicate="biolink:regulates",
               sources=[{"resource_id": "infores:source2", "resource_role": "primary_knowledge_source"}])

    query_genes = ["NCBIGene:348", "NCBIGene:351"]

    # Generate elements (edge_width parameter removed)
    elements = prepare_cytoscape_elements(
        G, query_genes, sizing_metric="gene_frequency"
    )

    # Verify _edge_width not in edge data (this was causing instability)
    for edge in elements["edges"]:
        assert "_edge_width" not in edge["data"], \
            "Edge width should not be in element data (causes collapse instability)"
        assert "_edge_font_size" not in edge["data"], \
            "Edge font size should not be in element data (causes collapse instability)"

    # Verify color attributes ARE present (needed for publication highlighting)
    for edge in elements["edges"]:
        assert "_line_color" in edge["data"], \
            "Line color should be in element data for publication highlighting"
        assert "_target_arrow_color" in edge["data"], \
            "Arrow color should be in element data for publication highlighting"

    # Verify basic edge structure
    assert len(elements["edges"]) == 2, "Should have 2 parallel edges"
    assert elements["edges"][0]["data"]["source"] == "NCBIGene:348"
    assert elements["edges"][0]["data"]["target"] == "NCBIGene:351"


def test_edge_styles_apply_width_dynamically():
    """Verify EdgeStyles use static width values.

    This test ensures that edge width is applied via EdgeStyle custom_styles
    as a static value, not read from element data attributes.
    """
    G = nx.MultiDiGraph()
    G.add_edge("A", "B", key=0, predicate="biolink:interacts_with")
    G.add_edge("A", "B", key=1, predicate="biolink:regulates")

    # Create styles with different widths
    styles_w2 = create_edge_styles(G, edge_width=2)
    styles_w5 = create_edge_styles(G, edge_width=5)

    # Should have 2 styles (one per predicate)
    assert len(styles_w2) == 2, "Should have 2 EdgeStyle objects (one per predicate)"
    assert len(styles_w5) == 2

    # Extract width from custom_styles
    # All predicates should have same width (applied uniformly)
    for style in styles_w2:
        assert style.custom_styles["width"] == 2, \
            "Width should be static value 2, not data() reference"

    for style in styles_w5:
        assert style.custom_styles["width"] == 5, \
            "Width should be static value 5, not data() reference"

    # Font size should scale with width
    # Font size formula: max(8, min(14, 6 + edge_width))
    for style in styles_w2:
        expected_font = max(8, min(14, 6 + 2))  # = 8
        assert style.custom_styles["font-size"] == expected_font

    for style in styles_w5:
        expected_font = max(8, min(14, 6 + 5))  # = 11
        assert style.custom_styles["font-size"] == expected_font

    # Verify color attributes use data() reference (dynamic per-edge)
    for style in styles_w2:
        assert style.custom_styles["line-color"] == "data(_line_color)", \
            "Line color should use data() for per-edge publication highlighting"
        assert style.custom_styles["target-arrow-color"] == "data(_target_arrow_color)", \
            "Arrow color should use data() for per-edge publication highlighting"


def test_publication_highlighting_preserved():
    """Verify _line_color still in edge data for pub filter.

    This test ensures that publication highlighting continues to work
    after removing edge_width from element data.
    """
    G = nx.MultiDiGraph()
    G.add_node("NCBIGene:348", category="Gene", label="APOE")
    G.add_node("NCBIGene:351", category="Gene", label="APP")

    # Add edge with publication highlight flag
    G.add_edge("NCBIGene:348", "NCBIGene:351", key=0,
               predicate="biolink:interacts_with",
               _has_filtered_pub=True)  # Marked as having filtered publication

    # Add edge without highlight
    G.add_edge("NCBIGene:348", "NCBIGene:351", key=1,
               predicate="biolink:regulates",
               _has_filtered_pub=False)

    query_genes = ["NCBIGene:348", "NCBIGene:351"]
    elements = prepare_cytoscape_elements(G, query_genes)

    # Verify color attributes are present for both edges
    highlighted_edge = elements["edges"][0]
    normal_edge = elements["edges"][1]

    assert "_line_color" in highlighted_edge["data"]
    assert "_target_arrow_color" in highlighted_edge["data"]

    # Highlighted edge should have teal color (#17BECF)
    assert highlighted_edge["data"]["_line_color"] == "#17BECF", \
        "Highlighted edge should have teal color"

    # Normal edge should have default gray color (#999999)
    assert normal_edge["data"]["_line_color"] == "#999999", \
        "Normal edge should have default gray color"

    # Verify flag is preserved
    assert highlighted_edge["data"]["_has_filtered_pub"] is True
    assert normal_edge["data"]["_has_filtered_pub"] is False


def test_node_sizing_stable_across_base_size_changes():
    """Verify node elements remain stable when base_node_size changes.

    This test ensures that node size factor (not absolute size) is stored
    in element data, keeping elements stable for collapse state preservation.
    """
    G = nx.MultiDiGraph()
    G.add_node("NCBIGene:348", category="Gene", label="APOE",
               gene_frequency=10, pagerank=0.15, betweenness=0.05, degree=5)
    G.add_node("NCBIGene:351", category="Gene", label="APP",
               gene_frequency=5, pagerank=0.10, betweenness=0.02, degree=3)

    query_genes = ["NCBIGene:348"]

    # Test with metric-based sizing
    elements_metric = prepare_cytoscape_elements(
        G, query_genes, sizing_metric="gene_frequency", use_metric_sizing=True
    )

    # Test with uniform sizing
    elements_uniform = prepare_cytoscape_elements(
        G, query_genes, sizing_metric="gene_frequency", use_metric_sizing=False
    )

    # Verify _size_factor is present (not _size)
    for node in elements_metric["nodes"]:
        assert "_size_factor" in node["data"], "Node should have _size_factor attribute"
        assert "_size" not in node["data"], "Node should NOT have _size attribute (causes instability)"
        assert "_font_size" not in node["data"], "Node should NOT have _font_size attribute (causes instability)"

    # Metric-based: nodes should have different size factors
    node1_metric_factor = elements_metric["nodes"][0]["data"]["_size_factor"]
    node2_metric_factor = elements_metric["nodes"][1]["data"]["_size_factor"]
    assert node1_metric_factor != node2_metric_factor, \
        "With metric sizing, nodes should have different size factors based on gene_frequency"

    # Uniform: nodes should have same size factor (1.0)
    node1_uniform_factor = elements_uniform["nodes"][0]["data"]["_size_factor"]
    node2_uniform_factor = elements_uniform["nodes"][1]["data"]["_size_factor"]
    # Both should be 1.0 (or close, accounting for query gene boost)
    assert 0.9 <= node1_uniform_factor <= 1.1, "Uniform sizing should have factor ~1.0"
    assert 0.9 <= node2_uniform_factor <= 1.1, "Uniform sizing should have factor ~1.0"


def test_node_styles_use_mapdata_for_sizing():
    """Verify NodeStyles use mapData to calculate size from _size_factor.

    This test ensures that node size is calculated dynamically via mapData,
    not read from absolute size values in element data.
    """
    from geneset_translator.ui.network_viz import create_node_styles

    G = nx.MultiDiGraph()
    G.add_node("A", category="Gene")
    G.add_node("B", category="Protein")

    # Create styles with different base_node_size values
    styles_size30 = create_node_styles(G, base_node_size=30, use_metric_sizing=True)
    styles_size50 = create_node_styles(G, base_node_size=50, use_metric_sizing=True)

    # Both should have 2 styles (one per category)
    assert len(styles_size30) == 2
    assert len(styles_size50) == 2

    # Verify mapData is used for width/height (not data() reference)
    for style in styles_size30:
        assert "mapData(_size_factor" in style.custom_styles["width"], \
            "Width should use mapData(_size_factor, ...) for dynamic sizing"
        assert "mapData(_size_factor" in style.custom_styles["height"], \
            "Height should use mapData(_size_factor, ...) for dynamic sizing"
        assert "mapData(_size_factor" in style.custom_styles["font-size"], \
            "Font-size should use mapData(_size_factor, ...) for dynamic sizing"

    # Verify different base_node_size produces different mapData ranges
    # Extract the mapData parameters (rough check)
    style30_width = styles_size30[0].custom_styles["width"]
    style50_width = styles_size50[0].custom_styles["width"]
    assert style30_width != style50_width, \
        "Different base_node_size should produce different mapData ranges"
