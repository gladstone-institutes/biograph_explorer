"""Tests for agent.display_ops: DisplayGraph holder + pure graph transforms."""

import networkx as nx

from geneset_translator.agent import display_ops as do


def _graph():
    """drug D connects to genes g1,g2; g3 is isolated-ish; query genes g1,g2,g3."""
    g = nx.MultiDiGraph()
    for n, cat in [
        ("NCBIGene:1", "Gene"),
        ("NCBIGene:2", "Gene"),
        ("NCBIGene:3", "Gene"),
        ("CHEBI:D", "ChemicalEntity"),
        ("CHEBI:E", "ChemicalEntity"),
    ]:
        g.add_node(n, category=cat, label=n, curie=n, is_query_gene=cat == "Gene")
    g.add_edge("CHEBI:D", "NCBIGene:1", key="e0", predicate="biolink:affects")
    g.add_edge("CHEBI:D", "NCBIGene:2", key="e1", predicate="biolink:affects")
    g.add_edge("CHEBI:E", "NCBIGene:3", key="e2", predicate="biolink:affects")
    return g


# -- DisplayGraph holder ---------------------------------------------------------------
def test_display_graph_lifecycle():
    dg = do.DisplayGraph()
    assert dg.is_empty() and dg.version == 0
    dg.replace(_graph())
    assert not dg.is_empty() and dg.version == 1
    snap = dg.snapshot()
    snap.add_node("X")  # snapshot is isolated
    assert "X" not in dg.snapshot()
    dg.apply(do.filter_connected_to, "CHEBI:D")
    assert dg.version == 2  # replace (1) + apply (2); snapshot does not bump
    dg.reset()
    assert dg.is_empty()


# -- pure transforms -------------------------------------------------------------------
def test_filter_connected_to():
    out = do.filter_connected_to(_graph(), "CHEBI:D")
    assert set(out.nodes) == {"CHEBI:D", "NCBIGene:1", "NCBIGene:2"}
    assert "NCBIGene:3" not in out  # connected to E, not D


def test_trim_to_top_degree_keeps_always_keep():
    g = _graph()
    # CHEBI:D has degree 2 (top). Keep top_n=1 but always keep the 3 query genes.
    out = do.trim_to_top_degree(g, top_n=1, always_keep=["NCBIGene:1", "NCBIGene:2", "NCBIGene:3"])
    assert "CHEBI:D" in out
    assert {"NCBIGene:1", "NCBIGene:2", "NCBIGene:3"} <= set(out.nodes)
    assert "CHEBI:E" not in out  # not top-1, not in always_keep


def test_filter_by_category():
    out = do.filter_by_category(_graph(), "ChemicalEntity", always_keep=["NCBIGene:1"])
    assert {"CHEBI:D", "CHEBI:E", "NCBIGene:1"} == set(out.nodes)


def test_merge_graphs_unions_nodes_and_edges():
    g1 = do.filter_connected_to(_graph(), "CHEBI:D")  # D,g1,g2
    g2 = do.filter_connected_to(_graph(), "CHEBI:E")  # E,g3
    out = do.merge_graphs(g1, g2)
    assert {"CHEBI:D", "CHEBI:E", "NCBIGene:1", "NCBIGene:2", "NCBIGene:3"} == set(out.nodes)
    assert out.number_of_edges() == 3


def test_add_disease_node_links_present_genes_only():
    out = do.add_disease_node(
        _graph(), "MONDO:0100096", "COVID-19",
        associated_gene_curies=["NCBIGene:1", "NCBIGene:999"],  # 999 not present -> skipped
    )
    assert out.nodes["MONDO:0100096"]["category"] == "Disease"
    assert out.nodes["MONDO:0100096"]["label"] == "COVID-19"
    assert out.has_edge("MONDO:0100096", "NCBIGene:1")
    assert not out.has_edge("MONDO:0100096", "NCBIGene:999")


def test_finalize_sets_gene_frequency_and_renders():
    out = do.finalize(_graph(), ["NCBIGene:1", "NCBIGene:2", "NCBIGene:3"])
    # CHEBI:D touches two query genes -> frequency 2
    assert out.nodes["CHEBI:D"]["gene_frequency"] == 2
    assert out.nodes["CHEBI:E"]["gene_frequency"] == 1
    # the edited graph still renders through the real renderer
    from geneset_translator.ui.network_viz import prepare_cytoscape_elements

    elements = prepare_cytoscape_elements(out, ["NCBIGene:1", "NCBIGene:2", "NCBIGene:3"])
    assert len(elements["nodes"]) == 5
