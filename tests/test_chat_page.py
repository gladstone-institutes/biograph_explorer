"""Guards for the chat page: pandas import + example datasets load cleanly.

Regression: the example loader used pd.read_csv but chat_page didn't import pandas; the broad
try/except hid it as an st.error, so it must be covered by an explicit test.
"""

import pandas as pd

from geneset_translator.ui import chat_page


def test_chat_page_imports_pandas():
    assert hasattr(chat_page, "pd"), "chat_page must import pandas (used by the example/CSV loaders)"


def test_example_datasets_load_and_have_gene_symbol_column():
    assert chat_page.EXAMPLES, "expected at least one example dataset"
    for _label, (path, disease, disease_name) in chat_page.EXAMPLES.items():
        df = pd.read_csv(path)
        assert "gene_symbol" in df.columns
        symbols = df["gene_symbol"].dropna().astype(str).str.strip().tolist()
        assert len(symbols) > 0
        assert disease.startswith("MONDO:")
        assert disease_name  # a human-readable name for the system prompt


def test_looks_like_curie():
    assert chat_page._looks_like_curie("MONDO:0005361")
    assert chat_page._looks_like_curie("NCBIGene:2322")
    assert not chat_page._looks_like_curie("acute myeloid leukemia")
    assert not chat_page._looks_like_curie("")


def test_graph_component_key_changes_on_display_edit_and_collapse():
    """The cytoscape key must change when the agent edits the display (version bump) or the collapse
    toggle flips, so the component remounts and visibly refreshes; otherwise it stays stable."""
    base = chat_page._graph_component_key("res_1", "dagre", 200, 0, True)
    assert base == chat_page._graph_component_key("res_1", "dagre", 200, 0, True)  # stable
    assert base != chat_page._graph_component_key("res_1", "dagre", 200, 1, True)  # version bump
    assert base != chat_page._graph_component_key("res_1", "dagre", 200, 0, False)  # collapse toggled
    assert base != chat_page._graph_component_key("res_2", "dagre", 200, 0, True)  # new result
    assert base != chat_page._graph_component_key("res_1", "cose", 200, 0, True)  # layout change
    # Metric-sizing mode changes node element data -> must remount.
    assert base != chat_page._graph_component_key("res_1", "dagre", 200, 0, True, True, "degree")
    assert base != chat_page._graph_component_key("res_1", "dagre", 200, 0, True, True, "gene_frequency")


def test_meta_edge_style_carries_width_so_slider_reaches_collapsed_edges():
    """Collapsed meta-edges are styled separately from EdgeStyles; the meta-edge style must carry the
    chosen edge width (and a matching font size) or the edge-width control won't affect them."""
    style = chat_page._meta_edge_style(7)
    assert style["width"] == 7
    assert isinstance(style["font-size"], int)  # scales with width via calculate_edge_font_size
    assert chat_page._meta_edge_style(2)["width"] != chat_page._meta_edge_style(9)["width"]


def test_graph_component_key_stable_across_node_and_edge_size():
    """Node size / edge width are styling-only: they must NOT change the key (positions preserved;
    the new styles are applied in place via the recomputed viz payload)."""
    # the key builder doesn't take size args at all -> identical key regardless of slider values
    k = chat_page._graph_component_key("res_1", "dagre", 200, 3, True, False, "degree")
    assert k == chat_page._graph_component_key("res_1", "dagre", 200, 3, True, False, "degree")


# Emoji blocks the project bans (kept narrow so scientific letters like alpha/beta are NOT flagged).
_EMOJI_RANGES = [
    (0x1F300, 0x1FAFF), (0x1F000, 0x1F0FF), (0x1F1E6, 0x1F1FF),
    (0x2600, 0x27BF), (0x2B00, 0x2BFF), (0xFE00, 0xFE0F),
]


def _has_emoji(text: str) -> bool:
    return any(any(lo <= ord(ch) <= hi for lo, hi in _EMOJI_RANGES) for ch in text)


def test_starter_questions_with_disease_and_geneset():
    symbols = ["IL13", "CCL26", "POSTN"]
    qs = chat_page._starter_questions(symbols, "Eosinophilic Esophagitis", "MONDO:0005361")
    # Headline is the translational disease + druggability question.
    assert "Eosinophilic Esophagitis" in qs[0]
    joined = " ".join(qs).lower()
    assert "drug" in qs[0].lower() and ("target" in qs[0].lower() or "associated" in qs[0].lower())
    # Breadth: network, pathfinder (both genes), cell-type/tissue, provenance.
    assert any("interact with each other" in q for q in qs)  # gene_network
    path = next(q for q in qs if "connect" in q)
    assert "IL13" in path and "CCL26" in path  # path_between names the first two genes
    assert "cell types or tissues" in joined  # cell_type_expression
    assert "knowledge sources" in joined  # data_sources / provenance
    assert not any(_has_emoji(q) for q in qs)  # project-wide emoji ban


def test_starter_questions_without_disease():
    symbols = ["FLT3", "BCL2", "NPM1"]
    qs = chat_page._starter_questions(symbols, None, "")
    # No disease token anywhere; item 1 is the druggability + evidence question.
    assert all("MONDO" not in q for q in qs)
    assert "drug" in qs[0].lower() and "evidence" in qs[0].lower()
    assert "associated with" not in qs[0].lower()  # no dangling disease anchor
    # Set-based features still present.
    assert any("interact with each other" in q for q in qs)
    assert any("connect" in q and "FLT3" in q and "BCL2" in q for q in qs)
    assert any("cell types or tissues" in q for q in qs)


def test_starter_questions_single_gene_omits_path_question():
    qs = chat_page._starter_questions(["FLT3"], "acute myeloid leukemia", "MONDO:0018874")
    # With <2 genes the path_between chip is omitted (no malformed "connect X and ").
    assert not any("connect" in q for q in qs)
    assert "acute myeloid leukemia" in qs[0]
    assert any("interact with each other" in q for q in qs)  # other set questions remain


def test_graph_to_cyjs_is_full_and_preserves_nested_edge_attrs():
    """The full-graph export must be valid Cytoscape.js JSON, report the true (uncapped) node/edge
    counts, and keep nested edge attributes (e.g. a publications list) intact for Cytoscape import."""
    import json

    import networkx as nx

    g = nx.MultiDiGraph()
    g.add_node("NCBIGene:2322", name="FLT3", category="biolink:Gene")
    g.add_node("CHEBI:1", name="midostaurin", category="biolink:Drug")
    g.add_edge(
        "NCBIGene:2322", "CHEBI:1",
        predicate="biolink:affected_by",
        publications=["PMID:1", "PMID:2"],
        primary_source="infores:drugbank",
    )
    cyjs, n_nodes, n_edges = chat_page._graph_to_cyjs(g)
    assert (n_nodes, n_edges) == (2, 1)

    data = json.loads(cyjs)
    assert "elements" in data and "nodes" in data["elements"] and "edges" in data["elements"]
    assert len(data["elements"]["nodes"]) == 2
    edge_data = data["elements"]["edges"][0]["data"]
    assert edge_data["publications"] == ["PMID:1", "PMID:2"]  # nested list survived
    assert edge_data["primary_source"] == "infores:drugbank"
