"""Tests for agent.tools: registry dispatch, schemas, and tool behaviors (fake gateway)."""

import types

import pandas as pd
import pytest

from geneset_translator.agent.tct_adapter import DictResultStash
from geneset_translator.agent import tools as T


class _FakeKG:
    def __init__(self, edges):
        self._edges = edges

    def __len__(self):
        return len(self._edges)

    def items(self):
        return list(self._edges.items())


class _FakeGateway:
    def resolve_genes(self, names, only_taxa=None, biolink_type=None):
        self.only_taxa = only_taxa
        self.biolink_type = biolink_type
        return {n: {"curie": f"NCBIGene:{i}", "types": ["biolink:Gene"]} for i, n in enumerate(names)}

    def neighborhood(self, gene_curie, target_categories=None):
        return types.SimpleNamespace(
            input_node_id=gene_curie,
            ranked=pd.DataFrame([{"output_node": "CHEBI:1", "Name": "drug"}]),
        )

    def path(self, n1, n2, intermediates):
        return types.SimpleNamespace(
            node1_id=n1, node2_id=n2, paths=pd.DataFrame([{"hop": 1}, {"hop": 2}])
        )

    def gene_network(self, curies):
        return _FakeKG(
            {"e0": {"subject": curies[0], "predicate": "biolink:interacts_with", "object": curies[1]}}
        )

    def edge_evidence(self, result, subject, object):
        return {"subject": subject, "object": object, "publications": ["PMID:1"], "edges_matched": 1}


def _ctx(query_genes=None):
    return T.ToolContext(
        tct=_FakeGateway(),
        stash=DictResultStash(),
        query_gene_curies=query_genes or [],
    )


def test_registry_schemas_cover_all_tools():
    reg = T.default_registry()
    names = {s["name"] for s in reg.schemas()}
    assert names == {
        "resolve_genes",
        "gene_neighborhood",
        "path_between",
        "gene_network",
        "edge_evidence",
        "cell_type_expression",
        "node_metadata",
        "data_sources",
        "filter_graph",
        "add_disease_node",
        "show_result",
        "cluster_graph",
    }
    for schema in reg.schemas():
        assert "description" in schema and "input_schema" in schema


def test_display_tools_are_not_parallel_safe_finders_are():
    reg = T.default_registry()
    assert reg.is_parallel_safe("gene_neighborhood") is True
    assert reg.is_parallel_safe("resolve_genes") is True
    # new read-only / network tools stay parallel-safe (and cacheable)
    assert reg.is_parallel_safe("node_metadata") is True
    assert reg.is_parallel_safe("data_sources") is True
    for name in ("filter_graph", "add_disease_node", "show_result", "cluster_graph"):
        assert reg.is_parallel_safe(name) is False


def test_cluster_graph_returns_digest_and_recolors_display():
    """cluster_graph summarizes the working graph and tags display nodes with their cluster so the
    viz can recolor. On a star (disease hub + gene leaves) it must report a hub-dominated topology
    with facets, not invent communities."""
    reg = T.default_registry()
    ctx = _ctx(query_genes=[])
    ctx.disease_curie = "MONDO:1"

    import networkx as nx

    g = nx.MultiDiGraph()
    g.add_node("MONDO:1", category="Disease", label="d")
    for i in range(30):
        n = f"NCBIGene:{i}"
        g.add_node(n, category="Gene", label=n)
        g.add_edge("MONDO:1", n, key=f"e{i}", predicate="biolink:associated_with")
    ctx.display.replace(g)

    result, is_error = reg.dispatch("cluster_graph", {}, ctx)
    assert is_error is False
    assert result["topology"] == "hub_dominated"
    assert "facet" in result["method"]
    assert result["category_facets"].get("Gene") == 30
    assert result["notes"]  # explains why community detection was skipped
    # display nodes are tagged with a cluster for recoloring
    snap = ctx.display.snapshot()
    assert any("cluster" in snap.nodes[n] for n in snap.nodes)


def test_cluster_graph_needs_a_graph():
    reg = T.default_registry()
    result, is_error = reg.dispatch("cluster_graph", {}, _ctx())
    assert is_error is True and "no graph" in result["error"]


# --------------------------------------------------------------------------------------
# node_metadata, data_sources, expression scope
# --------------------------------------------------------------------------------------
def test_node_metadata_compacts_go_and_gene_type(monkeypatch):
    class _FakeAnnotator:
        def __init__(self, *a, **k):
            pass

        def annotate_nodes(self, curies, fields="all"):
            return {c: types.SimpleNamespace(label=f"label-{c}") for c in curies if c != "NCBIGene:missing"}

        def _extract_annotation_features(self, node):
            return {
                "go_bp": [f"bp{i}" for i in range(15)],  # truncated to 10
                "go_mf": ["ATP binding"],
                "type_of_gene": "protein-coding",
                "alias": ["A", "B"],
            }

    monkeypatch.setattr("geneset_translator.core.node_annotator.NodeAnnotator", _FakeAnnotator)
    reg = T.default_registry()
    ctx = _ctx()
    result, is_error = reg.dispatch(
        "node_metadata", {"node_curies": ["NCBIGene:1", "NCBIGene:missing"]}, ctx
    )
    assert is_error is False
    assert result["n"] == 1  # the missing one is dropped
    meta = result["node_metadata"]["NCBIGene:1"]
    assert meta["label"] == "label-NCBIGene:1"
    assert meta["type_of_gene"] == "protein-coding"
    assert len(meta["go_bp"]) == 10  # truncated
    assert meta["aliases"] == ["A", "B"]


def test_data_sources_summarizes_infores(monkeypatch):
    edges = {
        "e0": {
            "subject": "NCBIGene:1",
            "object": "CHEBI:D",
            "predicate": "biolink:affects",
            "sources": [
                {"resource_id": "infores:drugbank", "resource_role": "primary_knowledge_source"},
                {"resource_id": "infores:biothings", "resource_role": "aggregator_knowledge_source"},
            ],
        }
    }
    fake_result = types.SimpleNamespace(
        knowledge_graph=types.SimpleNamespace(items=lambda: list(edges.items()))
    )
    reg = T.default_registry()
    ctx = _ctx()
    rid = ctx.stash.put(fake_result)

    monkeypatch.setattr(
        "geneset_translator.utils.infores_utils.download_infores_catalog",
        lambda cache_dir: {"raw": True},
    )
    monkeypatch.setattr(
        "geneset_translator.utils.infores_utils.parse_infores_catalog",
        lambda catalog: {
            "infores:drugbank": {"id": "infores:drugbank", "name": "DrugBank",
                                 "knowledge_level": "knowledge_assertion", "agent_type": "manual_agent"},
            "infores:biothings": {"id": "infores:biothings", "name": "BioThings",
                                  "knowledge_level": "knowledge_assertion", "agent_type": "automated_agent"},
        },
    )
    result, is_error = reg.dispatch("data_sources", {"result_id": rid}, ctx)
    assert is_error is False
    assert result["sources_in_result"] == 2
    assert result["total_sources"] == 2
    top_ids = {s["id"] for s in result["top_sources"]}
    assert "infores:drugbank" in top_ids


def test_data_sources_unknown_result_id_errors():
    reg = T.default_registry()
    result, is_error = reg.dispatch("data_sources", {"result_id": "nope"}, _ctx())
    assert is_error is True and "unknown result_id" in result["error"]


def test_compact_expression_scope_filters_hpa_keys():
    features = {
        "hpa_cell_type_specificity": "Cell type enhanced",
        "hpa_top_cell_types": [("T-cell", 9.0)],
        "hpa_tissue_specificity": "Tissue enhanced",
        "hpa_top_tissues": [("liver", 5.0)],
        "hpa_immune_cell_specificity": "Immune cell enhanced",
        "hpa_top_immune_cells": [("NK-cell", 3.0)],
        "go_bp": ["apoptosis"],
    }
    cell = T._compact_expression(features, "cell_type")
    assert set(cell) == {"hpa_cell_type_specificity", "hpa_top_cell_types"}  # no immune/tissue
    tissue = T._compact_expression(features, "tissue")
    assert set(tissue) == {"hpa_tissue_specificity", "hpa_top_tissues"}
    immune = T._compact_expression(features, "immune")
    assert set(immune) == {"hpa_immune_cell_specificity", "hpa_top_immune_cells"}
    allf = T._compact_expression(features, "all")
    assert "go_bp" in allf and "hpa_tissue_specificity" in allf and "hpa_immune_cell_specificity" in allf


def _seed_display(ctx):
    """Put a small normalized graph into ctx.display: drug D -> g1,g2; g3 separate."""
    import networkx as nx

    g = nx.MultiDiGraph()
    for n, cat in [("NCBIGene:1", "Gene"), ("NCBIGene:2", "Gene"), ("NCBIGene:3", "Gene"),
                   ("CHEBI:D", "ChemicalEntity")]:
        g.add_node(n, category=cat, label=n, curie=n, is_query_gene=cat == "Gene")
    g.add_edge("CHEBI:D", "NCBIGene:1", key="e0", predicate="biolink:affects")
    g.add_edge("CHEBI:D", "NCBIGene:2", key="e1", predicate="biolink:affects")
    ctx.display.replace(g)


def test_filter_graph_connected_to_mutates_display():
    reg = T.default_registry()
    ctx = _ctx(query_genes=["NCBIGene:1", "NCBIGene:2", "NCBIGene:3"])
    _seed_display(ctx)
    result, is_error = reg.dispatch("filter_graph", {"connected_to": "CHEBI:D"}, ctx)
    assert is_error is False
    assert set(ctx.display.snapshot().nodes) == {"CHEBI:D", "NCBIGene:1", "NCBIGene:2"}
    assert result["nodes_after"] == 3


def test_filter_graph_requires_a_criterion():
    reg = T.default_registry()
    ctx = _ctx(query_genes=["NCBIGene:1"])
    _seed_display(ctx)
    result, is_error = reg.dispatch("filter_graph", {}, ctx)
    assert is_error is True and "specify at least one" in result["error"]


def _seed_clustered_display(ctx):
    """Two clusters (C1, C2), each with a clear intra-cluster hub plus low-degree members."""
    import networkx as nx

    g = nx.MultiDiGraph()
    # cluster C1: hub NCBIGene:1 connects to 1a,1b,1c; 1d is a leaf
    # cluster C2: hub CHEBI:H connects to H1,H2; H3 is a leaf
    nodes = {
        "NCBIGene:1": ("Gene", "C1"), "1a": ("Gene", "C1"), "1b": ("Gene", "C1"),
        "1c": ("Gene", "C1"), "1d": ("Gene", "C1"),
        "CHEBI:H": ("ChemicalEntity", "C2"), "H1": ("ChemicalEntity", "C2"),
        "H2": ("ChemicalEntity", "C2"), "H3": ("ChemicalEntity", "C2"),
    }
    for n, (cat, cid) in nodes.items():
        g.add_node(n, category=cat, label=n, curie=n, cluster=cid, is_query_gene=n == "NCBIGene:1")
    for t in ("1a", "1b", "1c", "1d"):
        g.add_edge("NCBIGene:1", t, key=f"e_{t}", predicate="biolink:interacts_with")
    g.add_edge("1a", "1b", key="e_1a1b", predicate="biolink:interacts_with")  # 1a slightly higher degree
    for t in ("H1", "H2", "H3"):
        g.add_edge("CHEBI:H", t, key=f"e_{t}", predicate="biolink:affects")
    g.add_edge("H1", "H2", key="e_h1h2", predicate="biolink:affects")
    ctx.display.replace(g)


def test_filter_graph_top_per_cluster_is_balanced_across_clusters():
    reg = T.default_registry()
    ctx = _ctx(query_genes=["NCBIGene:1"])
    _seed_clustered_display(ctx)
    result, is_error = reg.dispatch("filter_graph", {"top_per_cluster": 2}, ctx)
    assert is_error is False
    kept = set(ctx.display.snapshot().nodes)
    c1 = {n for n in kept if n in {"NCBIGene:1", "1a", "1b", "1c", "1d"}}
    c2 = {n for n in kept if n in {"CHEBI:H", "H1", "H2", "H3"}}
    # both clusters represented (not gene-dominated), and the per-cluster hubs are kept
    assert "NCBIGene:1" in c1 and "CHEBI:H" in c2
    assert len(c1) >= 2 and len(c2) >= 2
    assert "top_2_per_cluster" in result["applied"]


def test_filter_graph_top_per_cluster_needs_clusters_first():
    reg = T.default_registry()
    ctx = _ctx(query_genes=["NCBIGene:1"])
    _seed_display(ctx)  # plain display, no 'cluster' attrs
    result, is_error = reg.dispatch("filter_graph", {"top_per_cluster": 3}, ctx)
    assert is_error is True and "cluster_graph first" in result["error"]


def test_show_result_sets_display_from_stash():
    reg = T.default_registry()
    ctx = _ctx(query_genes=["NCBIGene:1"])
    # gene_neighborhood stashes a fake nb whose to_networkx is exercised via tct_result_to_nx;
    # simpler: stash a fake KnowledgeGraph-like result with to_networkx.
    import networkx as nx
    import types

    tg = nx.MultiDiGraph()
    tg.add_node("NCBIGene:1", label="g1")
    tg.add_node("CHEBI:Z", label="z")
    tg.add_edge("NCBIGene:1", "CHEBI:Z", key="e0", predicate="biolink:affects",
                primary_sources=[], aggregator_sources=[], publications=[], supporting_text=[],
                confidence_scores={})
    fake_kg = types.SimpleNamespace(knowledge_graph=types.SimpleNamespace(
        to_networkx=lambda resolve_names=False, include_attributes=False: tg))
    rid = ctx.stash.put(fake_kg)
    result, is_error = reg.dispatch("show_result", {"result_id": rid}, ctx)
    assert is_error is False
    assert "CHEBI:Z" in ctx.display.snapshot().nodes


def test_add_disease_node_links_shown_genes():
    import types

    reg = T.default_registry()
    ctx = _ctx(query_genes=["NCBIGene:1", "NCBIGene:2", "NCBIGene:3"])
    ctx.disease_curie = "MONDO:0100096"
    ctx.disease_label = "COVID-19"
    _seed_display(ctx)  # NCBIGene:1,2,3 + CHEBI:D
    # disease->gene neighborhood returns two of the shown genes
    ctx.tct.neighborhood = lambda curie, cats=None: types.SimpleNamespace(
        ranked=pd.DataFrame([{"output_node": "NCBIGene:1"}, {"output_node": "NCBIGene:2"}])
    )
    result, is_error = reg.dispatch("add_disease_node", {}, ctx)
    assert is_error is False
    snap = ctx.display.snapshot()
    assert snap.nodes["MONDO:0100096"]["category"] == "Disease"
    assert snap.nodes["MONDO:0100096"]["label"] == "COVID-19"
    assert snap.has_edge("MONDO:0100096", "NCBIGene:1")
    assert result["linked_genes"] == 2


def test_finder_resets_display():
    reg = T.default_registry()
    ctx = _ctx(query_genes=["NCBIGene:1"])
    _seed_display(ctx)
    assert not ctx.display.is_empty()
    reg.dispatch("gene_neighborhood", {"gene_curie": "NCBIGene:1"}, ctx)
    assert ctx.display.is_empty()  # finder cleared the curated display


def test_dispatch_unknown_tool_is_error():
    reg = T.default_registry()
    result, is_error = reg.dispatch("does_not_exist", {}, _ctx())
    assert is_error is True
    assert "unknown tool" in result["error"]


def test_dispatch_tool_exception_is_error():
    reg = T.default_registry()
    # edge_evidence with a missing result_id raises ValueError -> surfaced as is_error
    result, is_error = reg.dispatch(
        "edge_evidence", {"result_id": "nope", "subject": "A", "object": "B"}, _ctx()
    )
    assert is_error is True
    assert "unknown result_id" in result["error"]


def test_resolve_genes_records_symbols_and_passes_taxon():
    reg = T.default_registry()
    ctx = _ctx()
    result, is_error = reg.dispatch("resolve_genes", {"names": ["FLT3", "BCL2"]}, ctx)
    assert is_error is False
    assert result["FLT3"]["curie"] == "NCBIGene:0"
    # curie -> symbol captured for friendly rendering
    assert ctx.curie_to_symbol["NCBIGene:0"] == "FLT3"
    assert ctx.tool_log == [{"name": "resolve_genes", "args": {"names": ["FLT3", "BCL2"]}}]


def test_gene_neighborhood_stashes_and_sets_latest_result_id():
    reg = T.default_registry()
    ctx = _ctx()
    result, is_error = reg.dispatch("gene_neighborhood", {"gene_curie": "NCBIGene:2322"}, ctx)
    assert is_error is False
    assert result["n_results"] == 1
    assert result["top_edges"] == [{"output_node": "CHEBI:1", "Name": "drug"}]
    # finder result stashed and marked as the latest to render
    assert ctx.latest_result_id == result["result_id"]
    assert ctx.stash.get(result["result_id"]) is not None


def test_gene_network_summarizes_edges():
    reg = T.default_registry()
    ctx = _ctx()
    result, is_error = reg.dispatch(
        "gene_network", {"gene_curies": ["NCBIGene:1", "NCBIGene:2"]}, ctx
    )
    assert is_error is False
    assert result["edge_count"] == 1
    assert result["top_edges"][0]["predicate"] == "biolink:interacts_with"


def test_dispatch_memoizes_identical_calls():
    reg = T.default_registry()
    ctx = _ctx()
    calls = {"n": 0}
    base = ctx.tct.resolve_genes

    def counting(names, only_taxa=None, biolink_type=None):
        calls["n"] += 1
        return base(names, only_taxa=only_taxa, biolink_type=biolink_type)

    ctx.tct.resolve_genes = counting
    r1, _ = reg.dispatch("resolve_genes", {"names": ["FLT3"]}, ctx)
    r2, _ = reg.dispatch("resolve_genes", {"names": ["FLT3"]}, ctx)
    assert calls["n"] == 1  # second call served from query_cache (no extra API call)
    assert r1 == r2
    reg.dispatch("resolve_genes", {"names": ["BCL2"]}, ctx)  # different args -> not cached
    assert calls["n"] == 2


def test_path_between_defaults_intermediates():
    reg = T.default_registry()
    ctx = _ctx()
    result, is_error = reg.dispatch(
        "path_between", {"node1_curie": "NCBIGene:596", "node2_curie": "CHEBI:1"}, ctx
    )
    assert is_error is False
    assert result["n_paths"] == 2
    assert result["node1_id"] == "NCBIGene:596"


# --------------------------------------------------------------------------------------
# Evidence enrichment (publication links + primary source on finder rows)
# --------------------------------------------------------------------------------------
def test_edge_publications_from_direct_key_and_attributes():
    # synthetic/recorded edge: publications live directly on the edge
    assert T._edge_publications({"publications": ["PMID:28644114", "28644114", "PMID:28644114"]}) == [
        "PMID:28644114"
    ]  # normalized + de-duplicated
    assert T._edge_publications({}) == []


def test_edge_primary_source_from_sources_and_synthetic():
    trapi = {"sources": [
        {"resource_id": "infores:aggA", "resource_role": "aggregator_knowledge_source"},
        {"resource_id": "infores:drugbank", "resource_role": "primary_knowledge_source"},
    ]}
    assert T._edge_primary_source(trapi) == "infores:drugbank"
    assert T._edge_primary_source({"primary_sources": ["infores:chembl"]}) == "infores:chembl"
    assert T._edge_primary_source({}) is None


def test_gene_neighborhood_attaches_publications_to_top_edges():
    """A neighborhood result with a knowledge graph carrying publications -> top_edges carry citable
    publication links + primary source for the input->output edge."""
    import types

    edge = {
        "subject": "NCBIGene:2322",
        "object": "CHEBI:1",
        "predicate": "biolink:affects",
        "publications": ["PMID:28644114"],
        "primary_sources": ["infores:drugbank"],
    }
    nb = types.SimpleNamespace(
        input_node_id="NCBIGene:2322",
        ranked=pd.DataFrame([{"output_node": "CHEBI:1", "Name": "drug"}]),
        knowledge_graph=_FakeKG({"e0": edge}),
    )

    class _GW(_FakeGateway):
        def neighborhood(self, gene_curie, target_categories=None):
            return nb

    reg = T.default_registry()
    ctx = T.ToolContext(tct=_GW(), stash=DictResultStash())
    result, is_error = reg.dispatch("gene_neighborhood", {"gene_curie": "NCBIGene:2322"}, ctx)
    assert is_error is False
    row = result["top_edges"][0]
    assert row["publications"] == ["PMID:28644114"]
    assert row["primary_source"] == "infores:drugbank"


def test_gene_network_attaches_publications_to_top_edges():
    edge = {
        "subject": "NCBIGene:1",
        "object": "NCBIGene:2",
        "predicate": "biolink:interacts_with",
        "publications": ["PMID:20705237"],
        "sources": [{"resource_id": "infores:biogrid", "resource_role": "primary_knowledge_source"}],
    }

    class _GW(_FakeGateway):
        def gene_network(self, curies):
            return _FakeKG({"e0": edge})

    reg = T.default_registry()
    ctx = T.ToolContext(tct=_GW(), stash=DictResultStash())
    result, is_error = reg.dispatch(
        "gene_network", {"gene_curies": ["NCBIGene:1", "NCBIGene:2"]}, ctx
    )
    assert is_error is False
    row = result["top_edges"][0]
    assert row["publications"] == ["PMID:20705237"]
    assert row["primary_source"] == "infores:biogrid"
