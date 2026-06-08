"""Tests for agent.tct_adapter: stash, gateway (mocked TCT), and result->nx conversion."""

import networkx as nx
import pytest

from geneset_translator.agent import tct_adapter as ta


# --------------------------------------------------------------------------------------
# ResultStash
# --------------------------------------------------------------------------------------
def test_dict_result_stash_put_get_roundtrip():
    store: dict = {}
    stash = ta.DictResultStash(store)
    rid1 = stash.put({"a": 1})
    rid2 = stash.put({"b": 2})
    assert rid1 != rid2
    assert stash.get(rid1) == {"a": 1}
    assert stash.get(rid2) == {"b": 2}
    # backing dict persists (chat page keeps this in session_state)
    assert set(store) == {rid1, rid2}
    assert stash.get("missing") is None


# --------------------------------------------------------------------------------------
# resolve_genes (human taxon by default)
# --------------------------------------------------------------------------------------
class _FakeNode:
    def __init__(self, curie, categories, label):
        self.curie = curie
        self.categories = categories
        self.label = label


def test_resolve_genes_passes_human_taxon_and_shapes_output(monkeypatch):
    calls = []

    def fake_batch_lookup(strings, **kwargs):
        calls.append({"strings": list(strings), "kwargs": dict(kwargs)})
        return {
            "FLT3": _FakeNode("NCBIGene:2322", ["biolink:Gene", "biolink:Protein"], "FLT3"),
            "NOPE": None,
        }

    monkeypatch.setattr(ta.name_resolver, "batch_lookup", fake_batch_lookup)
    gw = ta.TctGateway(resources=object())
    out = gw.resolve_genes(["FLT3", "NOPE"])

    # pass 1 carries the human taxon filter
    assert calls[0]["kwargs"].get("only_taxa") == ta.HUMAN_TAXON
    assert out["FLT3"]["curie"] == "NCBIGene:2322"
    assert out["FLT3"]["types"][:2] == ["biolink:Gene", "biolink:Protein"]
    assert out["NOPE"] is None


def test_resolve_genes_retries_unresolved_without_taxon(monkeypatch):
    """A drug (no taxon) fails the human-taxon pass, then resolves on the taxon-free retry."""
    calls = []

    def fake_batch_lookup(strings, **kwargs):
        calls.append({"strings": list(strings), "kwargs": dict(kwargs)})
        if kwargs.get("only_taxa"):  # pass 1: gene resolves, drug does not
            return {
                "BCL2": _FakeNode("NCBIGene:596", ["biolink:Gene"], "BCL2"),
                "venetoclax": None,
            }
        return {"venetoclax": _FakeNode("PUBCHEM.COMPOUND:49846579", ["biolink:SmallMolecule"], "venetoclax")}

    monkeypatch.setattr(ta.name_resolver, "batch_lookup", fake_batch_lookup)
    gw = ta.TctGateway(resources=object())
    out = gw.resolve_genes(["BCL2", "venetoclax"])

    assert len(calls) == 2
    assert calls[0]["kwargs"].get("only_taxa") == ta.HUMAN_TAXON
    assert calls[1]["strings"] == ["venetoclax"]  # only the unresolved name is retried
    assert calls[1]["kwargs"].get("only_taxa") is None  # taxon filter dropped on retry
    assert out["BCL2"]["curie"] == "NCBIGene:596"
    assert out["venetoclax"]["curie"] == "PUBCHEM.COMPOUND:49846579"


def test_resolve_disease_drops_taxon_and_passes_type(monkeypatch):
    """A disease name resolved with biolink_type skips the taxon filter (single typed lookup),
    so it is not mis-mapped to a same-named gene."""
    calls = []

    def fake_batch_lookup(strings, **kwargs):
        calls.append({"strings": list(strings), "kwargs": dict(kwargs)})
        return {"acute myeloid leukemia": _FakeNode("MONDO:0018874", ["biolink:Disease"], "AML")}

    monkeypatch.setattr(ta.name_resolver, "batch_lookup", fake_batch_lookup)
    gw = ta.TctGateway(resources=object())
    out = gw.resolve_genes(["acute myeloid leukemia"], biolink_type="biolink:Disease")

    assert len(calls) == 1  # no human-taxon two-pass for a non-gene type
    assert calls[0]["kwargs"].get("only_taxa") is None
    assert calls[0]["kwargs"].get("biolink_type") == "biolink:Disease"
    assert out["acute myeloid leukemia"]["curie"] == "MONDO:0018874"


def test_neighborhood_unnormalizable_input_raises_clean_error(monkeypatch):
    """TCT.Neighborhood_finder raises AttributeError on an unnormalizable node; the gateway
    converts it to an actionable ValueError."""
    class _BoomTCT:
        @staticmethod
        def Neighborhood_finder(input_node, node2_categories, resources):
            raise AttributeError("'NoneType' object has no attribute 'types'")

    monkeypatch.setattr(ta, "TCT", _BoomTCT)
    gw = ta.TctGateway(resources=object())
    with pytest.raises(ValueError, match="node normalizer"):
        gw.neighborhood("BOGUS:123", ["biolink:Gene"])


# --------------------------------------------------------------------------------------
# finders call through to TCT with the right arguments
# --------------------------------------------------------------------------------------
class _FakeTCT:
    @staticmethod
    def Neighborhood_finder(input_node, node2_categories, resources):
        return {"kind": "nb", "input": input_node, "cats": node2_categories, "res": resources}

    @staticmethod
    def Path_finder(input_node1, input_node2, intermediate_categories, resources):
        return {"kind": "path", "n1": input_node1, "n2": input_node2, "inter": intermediate_categories}


def test_neighborhood_defaults_to_drug_categories(monkeypatch):
    monkeypatch.setattr(ta, "TCT", _FakeTCT)
    gw = ta.TctGateway(resources="RES")
    r = gw.neighborhood("NCBIGene:2322")
    assert r["input"] == "NCBIGene:2322"
    assert r["cats"] == ta.DEFAULT_TARGET_CATEGORIES
    assert r["res"] == "RES"


def test_path_passes_intermediates(monkeypatch):
    monkeypatch.setattr(ta, "TCT", _FakeTCT)
    gw = ta.TctGateway(resources="RES")
    r = gw.path("NCBIGene:596", "CHEBI:1", ["biolink:Gene"])
    assert r["n1"] == "NCBIGene:596" and r["n2"] == "CHEBI:1"
    assert r["inter"] == ["biolink:Gene"]


# --------------------------------------------------------------------------------------
# edge_evidence
# --------------------------------------------------------------------------------------
class _FakeKG:
    def __init__(self, edges):
        self._edges = edges

    def items(self):
        return self._edges.items()

    def to_networkx(self, resolve_names=False, include_attributes=False):  # for _knowledge_graphs guard
        return nx.MultiDiGraph()


class _FakeResult:
    def __init__(self, kg):
        self.knowledge_graph = kg


def test_edge_evidence_aggregates_matching_edges(monkeypatch):
    edges = {
        "e0": {"subject": "A", "object": "B", "predicate": "biolink:treats", "attributes": [1]},
        "e1": {"subject": "A", "object": "C", "predicate": "biolink:affects", "attributes": [2]},
    }
    monkeypatch.setattr(
        ta,
        "extract_rich_edge_attributes",
        lambda attrs: {
            "publications": ["PMID:9"],
            "supporting_text": ["evidence text"],
            "confidence_scores": {"s": 0.9},
        },
    )
    gw = ta.TctGateway(resources=object())
    ev = gw.edge_evidence(_FakeResult(_FakeKG(edges)), "A", "B")
    assert ev["edges_matched"] == 1
    assert ev["predicates"] == ["biolink:treats"]
    assert ev["publications"] == ["PMID:9"]
    assert ev["confidence_scores"] == {"s": 0.9}


# --------------------------------------------------------------------------------------
# tct_result_to_nx
# --------------------------------------------------------------------------------------
class _FakeKGGraph:
    """A KnowledgeGraph that returns a prebuilt TCT-style nx graph."""

    def __init__(self, g):
        self._g = g

    def to_networkx(self, resolve_names=False, include_attributes=False):
        return self._g


def test_tct_result_to_nx_normalizes_neighborhood_result():
    tct_g = nx.MultiDiGraph()
    tct_g.add_node("NCBIGene:2322", label="FLT3")
    tct_g.add_node("PUBCHEM.COMPOUND:1", label="Drug")
    tct_g.add_edge(
        "NCBIGene:2322", "PUBCHEM.COMPOUND:1", key="e0",
        predicate="biolink:interacts_with", primary_sources=["infores:x"],
        aggregator_sources=[], publications=[], supporting_text=[], confidence_scores={},
    )
    result = _FakeResult(_FakeKGGraph(tct_g))
    out = ta.tct_result_to_nx(result, ["NCBIGene:2322"], curie_to_symbol={"NCBIGene:2322": "FLT3"})
    assert out.nodes["NCBIGene:2322"]["is_query_gene"] is True
    assert out.nodes["PUBCHEM.COMPOUND:1"]["category"] == "ChemicalEntity"
    assert out.nodes["PUBCHEM.COMPOUND:1"]["gene_frequency"] == 1
    (_, _, data) = list(out.edges(data=True))[0]
    assert data["predicate"] == "biolink:interacts_with"
