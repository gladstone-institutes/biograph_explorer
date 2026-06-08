"""Tests for agent.viz_normalizer: TCT graph -> network_viz attribute shape."""

import networkx as nx
import pytest

from geneset_translator.agent import viz_normalizer as vn


def _fake_tct_graph():
    """Mimic KnowledgeGraph.to_networkx(resolve_names=True, include_attributes=True)."""
    g = nx.MultiDiGraph()
    g.add_node("NCBIGene:2322", label="FLT3")
    g.add_node("PUBCHEM.COMPOUND:153999", label="Ruboxistaurin")
    g.add_node("MONDO:0018874", label="acute myeloid leukemia")
    g.add_edge(
        "NCBIGene:2322",
        "PUBCHEM.COMPOUND:153999",
        key="e0",
        predicate="biolink:physically_interacts_with",
        primary_sources=["infores:foo"],
        aggregator_sources=["infores:bar"],
        publications=["PMID:1"],
        supporting_text=["FLT3 interacts with the compound."],
        confidence_scores={"tmkp_confidence_score": 0.85},
    )
    return g


@pytest.mark.parametrize(
    "curie,expected",
    [
        ("NCBIGene:2322", "Gene"),
        ("HGNC:3765", "Gene"),
        ("UniProtKB:P36888", "Protein"),
        ("PUBCHEM.COMPOUND:153999", "ChemicalEntity"),
        ("CHEMBL.COMPOUND:CHEMBL1", "ChemicalEntity"),
        ("CHEBI:1234", "ChemicalEntity"),
        ("MONDO:0018874", "Disease"),
        ("GO:0006915", "BiologicalProcess"),
        ("REACT:R-HSA-1", "Pathway"),
        ("HP:0001234", "PhenotypicFeature"),
        ("UMLS:C0001", "Other"),
        ("no-colon", "Other"),
    ],
)
def test_classify_category(curie, expected):
    assert vn.classify_category(curie) == expected


def test_normalize_node_attributes():
    g = _fake_tct_graph()
    out = vn.normalize(g, ["NCBIGene:2322"], curie_to_symbol={"NCBIGene:2322": "FLT3"})

    gene = out.nodes["NCBIGene:2322"]
    assert gene["category"] == "Gene"
    assert gene["is_query_gene"] is True
    assert gene["original_symbol"] == "FLT3"
    assert gene["label"] == "FLT3"
    assert gene["curie"] == "NCBIGene:2322"

    drug = out.nodes["PUBCHEM.COMPOUND:153999"]
    assert drug["category"] == "ChemicalEntity"
    assert drug["is_query_gene"] is False
    # convergence: the drug touches one query gene
    assert drug["gene_frequency"] == 1


def test_normalize_edge_attributes_match_renderer_shape():
    g = _fake_tct_graph()
    out = vn.normalize(g, ["NCBIGene:2322"])
    # single edge
    (_, _, data) = list(out.edges(data=True))[0]
    assert data["predicate"] == "biolink:physically_interacts_with"
    # sources reshaped to the TRAPI dict shape network_viz formats
    assert {"resource_id": "infores:foo", "resource_role": "primary_knowledge_source",
            "upstream_resource_ids": []} in data["sources"]
    assert {"resource_id": "infores:bar", "resource_role": "aggregator_knowledge_source",
            "upstream_resource_ids": []} in data["sources"]
    assert data["publications"] == ["PMID:1"]
    assert data["sentences"] == ["FLT3 interacts with the compound."]
    assert data["confidence_scores"] == {"tmkp_confidence_score": 0.85}
    # keys present even when empty so the renderer's .get(...) is well-defined
    assert data["knowledge_level"] is None
    assert data["qualifiers"] == []


def test_publications_coerced_to_strings_and_renderable():
    """Regression: TCT returns some publications as ints; the renderer calls .startswith."""
    g = nx.MultiDiGraph()
    g.add_node("NCBIGene:2322", label="FLT3")
    g.add_node("CHEBI:1", label="drug")
    g.add_edge(
        "NCBIGene:2322", "CHEBI:1", key="e0",
        predicate="biolink:affects",
        primary_sources=["infores:x"], aggregator_sources=[],
        publications=[12345, "PMID:67890", "https://example.com/p"],
        supporting_text=[], confidence_scores={},
    )
    out = vn.normalize(g, ["NCBIGene:2322"])
    (_, _, data) = list(out.edges(data=True))[0]
    assert data["publications"] == ["PMID:12345", "PMID:67890", "https://example.com/p"]

    # the real renderer must not crash on the normalized graph
    from geneset_translator.ui.network_viz import prepare_cytoscape_elements

    elements = prepare_cytoscape_elements(out, ["NCBIGene:2322"])
    assert len(elements["nodes"]) == 2 and len(elements["edges"]) == 1


def test_compute_gene_frequency_counts_distinct_query_genes():
    g = nx.MultiDiGraph()
    g.add_edge("g1", "hub")
    g.add_edge("g2", "hub")
    g.add_edge("g1", "leaf")
    freq = vn.compute_gene_frequency(g, ["g1", "g2"])
    assert freq["hub"] == 2
    assert freq["leaf"] == 1
    assert freq["g1"] == 0
