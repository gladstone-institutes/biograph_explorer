"""Tests for agent.graph_clustering.summarize_graph: topology-aware structured summaries.

The key behaviors: community detection ONLY on a real mesh (after stripping universal hubs), and a
graceful fallback to category/predicate facets on hub-dominated/star graphs (the user's constraint:
modularity is meaningless when every edge funnels through one disease/gene hub)."""

import networkx as nx

from geneset_translator.agent.graph_clustering import (
    MAX_CLUSTERS_REPORTED,
    TOP_MEMBERS_PER_CLUSTER,
    summarize_graph,
)


def _mesh_two_communities():
    """Two dense gene cliques joined by a single bridge edge."""
    g = nx.MultiDiGraph()
    a = [f"NCBIGene:{i}" for i in range(6)]
    b = [f"NCBIGene:{100 + i}" for i in range(6)]
    for grp in (a, b):
        for i in range(len(grp)):
            for j in range(i + 1, len(grp)):
                g.add_edge(grp[i], grp[j], predicate="biolink:interacts_with")
    g.add_edge(a[0], b[0], predicate="biolink:interacts_with")
    for n in g.nodes:
        g.nodes[n]["category"] = "biolink:Gene"
    return g, a, b


def test_mesh_uses_community_detection():
    g, a, b = _mesh_two_communities()
    digest, node_cluster = summarize_graph(g, query_genes=[])
    assert digest["topology"] == "mesh"
    assert "louvain" in digest["method"]
    # the two cliques land in different communities
    non_hub = [n for n in a[1:] ]  # interior nodes of clique A (not the bridge hub a[0])
    clusters_of_a = {node_cluster[n] for n in non_hub if n in node_cluster}
    clusters_of_b = {node_cluster[n] for n in b[1:] if n in node_cluster}
    assert clusters_of_a and clusters_of_b and clusters_of_a.isdisjoint(clusters_of_b)


def test_star_is_hub_dominated_and_faceted_not_clustered():
    g = nx.MultiDiGraph()
    g.add_node("MONDO:1", category="biolink:Disease")
    for i in range(40):
        n = f"NCBIGene:{i}"
        g.add_node(n, category="biolink:Gene")
        g.add_edge("MONDO:1", n, predicate="biolink:associated_with")
    digest, _ = summarize_graph(g, query_genes=[], disease_curie="MONDO:1")
    assert digest["topology"] == "hub_dominated"
    assert "facet" in digest["method"]
    assert digest["notes"]  # explains the fallback
    assert digest["category_facets"].get("Gene") == 40
    assert digest["predicate_facets"].get("associated_with") == 40


def test_members_are_named_from_label_and_clusters_carry_predicates():
    """Non-query members (proteins, concepts, chemicals) must be reported by their resolved name
    (the ``label`` attr viz_normalizer sets), not bare CURIEs, and each cluster carries the predicates
    that connect it -- so the summary can describe ALL entity types usefully."""
    g = nx.MultiDiGraph()
    g.add_node("MONDO:1", category="biolink:Disease", label="acute myeloid leukemia")
    # one named protein and one named UMLS concept hanging off the disease hub
    g.add_node("UniProtKB:P31994", category="biolink:Protein", label="Fc receptor FCGR2B")
    g.add_node("UMLS:C0007634", category="biolink:Other", label="Cells")
    for i in range(20):
        n = f"NCBIGene:{i}"
        g.add_node(n, category="biolink:Gene", label=f"GENE{i}")
        g.add_edge("MONDO:1", n, predicate="biolink:associated_with")
    g.add_edge("MONDO:1", "UniProtKB:P31994", predicate="biolink:affects")
    g.add_edge("MONDO:1", "UMLS:C0007634", predicate="biolink:coexists_with")

    digest, _ = summarize_graph(g, query_genes=[], disease_curie="MONDO:1")
    all_members = [m for c in digest["clusters"] for m in c["top_members"]]
    assert "Fc receptor FCGR2B" in all_members  # named, not "UniProtKB:P31994"
    assert "Cells" in all_members               # named, not "UMLS:C0007634"
    assert not any(m.startswith(("UniProtKB:", "UMLS:", "NCBIGene:", "MONDO:")) for m in all_members)
    assert all("dominant_predicates" in c for c in digest["clusters"])
    assert any(c["dominant_predicates"] for c in digest["clusters"])


def test_query_gene_neighborhood_star_does_not_fabricate_modules():
    """A single query gene expanded to many drugs is a star around that gene -> facets, no communities."""
    g = nx.MultiDiGraph()
    g.add_node("NCBIGene:2322", category="biolink:Gene")
    for i in range(50):
        d = f"CHEBI:{i}"
        g.add_node(d, category="biolink:ChemicalEntity")
        g.add_edge("NCBIGene:2322", d, predicate="biolink:affects")
    digest, _ = summarize_graph(g, query_genes=["NCBIGene:2322"])
    assert digest["topology"] == "hub_dominated"
    assert digest["category_facets"].get("ChemicalEntity") == 50


def test_fragmented_graph_clusters_are_components():
    g = nx.MultiDiGraph()
    for k in range(5):
        a, b = f"X:{k}a", f"X:{k}b"
        g.add_edge(a, b, predicate="biolink:related_to")
        g.nodes[a]["category"] = "biolink:Gene"
        g.nodes[b]["category"] = "biolink:Gene"
    digest, node_cluster = summarize_graph(g, query_genes=[])
    assert digest["topology"] == "fragmented"
    assert digest["n_clusters"] == 5
    assert len(set(node_cluster.values())) == 5


def test_digest_is_compact_and_uses_symbols():
    g, a, b = _mesh_two_communities()
    c2s = {n: f"SYM{i}" for i, n in enumerate(a + b)}
    digest, _ = summarize_graph(g, query_genes=[], curie_to_symbol=c2s)
    assert len(digest["clusters"]) <= MAX_CLUSTERS_REPORTED + 1  # +1 for the 'hubs' group
    for c in digest["clusters"]:
        assert len(c["top_members"]) <= TOP_MEMBERS_PER_CLUSTER
        assert all(m.startswith("SYM") for m in c["top_members"])  # mapped to symbols


def test_clustering_is_deterministic():
    g, _, _ = _mesh_two_communities()
    d1, nc1 = summarize_graph(g, query_genes=[], seed=0)
    d2, nc2 = summarize_graph(g, query_genes=[], seed=0)
    assert nc1 == nc2
    assert [c["size"] for c in d1["clusters"]] == [c["size"] for c in d2["clusters"]]
