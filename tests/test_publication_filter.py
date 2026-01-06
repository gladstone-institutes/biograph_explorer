"""Tests for publication filter edge marking and styling."""
import pytest
import networkx as nx
from biograph_explorer.ui.network_viz import filter_graph_by_publication
from biograph_explorer.utils.biolink_predicates import get_predicate_depths


def create_test_graph_with_publications():
    """Create a MultiDiGraph with contrived publications for testing.

    Graph structure:
    - 2 query genes: GENE_A, GENE_B
    - 2 intermediates: PROTEIN_X, PROTEIN_Y
    - 1 disease: DISEASE_Z

    Publications:
    - PMID:11111 on edge GENE_A -> PROTEIN_X (predicate: physically_interacts_with)
    - PMID:22222 on edge GENE_B -> PROTEIN_Y (predicate: affects)
    - No pub on edge PROTEIN_X -> DISEASE_Z (predicate: associated_with)
    - PMID:11111 on edge PROTEIN_Y -> DISEASE_Z (predicate: biomarker_for)
    """
    G = nx.MultiDiGraph()

    G.add_node("GENE_A", category="Gene", is_query_gene=True, label="Gene A")
    G.add_node("GENE_B", category="Gene", is_query_gene=True, label="Gene B")
    G.add_node("PROTEIN_X", category="Protein", is_query_gene=False, label="Protein X")
    G.add_node("PROTEIN_Y", category="Protein", is_query_gene=False, label="Protein Y")
    G.add_node("DISEASE_Z", category="Disease", is_query_gene=False, label="Disease Z")

    G.add_edge("GENE_A", "PROTEIN_X", key="edge_0",
               predicate="biolink:physically_interacts_with",
               publications=["PMID:11111"])

    G.add_edge("GENE_B", "PROTEIN_Y", key="edge_1",
               predicate="biolink:affects",
               publications=["PMID:22222"])

    G.add_edge("PROTEIN_X", "DISEASE_Z", key="edge_2",
               predicate="biolink:associated_with",
               publications=[])

    G.add_edge("PROTEIN_Y", "DISEASE_Z", key="edge_3",
               predicate="biolink:biomarker_for",
               publications=["PMID:11111"])

    return G


class TestFilterGraphByPublication:
    """Tests for filter_graph_by_publication function."""

    def test_marks_edges_with_publication(self):
        """Edges with the filtered publication should have _has_filtered_pub=True."""
        G = create_test_graph_with_publications()
        query_genes = ["GENE_A", "GENE_B"]

        filtered = filter_graph_by_publication(G, "PMID:11111", query_genes, "DISEASE_Z")

        # Check edge_0 (GENE_A -> PROTEIN_X) - has PMID:11111
        assert filtered["GENE_A"]["PROTEIN_X"]["edge_0"].get("_has_filtered_pub") is True

        # Check edge_3 (PROTEIN_Y -> DISEASE_Z) - has PMID:11111
        assert filtered["PROTEIN_Y"]["DISEASE_Z"]["edge_3"].get("_has_filtered_pub") is True

    def test_marks_edges_without_publication(self):
        """Edges without the filtered publication should have _has_filtered_pub=False."""
        G = create_test_graph_with_publications()
        query_genes = ["GENE_A", "GENE_B"]

        filtered = filter_graph_by_publication(G, "PMID:11111", query_genes, "DISEASE_Z")

        # Check edge_1 (GENE_B -> PROTEIN_Y) - has PMID:22222, not PMID:11111
        if "GENE_B" in filtered and "PROTEIN_Y" in filtered["GENE_B"]:
            assert filtered["GENE_B"]["PROTEIN_Y"]["edge_1"].get("_has_filtered_pub") is False

        # Check edge_2 (PROTEIN_X -> DISEASE_Z) - no publications
        if "PROTEIN_X" in filtered and "DISEASE_Z" in filtered["PROTEIN_X"]:
            assert filtered["PROTEIN_X"]["DISEASE_Z"]["edge_2"].get("_has_filtered_pub") is False

    def test_pmid_case_normalization(self):
        """PMID case differences should match correctly."""
        G = nx.MultiDiGraph()
        G.add_node("A", is_query_gene=True, category="Gene")
        G.add_node("B", is_query_gene=False, category="Protein")
        G.add_edge("A", "B", key="e1", predicate="biolink:affects", publications=["pmid:12345"])

        # Filter with uppercase format
        filtered = filter_graph_by_publication(G, "PMID:12345", ["A"], None)

        # Should still match after case normalization
        assert filtered["A"]["B"]["e1"].get("_has_filtered_pub") is True


class TestPriorityEdgeLabelSelection:
    """Tests for priority_edge_label computation logic."""

    def compute_priority_label(self, graph, pub_filter=None):
        """Replicate the priority_label computation logic from app.py."""
        all_predicates = set()
        filtered_pub_predicates = set()

        for _, _, _, data in graph.edges(keys=True, data=True):
            pred = data.get("predicate", "").replace("biolink:", "")
            if pred:
                all_predicates.add(pred)
                if data.get('_has_filtered_pub'):
                    filtered_pub_predicates.add(pred)

        predicates_to_use = filtered_pub_predicates if (pub_filter and filtered_pub_predicates) else all_predicates

        if predicates_to_use:
            depths = get_predicate_depths()
            sorted_predicates = sorted(
                predicates_to_use,
                key=lambda p: (-depths.get(p.lower().replace(" ", "_"), 0), p)
            )
            return sorted_predicates[0]
        return None

    def test_without_pub_filter_uses_all_predicates(self):
        """Without publication filter, should use the most specific predicate overall."""
        G = create_test_graph_with_publications()

        # No publication filter
        priority = self.compute_priority_label(G, pub_filter=None)

        # Should pick from all predicates based on depth
        assert priority is not None
        # All predicates in the graph
        all_preds = {"physically_interacts_with", "affects", "associated_with", "biomarker_for"}
        assert priority in all_preds

    def test_with_pub_filter_uses_only_filtered_predicates(self):
        """With publication filter, should only use predicates from edges with that pub."""
        G = create_test_graph_with_publications()
        query_genes = ["GENE_A", "GENE_B"]

        filtered = filter_graph_by_publication(G, "PMID:11111", query_genes, "DISEASE_Z")
        priority = self.compute_priority_label(filtered, pub_filter="PMID:11111")

        # Priority should be one of the predicates on edges with PMID:11111
        valid_predicates = {"physically_interacts_with", "biomarker_for"}
        assert priority in valid_predicates, f"Expected one of {valid_predicates}, got {priority}"

    def test_pub_filter_with_no_matching_edges_falls_back(self):
        """If pub_filter is set but no edges match, should fall back to all predicates."""
        G = create_test_graph_with_publications()

        # Manually set no edges as having filtered pub
        for u, v, key in G.edges(keys=True):
            G[u][v][key]['_has_filtered_pub'] = False

        # Even with pub_filter set, should use all predicates since none match
        priority = self.compute_priority_label(G, pub_filter="PMID:99999")

        # Should return a predicate (not None)
        assert priority is not None


class TestHighlightedEdgeCount:
    """Tests for counting edges marked with _has_filtered_pub."""

    def test_count_matches_expected(self):
        """Highlighted count should match number of edges with the publication."""
        G = create_test_graph_with_publications()
        query_genes = ["GENE_A", "GENE_B"]

        filtered = filter_graph_by_publication(G, "PMID:11111", query_genes, "DISEASE_Z")

        highlighted = sum(1 for _, _, _, d in filtered.edges(keys=True, data=True)
                        if d.get('_has_filtered_pub'))

        # Should be 2 edges with PMID:11111
        assert highlighted == 2, f"Expected 2 highlighted edges, got {highlighted}"

    def test_total_edges_preserved(self):
        """All edges should have _has_filtered_pub attribute set."""
        G = create_test_graph_with_publications()
        query_genes = ["GENE_A", "GENE_B"]

        filtered = filter_graph_by_publication(G, "PMID:11111", query_genes, "DISEASE_Z")

        # All edges in filtered graph should have _has_filtered_pub attribute
        for u, v, key, data in filtered.edges(keys=True, data=True):
            assert '_has_filtered_pub' in data, f"Edge {u}->{v} ({key}) missing _has_filtered_pub"
