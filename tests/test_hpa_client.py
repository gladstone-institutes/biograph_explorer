"""Tests for Human Protein Atlas client functionality.

Tests for HPAClient and HPAAnnotation classes.
"""

import pytest
import networkx as nx
from unittest.mock import MagicMock, patch

from geneset_translator.core.hpa_client import HPAClient, HPAAnnotation


class TestHPAAnnotation:
    """Tests for HPAAnnotation dataclass."""

    def test_from_hpa_response_basic(self):
        """Test creating HPAAnnotation from HPA API response."""
        response = {
            "Gene": "APOE",
            "RNA single cell type specificity": "Cell type enhanced",
            "RNA single cell type distribution": "Detected in many",
            "RNA single cell type specific nCPM": {
                "Hepatocytes": 9962.2,
                "Melanocytes": 3427.4,
                "Macrophages": 1733.0,
            },
            "Disease involvement": ["Alzheimer disease", "Neurodegeneration"],
            "Protein class": ["Plasma proteins", "Disease related genes"],
            "RNA tissue specificity": "Tissue enhanced",
            "RNA tissue specific nTPM": {
                "liver": 6533.9,
                "brain": 2714.9,
            },
        }

        annotation = HPAAnnotation.from_hpa_response("ENSG00000130203", response)

        # Basic fields
        assert annotation.ensembl_id == "ENSG00000130203"
        assert annotation.gene_symbol == "APOE"

        # Single cell specificity
        assert annotation.cell_type_specificity == "Cell type enhanced"
        assert annotation.cell_type_distribution == "Detected in many"
        assert annotation.cell_type_ncpm["Hepatocytes"] == 9962.2
        assert len(annotation.top_cell_types) == 3
        assert annotation.top_cell_types[0][0] == "Hepatocytes"  # Highest first

        # Disease involvement
        assert "Alzheimer disease" in annotation.disease_involvement
        assert "Neurodegeneration" in annotation.disease_involvement

        # Protein class
        assert "Plasma proteins" in annotation.protein_class
        assert "Disease related genes" in annotation.protein_class

        # Tissue specificity
        assert annotation.tissue_specificity == "Tissue enhanced"
        assert annotation.tissue_ntpm["liver"] == 6533.9
        assert len(annotation.top_tissues) == 2
        assert annotation.top_tissues[0][0] == "liver"  # Highest first

    def test_from_hpa_response_immune_cell_data(self):
        """Test HPAAnnotation parses immune/blood cell specificity data."""
        response = {
            "Gene": "TNF",
            "RNA blood cell specificity": "Immune cell enhanced",
            "RNA blood cell distribution": "Detected in many",
            "RNA blood cell specific nTPM": {
                "non-classical monocyte": 73.4,
                "classical monocyte": 45.2,
                "NK-cell": 12.1,
            },
            "RNA blood lineage specificity": "Group enriched",
            "RNA blood lineage specific nTPM": {
                "monocytes": 73.4,
                "T-cells": 41.4,
                "granulocytes": 22.7,
            },
        }

        annotation = HPAAnnotation.from_hpa_response("ENSG00000232810", response)

        # Immune cell specificity
        assert annotation.immune_cell_specificity == "Immune cell enhanced"
        assert annotation.immune_cell_distribution == "Detected in many"
        assert annotation.immune_cell_ntpm["non-classical monocyte"] == 73.4
        assert len(annotation.top_immune_cells) == 3
        assert annotation.top_immune_cells[0][0] == "non-classical monocyte"  # Highest first

        # Immune lineage specificity
        assert annotation.immune_lineage_specificity == "Group enriched"
        assert len(annotation.top_immune_lineages) == 3
        assert annotation.top_immune_lineages[0][0] == "monocytes"  # Highest first

    def test_from_hpa_response_missing_data(self):
        """Test HPAAnnotation handles missing fields gracefully."""
        response = {"Gene": "TEST"}

        annotation = HPAAnnotation.from_hpa_response("ENSG00000000001", response)

        assert annotation.ensembl_id == "ENSG00000000001"
        assert annotation.gene_symbol == "TEST"
        assert annotation.cell_type_specificity is None
        assert annotation.cell_type_ncpm == {}
        assert annotation.top_cell_types == []
        assert annotation.disease_involvement == []
        assert annotation.protein_class == []
        assert annotation.tissue_specificity is None
        assert annotation.tissue_ntpm == {}
        assert annotation.top_tissues == []
        # Immune cell fields should also be empty/None
        assert annotation.immune_cell_specificity is None
        assert annotation.immune_cell_ntpm == {}
        assert annotation.top_immune_cells == []
        assert annotation.immune_lineage_specificity is None
        assert annotation.top_immune_lineages == []


class TestHPAClientEnsemblLookup:
    """Tests for HPAClient Ensembl ID lookup functionality."""

    def test_get_ensembl_ids_success(self):
        """Test Ensembl ID lookup returns valid results."""
        client = HPAClient()

        # Test with known genes
        result = client.get_ensembl_ids(["348", "7157"])  # APOE, TP53

        # Result is now Dict[str, List[str]] to handle genes with multiple Ensembl IDs
        assert "348" in result
        assert isinstance(result["348"], list)
        assert len(result["348"]) >= 1
        assert result["348"][0].startswith("ENSG")
        assert "7157" in result
        assert isinstance(result["7157"], list)
        assert result["7157"][0].startswith("ENSG")

    def test_get_ensembl_ids_empty_list(self):
        """Test empty input returns empty dict."""
        client = HPAClient()
        result = client.get_ensembl_ids([])
        assert result == {}

    def test_get_ensembl_ids_invalid_gene(self):
        """Test invalid gene IDs are handled gracefully."""
        client = HPAClient()
        result = client.get_ensembl_ids(["INVALID_123456789"])
        assert "INVALID_123456789" not in result


class TestHPAClientFetch:
    """Tests for HPAClient cell type specificity fetch."""

    def test_fetch_cell_type_specificity_success(self):
        """Test fetching HPA data for known genes."""
        client = HPAClient()

        # APOE - known to have cell type data
        result = client.fetch_cell_type_specificity(["ENSG00000130203"])

        assert "ENSG00000130203" in result
        annotation = result["ENSG00000130203"]
        assert annotation.gene_symbol == "APOE"
        assert annotation.cell_type_specificity is not None
        assert len(annotation.cell_type_ncpm) > 0

    def test_fetch_cell_type_specificity_empty_list(self):
        """Test empty input returns empty dict."""
        client = HPAClient()
        result = client.fetch_cell_type_specificity([])
        assert result == {}

    def test_fetch_cell_type_specificity_invalid_ensembl(self):
        """Test invalid Ensembl IDs are handled gracefully."""
        client = HPAClient()
        result = client.fetch_cell_type_specificity(["ENSG99999999999"])
        assert "ENSG99999999999" not in result


class TestHPAClientGraphAnnotation:
    """Tests for HPAClient graph annotation functionality."""

    def test_annotate_graph_adds_hpa_data(self):
        """Test that HPA annotations are added to graph nodes."""
        client = HPAClient()

        # Create a simple test graph with gene nodes
        G = nx.DiGraph()
        G.add_node("NCBIGene:348", label="APOE", category="Gene")
        G.add_node("NCBIGene:7157", label="TP53", category="Gene")
        G.add_node("GO:0006915", label="apoptotic process", category="BiologicalProcess")

        # Initialize annotation_features for gene nodes
        G.nodes["NCBIGene:348"]["annotation_features"] = {}
        G.nodes["NCBIGene:7157"]["annotation_features"] = {}

        # Annotate
        annotated_graph, metadata = client.annotate_graph(G)

        # Check metadata
        assert "hpa_annotated_count" in metadata
        assert metadata["hpa_total_genes"] == 2  # Only gene nodes

        # Check at least one gene was annotated
        assert metadata["hpa_annotated_count"] >= 1

        # Check APOE has HPA data (including expanded fields)
        apoe_features = annotated_graph.nodes["NCBIGene:348"].get("annotation_features", {})
        if "hpa_ensembl_id" in apoe_features:
            assert apoe_features["hpa_ensembl_id"].startswith("ENSG")
            assert "hpa_cell_type_specificity" in apoe_features
            # Check expanded fields
            assert "hpa_top_cell_types" in apoe_features
            assert len(apoe_features["hpa_top_cell_types"]) <= 5  # Top 5
            assert "hpa_disease_involvement" in apoe_features
            assert "hpa_protein_class" in apoe_features
            assert "hpa_tissue_specificity" in apoe_features
            assert "hpa_top_tissues" in apoe_features
            assert len(apoe_features["hpa_top_tissues"]) <= 3  # Top 3

    def test_annotate_graph_no_gene_nodes(self):
        """Test graph with no gene nodes returns zero annotations."""
        client = HPAClient()

        # Create graph with only non-gene nodes
        G = nx.DiGraph()
        G.add_node("GO:0006915", label="apoptotic process")
        G.add_node("MONDO:0004975", label="Alzheimer disease")

        annotated_graph, metadata = client.annotate_graph(G)

        assert metadata["hpa_annotated_count"] == 0
        assert metadata["hpa_total_genes"] == 0


class TestHPAClientUniProtKBMapping:
    """Tests for HPAClient UniProtKB to Ensembl ID mapping."""

    def test_get_ensembl_ids_from_uniprot_success(self):
        """Test UniProtKB ID lookup returns valid results."""
        client = HPAClient()

        # Test with known proteins (APP = P05067, CRYAA = P02489)
        result = client.get_ensembl_ids_from_uniprot(["P05067", "P02489"])

        assert "P05067" in result  # APP
        assert result["P05067"].startswith("ENSG")
        assert "P02489" in result  # CRYAA
        assert result["P02489"].startswith("ENSG")

    def test_get_ensembl_ids_from_uniprot_empty_list(self):
        """Test empty input returns empty dict."""
        client = HPAClient()
        result = client.get_ensembl_ids_from_uniprot([])
        assert result == {}

    def test_get_ensembl_ids_from_uniprot_invalid_id(self):
        """Test invalid UniProtKB IDs are handled gracefully."""
        client = HPAClient()
        result = client.get_ensembl_ids_from_uniprot(["INVALID_UNIPROT_999"])
        assert "INVALID_UNIPROT_999" not in result


class TestHPAClientProteinAnnotation:
    """Tests for HPAClient protein node annotation functionality."""

    def test_annotate_graph_protein_nodes(self):
        """Test that HPA annotations are added to protein nodes."""
        client = HPAClient()

        # Create a test graph with protein nodes
        # P05067 = APP (Amyloid precursor protein)
        G = nx.DiGraph()
        G.add_node("UniProtKB:P05067", label="APP", category="Protein")
        G.add_node("GO:0006915", label="apoptotic process", category="BiologicalProcess")

        # Initialize annotation_features
        G.nodes["UniProtKB:P05067"]["annotation_features"] = {}

        # Annotate
        annotated_graph, metadata = client.annotate_graph(G)

        # Check metadata includes protein counts
        assert "hpa_total_proteins" in metadata
        assert metadata["hpa_total_proteins"] == 1
        assert metadata["hpa_total_genes"] == 0

        # Check protein was annotated
        if metadata["hpa_proteins_annotated"] > 0:
            app_features = annotated_graph.nodes["UniProtKB:P05067"].get("annotation_features", {})
            assert "hpa_ensembl_id" in app_features
            assert app_features["hpa_ensembl_id"].startswith("ENSG")

    def test_annotate_graph_mixed_gene_and_protein_nodes(self):
        """Test that HPA annotations work with both gene and protein nodes."""
        client = HPAClient()

        # Create a test graph with both gene and protein nodes
        G = nx.DiGraph()
        G.add_node("NCBIGene:348", label="APOE", category="Gene")  # APOE gene
        G.add_node("UniProtKB:P05067", label="APP", category="Protein")  # APP protein
        G.add_node("GO:0006915", label="apoptotic process", category="BiologicalProcess")

        # Initialize annotation_features
        G.nodes["NCBIGene:348"]["annotation_features"] = {}
        G.nodes["UniProtKB:P05067"]["annotation_features"] = {}

        # Annotate
        annotated_graph, metadata = client.annotate_graph(G)

        # Check metadata
        assert metadata["hpa_total_genes"] == 1
        assert metadata["hpa_total_proteins"] == 1
        assert "hpa_genes_annotated" in metadata
        assert "hpa_proteins_annotated" in metadata

        # Total annotated should be sum of genes + proteins annotated
        assert metadata["hpa_annotated_count"] == (
            metadata["hpa_genes_annotated"] + metadata["hpa_proteins_annotated"]
        )
