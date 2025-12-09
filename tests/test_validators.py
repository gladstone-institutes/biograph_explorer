"""Tests for input validators."""

import pytest
from biograph_explorer.utils import validate_gene_list, validate_disease_curie, ValidationError
from biograph_explorer.utils.validators import validate_curie, validate_convergence_threshold


class TestGeneListValidation:
    """Test suite for gene list validation."""

    def test_validate_valid_gene_list(self):
        """Test validation of valid gene list."""
        result = validate_gene_list(["APOE", "APP", "PSEN1"])
        assert result == ["APOE", "APP", "PSEN1"]

    def test_validate_empty_gene_list(self):
        """Test that empty gene list raises ValidationError."""
        with pytest.raises(ValidationError):
            validate_gene_list([])

    def test_validate_too_many_genes(self):
        """Test that gene list exceeding max_genes raises ValidationError."""
        genes = [f"GENE{i}" for i in range(101)]
        with pytest.raises(ValidationError):
            validate_gene_list(genes, max_genes=100)

    def test_gene_list_trimming(self):
        """Test that gene symbols are trimmed and cleaned."""
        result = validate_gene_list([" APOE ", " APP "])
        assert result == ["APOE", "APP"]

    def test_gene_list_uppercasing(self):
        """Test that gene symbols are converted to uppercase."""
        result = validate_gene_list(["apoe", "App"])
        assert result == ["APOE", "APP"]

    def test_duplicate_genes(self):
        """Test handling of duplicate gene symbols."""
        result = validate_gene_list(["APOE", "APOE", "APP"])
        assert result == ["APOE", "APP"]


class TestDiseaseCURIEValidation:
    """Test suite for disease CURIE validation."""

    def test_validate_valid_mondo_curie(self):
        """Test validation of valid MONDO CURIE."""
        result = validate_disease_curie("MONDO:0004975")
        assert result == "MONDO:0004975"

    def test_validate_valid_doid_curie(self):
        """Test validation of valid DOID CURIE."""
        result = validate_disease_curie("DOID:10652")
        assert result == "DOID:10652"

    def test_validate_invalid_curie_format(self):
        """Test that invalid CURIE format raises ValidationError."""
        with pytest.raises(ValidationError):
            validate_disease_curie("not-a-curie")

    def test_validate_empty_curie(self):
        """Test that empty string raises ValidationError."""
        with pytest.raises(ValidationError):
            validate_disease_curie("")


class TestCURIEValidation:
    """Test suite for generic CURIE validation."""

    def test_valid_curie_formats(self):
        """Test various valid CURIE formats."""
        valid_curies = [
            "NCBIGene:1803",
            "MONDO:0004975",
            "DOID:10652",
            "HGNC:1234",
            "UniProtKB:P12345",
        ]
        for curie in valid_curies:
            assert validate_curie(curie) is True

    def test_invalid_curie_formats(self):
        """Test various invalid CURIE formats."""
        invalid_curies = [
            "not-a-curie",
            ":12345",
            "PREFIX:",
            "12345",
        ]
        for curie in invalid_curies:
            assert validate_curie(curie) is False


class TestParameterValidation:
    """Test suite for parameter validation."""

    def test_convergence_threshold_validation(self):
        """Test convergence threshold validation."""
        assert validate_convergence_threshold(2) == 2

    def test_negative_threshold(self):
        """Test that negative threshold raises ValidationError."""
        with pytest.raises(ValidationError):
            validate_convergence_threshold(0)

    def test_excessive_threshold(self):
        """Test that excessively high threshold raises ValidationError."""
        with pytest.raises(ValidationError):
            validate_convergence_threshold(51)
