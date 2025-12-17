"""Tests for TRAPI client functionality.

Tests for the new Disease → BiologicalProcess query pattern.
"""

import pytest
from pathlib import Path
from unittest.mock import MagicMock, patch

from biograph_explorer.utils.biolink_predicates import (
    DISEASE_BP_INFORMATIVE_PREDICATES,
    DISEASE_BP_NOISE_PREDICATES,
    DEFAULT_BP_INTERMEDIATE_CATEGORIES,
    filter_disease_bp_predicates,
)


class TestDiseaseBPPredicates:
    """Tests for Disease → BiologicalProcess predicate constants and filtering."""

    def test_informative_predicates_defined(self):
        """Informative predicates list is defined and non-empty."""
        assert DISEASE_BP_INFORMATIVE_PREDICATES is not None
        assert len(DISEASE_BP_INFORMATIVE_PREDICATES) > 0
        # Should include key predicates
        assert "biolink:affects" in DISEASE_BP_INFORMATIVE_PREDICATES
        assert "biolink:causes" in DISEASE_BP_INFORMATIVE_PREDICATES
        assert "biolink:disrupts" in DISEASE_BP_INFORMATIVE_PREDICATES

    def test_noise_predicates_defined(self):
        """Noise predicates list is defined and includes text mining."""
        assert DISEASE_BP_NOISE_PREDICATES is not None
        assert len(DISEASE_BP_NOISE_PREDICATES) > 0
        # Must include the main text mining predicate
        assert "biolink:occurs_together_in_literature_with" in DISEASE_BP_NOISE_PREDICATES

    def test_no_overlap_between_informative_and_noise(self):
        """Informative and noise predicate lists should not overlap."""
        informative_set = set(DISEASE_BP_INFORMATIVE_PREDICATES)
        noise_set = set(DISEASE_BP_NOISE_PREDICATES)
        overlap = informative_set & noise_set
        assert len(overlap) == 0, f"Overlap found: {overlap}"

    def test_filter_disease_bp_predicates_keeps_informative(self):
        """Filter should keep edges with informative predicates."""
        edges = [
            {"subject": "MONDO:123", "object": "GO:456", "predicate": "biolink:affects"},
            {"subject": "MONDO:123", "object": "GO:789", "predicate": "biolink:causes"},
        ]
        filtered = filter_disease_bp_predicates(edges, use_informative_only=True)
        assert len(filtered) == 2

    def test_filter_disease_bp_predicates_removes_noise(self):
        """Filter should remove edges with noise predicates."""
        edges = [
            {"subject": "MONDO:123", "object": "GO:456", "predicate": "biolink:affects"},
            {"subject": "MONDO:123", "object": "GO:789", "predicate": "biolink:occurs_together_in_literature_with"},
            {"subject": "MONDO:123", "object": "GO:012", "predicate": "biolink:coexists_with"},
        ]
        filtered = filter_disease_bp_predicates(edges, use_informative_only=True)
        # Only the 'affects' edge should remain
        assert len(filtered) == 1
        assert filtered[0]["predicate"] == "biolink:affects"

    def test_filter_disease_bp_predicates_disabled(self):
        """When filtering is disabled, all edges should pass through."""
        edges = [
            {"subject": "MONDO:123", "object": "GO:456", "predicate": "biolink:affects"},
            {"subject": "MONDO:123", "object": "GO:789", "predicate": "biolink:occurs_together_in_literature_with"},
        ]
        filtered = filter_disease_bp_predicates(edges, use_informative_only=False)
        assert len(filtered) == 2


class TestDefaultBPIntermediateCategories:
    """Tests for the default intermediate categories for BP queries."""

    def test_default_categories_defined(self):
        """Default intermediate categories list is defined."""
        assert DEFAULT_BP_INTERMEDIATE_CATEGORIES is not None
        assert len(DEFAULT_BP_INTERMEDIATE_CATEGORIES) > 0

    def test_default_categories_include_key_types(self):
        """Default categories should include key intermediate types."""
        # Based on notebook exploration, all these should be viable
        expected = [
            "biolink:ChemicalEntity",
            "biolink:Protein",
            "biolink:Gene",
            "biolink:Pathway",
        ]
        for cat in expected:
            assert cat in DEFAULT_BP_INTERMEDIATE_CATEGORIES, f"Missing: {cat}"

    def test_all_categories_have_biolink_prefix(self):
        """All categories should have biolink: prefix."""
        for cat in DEFAULT_BP_INTERMEDIATE_CATEGORIES:
            assert cat.startswith("biolink:"), f"Missing prefix: {cat}"


class TestTRAPIClientDiseaseBP:
    """Tests for TRAPIClient Disease → BiologicalProcess methods.

    Note: These are unit tests that don't make actual API calls.
    Integration tests would require mocking the TCT library.
    """

    @pytest.fixture
    def mock_client(self):
        """Create a mock TRAPIClient for testing."""
        # Import here to avoid issues if TCT is not installed
        try:
            from biograph_explorer.core.trapi_client import TRAPIClient

            with patch.object(TRAPIClient, '_load_translator_resources'):
                client = TRAPIClient(cache_dir=Path("/tmp/test_cache"))
                client.metaKG = MagicMock()
                client.APInames = {}
                client.Translator_KP_info = {}
                return client
        except ImportError:
            pytest.skip("TCT library not available")

    def test_cache_disease_bp_results(self, mock_client, tmp_path):
        """Test caching of Disease → BP discovery results."""
        mock_client.cache_dir = tmp_path

        bp_curies = ["GO:0006915", "GO:0006914"]
        metadata = {
            "disease_curie": "MONDO:0004975",
            "total_edges_before_filter": 100,
            "total_edges_after_filter": 20,
        }

        cache_file = mock_client._cache_disease_bp_results(
            "MONDO:0004975", bp_curies, metadata
        )

        assert cache_file.exists()
        assert "disease_bp_MONDO_0004975" in cache_file.name
        assert cache_file.suffix == ".json"

    def test_list_cached_disease_bp_results_empty(self, mock_client, tmp_path):
        """Test listing cached results when none exist."""
        mock_client.cache_dir = tmp_path

        results = mock_client.list_cached_disease_bp_results()
        assert results == []

    def test_list_cached_disease_bp_results(self, mock_client, tmp_path):
        """Test listing cached results after caching."""
        mock_client.cache_dir = tmp_path

        # Cache some results first
        bp_curies = ["GO:0006915"]
        metadata = {"disease_curie": "MONDO:0004975"}
        mock_client._cache_disease_bp_results("MONDO:0004975", bp_curies, metadata)

        results = mock_client.list_cached_disease_bp_results()
        assert len(results) == 1
        assert results[0]["disease_curie"] == "MONDO:0004975"
        assert results[0]["bp_count"] == 1

    def test_load_cached_disease_bp_results(self, mock_client, tmp_path):
        """Test loading cached results."""
        mock_client.cache_dir = tmp_path

        # Cache some results first
        bp_curies = ["GO:0006915", "GO:0006914"]
        metadata = {"disease_curie": "MONDO:0004975", "test_key": "test_value"}
        cache_file = mock_client._cache_disease_bp_results("MONDO:0004975", bp_curies, metadata)

        # Load them back
        loaded_curies, loaded_metadata = mock_client.load_cached_disease_bp_results(cache_file)

        assert loaded_curies == bp_curies
        assert loaded_metadata["test_key"] == "test_value"

    def test_load_cached_disease_bp_results_not_found(self, mock_client, tmp_path):
        """Test loading non-existent cache file raises error."""
        mock_client.cache_dir = tmp_path

        with pytest.raises(FileNotFoundError):
            mock_client.load_cached_disease_bp_results(tmp_path / "nonexistent.json")
