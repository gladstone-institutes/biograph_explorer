"""Tests for TRAPI client functionality.

Tests for the pathfinder query: Gene → [Intermediate] → Disease.
"""

import pytest
from pathlib import Path
from unittest.mock import MagicMock, patch

from geneset_translator.utils.biolink_predicates import (
    GRANULARITY_PRESETS,
    filter_predicates_by_granularity,
)


class TestGranularityPresets:
    """Tests for predicate granularity presets."""

    def test_presets_defined(self):
        """Granularity presets are defined."""
        assert GRANULARITY_PRESETS is not None
        assert len(GRANULARITY_PRESETS) > 0

    def test_preset_levels(self):
        """Preset levels are correctly defined."""
        assert "All Relationships" in GRANULARITY_PRESETS
        assert "Standard" in GRANULARITY_PRESETS
        assert "Specific Only" in GRANULARITY_PRESETS

    def test_preset_min_depths(self):
        """Preset min_depth values are correctly ordered."""
        assert GRANULARITY_PRESETS["All Relationships"]["min_depth"] == 1
        assert GRANULARITY_PRESETS["Standard"]["min_depth"] == 2
        assert GRANULARITY_PRESETS["Specific Only"]["min_depth"] == 3


class TestFilterPredicatesByGranularity:
    """Tests for predicate filtering by granularity level."""

    def test_filter_preserves_specific_predicates(self):
        """Filter should preserve specific predicates at any level."""
        predicates = ["biolink:causes", "biolink:affects", "biolink:treats"]
        filtered = filter_predicates_by_granularity(predicates, min_depth=2)
        # Specific predicates should pass through
        assert len(filtered) > 0

    def test_filter_excludes_literature_when_requested(self):
        """Filter should exclude literature predicates when exclude_literature=True."""
        predicates = [
            "biolink:causes",
            "biolink:occurs_together_in_literature_with",
        ]
        filtered = filter_predicates_by_granularity(
            predicates, min_depth=1, exclude_literature=True
        )
        assert "biolink:occurs_together_in_literature_with" not in filtered

    def test_filter_excludes_coexpression_when_requested(self):
        """Filter should exclude coexpression predicates when exclude_coexpression=True."""
        predicates = ["biolink:causes", "biolink:coexpressed_with"]
        filtered = filter_predicates_by_granularity(
            predicates, min_depth=1, exclude_coexpression=True
        )
        assert "biolink:coexpressed_with" not in filtered

    def test_filter_excludes_homology_when_requested(self):
        """Filter should exclude homology predicates when exclude_homology=True."""
        predicates = [
            "biolink:causes",
            "biolink:homologous_to",
            "biolink:orthologous_to",
        ]
        filtered = filter_predicates_by_granularity(
            predicates, min_depth=1, exclude_homology=True
        )
        assert "biolink:homologous_to" not in filtered
        assert "biolink:orthologous_to" not in filtered


class TestTRAPIClient:
    """Tests for TRAPIClient functionality.

    Note: These are unit tests that don't make actual API calls.
    Integration tests would require mocking the TCT library.
    """

    @pytest.fixture
    def mock_client(self):
        """Create a mock TRAPIClient for testing."""
        # Import here to avoid issues if TCT is not installed
        try:
            from geneset_translator.core.trapi_client import TRAPIClient

            with patch.object(TRAPIClient, "_load_translator_resources"):
                client = TRAPIClient(cache_dir=Path("/tmp/test_cache"))
                client.metaKG = MagicMock()
                client.APInames = {}
                client.Translator_KP_info = {}
                return client
        except ImportError:
            pytest.skip("TCT library not available")

    def test_client_initialization(self, mock_client):
        """Test client initializes correctly."""
        assert mock_client is not None
        assert mock_client.cache_dir is not None

    def test_cache_response(self, mock_client, tmp_path):
        """Test caching of query responses."""
        from geneset_translator.core.trapi_client import TRAPIResponse

        mock_client.cache_dir = tmp_path

        response = TRAPIResponse(
            query_id="test_query",
            input_genes=["NCBIGene:123"],
            target_disease="MONDO:0004975",
            edges=[],
            metadata={},
            apis_queried=10,
            apis_succeeded=5,
        )

        cache_file = mock_client._cache_response(response)

        assert cache_file.exists()
        assert cache_file.suffix == ".json"

    def test_list_cached_queries_empty(self, mock_client, tmp_path):
        """Test listing cached queries when none exist."""
        mock_client.cache_dir = tmp_path

        results = mock_client.list_cached_queries()
        assert results == []
