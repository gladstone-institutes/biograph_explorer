"""Tests for GeneSet Translator package.

Tests:
- Package initialization
- Version availability


"""

from geneset_translator import __version__


class TestPackage:
    """Test suite for package-level functionality."""

    def test_version_available(self):
        """Test that package version is accessible."""
        assert __version__ is not None
        assert isinstance(__version__, str)
