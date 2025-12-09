"""Tests for BioGraph Explorer package.

Tests:
- Package initialization
- Version availability

Phase 2 Status: Stub created
TODO: Add integration tests if needed
"""

from biograph_explorer import __version__


class TestPackage:
    """Test suite for package-level functionality."""

    def test_version_available(self):
        """Test that package version is accessible."""
        assert __version__ is not None
        assert isinstance(__version__, str)
