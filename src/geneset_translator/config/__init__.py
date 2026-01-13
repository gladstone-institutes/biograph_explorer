"""Configuration module for GeneSet Translator.

Handles application settings, logging configuration, and environment variables.

Phase 2 Status: Stub created
"""

from .settings import Settings, get_settings
from .logging_config import setup_logging

__all__ = ["Settings", "get_settings", "setup_logging"]
