"""Application settings and configuration management.

Manages:
- TRAPI API endpoints and rate limits
- Cache directories and TTL settings
- Visualization defaults (max nodes, edge sampling)
- Claude API configuration (future Phase 3)

Phase 2 Status: Stub created
"""

from pathlib import Path
from typing import Optional
from pydantic import BaseModel, Field


class Settings(BaseModel):
    """Application configuration settings."""

    # Data directories
    data_dir: Path = Field(default=Path("data"), description="Root data directory")
    cache_dir: Path = Field(default=Path("data/cache"), description="TRAPI response cache")
    sessions_dir: Path = Field(default=Path("data/sessions"), description="Pickled graphs")
    exports_dir: Path = Field(default=Path("data/exports"), description="HTML/PNG exports")
    logs_dir: Path = Field(default=Path("data/logs"), description="API query logs")

    # TRAPI settings
    trapi_timeout: int = Field(default=30, description="API timeout in seconds")
    trapi_rate_limit: int = Field(default=10, description="Max queries per second")
    trapi_max_workers: int = Field(default=5, description="Parallel query workers")
    trapi_cache_ttl: int = Field(default=86400 * 7, description="Cache TTL in seconds (default: 7 days)")

    # Disease-BiologicalProcess query settings
    disease_bp_timeout: int = Field(
        default=600, description="Timeout for Disease→BP queries in seconds (extended for 2-stage pattern)"
    )
    use_cached_bp_default: bool = Field(
        default=True, description="Default: use cached BiologicalProcesses if available"
    )

    # Graph settings
    max_viz_edges: int = Field(default=50, description="Max edges to visualize")
    min_edges_per_gene: int = Field(default=2, description="Min edges per query gene in viz")
    convergence_threshold: int = Field(default=2, description="Min gene frequency for convergence")

    # Clustering settings
    clustering_algorithm: str = Field(default="louvain", description="Community detection algorithm")
    centrality_metrics: list[str] = Field(
        default=["pagerank", "betweenness", "degree"], description="Centrality metrics to compute"
    )

    # Claude API settings (Phase 3)
    claude_api_key: Optional[str] = Field(default=None, description="Anthropic API key")
    claude_model: str = Field(default="claude-haiku-4", description="Claude model to use")
    claude_max_tokens: int = Field(default=1000, description="Max tokens per response")

    # UI settings
    streamlit_theme: str = Field(default="light", description="Streamlit theme")
    show_debug_info: bool = Field(default=False, description="Show debug information in UI")

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        # Ensure directories exist
        self.data_dir.mkdir(exist_ok=True)
        self.cache_dir.mkdir(exist_ok=True, parents=True)
        self.sessions_dir.mkdir(exist_ok=True, parents=True)
        self.exports_dir.mkdir(exist_ok=True, parents=True)
        self.logs_dir.mkdir(exist_ok=True, parents=True)


# Singleton instance
_settings: Optional[Settings] = None


def get_settings() -> Settings:
    """Get or create the application settings singleton.

    Returns:
        Settings instance with default or environment-configured values
    """
    global _settings
    if _settings is None:
        _settings = Settings()
    return _settings


def reset_settings() -> None:
    """Reset settings singleton (useful for testing)."""
    global _settings
    _settings = None
