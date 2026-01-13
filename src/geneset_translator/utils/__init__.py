"""Utility modules for GeneSet Translator.

Provides:
- Input validation for gene lists and disease CURIEs
- Data formatters for display
- Persistence (pickle NetworkX graphs, JSON caching)
- Biolink predicate hierarchy filtering
"""

from .validators import validate_gene_list, validate_disease_curie, ValidationError
from .formatters import format_node_label, format_edge_label
from .biolink_predicates import (
    GRANULARITY_PRESETS,
    filter_predicates_by_granularity,
    get_allowed_predicates_for_display,
    get_excluded_predicates_for_display,
    get_predicate_depths,
    get_predicate_info,
)
from .publication_utils import (
    get_publication_frequency,
    normalize_publication_id,
    format_publication_display,
    validate_publication_extraction,
)

__all__ = [
    "validate_gene_list",
    "validate_disease_curie",
    "ValidationError",
    "format_node_label",
    "format_edge_label",
    "GRANULARITY_PRESETS",
    "filter_predicates_by_granularity",
    "get_allowed_predicates_for_display",
    "get_excluded_predicates_for_display",
    "get_predicate_depths",
    "get_predicate_info",
    "get_publication_frequency",
    "normalize_publication_id",
    "format_publication_display",
    "validate_publication_extraction",
]
