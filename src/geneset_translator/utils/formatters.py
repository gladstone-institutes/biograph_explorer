"""Data formatting utilities for display.

Handles:
- Node label formatting (truncate long names, add metadata)
- Edge label formatting (predicates → human-readable)
- Export formatting (HTML, JSON)

Phase 2 Status: Stub created
TODO: Implement formatters
"""

from typing import Dict, List, Any, Optional
import networkx as nx


def format_node_label(
    node_id: str,
    graph: nx.DiGraph,
    max_length: int = 30,
    include_category: bool = False,
) -> str:
    """Format node label for display.

    Args:
        node_id: Node ID
        graph: NetworkX graph containing node
        max_length: Max label length before truncation
        include_category: Whether to include node category

    Returns:
        Formatted label string

    TODO: Implement node label formatting
    """
    raise NotImplementedError("TODO: Implement node label formatting")


def format_edge_label(predicate: str) -> str:
    """Format biolink predicate for human readability.

    Args:
        predicate: biolink predicate (e.g., "biolink:associated_with")

    Returns:
        Human-readable label (e.g., "associated with")
    """
    # Remove biolink prefix and convert underscores to spaces
    if predicate.startswith("biolink:"):
        predicate = predicate[8:]
    return predicate.replace("_", " ")



def truncate_string(s: str, max_length: int = 50) -> str:
    """Truncate string with ellipsis if too long.

    Args:
        s: String to truncate
        max_length: Maximum length

    Returns:
        Truncated string
    """
    if len(s) <= max_length:
        return s
    return s[: max_length - 3] + "..."


def format_number(n: float, decimals: int = 2) -> str:
    """Format number for display.

    Args:
        n: Number to format
        decimals: Number of decimal places

    Returns:
        Formatted number string
    """
    if isinstance(n, int):
        return f"{n:,}"
    return f"{n:,.{decimals}f}"
