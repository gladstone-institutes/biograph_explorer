"""Data formatting utilities for display.

Handles:
- Node label formatting (truncate long names, add metadata)
- Edge label formatting (predicates → human-readable)
- Export formatting (HTML, JSON)

Phase 2 Status: Stub created
TODO: Implement formatters
"""

import re
from typing import Dict, List, Any, Optional
import networkx as nx

# Unicode emoji / pictograph blocks only. Deliberately EXCLUDES Greek and other scientific
# letters (alpha, beta, mu, ...) and ordinary punctuation so gene/chemical text is preserved.
_EMOJI_RE = re.compile(
    "["
    "\U0001F300-\U0001FAFF"  # symbols & pictographs, emoticons, supplemental, extended-A
    "\U0001F000-\U0001F0FF"  # mahjong / dominoes / playing cards
    "\U0001F1E6-\U0001F1FF"  # regional indicator (flags)
    "\U00002600-\U000027BF"  # misc symbols + dingbats (warning, check, cross, ...)
    "\U00002B00-\U00002BFF"  # misc symbols and arrows (stars, etc.)
    "\U0000FE00-\U0000FE0F"  # variation selectors (emoji presentation)
    "]+",
    flags=re.UNICODE,
)


def strip_emoji(text: str) -> str:
    """Remove emoji / pictographic symbols from text, preserving scientific letters and markdown.

    Scoped to emoji unicode blocks so Greek letters (used in gene/protein names) and normal
    punctuation survive. Collapses any whitespace left by a removed emoji.
    """
    if not text:
        return text
    cleaned = _EMOJI_RE.sub("", text)
    # Tidy up spaces left dangling by removed symbols (e.g. "drug X" stays single-spaced; a removed
    # trailing emoji does not leave a trailing space). Per-line rstrip preserves markdown structure.
    cleaned = re.sub(r"[ \t]{2,}", " ", cleaned)
    cleaned = re.sub(r" +([.,;:!?])", r"\1", cleaned)
    cleaned = "\n".join(line.rstrip() for line in cleaned.split("\n"))
    return cleaned


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
