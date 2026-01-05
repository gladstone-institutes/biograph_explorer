"""RAG chat interface for Streamlit UI.

Features:
- Chat-style Q&A interface with Claude
- Structured citations with validation
- Click citation → expand PyVis subgraph
- Conversation history
- Example questions
- Token usage tracking

Phase 3 Status: Stub only - not implemented in Phase 2
TODO: Implement in Phase 3 with Anthropic SDK
"""

from typing import List, Dict, Any, Optional
import streamlit as st
import networkx as nx


def render_rag_chat(
    graph: nx.DiGraph,
) -> None:
    """Render RAG chat interface.

    Args:
        graph: NetworkX knowledge graph

    Example (Phase 3):
        >>> render_rag_chat(graph)
        # User can ask questions, get answers with citations

    Phase 3: Not implemented
    """
    raise NotImplementedError("Phase 3: RAG chat not implemented")


def render_chat_history(messages: List[Dict[str, Any]]) -> None:
    """Render chat message history.

    Args:
        messages: List of message dicts (role, content, citations)

    Phase 3: Not implemented
    """
    raise NotImplementedError("Phase 3: Chat history rendering not implemented")


def render_example_questions() -> Optional[str]:
    """Render example question buttons.

    Returns:
        Selected example question or None

    Phase 3: Not implemented
    """
    EXAMPLE_QUESTIONS = [
        "Why is BACE1 a high-priority target?",
        "What pathways connect APOE and TREM2?",
        "Which nodes are most central in the amyloid processing cluster?",
        "What evidence links CD33 to neuroinflammation?",
    ]
    raise NotImplementedError("Phase 3: Example questions not implemented")


def render_citation_subgraph(
    citation: Any,
    graph: nx.DiGraph,
) -> None:
    """Render PyVis subgraph for a citation.

    Args:
        citation: Citation object with node_ids
        graph: Full NetworkX graph

    Phase 3: Not implemented
    """
    raise NotImplementedError("Phase 3: Citation subgraph not implemented")


def render_token_usage(input_tokens: int, output_tokens: int, cost: float) -> None:
    """Render token usage and cost information.

    Args:
        input_tokens: Number of input tokens
        output_tokens: Number of output tokens
        cost: Estimated cost in USD

    Phase 3: Not implemented
    """
    raise NotImplementedError("Phase 3: Token usage display not implemented")
