"""RAG system for LLM-assisted graph exploration using Claude.

Implements:
- 3-layer context strategy (graph stats → relevant subgraph → node details)
- Claude Haiku 4 integration with tool use for structured citations
- Citation validation against actual graph
- Subgraph extraction for visualization


TODO: Implement in Phase 3 with Anthropic SDK
"""

from typing import List, Optional
import networkx as nx
from pydantic import BaseModel, Field


class Citation(BaseModel):
    """Structured citation from LLM response."""

    node_ids: List[str] = Field(description="Cited node IDs")
    metric_name: Optional[str] = Field(default=None, description="Metric referenced")
    metric_value: Optional[float] = Field(default=None, description="Metric value")
    reasoning: str = Field(description="Why this node is relevant")


class RAGResponse(BaseModel):
    """Response from RAG system."""

    answer: str = Field(description="Natural language answer")
    citations: List[Citation] = Field(default_factory=list, description="Structured citations")
    context_tokens: int = Field(description="Tokens used for context")
    subgraph: Optional[nx.DiGraph] = Field(default=None, description="Cited subgraph for visualization")

    class Config:
        arbitrary_types_allowed = True


class RAGSystem:
    """RAG system for graph exploration with Claude.

    Example (Phase 3):
        >>> rag = RAGSystem(api_key="sk-...", model="claude-haiku-4")
        >>> response = rag.ask_question(
        ...     question="Why is BACE1 a high-priority target?",
        ...     graph=graph
        ... )
        >>> print(response.answer)
        >>> # Visualize cited subgraph
        >>> viz.render(response.subgraph)
    """

    def __init__(
        self,
        api_key: str,
        model: str = "claude-haiku-4",
        max_tokens: int = 1000,
    ):
        """Initialize RAG system.

        Args:
            api_key: Anthropic API key
            model: Claude model to use
            max_tokens: Max tokens per response
        """
        raise NotImplementedError("Phase 3: Not implemented yet")

    def ask_question(
        self,
        question: str,
        graph: nx.DiGraph,
    ) -> RAGResponse:
        """Answer question about the knowledge graph.

        Args:
            question: User question
            graph: Full NetworkX knowledge graph

        Returns:
            RAGResponse with answer, citations, and subgraph

        Phase 3: Not implemented
        """
        raise NotImplementedError("Phase 3: RAG system not implemented")

    def _build_context(
        self,
        graph: nx.DiGraph,
        question: str,
    ) -> str:
        """Build context for LLM.

        Phase 3: Not implemented
        """
        raise NotImplementedError("Phase 3: Context building not implemented")

    def _validate_citations(
        self,
        citations: List[Citation],
        graph: nx.DiGraph,
    ) -> List[Citation]:
        """Validate citations against actual graph data.

        Phase 3: Not implemented
        """
        raise NotImplementedError("Phase 3: Citation validation not implemented")
