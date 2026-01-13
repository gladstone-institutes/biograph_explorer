"""Core modules for GeneSet Translator.

Contains the main business logic:
- TRAPI client for querying Translator APIs
- Graph builder for NetworkX construction
- Node annotator for enriching nodes with metadata
- RAG system for LLM-assisted exploration (Phase 3)

Phase 2 Status: Stubs created, implementation in progress
"""

from .trapi_client import TRAPIClient, TRAPIResponse
from .graph_builder import GraphBuilder, KnowledgeGraph
from .node_annotator import NodeAnnotator
from .translator_node import TranslatorNode, TranslatorAttribute, TranslatorEdge
from .hpa_client import HPAClient, HPAAnnotation

# Phase 3 (stub)
# from .rag_system import RAGSystem

__all__ = [
    "TRAPIClient",
    "TRAPIResponse",
    "GraphBuilder",
    "KnowledgeGraph",
    "NodeAnnotator",
    "TranslatorNode",
    "TranslatorAttribute",
    "TranslatorEdge",
    "HPAClient",
    "HPAAnnotation",
]
