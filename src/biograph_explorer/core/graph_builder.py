"""NetworkX graph builder from TRAPI responses.

Handles:
- Conversion from TRAPI edge format to NetworkX DiGraph
- Rich node attributes (labels, categories, is_query_gene flags)
- Edge attributes (predicates, qualifiers, knowledge sources)
- Gene frequency calculation (convergence metric)
- Subgraph extraction

"""

from typing import List, Dict, Any, Optional, Set
import networkx as nx
from pydantic import BaseModel, Field
import logging

try:
    from TCT import name_resolver
    TCT_AVAILABLE = True
except ImportError:
    TCT_AVAILABLE = False

logger = logging.getLogger(__name__)


class KnowledgeGraph(BaseModel):
    """Wrapper for NetworkX graph with metadata."""

    graph: Any = Field(description="NetworkX DiGraph (excluded from serialization)")
    num_nodes: int = Field(description="Number of nodes")
    num_edges: int = Field(description="Number of edges")
    query_genes: List[str] = Field(description="Input gene CURIEs")
    node_categories: Dict[str, int] = Field(default_factory=dict, description="Node category counts")
    metadata: Dict[str, Any] = Field(default_factory=dict, description="Graph metadata")

    class Config:
        arbitrary_types_allowed = True


class GraphBuilder:
    """Builds NetworkX graphs from TRAPI query results.

    Example:
        >>> builder = GraphBuilder()
        >>> edges = [{"subject": "NCBIGene:1803", "object": "MONDO:0100096", "predicate": "biolink:associated_with"}]
        >>> gene_curies = ["NCBIGene:1803"]
        >>> kg = builder.build_from_trapi_edges(edges, gene_curies)
        >>> print(f"Built graph: {kg.num_nodes} nodes, {kg.num_edges} edges")
    """

    def __init__(self):
        """Initialize graph builder."""
        if not TCT_AVAILABLE:
            raise ImportError("TCT library required for node name resolution")

    def build_from_trapi_edges(
        self,
        edges: List[Dict[str, Any]],
        query_gene_curies: List[str],
        curie_to_symbol: Optional[Dict[str, str]] = None,
        curie_to_name: Optional[Dict[str, str]] = None,
    ) -> KnowledgeGraph:
        """Build NetworkX DiGraph from TRAPI edge list.

        Extracted from notebook cells 12, 18.

        Args:
            edges: List of TRAPI edges (dicts with subject, object, predicate)
            query_gene_curies: List of input gene CURIEs for highlighting
            curie_to_symbol: Optional mapping of gene CURIEs to original symbols
            curie_to_name: Optional mapping of CURIEs to human-readable names.
                If provided, skips network lookup for names (uses cached names).

        Returns:
            KnowledgeGraph with NetworkX DiGraph and metadata
        """
        if curie_to_symbol is None:
            curie_to_symbol = {}
        if not edges:
            logger.warning("No edges provided - returning empty graph")
            return KnowledgeGraph(
                graph=nx.DiGraph(),
                num_nodes=0,
                num_edges=0,
                query_genes=query_gene_curies,
                node_categories={},
                metadata={},
            )

        logger.info(f"Building graph from {len(edges)} edges...")

        # Extract subject-predicate-object triples
        subjects = []
        objects = []
        predicates = []

        for edge in edges:
            if isinstance(edge, dict):
                subjects.append(edge.get("subject", ""))
                objects.append(edge.get("object", ""))
                predicates.append(edge.get("predicate", ""))

        # Get unique nodes for name lookup
        unique_nodes = list(set(subjects + objects))

        # Use cached names if provided, otherwise fetch from network
        if curie_to_name:
            curie_to_label = curie_to_name
            logger.info(f"Using {len(curie_to_label)} cached node names")
        else:
            logger.info(f"Looking up names for {len(unique_nodes)} unique nodes...")
            curie_to_label = self._lookup_node_names(unique_nodes)

        # Build NetworkX graph with CURIEs as node IDs
        graph = nx.DiGraph()

        for i, (subj, obj, pred) in enumerate(zip(subjects, objects, predicates)):
            if subj and obj:  # Skip empty entries
                # Preserve all edge attributes from TRAPI
                edge_attrs = edges[i]

                # Extract sources information (full detail)
                sources = edge_attrs.get('sources', [])

                # Extract attributes array (contains publications, sentences, confidence scores, etc.)
                attributes = edge_attrs.get('attributes', [])

                # Extract publications from all patterns (top-level + nested)
                publications = self._extract_publications_robust(attributes)

                # Extract supporting text/sentences from all patterns
                sentences = self._extract_supporting_text_robust(attributes)

                # Extract confidence scores from all patterns
                confidence_scores = self._extract_confidence_scores_robust(attributes)

                graph.add_edge(
                    subj,
                    obj,
                    predicate=pred,
                    sources=sources,  # Full source information
                    attributes=attributes,  # Complete attributes array
                    publications=publications,  # Extracted publications
                    sentences=sentences,  # Supporting sentences
                    confidence_scores=confidence_scores,  # Confidence scores
                    qualifiers=edge_attrs.get('qualifiers', []),
                    query_result_id=edge_attrs.get('query_result_id'),
                    # Keep legacy knowledge_source for backward compatibility
                    knowledge_source=edge_attrs.get('knowledge_source', []),
                )

        logger.info(f"Created graph: {graph.number_of_nodes()} nodes, {graph.number_of_edges()} edges")

        # Add node attributes
        self._add_node_attributes(graph, query_gene_curies, curie_to_label, curie_to_symbol)

        # Count categories
        node_categories = {}
        for node in graph.nodes():
            cat = graph.nodes[node].get("category", "Unknown")
            node_categories[cat] = node_categories.get(cat, 0) + 1

        logger.info(f"Node categories: {node_categories}")

        return KnowledgeGraph(
            graph=graph,
            num_nodes=graph.number_of_nodes(),
            num_edges=graph.number_of_edges(),
            query_genes=query_gene_curies,
            node_categories=node_categories,
            metadata={
                "density": nx.density(graph),
                "is_connected": nx.is_weakly_connected(graph),
            },
        )

    def _add_node_attributes(
        self,
        graph: nx.DiGraph,
        query_gene_curies: List[str],
        curie_to_label: Dict[str, str],
        curie_to_symbol: Dict[str, str],
    ) -> None:
        """Add rich attributes to graph nodes.

        Attributes:
            - label: Human-readable name
            - original_symbol: Original gene symbol for query genes
            - category: biolink category (Gene, Disease, Protein, etc.)
            - is_query_gene: Boolean flag for input genes
            - curie: Node identifier

        Extracted from notebook cell 18.

        Args:
            graph: NetworkX graph to annotate
            query_gene_curies: List of input gene CURIEs
            curie_to_label: Dictionary mapping CURIEs to labels
            curie_to_symbol: Dictionary mapping CURIEs to original gene symbols
        """
        for node in graph.nodes():
            # Add label
            graph.nodes[node]["label"] = curie_to_label.get(node, node)
            graph.nodes[node]["curie"] = node

            # Add original symbol for query genes
            if node in curie_to_symbol:
                graph.nodes[node]["original_symbol"] = curie_to_symbol[node]

            # Classify category based on CURIE prefix
            if node in query_gene_curies:
                graph.nodes[node]["category"] = "Gene"
                graph.nodes[node]["is_query_gene"] = True
            elif node.startswith("MONDO:") or node.startswith("DOID:"):
                graph.nodes[node]["category"] = "Disease"
                graph.nodes[node]["is_query_gene"] = False
            elif node.startswith("NCBIGene:") or node.startswith("HGNC:"):
                graph.nodes[node]["category"] = "Gene"
                graph.nodes[node]["is_query_gene"] = False
            elif node.startswith("UniProtKB:") or node.startswith("PR:"):
                graph.nodes[node]["category"] = "Protein"
                graph.nodes[node]["is_query_gene"] = False
            elif node.startswith("CHEBI:") or node.startswith("CHEMBL:"):
                graph.nodes[node]["category"] = "ChemicalEntity"
                graph.nodes[node]["is_query_gene"] = False
            elif node.startswith("GO:"):
                graph.nodes[node]["category"] = "BiologicalProcess"
                graph.nodes[node]["is_query_gene"] = False
            else:
                graph.nodes[node]["category"] = "Other"
                graph.nodes[node]["is_query_gene"] = False

    def calculate_gene_frequency(
        self,
        graph: nx.DiGraph,
        query_genes: List[str],
    ) -> Dict[str, int]:
        """Calculate gene frequency (convergence metric) for each node.

        Gene frequency = number of query genes that have a path to this node.
        High gene frequency indicates a convergent node.

        Args:
            graph: NetworkX DiGraph
            query_genes: List of query gene node IDs

        Returns:
            Dictionary mapping node ID to gene frequency
        """
        logger.info(f"Calculating gene frequency for {len(query_genes)} query genes...")

        gene_frequency = {node: 0 for node in graph.nodes()}

        # For each query gene, find all reachable nodes
        for gene in query_genes:
            if gene not in graph:
                logger.warning(f"Query gene {gene} not found in graph")
                continue

            # Find all nodes reachable from this gene (BFS/DFS)
            reachable = nx.descendants(graph, gene)
            reachable.add(gene)  # Include the gene itself

            # Increment frequency for all reachable nodes
            for node in reachable:
                gene_frequency[node] += 1

        # Also count reverse direction (nodes that can reach the gene)
        for gene in query_genes:
            if gene not in graph:
                continue

            # Find all nodes that can reach this gene
            ancestors = nx.ancestors(graph, gene)

            for node in ancestors:
                gene_frequency[node] += 1

        logger.info(f"Gene frequency calculated. Max frequency: {max(gene_frequency.values()) if gene_frequency else 0}")

        return gene_frequency

    def extract_subgraph(
        self,
        graph: nx.DiGraph,
        center_nodes: List[str],
        k_hops: int = 1,
    ) -> nx.DiGraph:
        """Extract k-hop subgraph around specified nodes.

        Args:
            graph: Source NetworkX graph
            center_nodes: Center node IDs
            k_hops: Number of hops to include (default: 1)

        Returns:
            Subgraph as NetworkX DiGraph

        TODO: Implement subgraph extraction for RAG citations (Phase 3)
        """
        raise NotImplementedError("TODO: Implement k-hop subgraph extraction")

    def _lookup_node_names(self, curies: List[str]) -> Dict[str, str]:
        """Look up human-readable names for CURIEs using TCT.

        Extracted from notebook cell 12.

        Args:
            curies: List of CURIEs

        Returns:
            Dictionary mapping CURIE to name
        """
        node_info_dict = name_resolver.batch_lookup(curies)

        curie_to_name = {}
        for curie in curies:
            info = node_info_dict.get(curie)
            if info and hasattr(info, "name") and info.name:
                curie_to_name[curie] = info.name
            else:
                # Use CURIE as fallback
                curie_to_name[curie] = curie

        return curie_to_name

    def _extract_publications_robust(self, attributes: List[Dict[str, Any]]) -> List[str]:
        """Extract publications from all TRAPI patterns (top-level + nested).

        Handles 3 patterns:
        1. Top-level with original_attribute_name="publications"
        2. Top-level with attribute_type_id="biolink:publications"
        3. Nested inside biolink:has_supporting_study_result

        Args:
            attributes: List of TRAPI attribute dictionaries

        Returns:
            Deduplicated list of publication IDs
        """
        pubs = []

        for attr in attributes:
            # Pattern 1 & 2: Top-level publications
            if attr.get('attribute_type_id') == 'biolink:publications':
                value = attr.get('value', [])
                # Handle both list and single string values
                if isinstance(value, list):
                    pubs.extend(value)
                elif value:
                    pubs.append(value)

            # Pattern 3: Nested in has_supporting_study_result
            if attr.get('attribute_type_id') == 'biolink:has_supporting_study_result':
                for nested in attr.get('attributes', []):
                    if nested.get('attribute_type_id') == 'biolink:publications':
                        value = nested.get('value', [])
                        if isinstance(value, list):
                            pubs.extend(value)
                        elif value:
                            pubs.append(value)

        # Deduplicate and filter out None/empty strings
        return list(set(filter(None, pubs)))

    def _extract_supporting_text_robust(self, attributes: List[Dict[str, Any]]) -> List[str]:
        """Extract supporting text/sentences from all TRAPI patterns.

        Handles 3 patterns:
        1. Legacy: original_attribute_name="sentences"
        2. Modern: attribute_type_id="biolink:supporting_text" (top-level)
        3. Nested: biolink:supporting_text inside has_supporting_study_result

        Args:
            attributes: List of TRAPI attribute dictionaries

        Returns:
            List of supporting text strings
        """
        texts = []

        for attr in attributes:
            # Pattern 1: Legacy sentences attribute
            if attr.get('original_attribute_name') == 'sentences':
                value = attr.get('value', '')
                if value:
                    texts.append(value)

            # Pattern 2: Modern top-level supporting_text
            if attr.get('attribute_type_id') == 'biolink:supporting_text':
                value = attr.get('value', '')
                if value:
                    texts.append(value)

            # Pattern 3: Nested in has_supporting_study_result
            if attr.get('attribute_type_id') == 'biolink:has_supporting_study_result':
                for nested in attr.get('attributes', []):
                    if nested.get('attribute_type_id') == 'biolink:supporting_text':
                        value = nested.get('value', '')
                        if value:
                            texts.append(value)

        # Filter out empty strings
        return [t for t in texts if t]

    def _extract_confidence_scores_robust(self, attributes: List[Dict[str, Any]]) -> Dict[str, float]:
        """Extract confidence scores from all TRAPI patterns.

        Handles multiple patterns:
        1. Legacy: original_attribute_name="tmkp_confidence_score"
        2. Modern top-level: attribute_type_id="biolink:extraction_confidence_score"
        3. STRING: original_attribute_name="Combined_score"
        4. Nested: extraction_confidence_score inside has_supporting_study_result

        Args:
            attributes: List of TRAPI attribute dictionaries

        Returns:
            Dictionary mapping score type to value
        """
        scores = {}

        for attr in attributes:
            orig_name = attr.get('original_attribute_name', '')
            attr_type = attr.get('attribute_type_id', '')
            value = attr.get('value')

            # Pattern 1: tmkp_confidence_score (legacy text mining)
            if orig_name == 'tmkp_confidence_score' and value is not None:
                try:
                    scores['tmkp_confidence_score'] = float(value)
                except (ValueError, TypeError):
                    pass

            # Pattern 2: extraction_confidence_score (modern)
            if attr_type == 'biolink:extraction_confidence_score' and value is not None:
                try:
                    scores['extraction_confidence_score'] = float(value)
                except (ValueError, TypeError):
                    pass

            # Pattern 3: STRING Combined_score
            if orig_name == 'Combined_score' and value is not None:
                try:
                    scores['combined_score'] = float(value)
                except (ValueError, TypeError):
                    pass

            # Pattern 4: Nested in has_supporting_study_result
            if attr_type == 'biolink:has_supporting_study_result':
                for nested in attr.get('attributes', []):
                    if nested.get('attribute_type_id') == 'biolink:extraction_confidence_score':
                        nested_value = nested.get('value')
                        if nested_value is not None:
                            try:
                                # Use different key for nested scores
                                scores['extraction_confidence_score_nested'] = float(nested_value)
                            except (ValueError, TypeError):
                                pass

        return scores
