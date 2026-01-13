"""LLM-assisted category summaries with citation graphs using Claude Haiku 4.

This module generates verifiable category-specific summaries with XML-based
citation extraction. LLM identifies relevant nodes by CURIE, and the system
extracts all edge/publication/sentence data directly from the graph.
Uses token-aware sampling to optimize costs while maintaining analytical depth.
"""

import hashlib
import json
import logging
import re
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple

import networkx as nx
from pydantic import BaseModel, Field

from ..utils.biolink_predicates import get_predicate_info

logger = logging.getLogger(__name__)


class CitationGraph(BaseModel):
    """Citation graph linking claims to nodes/edges/publications."""

    citation_id: int = Field(description="Unique citation ID")
    claim: str = Field(description="The specific claim being cited")
    node_ids: List[str] = Field(default_factory=list, description="Node CURIEs from LLM")
    edge_ids: List[str] = Field(default_factory=list, description="Edge IDs extracted from graph")
    context_node_ids: List[str] = Field(default_factory=list, description="Additional context nodes")
    publication_ids: List[str] = Field(default_factory=list, description="PMIDs extracted from graph")
    sentences: List[str] = Field(default_factory=list, description="Supporting text extracted from graph")
    confidence: str = Field(default="medium", description="Citation confidence: low/medium/high")


class SummaryData(BaseModel):
    """Complete summary with text and citation graphs."""

    category: str = Field(description="Node category (e.g., Protein, BiologicalProcess)")
    summary_text: str = Field(description="Summary text with [Citation N] markers")
    citations: List[CitationGraph] = Field(default_factory=list, description="Citation graphs")
    metadata: Dict[str, Any] = Field(default_factory=dict, description="Tokens, sampling strategy, etc.")


class StagedCategoryQuery(BaseModel):
    """Pre-computed context with ACTUAL token counts for a category summary.

    Used for accurate cost estimation before LLM execution.
    Input token count is from Anthropic's official token counting API.
    """

    category: str = Field(description="Node category")
    context: Dict[str, Any] = Field(description="Prepared JSON context for LLM")
    input_tokens: int = Field(description="Total input tokens from Anthropic token counting API")
    system_prompt_tokens: int = Field(description="System prompt tokens (estimated for display breakdown)")
    output_tokens_estimate: int = Field(default=2000, description="Estimated output tokens")
    nodes_total: int = Field(description="Total nodes in this category")
    nodes_sampled: int = Field(description="Nodes included after sampling")
    edges_sampled: int = Field(description="Edges in sampled subgraph")
    is_cached: bool = Field(default=False, description="Whether cached result exists")

    @property
    def total_input_tokens(self) -> int:
        """Total input tokens (already includes system prompt from API)."""
        return self.input_tokens

    @property
    def context_tokens(self) -> int:
        """Estimated context tokens (total minus system prompt estimate)."""
        return max(0, self.input_tokens - self.system_prompt_tokens)

    @property
    def estimated_cost(self) -> float:
        """Calculate estimated cost using Haiku 4.5 pricing.

        Pricing: $1.00/M input tokens, $5.00/M output tokens
        """
        input_cost = (self.input_tokens / 1_000_000) * 1.00
        output_cost = (self.output_tokens_estimate / 1_000_000) * 5.00
        return input_cost + output_cost


class LLMSummarizer:
    """Generate citation-based category summaries using Claude Haiku 4."""

    # Default prompt file location (relative to this module)
    DEFAULT_PROMPT_PATH = Path(__file__).parent.parent / "config" / "prompts" / "category_summary.txt"

    def __init__(
        self,
        model: str = "claude-haiku-4-5",
        cache_dir: Path = Path("data/cache"),
        prompt_path: Optional[Path] = None
    ):
        """Initialize LLM summarizer.

        Args:
            model: Claude model to use
            cache_dir: Directory for caching summaries and citations
            prompt_path: Optional path to custom prompt file. If not provided,
                        uses the default prompt at config/prompts/category_summary.txt.
                        The prompt file should contain {category} placeholders.

        Note:
            Requires ANTHROPIC_API_KEY environment variable to be set.
            The Anthropic SDK will automatically read it from the environment.
        """
        try:
            from anthropic import Anthropic
            from dotenv import load_dotenv

            load_dotenv()
            # Anthropic SDK automatically reads from ANTHROPIC_API_KEY env var
            self.client = Anthropic()
            logger.info(f"Anthropic client initialized successfully")
        except ImportError:
            raise ImportError("anthropic package required. Install with: poetry add anthropic")

        self.model = model
        self.prompt_path = prompt_path or self.DEFAULT_PROMPT_PATH
        self.summary_cache_dir = cache_dir / "summaries"
        self.citation_cache_dir = cache_dir / "citation_graphs"
        self.summary_cache_dir.mkdir(parents=True, exist_ok=True)
        self.citation_cache_dir.mkdir(parents=True, exist_ok=True)

    def generate_category_summary(
        self,
        graph: nx.MultiDiGraph,
        category: str,
        query_genes: List[str],
        disease_curie: str,
        infores_metadata: Optional[Dict] = None,
        max_nodes: int = 20,
        min_gene_frequency: int = 2
    ) -> SummaryData:
        """Generate citation-based summary for a specific category using JSON context.

        Args:
            graph: Full knowledge graph
            category: Node category to summarize (e.g., "Protein", "BiologicalProcess")
            query_genes: List of input gene CURIEs
            disease_curie: Target disease CURIE
            infores_metadata: Optional knowledge source metadata
            max_nodes: Maximum intermediate nodes to include (user-adjustable)
            min_gene_frequency: Minimum gene_frequency for intermediate nodes (user-adjustable)

        Returns:
            SummaryData with summary text and citations
        """
        # Check cache first
        cache_key = self._generate_cache_key(query_genes, disease_curie, category, graph, min_gene_frequency)
        cached_summary = self._load_from_cache(cache_key, category)
        if cached_summary:
            logger.info(f"Using cached summary for {category}")
            cached_summary.metadata['from_cache'] = True
            return cached_summary

        # Filter to category nodes
        category_nodes = [
            node for node, data in graph.nodes(data=True)
            if data.get('category') == category
        ]

        if not category_nodes:
            logger.warning(f"No nodes found for category {category}")
            return SummaryData(
                category=category,
                summary_text=f"No {category} nodes found in the knowledge graph.",
                citations=[],
                metadata={'error': 'no_nodes'}
            )

        logger.info(f"Generating summary for {category}: {len(category_nodes)} nodes")

        # Sample top intermediates by query gene connectivity (matches visualization approach)
        sampled_subgraph = self._sample_nodes_by_query_gene_connections(
            graph,
            category_nodes,
            query_genes,
            max_nodes,
            disease_curie=disease_curie,
            min_gene_frequency=min_gene_frequency
        )

        # Prepare JSON context
        context = self._prepare_json_context(
            sampled_subgraph,
            category_nodes,
            query_genes,
            disease_curie,
            category,
            infores_metadata
        )

        # Log token usage (use Anthropic API for accurate count)
        system_prompt = self._build_system_prompt(category)
        context_json_str = json.dumps(context, indent=2)
        token_count = self._count_tokens_via_api(system_prompt, context_json_str)
        logger.info(f"JSON context for {category}: {token_count} tokens, "
                   f"{len(context['nodes'])} nodes, {len(context['edges'])} edges")

        # Generate with citations
        summary_text, citations = self._generate_with_citations(context, category, graph)

        # Create summary data
        summary_data = SummaryData(
            category=category,
            summary_text=summary_text,
            citations=citations,
            metadata={
                'timestamp': datetime.now().isoformat(),
                'nodes_total': len(category_nodes),
                'nodes_sampled': sampled_subgraph.number_of_nodes(),
                'edges_sampled': sampled_subgraph.number_of_edges(),
                'max_nodes': max_nodes,
                'sampling_strategy': 'query_gene_connectivity',
                'token_count': token_count,
                'model': self.model,
                'format_version': 'json_v3',
                'from_cache': False
            }
        )

        # Cache the result
        self._save_to_cache(cache_key, category, summary_data)

        return summary_data

    def _generate_cache_key(
        self,
        query_genes: List[str],
        disease_curie: str,
        category: str,
        graph: nx.MultiDiGraph,
        min_gene_frequency: int = 2
    ) -> str:
        """Generate hash-based cache key including model, format version, and filter settings."""
        format_version = "json_v3"  # Increment when changing JSON structure or filters
        key_str = f"{','.join(sorted(query_genes))}|{disease_curie}|{category}|{graph.number_of_edges()}|{self.model}|{format_version}|minfreq{min_gene_frequency}"
        return hashlib.md5(key_str.encode()).hexdigest()

    def _load_from_cache(self, cache_key: str, category: str) -> Optional[SummaryData]:
        """Load summary from cache if available and fresh."""
        summary_file = self.summary_cache_dir / cache_key / f"summary_{category}.json"

        if summary_file.exists():
            cache_age = datetime.now() - datetime.fromtimestamp(summary_file.stat().st_mtime)
            if cache_age < timedelta(days=30):
                with open(summary_file, 'r') as f:
                    data = json.load(f)
                    return SummaryData(**data)

        return None

    def _save_to_cache(self, cache_key: str, category: str, summary_data: SummaryData):
        """Save summary to cache."""
        cache_subdir = self.summary_cache_dir / cache_key
        cache_subdir.mkdir(parents=True, exist_ok=True)

        summary_file = cache_subdir / f"summary_{category}.json"
        with open(summary_file, 'w') as f:
            json.dump(summary_data.model_dump(), f, indent=2)

        # Also save citations separately for citation viewer
        citation_subdir = self.citation_cache_dir / cache_key
        citation_subdir.mkdir(parents=True, exist_ok=True)

        citation_file = citation_subdir / f"citations_{category}.json"
        with open(citation_file, 'w') as f:
            json.dump([c.model_dump() for c in summary_data.citations], f, indent=2)

        logger.info(f"Cached summary and citations for {category}")

    def _count_tokens_via_api(self, system_prompt: str, user_content: str) -> int:
        """Count tokens using Anthropic's official token counting API.

        This provides EXACT token counts matching what the API will charge.
        Uses the messages.count_tokens endpoint.

        Args:
            system_prompt: System prompt text
            user_content: User message content (JSON context)

        Returns:
            Exact input token count from Anthropic API
        """
        try:
            response = self.client.messages.count_tokens(
                model=self.model,
                system=system_prompt,
                messages=[{
                    "role": "user",
                    "content": user_content
                }]
            )
            return response.input_tokens
        except Exception as e:
            logger.warning(f"Token counting API failed, falling back to character estimate: {e}")
            # Fallback: rough estimate of 4 chars per token
            return len(system_prompt + user_content) // 4

    def _estimate_system_prompt_tokens(self, system_prompt: str) -> int:
        """Estimate system prompt tokens for display breakdown.

        Uses the Anthropic API with an empty user message to get
        an accurate system prompt token count.

        Args:
            system_prompt: System prompt text

        Returns:
            Estimated system prompt token count
        """
        try:
            response = self.client.messages.count_tokens(
                model=self.model,
                system=system_prompt,
                messages=[{
                    "role": "user",
                    "content": ""
                }]
            )
            return response.input_tokens
        except Exception as e:
            logger.warning(f"System prompt token counting failed: {e}")
            # Fallback: rough estimate of 4 chars per token
            return len(system_prompt) // 4

    def stage_category_query(
        self,
        graph: nx.MultiDiGraph,
        category: str,
        query_genes: List[str],
        disease_curie: str,
        max_nodes: int = 20,
        min_gene_frequency: int = 2,
        infores_metadata: Optional[Dict] = None
    ) -> StagedCategoryQuery:
        """Stage a category query by preparing context and counting actual tokens.

        This prepares everything needed for an LLM call WITHOUT making the call,
        allowing accurate cost estimation before execution.

        Args:
            graph: Full knowledge graph
            category: Node category to summarize
            query_genes: List of input gene CURIEs
            disease_curie: Target disease CURIE
            max_nodes: Maximum category nodes to include (user-adjustable)
            min_gene_frequency: Minimum gene_frequency for intermediate nodes (user-adjustable)
            infores_metadata: Optional knowledge source metadata

        Returns:
            StagedCategoryQuery with context, actual token counts, and cost estimate
        """
        # Check if cached result exists
        cache_key = self._generate_cache_key(query_genes, disease_curie, category, graph, min_gene_frequency)
        cached_summary = self._load_from_cache(cache_key, category)
        is_cached = cached_summary is not None

        # Filter to category nodes
        category_nodes = [
            node for node, data in graph.nodes(data=True)
            if data.get('category') == category
        ]

        nodes_total = len(category_nodes)

        if not category_nodes:
            # No nodes for this category
            return StagedCategoryQuery(
                category=category,
                context={},
                input_tokens=0,
                system_prompt_tokens=0,
                output_tokens_estimate=0,
                nodes_total=0,
                nodes_sampled=0,
                edges_sampled=0,
                is_cached=is_cached
            )

        # Sample nodes by query gene connectivity
        sampled_subgraph = self._sample_nodes_by_query_gene_connections(
            graph,
            category_nodes,
            query_genes,
            max_nodes,
            disease_curie=disease_curie,
            min_gene_frequency=min_gene_frequency
        )

        # Prepare JSON context
        context = self._prepare_json_context(
            sampled_subgraph,
            category_nodes,
            query_genes,
            disease_curie,
            category,
            infores_metadata
        )

        # Build system prompt
        system_prompt = self._build_system_prompt(category)
        context_json_str = json.dumps(context, indent=2)

        # Count tokens using Anthropic's official token counting API
        # This gives exact token counts matching what the API will charge
        input_tokens = self._count_tokens_via_api(system_prompt, context_json_str)

        # For display purposes, estimate the breakdown (API gives total only)
        # Use Anthropic API with empty user message to get system prompt tokens
        system_prompt_tokens = self._estimate_system_prompt_tokens(system_prompt)

        # Estimate output tokens based on context complexity
        # Based on observed data:
        # - Small contexts (<10 nodes): ~500-800 output tokens
        # - Large contexts (>20 nodes): ~2000+ output tokens
        # Formula: base 500 + 75 tokens per node (up to max 2500)
        nodes_count = sampled_subgraph.number_of_nodes()
        output_tokens_estimate = min(500 + (nodes_count * 75), 2500)

        logger.info(
            f"Staged {category}: {input_tokens:,} input tokens (API), "
            f"{nodes_count} nodes, {sampled_subgraph.number_of_edges()} edges, "
            f"~{output_tokens_estimate:,} output tokens (est)"
            f"{' (cached)' if is_cached else ''}"
        )

        return StagedCategoryQuery(
            category=category,
            context=context,
            input_tokens=input_tokens,
            system_prompt_tokens=system_prompt_tokens,
            output_tokens_estimate=output_tokens_estimate,
            nodes_total=nodes_total,
            nodes_sampled=nodes_count,
            edges_sampled=sampled_subgraph.number_of_edges(),
            is_cached=is_cached
        )

    def stage_all_categories(
        self,
        graph: nx.MultiDiGraph,
        categories: List[str],
        query_genes: List[str],
        disease_curie: str,
        max_nodes: int = 20,
        min_gene_frequency: int = 2,
        infores_metadata: Optional[Dict] = None
    ) -> Tuple[List[StagedCategoryQuery], float]:
        """Stage all selected categories and compute total cost.

        Args:
            graph: Full knowledge graph
            categories: List of categories to stage
            query_genes: List of input gene CURIEs
            disease_curie: Target disease CURIE
            max_nodes: Maximum category nodes to include
            min_gene_frequency: Minimum gene_frequency for intermediate nodes
            infores_metadata: Optional knowledge source metadata

        Returns:
            Tuple of (list of staged queries, total estimated cost)
        """
        staged_queries = []
        total_cost = 0.0

        for category in categories:
            staged = self.stage_category_query(
                graph=graph,
                category=category,
                query_genes=query_genes,
                disease_curie=disease_curie,
                max_nodes=max_nodes,
                min_gene_frequency=min_gene_frequency,
                infores_metadata=infores_metadata
            )
            staged_queries.append(staged)

            # Only count cost for non-cached queries
            if not staged.is_cached:
                total_cost += staged.estimated_cost

        logger.info(
            f"Staged {len(categories)} categories: total cost ${total_cost:.4f} "
            f"({sum(1 for sq in staged_queries if sq.is_cached)} cached)"
        )

        return staged_queries, total_cost

    def _sample_nodes_by_query_gene_connections(
        self,
        graph: nx.MultiDiGraph,
        category_nodes: List[str],
        query_genes: List[str],
        max_nodes: int,
        disease_curie: Optional[str] = None,
        min_gene_frequency: int = 2
    ) -> nx.MultiDiGraph:
        """Sample category nodes by query gene connectivity (matches visualization approach).

        This sampling strategy mirrors the "Top Intermediates" approach used in the
        Network tab visualization, ensuring consistency between what users see and
        what the LLM summarizes.

        Strategy:
            1. Score each category node by number of direct edges to/from query genes
            2. Filter nodes by min_gene_frequency threshold
            3. Take top-N category nodes by score
            4. Always include query genes and disease node
            5. Return subgraph with selected nodes and ALL their connecting edges

        Args:
            graph: Full knowledge graph
            category_nodes: Nodes in the target category
            query_genes: Query gene CURIEs
            max_nodes: Maximum number of category nodes to include
            disease_curie: Target disease CURIE (always included if present)
            min_gene_frequency: Minimum gene_frequency to include in LLM context (user-adjustable)

        Returns:
            Subgraph with top category nodes by query gene connectivity
        """
        query_genes_set = set(query_genes)

        # Score each category node by query gene connections
        # Filter: only include nodes with gene_frequency >= min_gene_frequency (convergent nodes)
        node_scores = {}
        filtered_count = 0
        for node in category_nodes:
            if node not in graph.nodes():
                continue
            node_data = graph.nodes[node]
            gene_freq = node_data.get('gene_frequency', 0)

            # Skip nodes with gene_frequency < 2 (not convergent)
            if gene_freq < min_gene_frequency:
                filtered_count += 1
                continue

            score = 0
            for gene in query_genes_set:
                if gene in graph.nodes():
                    # Count edges in both directions
                    if graph.has_edge(gene, node):
                        score += 1
                    if graph.has_edge(node, gene):
                        score += 1
            node_scores[node] = score

        logger.info(f"Filtered {filtered_count} nodes with gene_frequency < {min_gene_frequency}")

        # Sort by score (descending) and take top-N
        top_nodes = sorted(node_scores.items(), key=lambda x: x[1], reverse=True)[:max_nodes]

        # Build node set: top category nodes + query genes (exempt from filter) + disease
        nodes_to_include = set(node for node, _ in top_nodes)
        # Query genes are always included regardless of gene_frequency
        nodes_to_include.update(g for g in query_genes if g in graph.nodes())

        # Always include disease node
        if disease_curie and disease_curie in graph.nodes():
            nodes_to_include.add(disease_curie)

        # Create subgraph with all edges between selected nodes
        subgraph = graph.subgraph(nodes_to_include).copy()

        # Log sampling statistics
        top_score = top_nodes[0][1] if top_nodes else 0
        min_score = top_nodes[-1][1] if top_nodes else 0
        logger.info(
            f"Sampled {len(top_nodes)} category nodes from {len(node_scores)} total "
            f"(max: {max_nodes}, score range: {min_score}-{top_score})"
        )
        logger.info(
            f"Subgraph: {subgraph.number_of_nodes()} nodes, {subgraph.number_of_edges()} edges"
        )

        return subgraph


    def _prepare_json_context(
        self,
        subgraph: nx.MultiDiGraph,
        all_category_nodes: List[str],
        query_genes: List[str],
        disease_curie: str,
        category: str,
        infores_metadata: Optional[Dict]
    ) -> Dict[str, Any]:
        """Prepare structured JSON context with complete edge information.

        Args:
            subgraph: Sampled subgraph with importance-scored edges
            all_category_nodes: All nodes in category (for stats)
            query_genes: Input gene CURIEs
            disease_curie: Target disease
            category: The category being summarized (e.g., "BiologicalProcess")
            infores_metadata: Knowledge source metadata

        Returns:
            Structured JSON context dictionary
        """
        # Query context
        query_context = {
            "query_genes": query_genes,
            "disease_curie": disease_curie,
            "category": category,  # Use the category parameter directly
            "total_category_nodes": len(all_category_nodes),
            "sampled_nodes": subgraph.number_of_nodes()
        }

        # Nodes
        nodes = []
        for node_id, data in subgraph.nodes(data=True):
            node_obj = {
                "curie": node_id,
                "label": data.get('label', node_id),
                "category": data.get('category', 'Unknown'),
                "gene_frequency": data.get('gene_frequency', 0)
            }

            # Mark disease node explicitly
            if node_id == disease_curie:
                node_obj['is_disease'] = True

            # Optional attributes
            if 'is_query_gene' in data:
                node_obj['is_query_gene'] = data['is_query_gene']
            if 'pagerank' in data:
                node_obj['pagerank'] = round(data['pagerank'], 4)
            if 'betweenness' in data:
                node_obj['betweenness'] = round(data['betweenness'], 4)

            # HPA annotation fields (skip GO annotations per user request)
            annotation_features = data.get('annotation_features', {})
            if annotation_features:
                # Tissue specificity
                if 'hpa_tissue_specificity' in annotation_features:
                    node_obj['tissue_specificity'] = annotation_features['hpa_tissue_specificity']
                if 'hpa_top_tissues' in annotation_features:
                    node_obj['top_tissues'] = annotation_features['hpa_top_tissues']

                # Cell type specificity
                if 'hpa_cell_type_specificity' in annotation_features:
                    node_obj['cell_type_specificity'] = annotation_features['hpa_cell_type_specificity']
                if 'hpa_top_cell_types' in annotation_features:
                    node_obj['top_cell_types'] = annotation_features['hpa_top_cell_types']

                # Immune cell specificity
                if 'hpa_immune_cell_specificity' in annotation_features:
                    node_obj['immune_cell_specificity'] = annotation_features['hpa_immune_cell_specificity']
                if 'hpa_top_immune_cells' in annotation_features:
                    node_obj['top_immune_cells'] = annotation_features['hpa_top_immune_cells']

                # Disease involvement and protein class
                if 'hpa_disease_involvement' in annotation_features:
                    node_obj['disease_involvement'] = annotation_features['hpa_disease_involvement']
                if 'hpa_protein_class' in annotation_features:
                    node_obj['protein_class'] = annotation_features['hpa_protein_class']

            nodes.append(node_obj)

        # Edges
        edges = []
        for u, v, key, data in subgraph.edges(keys=True, data=True):
            # Format: subject→predicate→object (semantic triple order)
            predicate_for_id = key.rsplit('_', 1)[0] if '_' in key else key
            edge_id = f"{u}→{predicate_for_id}→{v}"

            # Extract sources
            sources = data.get('sources', [])
            source_list = []
            for s in sources:
                if isinstance(s, dict):
                    source_list.append({
                        "resource_id": s.get('resource_id', ''),
                        "resource_role": s.get('resource_role', '')
                    })
                else:
                    source_list.append({"resource_id": str(s)})

            # Extract publications (no truncation, but limit display to 20 for readability)
            publications = data.get('publications', [])
            pub_list = publications[:20] if len(publications) > 20 else publications

            # Extract supporting text
            supporting_text = data.get('sentences', [])
            text_list = supporting_text[:5] if len(supporting_text) > 5 else supporting_text

            # Get predicate definition from Biolink model
            predicate = data.get('predicate', 'related_to')
            predicate_name = predicate.replace('biolink:', '') if predicate else 'related_to'
            predicate_info = get_predicate_info(predicate_name)
            predicate_definition = ""
            predicate_inverse = ""
            if predicate_info:
                predicate_definition = predicate_info.get('description', '')
                predicate_inverse = predicate_info.get('inverse', '')

            edge_obj = {
                "edge_id": edge_id,
                "subject": u,
                "object": v,
                "predicate": predicate,
                "predicate_definition": predicate_definition,
                "predicate_inverse": predicate_inverse,
                "sources": source_list,
                "publications": pub_list,
                "supporting_text": text_list,
                "confidence_scores": data.get('confidence_scores', {}),
                "importance_score": round(data.get('importance_score', 0.0), 2),
                "publication_count": data.get('publication_count', len(publications)),
                "connects_to_query_genes": data.get('connects_to_query_genes', [])
            }

            edges.append(edge_obj)

        # Knowledge sources summary
        knowledge_sources = {}
        if infores_metadata:
            summary = infores_metadata.get('summary', {})
            knowledge_sources = {
                "total_sources": summary.get('total_sources', 0),
                "top_sources": [
                    {
                        "name": src['name'],
                        "resource_id": src.get('id', ''),
                        "edge_count": src['edge_count']
                    }
                    for src in summary.get('top_sources', [])[:5]
                ]
            }

        return {
            "query_context": query_context,
            "nodes": nodes,
            "edges": edges,
            "knowledge_sources": knowledge_sources
        }

    def _generate_with_citations(
        self,
        context: Dict[str, Any],
        category: str,
        graph: nx.MultiDiGraph
    ) -> Tuple[str, List[CitationGraph]]:
        """Generate summary with XML-based citation extraction.

        The LLM outputs citations with node CURIEs in XML format, and the system
        extracts all edge/publication/sentence data directly from the graph.

        Args:
            context: Prepared JSON context dictionary
            category: Node category
            graph: Full graph for citation extraction

        Returns:
            Tuple of (summary_text, citations)
        """
        # Build the XML-structured prompt
        system_prompt = self._build_system_prompt(category)

        try:
            # Convert context to JSON string
            context_json_str = json.dumps(context, indent=2)

            # Single API call with XML-structured prompt
            response = self.client.messages.create(
                model=self.model,
                max_tokens=8000,  # Increased to ensure room for citations + summary
                temperature=0.1,  # Low temperature for consistent output
                system=system_prompt,
                messages=[{
                    "role": "user",
                    "content": context_json_str
                }]
            )

            # Extract response text
            response_text = ""
            for block in response.content:
                if block.type == "text":
                    response_text += block.text

            # Parse XML citations and extract graph data
            # Pass query context for complete pathway extraction
            query_genes = context.get('query_context', {}).get('query_genes', [])
            disease_curie = context.get('query_context', {}).get('disease_curie')
            citations = self._parse_xml_citations(
                response_text, graph, query_genes=query_genes, disease_curie=disease_curie
            )

            # Extract summary from <summary> tags (with robust fallbacks)
            summary_text = self._extract_summary(response_text, category=category)

            if not summary_text:
                logger.warning(f"No summary found in response for {category}")
                summary_text = f"Summary extraction failed for {category}. Check data/cache/llm_debug/ for raw response."

            # Log actual token usage from API response
            if hasattr(response, 'usage'):
                actual_input = response.usage.input_tokens
                actual_output = response.usage.output_tokens
                actual_cost = (actual_input / 1_000_000) * 1.00 + (actual_output / 1_000_000) * 5.00
                logger.info(
                    f"Generated summary for {category}: {len(citations)} citations, {len(summary_text)} chars | "
                    f"Actual usage: {actual_input:,} input + {actual_output:,} output = ${actual_cost:.4f}"
                )
            else:
                logger.info(f"Generated summary for {category}: {len(citations)} citations, {len(summary_text)} chars")

            return summary_text, citations

        except Exception as e:
            logger.error(f"LLM generation failed for {category}: {e}")
            return f"Failed to generate summary for {category}: {e}", []

    def _build_system_prompt(self, category: str) -> str:
        """Build the system prompt for category summaries from external file.

        Loads the prompt template from the configured prompt_path and substitutes
        {category} placeholders with the actual category name.

        Args:
            category: The node category being summarized (e.g., "Protein", "ChemicalEntity")

        Returns:
            Formatted system prompt string

        Raises:
            FileNotFoundError: If the prompt file doesn't exist
        """
        try:
            with open(self.prompt_path, 'r') as f:
                prompt_template = f.read()

            # Substitute {category} placeholders
            return prompt_template.format(category=category)

        except FileNotFoundError:
            logger.error(f"Prompt file not found: {self.prompt_path}")
            raise FileNotFoundError(
                f"Category summary prompt file not found at {self.prompt_path}. "
                f"Expected location: src/geneset_translator/config/prompts/category_summary.txt"
            )

    def _flatten_to_strings(self, data: Any) -> List[str]:
        """Flatten nested lists/structures and extract only string elements.

        Handles cases where LLM returns nested arrays like [["CURIE1", "CURIE2"]]
        or mixed types like ["CURIE1", 123, ["nested"]].

        Args:
            data: Parsed JSON data (could be list, string, or nested structure)

        Returns:
            Flat list of strings
        """
        result = []
        if isinstance(data, str):
            result.append(data)
        elif isinstance(data, list):
            for item in data:
                result.extend(self._flatten_to_strings(item))
        # Ignore non-string, non-list items (ints, dicts, None, etc.)
        return result

    def _parse_xml_citations(
        self,
        response_text: str,
        graph: nx.MultiDiGraph,
        query_genes: Optional[List[str]] = None,
        disease_curie: Optional[str] = None
    ) -> List[CitationGraph]:
        """Parse XML citations and extract full graph data for cited nodes.

        Args:
            response_text: LLM response containing XML citation blocks
            graph: Full graph for extracting edge/publication data
            query_genes: Query gene CURIEs for complete pathway extraction
            disease_curie: Disease CURIE for complete pathway extraction

        Returns:
            List of CitationGraph objects with extracted graph data
        """
        citations = []
        pattern = r'<save_citation_graph>(.*?)</save_citation_graph>'

        for match in re.finditer(pattern, response_text, re.DOTALL):
            block = match.group(1)

            # Extract citation_id
            citation_id_match = re.search(r'<citation_id>(\d+)</citation_id>', block)
            if not citation_id_match:
                logger.warning("Citation block missing citation_id, skipping")
                continue
            citation_id = int(citation_id_match.group(1))

            # Extract claim
            claim_match = re.search(r'<claim>(.*?)</claim>', block, re.DOTALL)
            claim = claim_match.group(1).strip() if claim_match else ""

            # Extract confidence
            confidence_match = re.search(r'<confidence>(.*?)</confidence>', block)
            confidence = confidence_match.group(1).strip() if confidence_match else "medium"

            # Extract node CURIEs from LLM output
            node_curies = []
            node_curies_match = re.search(r'<node_curies>(.*?)</node_curies>', block, re.DOTALL)
            if node_curies_match:
                try:
                    parsed = json.loads(node_curies_match.group(1).strip())
                    # Flatten and filter to only strings (handles nested lists, non-string items)
                    node_curies = self._flatten_to_strings(parsed)
                except json.JSONDecodeError:
                    logger.warning(f"Citation {citation_id}: Failed to parse node_curies JSON")
                    # Try to extract CURIEs using regex as fallback
                    curie_pattern = r'"([^"]+)"'
                    node_curies = re.findall(curie_pattern, node_curies_match.group(1))

            # Extract ALL edges, publications, and sentences from graph for cited nodes
            # Include query genes and disease for complete pathway traceability
            edge_ids, publication_ids, sentences, avg_importance = self._extract_edges_for_nodes(
                graph, node_curies, query_genes=query_genes, disease_curie=disease_curie
            )

            # Auto-adjust confidence based on extracted evidence
            # Override LLM confidence if evidence strongly disagrees
            pub_count = len(publication_ids)
            if avg_importance > 7.0 and pub_count > 5:
                adjusted_confidence = "high"
            elif avg_importance > 4.0 or pub_count > 2:
                adjusted_confidence = "medium"
            elif pub_count < 2:
                adjusted_confidence = "low"
            else:
                adjusted_confidence = confidence  # Keep LLM's assessment

            citations.append(CitationGraph(
                citation_id=citation_id,
                claim=claim,
                node_ids=node_curies,
                edge_ids=edge_ids,
                publication_ids=publication_ids,
                sentences=sentences,
                confidence=adjusted_confidence
            ))

        logger.info(f"Parsed {len(citations)} citations from XML response")
        return citations

    def _extract_edges_for_nodes(
        self,
        graph: nx.MultiDiGraph,
        node_curies: List[str],
        query_genes: Optional[List[str]] = None,
        disease_curie: Optional[str] = None
    ) -> Tuple[List[str], List[str], List[str], float]:
        """Extract all edges, publications, and sentences connecting the cited nodes.

        Only includes nodes explicitly cited by the LLM. Query genes and disease
        are NOT automatically added - they should only appear if the LLM included
        them in the citation's node_curies list.

        Args:
            graph: Full knowledge graph
            node_curies: List of node CURIEs from the citation
            query_genes: Unused, kept for API compatibility
            disease_curie: Unused, kept for API compatibility

        Returns:
            Tuple of (edge_ids, publication_ids, sentences, avg_importance_score)
        """
        # Only use nodes explicitly cited by the LLM
        node_set = set(node_curies)

        edge_ids = []
        all_publications = []
        all_sentences = []
        importance_scores = []

        # Find all edges where BOTH endpoints are in the cited node set
        for u, v, key, data in graph.edges(keys=True, data=True):
            if u in node_set and v in node_set:
                # Format: subject→predicate→object (semantic triple order)
                predicate_for_id = key.rsplit('_', 1)[0] if '_' in key else key
                edge_id = f"{u}→{predicate_for_id}→{v}"
                edge_ids.append(edge_id)

                # Collect publications from this edge
                pubs = data.get('publications', [])
                all_publications.extend(pubs)

                # Collect sentences/supporting_text from this edge
                sentences = data.get('sentences', [])
                if sentences:
                    all_sentences.extend(sentences)

                # Track importance scores for confidence adjustment
                importance = data.get('importance_score', 0.0)
                importance_scores.append(importance)

        # Deduplicate publications (these should be strings like "PMID:12345")
        unique_publications = list(set(p for p in all_publications if isinstance(p, str)))

        # Deduplicate sentences - handle potential non-string items
        seen_sentences = set()
        unique_sentences = []
        for s in all_sentences:
            if isinstance(s, str) and s not in seen_sentences:
                seen_sentences.add(s)
                unique_sentences.append(s)

        # Calculate average importance
        avg_importance = sum(importance_scores) / len(importance_scores) if importance_scores else 0.0

        return edge_ids, unique_publications, unique_sentences, avg_importance

    def _extract_summary(self, response_text: str, category: str = "") -> str:
        """Extract summary from LLM response using multiple strategies.

        Tries extraction patterns in order of priority:
        1. Exact <summary>...</summary> tags
        2. Case-insensitive summary tags
        3. Summary tags with attributes
        4. Markdown "## Summary" or "**Summary:**" sections
        5. Everything after the last </save_citation_graph> tag
        6. Paragraphs containing [Citation N] markers

        Args:
            response_text: Full LLM response text
            category: Category name for logging

        Returns:
            Extracted summary text (never empty - falls back to heuristic extraction)
        """
        # Log raw response for debugging
        self._log_raw_response(response_text, category)

        # Strategy 1: Exact <summary>...</summary> tags
        match = re.search(r'<summary>(.*?)</summary>', response_text, re.DOTALL)
        if match:
            logger.info(f"Summary extracted using Strategy 1 (exact tags)")
            return match.group(1).strip()

        # Strategy 2: Case-insensitive summary tags
        match = re.search(r'<summary>(.*?)</summary>', response_text, re.DOTALL | re.IGNORECASE)
        if match:
            logger.info(f"Summary extracted using Strategy 2 (case-insensitive)")
            return match.group(1).strip()

        # Strategy 3: Summary tags with attributes (e.g., <summary type="...">)
        match = re.search(r'<summary[^>]*>(.*?)</summary>', response_text, re.DOTALL | re.IGNORECASE)
        if match:
            logger.info(f"Summary extracted using Strategy 3 (tags with attributes)")
            return match.group(1).strip()

        # Strategy 4: Markdown section headers
        # Look for "## Summary", "### Summary", "**Summary:**", etc.
        patterns = [
            r'##\s*Summary\s*\n(.*?)(?=\n##|\n\*\*|$)',  # ## Summary header
            r'\*\*Summary:?\*\*\s*\n?(.*?)(?=\n\*\*|\n##|$)',  # **Summary:** bold
            r'Summary:\s*\n(.*?)(?=\n##|\n\*\*|$)',  # Plain "Summary:" label
        ]
        for i, pattern in enumerate(patterns, start=1):
            match = re.search(pattern, response_text, re.DOTALL | re.IGNORECASE)
            if match:
                text = match.group(1).strip()
                if len(text) > 100:  # Minimum viable summary length
                    logger.info(f"Summary extracted using Strategy 4.{i} (markdown section)")
                    return text

        # Strategy 5: Extract everything after the last citation block
        # This assumes the summary comes after all citations
        last_citation_end = response_text.rfind('</save_citation_graph>')
        if last_citation_end != -1:
            after_citations = response_text[last_citation_end + len('</save_citation_graph>'):].strip()
            # Clean up markdown/XML artifacts
            after_citations = re.sub(r'^#+\s*TURN\s*2.*?\n', '', after_citations, flags=re.IGNORECASE)
            after_citations = re.sub(r'^#+\s*GENERATE\s*SUMMARY.*?\n', '', after_citations, flags=re.IGNORECASE)
            after_citations = after_citations.strip()
            if len(after_citations) > 100:
                logger.info(f"Summary extracted using Strategy 5 (after last citation)")
                return after_citations

        # Strategy 6: If nothing else works, try to find the longest paragraph
        # that contains citation markers [Citation N]
        paragraphs = response_text.split('\n\n')
        citation_paragraphs = [p for p in paragraphs if re.search(r'\[Citation \d+\]', p)]
        if citation_paragraphs:
            # Join paragraphs with citations
            summary = '\n\n'.join(citation_paragraphs)
            logger.info(f"Summary extracted using Strategy 6 (citation paragraphs)")
            return summary.strip()

        logger.error(f"All extraction strategies failed for {category}")
        return ""

    def _log_raw_response(self, response_text: str, category: str) -> None:
        """Log raw LLM response to file for debugging.

        Args:
            response_text: Full LLM response text
            category: Category name for filename
        """
        try:
            from pathlib import Path
            from datetime import datetime as dt

            log_dir = Path("data/cache/llm_debug")
            log_dir.mkdir(parents=True, exist_ok=True)

            timestamp = dt.now().strftime("%Y%m%d_%H%M%S")
            log_file = log_dir / f"raw_response_{category}_{timestamp}.txt"

            with open(log_file, 'w') as f:
                f.write(f"Category: {category}\n")
                f.write(f"Timestamp: {timestamp}\n")
                f.write(f"Response length: {len(response_text)} chars\n")
                f.write("=" * 80 + "\n")
                f.write(response_text)

            logger.info(f"Raw LLM response logged to {log_file}")
        except Exception as e:
            logger.warning(f"Failed to log raw response: {e}")
