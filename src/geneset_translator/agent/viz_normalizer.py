"""Map a TCT ``to_networkx`` graph to the attribute shape ``ui.network_viz`` expects.

TCT's ``KnowledgeGraph.to_networkx(resolve_names=True, include_attributes=True)`` returns a
MultiDiGraph whose nodes carry only ``label`` and whose edges carry ``predicate``,
``primary_sources``, ``aggregator_sources``, ``publications``, ``supporting_text`` and
``confidence_scores`` (verified against commit f702893). The existing renderer
(``prepare_cytoscape_elements``) instead reads the richer "raw" shape that ``GraphBuilder``
produces and does the string formatting itself. This module is the single place that knows
how to translate between the two, plus it derives node ``category`` from the CURIE prefix
(TCT's graph has no category) and computes the ``gene_frequency`` convergence metric.

Pure NetworkX: no TCT, no Streamlit imports.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Union

import networkx as nx

# Short category names matching ui.network_viz.CATEGORY_COLORS / CATEGORY_ICONS keys.
# Ordered longest-prefix-first within each category is unnecessary because we match by
# explicit prefix membership below.
_PREFIX_TO_CATEGORY: Dict[str, str] = {
    # Genes
    "NCBIGene": "Gene",
    "HGNC": "Gene",
    "ENSEMBL": "Gene",
    "OMIM.GENE": "Gene",
    # Proteins
    "UniProtKB": "Protein",
    "PR": "Protein",
    # Chemicals / drugs
    "CHEBI": "ChemicalEntity",
    "CHEMBL.COMPOUND": "ChemicalEntity",
    "CHEMBL": "ChemicalEntity",
    "PUBCHEM.COMPOUND": "ChemicalEntity",
    "DRUGBANK": "ChemicalEntity",
    "UNII": "ChemicalEntity",
    "RXCUI": "ChemicalEntity",
    "KEGG": "ChemicalEntity",
    "GTOPDB": "ChemicalEntity",
    "HMDB": "ChemicalEntity",
    "INCHIKEY": "ChemicalEntity",
    "MESH": "ChemicalEntity",
    # Disease
    "MONDO": "Disease",
    "DOID": "Disease",
    "OMIM": "Disease",
    "ORPHANET": "Disease",
    "Orphanet": "Disease",
    "MEDDRA": "Disease",
    # Phenotype
    "HP": "PhenotypicFeature",
    "EFO": "PhenotypicFeature",
    # Process / pathway / function
    "GO": "BiologicalProcess",
    "REACT": "Pathway",
    "SMPDB": "Pathway",
    "PANTHER.PATHWAY": "Pathway",
    "WIKIPATHWAYS": "Pathway",
    # Anatomy / cell
    "UBERON": "AnatomicalEntity",
    "CL": "Cell",
    "CLO": "Cell",
}


def classify_category(curie: str) -> str:
    """Classify a CURIE into a short biolink-ish category from its prefix.

    Falls back to "Other" for unmapped prefixes (e.g. UMLS, which is ambiguous).
    """
    if not curie or ":" not in curie:
        return "Other"
    prefix = curie.split(":", 1)[0]
    return _PREFIX_TO_CATEGORY.get(prefix, "Other")


def _norm_publication(pub: Any) -> str:
    """Coerce a publication entry to a string the renderer can format.

    TCT returns some publications as bare integer PMIDs; the renderer calls ``.startswith`` on
    each, so non-strings crash it. Bare digits are prefixed ``PMID:`` so they render as links.
    """
    s = str(pub).strip()
    return f"PMID:{s}" if s.isdigit() else s


def _build_sources(primary: List[str], aggregator: List[str]) -> List[dict]:
    """Convert TCT primary/aggregator infores id lists to the renderer's source shape.

    ``network_viz`` expects each source as ``{resource_id, resource_role, upstream_resource_ids}``.
    """
    sources: List[dict] = []
    for rid in primary or []:
        sources.append(
            {
                "resource_id": rid,
                "resource_role": "primary_knowledge_source",
                "upstream_resource_ids": [],
            }
        )
    for rid in aggregator or []:
        sources.append(
            {
                "resource_id": rid,
                "resource_role": "aggregator_knowledge_source",
                "upstream_resource_ids": [],
            }
        )
    return sources


def compute_gene_frequency(
    graph: Union[nx.DiGraph, nx.MultiDiGraph], query_genes: List[str]
) -> Dict[str, int]:
    """Convergence metric: number of distinct query genes directly linked to each node.

    Mirrors ``GraphBuilder.calculate_gene_frequency`` but kept here so the agent path has
    no dependency back into ``core``.
    """
    query_set = set(query_genes)
    freq: Dict[str, int] = {}
    for node in graph.nodes():
        connected = set()
        for predecessor in graph.predecessors(node):
            if predecessor in query_set:
                connected.add(predecessor)
        for successor in graph.successors(node):
            if successor in query_set:
                connected.add(successor)
        freq[node] = len(connected)
    return freq


def normalize(
    tct_graph: Union[nx.DiGraph, nx.MultiDiGraph],
    query_gene_curies: List[str],
    curie_to_symbol: Optional[Dict[str, str]] = None,
    disease_curie: Optional[str] = None,
) -> nx.MultiDiGraph:
    """Return a MultiDiGraph carrying the node/edge attributes ``network_viz`` reads.

    Args:
        tct_graph: output of ``KnowledgeGraph.to_networkx(resolve_names=True, include_attributes=True)``.
        query_gene_curies: input gene CURIEs (marked as query genes / sized larger).
        curie_to_symbol: optional gene CURIE -> user symbol, used for ``original_symbol``.
        disease_curie: optional disease anchor (rendered as a triangle like the Classic UI).
    """
    curie_to_symbol = curie_to_symbol or {}
    query_set = set(query_gene_curies)

    out = nx.MultiDiGraph()

    for node, data in tct_graph.nodes(data=True):
        label = data.get("label") or node
        is_query = node in query_set
        category = "Gene" if is_query else classify_category(node)
        attrs = {
            "label": label,
            "curie": node,
            "category": category,
            "is_query_gene": is_query,
            "is_disease_associated_bp": bool(disease_curie) and node == disease_curie,
            "synonyms": [],
        }
        if node in curie_to_symbol:
            attrs["original_symbol"] = curie_to_symbol[node]
        out.add_node(node, **attrs)

    is_multi = tct_graph.is_multigraph()
    edge_iter = (
        tct_graph.edges(keys=True, data=True)
        if is_multi
        else ((u, v, None, d) for u, v, d in tct_graph.edges(data=True))
    )
    for i, (u, v, key, data) in enumerate(edge_iter):
        predicate = data.get("predicate", "") or ""
        edge_key = key if key is not None else f"{predicate.replace('biolink:', '') or 'edge'}_{i}"
        out.add_edge(
            u,
            v,
            key=edge_key,
            predicate=predicate,
            sources=_build_sources(
                data.get("primary_sources", []), data.get("aggregator_sources", [])
            ),
            publications=[_norm_publication(p) for p in (data.get("publications") or [])],
            sentences=[str(s) for s in (data.get("supporting_text") or []) if s],
            confidence_scores=dict(data.get("confidence_scores", {}) or {}),
            knowledge_level=None,
            agent_type=None,
            qualifiers=[],
            attributes=[],
        )

    for node, f in compute_gene_frequency(out, query_gene_curies).items():
        out.nodes[node]["gene_frequency"] = f

    return out
