"""Agent tools: thin wrappers over the TctGateway (plus the HPA/annotation extension).

Each tool is a small object implementing the ``Tool`` protocol (name, description,
input_schema, run). ``ToolRegistry`` exposes the Anthropic tool schemas and a uniform
``dispatch`` that catches errors into agent-visible feedback and records successful calls.
Adding a tool means adding one object to the registry; the agent loop never changes.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Protocol, runtime_checkable

from . import display_ops
from .display_ops import DisplayState
from .graph_clustering import summarize_graph
from .tct_adapter import (
    HUMAN_TAXON,
    ResultStash,
    TctGateway,
    _knowledge_graphs,
    tct_result_to_nx,
)

logger = logging.getLogger(__name__)


@dataclass
class ToolContext:
    """Injected dependencies + per-session mutable state a tool may read/write."""

    tct: TctGateway
    stash: ResultStash
    query_gene_curies: List[str] = field(default_factory=list)
    curie_to_symbol: Dict[str, str] = field(default_factory=dict)
    disease_curie: Optional[str] = None
    disease_label: Optional[str] = None
    annotations: Dict[str, Dict[str, Any]] = field(default_factory=dict)
    tool_log: List[Dict[str, Any]] = field(default_factory=list)
    latest_result_id: Optional[str] = None
    query_cache: Dict[str, Any] = field(default_factory=dict)  # (tool,args) -> result, saves API calls
    display: DisplayState = field(default_factory=display_ops.DisplayGraph)


@runtime_checkable
class Tool(Protocol):
    name: str
    description: str
    input_schema: Dict[str, Any]
    parallel_safe: bool  # True -> may run concurrently in a step; False -> mutates shared display state

    def run(self, args: Dict[str, Any], ctx: ToolContext) -> Dict[str, Any]: ...


def _df_head_records(df: Any, n: int) -> "tuple[List[dict], int]":
    """Compact a pandas DataFrame to head(n) records; never raise."""
    if df is None:
        return [], 0
    try:
        return df.head(n).to_dict("records"), int(len(df))
    except Exception:  # noqa: BLE001
        return [], 0


def _edge_publications(edge: Dict[str, Any]) -> List[str]:
    """Up to a few normalized publication links for an edge, from either the direct ``publications``
    key (synthetic/recorded edges) or the TRAPI ``attributes`` (live TCT edges)."""
    from geneset_translator.utils.publication_utils import normalize_publication_id

    raw = edge.get("publications")
    if not raw and edge.get("attributes"):
        try:
            from TCT.attribute_extraction import extract_rich_edge_attributes

            raw = extract_rich_edge_attributes(edge["attributes"]).get("publications")
        except Exception:  # noqa: BLE001
            raw = None
    out: List[str] = []
    for pub in raw or []:
        norm = normalize_publication_id(str(pub))
        if norm and norm.isdigit():  # TCT often returns bare-digit PMIDs -> make them linkable
            norm = f"PMID:{norm}"
        if norm and norm not in out:
            out.append(norm)
    return out


def _edge_primary_source(edge: Dict[str, Any]) -> Optional[str]:
    """The primary knowledge source infores id for an edge (TRAPI ``sources`` or the
    synthetic ``primary_sources`` key)."""
    for src in edge.get("sources", []) or []:
        if isinstance(src, dict) and src.get("resource_role") == "primary_knowledge_source":
            return src.get("resource_id")
    prim = edge.get("primary_sources") or []
    return prim[0] if prim else None


def _edge_evidence_index(result: Any, max_pubs: int = 3) -> Dict[tuple, Dict[str, Any]]:
    """Map (subject, object) -> {publications, primary_source} from a result's knowledge graph(s),
    so finder rows can carry citable evidence without an extra edge_evidence call."""
    from .tct_adapter import _knowledge_graphs

    index: Dict[tuple, Dict[str, Any]] = {}
    for kg in _knowledge_graphs(result):
        try:
            items = list(kg.items())
        except Exception:  # noqa: BLE001
            continue
        for _eid, edge in items:
            if not isinstance(edge, dict):
                continue
            key = (edge.get("subject"), edge.get("object"))
            entry = index.setdefault(key, {"publications": [], "primary_source": None})
            for pub in _edge_publications(edge):
                if pub not in entry["publications"]:
                    entry["publications"].append(pub)
            if not entry["primary_source"]:
                entry["primary_source"] = _edge_primary_source(edge)
    for entry in index.values():
        entry["publications"] = entry["publications"][:max_pubs]
    return index


# --------------------------------------------------------------------------------------
# Tools
# --------------------------------------------------------------------------------------
class ResolveGenesTool:
    name = "resolve_genes"
    description = (
        "Resolve gene symbols or entity names to CURIEs. ALWAYS call this FIRST. Every other tool "
        "requires CURIEs (e.g. NCBIGene:2322), never bare symbols. Resolves to human genes by "
        "default. When resolving a DISEASE or DRUG by name, pass biolink_type (e.g. "
        "'biolink:Disease', 'biolink:Drug') so the name is not mis-matched to a same-named gene "
        "(e.g. 'acute myeloid leukemia' would otherwise resolve to the gene RUNX1)."
    )
    input_schema = {
        "type": "object",
        "properties": {
            "names": {
                "type": "array",
                "items": {"type": "string"},
                "description": "Gene symbols or entity names, e.g. ['FLT3','BCL2','venetoclax'].",
            },
            "biolink_type": {
                "type": "string",
                "description": (
                    "Optional biolink type to constrain resolution for NON-gene names, e.g. "
                    "'biolink:Disease' for a disease name or 'biolink:Drug'/'biolink:SmallMolecule' "
                    "for a drug. Omit for genes (the human-gene default)."
                ),
            },
        },
        "required": ["names"],
    }

    def run(self, args: Dict[str, Any], ctx: ToolContext) -> Dict[str, Any]:
        resolved = ctx.tct.resolve_genes(
            args["names"],
            only_taxa=args.get("only_taxa", HUMAN_TAXON),
            biolink_type=args.get("biolink_type"),
        )
        # Remember CURIE -> user-facing symbol so rendered graphs show friendly names.
        for name, info in resolved.items():
            if info and info.get("curie"):
                ctx.curie_to_symbol.setdefault(info["curie"], name)
        return resolved


class GeneNeighborhoodTool:
    name = "gene_neighborhood"
    description = (
        "Find the 1-hop neighborhood of ONE entity, filtered to target categories. Works on ANY "
        "resolved CURIE, not just genes: gene -> drugs/chemicals (the default categories), "
        "disease -> genes (pass target_categories ['biolink:Gene']), drug -> targets (pass "
        "['biolink:Gene','biolink:Protein']). gene_curie MUST be a CURIE from resolve_genes."
    )
    input_schema = {
        "type": "object",
        "properties": {
            "gene_curie": {
                "type": "string",
                "description": (
                    "Any entity CURIE from resolve_genes (gene, disease, or drug), e.g. "
                    "NCBIGene:2322 or MONDO:0100096."
                ),
            },
            "target_categories": {
                "type": "array",
                "items": {"type": "string"},
                "description": (
                    "biolink categories for the neighbor nodes. Omit for the drug/chemical default "
                    "(good for gene -> drugs). For disease -> genes pass ['biolink:Gene']; for "
                    "drug -> targets pass ['biolink:Gene','biolink:Protein']."
                ),
            },
        },
        "required": ["gene_curie"],
    }

    def run(self, args: Dict[str, Any], ctx: ToolContext) -> Dict[str, Any]:
        nb = ctx.tct.neighborhood(args["gene_curie"], args.get("target_categories"))
        result_id = ctx.stash.put(nb)
        ctx.display.reset()  # new base graph; clears any prior display curation
        top, n = _df_head_records(getattr(nb, "ranked", None), 25)
        input_node = getattr(nb, "input_node_id", args["gene_curie"])
        index = _edge_evidence_index(nb)
        for row in top[:10]:  # attach citable evidence to the top edges (bounded for token budget)
            ev = index.get((input_node, row.get("output_node")))
            if ev and (ev["publications"] or ev["primary_source"]):
                row["publications"] = ev["publications"]
                row["primary_source"] = ev["primary_source"]
        return {
            "result_id": result_id,
            "input_node_id": input_node,
            "n_results": n,
            "top_edges": top,
        }


class PathBetweenTool:
    name = "path_between"
    description = (
        "Find connecting paths between TWO entities (gene, drug, disease, ...). Both inputs MUST "
        "be CURIEs from resolve_genes. intermediate_categories are the allowed middle-node types."
    )
    input_schema = {
        "type": "object",
        "properties": {
            "node1_curie": {"type": "string", "description": "First entity CURIE."},
            "node2_curie": {"type": "string", "description": "Second entity CURIE."},
            "intermediate_categories": {
                "type": "array",
                "items": {"type": "string"},
                "description": (
                    "Allowed intermediate node categories, e.g. "
                    "['biolink:Gene','biolink:Protein','biolink:ChemicalEntity']."
                ),
            },
        },
        "required": ["node1_curie", "node2_curie"],
    }

    def run(self, args: Dict[str, Any], ctx: ToolContext) -> Dict[str, Any]:
        intermediates = args.get("intermediate_categories") or [
            "biolink:Gene",
            "biolink:Protein",
            "biolink:ChemicalEntity",
        ]
        pr = ctx.tct.path(args["node1_curie"], args["node2_curie"], intermediates)
        result_id = ctx.stash.put(pr)
        ctx.display.reset()  # new base graph; clears any prior display curation
        top, n = _df_head_records(getattr(pr, "paths", None), 25)
        # A path spans several edges; surface a small pool of the path's publication links so the
        # model can cite the connection (use edge_evidence for a specific intermediate edge).
        index = _edge_evidence_index(pr, max_pubs=2)
        evidence_pubs: List[str] = []
        for entry in index.values():
            for pub in entry["publications"]:
                if pub not in evidence_pubs:
                    evidence_pubs.append(pub)
        return {
            "result_id": result_id,
            "node1_id": getattr(pr, "node1_id", args["node1_curie"]),
            "node2_id": getattr(pr, "node2_id", args["node2_curie"]),
            "n_paths": n,
            "top_paths": top,
            "evidence_publications": evidence_pubs[:8],
        }


class GeneNetworkTool:
    name = "gene_network"
    description = (
        "Find direct relationships AMONG a set of genes (a gene-gene interaction network). "
        "gene_curies MUST be CURIEs from resolve_genes. Returns top_hubs (the input genes ranked "
        "by their degree in this network) so you can pick the 'most connected' genes from the data."
    )
    input_schema = {
        "type": "object",
        "properties": {
            "gene_curies": {
                "type": "array",
                "items": {"type": "string"},
                "description": "Gene CURIEs, e.g. ['NCBIGene:2322','NCBIGene:596'].",
            }
        },
        "required": ["gene_curies"],
    }

    def run(self, args: Dict[str, Any], ctx: ToolContext) -> Dict[str, Any]:
        gene_curies = args["gene_curies"]
        kg = ctx.tct.gene_network(gene_curies)
        result_id = ctx.stash.put(kg)
        ctx.display.reset()  # new base graph; clears any prior display curation

        degree: Dict[str, int] = {g: 0 for g in gene_curies}
        top: List[dict] = []
        try:
            for _edge_id, edge in kg.items():
                if not isinstance(edge, dict):
                    continue
                subj, obj = edge.get("subject"), edge.get("object")
                if subj in degree:
                    degree[subj] += 1
                if obj in degree:
                    degree[obj] += 1
                if len(top) < 25:
                    row = {"subject": subj, "predicate": edge.get("predicate"), "object": obj}
                    if len(top) < 10:  # citable evidence for the top edges (bounded for tokens)
                        pubs = _edge_publications(edge)
                        prim = _edge_primary_source(edge)
                        if pubs:
                            row["publications"] = pubs
                        if prim:
                            row["primary_source"] = prim
                    top.append(row)
        except Exception:  # noqa: BLE001
            pass
        try:
            edge_count = len(kg)
        except Exception:  # noqa: BLE001
            edge_count = len(top)

        top_hubs = [
            {"gene_curie": g, "symbol": ctx.curie_to_symbol.get(g, g), "degree": d}
            for g, d in sorted(degree.items(), key=lambda kv: kv[1], reverse=True)
        ][:10]
        return {
            "result_id": result_id,
            "edge_count": edge_count,
            "top_hubs": top_hubs,
            "top_edges": top,
        }


class EdgeEvidenceTool:
    name = "edge_evidence"
    description = (
        "Get supporting publications, sentences, and confidence scores for the edge(s) between "
        "two nodes in a previous result. Pass the finder's result_id and the subject/object CURIEs."
    )
    input_schema = {
        "type": "object",
        "properties": {
            "result_id": {"type": "string", "description": "result_id from a finder tool."},
            "subject": {"type": "string", "description": "Subject node CURIE."},
            "object": {"type": "string", "description": "Object node CURIE."},
        },
        "required": ["result_id", "subject", "object"],
    }

    def run(self, args: Dict[str, Any], ctx: ToolContext) -> Dict[str, Any]:
        result = ctx.stash.get(args["result_id"])
        if result is None:
            raise ValueError(
                f"unknown result_id {args['result_id']!r}; run a finder tool first to create one"
            )
        return ctx.tct.edge_evidence(result, args["subject"], args["object"])


def _expression_key_filter(scope: str):
    """Return a predicate selecting which HPA/annotation feature keys to keep for a scope.

    HPA writes its features under ``hpa_*`` keys (see hpa_client._apply_hpa_annotations):
    hpa_cell_type_*, hpa_tissue_*/hpa_top_tissues, hpa_immune_*/hpa_top_immune_*.
    """
    if scope == "cell_type":
        return lambda k: "cell_type" in k and "immune" not in k
    if scope == "tissue":
        return lambda k: "tissue" in k
    if scope == "immune":
        return lambda k: "immune" in k
    # "all": every HPA expression key plus GO terms and protein class.
    return lambda k: k.startswith("hpa_") or k.startswith("go_") or "protein_class" in k


def _compact_expression(features: Dict[str, Any], scope: str = "all") -> Dict[str, Any]:
    """Pick the expression/annotation fields for the requested scope and truncate for context."""
    keep = _expression_key_filter(scope)
    out: Dict[str, Any] = {}
    for key, value in features.items():
        if keep(key.lower()):
            out[key] = value[:5] if isinstance(value, list) else value
    return out


class CellTypeExpressionTool:
    name = "cell_type_expression"
    description = (
        "Get Human Protein Atlas expression and Node Annotator metadata (e.g. GO terms) for the "
        "genes in a previous result. Returns cell-type, TISSUE, and IMMUNE-lineage specificity. "
        "Call this when the user asks about expression, cell type, tissue, or immune/blood-cell "
        "specificity. Use scope to narrow it (cell_type | tissue | immune | all). Pass the "
        "finder's result_id."
    )
    input_schema = {
        "type": "object",
        "properties": {
            "result_id": {"type": "string", "description": "result_id from a finder tool."},
            "gene_curies": {
                "type": "array",
                "items": {"type": "string"},
                "description": "Optional subset of gene CURIEs to summarize; defaults to the gene set.",
            },
            "scope": {
                "type": "string",
                "enum": ["cell_type", "tissue", "immune", "all"],
                "description": (
                    "Which expression facet to return: cell_type (single-cell), tissue, immune "
                    "(blood/immune lineage), or all (default)."
                ),
            },
        },
        "required": ["result_id"],
    }

    def run(self, args: Dict[str, Any], ctx: ToolContext) -> Dict[str, Any]:
        result = ctx.stash.get(args["result_id"])
        if result is None:
            raise ValueError(f"unknown result_id {args['result_id']!r}; run a finder tool first")

        from geneset_translator.core.hpa_client import HPAClient
        from geneset_translator.core.node_annotator import NodeAnnotator

        graph = tct_result_to_nx(
            result, ctx.query_gene_curies, ctx.curie_to_symbol, ctx.disease_curie
        )
        try:
            graph, _ = HPAClient().annotate_graph(graph)
        except Exception as e:  # noqa: BLE001
            logger.warning("HPA annotation failed: %s", e)
        try:
            graph, _ = NodeAnnotator().annotate_graph(graph)
        except Exception as e:  # noqa: BLE001
            logger.warning("Node annotation failed: %s", e)

        for node, data in graph.nodes(data=True):
            features = data.get("annotation_features")
            if features:
                ctx.annotations[node] = features

        scope = args.get("scope", "all")
        genes = args.get("gene_curies") or ctx.query_gene_curies
        expression = {
            gene: _compact_expression(ctx.annotations[gene], scope)
            for gene in genes
            if gene in ctx.annotations
        }
        return {
            "result_id": args["result_id"],
            "scope": scope,
            "annotated_nodes": len(ctx.annotations),
            "expression": expression,
        }


class NodeMetadataTool:
    name = "node_metadata"
    description = (
        "Get functional metadata for one or more entities: Gene Ontology terms (biological "
        "process / molecular function / cellular component), gene type, and aliases. Call this for "
        "'what does gene X do', 'what processes/functions is it involved in', or 'what is its gene "
        "type'. Standalone: pass CURIEs directly (no finder result needed first)."
    )
    input_schema = {
        "type": "object",
        "properties": {
            "node_curies": {
                "type": "array",
                "items": {"type": "string"},
                "description": "Entity CURIEs from resolve_genes, e.g. ['NCBIGene:2322'].",
            }
        },
        "required": ["node_curies"],
    }

    def run(self, args: Dict[str, Any], ctx: ToolContext) -> Dict[str, Any]:
        from geneset_translator.core.node_annotator import NodeAnnotator

        curies = args["node_curies"]
        annotator = NodeAnnotator()
        nodes = annotator.annotate_nodes(curies)
        out: Dict[str, Any] = {}
        for curie in curies:
            node = nodes.get(curie)
            if node is None:
                continue
            feats = annotator._extract_annotation_features(node)
            compact = {"label": node.label or curie}
            if feats.get("type_of_gene"):
                compact["type_of_gene"] = feats["type_of_gene"]
            for go_key in ("go_bp", "go_mf", "go_cc"):
                terms = feats.get(go_key)
                if terms:
                    compact[go_key] = terms[:10]
            aliases = feats.get("alias")
            if aliases:
                compact["aliases"] = aliases[:10]
            out[curie] = compact
        return {"node_metadata": out, "n": len(out)}


class DataSourcesTool:
    name = "data_sources"
    description = (
        "Report which knowledge sources back a previous result and how reliable they are: the "
        "primary/aggregator knowledge sources, their knowledge levels (curated vs text-mined vs "
        "predicted), and the top contributors by edge count. Call this for 'where does this come "
        "from', 'which sources', or 'how reliable is this'. Pass the finder's result_id."
    )
    input_schema = {
        "type": "object",
        "properties": {
            "result_id": {"type": "string", "description": "result_id from a finder tool."},
        },
        "required": ["result_id"],
    }

    def run(self, args: Dict[str, Any], ctx: ToolContext) -> Dict[str, Any]:
        import networkx as nx

        from geneset_translator.config.settings import get_settings
        from geneset_translator.utils.infores_utils import (
            download_infores_catalog,
            extract_unique_sources,
            generate_source_summary,
            parse_infores_catalog,
        )

        result = ctx.stash.get(args["result_id"])
        if result is None:
            raise ValueError(f"unknown result_id {args['result_id']!r}; run a finder tool first")

        # Build a lightweight graph that PRESERVES the raw TRAPI ``sources`` per edge, so the
        # infores helpers (which read edge['sources']) work. The normalized render graph drops it.
        graph = nx.MultiDiGraph()
        for kg in _knowledge_graphs(result):
            for edge_id, edge in kg.items():
                if not isinstance(edge, dict):
                    continue
                graph.add_edge(
                    edge.get("subject"),
                    edge.get("object"),
                    key=edge_id,
                    sources=edge.get("sources", []),
                )

        relevant = extract_unique_sources(graph)
        cache_dir = get_settings().cache_dir / "infores"
        catalog = parse_infores_catalog(download_infores_catalog(cache_dir))
        filtered = {k: v for k, v in catalog.items() if k in relevant}
        summary = generate_source_summary(filtered, graph)
        summary["sources_in_result"] = len(relevant)
        summary["result_id"] = args["result_id"]
        return summary


# --------------------------------------------------------------------------------------
# Display-control tools (edit the on-screen graph in memory; no new query unless noted)
# --------------------------------------------------------------------------------------
def _working_display_graph(ctx: ToolContext):
    """The graph display tools operate on: the curated display if any, else the latest result."""
    if not ctx.display.is_empty():
        return ctx.display.snapshot()
    if ctx.latest_result_id:
        result = ctx.stash.get(ctx.latest_result_id)
        if result is not None:
            return tct_result_to_nx(
                result, ctx.query_gene_curies, ctx.curie_to_symbol, ctx.disease_curie
            )
    return None


def _always_keep(ctx: ToolContext) -> set:
    keep = set(ctx.query_gene_curies)
    if ctx.disease_curie:
        keep.add(ctx.disease_curie)
    return keep


class FilterGraphTool:
    name = "filter_graph"
    parallel_safe = False
    description = (
        "Change WHAT IS SHOWN in the current on-screen graph, without running a new query. Call "
        "this for 'show only', 'trim', 'focus on', 'genes connected to X', or 'the most connected'. "
        "connected_to=<CURIE> keeps that node and its direct neighbors; category keeps one node "
        "type (e.g. 'ChemicalEntity'); top_n keeps the most-connected nodes by degree OVERALL; "
        "top_per_cluster keeps the most-connected nodes WITHIN EACH cluster (run cluster_graph first) "
        "-- use this, not top_n, for a balanced 'top nodes of each cluster' view. Query genes and the "
        "disease are always kept."
    )
    input_schema = {
        "type": "object",
        "properties": {
            "connected_to": {"type": "string", "description": "CURIE; keep this node and its neighbors."},
            "category": {"type": "string", "description": "Short category to keep, e.g. ChemicalEntity, Gene."},
            "top_n": {"type": "integer", "description": "Keep the top_n most-connected nodes by degree OVERALL."},
            "top_per_cluster": {
                "type": "integer",
                "description": "Keep the top-K most-connected nodes within EACH cluster (needs cluster_graph first).",
            },
        },
    }

    def run(self, args: Dict[str, Any], ctx: ToolContext) -> Dict[str, Any]:
        graph = _working_display_graph(ctx)
        if graph is None or graph.number_of_nodes() == 0:
            raise ValueError("no graph to filter yet; run a finder tool first")
        before = graph.number_of_nodes()
        keep = _always_keep(ctx)
        applied: List[str] = []
        if args.get("connected_to"):
            graph = display_ops.filter_connected_to(graph, args["connected_to"])
            applied.append(f"connected_to={args['connected_to']}")
        if args.get("category"):
            graph = display_ops.filter_by_category(graph, args["category"], keep)
            applied.append(f"category={args['category']}")
        if args.get("top_n"):
            graph = display_ops.trim_to_top_degree(graph, int(args["top_n"]), keep)
            applied.append(f"top_{int(args['top_n'])}_by_degree")
        if args.get("top_per_cluster"):
            if not any(graph.nodes[n].get("cluster") is not None for n in graph.nodes):
                raise ValueError(
                    "no clusters on the current graph; run cluster_graph first, then filter by top_per_cluster"
                )
            graph = display_ops.trim_to_top_per_cluster(graph, int(args["top_per_cluster"]), keep)
            applied.append(f"top_{int(args['top_per_cluster'])}_per_cluster")
        if not applied:
            raise ValueError("specify at least one of: connected_to, category, top_n, top_per_cluster")
        graph = display_ops.finalize(graph, ctx.query_gene_curies)
        ctx.display.replace(graph)
        return {
            "applied": applied,
            "nodes_before": before,
            "nodes_after": graph.number_of_nodes(),
            "edges_after": graph.number_of_edges(),
        }


class AddDiseaseNodeTool:
    name = "add_disease_node"
    parallel_safe = False
    description = (
        "Add the disease as a node in the current graph and connect it to the genes already shown. "
        "Use for 'add the disease to the graph/picture'. Runs one disease->gene association query, "
        "then links the disease to the on-screen genes."
    )
    input_schema = {
        "type": "object",
        "properties": {
            "disease_curie": {
                "type": "string",
                "description": "Disease CURIE; defaults to the session disease if omitted.",
            }
        },
    }

    def run(self, args: Dict[str, Any], ctx: ToolContext) -> Dict[str, Any]:
        disease = args.get("disease_curie") or ctx.disease_curie
        if not disease:
            raise ValueError("no disease set; pass disease_curie or set a disease in the sidebar")
        graph = _working_display_graph(ctx)
        if graph is None or graph.number_of_nodes() == 0:
            raise ValueError("no graph to add the disease to; run a finder tool first")
        nb = ctx.tct.neighborhood(disease, ["biolink:Gene"])
        associated: List[str] = []
        ranked = getattr(nb, "ranked", None)
        try:
            if ranked is not None and "output_node" in ranked.columns:
                associated = [str(x) for x in ranked["output_node"].tolist()]
        except Exception:  # noqa: BLE001
            pass
        graph = display_ops.add_disease_node(graph, disease, ctx.disease_label, associated)
        graph = display_ops.finalize(graph, ctx.query_gene_curies)
        ctx.display.replace(graph)
        linked = len(set(graph.successors(disease))) if disease in graph else 0
        return {
            "disease": disease,
            "label": ctx.disease_label or disease,
            "linked_genes": linked,
            "nodes_after": graph.number_of_nodes(),
        }


class ShowResultTool:
    name = "show_result"
    parallel_safe = False
    description = (
        "Set or merge a previously collected result (by result_id) into the displayed graph, "
        "instead of only ever showing the latest. merge=true adds it to what is shown; merge=false "
        "(default) replaces. Use to bring back an earlier result or combine several."
    )
    input_schema = {
        "type": "object",
        "properties": {
            "result_id": {"type": "string", "description": "A result_id from a previous finder tool."},
            "merge": {"type": "boolean", "description": "True = add to the current graph; False = replace."},
        },
        "required": ["result_id"],
    }

    def run(self, args: Dict[str, Any], ctx: ToolContext) -> Dict[str, Any]:
        result = ctx.stash.get(args["result_id"])
        if result is None:
            raise ValueError(f"unknown result_id {args['result_id']!r}")
        graph = tct_result_to_nx(result, ctx.query_gene_curies, ctx.curie_to_symbol, ctx.disease_curie)
        merge = bool(args.get("merge", False))
        if merge and not ctx.display.is_empty():
            graph = display_ops.merge_graphs(ctx.display.snapshot(), graph)
        graph = display_ops.finalize(graph, ctx.query_gene_curies)
        ctx.display.replace(graph)
        return {"result_id": args["result_id"], "merged": merge, "nodes_after": graph.number_of_nodes()}


class ClusterGraphTool:
    name = "cluster_graph"
    parallel_safe = False  # writes the display (tags nodes with their cluster for recoloring)
    description = (
        "Summarize the STRUCTURE of a large result so you can describe it accurately instead of "
        "guessing from the top edges. Returns a topology-aware digest: communities (for a real "
        "interaction mesh, each with size + top members + dominant category), or -- when the graph is "
        "hub-dominated (a star around one gene/disease) -- node-category and edge-predicate facets "
        "instead, plus connected components and overall counts. Also recolors the on-screen graph by "
        "cluster. Call this after a finder returns a large result (hundreds+ of nodes/edges); base "
        "your summary on what it returns and do not invent modules it did not find. Optional result_id "
        "selects which result (defaults to the current graph)."
    )
    input_schema = {
        "type": "object",
        "properties": {
            "result_id": {
                "type": "string",
                "description": "A result_id from a previous finder; defaults to the current graph.",
            }
        },
    }

    def run(self, args: Dict[str, Any], ctx: ToolContext) -> Dict[str, Any]:
        rid = args.get("result_id")
        if rid:
            result = ctx.stash.get(rid)
            if result is None:
                raise ValueError(f"unknown result_id {rid!r}")
            graph = tct_result_to_nx(result, ctx.query_gene_curies, ctx.curie_to_symbol, ctx.disease_curie)
        else:
            graph = _working_display_graph(ctx)
        if graph is None or graph.number_of_nodes() == 0:
            raise ValueError("no graph to cluster yet; run a finder tool first")

        digest, node_cluster = summarize_graph(
            graph, ctx.query_gene_curies, ctx.disease_curie, ctx.curie_to_symbol
        )
        # Tag nodes with their cluster so network_viz can recolor by community/facet.
        for node, cid in node_cluster.items():
            if node in graph:
                graph.nodes[node]["cluster"] = cid
        graph = display_ops.finalize(graph, ctx.query_gene_curies)
        ctx.display.replace(graph)
        return digest


# --------------------------------------------------------------------------------------
# Registry
# --------------------------------------------------------------------------------------
class ToolRegistry:
    """Holds tools, exposes their Anthropic schemas, and dispatches calls uniformly."""

    def __init__(self, tools: List[Tool]) -> None:
        self._tools: Dict[str, Tool] = {t.name: t for t in tools}

    def schemas(self) -> List[Dict[str, Any]]:
        return [
            {"name": t.name, "description": t.description, "input_schema": t.input_schema}
            for t in self._tools.values()
        ]

    def names(self) -> List[str]:
        return list(self._tools)

    def is_parallel_safe(self, name: str) -> bool:
        tool = self._tools.get(name)
        return bool(getattr(tool, "parallel_safe", True)) if tool else True

    def dispatch(self, name: str, args: Dict[str, Any], ctx: ToolContext):
        """Run a tool by name. Returns (result, is_error). Never raises.

        Identical (tool, args) calls to parallel-safe (idempotent) tools are memoized in
        ``ctx.query_cache`` so repeats reuse the stashed result instead of re-hitting the network.
        Display tools (``parallel_safe=False``) are stateful and always execute (never cached).
        """
        tool = self._tools.get(name)
        if tool is None:
            return {"error": f"unknown tool: {name}"}, True

        cacheable = getattr(tool, "parallel_safe", True)
        cache_key = f"{name}:{json.dumps(args or {}, sort_keys=True, default=str)}"
        if cacheable and cache_key in ctx.query_cache:
            logger.info("tool %s: cache hit (no API call)", name)
            result = ctx.query_cache[cache_key]
            ctx.tool_log.append({"name": name, "args": args or {}})
            if isinstance(result, dict) and result.get("result_id"):
                ctx.latest_result_id = result["result_id"]
            return result, False

        try:
            result = tool.run(args or {}, ctx)
        except Exception as e:  # noqa: BLE001 - surface as agent-visible error for self-debug
            logger.exception("tool %s failed", name)
            return {"error": f"{type(e).__name__}: {e}"}, True

        if cacheable:
            ctx.query_cache[cache_key] = result
        ctx.tool_log.append({"name": name, "args": args or {}})
        if isinstance(result, dict) and result.get("result_id"):
            ctx.latest_result_id = result["result_id"]
        return result, False


def default_registry() -> ToolRegistry:
    """The standard tool set for the chat agent."""
    return ToolRegistry(
        [
            ResolveGenesTool(),
            GeneNeighborhoodTool(),
            PathBetweenTool(),
            GeneNetworkTool(),
            EdgeEvidenceTool(),
            CellTypeExpressionTool(),
            NodeMetadataTool(),
            DataSourcesTool(),
            FilterGraphTool(),
            AddDiseaseNodeTool(),
            ShowResultTool(),
            ClusterGraphTool(),
        ]
    )
