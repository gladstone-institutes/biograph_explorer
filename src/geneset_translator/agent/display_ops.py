"""Display-graph state holder + pure NetworkX transforms for agent-driven graph control.

The agent edits what is *shown* (trim, focus, add the disease node, merge results) through these
operations rather than re-querying. ``DisplayGraph`` is the thread-safe, versioned state holder
(the chat page renders it); the module-level functions are pure transforms on the normalized graph
shape that ``ui.network_viz`` consumes (see ``viz_normalizer``). No TCT, no Streamlit imports.
"""

from __future__ import annotations

import threading
from typing import Any, Callable, Iterable, Optional, Protocol, Union, runtime_checkable

import networkx as nx

from . import viz_normalizer

Graph = Union[nx.DiGraph, nx.MultiDiGraph]


# --------------------------------------------------------------------------------------
# State holder
# --------------------------------------------------------------------------------------
@runtime_checkable
class DisplayState(Protocol):
    """Abstraction the display tools depend on (concrete impl: ``DisplayGraph``)."""

    version: int

    def is_empty(self) -> bool: ...
    def snapshot(self) -> nx.MultiDiGraph: ...
    def replace(self, graph: Graph) -> None: ...
    def reset(self) -> None: ...
    def apply(self, fn: Callable[..., nx.MultiDiGraph], *args: Any, **kwargs: Any) -> None: ...


class DisplayGraph:
    """Thread-safe, versioned holder for the curated display graph.

    Tool calls in one agent step may run concurrently, so every mutation is locked. ``version``
    bumps on every change and is used by the chat page to invalidate its element cache while
    keeping the streamlit-cytoscape component ``key`` stable (in-place update, preserved viewport).
    """

    def __init__(self) -> None:
        self._g: nx.MultiDiGraph = nx.MultiDiGraph()
        self._lock = threading.Lock()
        self.version = 0

    def is_empty(self) -> bool:
        with self._lock:
            return self._g.number_of_nodes() == 0

    def snapshot(self) -> nx.MultiDiGraph:
        with self._lock:
            return self._g.copy()

    def replace(self, graph: Graph) -> None:
        with self._lock:
            self._g = nx.MultiDiGraph(graph)
            self.version += 1

    def reset(self) -> None:
        with self._lock:
            self._g = nx.MultiDiGraph()
            self.version += 1

    def apply(self, fn: Callable[..., nx.MultiDiGraph], *args: Any, **kwargs: Any) -> None:
        """Replace the held graph with ``fn(current, *args, **kwargs)`` under the lock."""
        with self._lock:
            self._g = fn(self._g, *args, **kwargs)
            self.version += 1


# --------------------------------------------------------------------------------------
# Pure transforms (return a new graph; preserve the renderer's node/edge attributes)
# --------------------------------------------------------------------------------------
def _neighbors(graph: Graph, node: str) -> set:
    if node not in graph:
        return set()
    return set(graph.predecessors(node)) | set(graph.successors(node))


def filter_connected_to(graph: Graph, curie: str) -> nx.MultiDiGraph:
    """Keep ``curie`` and the nodes directly connected to it (e.g. 'genes connected to this drug')."""
    keep = {curie} | _neighbors(graph, curie)
    return nx.MultiDiGraph(graph.subgraph(keep).copy())


def trim_to_top_degree(
    graph: Graph, top_n: int, always_keep: Optional[Iterable[str]] = None
) -> nx.MultiDiGraph:
    """Keep the ``top_n`` highest-degree nodes plus ``always_keep`` (query genes / disease)."""
    always = set(always_keep or [])
    ranked = sorted(graph.degree, key=lambda nd: nd[1], reverse=True)
    keep = {n for n, _ in ranked[: max(top_n, 0)]} | (always & set(graph.nodes))
    return nx.MultiDiGraph(graph.subgraph(keep).copy())


def trim_to_top_per_cluster(
    graph: Graph, k: int, always_keep: Optional[Iterable[str]] = None
) -> nx.MultiDiGraph:
    """Keep the ``k`` highest-degree nodes WITHIN each cluster (node attr ``cluster``, set by the
    cluster_graph tool) plus ``always_keep`` -- a balanced per-cluster view, unlike trim_to_top_degree
    which keeps the highest-degree nodes overall (dominated by the largest cluster on hub graphs)."""
    always = set(always_keep or [])
    by_cluster: dict = {}
    for node, data in graph.nodes(data=True):
        cid = data.get("cluster")
        if cid is not None:
            by_cluster.setdefault(cid, []).append(node)
    keep = always & set(graph.nodes)
    for nodes in by_cluster.values():
        ranked = sorted(nodes, key=lambda n: graph.degree(n), reverse=True)
        keep.update(ranked[: max(k, 0)])
    return nx.MultiDiGraph(graph.subgraph(keep).copy())


def filter_by_category(
    graph: Graph, category: str, always_keep: Optional[Iterable[str]] = None
) -> nx.MultiDiGraph:
    """Keep nodes of ``category`` (short form, e.g. 'ChemicalEntity') plus ``always_keep``."""
    always = set(always_keep or [])
    keep = {
        n for n, d in graph.nodes(data=True) if d.get("category") == category
    } | (always & set(graph.nodes))
    return nx.MultiDiGraph(graph.subgraph(keep).copy())


def merge_graphs(g1: Graph, g2: Graph) -> nx.MultiDiGraph:
    """Union of two normalized graphs (nodes + edges); g2 attributes win on overlap."""
    return nx.compose(nx.MultiDiGraph(g1), nx.MultiDiGraph(g2))


def add_disease_node(
    graph: Graph,
    disease_curie: str,
    disease_label: Optional[str],
    associated_gene_curies: Iterable[str],
) -> nx.MultiDiGraph:
    """Add the disease as a node and link it to the associated genes already in the graph."""
    out = nx.MultiDiGraph(graph)
    out.add_node(
        disease_curie,
        label=disease_label or disease_curie,
        curie=disease_curie,
        category="Disease",
        is_query_gene=False,
        is_disease_associated_bp=False,
        synonyms=[],
    )
    for gene in associated_gene_curies:
        if gene in out and gene != disease_curie:
            out.add_edge(
                disease_curie,
                gene,
                key=f"associated_with_{gene}",
                predicate="biolink:associated_with",
                sources=[],
                publications=[],
                sentences=[],
                confidence_scores={},
                knowledge_level=None,
                agent_type=None,
                qualifiers=[],
                attributes=[],
            )
    return out


def drop_orphans(graph: Graph) -> nx.MultiDiGraph:
    """Return a copy with orphan nodes removed -- nodes with no edge to any OTHER node (degree 0 or
    self-loops only). Floating dots carry no information and clutter the view. Never mutates the input."""
    out = nx.MultiDiGraph(graph)
    orphans = [
        n for n in out.nodes
        if not ((set(out.predecessors(n)) | set(out.successors(n))) - {n})
    ]
    out.remove_nodes_from(orphans)
    return out


def finalize(graph: Graph, query_gene_curies: Iterable[str]) -> nx.MultiDiGraph:
    """Recompute gene_frequency and backfill category/is_query_gene so the result renders cleanly.
    Orphan nodes (no edge to any other node) are dropped so the agent never displays floating dots."""
    out = drop_orphans(graph)
    query_set = set(query_gene_curies)
    for node, data in out.nodes(data=True):
        data.setdefault("curie", node)
        if not data.get("category"):
            data["category"] = "Gene" if node in query_set else viz_normalizer.classify_category(node)
        data.setdefault("is_query_gene", node in query_set)
        data.setdefault("label", node)
    for node, freq in viz_normalizer.compute_gene_frequency(out, list(query_set)).items():
        out.nodes[node]["gene_frequency"] = freq
    return out
