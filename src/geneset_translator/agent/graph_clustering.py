"""Topology-aware structured summaries of a result graph, for the chat agent.

The agent never sees a large finder result in full (it gets only a top slice), so it cannot summarize
big graphs reliably. ``summarize_graph`` computes structure deterministically and returns a COMPACT
digest the agent can reason over, plus a ``node -> cluster`` map used to recolor the display.

Clustering is topology-aware because modularity-style community detection degenerates on highly
constrained graphs (a star / universal-hub graph -- e.g. one disease node every path ends at -- or a
near-bipartite layered graph). So:

- ``mesh``         (a real interaction network)  -> Louvain communities on the largest component with
                                                    universal hubs stripped, so peripheral modules show.
- ``hub_dominated`` (one/few nodes touch nearly all) -> NO modularity; facet by node category + edge
                                                    predicate, plus a k-core, and say so.
- ``fragmented``   (many components)              -> clusters are the connected components.

Pure networkx (no Streamlit/TCT) so it is unit-testable.
"""

from __future__ import annotations

from collections import Counter
from typing import Any, Dict, Iterable, List, Optional, Tuple

import networkx as nx

# Tuning knobs (module-level so tests and callers can reason about them).
HUB_DEGREE_FRACTION = 0.5     # a node touching >= this fraction of all others is a "universal hub"
LARGE_QUESTION_HINT = 150     # results at/above this are "large" (mirrors the agent prompt rule)
MAX_EDGES_FOR_CLUSTERING = 20000  # above this, cluster on a degree-thresholded subgraph (noted)
MAX_CLUSTERS_REPORTED = 12
TOP_MEMBERS_PER_CLUSTER = 6


def _simple_undirected(graph: nx.Graph) -> nx.Graph:
    """Collapse a (Multi)DiGraph to a simple undirected graph for structural analysis."""
    ug = nx.Graph()
    ug.add_nodes_from(graph.nodes(data=True))
    for u, v in graph.edges():
        if u != v:
            ug.add_edge(u, v)
    return ug


def _node_category(graph: nx.Graph, node: str) -> str:
    """Best-effort biolink category for a node (attribute first, else CURIE prefix)."""
    data = graph.nodes.get(node, {})
    cat = data.get("category") or data.get("node_category")
    if cat:
        return str(cat).replace("biolink:", "")
    prefix = str(node).split(":", 1)[0] if ":" in str(node) else ""
    return prefix or "Unknown"


def _label(node: str, curie_to_symbol: Optional[Dict[str, str]], graph: nx.Graph) -> str:
    if curie_to_symbol and node in curie_to_symbol:
        return curie_to_symbol[node]
    data = graph.nodes.get(node, {})
    # viz_normalizer stores the resolved name on ``label``; ``name`` is a defensive fallback. This is
    # what makes protein / chemical / concept members readable instead of bare CURIEs.
    name = data.get("label") or data.get("name")
    return str(name) if name else str(node)


def _predicate_facets(graph: nx.Graph) -> Dict[str, int]:
    counts: Counter = Counter()
    for _u, _v, data in graph.edges(data=True):
        pred = data.get("predicate")
        if pred:
            counts[str(pred).replace("biolink:", "")] += 1
    return dict(counts.most_common(12))


def _category_facets(graph: nx.Graph) -> Dict[str, int]:
    counts: Counter = Counter(_node_category(graph, n) for n in graph.nodes)
    return dict(counts.most_common(12))


def _cluster_dominant_predicates(graph: nx.Graph, member_set: set, top: int = 3) -> List[str]:
    """Top predicates on edges incident to a cluster's members -- says HOW the cluster connects."""
    counts: Counter = Counter()
    for u, v, data in graph.edges(data=True):
        if u in member_set or v in member_set:
            pred = data.get("predicate")
            if pred:
                counts[str(pred).replace("biolink:", "")] += 1
    return [p for p, _ in counts.most_common(top)]


def _cluster_record(
    ug: nx.Graph, members: Iterable[str], cluster_id: str,
    curie_to_symbol: Optional[Dict[str, str]], graph: nx.Graph,
) -> Dict[str, Any]:
    members = list(members)
    member_set = set(members)
    # Rank by OVERALL degree so each cluster surfaces its most-connected (most informative) named
    # members -- better than intra-cluster degree for category-facet clusters whose same-category
    # members barely interconnect (e.g. proteins/concepts that mostly link to genes).
    ranked = sorted(members, key=lambda n: ug.degree(n) if n in ug else 0, reverse=True)
    cats = Counter(_node_category(graph, n) for n in members)
    return {
        "id": cluster_id,
        "size": len(members),
        "top_members": [_label(n, curie_to_symbol, graph) for n in ranked[:TOP_MEMBERS_PER_CLUSTER]],
        "dominant_category": cats.most_common(1)[0][0] if cats else "Unknown",
        "dominant_predicates": _cluster_dominant_predicates(graph, member_set),
    }


def _universal_hubs(ug: nx.Graph, anchors: set) -> set:
    """Nodes that connect to (nearly) everything: explicit anchors (query genes / disease) plus any
    node whose degree exceeds HUB_DEGREE_FRACTION of all other nodes."""
    n = ug.number_of_nodes()
    hubs = {a for a in anchors if a in ug}
    if n > 2:
        thresh = HUB_DEGREE_FRACTION * (n - 1)
        hubs |= {node for node, deg in ug.degree() if deg >= thresh}
    return hubs


def _leaf_fraction(ug: nx.Graph, hubs: set) -> float:
    """Fraction of non-hub nodes that are leaves attached only to hubs (star-ness indicator)."""
    periphery = [n for n in ug.nodes if n not in hubs]
    if not periphery:
        return 0.0
    star_leaves = 0
    for n in periphery:
        nbrs = set(ug.neighbors(n))
        if nbrs and nbrs <= hubs:  # every neighbor is a hub
            star_leaves += 1
    return star_leaves / len(periphery)


def _classify_topology(ug: nx.Graph, hubs: set) -> str:
    components = list(nx.connected_components(ug))
    n = ug.number_of_nodes()
    if n == 0:
        return "fragmented"
    largest = max((len(c) for c in components), default=0)
    # Many components and no single dominant one -> fragmented (cluster = components).
    if len(components) > 1 and largest < 0.6 * n:
        return "fragmented"
    # Dominant hub(s) with most of the periphery hanging off them -> hub_dominated (star/bipartite).
    if hubs and _leaf_fraction(ug, hubs) >= 0.5:
        return "hub_dominated"
    return "mesh"


def _louvain(ug: nx.Graph, seed: int) -> List[set]:
    try:
        return list(nx.community.louvain_communities(ug, seed=seed))
    except Exception:  # noqa: BLE001 - fall back to a deterministic modularity method
        return [set(c) for c in nx.community.greedy_modularity_communities(ug)]


def summarize_graph(
    graph: nx.Graph,
    query_genes: Optional[Iterable[str]] = None,
    disease_curie: Optional[str] = None,
    curie_to_symbol: Optional[Dict[str, str]] = None,
    seed: int = 0,
) -> Tuple[Dict[str, Any], Dict[str, str]]:
    """Return ``(digest, node_cluster)``.

    ``digest`` is a compact, agent-facing structured summary; ``node_cluster`` maps each node to its
    cluster id (for recoloring the display)."""
    anchors = set(query_genes or [])
    if disease_curie:
        anchors.add(disease_curie)

    n_nodes = graph.number_of_nodes()
    n_edges = graph.number_of_edges()
    notes: List[str] = []

    ug = _simple_undirected(graph)

    # Guard very large graphs: cluster on the highest-degree core (noted, not silent).
    if ug.number_of_edges() > MAX_EDGES_FOR_CLUSTERING:
        keep = {n for n, _ in sorted(ug.degree(), key=lambda kv: kv[1], reverse=True)[:4000]} | (
            anchors & set(ug.nodes)
        )
        ug = ug.subgraph(keep).copy()
        notes.append(
            f"graph has >{MAX_EDGES_FOR_CLUSTERING} edges; clustered on the {ug.number_of_nodes()} "
            "highest-degree nodes (structure is approximate)."
        )

    hubs = _universal_hubs(ug, anchors)
    topology = _classify_topology(ug, hubs)

    components = sorted(
        (len(c) for c in nx.connected_components(ug)), reverse=True
    )
    node_cluster: Dict[str, str] = {}
    clusters: List[Dict[str, Any]] = []

    if topology == "mesh":
        method = "louvain_communities (universal hubs stripped)"
        core = ug.subgraph([n for n in ug.nodes if n not in hubs])
        comms = [c for c in _louvain(core, seed) if c]
        comms.sort(key=len, reverse=True)
        for i, members in enumerate(comms[:MAX_CLUSTERS_REPORTED]):
            cid = f"C{i + 1}"
            for m in members:
                node_cluster[m] = cid
            clusters.append(_cluster_record(ug, members, cid, curie_to_symbol, graph))
        for h in hubs:  # hubs form their own labelled group so they still get a color
            node_cluster[h] = "hubs"
        if hubs:
            clusters.append(_cluster_record(ug, hubs, "hubs", curie_to_symbol, graph))
        if len(comms) > MAX_CLUSTERS_REPORTED:
            notes.append(f"{len(comms)} communities found; reporting the {MAX_CLUSTERS_REPORTED} largest.")

    elif topology == "hub_dominated":
        method = "category/predicate facets (hub-dominated; modularity not meaningful)"
        notes.append(
            "network is hub-dominated (star/bipartite); grouped by node category and edge predicate "
            "rather than community."
        )
        # k-core surfaces any denser sub-core beyond the trivial star.
        try:
            core_k = nx.k_core(ug) if ug.number_of_edges() else nx.Graph()
            if core_k.number_of_nodes() and core_k.number_of_nodes() < ug.number_of_nodes():
                notes.append(f"k-core ({core_k.number_of_nodes()} nodes) is the densest sub-core.")
        except Exception:  # noqa: BLE001
            pass
        # clusters = node categories (so recolor groups by category)
        by_cat: Dict[str, List[str]] = {}
        for node in ug.nodes:
            by_cat.setdefault(_node_category(ug, node), []).append(node)
        for cat, members in sorted(by_cat.items(), key=lambda kv: len(kv[1]), reverse=True)[
            :MAX_CLUSTERS_REPORTED
        ]:
            for m in members:
                node_cluster[m] = cat
            clusters.append(_cluster_record(ug, members, cat, curie_to_symbol, graph))

    else:  # fragmented
        method = "connected components"
        comps = sorted(nx.connected_components(ug), key=len, reverse=True)
        for i, members in enumerate(comps[:MAX_CLUSTERS_REPORTED]):
            cid = f"component {i + 1}"
            for m in members:
                node_cluster[m] = cid
            clusters.append(_cluster_record(ug, members, cid, curie_to_symbol, graph))
        if len(comps) > MAX_CLUSTERS_REPORTED:
            notes.append(f"{len(comps)} components; reporting the {MAX_CLUSTERS_REPORTED} largest.")

    digest = {
        "n_nodes": n_nodes,
        "n_edges": n_edges,
        "topology": topology,
        "method": method,
        "n_clusters": len(clusters),
        "components": components[:10],
        "clusters": clusters,
        "category_facets": _category_facets(ug),
        "predicate_facets": _predicate_facets(graph),
        "notes": notes,
    }
    return digest, node_cluster
