"""The single boundary between the agent and TCT (Translator Component Toolkit).

Everything TCT-specific for the agent path lives here: loading resources, the finder calls,
gene resolution, edge-evidence extraction, and converting a TCT result into a NetworkX graph
for rendering. If the app later moves from this git branch to a published PyPI TCT, this is
the one file to change.

Also defines the ``ResultStash`` (full result objects kept server-side, keyed by ``result_id``)
so large graphs never go into the LLM context.
"""

from __future__ import annotations

import logging
import threading
import time
from typing import Any, Dict, List, Optional, Protocol, Union, runtime_checkable

import networkx as nx

from . import viz_normalizer

try:
    from TCT import TCT, name_resolver, translator_query
    from TCT.translator_resources import TranslatorResources
    from TCT.attribute_extraction import extract_rich_edge_attributes

    TCT_AVAILABLE = True
except ImportError:  # pragma: no cover - TCT optional at import time
    TCT_AVAILABLE = False

logger = logging.getLogger(__name__)

# Default to human so gene symbols resolve to the canonical human CURIE. Verified: plain
# batch_lookup returns non-human orthologs (e.g. FLT3 -> NCBIGene:100515445, 0 edges),
# while only_taxa="NCBITaxon:9606" returns FLT3 -> NCBIGene:2322.
HUMAN_TAXON = "NCBITaxon:9606"

DEFAULT_TARGET_CATEGORIES = [
    "biolink:Drug",
    "biolink:SmallMolecule",
    "biolink:ChemicalEntity",
]


def load_resources():
    """Build the (expensive) TranslatorResources once. Wrap with @st.cache_resource at the UI layer."""
    if not TCT_AVAILABLE:
        raise ImportError("TCT library not installed.")
    return TranslatorResources.load()


# --------------------------------------------------------------------------------------
# Result stash
# --------------------------------------------------------------------------------------
@runtime_checkable
class ResultStash(Protocol):
    def put(self, obj: Any) -> str: ...
    def get(self, result_id: str) -> Any: ...


class DictResultStash:
    """ResultStash backed by a plain dict (the chat page passes ``st.session_state``'s dict).

    Thread-safe: tool calls in one agent step run concurrently, so ``put`` must lock to avoid
    two threads computing the same ``result_id`` from ``len(store)``.
    """

    def __init__(self, store: Optional[Dict[str, Any]] = None) -> None:
        self._store: Dict[str, Any] = store if store is not None else {}
        self._lock = threading.Lock()

    def put(self, obj: Any) -> str:
        with self._lock:
            result_id = f"res_{len(self._store) + 1}"
            self._store[result_id] = obj
            return result_id

    def get(self, result_id: str) -> Any:
        return self._store.get(result_id)


# --------------------------------------------------------------------------------------
# Gateway
# --------------------------------------------------------------------------------------
class TctGateway:
    """Thin, injectable wrapper over the TCT result-class API.

    Methods return raw TCT objects / DataFrames; the tool layer compacts and stashes them.
    """

    def __init__(self, resources: Any) -> None:
        if not TCT_AVAILABLE:
            raise ImportError("TCT library not installed.")
        self.resources = resources

    # -- resolution --------------------------------------------------------------------
    def resolve_genes(
        self,
        names: List[str],
        only_taxa: Optional[str] = HUMAN_TAXON,
        biolink_type: Optional[str] = None,
    ) -> Dict[str, Optional[Dict[str, Any]]]:
        """Resolve names to CURIEs.

        For GENES (the default, ``biolink_type`` unset): two-pass when ``only_taxa`` is set --
        pass 1 looks up with the taxon filter so gene symbols pick the canonical human CURIE;
        pass 2 retries any unresolved name WITHOUT the taxon filter so chemicals/drugs (which
        carry no taxon) still resolve (e.g. ['BCL2', 'venetoclax'] -> NCBIGene:596 + a chemical).

        For NON-GENE entities (``biolink_type`` like 'biolink:Disease' / 'biolink:Drug'): drop the
        taxon filter entirely and constrain by type. This is required because the human-taxon
        filter mis-maps a disease/drug NAME onto a same-named gene -- e.g. with the taxon filter
        'acute myeloid leukemia' resolves to NCBIGene:861 (RUNX1), but without it (and typed
        Disease) it correctly resolves to MONDO:0018874.

        Returns ``{name: {curie, types, label}}`` (or ``{name: None}`` when unresolved).
        """
        is_gene_type = (
            biolink_type is None
            or "gene" in biolink_type.lower()
            or "protein" in biolink_type.lower()
        )
        if not is_gene_type:
            return self._lookup(names, only_taxa=None, biolink_type=biolink_type)

        out = self._lookup(names, only_taxa, biolink_type=biolink_type)
        unresolved = [n for n in names if out.get(n) is None]
        if unresolved and only_taxa:
            logger.info(
                "TCT resolve_genes: retrying %d unresolved name(s) without taxon filter", len(unresolved)
            )
            out.update(self._lookup(unresolved, None, biolink_type=biolink_type))
        return out

    def _lookup(
        self,
        names: List[str],
        only_taxa: Optional[str],
        biolink_type: Optional[str] = None,
    ) -> Dict[str, Optional[Dict[str, Any]]]:
        kwargs: Dict[str, Any] = {}
        if only_taxa:
            kwargs["only_taxa"] = only_taxa
        if biolink_type:
            kwargs["biolink_type"] = biolink_type
        logger.info(
            "TCT resolve_genes: %d name(s) (taxon=%s type=%s)", len(names), only_taxa, biolink_type
        )
        t0 = time.time()
        info = name_resolver.batch_lookup(strings=list(names), **kwargs)
        logger.info("TCT resolve_genes: done in %.1fs", time.time() - t0)
        out: Dict[str, Optional[Dict[str, Any]]] = {}
        for name in names:
            node = info.get(name)
            if node is None or not getattr(node, "curie", None):
                out[name] = None
                continue
            types = getattr(node, "categories", None) or getattr(node, "types", None) or []
            out[name] = {
                "curie": node.curie,
                "types": list(types)[:6],
                "label": getattr(node, "label", None),
            }
        return out

    # -- finders -----------------------------------------------------------------------
    def neighborhood(self, gene_curie: str, target_categories: Optional[List[str]] = None):
        cats = target_categories or list(DEFAULT_TARGET_CATEGORIES)
        logger.info("TCT Neighborhood_finder: %s -> %s", gene_curie, cats)
        t0 = time.time()
        try:
            nb = TCT.Neighborhood_finder(
                input_node=gene_curie, node2_categories=cats, resources=self.resources
            )
        except AttributeError as e:
            # TCT.Neighborhood_finder dereferences node_normalizer.get_normalized_nodes(...).types
            # without a None-guard; an unnormalizable input node surfaces as an opaque
            # AttributeError. Convert it to an actionable message the agent can recover from.
            raise ValueError(
                f"could not run a neighborhood for {gene_curie!r}: the Translator node normalizer "
                f"did not recognize it. Re-resolve the entity to a standard CURIE (for a disease or "
                f"drug, resolve by name with the matching biolink_type) and try again."
            ) from e
        logger.info("TCT Neighborhood_finder: done in %.1fs", time.time() - t0)
        return nb

    def path(self, node1_curie: str, node2_curie: str, intermediate_categories: List[str]):
        logger.info("TCT Path_finder: %s <-> %s via %s", node1_curie, node2_curie, intermediate_categories)
        t0 = time.time()
        pr = TCT.Path_finder(
            input_node1=node1_curie,
            input_node2=node2_curie,
            intermediate_categories=intermediate_categories,
            resources=self.resources,
        )
        logger.info("TCT Path_finder: done in %.1fs", time.time() - t0)
        return pr

    def gene_network(self, gene_curies: List[str]):
        sub = ["biolink:Gene"]
        obj = ["biolink:Gene"]
        sele_predicates = list(
            set(TCT.select_concept(sub_list=sub, obj_list=obj, metaKG=self.resources.meta_kg))
        )
        sele_apis = list(TCT.select_API(sub_list=sub, obj_list=obj, metaKG=self.resources.meta_kg))
        logger.info(
            "TCT gene_network: %d gene(s) across %d API(s), %d predicate(s)",
            len(gene_curies), len(sele_apis), len(sele_predicates),
        )
        query_json = TCT.format_query_json(
            list(gene_curies), list(gene_curies), sub, obj, sele_predicates
        )
        t0 = time.time()
        kg = translator_query.parallel_api_query(
            query_json=query_json,
            select_APIs=sele_apis,
            resources=self.resources,
            max_workers=max(1, len(sele_apis)),
        )
        logger.info("TCT gene_network: done in %.1fs", time.time() - t0)
        return kg

    # -- evidence ----------------------------------------------------------------------
    def edge_evidence(self, result_obj: Any, subject: str, object: str) -> Dict[str, Any]:
        """Aggregate publications/supporting text/scores for edges between two nodes."""
        pubs: List[str] = []
        texts: List[str] = []
        scores: Dict[str, Any] = {}
        predicates: List[str] = []
        matched = 0
        for kg in _knowledge_graphs(result_obj):
            for _edge_id, edge in kg.items():
                if not isinstance(edge, dict):
                    continue
                if edge.get("subject") == subject and edge.get("object") == object:
                    matched += 1
                    if edge.get("predicate"):
                        predicates.append(edge["predicate"])
                    rich = extract_rich_edge_attributes(edge.get("attributes", []))
                    pubs.extend(rich.get("publications", []) or [])
                    texts.extend(rich.get("supporting_text", []) or [])
                    cs = rich.get("confidence_scores", {}) or {}
                    if isinstance(cs, dict):
                        scores.update(cs)
        return {
            "subject": subject,
            "object": object,
            "edges_matched": matched,
            "predicates": sorted(set(predicates)),
            "publications": sorted(set(filter(None, pubs))),
            "supporting_text": [t for t in texts if t][:20],
            "confidence_scores": scores,
        }


# --------------------------------------------------------------------------------------
# TCT result -> NetworkX (for rendering via ui.network_viz)
# --------------------------------------------------------------------------------------
def _knowledge_graphs(result: Any) -> List[Any]:
    """Collect the dict-like KnowledgeGraph(s) attached to a TCT result (or the result itself)."""
    kgs: List[Any] = []
    for attr in ("knowledge_graph", "knowledge_graph1", "knowledge_graph2"):
        kg = getattr(result, attr, None)
        if kg is not None:
            kgs.append(kg)
    if not kgs and hasattr(result, "items") and hasattr(result, "to_networkx"):
        kgs.append(result)  # the result is itself a KnowledgeGraph
    return kgs


def _result_to_tct_graph(result: Any) -> nx.MultiDiGraph:
    """Build a TCT MultiDiGraph (with attributes) from any finder result or KnowledgeGraph."""
    kg = getattr(result, "knowledge_graph", None)
    if kg is not None:
        return kg.to_networkx(resolve_names=True, include_attributes=True)

    if getattr(result, "knowledge_graph1", None) is not None:
        graphs = [
            getattr(result, a).to_networkx(resolve_names=True, include_attributes=True)
            for a in ("knowledge_graph1", "knowledge_graph2")
            if getattr(result, a, None) is not None
        ]
        return nx.compose_all(graphs) if graphs else nx.MultiDiGraph()

    if hasattr(result, "items") and hasattr(result, "to_networkx"):
        return result.to_networkx(resolve_names=True, include_attributes=True)

    if hasattr(result, "to_networkx"):
        return result.to_networkx(resolve_names=True)

    raise ValueError(f"Unsupported TCT result type for graph conversion: {type(result)!r}")


def tct_result_to_nx(
    result: Any,
    query_gene_curies: List[str],
    curie_to_symbol: Optional[Dict[str, str]] = None,
    disease_curie: Optional[str] = None,
) -> nx.MultiDiGraph:
    """Convert a TCT finder result into a render-ready graph for ``ui.network_viz``."""
    tct_graph = _result_to_tct_graph(result)
    return viz_normalizer.normalize(
        tct_graph, query_gene_curies, curie_to_symbol=curie_to_symbol, disease_curie=disease_curie
    )
