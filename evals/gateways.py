"""TctGateway stand-ins for evals.

- ``RecordedGateway`` returns deterministic, synthetic-but-realistic results offline (no network),
  including a ``knowledge_graph`` that supports ``.to_networkx(...)`` so the expression and display
  tools (which call ``tct_result_to_nx``) actually run during evals.
- ``CassetteGateway`` replays REAL TCT results recorded once to JSON (see ``record_cassettes.py``),
  so answer-quality / faithfulness can be judged against realistic evidence, deterministically.

Both return the same result shape via ``RecordedKnowledgeGraph``; only the edge data differs
(synthetic vs recorded-from-live).
"""

from __future__ import annotations

import json
import types
from pathlib import Path
from typing import Any, Dict, List, Optional

import networkx as nx
import pandas as pd

CASSETTE_DIR = Path(__file__).parent / "cassettes"


class RecordedKnowledgeGraph:
    """Dict-like KnowledgeGraph stand-in: ``.items()`` for edge_evidence/degree, ``.to_networkx()``
    for rendering/display tools. ``edges`` maps edge_id -> a TRAPI-ish edge dict."""

    def __init__(self, edges: Dict[str, dict]):
        self._edges = dict(edges)

    def __len__(self) -> int:
        return len(self._edges)

    def items(self):
        return list(self._edges.items())

    def to_networkx(self, resolve_names: bool = False, include_attributes: bool = False) -> nx.MultiDiGraph:
        g = nx.MultiDiGraph()
        for edge_id, e in self._edges.items():
            s, o = e["subject"], e["object"]
            g.add_node(s, label=e.get("subject_name", s))
            g.add_node(o, label=e.get("object_name", o))
            g.add_edge(
                s,
                o,
                key=edge_id,
                predicate=e.get("predicate", "biolink:related_to"),
                primary_sources=e.get("primary_sources", []),
                aggregator_sources=e.get("aggregator_sources", []),
                publications=e.get("publications", []),
                supporting_text=e.get("supporting_text", []),
                confidence_scores=e.get("confidence_scores", {}),
            )
        return g


def _neighborhood_result(gene_curie: str, drug_edges: List[dict]):
    edges = {f"e{i}": e for i, e in enumerate(drug_edges)}
    ranked = pd.DataFrame(
        [
            {
                "output_node": e["object"],
                "Name": e.get("object_name", e["object"]),
                "Num_of_primary_infores": len(e.get("primary_sources", [])) or 1,
            }
            for e in drug_edges
        ]
    )
    return types.SimpleNamespace(
        input_node_id=gene_curie, ranked=ranked, knowledge_graph=RecordedKnowledgeGraph(edges)
    )


def _path_result(node1: str, node2: str, kg1_edges: List[dict], kg2_edges: List[dict], paths: List[dict]):
    return types.SimpleNamespace(
        node1_id=node1,
        node2_id=node2,
        paths=pd.DataFrame(paths),
        knowledge_graph1=RecordedKnowledgeGraph({f"a{i}": e for i, e in enumerate(kg1_edges)}),
        knowledge_graph2=RecordedKnowledgeGraph({f"b{i}": e for i, e in enumerate(kg2_edges)}),
    )


def _edge_evidence_from_graphs(result: Any, subject: str, object: str) -> Dict[str, Any]:
    """Look up the matching edge(s) across a result's knowledge graph(s)."""
    pubs, texts, scores, preds, matched = [], [], {}, [], 0
    kgs = [getattr(result, a, None) for a in ("knowledge_graph", "knowledge_graph1", "knowledge_graph2")]
    if not any(kgs) and hasattr(result, "items"):
        kgs = [result]
    for kg in [k for k in kgs if k is not None]:
        for _eid, e in kg.items():
            if isinstance(e, dict) and e.get("subject") == subject and e.get("object") == object:
                matched += 1
                preds.append(e.get("predicate"))
                pubs.extend(e.get("publications", []))
                texts.extend(e.get("supporting_text", []))
                scores.update(e.get("confidence_scores", {}))
    return {
        "subject": subject,
        "object": object,
        "edges_matched": matched,
        "predicates": sorted(set(filter(None, preds))),
        "publications": sorted(set(pubs)),
        "supporting_text": texts[:20],
        "confidence_scores": scores,
    }


# Synthetic-but-plausible drugs reused for any neighborhood (offline tier).
_SYNTH_DRUGS = [
    ("PUBCHEM.COMPOUND:9829523", "midostaurin", ["infores:drugbank"], ["PMID:28644114"]),
    ("PUBCHEM.COMPOUND:49803313", "gilteritinib", ["infores:drugcentral"], ["PMID:31665578"]),
    ("CHEBI:63631", "sorafenib", ["infores:chembl"], ["PMID:18950845"]),
]


class RecordedGateway:
    """Deterministic, offline stand-in for TctGateway (synthetic but to_networkx-capable)."""

    def resolve_genes(self, names, only_taxa=None, biolink_type=None):
        return {
            name: {"curie": f"NCBIGene:{1000 + i}", "types": ["biolink:Gene"]}
            for i, name in enumerate(names)
        }

    def neighborhood(self, gene_curie, target_categories=None):
        edges = [
            {
                "subject": gene_curie,
                "object": curie,
                "predicate": "biolink:affects",
                "subject_name": gene_curie,
                "object_name": name,
                "primary_sources": srcs,
                "sources": [{"resource_id": s, "resource_role": "primary_knowledge_source"} for s in srcs]
                + [{"resource_id": "infores:automat-robokop", "resource_role": "aggregator_knowledge_source"}],
                "publications": pubs,
                "supporting_text": [f"{name} affects the gene product."],
                "confidence_scores": {"tmkp_confidence_score": 0.8},
            }
            for curie, name, srcs, pubs in _SYNTH_DRUGS
        ]
        return _neighborhood_result(gene_curie, edges)

    def path(self, node1_curie, node2_curie, intermediate_categories):
        inter = "NCBIGene:5290"
        kg1 = [{"subject": node1_curie, "object": inter, "predicate": "biolink:related_to",
                "object_name": "PIK3CA", "publications": ["PMID:20176728"]}]
        kg2 = [{"subject": inter, "object": node2_curie, "predicate": "biolink:related_to",
                "object_name": node2_curie, "publications": ["PMID:25409260"]}]
        paths = [{"n1": node1_curie, "intermediate": inter, "n2": node2_curie}]
        return _path_result(node1_curie, node2_curie, kg1, kg2, paths)

    def gene_network(self, gene_curies):
        edges = {
            f"e{i}": {
                "subject": gene_curies[i],
                "object": gene_curies[i + 1],
                "predicate": "biolink:interacts_with",
                "publications": ["PMID:20705237"],
            }
            for i in range(len(gene_curies) - 1)
        }
        return RecordedKnowledgeGraph(edges)

    def edge_evidence(self, result, subject, object):
        ev = _edge_evidence_from_graphs(result, subject, object)
        if ev["edges_matched"] == 0:  # offline fallback so the tool always has something
            ev.update({"edges_matched": 1, "publications": ["PMID:12345"],
                       "supporting_text": ["Experimental evidence of interaction."],
                       "confidence_scores": {"tmkp_confidence_score": 0.87}})
        return ev


class CassetteGateway:
    """Replays REAL TCT results recorded to JSON by record_cassettes.py (deterministic, real data)."""

    def __init__(self, cassette_dir: Path = CASSETTE_DIR):
        path = Path(cassette_dir) / "cassettes.json"
        if not path.exists():
            raise FileNotFoundError(
                f"No cassettes at {path}. Record them first: python -m evals.record_cassettes"
            )
        self._data = json.loads(path.read_text())

    def _get(self, kind: str, key: str) -> Optional[dict]:
        return self._data.get(kind, {}).get(key)

    def resolve_genes(self, names, only_taxa=None, biolink_type=None):
        rec = self._data.get("resolve", {})
        return {n: rec.get(n) for n in names}

    def neighborhood(self, gene_curie, target_categories=None):
        rec = self._get("neighborhood", gene_curie) or {"edges": [], "ranked": []}
        return types.SimpleNamespace(
            input_node_id=gene_curie,
            ranked=pd.DataFrame(rec["ranked"]),
            knowledge_graph=RecordedKnowledgeGraph({f"e{i}": e for i, e in enumerate(rec["edges"])}),
        )

    def path(self, node1_curie, node2_curie, intermediate_categories):
        # Recorded under one ordering; the agent may pass the two CURIEs in either order.
        rec = (
            self._get("path", f"{node1_curie}|{node2_curie}")
            or self._get("path", f"{node2_curie}|{node1_curie}")
            or {"kg1": [], "kg2": [], "paths": []}
        )
        return _path_result(node1_curie, node2_curie, rec["kg1"], rec["kg2"], rec["paths"])

    def gene_network(self, gene_curies):
        rec = self._get("gene_network", "|".join(sorted(gene_curies))) or {"edges": []}
        return RecordedKnowledgeGraph({f"e{i}": e for i, e in enumerate(rec["edges"])})

    def edge_evidence(self, result, subject, object):
        return _edge_evidence_from_graphs(result, subject, object)
