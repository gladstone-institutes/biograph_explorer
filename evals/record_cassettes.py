"""One-time: record REAL TCT results for the eval gene sets to evals/cassettes/cassettes.json.

Hits live TCT (slow, ~5-10 min) for the finder calls the eval cases need, so the CassetteGateway
can replay real data deterministically (and free) for faithfulness-on-real-data sweeps:

    python -m evals.record_cassettes
    python -m evals.runner --gateway cassette --model claude-haiku-4-5
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List

from dotenv import load_dotenv

from .cases import GENESET
from .gateways import CASSETTE_DIR

EXTRA_ENTITIES = ["venetoclax"]  # needed by path_bcl2_venetoclax


def _kg_to_edges(kg, name_map: Dict[str, str], extract) -> List[dict]:
    out: List[dict] = []
    try:
        items = list(kg.items())
    except Exception:  # noqa: BLE001
        return out
    for _eid, e in items[:200]:
        if not isinstance(e, dict):
            continue
        s, o = e.get("subject"), e.get("object")
        rich = extract(e.get("attributes", [])) if e.get("attributes") else {}
        out.append(
            {
                "subject": s,
                "object": o,
                "predicate": e.get("predicate", "biolink:related_to"),
                "subject_name": name_map.get(s, s),
                "object_name": name_map.get(o, o),
                "primary_sources": [],
                "sources": e.get("sources", []),  # captured for the data_sources tool on re-record
                "publications": (rich.get("publications") or [])[:5],
                "supporting_text": (rich.get("supporting_text") or [])[:2],
                "confidence_scores": rich.get("confidence_scores") or {},
            }
        )
    return out


def main() -> int:
    load_dotenv()
    from TCT import name_resolver
    from TCT.attribute_extraction import extract_rich_edge_attributes

    from geneset_translator.agent.tct_adapter import TctGateway, load_resources

    gw = TctGateway(load_resources())
    data: Dict[str, dict] = {"resolve": {}, "neighborhood": {}, "path": {}, "gene_network": {}}

    # 1. resolve
    names = GENESET + EXTRA_ENTITIES
    resolved = gw.resolve_genes(names)
    data["resolve"] = {n: resolved.get(n) for n in names}
    gene_curies = [resolved[g]["curie"] for g in GENESET if resolved.get(g)]
    print("resolved:", {n: (v and v["curie"]) for n, v in resolved.items()})

    def _names(curies):
        try:
            info = name_resolver.batch_lookup(strings=list(set(curies)))
            return {c: getattr(info.get(c), "label", c) or c for c in curies}
        except Exception:  # noqa: BLE001
            return {}

    # 2. neighborhoods (one per gene)
    for g in GENESET:
        if not resolved.get(g):
            continue
        curie = resolved[g]["curie"]
        print(f"neighborhood {g} ({curie})...", flush=True)
        nb = gw.neighborhood(curie)
        nm = _names([curie] + [str(x) for x in getattr(nb, "ranked").get("output_node", [])])
        data["neighborhood"][curie] = {
            "edges": _kg_to_edges(nb.knowledge_graph, nm, extract_rich_edge_attributes),
            "ranked": nb.ranked.head(25).to_dict("records"),
        }

    # 3. path BCL2 <-> venetoclax
    if resolved.get("BCL2") and resolved.get("venetoclax"):
        n1, n2 = resolved["BCL2"]["curie"], resolved["venetoclax"]["curie"]
        print(f"path {n1} <-> {n2}...", flush=True)
        pr = gw.path(n1, n2, ["biolink:Gene", "biolink:Protein", "biolink:ChemicalEntity"])
        data["path"][f"{n1}|{n2}"] = {
            "kg1": _kg_to_edges(getattr(pr, "knowledge_graph1", None) or _Empty(), {}, extract_rich_edge_attributes),
            "kg2": _kg_to_edges(getattr(pr, "knowledge_graph2", None) or _Empty(), {}, extract_rich_edge_attributes),
            "paths": pr.paths.head(25).to_dict("records"),
        }

    # 4. gene_network among all genes
    print("gene_network...", flush=True)
    kg = gw.gene_network(gene_curies)
    data["gene_network"]["|".join(sorted(gene_curies))] = {
        "edges": _kg_to_edges(kg, _names(gene_curies), extract_rich_edge_attributes)
    }

    CASSETTE_DIR.mkdir(parents=True, exist_ok=True)
    out = CASSETTE_DIR / "cassettes.json"
    out.write_text(json.dumps(data, default=str, indent=2))
    print(f"\nWrote {out}")
    return 0


class _Empty:
    def items(self):
        return []


if __name__ == "__main__":
    raise SystemExit(main())
