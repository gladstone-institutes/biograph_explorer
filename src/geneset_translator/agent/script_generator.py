"""Emit a standalone Python script that reproduces the TCT graphs an agent collected.

Input is the session tool log (list of ``{"name", "args"}`` from successful tool calls). Output is
a runnable, TCT-only script (no ``geneset_translator`` dependency) that re-runs each finder, builds a
NetworkX graph WITH edge attributes (predicates, publications, knowledge sources) via TCT's
``to_networkx``, and exports each graph to Cytoscape.js JSON (``.cyjs``) for import into Cytoscape
Desktop. Identical finder calls are de-duplicated so each unique (slow) query runs once.
"""

from __future__ import annotations

import json
from typing import Any, Dict, List, Optional

HUMAN_TAXON = "NCBITaxon:9606"
_DEFAULT_TARGETS = ["biolink:Drug", "biolink:SmallMolecule", "biolink:ChemicalEntity"]
_DEFAULT_INTERMEDIATES = ["biolink:Gene", "biolink:Protein", "biolink:ChemicalEntity"]

_HEADER = '''#!/usr/bin/env python
"""Reproduce the TCT graphs collected in a GeneSet Translator chat session.

Auto-generated, TCT-only (no geneset_translator dependency). Install TCT first, e.g.:
  pip install "TCT @ git+https://github.com/natalie-23-gill/Translator_component_toolkit.git@f702893"

Each finder is re-run live against the NCATS Translator network: expect a ~30s one-time resource
load and roughly ~60s per neighborhood / path / gene-network call. Every result is converted to a
NetworkX graph (with edge predicates, publications, and knowledge sources) and written to
./tct_repro_output/*.cyjs -- in Cytoscape Desktop use File > Import > Network from File to open them.

These .cyjs files are the COMPLETE graph (all nodes and edges). The in-app view samples down to the
"Max graph nodes" setting (~200) for rendering speed; this export is not capped (see the node/edge
counts printed for each graph below).
'''

_PREAMBLE = '''import json
from pathlib import Path

import networkx as nx
from TCT import TCT, name_resolver, translator_query
from TCT.translator_resources import TranslatorResources

resources = TranslatorResources.load()
OUT = Path("tct_repro_output")
OUT.mkdir(exist_ok=True)


def to_graph(result):
    """TCT result / KnowledgeGraph -> nx.MultiDiGraph with edge attributes (predicates,
    publications, knowledge sources). TCT-only; mirrors how the app converts finder results."""
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
    return result.to_networkx(resolve_names=True)


def export_cyjs(graph, path):
    """Write a NetworkX graph to Cytoscape.js JSON, importable into Cytoscape Desktop."""
    data = nx.cytoscape_data(graph)
    Path(path).write_text(json.dumps(data, indent=2, default=str))
    print(f"  wrote {graph.number_of_nodes()} nodes / {graph.number_of_edges()} edges -> {path}")
'''

_GENE_NETWORK_HELPER = '''

def run_gene_network(genes):
    """Among-the-set gene-gene KnowledgeGraph (mirrors the app's gene_network tool)."""
    sub = obj = ["biolink:Gene"]
    predicates = list(set(TCT.select_concept(sub_list=sub, obj_list=obj, metaKG=resources.meta_kg)))
    apis = list(TCT.select_API(sub_list=sub, obj_list=obj, metaKG=resources.meta_kg))
    query = TCT.format_query_json(list(genes), list(genes), sub, obj, predicates)
    return translator_query.parallel_api_query(
        query_json=query, select_APIs=apis, resources=resources, max_workers=max(1, len(apis))
    )
'''


def _py(value: Any) -> str:
    return repr(value)


def _emit_resolve(args: Dict[str, Any]) -> List[str]:
    """Faithfully reproduce TctGateway.resolve_genes: non-gene types resolve taxon-free with a
    biolink_type filter; genes use the taxon -> taxon-free two-pass fallback."""
    names = args.get("names", [])
    biolink_type = args.get("biolink_type")
    is_gene = (
        biolink_type is None
        or "gene" in biolink_type.lower()
        or "protein" in biolink_type.lower()
    )
    out = [f"# resolve_genes: {names}" + (f" (biolink_type={biolink_type})" if biolink_type else "")]
    if not is_gene:
        out.append(
            f"resolved = name_resolver.batch_lookup(strings={_py(names)}, "
            f"biolink_type={_py(biolink_type)})"
        )
    else:
        bt = f", biolink_type={_py(biolink_type)}" if biolink_type else ""
        out.append(
            f"resolved = name_resolver.batch_lookup(strings={_py(names)}, "
            f"only_taxa={_py(HUMAN_TAXON)}{bt})"
        )
        out.append("_unresolved = [n for n in resolved if resolved.get(n) is None]")
        out.append("if _unresolved:  # pass 2: chemicals/drugs carry no taxon")
        out.append(f"    resolved.update(name_resolver.batch_lookup(strings=_unresolved{bt}))")
    out.append("print({k: getattr(v, 'curie', None) for k, v in resolved.items()})")
    out.append("")
    return out


def generate_reproduction_script(
    tool_log: List[Dict[str, Any]],
    questions: Optional[List[str]] = None,
) -> str:
    """Return a standalone .py string reproducing the finder calls in ``tool_log`` as NetworkX
    graphs exported to Cytoscape.js JSON. Identical finder calls are emitted once."""
    lines: List[str] = [_HEADER]
    if questions:
        lines.append("Questions asked:")
        for q in questions:
            lines.append(f"  - {q}")
    lines.append('"""')
    lines.append("")
    lines.append(_PREAMBLE)
    if any(c.get("name") == "gene_network" for c in tool_log):
        lines.append(_GENE_NETWORK_HELPER)

    counters = {"nb": 0, "pr": 0, "kg": 0}
    seen: Dict[str, str] = {}  # call signature -> variable, to de-duplicate repeated calls

    for call in tool_log:
        name = call.get("name")
        args = call.get("args", {}) or {}

        if name in ("resolve_genes", "gene_neighborhood", "path_between", "gene_network"):
            sig = name + "|" + json.dumps(args, sort_keys=True, default=str)
            if sig in seen:
                lines.append(f"# (duplicate of {seen[sig]}; skipped to avoid re-running)")
                lines.append("")
                continue
        else:
            sig = None

        if name == "resolve_genes":
            seen[sig] = "resolved"
            lines.extend(_emit_resolve(args))

        elif name == "gene_neighborhood":
            counters["nb"] += 1
            var = f"nb{counters['nb']}"
            seen[sig] = var
            curie = args.get("gene_curie")
            targets = args.get("target_categories") or _DEFAULT_TARGETS
            safe = str(curie).replace(":", "_")
            lines.append(f"# gene_neighborhood: {curie} -> {targets}")
            lines.append(
                f"{var} = TCT.Neighborhood_finder(input_node={_py(curie)}, "
                f"node2_categories={_py(targets)}, resources=resources)"
            )
            lines.append(f"print({var}.ranked.head(25))")
            lines.append(f"export_cyjs(to_graph({var}), OUT / {_py(f'{var}_{safe}.cyjs')})")
            lines.append("")

        elif name == "path_between":
            counters["pr"] += 1
            var = f"pr{counters['pr']}"
            seen[sig] = var
            n1, n2 = args.get("node1_curie"), args.get("node2_curie")
            intermediates = args.get("intermediate_categories") or _DEFAULT_INTERMEDIATES
            lines.append(f"# path_between: {n1} <-> {n2}")
            lines.append(
                f"{var} = TCT.Path_finder(input_node1={_py(n1)}, input_node2={_py(n2)}, "
                f"intermediate_categories={_py(intermediates)}, resources=resources)"
            )
            lines.append(f"print({var}.paths.head(25))")
            lines.append(f"export_cyjs(to_graph({var}), OUT / {_py(f'{var}.cyjs')})")
            lines.append("")

        elif name == "gene_network":
            counters["kg"] += 1
            var = f"kg{counters['kg']}"
            seen[sig] = var
            curies = args.get("gene_curies", [])
            lines.append(f"# gene_network among {len(curies)} genes: {curies}")
            lines.append(f"{var} = run_gene_network({_py(curies)})")
            lines.append(f"print('edges:', len({var}))")
            lines.append(f"export_cyjs(to_graph({var}), OUT / {_py(f'{var}_network.cyjs')})")
            lines.append("")

        elif name in ("edge_evidence", "cell_type_expression", "node_metadata", "data_sources"):
            # Reproducing these needs geneset_translator (NodeAnnotator / HPA / infores catalog) or
            # the already-collected result, so a TCT-only script can only note them.
            rid = args.get("result_id")
            ref = f" on {rid}" if rid else ""
            lines.append(
                f"# app-only step not reproduced in this TCT-only script "
                f"(needs geneset_translator): {name}{ref} {args}"
            )
            lines.append("")

        elif name in ("filter_graph", "add_disease_node", "show_result"):
            # Display curation edits the in-app graph object; not a standalone TCT operation.
            lines.append(
                f"# display-only step (edits the in-app graph, not reproduced): {name} {args}"
            )
            lines.append("")

    return "\n".join(lines)
