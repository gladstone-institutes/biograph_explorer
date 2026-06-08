"""Emit a standalone Python script that reproduces the TCT results an agent collected.

Input is the session tool log (list of ``{"name", "args"}`` from successful tool calls). Output
is a runnable script using ``TranslatorResources.load()`` and the TCT finders, so a user can
reproduce (and audit) exactly what the agent ran.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

HUMAN_TAXON = "NCBITaxon:9606"
_DEFAULT_TARGETS = ["biolink:Drug", "biolink:SmallMolecule", "biolink:ChemicalEntity"]
_DEFAULT_INTERMEDIATES = ["biolink:Gene", "biolink:Protein", "biolink:ChemicalEntity"]

_HEADER = '''#!/usr/bin/env python
"""Reproduce the TCT results collected in a GeneSet Translator chat session.

Auto-generated. Install TCT first, e.g.:
  pip install "TCT @ git+https://github.com/natalie-23-gill/Translator_component_toolkit.git@f702893"
'''

_IMPORTS = '''from TCT import TCT, name_resolver, translator_query
from TCT.translator_resources import TranslatorResources

resources = TranslatorResources.load()
'''


def _py(value: Any) -> str:
    return repr(value)


def generate_reproduction_script(
    tool_log: List[Dict[str, Any]],
    questions: Optional[List[str]] = None,
) -> str:
    """Return a standalone .py string reproducing the finder calls in ``tool_log``."""
    lines: List[str] = [_HEADER]
    if questions:
        lines.append("Questions asked:")
        for q in questions:
            lines.append(f"  - {q}")
    lines.append('"""')
    lines.append("")
    lines.append(_IMPORTS)

    counters = {"nb": 0, "pr": 0, "kg": 0}

    for call in tool_log:
        name = call.get("name")
        args = call.get("args", {}) or {}

        if name == "resolve_genes":
            names = args.get("names", [])
            lines.append(f"# resolve_genes: {names}")
            lines.append(
                f"resolved = name_resolver.batch_lookup(strings={_py(names)}, only_taxa={_py(HUMAN_TAXON)})"
            )
            lines.append("print({k: getattr(v, 'curie', None) for k, v in resolved.items()})")
            lines.append("")

        elif name == "gene_neighborhood":
            counters["nb"] += 1
            var = f"nb{counters['nb']}"
            targets = args.get("target_categories") or _DEFAULT_TARGETS
            lines.append(f"# gene_neighborhood: {args.get('gene_curie')}")
            lines.append(
                f"{var} = TCT.Neighborhood_finder(input_node={_py(args.get('gene_curie'))}, "
                f"node2_categories={_py(targets)}, resources=resources)"
            )
            lines.append(f"print({var}.ranked.head(25))")
            lines.append("")

        elif name == "path_between":
            counters["pr"] += 1
            var = f"pr{counters['pr']}"
            intermediates = args.get("intermediate_categories") or _DEFAULT_INTERMEDIATES
            lines.append(
                f"# path_between: {args.get('node1_curie')} <-> {args.get('node2_curie')}"
            )
            lines.append(
                f"{var} = TCT.Path_finder(input_node1={_py(args.get('node1_curie'))}, "
                f"input_node2={_py(args.get('node2_curie'))}, "
                f"intermediate_categories={_py(intermediates)}, resources=resources)"
            )
            lines.append(f"print({var}.paths.head(25))")
            lines.append("")

        elif name == "gene_network":
            counters["kg"] += 1
            var = f"kg{counters['kg']}"
            curies = args.get("gene_curies", [])
            lines.append(f"# gene_network among: {curies}")
            lines.append(f"_genes = {_py(curies)}")
            lines.append(
                "_predicates = list(set(TCT.select_concept("
                "sub_list=['biolink:Gene'], obj_list=['biolink:Gene'], metaKG=resources.meta_kg)))"
            )
            lines.append(
                "_apis = list(TCT.select_API("
                "sub_list=['biolink:Gene'], obj_list=['biolink:Gene'], metaKG=resources.meta_kg))"
            )
            lines.append(
                "_query = TCT.format_query_json(_genes, _genes, ['biolink:Gene'], "
                "['biolink:Gene'], _predicates)"
            )
            lines.append(
                f"{var} = translator_query.parallel_api_query(query_json=_query, "
                f"select_APIs=_apis, resources=resources, max_workers=max(1, len(_apis)))"
            )
            lines.append(f"print('edges:', len({var}))")
            lines.append("")

        elif name in (
            "edge_evidence",
            "cell_type_expression",
            "node_metadata",
            "data_sources",
            "filter_graph",
            "add_disease_node",
            "show_result",
        ):
            # Not a TCT finder call (operates on an already-collected result or external API);
            # note it for context.
            lines.append(f"# (skipped non-finder step: {name} {args})")
            lines.append("")

    return "\n".join(lines)
