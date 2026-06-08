"""Offline smoke test for the evals harness: stub LLM + recorded gateway, no network."""

import types

from evals.cases import CASES, EvalCase
from evals.grader import deterministic_grade
from evals.report import format_report
from evals.runner import RecordedGateway, run_case


def _text(text):
    return types.SimpleNamespace(type="text", text=text)


def _tool(name, tid, inp):
    return types.SimpleNamespace(type="tool_use", name=name, id=tid, input=inp)


def _resp(content, stop_reason):
    return types.SimpleNamespace(
        content=content,
        stop_reason=stop_reason,
        usage=types.SimpleNamespace(
            input_tokens=10, output_tokens=5,
            cache_read_input_tokens=0, cache_creation_input_tokens=0,
        ),
    )


class _ScriptedLLM:
    """Returns scripted responses; a list of ('tool', name, input) or ('end', text)."""

    def __init__(self, script):
        self._script = list(script)
        self._i = 0

    def create_message(self, **kwargs):
        step = self._script[min(self._i, len(self._script) - 1)]
        self._i += 1
        if step[0] == "tool":
            return _resp([_tool(step[1], f"t{self._i}", step[2])], "tool_use")
        return _resp([_text(step[1])], "end_turn")


def _case(case_id) -> EvalCase:
    return next(c for c in CASES if c.id == case_id)


def test_recorded_gateway_is_offline_and_shaped():
    gw = RecordedGateway()
    assert gw.resolve_genes(["FLT3", "BCL2"])["FLT3"]["curie"].startswith("NCBIGene:")
    nb = gw.neighborhood("NCBIGene:2322")
    assert list(nb.ranked.columns)  # has a ranked DataFrame
    kg = gw.gene_network(["NCBIGene:1", "NCBIGene:2"])
    assert len(kg) == 1


def test_run_case_captures_correct_tool_sequence():
    case = _case("drugs_target_flt3")
    llm = _ScriptedLLM(
        [
            ("tool", "resolve_genes", {"names": ["FLT3"]}),
            ("tool", "gene_neighborhood", {"gene_curie": "NCBIGene:1000"}),
            ("end", "Drugs targeting FLT3 include midostaurin and gilteritinib."),
        ]
    )
    rec = run_case(case, llm, model="claude-sonnet-4-6")
    assert rec["tool_names"] == ["resolve_genes", "gene_neighborhood"]
    grade = deterministic_grade(case, rec["tool_names"])
    assert grade["passed"] is True


def test_deterministic_grade_flags_wrong_tools():
    case = _case("drugs_target_flt3")  # forbids gene_network, expects gene_neighborhood
    grade = deterministic_grade(case, ["resolve_genes", "gene_network"])
    assert grade["passed"] is False
    assert "gene_neighborhood" in grade["missing_tools"]
    assert "gene_network" in grade["used_forbidden"]


def test_deterministic_grade_requires_resolve_before_finder():
    case = _case("drugs_target_flt3")
    grade = deterministic_grade(case, ["gene_neighborhood", "resolve_genes"])
    assert grade["resolved_first"] is False
    assert grade["passed"] is False


def test_offtopic_case_passes_when_no_tools_called():
    case = _case("offtopic_rejected")
    llm = _ScriptedLLM([("end", "I can only help with questions about your gene set.")])
    rec = run_case(case, llm, model="claude-sonnet-4-6")
    assert rec["tool_names"] == []
    grade = deterministic_grade(case, rec["tool_names"])
    assert grade["passed"] is True


def test_deterministic_grade_efficiency_metrics():
    case = _case("drugs_target_flt3")  # max_calls=2, forbids gene_network
    names = ["resolve_genes", "gene_neighborhood", "gene_neighborhood"]
    calls = [
        {"name": "resolve_genes", "args": {"names": ["FLT3"]}},
        {"name": "gene_neighborhood", "args": {"gene_curie": "NCBIGene:2322"}},
        {"name": "gene_neighborhood", "args": {"gene_curie": "NCBIGene:2322"}},  # duplicate
    ]
    g = deterministic_grade(case, names, calls)
    assert g["n_calls"] == 3
    assert g["redundant_calls"] == 1
    assert g["within_budget"] is False
    assert g["passed"] is False  # over budget + redundant

    g2 = deterministic_grade(
        case,
        ["resolve_genes", "gene_neighborhood"],
        [
            {"name": "resolve_genes", "args": {"names": ["FLT3"]}},
            {"name": "gene_neighborhood", "args": {"gene_curie": "NCBIGene:2322"}},
        ],
    )
    assert g2["within_budget"] and g2["redundant_calls"] == 0 and g2["passed"]


def test_format_report_runs():
    case = _case("drugs_target_flt3")
    results = [
        {
            "case": case,
            "tool_names": ["resolve_genes", "gene_neighborhood"],
            "final_text": "ok",
            "iterations": 2,
            "cost_usd": 0.001,
            "deterministic": deterministic_grade(case, ["resolve_genes", "gene_neighborhood"]),
        }
    ]
    report = format_report(results)
    assert "deterministic pass: 1/1" in report
    assert "drugs_target_flt3" in report
