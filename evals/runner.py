"""Run tool-calling eval cases through the real AgentLoop with a recorded (offline) gateway.

The recorded gateway returns canned, deterministic finder results so we exercise the model's
tool SELECTION without hitting the slow/flaky Translator network. The LLM still makes real calls
(that's the point), so this is paid; run it explicitly:

    python -m evals.runner [--model claude-sonnet-4-6] [--grader-model claude-sonnet-4-6]
                           [--no-judge] [--cases id1,id2]
"""

from __future__ import annotations

import sys
from typing import Any, Dict, List, Optional

from geneset_translator.agent.agent_loop import AgentLoop
from geneset_translator.agent.cost import CostTracker
from geneset_translator.agent.tct_adapter import DictResultStash
from geneset_translator.agent.tools import ToolContext, default_registry

from .cases import CASES, EvalCase
from .gateways import CassetteGateway, RecordedGateway  # noqa: F401 (RecordedGateway re-exported)
from .grader import ModelGrader, deterministic_grade


def run_case(
    case: EvalCase,
    llm_client: Any,
    model: str,
    *,
    gateway: Any = None,
    max_iterations: int = 8,
    cost_cap: float = 1.0,
) -> Dict[str, Any]:
    """Execute one case; capture the ordered tool calls, final answer, cost."""
    from geneset_translator.agent.agent_loop import build_system_prompt

    gateway = gateway if gateway is not None else RecordedGateway()
    registry = default_registry()
    ctx = ToolContext(
        tct=gateway,
        stash=DictResultStash(),
        query_gene_curies=[],
        disease_curie=case.disease,
        disease_label=case.disease,  # offline: label == curie (cases use CURIEs)
    )
    loop = AgentLoop(
        llm_client,
        registry,
        model=model,
        system_prompt=build_system_prompt(case.gene_symbols, case.disease),
        max_iterations=max_iterations,
        effort="medium",
    )
    cost = CostTracker()
    result = loop.run(
        [{"role": "user", "content": case.question}], ctx, cost, cost_cap_usd=cost_cap
    )
    tool_names = [c["name"] for c in ctx.tool_log]
    return {
        "case": case,
        "tool_names": tool_names,
        "tool_calls": list(ctx.tool_log),
        "final_text": result.final_text,
        "iterations": result.iterations,
        "cost_usd": cost.total_usd,
    }


def run_all(
    cases: List[EvalCase],
    llm_client: Any,
    model: str,
    grader: Optional[ModelGrader] = None,
    gateway: Any = None,
) -> List[Dict[str, Any]]:
    results = []
    for case in cases:
        rec = run_case(case, llm_client, model, gateway=gateway)
        rec["deterministic"] = deterministic_grade(case, rec["tool_names"], rec["tool_calls"])
        if grader is not None:
            try:
                rec["judge"] = grader.grade(case, rec["tool_names"], rec["final_text"]).model_dump()
            except Exception as e:  # noqa: BLE001
                rec["judge"] = {"error": str(e)}
        results.append(rec)
    return results


def _parse_args(argv: List[str]) -> Dict[str, Any]:
    opts: Dict[str, Any] = {
        "model": "claude-sonnet-4-6",
        "grader_model": "claude-sonnet-4-6",
        "judge": True,
        "cases": None,
        "gateway": "recorded",  # recorded (offline synthetic) | cassette (recorded-real) | live
    }
    it = iter(argv)
    for a in it:
        if a == "--model":
            opts["model"] = next(it)
        elif a == "--grader-model":
            opts["grader_model"] = next(it)
        elif a == "--no-judge":
            opts["judge"] = False
        elif a == "--cases":
            opts["cases"] = set(next(it).split(","))
        elif a == "--gateway":
            opts["gateway"] = next(it)
    return opts


def _make_gateway(kind: str) -> Any:
    if kind == "cassette":
        return CassetteGateway()
    if kind == "live":
        from geneset_translator.agent.tct_adapter import TctGateway, load_resources

        return TctGateway(load_resources())
    return RecordedGateway()


def main(argv: Optional[List[str]] = None) -> int:
    import os

    from dotenv import load_dotenv

    from geneset_translator.agent.agent_loop import AnthropicLLMClient

    from .report import format_report

    load_dotenv()  # pick up ANTHROPIC_API_KEY from .env, like the app does
    opts = _parse_args(argv if argv is not None else sys.argv[1:])
    if not os.environ.get("ANTHROPIC_API_KEY"):
        print("ANTHROPIC_API_KEY is required to run the live evals.", file=sys.stderr)
        return 2

    cases = CASES if not opts["cases"] else [c for c in CASES if c.id in opts["cases"]]
    client = AnthropicLLMClient()
    grader = ModelGrader(model=opts["grader_model"]) if opts["judge"] else None
    gateway = _make_gateway(opts["gateway"])

    print(f"Running {len(cases)} cases on {opts['model']} [gateway={opts['gateway']}]"
          f"{' (judge: ' + opts['grader_model'] + ')' if grader else ' (no judge)'}...\n")
    results = run_all(cases, client, opts["model"], grader, gateway=gateway)
    print(format_report(results))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
