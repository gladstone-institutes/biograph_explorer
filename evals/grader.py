"""Grading for tool-calling evals: a deterministic check plus an LLM judge.

deterministic_grade is pluggable and runs without any network. ModelGrader uses an LLM judge
(messages.parse structured output) for answer quality; it is the primary 'model-based grading'
mode requested. Both consume the captured tool-call sequence + final answer.
"""

from __future__ import annotations

import json
from typing import Any, Dict, List, Optional

from pydantic import BaseModel

from .cases import EvalCase

# Finder tools that require a CURIE (so resolve_genes should precede them).
FINDER_TOOLS = {"gene_neighborhood", "path_between", "gene_network"}


class Grade(BaseModel):
    """LLM-judge verdict for one transcript."""

    tool_selection_correct: bool
    resolved_first: bool
    order_ok: bool
    faithful: bool  # every factual claim (drugs, disease name, relationships) traces to a tool result
    answer_score: int  # 1 (poor) .. 5 (excellent)
    notes: str


def _is_subsequence(prefix: List[str], seq: List[str]) -> bool:
    i = 0
    for item in seq:
        if i < len(prefix) and item == prefix[i]:
            i += 1
    return i == len(prefix)


def deterministic_grade(
    case: EvalCase,
    tool_names: List[str],
    tool_calls: Optional[List[Dict[str, Any]]] = None,
) -> Dict[str, Any]:
    """Rule-based grade from the ordered tool-call names; no model required.

    ``tool_calls`` (list of {name, args}) enables redundancy detection (identical repeated calls).
    """
    called = set(tool_names)
    expected = set(case.expected_tools)
    forbidden = set(case.forbidden_tools)

    missing = sorted(expected - called)
    used_forbidden = sorted(forbidden & called)
    tool_selection_correct = not missing and not used_forbidden

    # resolve_genes must precede the first finder call (when required and a finder was used)
    resolved_first = True
    if case.must_resolve_first:
        finder_idx = next((i for i, t in enumerate(tool_names) if t in FINDER_TOOLS), None)
        if finder_idx is not None:
            resolve_idx = next((i for i, t in enumerate(tool_names) if t == "resolve_genes"), None)
            resolved_first = resolve_idx is not None and resolve_idx < finder_idx

    order_ok = _is_subsequence(case.expected_order_prefix, tool_names)

    # Efficiency: total calls, redundant (identical name+args) calls, and budget compliance.
    n_calls = len(tool_names)
    redundant_calls = 0
    if tool_calls:
        seen = set()
        for c in tool_calls:
            key = (c.get("name"), json.dumps(c.get("args", {}), sort_keys=True, default=str))
            if key in seen:
                redundant_calls += 1
            seen.add(key)
    within_budget = case.max_calls is None or n_calls <= case.max_calls

    return {
        "tool_selection_correct": tool_selection_correct,
        "resolved_first": resolved_first,
        "order_ok": order_ok,
        "missing_tools": missing,
        "used_forbidden": used_forbidden,
        "n_calls": n_calls,
        "redundant_calls": redundant_calls,
        "within_budget": within_budget,
        "passed": tool_selection_correct
        and resolved_first
        and order_ok
        and within_budget
        and redundant_calls == 0,
    }


_JUDGE_SYSTEM = (
    "You are grading a biomedical agent's tool use and answer. The agent answers gene-set "
    "questions by calling tools (resolve_genes, gene_neighborhood, path_between, gene_network, "
    "edge_evidence, cell_type_expression). resolve_genes must run before any finder (finders need "
    "CURIEs). Judge strictly and concisely against the rubric."
)


class ModelGrader:
    """LLM-as-judge. Default model is Sonnet 4.6; use Opus for final optimization passes."""

    def __init__(self, model: str = "claude-sonnet-4-6", api_key: Optional[str] = None) -> None:
        import anthropic

        self._client = anthropic.Anthropic(api_key=api_key) if api_key else anthropic.Anthropic()
        self.model = model

    def grade(self, case: EvalCase, tool_names: List[str], final_text: str) -> Grade:
        disease_note = f"\nDisease context: {case.disease}" if case.disease else ""
        user = (
            f"Question: {case.question}\n"
            f"Gene set: {', '.join(case.gene_symbols)}{disease_note}\n"
            f"Expected tools: {case.expected_tools or '(none; off-topic should call no tools)'}\n"
            f"Forbidden tools: {case.forbidden_tools or '(none)'}\n"
            f"Tools the agent actually called, in order: {tool_names or '(none)'}\n\n"
            f"Rubric: {case.answer_rubric}\n\n"
            f"Agent's final answer:\n{final_text}\n\n"
            "Grade strictly. faithful = every factual claim (drug names, the disease name, gene "
            "relationships) is grounded in the tool results and nothing is invented (e.g. guessing "
            "a disease name from a CURIE counts as NOT faithful). "
            "Return a JSON grade: tool_selection_correct, resolved_first, order_ok, faithful "
            "(booleans), answer_score (1-5 integer), and a one-sentence notes."
        )
        resp = self._client.messages.parse(
            model=self.model,
            max_tokens=1024,
            system=_JUDGE_SYSTEM,
            messages=[{"role": "user", "content": user}],
            output_format=Grade,
        )
        return resp.parsed_output
