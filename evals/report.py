"""Format eval results into a readable report table."""

from __future__ import annotations

from typing import Any, Dict, List


def format_report(results: List[Dict[str, Any]]) -> str:
    lines: List[str] = []
    header = f"{'case':<22} {'det':<5} {'calls':<6} {'faith':<6} {'score':<6} notes"
    lines.append(header)
    lines.append("-" * len(header))

    det_pass = 0
    faithful_count = 0
    judged = 0
    scores: List[int] = []
    for rec in results:
        case = rec["case"]
        det = rec["deterministic"]
        passed = det["passed"]
        det_pass += int(passed)
        judge = rec.get("judge") or {}
        score = judge.get("answer_score")
        if isinstance(score, int):
            scores.append(score)
        faith = judge.get("faithful")
        if faith is not None:
            judged += 1
            faithful_count += int(bool(faith))
        note = judge.get("notes") or judge.get("error") or ""
        flag = "PASS" if passed else "FAIL"
        if not passed:
            extra = []
            if det["missing_tools"]:
                extra.append(f"missing={det['missing_tools']}")
            if det["used_forbidden"]:
                extra.append(f"forbidden={det['used_forbidden']}")
            if not det["resolved_first"]:
                extra.append("not-resolved-first")
            if not det["order_ok"]:
                extra.append("bad-order")
            if not det.get("within_budget", True):
                extra.append(f"over-budget({det['n_calls']}>{case.max_calls})")
            if det.get("redundant_calls"):
                extra.append(f"redundant={det['redundant_calls']}")
            note = (note + " | " + "; ".join(extra)).strip(" |")
        calls = f"{det['n_calls']}" + (f"+{det['redundant_calls']}r" if det.get("redundant_calls") else "")
        faith_s = "-" if faith is None else ("yes" if faith else "NO")
        score_s = str(score) if score is not None else "-"
        lines.append(f"{case.id:<22} {flag:<5} {calls:<6} {faith_s:<6} {score_s:<6} {note[:70]}")

    lines.append("-" * len(header))
    avg = sum(scores) / len(scores) if scores else 0.0
    total_cost = sum(rec.get("cost_usd", 0.0) for rec in results)
    total_calls = sum(rec["deterministic"]["n_calls"] for rec in results)
    summary = f"deterministic pass: {det_pass}/{len(results)}"
    if scores:
        summary += f"   avg judge score: {avg:.2f}/5 (n={len(scores)})"
    if judged:
        summary += f"   faithful: {faithful_count}/{judged}"
    summary += f"   tool calls: {total_calls}   agent spend: ${total_cost:.4f}"
    lines.append(summary)
    return "\n".join(lines)
