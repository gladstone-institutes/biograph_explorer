"""Diagnostic: reproduce the EoE example query and break down per-API HTTP outcomes.

Mirrors the app defaults for the "Eosinophilic Esophagitis (10 genes)" example so we
can see exactly which Translator APIs fail, with what error_type, and what the actual
error message / HTTP status is.

Run:  uv run python scripts/diagnose_eoe_http.py
"""
import logging
import sys
from collections import Counter, defaultdict
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    stream=sys.stdout,
)

from geneset_translator.core.trapi_client import TRAPIClient

# --- EoE example, matching app defaults -------------------------------------
GENES = ["CCL26", "CAPN14", "STAT6", "IL13", "TSLP",
         "POSTN", "ALOX15", "CLC", "DSG1", "FLG"]
DISEASE = "MONDO:0005361"
INTERMEDIATE = ["biolink:Protein", "biolink:ChemicalEntity", "biolink:Gene"]
PREDICATE_MIN_DEPTH = 1  # "All Relationships" preset (app default)

# Bound the run: real default is 600s, but HTTP errors (429/5xx) return fast.
# 120s still distinguishes genuine HTTP errors from genuinely slow APIs.
TIMEOUT = int(sys.argv[1]) if len(sys.argv) > 1 else 120


def main() -> None:
    client = TRAPIClient(cache_dir=Path("data/cache"), timeout=TIMEOUT)
    response = client.query_gene_neighborhood(
        GENES,
        disease_curie=DISEASE,
        intermediate_categories=INTERMEDIATE,
        predicate_min_depth=PREDICATE_MIN_DEPTH,
        exclude_literature=True,
        exclude_coexpression=True,
        exclude_homology=True,
        progress_callback=lambda m: print(f"[progress] {m}"),
    )

    timings = response.metadata.get("api_timings", [])

    print("\n" + "=" * 78)
    print(f"EoE query: {response.apis_succeeded}/{response.apis_queried} APIs succeeded, "
          f"{len(response.edges)} edges, timeout={TIMEOUT}s")
    print("=" * 78)

    # Breakdown by error_type
    by_type = Counter(t.get("error_type") or ("ok" if t["success"] else "other")
                      for t in timings)
    print("\nOutcome counts by error_type:")
    for etype, count in by_type.most_common():
        print(f"  {etype:18s} {count}")

    # Group the failures and show the actual messages / statuses
    failures = defaultdict(list)
    for t in timings:
        if not t["success"]:
            failures[t.get("error_type") or "other"].append(
                (t["api_name"], t.get("error"), round(t.get("duration_seconds", 0), 1))
            )

    print("\nFailure detail (api | error | duration):")
    for etype in sorted(failures):
        print(f"\n  --- {etype} ({len(failures[etype])}) ---")
        for name, err, dur in sorted(failures[etype], key=lambda x: x[0]):
            print(f"    {name:45s} {dur:6.1f}s  {err}")

    # Slowest APIs overall (timeout candidates)
    print("\nSlowest 10 APIs:")
    for t in sorted(timings, key=lambda x: x.get("duration_seconds", 0), reverse=True)[:10]:
        flag = "OK " if t["success"] else (t.get("error_type") or "FAIL")
        print(f"    {t['duration_seconds']:6.1f}s  {flag:16s} {t['api_name']}")


if __name__ == "__main__":
    main()
