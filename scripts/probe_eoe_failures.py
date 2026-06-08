"""Probe the failing EoE APIs: print the actual endpoint URL + response body.

The app records only "HTTP {status}" and discards the body, which is where the
real reason lives. This re-sends the same optimized query to the failing APIs and
dumps URL + body so we can classify each failure (unsupported query shape vs stale
URL vs server-down).

Run:  uv run python scripts/probe_eoe_failures.py
"""
import logging
import sys
from copy import deepcopy
from pathlib import Path

import requests

logging.basicConfig(level=logging.WARNING, format="%(message)s", stream=sys.stdout)

from geneset_translator.core.trapi_client import TRAPIClient
import TCT

GENES = ["CCL26", "CAPN14", "STAT6", "IL13", "TSLP",
         "POSTN", "ALOX15", "CLC", "DSG1", "FLG"]
DISEASE = "MONDO:0005361"
INTERMEDIATE = ["biolink:Protein", "biolink:ChemicalEntity", "biolink:Gene"]

FAILING = [
    "RTX KG2 - TRAPI 1.5.0",                # 400
    "Microbiome KP - TRAPI 1.5.0",          # 400
    "CATRAX Pharmacogenomics KP - TRAPI 1.5.0",  # 400
    "MolePro",                              # 404
    "Multiomics KP - TRAPI 1.5.0",          # 404
    "Retriever",                            # 405
    "Connections Hypothesis Provider API",  # 500
]


def build_query(client):
    gene_map = client.normalize_genes(GENES)
    curies = list(gene_map.values())
    hop1 = list(set(TCT.select_concept(["biolink:Gene"], INTERMEDIATE, client.metaKG)))
    hop2 = list(set(TCT.select_concept(INTERMEDIATE, ["biolink:Disease"], client.metaKG)))
    return {
        "message": {"query_graph": {
            "nodes": {
                "n00": {"ids": curies, "categories": ["biolink:Gene"]},
                "n01": {"categories": INTERMEDIATE},
                "n02": {"ids": [DISEASE], "categories": ["biolink:Disease"]},
            },
            "edges": {
                "e00": {"subject": "n00", "object": "n01", "predicates": hop1},
                "e01": {"subject": "n01", "object": "n02", "predicates": hop2},
            },
        }}
    }


def main():
    client = TRAPIClient(cache_dir=Path("data/cache"), timeout=30)
    client._load_translator_resources()

    api_predicates = {}
    for api in set(client.metaKG["API"]):
        api_predicates[api] = list(set(client.metaKG[client.metaKG["API"] == api]["Predicate"]))

    query = build_query(client)

    for api in FAILING:
        url = client.APInames.get(api, "<NOT IN APInames>")
        opt = client._optimize_query_json(deepcopy(query), api, api_predicates)
        e00_preds = opt["message"]["query_graph"]["edges"]["e00"]["predicates"]
        print("\n" + "=" * 78)
        print(f"API : {api}")
        print(f"URL : {url}")
        print(f"e00 predicates sent: {len(e00_preds)}  |  edges in query graph: "
              f"{len(opt['message']['query_graph']['edges'])}")
        try:
            r = requests.post(url, json=opt, timeout=30)
            print(f"HTTP {r.status_code}")
            body = r.text.strip().replace("\n", " ")
            print(f"BODY: {body[:600]}")
        except Exception as e:
            print(f"EXC : {type(e).__name__}: {e}")


if __name__ == "__main__":
    main()
