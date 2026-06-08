"""Eval cases for tool-calling. Each case anchors a question to a gene set and states the
expected tool behavior plus an answer rubric for the LLM judge."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional


@dataclass
class EvalCase:
    id: str
    gene_symbols: List[str]
    question: str
    expected_tools: List[str]              # tools that SHOULD be called (set membership)
    expected_order_prefix: List[str] = field(default_factory=list)  # required leading order
    must_resolve_first: bool = True        # resolve_genes must precede any finder
    forbidden_tools: List[str] = field(default_factory=list)
    answer_rubric: str = ""
    disease: Optional[str] = None
    max_calls: Optional[int] = None        # efficiency budget: fail if more tool calls than this


# Note: node_metadata, data_sources, and cell_type_expression call LIVE external APIs
# (NodeAnnotator, the infores catalog, HPA) that are NOT cassetted and need REAL CURIEs, so judge
# those cases on tool SELECTION (the deterministic grade). Their answer faithfulness is only
# meaningful with --gateway live (real resolve -> real CURIEs -> real annotations / raw edge sources).
GENESET = ["FLT3", "NPM1", "NRAS", "BCL2"]

CASES: List[EvalCase] = [
    EvalCase(
        id="drugs_target_flt3",
        gene_symbols=GENESET,
        question="What drugs or chemicals target FLT3?",
        expected_tools=["resolve_genes", "gene_neighborhood"],
        expected_order_prefix=["resolve_genes", "gene_neighborhood"],
        forbidden_tools=["gene_network"],
        max_calls=2,
        answer_rubric=(
            "A good answer names specific drugs/chemicals returned for FLT3 and is grounded in the "
            "neighborhood result, not invented. It should not claim a gene-gene network was used."
        ),
    ),
    EvalCase(
        id="path_bcl2_venetoclax",
        gene_symbols=GENESET,
        question="Is there a path connecting BCL2 and venetoclax?",
        expected_tools=["resolve_genes", "path_between"],
        expected_order_prefix=["resolve_genes", "path_between"],
        answer_rubric=(
            "A good answer resolves BOTH BCL2 and venetoclax to CURIEs, runs path_between, and "
            "describes whether/how they connect via intermediates."
        ),
    ),
    EvalCase(
        id="gene_gene_network",
        gene_symbols=GENESET,
        question="Show how these genes interact with each other.",
        expected_tools=["resolve_genes", "gene_network"],
        expected_order_prefix=["resolve_genes", "gene_network"],
        forbidden_tools=["gene_neighborhood"],
        max_calls=2,
        answer_rubric=(
            "A good answer builds a gene-gene network among the four genes and summarizes the direct "
            "interactions found."
        ),
    ),
    EvalCase(
        id="edge_evidence_flt3",
        gene_symbols=GENESET,
        question="What drugs target FLT3, and what is the evidence for the top one?",
        expected_tools=["resolve_genes", "gene_neighborhood", "edge_evidence"],
        expected_order_prefix=["resolve_genes", "gene_neighborhood"],
        answer_rubric=(
            "A good answer finds FLT3 drugs, then calls edge_evidence for a specific FLT3-drug edge "
            "and reports publications/scores rather than inventing them."
        ),
    ),
    EvalCase(
        id="expression_microglia",
        gene_symbols=GENESET,
        question="Which of these genes are specifically expressed in particular cell types or tissues?",
        expected_tools=["cell_type_expression"],
        answer_rubric=(
            "A good answer uses cell_type_expression (HPA) and reports per-gene cell-type/tissue "
            "specificity. It may first run a finder to have a result to annotate."
        ),
    ),
    EvalCase(
        id="gene_function",
        gene_symbols=GENESET,
        question="What biological processes and molecular functions is FLT3 involved in?",
        expected_tools=["resolve_genes", "node_metadata"],
        expected_order_prefix=["resolve_genes", "node_metadata"],
        forbidden_tools=["gene_network", "path_between"],
        max_calls=2,
        answer_rubric=(
            "A good answer resolves FLT3, calls node_metadata, and reports the GO biological-process / "
            "molecular-function terms and gene type it returns, without inventing functions."
        ),
    ),
    EvalCase(
        id="disease_genes",
        gene_symbols=GENESET,
        question="What genes are associated with acute myeloid leukemia?",
        expected_tools=["resolve_genes", "gene_neighborhood"],
        expected_order_prefix=["resolve_genes", "gene_neighborhood"],
        forbidden_tools=["gene_network"],
        max_calls=2,
        answer_rubric=(
            "A good answer resolves the disease to a CURIE and uses gene_neighborhood (entity-agnostic, "
            "target_categories Gene) to list associated genes, grounded in the result."
        ),
    ),
    EvalCase(
        id="data_provenance",
        gene_symbols=GENESET,
        question=(
            "Find drugs targeting FLT3, then tell me which knowledge sources support those "
            "associations and how reliable they are."
        ),
        expected_tools=["resolve_genes", "gene_neighborhood", "data_sources"],
        expected_order_prefix=["resolve_genes", "gene_neighborhood"],
        answer_rubric=(
            "A good answer runs a finder, then calls data_sources and reports the actual knowledge "
            "sources / knowledge levels it returns, not invented provenance."
        ),
    ),
    EvalCase(
        id="tissue_expression",
        gene_symbols=GENESET,
        question="Which tissues most strongly express these genes?",
        expected_tools=["cell_type_expression"],
        answer_rubric=(
            "A good answer uses cell_type_expression (scope tissue) and reports per-gene tissue "
            "specificity from HPA. It may first run a finder to have a result to annotate."
        ),
    ),
    EvalCase(
        id="offtopic_rejected",
        gene_symbols=GENESET,
        question="What is the capital of France?",
        expected_tools=[],
        forbidden_tools=["resolve_genes", "gene_neighborhood", "path_between", "gene_network"],
        must_resolve_first=False,
        answer_rubric=(
            "A good answer politely declines / redirects to the gene set and does NOT call any "
            "Translator tools for an off-topic question."
        ),
    ),
]
