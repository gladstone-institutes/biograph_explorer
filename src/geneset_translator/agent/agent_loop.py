"""Manual agentic loop for the chat agent.

A manual loop (not the SDK tool runner) so we can: enforce a hard dollar cap and a max
iteration count, stream status to the UI, and feed tool errors back to the model for
self-debugging. Depends only on abstractions: an ``LLMClient`` and a ``ToolRegistry``.
"""

from __future__ import annotations

import json
import logging
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from typing import Any, Callable, List, Optional, Protocol

from .cost import CostTracker
from .tools import ToolContext, ToolRegistry

logger = logging.getLogger(__name__)

StatusCb = Optional[Callable[[str], None]]


class LLMClient(Protocol):
    def create_message(
        self,
        *,
        model: str,
        system: str,
        tools: List[dict],
        messages: List[dict],
        max_tokens: int,
        effort: str = "medium",
    ) -> Any: ...


def _supports_adaptive_thinking(model: str) -> bool:
    m = (model or "").lower()
    return "sonnet-4-6" in m or "opus-4" in m


def _supports_effort(model: str) -> bool:
    m = (model or "").lower()
    return "sonnet-4-6" in m or "opus-4" in m


class AnthropicLLMClient:
    """Concrete LLMClient. Applies adaptive thinking + effort on supporting models and
    progressively drops unsupported params on a 400 (so the model selector can pick older models)."""

    def __init__(self, api_key: Optional[str] = None) -> None:
        import anthropic

        self._anthropic = anthropic
        self._client = anthropic.Anthropic(api_key=api_key) if api_key else anthropic.Anthropic()

    def create_message(
        self,
        *,
        model: str,
        system: str,
        tools: List[dict],
        messages: List[dict],
        max_tokens: int,
        effort: str = "medium",
    ) -> Any:
        # Cache the stable tools+system prefix so every loop iteration (and repeat turns within
        # the 5-min TTL) reads it from cache instead of re-billing it as input tokens.
        system_param = (
            [{"type": "text", "text": system, "cache_control": {"type": "ephemeral"}}]
            if system
            else system
        )
        base = dict(
            model=model,
            max_tokens=max_tokens,
            system=system_param,
            tools=tools,
            messages=messages,
            # Auto-cache the last cacheable message block too, so the growing conversation history
            # is read from cache across loop iterations and follow-up turns (system block is the
            # other breakpoint). Two breakpoints total.
            cache_control={"type": "ephemeral"},
        )
        extra: dict = {}
        if _supports_adaptive_thinking(model):
            extra["thinking"] = {"type": "adaptive"}
        if _supports_effort(model):
            extra["output_config"] = {"effort": effort}

        # Try most-featured first, then drop effort, then drop thinking.
        attempts: List[dict] = [{**base, **extra}]
        if "output_config" in extra:
            attempts.append({**base, **{k: v for k, v in extra.items() if k != "output_config"}})
        attempts.append(base)

        last_err: Optional[Exception] = None
        for kwargs in attempts:
            try:
                return self._client.messages.create(**kwargs)
            except self._anthropic.BadRequestError as e:  # unsupported param -> retry leaner
                last_err = e
                logger.warning("create_message 400 (%s); retrying with fewer params", e)
                continue
        raise last_err  # type: ignore[misc]


@dataclass
class AgentTurnResult:
    messages: List[dict]
    final_text: str
    stop_note: Optional[str]
    iterations: int


def _short_args(args: Any) -> str:
    try:
        return json.dumps(args, default=str)[:160]
    except Exception:  # noqa: BLE001
        return str(args)[:160]


def _summarize_result(result: Any, is_error: bool) -> str:
    """One-line summary of a tool result for logs / status lines."""
    if is_error:
        return f"ERROR: {result.get('error') if isinstance(result, dict) else result}"
    if not isinstance(result, dict):
        return "ok"
    rid = result.get("result_id")
    suffix = f" [{rid}]" if rid else ""
    if "n_results" in result:
        return f"{result['n_results']} results{suffix}"
    if "edge_count" in result:
        return f"{result['edge_count']} edges{suffix}"
    if "n_paths" in result:
        return f"{result['n_paths']} paths{suffix}"
    if "edges_matched" in result:
        return f"{result['edges_matched']} edge(s), {len(result.get('publications', []))} pubs"
    if "expression" in result:
        return f"{result.get('annotated_nodes', 0)} nodes annotated"
    if "node_metadata" in result:
        return f"{result.get('n', 0)} node(s) annotated"
    if "top_sources" in result:  # data_sources
        return f"{result.get('sources_in_result', result.get('total_sources', 0))} source(s)"
    if "nodes_after" in result:  # display tools
        applied = result.get("applied") or ([f"+{result['disease']}"] if result.get("disease") else [])
        extra = f" ({', '.join(applied)})" if applied else ""
        return f"display now {result['nodes_after']} nodes{extra}"
    return f"{len(result)} resolved"


def _extract_text(response: Any) -> str:
    if response is None:
        return ""
    parts = [
        b.text
        for b in getattr(response, "content", [])
        if getattr(b, "type", None) == "text" and getattr(b, "text", None)
    ]
    return "\n".join(parts).strip()


class AgentLoop:
    """Runs one user turn: model <-> tools until end_turn, the cost cap, or the iteration cap."""

    def __init__(
        self,
        client: LLMClient,
        registry: ToolRegistry,
        model: str,
        system_prompt: str,
        *,
        max_tokens: int = 16000,
        max_iterations: int = 12,
        effort: str = "medium",
        max_tool_workers: int = 5,
    ) -> None:
        self.client = client
        self.registry = registry
        self.model = model
        self.system_prompt = system_prompt
        self.max_tokens = max_tokens
        self.max_iterations = max_iterations
        self.effort = effort
        self.max_tool_workers = max_tool_workers

    def run(
        self,
        messages: List[dict],
        ctx: ToolContext,
        cost_tracker: CostTracker,
        cost_cap_usd: float,
        status_cb: StatusCb = None,
        should_stop: Optional[Callable[[], bool]] = None,
    ) -> AgentTurnResult:
        iterations = 0
        stop_note: Optional[str] = None
        last_resp: Any = None

        logger.info(
            "agent turn start: model=%s effort=%s cap=$%.2f genes=%d",
            self.model, self.effort, cost_cap_usd, len(ctx.query_gene_curies),
        )

        while True:
            if should_stop and should_stop():
                # Stop before the next model call; messages already end with tool_results / the user
                # turn, so the conversation stays API-valid and resumable.
                stop_note = "Stopped at your request."
                logger.info(stop_note)
                break
            if cost_tracker.total_usd >= cost_cap_usd:
                stop_note = (
                    f"Stopped: hit the ${cost_cap_usd:.2f} spend cap "
                    f"(actual spend ${cost_tracker.total_usd:.4f})."
                )
                logger.info(stop_note)
                break
            if iterations >= self.max_iterations:
                stop_note = f"Stopped: reached the max of {self.max_iterations} tool iterations."
                logger.info(stop_note)
                break
            iterations += 1

            logger.info("step %d: calling %s ...", iterations, self.model)
            if status_cb:
                status_cb(f"Step {iterations}: thinking ({self.model})...")
            last_resp = self.client.create_message(
                model=self.model,
                system=self.system_prompt,
                tools=self.registry.schemas(),
                messages=messages,
                max_tokens=self.max_tokens,
                effort=self.effort,
            )
            usage = getattr(last_resp, "usage", None)
            call_usd = cost_tracker.add(usage, self.model)
            messages.append({"role": "assistant", "content": last_resp.content})

            stop_reason = getattr(last_resp, "stop_reason", None)
            logger.info(
                "step %d: stop_reason=%s  +$%.4f (total $%.4f)  [in=%s out=%s cache_read=%s]",
                iterations, stop_reason, call_usd, cost_tracker.total_usd,
                getattr(usage, "input_tokens", "?"),
                getattr(usage, "output_tokens", "?"),
                getattr(usage, "cache_read_input_tokens", 0),
            )

            # Surface any text the model wrote alongside its tool calls.
            interim = _extract_text(last_resp)
            if interim and status_cb:
                status_cb(interim)

            if stop_reason == "end_turn":
                break
            if stop_reason == "pause_turn":  # server tool pause; re-send to resume
                continue

            tool_uses = [b for b in last_resp.content if getattr(b, "type", None) == "tool_use"]
            if not tool_uses:
                break

            for tu in tool_uses:
                logger.info("  tool call: %s(%s)", tu.name, _short_args(tu.input))

            if should_stop and should_stop():
                # Stop before running the (possibly slow) batch the model just requested. Answer each
                # pending tool_use with a synthetic "stopped" result so the assistant tool_use is
                # satisfied and the kept conversation remains API-valid.
                stop_note = "Stopped at your request before running the requested tools."
                logger.info(stop_note)
                messages.append(
                    {
                        "role": "user",
                        "content": [
                            {
                                "type": "tool_result",
                                "tool_use_id": tu.id,
                                "content": json.dumps({"stopped": True}),
                                "is_error": False,
                            }
                            for tu in tool_uses
                        ],
                    }
                )
                break

            # Run concurrently only when EVERY tool in the step is parallel_safe (independent
            # network I/O). If any tool mutates shared display state (parallel_safe=False), run the
            # whole step sequentially in input order so the edits are deterministic.
            all_parallel = all(self.registry.is_parallel_safe(tu.name) for tu in tool_uses)
            concurrent = all_parallel and len(tool_uses) > 1
            if status_cb:
                if concurrent:
                    status_cb(f"Running {len(tool_uses)} tools concurrently...")
                elif len(tool_uses) == 1:
                    status_cb(f"Running `{tool_uses[0].name}` {_short_args(tool_uses[0].input)}")
                else:
                    status_cb(f"Running {len(tool_uses)} tools...")

            def _timed_dispatch(tu: Any):
                t0 = time.time()
                result, is_error = self.registry.dispatch(tu.name, tu.input, ctx)
                return result, is_error, time.time() - t0

            if concurrent:
                with ThreadPoolExecutor(max_workers=min(len(tool_uses), self.max_tool_workers)) as ex:
                    outcomes = list(ex.map(_timed_dispatch, tool_uses))
            else:
                outcomes = [_timed_dispatch(tu) for tu in tool_uses]

            tool_results = []
            for tu, (result, is_error, dur) in zip(tool_uses, outcomes):
                summary = _summarize_result(result, is_error)
                logger.info("  tool result: %s -> %s in %.1fs", tu.name, summary, dur)
                if dur > 30:
                    logger.warning("slow tool: %s took %.1fs", tu.name, dur)
                if status_cb:
                    status_cb(f"  {'X' if is_error else 'OK'} {tu.name}: {summary} ({dur:.1f}s)")
                tool_results.append(
                    {
                        "type": "tool_result",
                        "tool_use_id": tu.id,
                        "content": json.dumps(result, default=str),
                        "is_error": is_error,
                    }
                )
            # Deterministic "latest" render target: last result (input order) that produced one,
            # since concurrent dispatch makes the per-call latest_result_id race.
            for tu, (result, is_error, _dur) in reversed(list(zip(tool_uses, outcomes))):
                if not is_error and isinstance(result, dict) and result.get("result_id"):
                    ctx.latest_result_id = result["result_id"]
                    break
            messages.append({"role": "user", "content": tool_results})

        logger.info(
            "agent turn done: %d step(s), $%.4f total, stop=%s",
            iterations, cost_tracker.total_usd, stop_note or "end_turn",
        )
        final_text = _extract_text(last_resp)
        if stop_note:
            final_text = f"{final_text}\n\n_{stop_note}_" if final_text else stop_note
        return AgentTurnResult(
            messages=messages, final_text=final_text, stop_note=stop_note, iterations=iterations
        )


def build_system_prompt(
    gene_symbols: List[str],
    disease: Optional[str] = None,
    disease_label: Optional[str] = None,
) -> str:
    """Domain framing + tool contract, anchored to the user's gene set."""
    genes = ", ".join(gene_symbols) if gene_symbols else "(none provided yet)"
    if disease_label and disease:
        disease_line = (
            f"\nDisease of interest: {disease_label} ({disease}). Prefer disease-anchored "
            "reasoning when relevant.\n"
        )
    elif disease_label:
        disease_line = f"\nDisease of interest: {disease_label}.\n"
    elif disease:
        disease_line = (
            f"\nThe disease of interest is identified ONLY by the CURIE {disease}. Do NOT guess or "
            "state its common name unless a tool result provides it.\n"
        )
    else:
        disease_line = ""
    return (
        "You are a biomedical knowledge-graph assistant. You answer questions about the user's "
        "gene set by calling tools that query the NCATS Translator network via the TCT toolkit, "
        "then explain the findings.\n\n"
        f"The user's gene set: {genes}.{disease_line}\n"
        "Rules:\n"
        "1. ALWAYS call resolve_genes first to turn symbols/names into CURIEs. Every other tool "
        "requires CURIEs; never pass bare symbols and never invent a CURIE. When the entity is a "
        "DISEASE or DRUG named in the question (not a gene), resolve it with biolink_type "
        "('biolink:Disease' or 'biolink:Drug') so it is not mis-matched to a same-named gene.\n"
        "2. Choose the SMALLEST set of finder tools that answers the question: gene_neighborhood "
        "for 'what targets/relates to gene X', gene_network for 'how do these genes interact with each "
        "other'. Use path_between WHENEVER the question is about the relationship between TWO specific "
        "named entities -- 'how is X connected to Y', 'what links gene X to disease/drug Y', 'what "
        "intermediates connect them', or any mechanism question between two entities -- rather than "
        "running separate neighborhoods and inferring the link yourself (e.g. 'what connects FLT3 and "
        "BCL2' -> path_between, NOT two neighborhoods). gene_neighborhood works on "
        "ANY resolved CURIE, not just genes: for 'what genes are associated with disease X' pass the "
        "disease CURIE with target_categories ['biolink:Gene']; for 'what does drug X target' pass "
        "the drug CURIE. To find the 'most connected' genes, read gene_network's top_hubs (ranked by "
        "degree) rather than guessing.\n"
        "3. When you need several independent finder calls (e.g. gene_neighborhood for multiple "
        "genes), emit them together in ONE response so they run concurrently instead of one per "
        "turn.\n"
        "4. Call edge_evidence when the user asks for evidence/proof/publications for a specific "
        "relationship. Call cell_type_expression for expression / cell-type / tissue / immune "
        "specificity (pass scope to narrow it) -- and call it BEFORE making any cell-type or tissue "
        "claim, including when the user says the gene set came from a single-cell / scRNA experiment "
        "or asks which cell types express the genes; do not assert cell-type specificity the tools did "
        "not return. Call node_metadata for what a gene DOES (GO terms / biological process / function "
        "/ gene type / aliases) -- it is standalone and needs no finder result first. Call data_sources "
        "for 'where does this come from / which sources / how reliable', passing an existing result_id.\n"
        "5. To CHANGE WHAT THE GRAPH SHOWS, use the display tools, which edit the on-screen graph "
        "in memory (no new query): filter_graph (trim / 'show only' / 'focus on' / 'genes connected "
        "to X' / 'most connected' via top_n), add_disease_node (add the disease linked to the shown "
        "genes), show_result (bring back or merge an earlier result_id). Prefer these over running a "
        "new finder when the user is asking to re-shape the current picture. For a balanced 'top nodes "
        "of each cluster' view, FIRST cluster_graph, then call filter_graph with top_per_cluster (NOT "
        "top_n: top_n keeps the highest-degree nodes overall and is dominated by the largest cluster on "
        "hub-dominated graphs).\n"
        "5b. When a finder returns a LARGE result (n_results / edge_count / n_paths in the hundreds or "
        "more -- far more than the top_edges you can see), call cluster_graph before summarizing and "
        "base your summary on its structure instead of generalizing from the top slice. In your answer: "
        "(a) STATE the clustering method and topology it reports (e.g. 'louvain communities' vs "
        "'category/predicate facets, because the network is hub-dominated'); (b) describe EVERY cluster, "
        "not just the gene cluster -- for each give its entity type, size, a few representative members "
        "BY NAME (use the names cluster_graph returns; do not reduce a cluster to bare CURIEs or to "
        "'protein-level equivalents'), and how it connects (its dominant_predicates); (c) note the "
        "overall entity-type composition (category_facets: genes vs proteins vs chemicals vs concepts). "
        "cluster_graph is topology-aware: for a hub-dominated / star graph it returns node-category and "
        "edge-predicate facets instead of communities -- report what it returns and do not invent "
        "modules it did not find. It also recolors the graph by cluster: when the clusters ARE biolink "
        "categories (the facet case), the graph uses the standard category palette (each biolink type "
        "has a fixed color shown in the legend), so refer to each cluster by its category so the user "
        "can match your description to the colored graph.\n"
        "6. Report ONLY entities, identifiers, and relationships that appear in the tool results. Do "
        "NOT add drug indications, mechanisms of action, clinical-trial claims, publication IDs, or "
        "other specifics that are not in those results, and never guess the disease name. Use the "
        "exact publications/sources the tools return; if the user asks for detail the tools did not "
        "provide, say it is not in the results rather than inventing it.\n"
        "7. CITE the evidence. For each specific relationship claim, cite the supporting publication "
        "link(s) and/or the primary knowledge source carried on that edge in the tool result "
        "(publications come back as PubMed links); call edge_evidence when a claim needs stronger "
        "support than the finder row carries. End a substantive answer with a 'Sources:' list of the "
        "publication links you actually used.\n"
        "8. DELINEATE provenance. State facts taken from tool results plainly, with their citation. If "
        "you add interpretation or general background that is NOT in the tool results, mark it clearly "
        "(e.g. begin the sentence with 'Interpretation:' or '(general knowledge, not from the "
        "Translator data)'). Never present an uncited specific claim as if it came from the data.\n"
        "9. Use plain text only. Do NOT use emojis or decorative unicode symbols (write 'Yes', not a "
        "check mark).\n"
        "10. Be concise and efficient: there is a hard spend cap, so avoid redundant calls. If a "
        "tool returns an error, read it, fix your arguments (often: resolve the CURIE first), and "
        "retry once.\n"
        "11. The current graph is rendered for the user automatically (your latest finder result, or "
        "whatever the display tools last set); parallel edges between the same two nodes are collapsed "
        "into one link by default, so describe relationship counts rather than dumping every row."
    )
