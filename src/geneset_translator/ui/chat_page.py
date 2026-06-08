"""Chat Explorer page: composition root for the agent UI.

Wires the cached TranslatorResources -> TctGateway -> ToolRegistry -> AgentLoop, runs a user
turn, and renders the latest finder result with the existing network_viz + streamlit-cytoscape
stack. All chat state is namespaced ``chat_*`` in st.session_state so it never collides with the
Classic page.
"""

from __future__ import annotations

import logging
import os
import re
from typing import Any, Optional

import pandas as pd
import streamlit as st

from geneset_translator.agent.agent_loop import (
    AgentLoop,
    AnthropicLLMClient,
    build_system_prompt,
)
from geneset_translator.agent.cost import CostTracker
from geneset_translator.agent.display_ops import DisplayGraph
from geneset_translator.agent.script_generator import generate_reproduction_script
from geneset_translator.agent.tct_adapter import (
    DictResultStash,
    TctGateway,
    load_resources,
    tct_result_to_nx,
)
from geneset_translator.agent.tools import ToolContext, default_registry
from geneset_translator.ui.network_viz import render_network_visualization
from geneset_translator.utils.model_utils import DEFAULT_MODEL_ID, fetch_available_models
from geneset_translator.utils.validators import ValidationError, validate_gene_list

logger = logging.getLogger(__name__)

AGENT_DEFAULT_MODEL = "claude-sonnet-4-6"
LAYOUTS = ["dagre", "cose", "fcose", "cola", "breadthfirst", "circle", "grid", "concentric"]

# Example datasets (CSV path, disease CURIE, disease name) reused from the Classic UI for testing.
EXAMPLES = {
    "Eosinophilic Esophagitis (10 genes)": (
        "data/test_genes/eosinophilic_esophagitis_genes.csv",
        "MONDO:0005361",
        "Eosinophilic Esophagitis",
    ),
    "COVID-19 (10 genes)": (
        "data/test_genes/covid19_genes.csv",
        "MONDO:0100096",
        "COVID-19",
    ),
}


def _looks_like_curie(text: str) -> bool:
    return bool(re.match(r"^[A-Za-z][A-Za-z0-9.]*:\w+$", (text or "").strip()))


def _api_key() -> Optional[str]:
    if os.environ.get("ANTHROPIC_API_KEY"):
        return os.environ["ANTHROPIC_API_KEY"]
    try:
        from geneset_translator.config.settings import get_settings

        return get_settings().claude_api_key
    except Exception:  # noqa: BLE001
        return None


@st.cache_resource(show_spinner="Loading Translator resources (one-time, ~30s)...")
def _cached_resources():
    return load_resources()


def _get_gateway() -> TctGateway:
    return TctGateway(_cached_resources())


def _init_state() -> None:
    ss = st.session_state
    ss.setdefault("chat_messages", [])           # raw Anthropic message history
    ss.setdefault("chat_display", [])            # UI transcript [{role, content}]
    ss.setdefault("chat_stash_store", {})        # result_id -> full TCT result object
    ss.setdefault("chat_tool_log", [])           # successful {name, args}
    ss.setdefault("chat_latest_result_id", None)
    ss.setdefault("chat_query_gene_curies", [])
    ss.setdefault("chat_gene_symbols", [])
    ss.setdefault("chat_curie_to_symbol", {})
    ss.setdefault("chat_annotations", {})        # curie -> annotation_features
    ss.setdefault("chat_node_summaries", {})     # curie -> formatted AI summary text
    ss.setdefault("chat_query_cache", {})        # (tool,args) -> result, saves repeat API calls
    ss.setdefault("chat_disease", "")
    ss.setdefault("chat_disease_label", None)    # human disease name for the prompt (avoids hallucination)
    ss.setdefault("chat_display_graph", DisplayGraph())  # agent-curated display graph (display tools edit this)
    ss.setdefault("chat_cost", CostTracker())


def _build_context(gateway: TctGateway) -> ToolContext:
    """Construct a ToolContext that reads/writes the persisted session-state objects."""
    ss = st.session_state
    return ToolContext(
        tct=gateway,
        stash=DictResultStash(ss.chat_stash_store),
        query_gene_curies=ss.chat_query_gene_curies,
        curie_to_symbol=ss.chat_curie_to_symbol,
        disease_curie=ss.chat_disease or None,
        disease_label=ss.chat_disease_label,
        annotations=ss.chat_annotations,
        tool_log=ss.chat_tool_log,
        latest_result_id=ss.chat_latest_result_id,
        query_cache=ss.chat_query_cache,
        display=ss.chat_display_graph,
    )


# --------------------------------------------------------------------------------------
# Sidebar
# --------------------------------------------------------------------------------------
def _sidebar(api_key: str) -> None:
    ss = st.session_state
    with st.sidebar:
        st.header(":material/genetics: Gene set")
        method = st.radio(
            "Input method",
            ["Example dataset", "Upload CSV", "Paste symbols"],
            key="chat_input_method",
        )

        pending: list = []
        example_disease = ""
        example_disease_name = None

        if method == "Example dataset":
            choice = st.selectbox("Example", list(EXAMPLES), key="chat_example_choice")
            path, example_disease, example_disease_name = EXAMPLES[choice]
            try:
                df = pd.read_csv(path)
                pending = df["gene_symbol"].dropna().astype(str).str.strip().tolist()
                st.caption(f"{len(pending)} genes - {example_disease_name} ({example_disease})")
            except Exception as e:  # noqa: BLE001
                st.error(f"Could not load example: {e}")

        elif method == "Upload CSV":
            up = st.file_uploader(
                "Gene CSV (e.g. a filtered DEG table)",
                type=["csv"],
                key="chat_csv",
                help="Any CSV; pick which column holds the gene symbols below.",
            )
            if up is not None:
                try:
                    df = pd.read_csv(up)
                    cols = list(df.columns)
                    if cols:
                        default_col = "gene_symbol" if "gene_symbol" in cols else cols[0]
                        col = st.selectbox(
                            "Gene name column", cols, index=cols.index(default_col), key="chat_csv_col"
                        )
                        pending = df[col].dropna().astype(str).str.strip().tolist()
                        st.caption(f"{len(pending)} values in '{col}'")
                        with st.expander("Preview CSV"):
                            st.dataframe(df.head(20), width="stretch")
                    else:
                        st.error("That CSV has no columns.")
                except Exception as e:  # noqa: BLE001
                    st.error(f"Could not read CSV: {e}")

        else:  # Paste symbols
            text = st.text_area(
                "Gene symbols",
                placeholder="FLT3, NPM1, NRAS, BCL2",
                help="Comma- or newline-separated HUGO symbols.",
                key="chat_gene_text",
            )
            pending = [g.strip() for g in (text or "").replace(",", "\n").splitlines() if g.strip()]

        disease_val = st.text_input(
            "Disease context (optional)",
            value=ss.chat_disease or example_disease,
            placeholder="e.g. acute myeloid leukemia or MONDO:0018874",
            help="Steers the agent toward disease-anchored reasoning.",
        )
        if st.button("Set gene list", type="primary", use_container_width=True, disabled=not pending):
            _load_gene_set(pending, disease_val, example_disease_name)

        if ss.chat_gene_symbols:
            resolved = len(ss.chat_query_gene_curies)
            st.caption(f"Active: {len(ss.chat_gene_symbols)} genes ({resolved} resolved to CURIEs).")

        st.divider()
        st.header(":material/tune: Settings")

        models = fetch_available_models(api_key)
        ids = [m["id"] for m in models]
        # GENESET_AGENT_MODEL lets a deployment / test pin the default (e.g. a cheap model).
        preferred = os.environ.get("GENESET_AGENT_MODEL", AGENT_DEFAULT_MODEL)
        default = (
            preferred if preferred in ids
            else AGENT_DEFAULT_MODEL if AGENT_DEFAULT_MODEL in ids
            else (ids[0] if ids else DEFAULT_MODEL_ID)
        )
        st.selectbox(
            "Claude model",
            options=ids or [DEFAULT_MODEL_ID],
            index=(ids.index(default) if default in ids else 0),
            key="selected_model",
            format_func=lambda mid: next((m["display_name"] for m in models if m["id"] == mid), mid),
        )
        st.select_slider("Effort", options=["low", "medium", "high"], value="medium", key="chat_effort")
        st.slider("Max API spend ($)", 0.05, 5.0, 0.50, 0.05, key="chat_cost_cap")
        st.slider("Max graph nodes", 50, 400, 200, 10, key="chat_max_nodes")
        st.selectbox("Graph layout", LAYOUTS, index=0, key="chat_layout")

        st.divider()
        cost: CostTracker = ss.chat_cost
        st.metric("Actual spend this session", f"${cost.total_usd:.4f}", help="From response.usage; no estimates.")
        st.caption(
            f"{cost.input_tokens:,} in / {cost.output_tokens:,} out tokens over {len(cost.calls)} calls"
        )
        try:
            from geneset_translator.utils.logging_setup import get_log_file

            st.caption(f":material/description: Session log: `{get_log_file()}`")
        except Exception:  # noqa: BLE001
            pass
        if st.button("Reset conversation", use_container_width=True):
            _reset_conversation()


def _load_gene_set(symbols_raw: list, disease: str, disease_label: Optional[str] = None) -> None:
    ss = st.session_state
    try:
        symbols = validate_gene_list([str(s) for s in symbols_raw])
    except ValidationError as e:
        st.error(f"Invalid gene list: {e}")
        return
    with st.spinner("Resolving genes to CURIEs..."):
        gateway = _get_gateway()
        resolved = gateway.resolve_genes(symbols)
    curies, c2s = [], {}
    unresolved = []
    for sym in symbols:
        info = resolved.get(sym)
        if info and info.get("curie"):
            curies.append(info["curie"])
            c2s[info["curie"]] = sym
        else:
            unresolved.append(sym)
    disease = (disease or "").strip()
    label = disease_label
    if not label and disease and not _looks_like_curie(disease):
        label = disease  # user typed a disease name, not a CURIE
    # New gene set -> fresh conversation/results/caches
    ss.chat_gene_symbols = symbols
    ss.chat_query_gene_curies = curies
    ss.chat_curie_to_symbol = c2s
    ss.chat_disease = disease
    ss.chat_disease_label = label
    _clear_session_results()
    if unresolved:
        st.warning(f"Could not resolve: {', '.join(unresolved)}")
    st.rerun()


def _clear_session_results() -> None:
    """Drop conversation, stashed results, and caches (stash + query_cache must clear together,
    or a cache hit could point at an evicted result_id)."""
    ss = st.session_state
    ss.chat_messages = []
    ss.chat_display = []
    ss.chat_stash_store = {}
    ss.chat_tool_log = []
    ss.chat_latest_result_id = None
    ss.chat_annotations = {}
    ss.chat_query_cache = {}
    ss.chat_display_graph.reset()
    ss.pop("_chat_graph_cache", None)


def _reset_conversation() -> None:
    _clear_session_results()
    st.session_state.chat_cost = CostTracker()
    st.rerun()


# --------------------------------------------------------------------------------------
# Graph panel
# --------------------------------------------------------------------------------------
def _build_current_graph() -> Optional[Any]:
    """Return the graph to render: the agent-curated display graph if non-empty, else the latest
    finder result. Enrichment is merged on top. Cached so layout/AI-summary reruns don't rebuild;
    the cache key tracks the display ``version`` (display tools bump it) and annotation count.
    """
    ss = st.session_state
    display_active = not ss.chat_display_graph.is_empty()
    rid = ss.chat_latest_result_id
    if not display_active and not rid:
        return None

    if display_active:
        cache_key = ("display", ss.chat_display_graph.version, len(ss.chat_annotations))
    else:
        cache_key = ("result", rid, len(ss.chat_annotations))
    cached = ss.get("_chat_graph_cache")
    if cached and cached[0] == cache_key:
        return cached[1]

    if display_active:
        graph = ss.chat_display_graph.snapshot()
    else:
        result = ss.chat_stash_store.get(rid)
        if result is None:
            return None
        graph = tct_result_to_nx(
            result, ss.chat_query_gene_curies, ss.chat_curie_to_symbol, ss.chat_disease or None
        )
    for node, features in ss.chat_annotations.items():  # on-request HPA / Node Annotator enrichment
        if node in graph:
            graph.nodes[node]["annotation_features"] = features
    ss._chat_graph_cache = (cache_key, graph)
    return graph


def _render_graph_panel() -> None:
    from streamlit_cytoscape import InfopanelAction, streamlit_cytoscape

    ss = st.session_state
    if not ss.chat_latest_result_id and ss.chat_display_graph.is_empty():
        st.info("Ask a question below and the resulting graph appears here.", icon=":material/hub:")
        return
    try:
        with st.spinner("Building graph..."):
            graph = _build_current_graph()
            if graph is None:
                return
            viz = render_network_visualization(
                graph,
                query_genes=ss.chat_query_gene_curies,
                sizing_metric="gene_frequency",
                layout=ss.chat_layout,
                max_intermediates=ss.chat_max_nodes,
                node_summaries=ss.chat_node_summaries,
            )
    except Exception as e:  # noqa: BLE001
        logger.exception("graph render failed")
        st.warning(f"Could not render this result as a graph: {e}")
        return

    n_nodes = len(viz["elements"]["nodes"])
    n_edges = len(viz["elements"]["edges"])
    st.caption(f"Showing {n_nodes} nodes / {n_edges} edges (capped at {ss.chat_max_nodes}).")
    key = f"chat_cyto_{ss.chat_latest_result_id}_{ss.chat_layout}_{ss.chat_max_nodes}"
    streamlit_cytoscape(
        viz["elements"],
        layout=viz["layout"],
        node_styles=viz["node_styles"],
        edge_styles=viz["edge_styles"],
        key=key,
        hide_underscore_attrs=True,
        edge_actions=["collapse", "expand"],
        infopanel_actions=[InfopanelAction("ai_summary", "AI Summary", icon="science", spinner=True)],
    )
    _handle_ai_summary(key, graph)


def _handle_ai_summary(component_key: str, graph: Any) -> None:
    """If the infopanel AI Summary action fired for a new node, generate and cache it."""
    ss = st.session_state
    val = ss.get(component_key)
    if not isinstance(val, dict) or val.get("action") != "ai_summary":
        return
    data = val.get("data") or {}
    if data.get("element_group") != "nodes":
        return
    node_id = data.get("element_id")
    if not node_id or node_id in ss.chat_node_summaries:
        return  # already summarized -> avoids a rerun loop
    try:
        with st.spinner(f"Summarizing {node_id}..."):
            ss.chat_node_summaries[node_id] = _node_summary_text(graph, node_id)
    except Exception as e:  # noqa: BLE001
        logger.exception("node summary failed")
        st.warning(f"Could not summarize node: {e}")
        return
    st.rerun()


def _node_summary_text(graph: Any, node_id: str) -> str:
    """Generate a node summary via LLMSummarizer and format it for the infopanel."""
    from geneset_translator.core.llm_summarizer import LLMSummarizer
    from geneset_translator.ui.summary_tab import _format_publication_link

    ss = st.session_state
    summarizer = LLMSummarizer(model=ss.get("selected_model", AGENT_DEFAULT_MODEL))
    summary = summarizer.generate_node_summary(
        graph, node_id, ss.chat_disease or None, ss.chat_query_gene_curies
    )
    text = re.sub(r"\[Citation \d+\]", "", summary.summary_text)
    text = re.sub(r"\s+([.,;:])", r"\1", text)
    text = re.sub(r"\s{2,}", " ", text).strip()

    seen, pmids = set(), []
    for citation in summary.citations:
        for pub in citation.publication_ids:
            if pub not in seen:
                seen.add(pub)
                pmids.append(pub)
    if pmids:
        links = " | ".join(_format_publication_link(p) for p in pmids[:10])
        text = f"{text}\n\nSources: {links}"
    return text


def _render_tools_panel() -> None:
    """Enrichment button + reproduction-script / CSV downloads for the latest result."""
    ss = st.session_state
    rid = ss.chat_latest_result_id
    if not rid and not ss.chat_tool_log:
        return

    col1, col2, col3 = st.columns(3)

    with col1:
        if rid and st.button(
            ":material/biotech: Annotate graph (HPA / GO)",
            use_container_width=True,
            help="Add Human Protein Atlas expression and Node Annotator metadata to the graph.",
        ):
            _annotate_latest()

    with col2:
        if ss.chat_tool_log:
            questions = [m["content"] for m in ss.chat_display if m["role"] == "user"]
            script = generate_reproduction_script(ss.chat_tool_log, questions=questions)
            st.download_button(
                ":material/download: Reproduction script",
                data=script,
                file_name="reproduce_tct_results.py",
                mime="text/x-python",
                use_container_width=True,
            )

    with col3:
        csv = _latest_result_csv()
        if csv is not None:
            st.download_button(
                ":material/table: Ranked results (CSV)",
                data=csv,
                file_name="tct_ranked_results.csv",
                mime="text/csv",
                use_container_width=True,
            )


def _annotate_latest() -> None:
    ss = st.session_state
    gateway = _get_gateway()
    ctx = _build_context(gateway)
    registry = default_registry()
    with st.spinner("Fetching HPA / annotation data..."):
        result, is_error = registry.dispatch(
            "cell_type_expression", {"result_id": ss.chat_latest_result_id}, ctx
        )
    ss.chat_latest_result_id = ctx.latest_result_id
    ss.pop("_chat_graph_cache", None)  # force rebuild with annotations
    if is_error:
        st.warning(f"Annotation failed: {result.get('error')}")
    else:
        st.toast(f"Annotated {result.get('annotated_nodes', 0)} nodes", icon=":material/check:")
    st.rerun()


def _latest_result_csv() -> Optional[str]:
    ss = st.session_state
    result = ss.chat_stash_store.get(ss.chat_latest_result_id) if ss.chat_latest_result_id else None
    if result is None:
        return None
    df = getattr(result, "ranked", None)
    if df is None:
        df = getattr(result, "paths", None)
    if df is None:
        return None
    try:
        return df.to_csv(index=False)
    except Exception:  # noqa: BLE001
        return None


# --------------------------------------------------------------------------------------
# Chat
# --------------------------------------------------------------------------------------
def _render_history() -> None:
    for msg in st.session_state.chat_display:
        with st.chat_message(msg["role"]):
            st.markdown(msg["content"])


def _handle_prompt(prompt: str, api_key: str) -> None:
    ss = st.session_state
    logger.info("USER question: %s", prompt)
    ss.chat_display.append({"role": "user", "content": prompt})
    with st.chat_message("user"):
        st.markdown(prompt)

    with st.chat_message("assistant"):
        status = st.status("Working...", expanded=True)

        def status_cb(message: str) -> None:
            status.write(message)

        gateway = _get_gateway()
        ctx = _build_context(gateway)
        client = AnthropicLLMClient(api_key=api_key)
        registry = default_registry()
        system_prompt = build_system_prompt(
            ss.chat_gene_symbols, ss.chat_disease or None, ss.chat_disease_label
        )
        loop = AgentLoop(
            client,
            registry,
            model=ss.get("selected_model", AGENT_DEFAULT_MODEL),
            system_prompt=system_prompt,
            effort=ss.get("chat_effort", "medium"),
        )

        ss.chat_messages.append({"role": "user", "content": prompt})
        try:
            result = loop.run(
                ss.chat_messages,
                ctx,
                ss.chat_cost,
                cost_cap_usd=ss.get("chat_cost_cap", 0.50),
                status_cb=status_cb,
            )
            ss.chat_latest_result_id = ctx.latest_result_id
            status.update(label=f"Done ({result.iterations} step(s), ${ss.chat_cost.total_usd:.4f})", state="complete")
            answer = result.final_text or "(no answer)"
        except Exception as e:  # noqa: BLE001
            logger.exception("agent turn failed")
            status.update(label="Error", state="error")
            answer = f"Something went wrong: {e}"
        st.markdown(answer)

    logger.info("ASSISTANT answer: %s", answer.replace("\n", " ")[:2000])
    ss.chat_display.append({"role": "assistant", "content": answer})
    st.rerun()


# --------------------------------------------------------------------------------------
# Entry
# --------------------------------------------------------------------------------------
def render() -> None:
    _init_state()
    api_key = _api_key()
    if not api_key:
        st.warning("No Anthropic API key found. See the 'Enable AI Chat' page.", icon=":material/key:")
        return

    _sidebar(api_key)

    st.title(":material/forum: Chat Explorer")
    st.caption(
        "Ask biological questions about your gene set. A Claude agent runs the relevant "
        "Translator (TCT) queries and renders the result below."
    )

    has_genes = bool(st.session_state.chat_query_gene_curies)
    if not has_genes:
        st.info("Enter a gene set in the sidebar and click **Set gene list** to begin.", icon=":material/arrow_back:")

    _render_graph_panel()
    _render_tools_panel()
    st.divider()
    _render_history()

    prompt = st.chat_input(
        "e.g. What drugs target FLT3?" if has_genes else "Set a gene list first...",
        disabled=not has_genes,
    )
    if prompt:
        _handle_prompt(prompt, api_key)
