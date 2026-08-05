"""Chat Explorer page: composition root for the agent UI.

Wires the cached TranslatorResources -> TctGateway -> ToolRegistry -> AgentLoop, runs a user
turn, and renders the latest finder result with the existing network_viz + streamlit-cytoscape
stack. All chat state is namespaced ``chat_*`` in st.session_state so it never collides with the
Classic page.
"""

from __future__ import annotations

import copy
import json
import logging
import os
import queue
import re
import threading
from concurrent.futures import ThreadPoolExecutor
from typing import Any, Optional

import networkx as nx
import pandas as pd
import streamlit as st

from geneset_translator.agent.agent_loop import (
    AgentLoop,
    AnthropicLLMClient,
    build_system_prompt,
)
from geneset_translator.agent.cost import CostTracker
from geneset_translator.agent.display_ops import DisplayGraph, drop_orphans, finalize
from geneset_translator.agent.script_generator import generate_reproduction_script
from geneset_translator.agent.tct_adapter import (
    DictResultStash,
    TctGateway,
    load_resources,
    tct_result_to_nx,
)
from geneset_translator.agent.tools import ToolContext, default_registry
from geneset_translator.ui.network_viz import (
    calculate_edge_font_size,
    cluster_color_map,
    render_network_visualization,
)
from geneset_translator.utils.formatters import strip_emoji
from geneset_translator.utils.model_utils import DEFAULT_MODEL_ID, fetch_available_models
from geneset_translator.utils.validators import ValidationError, validate_gene_list

logger = logging.getLogger(__name__)

AGENT_DEFAULT_MODEL = "claude-sonnet-4-6"
LAYOUTS = ["dagre", "cose", "fcose", "cola", "breadthfirst", "circle", "grid", "concentric"]
# Metrics meaningful on a chat graph (degree is intrinsic; gene_frequency is computed by viz_normalizer).
SIZING_METRICS = ["degree", "gene_frequency"]

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


def _starter_questions(
    symbols: list, disease_label: Optional[str] = None, disease_curie: str = ""
) -> list:
    """Ordered, context-aware starter questions for a loaded gene set.

    Pure (no Streamlit) so it is unit-testable. Each question is worded to match the agent's routing
    rules in agent_loop.build_system_prompt and collectively demonstrates the app's distinctive
    features for a DEG workflow: disease->gene association + druggability (gene_neighborhood),
    gene-gene interactions (gene_network), shared intermediates (path_between), cell-type/tissue
    expression (cell_type_expression), and provenance (data_sources). The first item is the headline
    translational question (also reused as the chat-box placeholder).
    """
    dz = (disease_label or disease_curie or "").strip()
    questions: list = []
    if dz:
        questions.append(
            f"Which of these genes are most associated with {dz}, "
            "and what drugs or chemicals target them?"
        )
    else:
        questions.append("What drugs or chemicals target these genes, and what is the evidence?")
    questions.append("How do these genes interact with each other?")
    if len(symbols) >= 2:
        questions.append(
            f"What intermediate genes or proteins connect {symbols[0]} and {symbols[1]}?"
        )
    questions.append("Which cell types or tissues most express these genes?")
    questions.append("What knowledge sources support these associations, and how reliable are they?")
    return questions


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


@st.cache_resource
def _agent_executor() -> ThreadPoolExecutor:
    """Shared pool that runs an agent turn off the main thread, so the app stays interactive (and the
    Stop button responsive) while the turn runs. Mirrors classic_explorer's _node_summary_executor."""
    return ThreadPoolExecutor(max_workers=2)


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
    # In-flight turn state (the agent runs in a worker thread so the Stop button stays responsive).
    ss.setdefault("chat_running", False)
    ss.setdefault("chat_turn_future", None)
    ss.setdefault("chat_turn_stop", None)        # threading.Event
    ss.setdefault("chat_turn_status_q", None)    # queue.Queue[str] (worker -> main status)
    ss.setdefault("chat_turn_status_lines", [])  # main-thread-owned accumulation of drained status
    ss.setdefault("chat_turn_ctx", None)         # isolated ToolContext, committed on completion
    ss.setdefault("chat_turn_store", None)        # isolated stash dict, committed on completion
    ss.setdefault("chat_turn_cost", None)        # isolated CostTracker copy, committed on completion
    ss.setdefault("chat_turn_display_base", None)  # display.version baseline (detect a real display edit)


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
        st.toggle(
            "Collapse parallel edges",
            value=True,
            key="chat_collapse_edges",
            help="Merge multiple relationships between the same two nodes into one link. "
            "Click a collapsed link's 'expand' action to see the individual edges.",
        )
        st.slider("Node size", 10, 100, 30, 5, key="chat_node_size", help="Base node size in pixels.")
        st.slider("Edge width", 1, 10, 2, 1, key="chat_edge_width", help="Edge line width in pixels.")
        st.toggle(
            "Size nodes by metric",
            value=False,
            key="chat_use_metric_sizing",
            help="Scale node size by the selected metric instead of uniform sizing.",
        )
        st.selectbox(
            "Size metric", SIZING_METRICS, index=0, key="chat_sizing_metric",
            help="Metric used when 'Size nodes by metric' is on.",
        )

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
    ss.pop("_chat_viz_cache", None)
    ss.pop("_chat_full_cyjs", None)
    if ss.get("chat_turn_stop") is not None:  # signal any in-flight worker to wind down
        ss.chat_turn_stop.set()
    _clear_turn_state()


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
    graph = drop_orphans(graph)  # never display floating dots (covers finder results too)
    for node, features in ss.chat_annotations.items():  # on-request HPA / Node Annotator enrichment
        if node in graph:
            graph.nodes[node]["annotation_features"] = features
    ss._chat_graph_cache = (cache_key, graph)
    return graph


def _graph_to_cyjs(graph: Any) -> tuple:
    """Serialize a NetworkX graph to Cytoscape.js JSON. Pure (no session state) so it is testable.
    Returns ``(json_str, n_nodes, n_edges)``. ``default=str`` keeps non-JSON edge attrs serializable."""
    return (
        json.dumps(nx.cytoscape_data(graph), default=str),
        graph.number_of_nodes(),
        graph.number_of_edges(),
    )


def _full_graph_cyjs() -> Optional[tuple]:
    """Cytoscape.js JSON for the FULL (uncapped) current graph, or None if there is no graph.

    Uses ``_build_current_graph()`` (the complete, pre-sampling graph the renderer starts from, with
    annotations merged) so this needs no finder re-run. Memoized by graph identity + size so a large
    graph is not re-encoded on every Streamlit rerun."""
    graph = _build_current_graph()
    if graph is None or graph.number_of_nodes() == 0:
        return None
    ss = st.session_state
    key = (id(graph), graph.number_of_nodes(), graph.number_of_edges())
    cached = ss.get("_chat_full_cyjs")
    if cached and cached[0] == key:
        return cached[1]
    payload = _graph_to_cyjs(graph)
    ss._chat_full_cyjs = (key, payload)
    return payload


def _graph_component_key(
    rid: Optional[str], layout: str, max_nodes: int, display_version: int, collapse: bool,
    use_metric_sizing: bool = False, sizing_metric: str = "degree",
) -> str:
    """Stable cytoscape component key. Changing it REMOUNTS the component (guaranteeing a refresh);
    keeping it lets the component update in place. Includes display_version so agent display edits
    refresh, collapse so toggling re-inits, and the metric-sizing mode (which changes node element
    data). Node size / edge width are deliberately EXCLUDED -- they are styling-only and update in
    place so the layout/positions are preserved. node_summaries are excluded too (an AI Summary must
    not remount the graph)."""
    return (
        f"chat_cyto_{rid}_{layout}_{max_nodes}_v{display_version}_c{int(collapse)}"
        f"_m{int(use_metric_sizing)}_{sizing_metric}"
    )


def _meta_edge_style(edge_width: int) -> dict:
    """Cytoscape style for collapsed parallel ("meta") edges. The component styles meta-edges
    SEPARATELY from the per-predicate EdgeStyles, so the edge-width control only reaches them if we
    set width here too. Mirrors the Classic UI's meta-edge style (preserved colors + autorotate label)
    and adds width so the slider applies to collapsed links as well."""
    return {
        "text-rotation": "autorotate",
        "width": edge_width,
        "font-size": calculate_edge_font_size(edge_width),
        "line-color": "data(_preservedLineColor)",
        "target-arrow-color": "data(_preservedArrowColor)",
    }


def _render_cluster_legend(node_elements: list) -> None:
    """When the graph is colored by cluster (cluster_graph tool), show a compact id -> color legend.
    Uses the same cluster_color_map as network_viz so swatch colors match the rendered nodes."""
    present = []
    for n in node_elements:
        cid = n.get("data", {}).get("cluster")
        if cid and cid not in present:
            present.append(cid)
    if not present:
        return
    cmap = cluster_color_map(present)
    swatches = " ".join(
        f"<span style='display:inline-block;width:10px;height:10px;background:{cmap[c]};"
        f"border-radius:2px;margin:0 3px 0 8px;vertical-align:middle'></span>{c}"
        for c in sorted(present)
    )
    st.markdown(f"<small><b>Clusters:</b>{swatches}</small>", unsafe_allow_html=True)


def _render_graph_panel() -> None:
    from streamlit_cytoscape import InfopanelAction, streamlit_cytoscape

    ss = st.session_state
    if not ss.chat_latest_result_id and ss.chat_display_graph.is_empty():
        st.info("Ask a question below and the resulting graph appears here.", icon=":material/hub:")
        return

    collapse = ss.get("chat_collapse_edges", True)
    node_size = ss.get("chat_node_size", 30)
    edge_width = ss.get("chat_edge_width", 2)
    use_metric = ss.get("chat_use_metric_sizing", False)
    sizing_metric = ss.get("chat_sizing_metric", SIZING_METRICS[0])
    display_version = ss.chat_display_graph.version
    # Memoize the rendered payload so layout/AI-summary reruns don't re-sample the graph. Keyed by
    # everything that changes the payload (graph identity, layout, cap, collapse, #summaries, sizing).
    viz_key = (
        ss.chat_latest_result_id, display_version, ss.chat_layout, ss.chat_max_nodes,
        collapse, len(ss.chat_node_summaries), len(ss.chat_annotations),
        node_size, edge_width, use_metric, sizing_metric,
    )
    cached_viz = ss.get("_chat_viz_cache")
    try:
        graph = _build_current_graph()
        if graph is None:
            return
        if cached_viz and cached_viz[0] == viz_key:
            viz = cached_viz[1]
        else:
            with st.spinner("Building graph..."):
                viz = render_network_visualization(
                    graph,
                    query_genes=ss.chat_query_gene_curies,
                    layout=ss.chat_layout,
                    max_intermediates=ss.chat_max_nodes,
                    node_summaries=ss.chat_node_summaries,
                    base_node_size=node_size,
                    edge_width=edge_width,
                    use_metric_sizing=use_metric,
                    sizing_metric=sizing_metric,
                )
            ss._chat_viz_cache = (viz_key, viz)
    except Exception as e:  # noqa: BLE001
        logger.exception("graph render failed")
        st.warning(f"Could not render this result as a graph: {e}")
        return

    n_nodes = len(viz["elements"]["nodes"])
    n_edges = len(viz["elements"]["edges"])
    st.caption(
        f"Showing {n_nodes} nodes / {n_edges} edges (capped at {ss.chat_max_nodes}). "
        "Select a node and use its Remove action to declutter the view."
    )
    _render_cluster_legend(viz["elements"]["nodes"])
    key = _graph_component_key(
        ss.chat_latest_result_id, ss.chat_layout, ss.chat_max_nodes, display_version, collapse,
        use_metric, sizing_metric,
    )
    streamlit_cytoscape(
        viz["elements"],
        layout=viz["layout"],
        node_styles=viz["node_styles"],
        edge_styles=viz["edge_styles"],
        key=key,
        hide_underscore_attrs=True,
        node_actions=["remove"],  # built-in per-node Remove so users can delete cluttering nodes
        edge_actions=["collapse", "expand"],
        collapse_parallel_edges=collapse,
        meta_edge_style=_meta_edge_style(edge_width),  # so the edge-width control reaches collapsed edges
        infopanel_actions=[InfopanelAction("ai_summary", "AI Summary", icon="science", spinner=True)],
    )
    _handle_node_removal(key, graph)
    _handle_ai_summary(key, graph)


def _handle_node_removal(component_key: str, graph: Any) -> None:
    """Persist the component's built-in 'remove' node action. The fork removes the node(s) client-side
    and returns their ids; we drop them from the curated display graph so the removal sticks across
    reruns (and re-finalize, which also clears any orphans the removal creates). Removing edits the
    display, which bumps its version -> a new component key, so the stale event is not reprocessed."""
    ss = st.session_state
    val = ss.get(component_key)
    if not isinstance(val, dict) or val.get("action") != "remove":
        return
    node_ids = (val.get("data") or {}).get("node_ids") or []
    present = [n for n in node_ids if graph.has_node(n)]
    if not present:
        return  # nothing to remove (stale/empty event) -> avoids a rerun loop
    trimmed = nx.MultiDiGraph(graph)
    trimmed.remove_nodes_from(present)
    trimmed = finalize(trimmed, ss.chat_query_gene_curies)
    ss.chat_display_graph.replace(trimmed)
    ss.pop("_chat_graph_cache", None)
    ss.pop("_chat_viz_cache", None)
    ss.pop("_chat_full_cyjs", None)
    logger.info("user removed %d node(s) from the view: %s", len(present), present)
    st.rerun()


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
    # Bill the node-summary Claude call to the same session tracker the panel shows, so the panel
    # matches the log (this call is otherwise logged but never counted in the sidebar metric).
    usage = getattr(summarizer, "last_usage", None)
    if usage is not None:
        added = ss.chat_cost.add(usage, summarizer.model)
        logger.info("session cost +$%.4f (node summary %s); total $%.4f",
                    added, node_id, ss.chat_cost.total_usd)
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
    return strip_emoji(text)


def _render_tools_panel() -> None:
    """Enrichment button + reproduction-script / CSV downloads for the latest result."""
    ss = st.session_state
    rid = ss.chat_latest_result_id
    if not rid and not ss.chat_tool_log:
        return

    col1, col2, col3, col4 = st.columns(4)

    with col1:
        if rid and st.button(
            ":material/biotech: Annotate graph (HPA / GO)",
            use_container_width=True,
            help="Add Human Protein Atlas expression and Node Annotator metadata to the graph.",
        ):
            _annotate_latest()

    with col2:
        full = _full_graph_cyjs()
        if full is not None:
            cyjs, n_nodes, n_edges = full
            st.download_button(
                ":material/download: Full graph (.cyjs)",
                data=cyjs,
                file_name="tct_full_graph.cyjs",
                mime="application/json",
                use_container_width=True,
                help=f"The complete {n_nodes}-node / {n_edges}-edge graph (the on-screen view is "
                "sampled to the 'Max graph nodes' setting). Import into Cytoscape Desktop.",
            )

    with col3:
        if ss.chat_tool_log:
            questions = [m["content"] for m in ss.chat_display if m["role"] == "user"]
            script = generate_reproduction_script(ss.chat_tool_log, questions=questions)
            st.download_button(
                ":material/download: Reproduction script",
                data=script,
                file_name="reproduce_tct_results.py",
                mime="text/x-python",
                use_container_width=True,
                help="Standalone TCT-only Python that re-runs each finder, rebuilds the NetworkX "
                "graphs (with edge predicates, publications, and sources), and exports them to "
                "Cytoscape.js JSON (.cyjs) for import into Cytoscape Desktop.",
            )

    with col4:
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


def _build_turn_context(gateway: TctGateway, turn_store: dict) -> ToolContext:
    """An ISOLATED ToolContext for one turn: every mutable piece is copied from session state so the
    worker thread never mutates shared objects. Committed back on the main thread when the turn ends."""
    ss = st.session_state
    display = DisplayGraph()
    if not ss.chat_display_graph.is_empty():
        display.replace(ss.chat_display_graph.snapshot())
    return ToolContext(
        tct=gateway,
        stash=DictResultStash(turn_store),
        query_gene_curies=list(ss.chat_query_gene_curies),
        curie_to_symbol=dict(ss.chat_curie_to_symbol),
        disease_curie=ss.chat_disease or None,
        disease_label=ss.chat_disease_label,
        annotations=dict(ss.chat_annotations),
        tool_log=list(ss.chat_tool_log),
        latest_result_id=ss.chat_latest_result_id,
        query_cache=dict(ss.chat_query_cache),
        display=display,
    )


def _commit_turn(
    ss: Any, turn_ctx: ToolContext, turn_store: dict, turn_cost: CostTracker, result: Any,
    display_base_version: int,
) -> str:
    """Copy a finished turn's isolated state back into session state (main thread only) and return the
    display answer. Pure w.r.t. Streamlit (``ss`` is any attribute-bearing object) so it is testable.

    The turn's display edits are committed ONTO the persistent ``ss.chat_display_graph`` (rather than
    replacing the object) so its ``version`` stays monotonic across turns -- the graph/viz caches and
    the cytoscape component key key on that version to detect change, and a per-turn fresh object would
    reset the counter and stale-cache the graph. Only commit when the turn actually changed the display
    (version moved past the seeded baseline), so a text-only turn doesn't needlessly remount the graph.
    """
    ss.chat_messages = result.messages
    ss.chat_stash_store = turn_store
    ss.chat_tool_log = turn_ctx.tool_log
    ss.chat_query_cache = turn_ctx.query_cache
    ss.chat_annotations = turn_ctx.annotations
    ss.chat_latest_result_id = turn_ctx.latest_result_id
    if turn_ctx.display.version != display_base_version:
        if turn_ctx.display.is_empty():
            ss.chat_display_graph.reset()
        else:
            ss.chat_display_graph.replace(turn_ctx.display.snapshot())
    ss.chat_cost = turn_cost
    return strip_emoji(result.final_text or "(no answer)")


def _clear_turn_state() -> None:
    ss = st.session_state
    ss.chat_running = False
    for k in ("chat_turn_future", "chat_turn_stop", "chat_turn_status_q",
              "chat_turn_ctx", "chat_turn_store", "chat_turn_cost", "chat_turn_display_base"):
        ss[k] = None
    ss.chat_turn_status_lines = []


def _handle_prompt(prompt: str, api_key: str) -> None:
    """Start an agent turn in a worker thread and hand control back to Streamlit so the Stop button
    and status stay live. The worker runs on an ISOLATED context/cost (copy-in); results are committed
    back by ``_poll_agent_turn`` when the future completes. The worker never touches st.session_state."""
    ss = st.session_state
    logger.info("USER question: %s", prompt)
    ss.chat_display.append({"role": "user", "content": prompt})

    gateway = _get_gateway()  # cached resource; must be built on the main thread
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

    turn_store = dict(ss.chat_stash_store)
    turn_ctx = _build_turn_context(gateway, turn_store)
    display_base_version = turn_ctx.display.version  # baseline; the turn changed the display iff this moves
    turn_cost = copy.deepcopy(ss.chat_cost)  # cumulative copy keeps the existing spend-cap semantics
    messages_copy = list(ss.chat_messages) + [{"role": "user", "content": prompt}]

    stop_event = threading.Event()
    status_q: "queue.Queue[str]" = queue.Queue()
    future = _agent_executor().submit(
        loop.run,
        messages_copy,
        turn_ctx,
        turn_cost,
        ss.get("chat_cost_cap", 0.50),
        status_cb=status_q.put,
        should_stop=stop_event.is_set,
    )

    ss.chat_running = True
    ss.chat_turn_future = future
    ss.chat_turn_stop = stop_event
    ss.chat_turn_status_q = status_q
    ss.chat_turn_status_lines = []
    ss.chat_turn_ctx = turn_ctx
    ss.chat_turn_store = turn_store
    ss.chat_turn_cost = turn_cost
    ss.chat_turn_display_base = display_base_version
    st.rerun()


@st.fragment(run_every=1.0)
def _poll_agent_turn() -> None:
    """Poll the in-flight turn: stream status, expose Stop, and commit results when done. While the
    turn runs only this fragment reruns (the graph stays put); completion triggers a full app rerun."""
    ss = st.session_state
    future = ss.get("chat_turn_future")
    if future is None:
        return

    q = ss.get("chat_turn_status_q")
    if q is not None:
        while True:
            try:
                ss.chat_turn_status_lines.append(q.get_nowait())
            except queue.Empty:
                break

    status = st.status("Working...", expanded=True)
    for line in ss.chat_turn_status_lines[-12:]:
        status.write(line)

    if not future.done():
        if ss.chat_turn_stop is not None and ss.chat_turn_stop.is_set():
            status.update(label="Stopping after the current step...", state="running")
            st.button("Stopping...", disabled=True, use_container_width=True, key="chat_stop_pending")
        elif st.button(":material/stop_circle: Stop", use_container_width=True, key="chat_stop"):
            ss.chat_turn_stop.set()
            logger.info("user requested stop")
            st.rerun(scope="fragment")
        return

    try:
        result = future.result()
        answer = _commit_turn(
            ss, ss.chat_turn_ctx, ss.chat_turn_store, ss.chat_turn_cost, result,
            ss.chat_turn_display_base,
        )
        status.update(
            label=f"Done ({result.iterations} step(s), ${ss.chat_cost.total_usd:.4f})", state="complete"
        )
    except Exception as e:  # noqa: BLE001
        logger.exception("agent turn failed")
        answer = f"Something went wrong: {e}"
        status.update(label="Error", state="error")
    logger.info("ASSISTANT answer: %s", answer.replace("\n", " ")[:2000])
    ss.chat_display.append({"role": "assistant", "content": answer})
    _clear_turn_state()
    st.rerun(scope="app")


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

    ss = st.session_state

    # While a turn runs it lives in a worker thread; show the live status + Stop button (this fragment
    # polls every second) instead of the chat input, so the user can interrupt and there is no
    # double-submit.
    if ss.chat_running:
        _poll_agent_turn()
        return

    starters = (
        _starter_questions(ss.chat_gene_symbols, ss.chat_disease_label, ss.chat_disease)
        if has_genes
        else []
    )

    # Clickable starter questions: shown only before the conversation begins, to demonstrate the
    # app's distinctive features for a DEG + disease workflow without cluttering an active chat.
    clicked_q = None
    if starters and not ss.chat_display:
        st.caption("Try one of these to get started:")
        cols = st.columns(2)
        for i, q in enumerate(starters):
            if cols[i % 2].button(q, key=f"chat_starter_{i}", use_container_width=True):
                clicked_q = q

    placeholder = ("e.g. " + starters[0]) if starters else "Set a gene list first..."
    prompt = st.chat_input(placeholder, disabled=not has_genes)
    submitted = prompt or clicked_q
    if submitted:
        _handle_prompt(submitted, api_key)
