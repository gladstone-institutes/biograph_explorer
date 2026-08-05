# The Agent Loop

How the Chat Explorer works: the manual Claude agent loop, the tool layer it drives, and the
state and threading rules that keep it safe inside Streamlit.

All code referenced here lives in `src/geneset_translator/agent/`, plus
`src/geneset_translator/ui/chat_page.py` (the composition root).

## Why a manual loop

The loop in `agent_loop.py` is hand-written rather than delegated to the Anthropic SDK's tool
runner, because the app needs four things the runner does not give us:

1. **A hard dollar cap.** Checked before every model call, from actual reported usage, not estimates.
2. **A hard iteration cap.** A runaway plan stops at a known number of steps.
3. **Live status streaming.** Each step, tool call, and result summary is pushed to the UI while
   the turn is still running.
4. **Interruption.** A Stop button that leaves the conversation in an API-valid, resumable state.

The loop depends only on two abstractions, an `LLMClient` protocol and a `ToolRegistry`, so it can
be exercised in tests with neither Anthropic nor TCT present.

## Component map

```
ui/chat_page.py                  composition root: state, threading, rendering
  |
  +-- AgentLoop            (agent_loop.py)     model <-> tools until end_turn / cap / stop
  |     +-- AnthropicLLMClient                 prompt caching, model-param fallback
  |     +-- CostTracker    (cost.py)           actual-usage dollar accounting
  |     +-- ToolRegistry   (tools.py)          schemas + uniform dispatch
  |           +-- 12 Tool objects
  |                 +-- TctGateway   (tct_adapter.py)   the ONLY boundary to TCT
  |                 +-- ResultStash  (tct_adapter.py)   full results kept server-side
  |                 +-- DisplayGraph (display_ops.py)   the curated on-screen graph
  |                 +-- summarize_graph (graph_clustering.py)
  |
  +-- tct_result_to_nx -> viz_normalizer -> network_viz -> streamlit-cytoscape
```

The dependency direction is strict: `agent_loop` knows nothing about TCT, Streamlit, or graphs.
`tools` knows about graphs but reaches TCT only through `TctGateway`. `graph_clustering`,
`display_ops`, and `viz_normalizer` are pure NetworkX, which is what makes them unit-testable
without network access.

## One turn, end to end

1. The user submits a question. `chat_page._handle_prompt` builds the system prompt
   (`build_system_prompt`), constructs an **isolated** `ToolContext` (see
   [Threading](#threading-inside-streamlit)), and submits `AgentLoop.run` to a worker thread.
2. `AgentLoop.run` enters its `while True` loop. Before each model call it checks, in order:
   the stop flag, the cost cap, the iteration cap. Any hit sets `stop_note` and breaks.
3. It calls the model with the full tool schema set and the running message history.
4. `CostTracker.add(response.usage, model)` records the actual spend of that call.
5. The assistant message is appended to `messages` verbatim (content blocks, not text) so the
   conversation stays valid for the next call.
6. `stop_reason` decides what happens next:
   - `end_turn` breaks out; the turn is done.
   - `pause_turn` continues without doing anything, resending to resume.
   - anything else with `tool_use` blocks proceeds to dispatch.
   - anything else with no `tool_use` blocks breaks (nothing left to do).
7. Tools are dispatched (concurrently when eligible, see [Concurrency](#concurrency-and-parallel_safe)),
   each result is JSON-serialized into a `tool_result` block, and the batch is appended as a single
   user message.
8. `ctx.latest_result_id` is set to the **last** result in input order that produced one, so the
   render target is deterministic even though dispatch was concurrent.
9. Loop back to step 2.

On exit, `AgentTurnResult` carries the mutated `messages`, the final text (with `stop_note`
appended in italics if the turn was cut short), and the iteration count.

## The tool layer

### The `Tool` protocol

A tool is any object with four attributes and a `run` method:

```python
name: str
description: str
input_schema: dict     # JSON Schema, sent to the model as-is
parallel_safe: bool    # optional; absent means True
def run(self, args: dict, ctx: ToolContext) -> dict: ...
```

There is no base class and no registration decorator. Adding a tool means writing one object and
adding it to `default_registry()`. The agent loop never changes.

### `ToolContext`

The per-session mutable state a tool may read or write: the gateway, the stash, the resolved gene
CURIEs and their symbol map, the disease, accumulated annotations, the tool log, the query cache,
and the display graph. It is a plain dataclass, injected into every `run` call.

### `ToolRegistry.dispatch`

The uniform entry point. It **never raises**, returning `(result, is_error)` instead, and does four
things beyond calling `run`:

- **Memoization.** Identical `(name, args)` calls to parallel-safe tools are cached in
  `ctx.query_cache` and never re-hit the network. Translator finder calls take roughly a minute
  each, so this matters. Display tools (`parallel_safe=False`) are stateful and always execute.
- **Error surfacing.** An exception becomes `{"error": "TypeError: ..."}` with `is_error=True`,
  which the model sees and can self-correct from. The system prompt tells it to read the error,
  fix its arguments, and retry once.
- **Tool logging.** Successful calls are appended to `ctx.tool_log`, which is what
  `script_generator.py` later replays into a standalone reproduction script.
- **Result tracking.** Any result carrying a `result_id` updates `ctx.latest_result_id`.

### The `result_id` pattern

`ResultStash` (`tct_adapter.py`) keeps full TCT result objects server-side, keyed by a short id
(`res_1`, `res_2`, ...). Finder tools return a compact summary plus that id; follow-up tools take
the id back. A large graph is therefore never serialized into the model's context.

Concretely, `gene_neighborhood` returns the top 25 ranked rows with publications attached to at
most the top 10, not the several thousand edges behind them. `edge_evidence`,
`cell_type_expression`, `data_sources`, `cluster_graph`, and `show_result` all address that full
result by id.

`DictResultStash` is backed by a plain dict (the chat page hands it `st.session_state`'s store) and
locks on `put`, because concurrent tool calls would otherwise race on `len(store)` when computing
the next id.

### Tool catalog

| Tool | What it does | Network I/O | `parallel_safe` |
| --- | --- | --- | --- |
| `resolve_genes` | Names to CURIEs, human-taxon default, `biolink_type` for non-genes | yes | yes |
| `gene_neighborhood` | TCT `Neighborhood_finder` on any CURIE, filtered by target categories | yes | yes |
| `path_between` | TCT `Path_finder` between two entities via intermediate categories | yes | yes |
| `gene_network` | Gene-gene interactions among the set, plus `top_hubs` by degree | yes | yes |
| `edge_evidence` | Publications, supporting text, and scores for one subject/object pair | no (reads stash) | yes |
| `cell_type_expression` | HPA expression plus Node Annotator metadata, scoped | yes | yes |
| `node_metadata` | GO terms, gene type, aliases; standalone, needs no prior result | yes | yes |
| `data_sources` | Primary/aggregator sources, knowledge levels, top contributors | yes (cached catalog) | yes |
| `filter_graph` | Trim the display by neighbor, category, top-N, or top-per-cluster | no | **no** |
| `add_disease_node` | Add the disease and link it to the shown genes | yes (one query) | **no** |
| `show_result` | Set or merge an earlier `result_id` into the display | no | **no** |
| `cluster_graph` | Topology-aware structural digest, and recolor by cluster | no | **no** |

The first eight answer questions. The last four change **what is shown** without running a new
query, which is both faster and cheaper than re-querying to reshape a picture.

## Concurrency and `parallel_safe`

When a model response contains several `tool_use` blocks, the loop runs them concurrently **only if
every tool in that step is parallel-safe**, on a `ThreadPoolExecutor` capped at `max_tool_workers`
(default 5). If any tool in the step mutates shared display state, the whole step runs sequentially
in input order so the edits are deterministic.

`parallel_safe` is read with `getattr(tool, "parallel_safe", True)`, so **the default is safe-to-
parallelize**. Only the four display tools declare `parallel_safe = False`. A new tool that mutates
shared state must declare it explicitly.

One subtlety: the three finder tools do touch the display, calling `ctx.display.reset()` to clear
prior curation when a new base graph arrives. That stays parallel-safe because `reset` is
idempotent and taken under `DisplayGraph`'s lock, so concurrent finders cannot interleave into a
bad state.

The system prompt actively asks the model to batch independent finder calls into one response,
since that is what unlocks the concurrent path.

## Stopping, caps, and cost

### Two stop checkpoints

`should_stop` is polled twice per iteration, and each site leaves a valid conversation:

- **Before the model call.** `messages` already ends with tool results or the user turn, so it is
  valid as-is and the turn is resumable.
- **After the model returns tool calls, before running them.** This is the one that matters for a
  slow batch. Every pending `tool_use` is answered with a synthetic
  `{"stopped": true}` tool result, because the API requires every `tool_use` to be satisfied.
  Skipping this would corrupt the history and break the next turn.

### The spend cap

`CostTracker` (`cost.py`) accumulates **actual** usage only, from each `response.usage`. There is no
estimation anywhere. Prices come from `utils.model_utils.MODEL_PRICING` (overridable via the
`GENESET_MODEL_PRICING` env var) because the API does not return pricing. Cache tokens are priced
at the documented multipliers: reads at 0.1x input, writes at 1.25x input.

The cap is checked before each model call, so the effective ceiling is the cap plus one call. The
UI exposes it as a "Max API spend" slider (default $0.50) and shows running spend in the sidebar.
Node summaries triggered from the graph infopanel are billed to the same tracker so the sidebar
metric matches the log.

### The iteration cap

`max_iterations` defaults to 12. Hitting it produces a `stop_note`, not an exception.

## The LLM client

`AnthropicLLMClient` handles two things beyond the raw SDK call.

**Prompt caching.** Two `ephemeral` breakpoints: one on the system block (the stable
system-prompt-plus-tools prefix) and one auto-applied to the last cacheable message. Together these
mean each loop iteration and each follow-up turn within the cache TTL reads the prefix and the
grown conversation from cache rather than re-billing it as input.

**Model-param fallback.** Adaptive thinking and effort are only sent to models that support them
(`_supports_adaptive_thinking` / `_supports_effort`). If the API still returns a 400, the client
retries progressively leaner: full, then without `output_config`, then base. This is what lets the
sidebar model selector offer older models without the loop breaking.

## Threading inside Streamlit

An agent turn takes tens of seconds to minutes. Running it inline would freeze the app and make the
Stop button unreachable, so `chat_page.py` uses a copy-in / commit-out pattern:

1. **Cached executor.** `_agent_executor()` is an `@st.cache_resource` `ThreadPoolExecutor`, so the
   pool survives reruns.
2. **Copy in.** `_build_turn_context` deep-copies every mutable piece of session state into a fresh
   `ToolContext`: the stash store, gene CURIEs, symbol map, annotations, tool log, query cache, and a
   fresh `DisplayGraph` seeded from the current one. `CostTracker` is deep-copied too, so the
   cumulative spend-cap semantics survive.
3. **The worker never touches `st.*`.** It communicates through a `queue.Queue` for status lines and
   a `threading.Event` for stop.
4. **Poll from a fragment.** `_poll_agent_turn` is `@st.fragment(run_every=1.0)`. While the turn
   runs, only that fragment reruns, so the graph does not remount and the viewport is preserved.
5. **Commit out.** On completion, `_commit_turn` copies the isolated state back on the main thread.

The display commit has a wrinkle worth knowing: edits are committed **onto** the persistent
`chat_display_graph` object rather than replacing it, so its `version` counter stays monotonic
across turns. The graph cache, the viz cache, and the cytoscape component key all key on that
version to detect change; a fresh per-turn object would reset the counter and serve a stale graph.
The commit only fires when the turn actually moved the version past its seeded baseline, so a
text-only turn does not needlessly remount the graph.

## From TCT result to rendered graph

```
TCT result object
  -> _result_to_tct_graph()      handles NeighborhoodResult / PathResult / bare KnowledgeGraph
  -> viz_normalizer.normalize()  TCT attribute shape -> network_viz shape;
                                 derives category from CURIE prefix, computes gene_frequency
  -> display_ops transforms      optional agent curation (filter / merge / add disease / cluster)
  -> display_ops.finalize()      drop orphans, backfill category + is_query_gene, recompute frequency
  -> network_viz + streamlit-cytoscape
```

`viz_normalizer` exists because TCT's `to_networkx` and the app's original `GraphBuilder` produce
different attribute shapes, and the renderer was written against the latter. It is the single place
that knows the translation, and it also supplies the node `category` that TCT's graph has no
concept of.

`finalize` runs after every display edit. Dropping orphans matters: trimming a graph regularly
strands nodes, and floating dots carry no information.

## `cluster_graph` and structural summaries

The agent only ever sees a top slice of a finder result, so it cannot summarize a large graph
reliably from what is in its context. `graph_clustering.summarize_graph` computes structure
deterministically and returns a compact digest instead.

It is topology-aware because modularity-style community detection degenerates on the graph shapes
Translator actually returns:

- **`mesh`** (a real interaction network): Louvain communities on the largest component with
  universal hubs stripped, so peripheral modules become visible instead of being swallowed.
- **`hub_dominated`** (one or a few nodes touch nearly everything, for example a star around a
  single disease): no modularity at all. Facet by node category and edge predicate, report a k-core
  if a denser sub-core exists, and say explicitly that this is what happened.
- **`fragmented`** (many components, none dominant): clusters are the connected components.

Classification uses a hub-degree threshold (`HUB_DEGREE_FRACTION`, half of all other nodes) and a
leaf fraction, the share of non-hub nodes whose every neighbor is a hub. Graphs above
`MAX_EDGES_FOR_CLUSTERING` are clustered on their highest-degree core, and the digest says so in
`notes` rather than silently approximating.

The digest returns counts, topology, method, per-cluster records (size, top members by name,
dominant category, dominant predicates), category and predicate facets, and components. It also
returns a node-to-cluster map used to recolor the display, which is why the prompt tells the model
to describe clusters in terms the user can match against the colored graph.

## The system prompt

`build_system_prompt` is a single function producing the domain framing plus an explicit tool
contract, anchored to the user's actual gene set and disease. The rules that carry the most weight:

- Always `resolve_genes` first; never pass a bare symbol or invent a CURIE.
- Route by question shape: `path_between` for any "how is X connected to Y" question rather than two
  neighborhoods plus inference.
- Batch independent finder calls into one response.
- Call `cluster_graph` before summarizing a large result, and report what it found rather than
  generalizing from the visible top slice.
- Report only what the tools returned, cite the publications and sources carried on the edges, and
  mark any added interpretation as such.

When the disease is known only as a CURIE, the prompt explicitly forbids guessing its common name.
That is a hallucination guard, not a style preference.

## Reproduction scripts

`script_generator.generate_reproduction_script` turns `ctx.tool_log` into a standalone, TCT-only
Python script: it re-runs each finder, rebuilds the NetworkX graph with edge attributes, and exports
each to Cytoscape.js JSON. Identical calls are de-duplicated so each slow query runs once. It
faithfully reproduces `resolve_genes`'s two-pass taxon logic rather than emitting a naive
`batch_lookup`. Steps that need `geneset_translator` (annotation, provenance) or that only edit the
in-app display are emitted as comments, not silently dropped.

## Adding a tool

1. Write an object with `name`, `description`, `input_schema`, and `run(args, ctx)`.
2. Set `parallel_safe = False` if it mutates the display or any other shared state.
3. Add it to `default_registry()`.
4. If it returns a stashable result, `ctx.stash.put(...)` it and return the `result_id` plus a
   **bounded** summary. Cap rows, cap publications, and include a total count.
5. Add a routing rule to `build_system_prompt` if the model would not otherwise know when to reach
   for it.
6. If it is reproducible in pure TCT, add a branch to `script_generator`.

On descriptions: write them as contracts, not docstrings. The existing ones state when to call the
tool, what the arguments must be, and what goes wrong otherwise ("ALWAYS call this FIRST", "works on
ANY resolved CURIE, not just genes"). That text is the main lever on routing behavior.

## What is deliberately not in the loop

- **Model-level retries.** Only unsupported-parameter 400s are retried, and only by dropping
  parameters. A genuine API failure surfaces to the user.
- **Conversation persistence.** History lives in `st.session_state` for the session only.
- **Result eviction.** The stash grows for the session. `_clear_session_results` drops the stash and
  the query cache **together**, since a surviving cache entry could otherwise point at an evicted
  `result_id`.
- **Estimated costs.** Every figure comes from reported usage.

## Tests

| Area | File |
| --- | --- |
| Registry dispatch, caching, error handling | `tests/test_tools_registry.py` |
| Gateway, stash, result conversion | `tests/test_tct_adapter.py` |
| Cost accounting and the cap | `tests/test_cost_cap.py` |
| Display transforms | `tests/test_display_ops.py` |
| Topology classification and clustering | `tests/test_graph_clustering.py` |
| Attribute-shape normalization | `tests/test_viz_normalizer.py` |
| Script generation | `tests/test_script_generator.py` |
| Chat page state and threading helpers | `tests/test_chat_page.py` |

Model-graded tool-calling evals live in `evals/` and run live against the API:

```bash
uv run python -m evals.runner
```
