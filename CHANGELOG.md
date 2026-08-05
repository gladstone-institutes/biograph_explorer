# Changelog

## v2.1.0 (05/08/2026)

- Remove nodes from the Chat Explorer graph using the node's built-in Remove action; removals
  persist across reruns and any orphaned nodes they create are cleared automatically
- Added a caption hint pointing at the Remove action for decluttering the view
- Added `docs/agent_loop.md`, an architecture write-up of the chat agent: the manual loop, the tool
  layer and `result_id` pattern, concurrency and stop/cost caps, the Streamlit threading model, and
  how to add a tool

## v2.0.0 (09/06/2026)

- **New AI Chat Explorer**: ask questions about your gene set in plain English and a Claude agent
  runs the relevant NCATS Translator (TCT) queries and renders the result
- The app is now a two-page navigation shell; the original form-based UI is preserved unchanged as
  the Classic Explorer and still requires no API key
- Agent tools for gene neighborhoods, gene-to-gene interaction networks, pathfinding between two
  entities, edge evidence, cell-type and tissue expression, node metadata, and data provenance
- Topology-aware clustering of large results into modules, or biolink-category facets when the
  network is hub-dominated, with the graph recolored by cluster
- Display-control tools so the agent can trim, focus, merge, and annotate the on-screen graph
  without re-running a query
- Actual-usage cost tracking with a per-session spend cap and a Stop button to interrupt a run
- Exports: the full uncapped graph as Cytoscape.js JSON, ranked results as CSV, and a standalone
  TCT-only reproduction script
- Per-node AI summaries with citations in the graph infopanel
- Model selector backed by the live Anthropic model list, with graceful fallback on older models
- Migrated to the TCT result-class API (`TranslatorResources.load()`, `Neighborhood_finder`,
  `Path_finder`)
- Rotating file logging at `data/logs/chat_agent.log`
- Added a model-graded eval harness for agent tool-calling (`uv run python -m evals.runner`)

## v1.1.0 (18/03/2026)

- Removed gene input limit — any size gene set is now accepted
- Added version number to app UI for easier issue reporting

## v1.0.0 (13/01/2026)

- First release of `geneset_translator`!