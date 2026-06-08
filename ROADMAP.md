# GeneSet Translator Roadmap

Planned enhancements that build on the recent dependency refresh (TCT 0.1.6,
streamlit-cytoscape 0.2.x, anthropic 0.105), the dynamic Claude model selector,
the results-parsing refresh, and the click-a-node AI summary feature.

## Priority 1: adopt TCT 0.1.6 Finder classes

Replace the hand-rolled query construction in `trapi_client.py`
(`TCT.select_concept`, `TCT.select_API`, `TCT.format_query_json`, the manual
2-hop query dict, and the parallel `ThreadPoolExecutor` query loop) with
`NeighborhoodFinder` (1-hop), `PathFinder` (gene to disease), and
`NetworkFinder`.

Goal: less custom code, fewer API-selection/predicate bugs, and the performance
improvements baked into 0.1.6. Do this behind the existing `TCT_AVAILABLE` flag
and keep the current code path as a fallback until parity is verified against
`tests/fixtures/alzheimers_test_case.json`. This also overlaps with investigating
API selection/health (only a small fraction of queried APIs currently succeed).

## Priority 2: Claude-native "Ask in plain English" query box

A text box where the user types a question; Claude (reusing the `Anthropic`
client and the model selector) emits structured query parameters (hop pattern,
subject/object biolink categories, predicates, granularity, disease) using a
structured-output schema. These pre-fill the existing sidebar controls for
review-before-run, then flow into the unchanged `TRAPIClient` pipeline.
Gene/disease identifiers are still grounded through TCT `name_resolver`, not the
LLM. Deliberately NOT using TCT's OpenAI-based NL feature, to avoid a second
vendor and API key and to stay Claude-only.

Scope and flow constraints:

- **Queries must be directly about the loaded gene set.** The translation prompt
  is constrained so the question is interpreted only as "how do these genes
  relate to X." Off-topic or open-ended questions (general biology, anything not
  anchored to the gene set) are rejected with a short nudge to rephrase, rather
  than translated. Validate the structured output against the allowed
  category/predicate/hop vocabulary before running so the LLM cannot emit an
  arbitrary or unsupported query.
- **Pipeline:** plain-language query -> translate to TCT/TRAPI query parameters
  -> run -> cache (reuse the existing TRAPI response cache) -> graph the top
  connections (rank by gene frequency / centrality as the network tab already
  does, rather than dumping the full result).
- **Encourage disease grounding.** The input helper text and the prompt steer
  users to ground the question in known disease pathology or to name the disease
  directly in the query, so a disease CURIE is resolvable and the result can be
  filtered/ranked meaningfully (the 2-hop gene -> intermediate -> disease path
  and disease-relevance ranking only work when a disease anchor exists). When no
  disease is supplied, fall back to 1-hop neighborhood discovery and surface a
  hint that adding a disease enables better filtering.

## Backlog

- **Rich citation panel for node summaries:** optional dialog/panel below the
  graph rendering clickable citation cards and the citation-graph subgraph
  (reuse `summary_tab._convert_citations_to_links`,
  `summary_tab.render_citation_card`, and `citation_viewer`) for users who want
  more than the infopanel text.
- **Performance and caching:** unify the JSON cache layout, tune
  `trapi_max_workers`/timeouts, add `node_actions=["expand","remove"]` for
  incremental graph exploration, and profile large-graph rendering.
- **Anthropic 0.105 follow-ups:** optionally surface model
  `capabilities`/context-window from `models.list()` in the selector, and
  consider prompt caching for repeated summary context.
