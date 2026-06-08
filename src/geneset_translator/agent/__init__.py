"""Agent package: Claude-driven chat exploration over TCT (Translator Component Toolkit).

Layering (SOLID):
- ``tct_adapter``  - the single boundary to TCT (TctGateway), the result stash, and
  TCT-result -> NetworkX conversion. Nothing else in the app imports TCT for the agent path.
- ``viz_normalizer`` - pure-NetworkX mapping of a TCT graph to the attribute shape the
  existing ``ui.network_viz`` renderer expects.
- ``tools``       - the Tool protocol, the agent tools, and the ToolRegistry.
- ``cost``        - actual-usage cost accounting (no estimation).
- ``agent_loop``  - the manual agentic loop (LLMClient protocol + AgentLoop).
- ``script_generator`` - emit a standalone TCT reproduction script.
"""
