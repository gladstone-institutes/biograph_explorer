"""Model-graded tool-calling evals for the GeneSet Translator chat agent.

Run live (makes paid LLM calls) with:
    python -m evals.runner

The runner drives the real AgentLoop against a *recorded* TctGateway (canned tool results,
no slow/flaky Translator network) so we measure the model's TOOL SELECTION, then grades each
transcript both deterministically (did it call the expected tools, resolve first, in order)
and with an LLM judge (answer quality). See cases.py for the eval set.
"""
