"""Tests for agent.cost.CostTracker and the AgentLoop hard cost/iteration caps."""

import types

import pytest

from geneset_translator.agent.agent_loop import AgentLoop
from geneset_translator.agent.cost import CostTracker
from geneset_translator.agent.tct_adapter import DictResultStash
from geneset_translator.agent.tools import ToolContext, ToolRegistry


def _usage(input_tokens=0, output_tokens=0, cache_read=0, cache_create=0):
    return types.SimpleNamespace(
        input_tokens=input_tokens,
        output_tokens=output_tokens,
        cache_read_input_tokens=cache_read,
        cache_creation_input_tokens=cache_create,
    )


def _text_block(text):
    return types.SimpleNamespace(type="text", text=text)


def _tool_block(name, tool_id, inp):
    return types.SimpleNamespace(type="tool_use", name=name, id=tool_id, input=inp)


def _resp(content, stop_reason, usage=None):
    return types.SimpleNamespace(content=content, stop_reason=stop_reason, usage=usage)


class _StubLLM:
    """Returns scripted responses in order; repeats the last one forever."""

    def __init__(self, responses):
        self._responses = responses
        self.calls = 0

    def create_message(self, **kwargs):
        i = min(self.calls, len(self._responses) - 1)
        self.calls += 1
        return self._responses[i]


class _NoopTool:
    name = "noop"
    description = "no-op"
    input_schema = {"type": "object", "properties": {}}

    def run(self, args, ctx):
        return {"ok": True}


def _ctx():
    return ToolContext(tct=object(), stash=DictResultStash())


# -- CostTracker -----------------------------------------------------------------------
def test_cost_tracker_uses_corrected_pricing():
    ct = CostTracker()
    # sonnet 4.6 = $3 / $15 per 1M
    ct.add(_usage(input_tokens=1_000_000, output_tokens=1_000_000), "claude-sonnet-4-6")
    assert ct.total_usd == pytest.approx(18.0)
    assert ct.input_tokens == 1_000_000 and ct.output_tokens == 1_000_000

    ct2 = CostTracker()
    # opus 4.8 corrected to $5 / $25 per 1M
    ct2.add(_usage(input_tokens=1_000_000, output_tokens=0), "claude-opus-4-8")
    assert ct2.total_usd == pytest.approx(5.0)


def test_cost_tracker_prices_cache_tokens_at_multipliers():
    ct = CostTracker()
    # haiku $1/$5; cache read at 0.1x input, cache write at 1.25x input
    ct.add(_usage(cache_read=1_000_000, cache_create=1_000_000), "claude-haiku-4-5")
    assert ct.total_usd == pytest.approx(0.1 * 1.0 + 1.25 * 1.0)


def test_cost_tracker_handles_none_usage():
    ct = CostTracker()
    assert ct.add(None, "claude-sonnet-4-6") == 0.0
    assert ct.total_usd == 0.0


# -- AgentLoop caps --------------------------------------------------------------------
def _loop(stub, **kw):
    registry = ToolRegistry([_NoopTool()])
    return AgentLoop(stub, registry, "claude-sonnet-4-6", "system", **kw)


def test_loop_hard_stops_at_cost_cap():
    # Each call costs 0.6 (200k input @ $3/1M); cap is 0.5 -> exactly one call, then stop.
    tool_resp = _resp([_tool_block("noop", "t1", {})], "tool_use", _usage(input_tokens=200_000))
    stub = _StubLLM([tool_resp])
    loop = _loop(stub, max_iterations=99)
    ct = CostTracker()
    result = loop.run([{"role": "user", "content": "hi"}], _ctx(), ct, cost_cap_usd=0.5)

    assert stub.calls == 1
    assert result.iterations == 1
    assert ct.total_usd == pytest.approx(0.6)
    assert "spend cap" in result.stop_note


def test_loop_hard_stops_at_max_iterations():
    tool_resp = _resp([_tool_block("noop", "t1", {})], "tool_use", _usage(input_tokens=1))
    stub = _StubLLM([tool_resp])
    loop = _loop(stub, max_iterations=3)
    ct = CostTracker()
    result = loop.run([{"role": "user", "content": "hi"}], _ctx(), ct, cost_cap_usd=1000.0)

    assert result.iterations == 3
    assert "max of 3" in result.stop_note


def test_loop_ends_on_end_turn_with_text():
    end = _resp([_text_block("Here is the answer.")], "end_turn", _usage(input_tokens=10, output_tokens=5))
    stub = _StubLLM([end])
    loop = _loop(stub)
    ct = CostTracker()
    result = loop.run([{"role": "user", "content": "hi"}], _ctx(), ct, cost_cap_usd=1.0)

    assert result.stop_note is None
    assert result.final_text == "Here is the answer."
    assert result.iterations == 1


def test_anthropic_client_caches_system_prefix(monkeypatch):
    from geneset_translator.agent.agent_loop import AnthropicLLMClient

    client = AnthropicLLMClient(api_key="x")
    captured = {}

    def fake_create(**kwargs):
        captured.update(kwargs)
        return _resp([_text_block("ok")], "end_turn", _usage())

    monkeypatch.setattr(client._client.messages, "create", fake_create)
    client.create_message(
        model="claude-sonnet-4-6", system="SYS", tools=[],
        messages=[{"role": "user", "content": "hi"}], max_tokens=100,
    )
    assert isinstance(captured["system"], list)
    assert captured["system"][0]["text"] == "SYS"
    assert captured["system"][0]["cache_control"] == {"type": "ephemeral"}


class _SlowStashTool:
    name = "slow"
    description = "sleeps then stashes"
    input_schema = {"type": "object", "properties": {"i": {"type": "integer"}}}

    def run(self, args, ctx):
        import time

        time.sleep(0.2)
        rid = ctx.stash.put({"i": args.get("i")})
        return {"result_id": rid, "i": args.get("i")}


def test_concurrent_dispatch_is_parallel_ordered_and_collision_free():
    import time as _time

    from geneset_translator.agent.tct_adapter import DictResultStash

    step1 = _resp(
        [_tool_block("slow", "t0", {"i": 0}), _tool_block("slow", "t1", {"i": 1}),
         _tool_block("slow", "t2", {"i": 2})],
        "tool_use", _usage(input_tokens=1),
    )
    step2 = _resp([_text_block("done")], "end_turn", _usage(input_tokens=1))
    stub = _StubLLM([step1, step2])
    registry = ToolRegistry([_SlowStashTool()])
    loop = AgentLoop(stub, registry, "claude-sonnet-4-6", "sys", max_tool_workers=4)
    ctx = ToolContext(tct=object(), stash=DictResultStash())

    t0 = _time.time()
    result = loop.run([{"role": "user", "content": "go"}], ctx, CostTracker(), cost_cap_usd=10.0)
    elapsed = _time.time() - t0

    assert elapsed < 0.5  # 3 x 0.2s serial would exceed 0.6s; concurrent ~0.2s
    # three distinct result_ids stashed (no len()-race collision)
    assert len([k for k in ctx.stash._store if k.startswith("res_")]) == 3
    # tool_results preserved in input order
    tr_turn = next(
        m for m in result.messages
        if m["role"] == "user" and isinstance(m["content"], list)
        and m["content"] and m["content"][0].get("type") == "tool_result"
    )
    assert [b["tool_use_id"] for b in tr_turn["content"]] == ["t0", "t1", "t2"]
    # latest render target is the LAST input-order tool (i=2), deterministically
    assert ctx.stash.get(ctx.latest_result_id)["i"] == 2


class _RecorderTool:
    def __init__(self, name, parallel_safe, log):
        self.name = name
        self.parallel_safe = parallel_safe
        self.description = "recorder"
        self.input_schema = {"type": "object", "properties": {}}
        self._log = log

    def run(self, args, ctx):
        self._log.append(self.name)
        return {"ok": self.name}


def test_mixed_step_runs_sequentially_in_input_order():
    from geneset_translator.agent.tct_adapter import DictResultStash

    order = []
    registry = ToolRegistry(
        [
            _RecorderTool("a", True, order),
            _RecorderTool("b", False, order),  # not parallel_safe -> whole step is sequential
            _RecorderTool("c", True, order),
        ]
    )
    step1 = _resp(
        [_tool_block("a", "t0", {}), _tool_block("b", "t1", {}), _tool_block("c", "t2", {})],
        "tool_use", _usage(input_tokens=1),
    )
    step2 = _resp([_text_block("done")], "end_turn", _usage(input_tokens=1))
    loop = AgentLoop(_StubLLM([step1, step2]), registry, "claude-sonnet-4-6", "sys")
    ctx = ToolContext(tct=object(), stash=DictResultStash())
    loop.run([{"role": "user", "content": "go"}], ctx, CostTracker(), cost_cap_usd=10.0)
    assert order == ["a", "b", "c"]  # deterministic input order, not completion order


def test_system_prompt_uses_disease_label_not_curie():
    from geneset_translator.agent.agent_loop import build_system_prompt

    p = build_system_prompt(["FLT3"], disease="MONDO:0005361", disease_label="Eosinophilic Esophagitis")
    assert "Eosinophilic Esophagitis" in p and "MONDO:0005361" in p
    p2 = build_system_prompt(["FLT3"], disease="MONDO:0005361")  # bare CURIE, no label
    assert "do not guess" in p2.lower()


def test_loop_appends_tool_results_then_finishes():
    # First a tool call, then end_turn -> the user turn with tool_result must be in messages.
    r1 = _resp([_tool_block("noop", "t1", {})], "tool_use", _usage(input_tokens=10))
    r2 = _resp([_text_block("done")], "end_turn", _usage(input_tokens=10, output_tokens=2))
    stub = _StubLLM([r1, r2])
    loop = _loop(stub)
    ctx = _ctx()
    result = loop.run([{"role": "user", "content": "hi"}], ctx, CostTracker(), cost_cap_usd=1.0)

    assert result.final_text == "done"
    # one successful tool call logged
    assert ctx.tool_log == [{"name": "noop", "args": {}}]
    tool_result_turns = [
        m for m in result.messages
        if m["role"] == "user" and isinstance(m["content"], list)
        and m["content"] and m["content"][0].get("type") == "tool_result"
    ]
    assert len(tool_result_turns) == 1
    assert tool_result_turns[0]["content"][0]["tool_use_id"] == "t1"
