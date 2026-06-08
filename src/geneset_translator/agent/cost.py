"""Actual-usage cost accounting for the agent loop.

No estimation: every dollar figure is derived from the actual ``response.usage`` returned by
each API call. Token prices come from the single ``model_utils.MODEL_PRICING`` table (the
Anthropic API does not return pricing); cache tokens are priced at the documented multipliers
(reads ~0.1x input, writes ~1.25x input).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List

from geneset_translator.utils.model_utils import get_model_pricing

CACHE_READ_MULTIPLIER = 0.1
CACHE_WRITE_MULTIPLIER = 1.25


@dataclass
class CostTracker:
    """Accumulates actual token usage and the derived dollar cost across a session."""

    input_tokens: int = 0
    output_tokens: int = 0
    cache_read_tokens: int = 0
    cache_creation_tokens: int = 0
    total_usd: float = 0.0
    calls: List[Dict[str, Any]] = field(default_factory=list)

    def add(self, usage: Any, model: str) -> float:
        """Record one call's actual usage; returns the dollar cost of that call."""
        if usage is None:
            return 0.0
        in_tok = int(getattr(usage, "input_tokens", 0) or 0)
        out_tok = int(getattr(usage, "output_tokens", 0) or 0)
        cache_read = int(getattr(usage, "cache_read_input_tokens", 0) or 0)
        cache_create = int(getattr(usage, "cache_creation_input_tokens", 0) or 0)

        in_rate, out_rate = get_model_pricing(model)
        usd = (
            in_tok / 1_000_000 * in_rate
            + out_tok / 1_000_000 * out_rate
            + cache_read / 1_000_000 * in_rate * CACHE_READ_MULTIPLIER
            + cache_create / 1_000_000 * in_rate * CACHE_WRITE_MULTIPLIER
        )

        self.input_tokens += in_tok
        self.output_tokens += out_tok
        self.cache_read_tokens += cache_read
        self.cache_creation_tokens += cache_create
        self.total_usd += usd
        self.calls.append(
            {
                "model": model,
                "input_tokens": in_tok,
                "output_tokens": out_tok,
                "cache_read_tokens": cache_read,
                "cache_creation_tokens": cache_create,
                "usd": usd,
            }
        )
        return usd

    @property
    def total_tokens(self) -> int:
        return (
            self.input_tokens
            + self.output_tokens
            + self.cache_read_tokens
            + self.cache_creation_tokens
        )

    def summary(self) -> Dict[str, Any]:
        return {
            "input_tokens": self.input_tokens,
            "output_tokens": self.output_tokens,
            "cache_read_tokens": self.cache_read_tokens,
            "cache_creation_tokens": self.cache_creation_tokens,
            "total_tokens": self.total_tokens,
            "total_usd": round(self.total_usd, 6),
            "n_calls": len(self.calls),
        }
