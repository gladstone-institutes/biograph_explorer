"""Claude model discovery and pricing helpers.

Kept out of ``LLMSummarizer`` because model discovery runs in the sidebar before
any summarizer exists and must degrade gracefully when the API key is missing or
expired. Streamlit-level caching is applied at the call site, not here, so these
functions stay pure and testable.
"""

import logging
from typing import List, Dict, Optional, Tuple

logger = logging.getLogger(__name__)

# Default model selected in the UI when available.
DEFAULT_MODEL_ID = "claude-haiku-4-5"

# Known-good fallback list shown when the API key is missing/expired so the model
# selector still works. Ordered best-default-first. Keep roughly in sync with the
# models the project supports.
FALLBACK_MODELS: List[Dict[str, str]] = [
    {"id": "claude-haiku-4-5", "display_name": "Claude Haiku 4.5"},
    {"id": "claude-sonnet-4-5", "display_name": "Claude Sonnet 4.5"},
    {"id": "claude-opus-4-1", "display_name": "Claude Opus 4.1"},
]

# Pricing (USD per 1M tokens) keyed by a substring of the model id. Pricing is not
# returned by the models API, so a local map is required. Matched by longest key so
# a specific tier wins. Values are (input_per_million, output_per_million).
MODEL_PRICING: Dict[str, Tuple[float, float]] = {
    "haiku": (1.00, 5.00),
    "sonnet": (3.00, 15.00),
    "opus": (15.00, 75.00),
}

# Fallback pricing (Haiku-class) for an unknown model id; never raise on costing.
DEFAULT_PRICING: Tuple[float, float] = (1.00, 5.00)


def get_model_pricing(model_id: str) -> Tuple[float, float]:
    """Return ``(input_per_million, output_per_million)`` for a model id.

    Matches the longest pricing key contained in the id (e.g. "haiku" in
    "claude-haiku-4-5"). Falls back to ``DEFAULT_PRICING`` for unknown ids and
    never raises.
    """
    mid = (model_id or "").lower()
    best_key = None
    for key in MODEL_PRICING:
        if key in mid and (best_key is None or len(key) > len(best_key)):
            best_key = key
    return MODEL_PRICING[best_key] if best_key else DEFAULT_PRICING


def fetch_available_models(api_key: Optional[str]) -> List[Dict[str, str]]:
    """Return available Claude models as ``[{"id", "display_name"}, ...]``.

    Calls ``client.models.list()`` when an API key is present; on any failure
    (missing/expired key, network error, SDK change) logs a warning and returns
    ``FALLBACK_MODELS``. The default model is sorted first when present.
    """
    if not api_key:
        return list(FALLBACK_MODELS)

    try:
        from anthropic import Anthropic

        client = Anthropic(api_key=api_key)
        models: List[Dict[str, str]] = []
        for m in client.models.list(limit=100):
            model_id = getattr(m, "id", None)
            if not model_id or not model_id.startswith("claude-"):
                continue
            display = getattr(m, "display_name", None) or model_id
            models.append({"id": model_id, "display_name": display})

        if not models:
            logger.warning("models.list() returned no Claude models; using fallback list")
            return list(FALLBACK_MODELS)

        # Default model first, then preserve API order (most recent first).
        models.sort(key=lambda d: d["id"] != DEFAULT_MODEL_ID)
        return models
    except Exception as e:  # noqa: BLE001 - any failure must fall back gracefully
        logger.warning(f"Could not fetch Claude models ({e}); using fallback list")
        return list(FALLBACK_MODELS)
