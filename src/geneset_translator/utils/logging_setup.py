"""Idempotent logging setup for the app.

Adds a console handler plus a rotating file handler at ``data/logs/chat_agent.log`` so the agent's
progress (see ``agent/agent_loop.py`` and ``agent/tct_adapter.py``) is captured both in the
terminal and in a persistent file you can ``tail -f``. Safe to call on every Streamlit rerun: a
sentinel on the root logger prevents duplicate handlers.
"""

from __future__ import annotations

import logging
from logging.handlers import RotatingFileHandler
from pathlib import Path
from typing import Optional, Union

_FMT = logging.Formatter("%(asctime)s %(levelname)s %(name)s: %(message)s", datefmt="%H:%M:%S")
_QUIET_LIBS = ("httpx", "httpcore", "urllib3", "anthropic")
_LOG_FILENAME = "chat_agent.log"
_SENTINEL = "_gst_handlers_added"


def get_log_file(log_dir: Optional[Union[str, Path]] = None) -> Path:
    """Resolve the agent log file path (no side effects)."""
    if log_dir is None:
        from geneset_translator.config.settings import get_settings

        log_dir = get_settings().logs_dir
    return (Path(log_dir) / _LOG_FILENAME).resolve()


def configure_logging(log_dir: Optional[Union[str, Path]] = None) -> Path:
    """Configure console + rotating-file logging at INFO. Idempotent. Returns the log file path."""
    root = logging.getLogger()
    root.setLevel(logging.INFO)

    if log_dir is None:
        from geneset_translator.config.settings import get_settings

        log_dir = get_settings().logs_dir
    log_dir = Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    logfile = (log_dir / _LOG_FILENAME).resolve()

    if not getattr(root, _SENTINEL, False):
        console = logging.StreamHandler()
        console.setFormatter(_FMT)
        root.addHandler(console)

        file_handler = RotatingFileHandler(
            str(logfile), maxBytes=5_000_000, backupCount=3, encoding="utf-8"
        )
        file_handler.setFormatter(_FMT)
        root.addHandler(file_handler)

        setattr(root, _SENTINEL, True)
        logging.getLogger(__name__).info("Agent logging to %s", logfile)

    for name in _QUIET_LIBS:
        logging.getLogger(name).setLevel(logging.WARNING)

    return logfile
