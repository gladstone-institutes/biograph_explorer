"""Tests for utils.logging_setup: path resolution + idempotent handler setup + file capture."""

import logging

from geneset_translator.utils import logging_setup


def test_get_log_file_resolves_under_dir(tmp_path):
    assert logging_setup.get_log_file(tmp_path) == (tmp_path / "chat_agent.log").resolve()


def test_configure_logging_idempotent_and_writes(tmp_path):
    root = logging.getLogger()
    saved_handlers = list(root.handlers)
    had_sentinel = hasattr(root, "_gst_handlers_added")
    if had_sentinel:
        delattr(root, "_gst_handlers_added")
    try:
        p1 = logging_setup.configure_logging(tmp_path)
        n_after_first = len(root.handlers)
        p2 = logging_setup.configure_logging(tmp_path)  # second call must be a no-op
        n_after_second = len(root.handlers)

        assert p1 == p2 == (tmp_path / "chat_agent.log").resolve()
        assert n_after_second == n_after_first  # no duplicate handlers
        assert n_after_first >= len(saved_handlers) + 2  # console + rotating file added

        logging.getLogger("geneset_translator.test").info("hello-log-marker")
        for h in root.handlers:
            h.flush()
        assert "hello-log-marker" in p1.read_text(encoding="utf-8")
    finally:
        # Restore global root-logger state so we don't pollute other tests.
        for h in list(root.handlers):
            if h not in saved_handlers:
                h.close()
                root.removeHandler(h)
        if not had_sentinel and hasattr(root, "_gst_handlers_added"):
            delattr(root, "_gst_handlers_added")
