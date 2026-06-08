"""GeneSet Translator - navigation shell.

Single Streamlit entry point. Routes between the AI Chat Explorer (default when an
Anthropic API key is present) and the Classic Explorer (always available; the only
page when no key is set). This shell owns ``st.set_page_config``; the page modules
must not call it.
"""

import os

import streamlit as st
from dotenv import load_dotenv

load_dotenv()

# Console + rotating file logging at data/logs/chat_agent.log. Idempotent across reruns; the agent
# path logs its progress to both the terminal and that file (tail -f data/logs/chat_agent.log).
from geneset_translator.utils.logging_setup import configure_logging

configure_logging()

st.set_page_config(
    page_title="GeneSet Translator",
    page_icon=":material/biotech:",
    layout="wide",
    initial_sidebar_state="expanded",
)


def _has_api_key() -> bool:
    """True if an Anthropic API key is available via env or settings."""
    if os.environ.get("ANTHROPIC_API_KEY"):
        return True
    try:
        from geneset_translator.config.settings import get_settings

        return bool(get_settings().claude_api_key)
    except Exception:
        return False


def _chat_page() -> None:
    from geneset_translator.ui.chat_page import render

    render()


def _enable_chat_page() -> None:
    from geneset_translator.ui.no_key_help import render

    render()


_CLASSIC_PATH = "src/geneset_translator/ui/classic_explorer.py"

if _has_api_key():
    pages = [
        st.Page(
            _chat_page,
            title="Chat Explorer",
            icon=":material/forum:",
            url_path="chat",
            default=True,
        ),
        st.Page(
            _CLASSIC_PATH,
            title="Classic Explorer",
            icon=":material/hub:",
            url_path="classic",
        ),
    ]
else:
    pages = [
        st.Page(
            _CLASSIC_PATH,
            title="Classic Explorer",
            icon=":material/hub:",
            url_path="classic",
            default=True,
        ),
        st.Page(
            _enable_chat_page,
            title="Enable AI Chat",
            icon=":material/key:",
            url_path="enable_chat",
        ),
    ]

st.navigation(pages).run()
