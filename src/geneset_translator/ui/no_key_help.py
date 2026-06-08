"""Help page shown (in place of the Chat Explorer) when no Anthropic API key is set."""

import streamlit as st


def render() -> None:
    st.title(":material/key: Enable AI Chat")
    st.markdown(
        """
The **Chat Explorer** lets you ask biological questions about your gene set in plain
English and have a Claude agent run the relevant Translator (TCT) queries for you.

It needs an Anthropic API key, which is not currently set. To enable it:

1. Get a key from the [Anthropic Console](https://console.anthropic.com/).
2. Add it to a `.env` file in the project root:

   ```
   ANTHROPIC_API_KEY="sk-ant-..."
   ```

3. Restart the app.

Until then, use the **Classic Explorer** (in the sidebar) to run gene-neighborhood
queries with the form-based interface.
        """
    )
    st.info(
        "No key required for the Classic Explorer. The AI Chat Explorer will appear "
        "here automatically once `ANTHROPIC_API_KEY` is set.",
        icon=":material/info:",
    )
