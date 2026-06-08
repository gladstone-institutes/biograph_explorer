"""Guards for the chat page: pandas import + example datasets load cleanly.

Regression: the example loader used pd.read_csv but chat_page didn't import pandas; the broad
try/except hid it as an st.error, so it must be covered by an explicit test.
"""

import pandas as pd

from geneset_translator.ui import chat_page


def test_chat_page_imports_pandas():
    assert hasattr(chat_page, "pd"), "chat_page must import pandas (used by the example/CSV loaders)"


def test_example_datasets_load_and_have_gene_symbol_column():
    assert chat_page.EXAMPLES, "expected at least one example dataset"
    for _label, (path, disease, disease_name) in chat_page.EXAMPLES.items():
        df = pd.read_csv(path)
        assert "gene_symbol" in df.columns
        symbols = df["gene_symbol"].dropna().astype(str).str.strip().tolist()
        assert len(symbols) > 0
        assert disease.startswith("MONDO:")
        assert disease_name  # a human-readable name for the system prompt


def test_looks_like_curie():
    assert chat_page._looks_like_curie("MONDO:0005361")
    assert chat_page._looks_like_curie("NCBIGene:2322")
    assert not chat_page._looks_like_curie("acute myeloid leukemia")
    assert not chat_page._looks_like_curie("")
