"""Tests for utils.formatters, focused on strip_emoji (the emoji ban helper)."""

from geneset_translator.utils.formatters import strip_emoji


def test_strip_emoji_removes_common_emoji():
    assert strip_emoji("Yes ✅ done ❌") == "Yes done"  # check + cross
    assert "\U0001f9ec" not in strip_emoji("Genes \U0001f9ec are interesting")  # dna emoji
    assert strip_emoji("Top hit \U0001f48a midostaurin") == "Top hit midostaurin"  # pill
    assert strip_emoji("warning ⚠️ here") == "warning here"  # warning + variation selector


def test_strip_emoji_preserves_scientific_and_normal_text():
    # Greek letters used in gene/protein names must survive.
    assert strip_emoji("TNF-α and IFN-γ signaling") == "TNF-α and IFN-γ signaling"
    assert strip_emoji("p ≤ 0.05, µM dose") == "p ≤ 0.05, µM dose"
    # Plain markdown / identifiers untouched.
    text = "**FLT3** affects [PMID:28644114](https://pubmed.ncbi.nlm.nih.gov/28644114)."
    assert strip_emoji(text) == text


def test_strip_emoji_handles_empty_and_none():
    assert strip_emoji("") == ""
    assert strip_emoji(None) is None


def test_strip_emoji_tidies_whitespace_left_behind():
    assert strip_emoji("drug \U0001f48a , gene") == "drug, gene"
