"""Tests for agent.script_generator: emit a runnable TCT reproduction script."""

from geneset_translator.agent.script_generator import generate_reproduction_script


def _log():
    return [
        {"name": "resolve_genes", "args": {"names": ["FLT3", "BCL2"]}},
        {"name": "gene_neighborhood", "args": {"gene_curie": "NCBIGene:2322"}},
        {
            "name": "path_between",
            "args": {"node1_curie": "NCBIGene:596", "node2_curie": "CHEBI:1", "intermediate_categories": ["biolink:Gene"]},
        },
        {"name": "gene_network", "args": {"gene_curies": ["NCBIGene:2322", "NCBIGene:596"]}},
        {"name": "edge_evidence", "args": {"result_id": "res_1", "subject": "A", "object": "B"}},
    ]


def test_script_has_imports_and_resource_load():
    script = generate_reproduction_script(_log())
    assert "from TCT.translator_resources import TranslatorResources" in script
    assert "from TCT import TCT, name_resolver, translator_query" in script
    assert "TranslatorResources.load()" in script


def test_script_contains_each_finder_call_and_curies():
    script = generate_reproduction_script(_log(), questions=["What drugs target FLT3?"])
    assert "name_resolver.batch_lookup(strings=['FLT3', 'BCL2']" in script
    assert "only_taxa='NCBITaxon:9606'" in script
    assert "TCT.Neighborhood_finder(input_node='NCBIGene:2322'" in script
    assert "TCT.Path_finder(input_node1='NCBIGene:596'" in script
    assert "translator_query.parallel_api_query(" in script
    # question echoed into the header
    assert "What drugs target FLT3?" in script
    # non-finder step noted, not executed
    assert "skipped non-finder step: edge_evidence" in script


def test_generated_script_is_valid_python():
    script = generate_reproduction_script(_log())
    compile(script, "<reproduction>", "exec")  # raises SyntaxError if malformed


def test_empty_log_still_valid():
    script = generate_reproduction_script([])
    compile(script, "<reproduction>", "exec")
    assert "TranslatorResources.load()" in script
