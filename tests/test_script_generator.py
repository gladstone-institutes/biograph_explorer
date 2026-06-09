"""Tests for agent.script_generator: emit a runnable, TCT-only reproduction script that builds
NetworkX graphs and exports them to Cytoscape.js JSON."""

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
    assert "import networkx as nx" in script
    assert "TranslatorResources.load()" in script


def test_script_contains_each_finder_call_and_curies():
    script = generate_reproduction_script(_log(), questions=["What drugs target FLT3?"])
    assert "name_resolver.batch_lookup(strings=['FLT3', 'BCL2']" in script
    assert "only_taxa='NCBITaxon:9606'" in script
    assert "TCT.Neighborhood_finder(input_node='NCBIGene:2322'" in script
    assert "TCT.Path_finder(input_node1='NCBIGene:596'" in script
    assert "run_gene_network(['NCBIGene:2322', 'NCBIGene:596'])" in script
    assert "translator_query.parallel_api_query(" in script  # inside the gene_network helper
    # question echoed into the header
    assert "What drugs target FLT3?" in script
    # app-only step noted, not executed
    assert "app-only step not reproduced" in script and "edge_evidence" in script


def test_script_builds_and_exports_networkx_graphs():
    """Every finder result is converted to a NetworkX graph and exported to Cytoscape.js JSON."""
    script = generate_reproduction_script(_log())
    assert "def to_graph(result):" in script
    assert "def export_cyjs(graph, path):" in script
    assert "nx.cytoscape_data(graph)" in script
    assert "to_networkx(resolve_names=True, include_attributes=True)" in script
    # header makes clear the export is the full graph, not the sampled in-app view
    assert "COMPLETE graph" in script and "not capped" in script
    # each finder exports a .cyjs under the output dir
    assert "export_cyjs(to_graph(nb1), OUT / 'nb1_NCBIGene_2322.cyjs')" in script
    assert "export_cyjs(to_graph(pr1), OUT / 'pr1.cyjs')" in script
    assert "export_cyjs(to_graph(kg1), OUT / 'kg1_network.cyjs')" in script


def test_duplicate_finder_calls_are_deduplicated():
    """Repeated identical finder calls (common across turns) run once; later ones are skipped."""
    log = [
        {"name": "gene_neighborhood", "args": {"gene_curie": "NCBIGene:2322"}},
        {"name": "gene_neighborhood", "args": {"gene_curie": "NCBIGene:2322"}},  # exact dup
        {"name": "gene_neighborhood", "args": {"gene_curie": "NCBIGene:596"}},   # different
    ]
    script = generate_reproduction_script(log)
    assert script.count("TCT.Neighborhood_finder(input_node='NCBIGene:2322'") == 1
    assert "duplicate of nb1; skipped" in script
    assert "TCT.Neighborhood_finder(input_node='NCBIGene:596'" in script


def test_resolve_genes_is_faithful_to_gateway():
    """Gene resolution uses the taxon -> taxon-free two-pass; a typed disease/drug drops the taxon
    filter and constrains by biolink_type (mirrors TctGateway.resolve_genes)."""
    gene_script = generate_reproduction_script(
        [{"name": "resolve_genes", "args": {"names": ["FLT3"]}}]
    )
    assert "only_taxa='NCBITaxon:9606'" in gene_script
    assert "if _unresolved:" in gene_script  # two-pass fallback

    disease_script = generate_reproduction_script(
        [{"name": "resolve_genes", "args": {"names": ["acute myeloid leukemia"], "biolink_type": "biolink:Disease"}}]
    )
    assert "biolink_type='biolink:Disease'" in disease_script
    assert "only_taxa" not in disease_script  # typed non-gene resolution is taxon-free


def test_display_only_steps_are_noted():
    log = [
        {"name": "filter_graph", "args": {"top_n": 12}},
        {"name": "add_disease_node", "args": {"disease_curie": "MONDO:0100096"}},
        {"name": "show_result", "args": {"result_id": "res_1", "merge": False}},
    ]
    script = generate_reproduction_script(log)
    assert "display-only step" in script
    assert "filter_graph" in script and "add_disease_node" in script and "show_result" in script


def test_generated_script_is_valid_python():
    script = generate_reproduction_script(_log())
    compile(script, "<reproduction>", "exec")  # raises SyntaxError if malformed


def test_empty_log_still_valid():
    script = generate_reproduction_script([])
    compile(script, "<reproduction>", "exec")
    assert "TranslatorResources.load()" in script
