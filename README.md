# BioGraph Explorer 🧬

Streamlit application for multi-gene TRAPI query integration with NetworkX clustering and LLM-assisted exploration using Cytoscape.js visualization.

## Overview

BioGraph Explorer queries the NCATS Translator knowledge graph system to find connections between genes and diseases, then uses network analysis to identify convergent pathways and therapeutic targets.

**Features:**
- 🔍 Gene normalization using TCT library
- 🌐 Parallel TRAPI queries across 15+ Translator APIs
- 📊 NetworkX graph construction with rich node attributes
- 🎯 Louvain community detection
- 📈 Centrality analysis (PageRank, betweenness, degree)
- 🎨 Interactive Cytoscape.js visualization with multiple layout algorithms
- 💾 CSV import/export for gene lists
- 🔬 Material Icons for biological entity types

## Installation

### Prerequisites
- Python >=3.11,<4.0
- [Poetry](https://python-poetry.org/docs/#installation) for dependency management

### Install Dependencies

```bash
# Clone repository
git clone <repository-url>
cd biograph_explorer

# Install with Poetry
poetry install

# Activate environment
poetry shell
```

## Quick Start

### 1. Run the Streamlit App

```bash
streamlit run app.py
```

The app will open in your browser at `http://localhost:8501`

### 2. Use Example Datasets

**Option A: Via UI**
1. In the sidebar, select "Example Dataset"
2. Choose either:
   - **Alzheimer's Disease (15 genes)**: APOE, APP, PSEN1, PSEN2, MAPT, TREM2, etc.
   - **COVID-19 (10 genes)**: CD6, IFITM3, IFITM2, STAT5A, KLRG1, DPP4, etc.
3. Click **Run Query**

**Option B: Load from CSV**
```bash
# Example files included
data/test_genes/alzheimers_genes.csv
data/test_genes/covid19_genes.csv
```

### 3. Upload Your Own Data

Create a CSV file with a `gene_symbol` column:

```csv
gene_symbol,description
APOE,Lipid transport
APP,Amyloid precursor
PSEN1,γ-secretase component
```

Then upload via the sidebar.

## Workflow

1. **Input** → Gene symbols (manual, CSV, or example)
2. **Normalize** → TCT converts symbols to CURIEs (e.g., APOE → NCBIGene:348)
3. **Query** → Parallel TRAPI queries to 15 Translator APIs
4. **Build** → NetworkX graph with 400-800 nodes, 700+ edges
5. **Analyze** → Louvain communities, PageRank, convergence metrics
6. **Visualize** → Interactive results dashboard

## Expected Results

### Alzheimer's Test Case (15 genes)
- **Edges**: ~1,500 from 6-10 successful APIs
- **Nodes**: ~500
- **Communities**: 3-4 (amyloid processing, lipid metabolism, neuroinflammation)
- **Convergent nodes**: BACE1, APOE, amyloid-β
- **Top targets**: BACE1, APOE, TREM2

### COVID-19 Test Case (10 genes)
- **Edges**: ~800
- **Nodes**: ~400
- **Communities**: 2-3
- **API success rate**: 40-60% (normal)

## Project Structure

```
biograph_explorer/
├── app.py                       # Main Streamlit application
├── data/
│   ├── test_genes/
│   │   ├── alzheimers_genes.csv   # Alzheimer's test dataset
│   │   └── covid19_genes.csv      # COVID-19 test dataset
│   ├── cache/                     # TRAPI response cache
│   └── sessions/                  # Saved analysis sessions
├── src/biograph_explorer/
│   ├── core/
│   │   ├── trapi_client.py        # TRAPI query client (TCT integration)
│   │   ├── graph_builder.py       # NetworkX graph construction
│   │   ├── clustering_engine.py   # Louvain + centrality analysis
│   │   └── rag_system.py          # LLM integration (Phase 3)
│   ├── utils/
│   │   ├── validators.py          # Input validation
│   │   ├── formatters.py          # Data formatting
│   │   └── persistence.py         # Graph serialization
│   └── ui/
│       ├── input_panel.py         # Gene/disease input (Streamlit)
│       ├── query_status.py        # Progress tracking
│       ├── results_overview.py    # Dashboard
│       ├── network_viz.py         # Cytoscape.js visualization
│       └── convergence_view.py    # Convergent nodes table
└── tests/
    └── fixtures/
        └── alzheimers_test_case.json  # Test case specification
```

## Development

### Run Tests

```bash
poetry run pytest
poetry run pytest --cov=src/biograph_explorer
```

### Build Documentation

```bash
cd docs
make html
```

### Add Dependencies

```bash
poetry add <package-name>
```

## Key Technologies

- **[TCT (Translator Clinical Tools)](https://github.com/NCATSTranslator/Translator-All/wiki/TCT)**: TRAPI query library
- **[NetworkX](https://networkx.org/)**: Graph analysis
- **[python-louvain](https://github.com/taynaud/python-louvain)**: Community detection
- **[Streamlit](https://streamlit.io/)**: Web UI
- **[st-link-analysis](https://github.com/AlrasheedA/st-link-analysis)**: Interactive Cytoscape.js visualization
- **[Pydantic](https://docs.pydantic.dev/)**: Data validation

## API Notes

### TRAPI Query Behavior

- **Normal API success rate**: 40-60% (6-10 out of 15 APIs)
- **Failed APIs**: `NoneType is not iterable` = no data for this query (expected)
- **Timeout handling**: Automatic fallback if Plover APIs unreachable
- **Caching**: All responses cached to `data/cache/` for reuse

### Gene Normalization

- Uses TCT's `name_resolver.batch_lookup()`
- Accepts HUGO symbols (APOE, APP, etc.)
- Returns NCBIGene CURIEs
- Success rate: 90-100% for valid human genes

## Troubleshooting

### TCT Not Available
```bash
poetry add tct
```

### python-louvain Not Found
```bash
poetry add python-louvain
```

### TRAPI Queries Timing Out
- **Normal behavior**: Some APIs timeout (expected)
- **Check**: At least 5-6 APIs should succeed
- **Fallback**: Uses cached responses if available

### Empty Graph
- Check gene normalization succeeded
- Verify disease CURIE format (MONDO:, DOID:, etc.)
- Try different gene set (some genes have sparse data)

## Performance

| Operation | Target | Notes |
|-----------|--------|-------|
| Gene normalization | <5s | 15 genes |
| TRAPI batch query | <2min | With caching/parallelization |
| Graph construction | <10s | 500 nodes |
| Community detection | <5s | Louvain algorithm |
| Centrality metrics | <10s | PageRank + betweenness |

## Phase 2 Status ✅

**Implemented:**
- ✅ TRAPI client with parallel queries
- ✅ NetworkX graph builder
- ✅ Louvain community detection
- ✅ Centrality metrics (PageRank, betweenness, degree)
- ✅ Streamlit UI with CSV import
- ✅ Example datasets (Alzheimer's, COVID-19)
- ✅ Progress tracking
- ✅ Session state management

**Phase 3 (In Progress):**
- ✅ Cytoscape.js network visualization with st-link-analysis
- ✅ Multiple layout algorithms (cose, fcose, circle, grid, etc.)
- ✅ Material Icons for entity types
- 🔲 Claude RAG system integration
- 🔲 Citation validation with subgraph extraction
- 🔲 Interactive chat interface

## Contributing

Interested in contributing? Check out the [contributing guidelines](CONTRIBUTING.md). Please note that this project is released with a [Code of Conduct](CONDUCT.md). By contributing to this project, you agree to abide by its terms.

## License

`biograph_explorer` was created by Natalie Gill. It is licensed under the terms of the MIT license.

## Credits

`biograph_explorer` was created with [`cookiecutter`](https://cookiecutter.readthedocs.io/en/latest/) and the `py-pkgs-cookiecutter` [template](https://github.com/py-pkgs/py-pkgs-cookiecutter).

## Citation

If you use BioGraph Explorer in your research, please cite:

```
@software{biograph_explorer,
  author = {Gill, Natalie},
  title = {BioGraph Explorer: Multi-gene TRAPI Query Integration with Network Analysis},
  year = {2025},
  url = {https://github.com/yourusername/biograph_explorer}
}
```
