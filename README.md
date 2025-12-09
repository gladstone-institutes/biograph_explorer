# BioGraph Explorer 🧬

**Status**: Phase 2 Complete ✅ | Phase 3 In Progress 📋

Streamlit application for exploring biomedical knowledge graphs through multi-gene TRAPI queries, network clustering, and interactive visualization.

## Overview

BioGraph Explorer is a research tool for exploring different approaches to high-dimensional biomedical queries using the NCATS Translator knowledge graph system:

1. **TRAPI Query Execution**: Executes multi-gene pathfinder queries from NCATS Translator APIs using TCT (Translator Component Toolkit)
2. **Result Augmentation**: Enriches intermediate nodes with additional information from external data sources beyond what Translator returns (not implemented yet)
3. **Result Distillation**: Applies clustering and other methods to identify patterns and reduce complexity
4. **LLM-Assisted Exploration**: Uses AI summarization and natural language Q&A to extract insights from augmented query results (not implemented yet)

The goal is to discover effective workflows for transforming complex high dimensional Translator queries into interpretable biological insights for results from transcriptomics experiments.

## Key Technologies

- **[TCT (Translator Component Toolkit)](https://github.com/gloriachin/Translator_component_toolkit)**: TRAPI query library
- **[NetworkX](https://networkx.org/)**: Graph analysis
- **[python-louvain](https://github.com/taynaud/python-louvain)**: Community detection
- **[Streamlit](https://streamlit.io/)**: Web UI
- **[streamlit-cytoscape](https://github.com/natalie-23-gill/streamlit-cytoscape)**: Interactive Cytoscape.js visualization based on [st-link-analysis](https://github.com/AlrasheedA/st-link-analysis)
- **[Pydantic](https://docs.pydantic.dev/)**: Data validation

**Current Features (Phase 1 and 2):**
- ✅ Gene normalization using NCATS Translator
- ✅ Parallel TRAPI queries across 15+ Translator APIs
- ✅ NetworkX graph construction
- ✅ Louvain community detection & centrality analysis
- ✅ Interactive Cytoscape.js visualization with 9 layout algorithms
- ✅ CSV import and example datasets (Alzheimer's, COVID-19)
- ✅ Query caching for development

**Planned Features (Phase 3):**
- 📋 Highlight the query genes in the resulting networks
- 📋 Node resizing in network view - currently broken
- 📋 Better clustering algorithms - current louvain clustering is not ideal for a highly constrained network (only 1 disease node)
- 📋 Filter results by predicate types
- 📋 Session management (explicitly save/load sessions, currently just cached in a local directory)
- 📋 Built in disease curie lookup
- 📋 Export to Cytoscape button
- 📋 Improved query builder (3-hop queries, gene -> [intermediate] -> disease associated phenotypes)
- 📋 Augmenting the query results with external data sources
- 📋 LLM cluster summarization
- 📋 RAG-powered chat interface with Claude AI with knowledge provenance subgraph display


> **Note**: This is research software in active development. Some features are incomplete and subject to change.

## Installation

### Prerequisites
- Python >=3.11,<4.0
- [Poetry](https://python-poetry.org/docs/#installation) for dependency management

### Install Dependencies

```bash
# Clone repository
git clone https://github.com/gladstone-institutes/biograph_explorer.git
cd biograph_explorer

# Install with Poetry
poetry install

```

## Quick Start

### 1. Run the App

```bash
streamlit run app.py
```


![BioGraph Explorer Main Interface](static/images/main_page.png)   


### 2. Try an Example Query

1. In the sidebar, select **"Example Dataset"**
2. Choose **"COVID-19 (10 genes)"** (default)
3. Select query pattern
4. Select intermediate node types
5. Click **"Run Query"** (takes ~3-5 minutes)

### 3. Explore Results

The app displays results in three tabs:

- **Overview**: Graph statistics, node categories, top convergent nodes
- **Network**: Interactive Cytoscape.js visualization (drag to pan, scroll to zoom)
- **Communities**: Detected biological modules with hub nodes

### 4. Try Your Own Genes

- Create a CSV with a `gene_symbol` column (see [data/test_genes/covid19_genes.csv](data/test_genes/covid19_genes.csv) for format) or enter genes manually in the sidebar. 

## Understanding Your Results

**Query Patterns:**
- **1-hop** (`Gene → Any`): Finds all connections to your genes (broad discovery)
- **2-hop** (`Gene → Intermediate → Disease`): Finds intermediate nodes of slected types between genes and a disease (requires disease CURIE like `MONDO:0100096` for COVID-19)

**Key Metrics:**
- **Gene Frequency**: How many query genes connect to each node (convergence indicator)
- **PageRank**: Relative importance in the network
- **Communities**: Biological modules detected by Louvain clustering


## Current Limitations

- **Large graphs**: Network nodes are auto-sampled by top intermediates with the most connections to query genes for visualization performance
- **Session persistence**: Manual query re-run required or loading cached queries (better session management coming in Phase 3)


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


## Current Limitations

- **Large graphs**: Networks >200 nodes are auto-sampled for visualization performance
- **Session persistence**: Manual query re-run required or loading cached queries (auto-save coming in Phase 3)

## Troubleshooting

**"Query returned no results"**
- Normal for some API combinations - at least 5-6 APIs should succeed
- Check gene normalization succeeded (view progress messages)
- Try different genes (some have sparse data in Translator)

**"Graph is empty"**
- Verify disease CURIE format for 2-hop queries (`MONDO:0100096` for COVID-19, `MONDO:0004975` for Alzheimer's)
- Check that genes are valid HUGO symbols (uppercase)

**"Visualization is slow"**
- Use simpler layout algorithms (dagre, breadthfirst vs. fcose)
- Reduce `max_intermediates` slider to sample fewer nodes
- View individual communities instead of full graph


## Acknowledgements   

`biograph_explorer` was created with [`cookiecutter`](https://cookiecutter.readthedocs.io/en/latest/) and the `py-pkgs-cookiecutter` [template](https://github.com/py-pkgs/py-pkgs-cookiecutter). The streamlit theme was created by github user [jmedia65](https://github.com/jmedia65/awesome-streamlit-themes).

## AI Disclosure Statement

Generative AI tools (Claude Code, Anthropic) were used as coding assistants during the development of this package. The author maintains full responsibility for the accuracy, reproducibility, and scientific validity of all code. AI-assisted outputs were reviewed and validated against expected behavior before integration. The research questions, analytical approaches, parameter selections, and scientific interpretations were determined independently by the author without AI input.

## License & Citation

MIT License - Created by Natalie Gill

If you use BioGraph Explorer in your research:

```bibtex
@software{biograph_explorer,
  author = {Gill, Natalie},
  title = {BioGraph Explorer: Multi-gene TRAPI Query Integration with Network Analysis},
  year = {2025},
  url = {https://github.com/yourusername/biograph_explorer}
}
```
