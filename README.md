# BioGraph Explorer 🧬

**Status**: Phase 2 Complete ✅ | Phase 3 In Progress 📋

Streamlit application for exploring biomedical knowledge graphs through multi-gene TRAPI queries, network clustering, and interactive visualization.

## Overview

BioGraph Explorer queries the NCATS Translator knowledge graph system to find connections between genes and diseases, then identifies convergent pathways and potential therapeutic targets using network analysis.

**Current Features (Phase 1 and 2):**
- ✅ Gene normalization using NCATS Translator
- ✅ Parallel TRAPI queries across 15+ Translator APIs
- ✅ NetworkX graph construction
- ✅ Louvain community detection & centrality analysis
- ✅ Interactive Cytoscape.js visualization with 9 layout algorithms
- ✅ CSV import and example datasets (Alzheimer's, COVID-19)
- ✅ Query caching for performance

**Planned Features (Phase 3):**
- 📋 Improved query builder (3-hop queries, gene -> [intermediate] -> disease associated phenotypes)
- 📋 Better clustering algorithms suited for each query type
- 📋 LLM cluster summarization
- 📋 RAG-powered chat interface with Claude AI with knowledge provenance subgraph display
- 📋 Session management (save/load analysis sessions)

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


### 2. Try an Example Query

1. In the sidebar, select **"Example Dataset"**
2. Choose **"Alzheimer's Disease (15 genes)"**
3. Select query pattern: **"1-hop: Gene → Any Connection"** (recommended for first run)
4. Click **"Run Query"** (takes ~3-5 minutes)

### 3. Explore Results

The app displays results in three tabs:

- **Overview**: Graph statistics, node categories, top convergent nodes
- **Network**: Interactive Cytoscape.js visualization (drag to pan, scroll to zoom)
- **Communities**: Detected biological modules with hub nodes

### 4. Try Your Own Genes

Create a CSV with a `gene_symbol` column (see [data/test_genes/alzheimers_genes.csv](data/test_genes/alzheimers_genes.csv) for format) or enter genes manually in the sidebar.

## Understanding Your Results

**Query Patterns:**
- **1-hop** (`Gene → Any`): Finds all connections to your genes (broad discovery)
- **2-hop** (`Gene → Intermediate → Disease`): Finds therapeutic targets between genes and a disease (requires disease CURIE like `MONDO:0004975`)

**Key Metrics:**
- **Gene Frequency**: How many query genes connect to each node (convergence indicator)
- **PageRank**: Relative importance in the network
- **Communities**: Biological modules detected by Louvain clustering



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

## Current Limitations

- **Knowledge provenance incomplete**: Edge sources/publications not fully displayed (in progress)
- **API success rate**: Expect 40-60% of TRAPI APIs to succeed (6-10 out of 15) - this is normal
- **Large graphs**: Networks >200 nodes are auto-sampled for visualization performance
- **Session persistence**: Manual query re-run required (auto-save coming in Phase 3)

## Troubleshooting

**"Query returned no results"**
- Normal for some API combinations - at least 5-6 APIs should succeed
- Check gene normalization succeeded (view progress messages)
- Try different genes (some have sparse data in Translator)

**"Graph is empty"**
- Verify disease CURIE format for 2-hop queries (`MONDO:0004975` for Alzheimer's)
- Check that genes are valid HUGO symbols (uppercase)

**"Visualization is slow"**
- Use simpler layout algorithms (dagre, breadthfirst vs. fcose)
- Reduce `max_intermediates` slider to sample fewer nodes
- View individual communities instead of full graph


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
