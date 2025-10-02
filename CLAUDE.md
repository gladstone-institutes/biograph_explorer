# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This project aims to build a multi-gene pathfinder query system using TCT (Translator Component Toolkit) to query the NCATS Translator ecosystem. The system takes differentially expressed genes and a disease/condition, performs multiple Translator queries, clusters the resulting knowledge graphs, and presents them in a digestible format.

## Key Architecture Components

### Data Flow
1. **Input**: Set of DE genes (Ensembl gene ID, Ensembl transcript ID, or HUGO gene IDs) + single disease/condition
2. **Query**: Multiple TRAPI queries via TCT to Translator KGs
3. **Processing**: Collect and cluster resulting knowledge graphs
4. **Output**: Clustered graphs with semantic labels (using small open-source LLM)
5. **Storage**: Neo4j database for graph storage and querying
6. **Interface**: Streamlit RAG interface for "chatting" with networks using small models (e.g., Qwen2.5 14B)

### TRAPI (Translator Reasoner API) Standards
- All queries and responses must conform to TRAPI specification
- Uses Biolink model for semantic entity/relationship descriptions
- Query structure: query_graph → knowledge_graph → results (with bindings)
- Node keys in knowledge graphs are CURIEs (e.g., "MONDO:0005148", "CHEBI:6801")
- Edge predicates use Biolink predicates (e.g., "biolink:treats")

### TCT Integration
- Install: `pip install TCT`
- Key utilities:
  - PathFinder: For multi-gene pathfinding queries
  - ConnectionFinder: Finding connections between entities
  - NetworkFinder: Building networks
- Development install: Clone TCT repo and use `uv sync` or `pip install -e .`

## Test Dataset

Default test case from PMC11255397:
- **Genes**: CD6, IFITM3, IFITM2, STAT5A, KLRG1, DPP4, IL32, PIK3AP1, FYN, IL4R
- **Disease**: COVID-19

## Python Environment

The repository uses a conda environment located in `.conda/`. No requirements.txt or environment.yml file currently exists in the root.

## Development Setup

1. Ensure conda environment is activated (environment located in `.conda/`)
2. Install TCT: `pip install TCT`
3. For development mode TCT: Clone TCT repo and run `uv sync` or `pip install -e .`
4. Future: Neo4j database setup required for graph storage
5. Future: Streamlit app for RAG interface

## Design Principles

- **Minimal Viable System**: Focus on core functionality first
- **TRAPI Compliance**: All interactions must conform to TRAPI standards
- **Modularity**: Separate query, clustering, labeling, storage, and interface components
- **Local-First**: Target small models that can run on laptops (Qwen2.5 14B)

## Reference Documentation

- TRAPI spec: `relevant_documentation/TRAPI_documentation.txt`
- TCT docs: `relevant_documentation/TCT_documentation.txt`
- Online TCT docs: https://ncatstranslator.github.io/Translator_component_toolkit/
