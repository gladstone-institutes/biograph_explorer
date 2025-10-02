# Project Roadmap: Multi-Gene Pathfinder

## Phase 1: Foundation & POC ✅ COMPLETED

### Accomplishments
- **Gene Normalization**: Successfully normalized 10/10 COVID-19 DE genes using TCT's `name_resolver.batch_lookup()`
- **TCT Integration**: Implemented neighborhood discovery query pattern (empty target list strategy)
- **API Querying**: Queried 15 Translator APIs with graceful degradation (6/15 succeeded = 40% success rate)
- **Data Collection**: Retrieved 1,537 edges from knowledge graphs, saved to JSON
- **Interactive Visualization**: ipycytoscape graph with 3-tier color coding (DE genes, disease, other nodes)
- **Robust Error Handling**: Fallback for API timeouts, informative error messages
- **Connectivity Preservation**: Per-gene sampling guarantees all query genes appear in visualization

### Key Technical Decisions
1. **Single Notebook Architecture**: All functionality contained in `multi_gene_pathfinder.ipynb`
2. **Neighborhood Discovery**: Empty target list finds all connections rather than forcing direct gene→disease paths
3. **Partial Success Acceptable**: 6/15 APIs sufficient for POC validation
4. **Smart Sampling**: Guarantee ≥2 edges per query gene to prevent missing nodes
5. **IOPub Optimization**: Reduced to 50 edges with circle layout to prevent Jupyter crashes

### Current Data Files
- `data/raw/tct_results_*.json`: Full query results (1,537 edges, 2.6MB)
- `multi_gene_pathfinder.ipynb`: Complete working notebook

---

## Phase 2: Graph Processing & Clustering 🎯 NEXT

### Objectives
Convert raw edge data into structured knowledge graphs and identify functional clusters.

### Tasks
1. **Knowledge Graph Construction**
   - Build NetworkX graph from edge DataFrame
   - Add node attributes (labels, categories, types)
   - Add edge attributes (predicates, qualifiers, publications, knowledge_source)
   - Preserve provenance metadata

2. **Graph Clustering**
   - Implement community detection algorithms:
     - Louvain method (modularity optimization)
     - Label propagation (fast, scalable)
     - Spectral clustering (quality clusters)
   - Compare clustering results across algorithms
   - Evaluate clustering quality (modularity, silhouette scores)

3. **Cluster Analysis**
   - Extract subgraphs for each cluster
   - Compute cluster statistics (size, density, centrality measures)
   - Identify hub nodes within clusters
   - Export cluster assignments to DataFrame

4. **Visualization Enhancements**
   - Color nodes by cluster assignment
   - Add cluster boundaries/labels
   - Implement cluster-level graph view (meta-graph)

### Expected Outputs
- `edges_clustered.csv`: Edges with cluster assignments
- `cluster_statistics.json`: Cluster metrics and characteristics
- Interactive visualization showing cluster organization

---

## Phase 3: Semantic Labeling 🔮 FUTURE

### Objectives
Automatically generate human-readable labels for clusters using lightweight LLM.

### Tasks
1. **Cluster Context Extraction**
   - For each cluster, collect:
     - Node labels and categories
     - Predicate types
     - Qualifier patterns
     - Sample publications/evidence
   - Format as structured prompt for LLM

2. **LLM Integration**
   - Test small models: Qwen2.5 3B/7B (runs on laptop)
   - Design prompt template: "Summarize this biological pathway cluster in 3-5 words"
   - Implement batch processing for all clusters
   - Add fallback: rule-based labels if LLM unavailable

3. **Label Refinement**
   - Validate labels don't overlap
   - Ensure labels reflect dominant biological themes
   - Allow manual override in notebook interface

### Expected Outputs
- `cluster_labels.json`: Cluster ID → semantic label mapping
- Labeled visualization with cluster names

---

## Phase 4: Neo4j Integration 📊 FUTURE

### Objectives
Load knowledge graphs into Neo4j database for persistent querying and augmentation.

### Tasks
1. **Schema Design**
   - Map Biolink categories to Neo4j node labels
   - Map predicates to relationship types
   - Design property model for qualifiers and provenance
   - Plan indexes (CURIE, label, cluster_id)

2. **Data Import**
   - Convert NetworkX → Neo4j using py2neo or neo4j-python-driver
   - Batch import for performance (use UNWIND)
   - Add cluster assignments as node properties
   - Import semantic labels

3. **Query Interface**
   - Implement Cypher query templates:
     - "Find paths between gene X and disease Y"
     - "Show cluster containing gene X"
     - "What genes connect to X via predicate P?"
   - Create reusable query functions in notebook

4. **Optional: External Data Augmentation**
   - Integrate Human Protein Atlas (tissue expression)
   - Add protein-protein interactions (STRING)
   - Link to clinical trials or drug databases

### Expected Outputs
- Neo4j database with full knowledge graph
- Cypher query library
- Updated notebook with Neo4j query examples

---

## Phase 5: Streamlit RAG Interface 💬 FUTURE

### Objectives
Build conversational interface to "chat" with knowledge graphs using RAG pattern.

### Tasks
1. **RAG Pipeline Design**
   - Vector embeddings for nodes/clusters (sentence-transformers)
   - Vector store: FAISS or ChromaDB
   - Retrieval: find relevant graph context for user questions
   - Generation: use Qwen2.5 14B to answer with context

2. **Streamlit App Development**
   - Chat interface with message history
   - Query suggestions (e.g., "How does FYN relate to COVID-19?")
   - Visualize retrieved graph context inline
   - Export conversation to markdown

3. **LLM Integration**
   - Local inference: llama-cpp-python or vLLM
   - Prompt engineering: format graph context effectively
   - Response grounding: cite specific edges/publications
   - Handle follow-up questions with conversation memory

4. **Deployment**
   - Docker container with Neo4j + Streamlit + model
   - Resource optimization for laptop deployment
   - Documentation for setup and usage

### Expected Outputs
- `streamlit_app/` directory with complete application
- Docker Compose configuration
- User guide for deployment

---

## Phase 6: Polish & Documentation 📝 FINAL

### Tasks
1. **Code Cleanup**
   - Refactor notebook into reusable functions
   - Add comprehensive docstrings
   - Type hints for all functions
   - Unit tests for core utilities

2. **Performance Optimization**
   - Profile slow operations (clustering, visualization)
   - Optimize graph algorithms
   - Cache expensive computations
   - Increase visualization edge limit if IOPub allows

3. **Documentation**
   - Complete README with setup instructions
   - Tutorial notebook with walkthrough
   - API reference for key functions
   - Architecture diagrams

4. **Validation**
   - Test with additional gene sets
   - Test with different diseases
   - Compare results with literature
   - Gather user feedback

---

## Technical Debt & Improvements

### Known Limitations
1. **Visualization Constraints**: Limited to 50 edges due to IOPub rate limit
   - Consider alternative: Plotly Dash (server-side rendering)
   - Or increase Jupyter `iopub_msg_rate_limit` config

2. **API Coverage**: Only 6/15 APIs succeeded
   - Investigate failures for remaining 9 APIs
   - Add retry logic with exponential backoff
   - Consider alternative APIs (e.g., SPOKE, RTX-KG2)

3. **Missing Connections**: 2/10 genes had no edges in results
   - May need multi-hop queries (2-hop or 3-hop paths)
   - Consider expanding predicates beyond current set

4. **Scalability**: Current approach tested with 10 genes
   - Test with 50-100 gene lists (RNA-seq typical output)
   - May need batching strategy or distributed queries

### Potential Enhancements
1. **Query Strategies**: Add support for gene→gene queries (find common pathways)
2. **Edge Filtering**: Filter by evidence strength, publication count, KP trustworthiness
3. **Differential Analysis**: Compare graphs for different diseases or conditions
4. **Export Options**: GraphML, Cytoscape JSON, CSV for external tools
5. **Provenance Tracking**: Detailed logs of which APIs contributed which edges

---

## Resources & References

### Documentation
- TRAPI Specification: https://github.com/NCATSTranslator/ReasonerAPI
- Biolink Model: https://biolink.github.io/biolink-model/
- TCT Repository: https://github.com/gloriachin/Translator_component_toolkit

### Key Notebooks
- PathFinder Example: `notebooks/Path_finder.ipynb` (reference implementation)
- ConnectionFinder Example: `notebooks/Connection_finder.ipynb`

### Data Sources
- Test Gene List: PMC11255397 (COVID-19 differentially expressed genes)
- Successful KP APIs: Automat-ubergraph, Automat-robokop, Automat-pharos, Automat-hetionet, Service Provider TRAPI, BioThings Explorer

### Tools & Libraries
- TCT: Translator Component Toolkit
- NetworkX: Graph algorithms
- ipycytoscape: Interactive graph visualization
- Neo4j: Graph database
- Streamlit: Web app framework
- Qwen2.5: Lightweight LLM for semantic labeling

---

## Getting Started with Phase 2

To continue development, open `multi_gene_pathfinder.ipynb` and add a new section:

```python
# =============================================================================
# PHASE 2: GRAPH CLUSTERING
# =============================================================================

# 1. Build full NetworkX graph from edges_df
G = nx.from_pandas_edgelist(
    edges_df,
    source='SubjectCURIE',
    target='ObjectCURIE',
    edge_attr=['Predicate', 'PublicationCount'],
    create_using=nx.DiGraph()
)

# 2. Add node attributes (labels, categories)
for _, row in nodes_df.iterrows():
    G.nodes[row['CURIE']]['label'] = row['Label']
    G.nodes[row['CURIE']]['category'] = row['Category']

# 3. Apply community detection
import community  # python-louvain
communities = community.best_partition(G.to_undirected())

# 4. Analyze clusters
for cluster_id in set(communities.values()):
    cluster_nodes = [n for n, c in communities.items() if c == cluster_id]
    print(f"Cluster {cluster_id}: {len(cluster_nodes)} nodes")
```

Continue iteratively building on the foundation established in Phase 1.
