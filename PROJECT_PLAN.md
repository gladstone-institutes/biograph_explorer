# BioGraph Explorer: Project Plan

**Goal**: Streamlit application for multi-gene TRAPI query integration with NetworkX clustering and LLM-assisted exploration using Cytoscape.js visualization

**Current Status**: Phase 2 (Clustering & UI) ✅ COMPLETED | Phase 3 (RAG) 🎯 IN PROGRESS

**Timeline**: 4 weeks MVP
**Team**: 1-2 developers
**Stack**: Python, Streamlit, NetworkX, st-link-analysis (Cytoscape.js), Claude API, TRAPI

**Working Application**: `app.py` (full Streamlit app with interactive visualization)

---

## Architecture

```
biograph_explorer/
├── config/
│   ├── settings.py              # Configuration management
│   └── logging_config.py
├── core/
│   ├── trapi_client.py          # TRAPI query & caching ✅
│   ├── graph_builder.py         # TRAPI → NetworkX conversion ✅
│   ├── clustering_engine.py     # NetworkX analysis (centrality, communities) ✅
│   └── rag_system.py            # Claude + NetworkX RAG (Phase 3)
├── utils/
│   ├── validators.py            # Input validation ✅
│   ├── formatters.py            # Data formatting ✅
│   └── persistence.py           # Pickle NetworkX graphs ✅
├── ui/
│   ├── input_panel.py           # Gene/disease input (integrated in app.py) ✅
│   ├── query_status.py          # Progress tracking (integrated in app.py) ✅
│   ├── results_overview.py      # Summary dashboard (integrated in app.py) ✅
│   ├── convergence_view.py      # Convergent nodes table (integrated in app.py) ✅
│   ├── network_viz.py           # Cytoscape.js rendering via st-link-analysis ✅
│   └── rag_chat.py              # Visual RAG chat interface (Phase 3)
├── data/
│   ├── cache/                   # TRAPI response cache ✅
│   ├── sessions/                # Session state (future: pickled graphs)
│   ├── exports/                 # Visualization exports
│   └── test_genes/              # Example datasets (Alzheimer's, COVID-19) ✅
└── tests/
    ├── test_trapi_client.py     ✅
    ├── test_graph_builder.py    ✅
    ├── test_clustering.py       ✅
    ├── test_validators.py       ✅
    └── fixtures/
        └── alzheimers_test_case.json
```

---

### Test Gene List (15 genes)
| Gene Symbol | Gene ID | Known Role | Expected Finding |
|-------------|---------|------------|------------------|
| APOE | HGNC:613 | Lipid transport | Central hub node |
| APP | HGNC:620 | Amyloid precursor | Converges on BACE1 |
| PSEN1 | HGNC:9508 | γ-secretase component | Converges on Aβ production |
| PSEN2 | HGNC:9509 | γ-secretase component | Converges on Aβ production |
| MAPT | HGNC:6893 | Tau protein | Separate cluster from amyloid |
| TREM2 | HGNC:17761 | Microglial receptor | Inflammation pathway |
| CLU | HGNC:2095 | Clusterin/apoJ | Converges with APOE on lipid |
| CR1 | HGNC:2328 | Complement receptor | Inflammation cluster |
| BIN1 | HGNC:1052 | Endocytosis | Converges on APP trafficking |
| PICALM | HGNC:8301 | Clathrin assembly | Converges on endocytosis |
| CD33 | HGNC:1659 | Myeloid receptor | Inflammation pathway |
| MS4A6A | HGNC:13378 | Membrane protein | Immune cluster |
| ABCA7 | HGNC:29 | ATP-binding cassette | Lipid metabolism with APOE |
| SORL1 | HGNC:11140 | Sorting receptor | APP trafficking |
| BACE1 | HGNC:933 | β-secretase | Central convergence point |

**Disease**: Alzheimer's Disease (MONDO:0004975)

### Expected Outcomes (Validation Criteria)

1. **Convergent nodes**: BACE1, β-amyloid, cholesterol metabolism proteins
2. **Communities**:
   - Cluster 1: Amyloid processing (APP, PSEN1/2, BACE1, BIN1)
   - Cluster 2: Lipid metabolism (APOE, CLU, ABCA7)
   - Cluster 3: Neuroinflammation (TREM2, CD33, CR1, MS4A6A)
3. **Top targets**: BACE1, APOE, TREM2 (all have drugs in development)
4. **Regulatory circuits**: APOE-lipid homeostasis feedback loop

## Core Components

### 1. TRAPI Client (`core/trapi_client.py`) ✅ COMPLETED
- Batch query genes to disease
- Support for 1-hop (neighborhood discovery) and 2-hop (gene → intermediate → disease) patterns
- Intermediate entity type filtering (Protein, ChemicalEntity, Gene, etc.)
- Response caching in `data/cache/`
- Retry logic with exponential backoff
- Progress callbacks for Streamlit UI
- Graceful degradation (accepts 40-60% API success rate)

### 2. Graph Builder (`core/graph_builder.py`) ✅ COMPLETED
- Convert TRAPI responses → NetworkX DiGraph
- Track source gene → node mappings (CURIE to symbol)
- Calculate gene_frequency (convergence metric)
- Rich node attributes: label, category, original_symbol, is_query_gene
- Preserve edge attributes: predicate, knowledge_source, publications
- Export to pickle format

### 3. Clustering Engine (`core/clustering_engine.py`) ✅ COMPLETED
**Algorithms**:
- Centrality: PageRank, betweenness, degree
- Convergent nodes: Filter by gene_frequency ≥ threshold
- Communities: Louvain algorithm (python-louvain)
- Graph statistics: density, average_degree, connected_components

**Output**: ClusteringResults with:
- List of communities with size, density, top nodes
- Modularity score
- Graph statistics
- Ranked convergent nodes

### 4. Cytoscape.js Visualization (`ui/network_viz.py`) ✅ COMPLETED
**Features**:
- st-link-analysis component (Streamlit-native Cytoscape.js)
- Multiple layout algorithms: cose, fcose, circle, grid, breadthfirst, concentric
- Material Icons for node categories (biotech, local_hospital, science, medication, nature, etc.)
- Node sizing by: gene_frequency, pagerank, betweenness, degree
- Color scheme by category:
  - Gene: Red (#E74C3C)
  - Disease: Purple (#9B59B6)
  - Protein: Cyan (#4ECDC4)
  - ChemicalEntity: Orange (#F39C12)
  - BiologicalProcess: Green (#2ECC71)
  - Cluster: Teal (#16A085)
  - Other: Blue (#3498DB)
- Smart graph sampling (guarantees ≥2 edges per query gene)
- Built-in fullscreen and JSON export

**Configuration**:
- <200 nodes: Full rendering
- \>200 nodes: Auto-sample with warning

### 5. RAG System (`core/rag_system.py`) 📋 PLANNED (Phase 3)
**Context Strategy** (3-layer):
- Layer 1: Graph statistics + top convergent nodes (2K tokens)
- Layer 2: Question-relevant subgraph (3-8K tokens)
- Layer 3: Detailed node properties (1-2K tokens per node)

**Citation Mechanism**:
- Claude Haiku 4 with tool use (structured output)
- Tool: `cite_graph_evidence` returns node_ids, metric_name, metric_value
- Validation: Check citations against actual graph
- Visualization: Extract subgraph → render with Cytoscape.js

---

## Test Case: Alzheimer's Disease (15 genes)

**Genes**: APOE, APP, PSEN1, PSEN2, MAPT, TREM2, CLU, CR1, BIN1, PICALM, CD33, MS4A6A, ABCA7, SORL1, BACE1
**Disease**: MONDO:0004975

**Expected Results**:
- Convergent nodes: BACE1, amyloid-β, cholesterol pathway proteins
- Communities: 3-4 clusters (amyloid processing, lipid metabolism, neuroinflammation)
- Top targets: BACE1, APOE, TREM2

**Validation**: Results match known biology (80%+ accuracy)

---

## Data Persistence Strategy

### Current Implementation

| Stage | Location | Status | Purpose |
|-------|----------|--------|---------|
| TRAPI Responses | `data/cache/` | ✅ | Individual API responses cached by query hash |
| Session State | Streamlit session_state | ✅ | In-memory graph, clustering results, query genes |
| Exports | `data/exports/` | ✅ | Standalone HTML visualizations |

### Planned: Load Previous Query Feature 📋 NEW

**Goal**: Allow users to reload previous query results from cached data without re-querying APIs

**Implementation**:
1. **Session Folder Structure**:
   ```
   data/sessions/{session_id}/
   ├── input_config.json          # Gene list, disease CURIE, query pattern
   ├── trapi_response.pkl         # Full TRAPIResponse object
   ├── graph.pkl                  # NetworkX graph with all attributes
   ├── clustering_results.pkl     # ClusteringResults object
   └── metadata.json              # Timestamp, query stats, API success rate
   ```

2. **UI Integration**:
   - Sidebar: "Load Previous Query" expander
   - List available sessions from `data/sessions/` with metadata (date, gene count, disease)
   - Select session → populate session state → navigate to results view
   - Display "⚠️ Loaded from cache" banner

3. **Save Session**:
   - Button in results view: "💾 Save Session"
   - Generate session_id from timestamp + gene list hash
   - Serialize all session state to folder
   - Show success message with session_id

4. **Auto-save**:
   - Optional: Auto-save after successful query completion
   - Configurable in settings (default: enabled)

**Benefits**:
- No re-querying when exploring different visualizations
- Share sessions via folder export
- Track analysis history
- Resume interrupted sessions

---

## Development Phases

### Phase 1: Foundation & POC ✅ COMPLETED
**Status**: All functionality working in `notebooks/multi_gene_pathfinder.ipynb`

**Accomplishments**:
- ✅ Gene normalization using TCT's `name_resolver.batch_lookup()`
- ✅ TRAPI integration with neighborhood discovery pattern (empty target list)
- ✅ Parallel API querying (15 Translator APIs, 6/15 success rate = normal)
- ✅ NetworkX graph construction from TRAPI responses
- ✅ Rich node attributes (labels, categories, is_query_gene flags)
- ✅ Interactive ipycytoscape visualization with 3-tier color coding
- ✅ Smart edge sampling guaranteeing all query genes appear
- ✅ JSON persistence of raw TRAPI results in `data/raw/`
- ✅ Robust error handling and graceful API degradation

**Test Results**:
- Successfully normalized 10/10 COVID-19 DE genes
- Retrieved 1,537 edges (812 unique after dedup) from knowledge graphs
- Built graph with 476 nodes, 723 edges, single connected component

---

### Phase 2: Graph Processing, Clustering & UI ✅ COMPLETED
**Goal**: Convert notebook prototype to modular package with clustering and Streamlit UI

**Accomplishments**:
1. **Modular Package Structure** ✅
   - ✅ `core/trapi_client.py` - TRAPI querying with caching
   - ✅ `core/graph_builder.py` - NetworkX graph construction
   - ✅ `core/clustering_engine.py` - Louvain community detection + centrality
   - ✅ `utils/persistence.py` - Session persistence utilities
   - ✅ `utils/validators.py` - Input validation for genes and disease CURIEs
   - ✅ `utils/formatters.py` - Data formatting utilities
   - ✅ `ui/network_viz.py` - Cytoscape.js visualization

2. **Clustering Implementation** ✅
   - ✅ Louvain community detection
   - ✅ Centrality metrics: PageRank, betweenness, degree
   - ✅ Hub node identification within clusters
   - ✅ Graph statistics: density, modularity, connected components

3. **Streamlit UI** ✅
   - ✅ Input panel with three methods: Example datasets, CSV upload, Manual entry
   - ✅ Query pattern selector: 1-hop vs 2-hop with intermediate type filtering
   - ✅ Progress tracking during query execution
   - ✅ Three-tab results view:
     - Overview: Graph metrics, statistics, node categories
     - Network: Interactive Cytoscape.js visualization with controls
     - Communities: Community detection results with top nodes
   - ✅ Node inspection panel with detailed metrics and edges
   - ✅ Layout algorithm selector (6 options)
   - ✅ Node sizing metric selector
   - ✅ Graph sampling slider
   - ✅ Cluster view toggle (meta-graph)

4. **Visualization Migration** ✅
   - ✅ Migrated from PyVis to st-link-analysis (Cytoscape.js)
   - ✅ Material Icons for biological entity types
   - ✅ Multiple layout algorithms
   - ✅ Better Streamlit integration (no iframe)
   - ✅ Built-in fullscreen and JSON export

**Outcomes**:
- ✅ Modular codebase following architecture plan
- ✅ 3-4 detected communities for Alzheimer's test case
- ✅ Full-featured Streamlit app for running queries
- ✅ Example datasets (Alzheimer's, COVID-19) included
- ✅ Test suite with >80% coverage

---

### Phase 3: RAG + Session Management 🎯 IN PROGRESS
**Goal**: LLM-assisted exploration with Claude API integration and session persistence

**Tasks**:
1. **Session Management** 📋 NEW
   - [ ] Implement session folder structure (`data/sessions/{id}/`)
   - [ ] Create session save/load functions in `utils/persistence.py`
   - [ ] Add "Load Previous Query" UI in sidebar
   - [ ] Add "Save Session" button in results view
   - [ ] Display session metadata (date, genes, disease, stats)
   - [ ] Optional auto-save after query completion

2. **RAG System**
   - [ ] Implement 3-layer context strategy
   - [ ] Claude Haiku 4 integration with tool use
   - [ ] Citation validation logic
   - [ ] Subgraph extraction for citations

3. **Interactive Chat UI**
   - [ ] Streamlit chat interface in new tab
   - [ ] Click citation → highlight nodes in visualization
   - [ ] Question answering about convergent nodes
   - [ ] Explanation of community structure
   - [ ] Drug target recommendations

4. **Multiple Predicates Handling** 📋 NEW
   - [ ] Switch from `DiGraph` to `MultiDiGraph` in `graph_builder.py`
   - [ ] Update edge processing in `network_viz.py` to aggregate predicates for visualization
   - [ ] Display multiple predicates as comma-separated or count (e.g., "3 relationships")
   - [ ] Store full predicate list in edge tooltip data
   - [ ] Update edge styles to handle aggregated predicate labels
   - [ ] Add style for "multiple_predicates" edges

**Expected Outcomes**:
- [ ] Users can save and reload query sessions
- [ ] RAG system answers questions about graph
- [ ] Citations link to visual evidence in Cytoscape.js
- [ ] Chat interface integrated in results view
- [ ] All TRAPI relationships preserved without edge overwrites

---

### Phase 4: UI Polish & Production 📋 PLANNED
**Goal**: Production-ready application with deployment

**Tasks**:
1. **Enhanced UI**
   - [ ] Convergent nodes table with sorting/filtering
   - [ ] Session history viewer
   - [ ] Batch session export/import
   - [ ] Export to PNG, PDF reports

2. **Testing & Documentation**
   - [ ] Expand test suite to >90% coverage
   - [ ] User guide with screenshots
   - [ ] API documentation (ReadTheDocs)
   - [ ] Video tutorials



---

## Dependencies

**Current** (`pyproject.toml`):
```toml
python = ">=3.11,<4.0"
tct = "^0.1.4"
streamlit = "^1.50.0"
st-link-analysis = "^0.4.0"  # ✅ NEW (replaced pyvis)
pydantic = "^2.11.9"
networkx = "^3.5"
pandas = "^2.2.1"
python-louvain = "^0.16"
```

**Planned** (Phase 3):
```toml
anthropic = "^0.40.0"  # Claude API
```

---

## Cytoscape.js Configuration

**Layout Algorithms** (st-link-analysis):
- **cose**: Physics-based, good for general graphs (default)
- **fcose**: Fast force-directed, better for large graphs
- **circle**: Circular layout
- **grid**: Grid layout
- **breadthfirst**: Hierarchical tree layout
- **concentric**: Concentric circles by centrality

**Node Styling**:
- Size: `15px + (normalized_metric * 30px)` (range: 15-45px)
- Query genes: minimum 25px
- Color: Category-based (see section 4 above)
- Icon: Material Icons by category

**Edge Styling**:
- Color: Gray (#808080)
- Caption: Predicate label (cleaned biolink term)
- Directed: Yes (arrows shown)
- Curve: Bezier

**Interaction**:
- Drag to pan
- Scroll to zoom
- Click node to select (future: show in details panel)
- Fullscreen button
- JSON export button

---

## Performance Targets

| Operation | Target | Notes |
|-----------|--------|-------|
| 15-gene TRAPI batch | <5 min | With caching/parallelization |
| Clustering analysis | <30 sec | NetworkX in-memory |
| Cytoscape.js render (<100 nodes) | <2 sec | Fast layout algorithms |
| Cytoscape.js render (100-200 nodes) | <5 sec | cose/fcose layouts |
| LLM response | <10 sec | Claude Haiku (Phase 3) |
| Page load | <2 sec | Any view |
| Session save | <5 sec | Pickle serialization |
| Session load | <3 sec | Pickle deserialization |

---

## Model Selection: Claude Haiku 4 (Phase 3)

**Why Haiku**:
- Task is information extraction, not complex reasoning
- Graph analysis already done by NetworkX
- Tool use for structured citations
- Cost: $0.02 per session (5 questions × 12K input + 1K output)

---

## Success Metrics

**Phase 1 (Completed)** ✅:
- ✅ Process 10 genes in <5 minutes (achieved: 1,537 raw edges, ~2 min)
- ✅ Build NetworkX graph from TRAPI responses (476 nodes, 723 edges)
- ✅ Render interactive visualization (ipycytoscape, 50-edge sampling)
- ✅ Gene normalization with 100% success rate (10/10 COVID genes)
- ✅ Graceful API degradation (40% success rate = sufficient)

**Phase 2 (Completed)** ✅:
- ✅ Detect 3-4 communities in Alzheimer's test case
- ✅ Identify convergent nodes (gene_frequency metric)
- ✅ Modular package structure following architecture plan
- ✅ Full-featured Streamlit UI with three-tab results view
- ✅ Interactive Cytoscape.js visualization with 6 layout algorithms
- ✅ Test suite with >80% coverage

**Phase 3 (In Progress)** 🎯:
- [ ] Session save/load functionality working
- [ ] Users can browse and reload previous queries
- [ ] RAG system with Claude Haiku 4
- [ ] Citation accuracy >90% (metrics match graph)
- [ ] Chat interface integrated in results view

**Phase 4 (Planned)** 📋:
- [ ] Non-expert can run analysis in <10 clicks
- [ ] Complete documentation (README, user guide, API docs)
- [ ] Streamlit Cloud deployment
- [ ] PyPI package release

---

## Risks & Mitigations

| Risk | Mitigation |
|------|------------|
| TRAPI downtime | ✅ Cache all responses, retry logic |
| Large graphs (>200 nodes) | ✅ Auto-sample with warning, cluster view |
| Citation parsing errors | Use structured output (tool use) - Phase 3 |
| Browser compatibility | Test Chrome/Firefox/Safari |
| Session file corruption | Validate on load, keep backups |
| Large session files | Compress pickles, set retention policy |

---

## Future Enhancements (Post-MVP)

- Multi-disease comparison view
- Augmenting TRAPI graphs with additional data sources:
  - WikiPathways
  - PFOCR
  - Human Protein Atlas
  - STRING protein interactions
- Differential analysis (compare gene sets)
- Export to Cytoscape desktop format
- API endpoint for programmatic access
- Collaborative sessions (multi-user)
- Integration with notebook environments (Jupyter widget)

---

## Deliverables

**Phase 2 (Completed)** ✅:
1. ✅ Working Streamlit app (`app.py`)
2. ✅ Modular package structure
3. ✅ Interactive Cytoscape.js visualization
4. ✅ Example datasets (Alzheimer's, COVID-19)
5. ✅ Test suite with >80% coverage
6. ✅ Documentation (README, CLAUDE.md)

**Phase 3 (In Progress)** 🎯:
1. Session management system
2. RAG chat interface
3. Citation validation
4. Enhanced user guide

**Phase 4 (Planned)** 📋:
1. Streamlit Cloud deployment
2. ReadTheDocs documentation
3. PyPI package
4. Video tutorials
5. Example session exports
