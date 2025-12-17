"""BioGraph Explorer - Streamlit Application

Main Streamlit app for multi-gene TRAPI query integration with NetworkX clustering.
"""

import streamlit as st
from pathlib import Path
import pandas as pd
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)

from biograph_explorer.core import TRAPIClient, GraphBuilder, ClusteringEngine
from biograph_explorer.utils import validate_gene_list, validate_disease_curie, ValidationError
from biograph_explorer.utils.biolink_predicates import (
    GRANULARITY_PRESETS,
    get_allowed_predicates_for_display,
)

# Page config
st.set_page_config(
    page_title="BioGraph Explorer",
    page_icon=":material/biotech:",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Initialize session state
if 'graph' not in st.session_state:
    st.session_state.graph = None
if 'knowledge_graph' not in st.session_state:
    st.session_state.knowledge_graph = None  # Store KnowledgeGraph object
if 'clustering_results' not in st.session_state:
    st.session_state.clustering_results = None
if 'query_genes' not in st.session_state:
    st.session_state.query_genes = []
if 'response' not in st.session_state:
    st.session_state.response = None
if 'selected_node' not in st.session_state:
    st.session_state.selected_node = None
if 'disease_curie' not in st.session_state:
    st.session_state.disease_curie = None

# Title
st.title(":material/biotech: BioGraph Explorer")
st.markdown("Multi-gene TRAPI query integration with NetworkX clustering and visualization")

# Sidebar: Input Configuration
st.sidebar.header(":material/input: Input Configuration")

# Input method selection
input_method = st.sidebar.radio(
    "Input Method",
    ["Example Dataset", "Upload CSV", "Manual Entry", "Load Cached Query"],
    index=0
)

genes = []
disease_curie = "MONDO:0100096"  # Default: COVID-19

if input_method == "Example Dataset":
    dataset_choice = st.sidebar.selectbox(
        "Select Example",
        ["COVID-19 (10 genes)", "Alzheimer's Disease (15 genes)"]
    )
    
    if dataset_choice == "Alzheimer's Disease (15 genes)":
        csv_path = "data/test_genes/alzheimers_genes.csv"
        disease_curie = "MONDO:0004975"
    else:
        csv_path = "data/test_genes/covid19_genes.csv"
        disease_curie = "MONDO:0100096"

    # Update session state when dataset changes
    st.session_state.disease_selected_curie = disease_curie

    try:
        df = pd.read_csv(csv_path)
        genes = df['gene_symbol'].tolist()
        st.sidebar.success(f":material/check_circle: Loaded {len(genes)} genes from example dataset")
        with st.sidebar.expander("View genes"):
            st.dataframe(df, width='stretch')
    except Exception as e:
        st.sidebar.error(f"Error loading example: {e}")

elif input_method == "Upload CSV":
    uploaded_file = st.sidebar.file_uploader(
        "Upload Gene List CSV",
        type=["csv"],
        help="CSV file with 'gene_symbol' column"
    )
    
    if uploaded_file:
        try:
            df = pd.read_csv(uploaded_file)
            if 'gene_symbol' in df.columns:
                genes = df['gene_symbol'].tolist()
                st.sidebar.success(f":material/check_circle: Loaded {len(genes)} genes from CSV")
                with st.sidebar.expander("View genes"):
                    st.dataframe(df, width='stretch')
            else:
                st.sidebar.error("CSV must have 'gene_symbol' column")
        except Exception as e:
            st.sidebar.error(f"Error reading CSV: {e}")

elif input_method == "Manual Entry":
    gene_text = st.sidebar.text_area(
        "Gene Symbols",
        height=150,
        help="Enter gene symbols, one per line or comma-separated",
        placeholder="APOE\nAPP\nPSEN1\n..."
    )

    if gene_text:
        # Parse input
        genes = [g.strip().upper() for g in gene_text.replace(',', '\n').split('\n') if g.strip()]
        st.sidebar.info(f":material/edit_note: {len(genes)} genes entered")

else:  # Load Cached Query
    # List available cached files using TRAPIClient
    cache_dir = Path("data/cache")
    client = TRAPIClient(cache_dir=cache_dir)
    cached_files = client.list_cached_queries()

    if cached_files:
        selected_cache = st.sidebar.selectbox(
            "Select Cached Query",
            options=range(len(cached_files)),
            format_func=lambda i: cached_files[i]["label"],
            help="Load a previously cached TRAPI query result"
        )

        if st.sidebar.button(":material/folder_open: Load Cache", type="primary"):
            selected_file = cached_files[selected_cache]["path"]

            try:
                cached_response = client._load_cached_response(selected_file)

                if cached_response:
                    st.sidebar.success(f":material/check_circle: Loaded {len(cached_response.edges)} edges from cache")

                    # Build graph from cached response
                    progress_bar = st.progress(0)
                    status_text = st.empty()

                    status_text.text("Building graph from cached data...")
                    progress_bar.progress(30)

                    builder = GraphBuilder()
                    engine = ClusteringEngine()

                    curie_to_symbol = cached_response.metadata.get("curie_to_symbol", {})
                    curie_to_name = cached_response.metadata.get("curie_to_name", {})

                    # Get disease-associated BP CURIEs for triangle rendering (if this was a Disease BioProcesses query)
                    disease_bp_curies = cached_response.metadata.get("bioprocess_endpoints", None)

                    kg = builder.build_from_trapi_edges(
                        cached_response.edges,
                        cached_response.input_genes,
                        curie_to_symbol,
                        curie_to_name,
                        disease_bp_curies=disease_bp_curies,
                    )

                    progress_bar.progress(60)

                    # Calculate gene frequency
                    gene_freq = builder.calculate_gene_frequency(kg.graph, cached_response.input_genes)
                    for node, freq in gene_freq.items():
                        kg.graph.nodes[node]['gene_frequency'] = freq

                    progress_bar.progress(70)

                    # Cluster
                    status_text.text("Detecting communities...")
                    results = engine.analyze_graph(kg.graph, cached_response.input_genes)

                    progress_bar.progress(90)

                    # Save to session state
                    st.session_state.graph = kg.graph
                    st.session_state.knowledge_graph = kg  # Store KnowledgeGraph object
                    st.session_state.clustering_results = results
                    st.session_state.query_genes = cached_response.input_genes
                    st.session_state.response = cached_response
                    st.session_state.disease_curie = cached_response.target_disease

                    progress_bar.progress(100)
                    status_text.text(":material/check_circle: Cache loaded successfully!")

                    st.rerun()

                else:
                    st.sidebar.error("Failed to load cached file")

            except Exception as e:
                st.sidebar.error(f"Error loading cache: {e}")
                import traceback
                st.sidebar.code(traceback.format_exc())
    else:
        st.sidebar.warning("No cached queries found. Run a query first to create cached results.")

# Disease CURIE lookup
st.sidebar.subheader(":material/search: Disease Lookup")

# Initialize session state for disease lookup
if 'disease_search_results' not in st.session_state:
    st.session_state.disease_search_results = []
if 'disease_selected_curie' not in st.session_state:
    st.session_state.disease_selected_curie = disease_curie

# Input mode toggle
disease_input_mode = st.sidebar.radio(
    "Disease Input Mode",
    ["Search by Name", "Enter CURIE Directly"],
    index=0,
    horizontal=True,
    label_visibility="collapsed"
)

if disease_input_mode == "Search by Name":
    disease_search_query = st.sidebar.text_input(
        "Search Disease",
        placeholder="e.g., Alzheimer, COVID-19, diabetes",
        help="Type a disease name to search. Results will appear below."
    )

    if disease_search_query and len(disease_search_query) >= 2:
        # Perform lookup
        client = TRAPIClient(cache_dir=Path("data/cache"))
        with st.sidebar.status("Searching...", expanded=False):
            search_results = client.lookup_disease(disease_search_query, max_results=10)
        st.session_state.disease_search_results = search_results

        if search_results:
            # Build options for selectbox
            options = {
                f"{r['label']} ({r['curie']})": r
                for r in search_results
            }

            selected_option = st.sidebar.selectbox(
                "Select Disease",
                options=list(options.keys()),
                help="Choose from the search results"
            )

            if selected_option:
                selected_disease = options[selected_option]
                disease_curie = selected_disease['curie']
                st.session_state.disease_selected_curie = disease_curie

                # Show confirmation details
                with st.sidebar.expander(":material/info: Disease Details", expanded=False):
                    st.markdown(f"**CURIE:** `{selected_disease['curie']}`")
                    st.markdown(f"**Label:** {selected_disease['label']}")
                    if selected_disease.get('types'):
                        types_str = ", ".join([t.replace("biolink:", "") for t in selected_disease['types'][:3]])
                        st.markdown(f"**Types:** {types_str}")
                    if selected_disease.get('synonyms'):
                        synonyms_preview = selected_disease['synonyms'][:5]
                        st.markdown(f"**Synonyms:** {', '.join(synonyms_preview)}")
                        if len(selected_disease['synonyms']) > 5:
                            st.caption(f"...and {len(selected_disease['synonyms']) - 5} more")
        else:
            st.sidebar.warning(f"No diseases found for '{disease_search_query}'")
            disease_curie = st.session_state.disease_selected_curie
    else:
        disease_curie = st.session_state.disease_selected_curie
        if disease_curie:
            st.sidebar.info(f":material/check_circle: Using: `{disease_curie}`")

else:
    # Direct CURIE input mode
    disease_curie = st.sidebar.text_input(
        "Disease CURIE",
        value=st.session_state.disease_selected_curie or disease_curie,
        help="Disease CURIE ID (e.g., MONDO:0004975 for Alzheimer's)"
    )
    st.session_state.disease_selected_curie = disease_curie

# Intermediate entity type selector
st.sidebar.markdown("---")
st.sidebar.subheader(":material/tune: Query Configuration")

# Query pattern selector
query_pattern = st.sidebar.radio(
    "Query Pattern",
    options=[
        "1-hop (Neighborhood Discovery)",
        "2-hop (Gene → Intermediate → Disease)",
        "2-hop (Gene → Intermediate → Disease BioProcesses)",
    ],
    index=1,  # Default to 2-hop Disease
    help="1-hop: direct connections. 2-hop to Disease: targets between genes and disease. 2-hop to BioProcesses: targets between genes and disease-associated biological processes."
)

# Intermediate type options for the query patterns
INTERMEDIATE_OPTIONS = [
    "Protein",
    "ChemicalEntity (Drugs/Metabolites)",
    "Gene",
    "PhenotypicFeature",
    "Pathway",
    "BiologicalProcess",
    "AnatomicalEntity",
    "MolecularActivity",
    "CellularComponent",
]

if query_pattern == "2-hop (Gene → Intermediate → Disease)":
    st.sidebar.markdown("**Path:** Gene → **[Intermediate]** → Disease")

    intermediate_types = st.sidebar.multiselect(
        "Intermediate Entity Types",
        options=INTERMEDIATE_OPTIONS,
        default=["Protein", "ChemicalEntity (Drugs/Metabolites)", "Gene"],
        help="Select intermediate entity types for 2-hop query. These are potential therapeutic targets."
    )

    if not disease_curie:
        st.sidebar.warning(":material/warning: Disease CURIE required for 2-hop queries")

elif query_pattern == "2-hop (Gene → Intermediate → Disease BioProcesses)":
    st.sidebar.markdown("**Path:** Gene → **[Intermediate]** → Disease-associated BiologicalProcess")

    intermediate_types = st.sidebar.multiselect(
        "Intermediate Entity Types",
        options=INTERMEDIATE_OPTIONS,
        default=["Protein", "ChemicalEntity (Drugs/Metabolites)", "Gene", "Pathway", "MolecularActivity"],
        help="Select intermediate entity types for 2-hop query. These connect genes to disease-associated biological processes."
    )

    if not disease_curie:
        st.sidebar.warning(":material/warning: Disease CURIE required for BiologicalProcess queries")

    # BiologicalProcess Discovery settings
    with st.sidebar.expander("BiologicalProcess Discovery Settings", expanded=False):
        filter_bp_predicates = st.checkbox(
            "Filter to informative predicates",
            value=True,
            help="Exclude text mining noise like 'occurs_together_in_literature_with' (removes ~85% of noisy edges)"
        )

        use_cached_bp = st.checkbox(
            "Use cached BiologicalProcesses if available",
            value=True,
            help="Skip Stage 1 discovery if we have cached BiologicalProcesses for this disease"
        )

else:
    st.sidebar.markdown("**Path:** Gene → **[Any Connection]**")
    intermediate_types = []  # Empty for 1-hop
    filter_bp_predicates = True  # Default
    use_cached_bp = True  # Default

# Ensure variables are defined for non-BP query patterns
if query_pattern != "2-hop (Gene → Intermediate → Disease BioProcesses)":
    filter_bp_predicates = True
    use_cached_bp = True

# Predicate Granularity Section
st.sidebar.markdown("### Predicate Filtering")

predicate_preset = st.sidebar.selectbox(
    "Relationship Granularity",
    options=["All Relationships", "Standard", "Specific Only"],
    index=2,  # Default to "Specific Only"
    help="Filter out vague relationship types. Standard excludes generic predicates like 'related_to'."
)

exclude_literature = st.sidebar.checkbox(
    "Exclude literature co-occurrence",
    value=True,
    help="Exclude 'occurs_together_in_literature_with' predicate"
)

exclude_coexpression = st.sidebar.checkbox(
    "Exclude coexpression",
    value=False,
    help="Exclude 'coexpressed_with' predicate"
)

exclude_homology = st.sidebar.checkbox(
    "Exclude homology",
    value=True,
    help="Exclude 'homologous_to', 'orthologous_to', 'paralogous_to', 'xenologous_to' predicates"
)

# Show included predicates in expander
with st.sidebar.expander("View included predicates"):
    min_depth = GRANULARITY_PRESETS[predicate_preset]["min_depth"]
    allowed = get_allowed_predicates_for_display(min_depth, exclude_literature, exclude_coexpression, exclude_homology)
    st.write(f"**{len(allowed)}** predicates included:")
    # Display all predicates in a scrollable text area
    predicate_display = "\n".join(sorted(allowed))
    st.text_area(
        "Predicates",
        value=predicate_display,
        height=200,
        disabled=True,
        label_visibility="collapsed"
    )

st.sidebar.markdown("---")

# Query button
run_query = st.sidebar.button(
    ":material/rocket_launch: Run Query",
    type="primary",
    disabled=len(genes) == 0
)

if len(genes) > 0:
    st.sidebar.markdown(f"**Ready to query:** {len(genes)} genes")

# Main area
if not run_query and not st.session_state.graph:
    # Welcome screen
    st.info(":material/arrow_back: Configure your query in the sidebar and click **Run Query** to start")

    col1, col2 = st.columns(2)
    with col1:
        st.subheader(":material/info: What it does")
        st.markdown("""
        - **Normalize genes** using TCT name resolver
        - **Query TRAPI** APIs for knowledge graph edges
        - **Build graph** with NetworkX
        - **Detect communities** using Louvain algorithm
        - **Visualize** with interactive Cytoscape.js graphs
        """)

    with col2:
        st.subheader(":material/dataset: Example Datasets")
        st.markdown("""
        **Alzheimer's Disease** (15 genes)
        - APOE, APP, PSEN1, PSEN2, MAPT, TREM2, CLU, CR1, BIN1, PICALM, CD33, MS4A6A, ABCA7, SORL1, BACE1
        
        **COVID-19** (10 genes)
        - CD6, IFITM3, IFITM2, STAT5A, KLRG1, DPP4, IL32, PIK3AP1, FYN, IL4R
        """)

# Execute query
if run_query:
    try:
        # Validate input
        validated_genes = validate_gene_list(genes, min_genes=1, max_genes=50)
        
        if disease_curie:
            disease_curie = validate_disease_curie(disease_curie)
        
        # Progress tracking
        progress_bar = st.progress(0)
        status_text = st.empty()
        
        # Step 1: Initialize clients
        status_text.text("Initializing TRAPI client...")
        progress_bar.progress(10)

        client = TRAPIClient(cache_dir=Path("data/cache"))
        builder = GraphBuilder()
        engine = ClusteringEngine()

        # Convert intermediate type selections to biolink categories
        type_mapping = {
            "Protein": "biolink:Protein",
            "ChemicalEntity (Drugs/Metabolites)": "biolink:ChemicalEntity",
            "PhenotypicFeature": "biolink:PhenotypicFeature",
            "Pathway": "biolink:Pathway",
            "BiologicalProcess": "biolink:BiologicalProcess",
            "Gene": "biolink:Gene",
            "AnatomicalEntity": "biolink:AnatomicalEntity",
            "MolecularActivity": "biolink:MolecularActivity",
            "CellularComponent": "biolink:CellularComponent",
        }

        def progress_callback(msg):
            status_text.text(msg)

        # Get predicate filtering settings
        predicate_min_depth = GRANULARITY_PRESETS[predicate_preset]["min_depth"]

        # Handle different query patterns
        if query_pattern == "2-hop (Gene → Intermediate → Disease BioProcesses)" and intermediate_types:
            # New: 2-hop query to disease-associated BiologicalProcesses
            intermediate_categories = [type_mapping[t] for t in intermediate_types if t in type_mapping]

            if not disease_curie:
                raise ValidationError("Disease CURIE is required for BiologicalProcess queries")

            status_text.text("Stage 1: Discovering disease-associated BiologicalProcesses...")
            progress_bar.progress(15)

            # Check for cached BP results if use_cached_bp is enabled
            bp_curies = None
            bp_metadata = None

            if use_cached_bp:
                cached_bp_results = client.list_cached_disease_bp_results(disease_curie)
                if cached_bp_results:
                    # Use the most recent cache
                    most_recent = cached_bp_results[0]
                    try:
                        bp_curies, bp_metadata = client.load_cached_disease_bp_results(most_recent['path'])
                        st.info(f":material/cached: Using cached BiologicalProcesses: {len(bp_curies)} BPs from {most_recent['timestamp_str']}")
                    except Exception as e:
                        st.warning(f":material/warning: Could not load cached BPs: {e}. Running fresh discovery...")
                        bp_curies = None
                        bp_metadata = None

            # Run the full 2-stage query
            response = client.query_gene_to_bioprocesses(
                validated_genes,
                disease_curie=disease_curie,
                intermediate_categories=intermediate_categories,
                bp_curies=bp_curies,
                bp_metadata=bp_metadata,
                filter_disease_bp=filter_bp_predicates,
                predicate_min_depth=predicate_min_depth,
                exclude_literature=exclude_literature,
                exclude_coexpression=exclude_coexpression,
                exclude_homology=exclude_homology,
                timeout_override=600,  # Extended timeout for 2-stage query
                progress_callback=progress_callback,
            )

            progress_bar.progress(60)

            # Show Stage 1 summary
            stage1_meta = response.metadata.get("stage1_metadata", {})
            bp_count = response.metadata.get("bioprocess_count", 0)
            edges_before = stage1_meta.get("total_edges_before_filter", 0)
            edges_after = stage1_meta.get("total_edges_after_filter", 0)

            if edges_before > 0 and filter_bp_predicates:
                st.info(f":material/filter_alt: Stage 1: Found {bp_count} BiologicalProcesses (filtered {edges_before} → {edges_after} edges)")

            edges_removed = response.metadata.get("edges_removed", 0)
            filter_info = f" (filtered {edges_removed} vague relationships)" if edges_removed > 0 else ""
            intermediate_cats = response.metadata.get("intermediate_categories", [])
            intermediate_str = ", ".join([c.replace("biolink:", "") for c in intermediate_cats])

            st.success(f":material/check_circle: 2-hop query: Found {len(response.edges)} edges from {response.apis_succeeded}/{response.apis_queried} APIs{filter_info}")
            st.info(f":material/timeline: Query: Gene → [{intermediate_str}] → {bp_count} Disease BiologicalProcesses")

        elif query_pattern == "2-hop (Gene → Intermediate → Disease)" and intermediate_types:
            # Standard 2-hop query
            intermediate_categories = [type_mapping[t] for t in intermediate_types if t in type_mapping]

            if not disease_curie:
                raise ValidationError("Disease CURIE is required for 2-hop queries")

            status_text.text(f"Querying Translator APIs for {len(validated_genes)} genes...")
            progress_bar.progress(20)

            response = client.query_gene_neighborhood(
                validated_genes,
                disease_curie=disease_curie,
                intermediate_categories=intermediate_categories,
                predicate_min_depth=predicate_min_depth,
                exclude_literature=exclude_literature,
                exclude_coexpression=exclude_coexpression,
                exclude_homology=exclude_homology,
                progress_callback=progress_callback
            )

            progress_bar.progress(60)

            edges_removed = response.metadata.get("edges_removed", 0)
            filter_info = f" (filtered {edges_removed} vague relationships)" if edges_removed > 0 else ""
            intermediate_cats = response.metadata.get("intermediate_categories", [])
            intermediate_str = ", ".join([c.replace("biolink:", "") for c in intermediate_cats])
            st.success(f":material/check_circle: 2-hop query: Found {len(response.edges)} edges from {response.apis_succeeded}/{response.apis_queried} APIs{filter_info}")
            st.info(f":material/timeline: Query: Gene → [{intermediate_str}] → {disease_curie}")

        else:
            # 1-hop neighborhood discovery
            status_text.text(f"Querying Translator APIs for {len(validated_genes)} genes...")
            progress_bar.progress(20)

            response = client.query_gene_neighborhood(
                validated_genes,
                disease_curie=disease_curie,
                intermediate_categories=None,
                predicate_min_depth=predicate_min_depth,
                exclude_literature=exclude_literature,
                exclude_coexpression=exclude_coexpression,
                exclude_homology=exclude_homology,
                progress_callback=progress_callback
            )

            progress_bar.progress(60)

            edges_removed = response.metadata.get("edges_removed", 0)
            filter_info = f" (filtered {edges_removed} vague relationships)" if edges_removed > 0 else ""
            st.success(f":material/check_circle: 1-hop query: Found {len(response.edges)} edges from {response.apis_succeeded}/{response.apis_queried} APIs{filter_info}")
        
        # Step 3: Build graph
        status_text.text("Building knowledge graph...")
        progress_bar.progress(70)

        # Get curie_to_symbol and curie_to_name mappings from response metadata
        curie_to_symbol = response.metadata.get("curie_to_symbol", {})
        curie_to_name = response.metadata.get("curie_to_name", {})

        # Get disease-associated BP CURIEs for triangle rendering (only for Disease BioProcesses pattern)
        disease_bp_curies = None
        if query_pattern == "2-hop (Gene → Intermediate → Disease BioProcesses)":
            disease_bp_curies = response.metadata.get("bioprocess_endpoints", [])

        kg = builder.build_from_trapi_edges(
            response.edges,
            response.input_genes,
            curie_to_symbol,
            curie_to_name,
            disease_bp_curies=disease_bp_curies,
        )
        
        # Calculate gene frequency
        gene_freq = builder.calculate_gene_frequency(kg.graph, response.input_genes)
        for node, freq in gene_freq.items():
            kg.graph.nodes[node]['gene_frequency'] = freq
        
        progress_bar.progress(80)
        
        # Step 4: Cluster
        status_text.text("Detecting communities...")
        results = engine.analyze_graph(kg.graph, response.input_genes)
        
        progress_bar.progress(90)

        # Save to session state
        st.session_state.graph = kg.graph
        st.session_state.knowledge_graph = kg  # Store KnowledgeGraph object
        st.session_state.clustering_results = results
        st.session_state.query_genes = response.input_genes
        st.session_state.response = response
        st.session_state.disease_curie = disease_curie  # Store for visualization sampling
        
        progress_bar.progress(100)
        status_text.markdown(":material/check_circle: Analysis complete!")
        
    except ValidationError as e:
        st.error(f"Validation error: {e}")
    except Exception as e:
        st.error(f"Error during query: {e}")
        import traceback
        st.code(traceback.format_exc())

# Display results
if st.session_state.graph:
    st.divider()
    
    # Tabs for different views
    tab_network, tab_overview, tab_communities = st.tabs([":material/hub: Network", ":material/analytics: Overview", ":material/group_work: Communities"])

    with tab_overview:
        st.header("Analysis Overview")
        
        # Key metrics
        col1, col2, col3, col4 = st.columns(4)

        with col1:
            st.metric("Nodes", st.session_state.graph.number_of_nodes(), border=True)

        with col2:
            st.metric("Edges", st.session_state.graph.number_of_edges(), border=True)

        with col3:
            st.metric("Communities", st.session_state.clustering_results.num_communities, border=True)

        with col4:
            import networkx as nx
            st.metric("Density", f"{nx.density(st.session_state.graph):.4f}", border=True)
        
        # Graph statistics
        st.subheader("Graph Statistics")
        stats_col1, stats_col2 = st.columns(2)
        
        with stats_col1:
            st.json(st.session_state.clustering_results.graph_stats)
        
        with stats_col2:
            st.metric("Modularity", f"{st.session_state.clustering_results.modularity:.3f}", border=True)
            st.caption("Higher modularity indicates better community structure")
        
        # Node categories
        st.subheader("Node Categories")
        if st.session_state.knowledge_graph and st.session_state.knowledge_graph.node_categories:
            category_df = pd.DataFrame([
                {"Category": cat, "Count": count}
                for cat, count in sorted(
                    st.session_state.knowledge_graph.node_categories.items(),
                    key=lambda x: x[1],
                    reverse=True
                )
            ])
            st.dataframe(category_df, width='stretch')
        else:
            st.info("No node category data available")
    
    with tab_network:
        st.header("Knowledge Graph Visualization")

        # Visualization controls - Row 1
        col1, col2, col3, col4 = st.columns([2, 2, 2, 2])

        with col1:
            sizing_metric = st.selectbox(
                "Node Size By",
                ["gene_frequency", "pagerank", "betweenness", "degree"],
                index=0,
                help="Metric used to determine node size"
            )

        with col2:
            layout_options = [
                ("dagre", "Dagre (Hierarchical) - DEFAULT"),
                ("fcose", "fCoSE (Force-Directed)"),
                ("cola", "Cola (Constrained Force)"),
                ("cose-bilkent", "CoSE-Bilkent (Advanced Force)"),
                ("breadthfirst", "Breadth-First (Radial)"),
                ("cose", "CoSE (Basic Force)"),
                ("circle", "Circle"),
                ("grid", "Grid"),
                ("concentric", "Concentric"),
            ]

            layout = st.selectbox(
                "Layout Algorithm",
                options=[opt[0] for opt in layout_options],
                format_func=lambda x: next(opt[1] for opt in layout_options if opt[0] == x),
                index=0,
                help="Dagre creates hierarchical layouts ideal for gene → intermediate → disease flows"
            )

        with col3:
            # Calculate total intermediate nodes (non-query genes)
            if st.session_state.graph:
                query_gene_set = set(st.session_state.query_genes)
                intermediate_count = sum(
                    1 for node in st.session_state.graph.nodes()
                    if node not in query_gene_set and node != st.session_state.disease_curie
                )
                max_slider_value = max(intermediate_count, 10)  # At least 10
                default_value = min(20, intermediate_count)  # Default to 200 or total if less
            else:
                max_slider_value = 500
                default_value = 20

            max_intermediates = st.slider(
                "Top Intermediates",
                min_value=10,
                max_value=max_slider_value,
                value=default_value,
                step=10,
                help="Number of top intermediate nodes to display, ranked by connections to query genes"
            )

        with col4:
            # Build cluster options
            cluster_options = ["All Nodes"]
            if st.session_state.clustering_results:
                for comm in st.session_state.clustering_results.communities:
                    cluster_options.append(f"Cluster {comm.community_id} ({comm.size} nodes)")

            selected_cluster = st.selectbox(
                "Filter by Cluster",
                options=cluster_options,
                index=0,
                help="View all nodes or filter to a specific community"
            )

        # Visualization controls - Row 2 (Node sizing)
        col5, col6, col7 = st.columns([2, 2, 4])

        with col5:
            base_node_size = st.slider(
                "Base Node Size",
                min_value=10,
                max_value=100,
                value=30,
                step=5,
                help="Base size for all nodes in pixels"
            )

        with col6:
            use_metric_sizing = st.checkbox(
                "Use Metric-Based Sizing",
                value=True,
                help="When enabled, node sizes vary by the selected metric. When disabled, all nodes are the same size."
            )

        with col7:
            edge_width = st.slider(
                "Edge Width",
                min_value=1,
                max_value=10,
                value=2,
                step=1,
                help="Width of edges in pixels"
            )

        # Render streamlit-cytoscape visualization
        from biograph_explorer.ui.network_viz import render_network_visualization
        from streamlit_cytoscape import streamlit_cytoscape

        # Determine which graph to visualize based on cluster selection
        if selected_cluster != "All Nodes" and st.session_state.clustering_results:
            # Extract cluster ID from selection (e.g., "Cluster 0 (45 nodes)" -> 0)
            cluster_id = int(selected_cluster.split()[1])

            # Find the selected community
            selected_community = None
            for comm in st.session_state.clustering_results.communities:
                if comm.community_id == cluster_id:
                    selected_community = comm
                    break

            if selected_community:
                # Get cluster nodes
                cluster_nodes = set(selected_community.nodes)

                # Include ALL nodes connected to ANY cluster node
                # This captures both upstream (query genes) and downstream (disease BPs) connections
                nodes_to_include = set(cluster_nodes)
                for node in cluster_nodes:
                    # Add all neighbors (both directions)
                    nodes_to_include.update(st.session_state.graph.predecessors(node))
                    nodes_to_include.update(st.session_state.graph.successors(node))

                # Create subgraph with expanded node set
                display_graph = st.session_state.graph.subgraph(nodes_to_include).copy()

                # Mark which nodes are "native" to cluster vs "connected" (for opacity styling)
                for node in display_graph.nodes():
                    display_graph.nodes[node]['in_cluster'] = node in cluster_nodes

                # Count nodes
                native_count = len(cluster_nodes & set(display_graph.nodes()))
                connected_count = display_graph.number_of_nodes() - native_count

                if connected_count > 0:
                    st.info(f":material/info: Showing Cluster {cluster_id}: {native_count} cluster nodes + {connected_count} connected nodes, {display_graph.number_of_edges()} edges")
                else:
                    st.info(f":material/info: Showing Cluster {cluster_id}: {native_count} nodes, {display_graph.number_of_edges()} edges")
            else:
                display_graph = st.session_state.graph
        else:
            display_graph = st.session_state.graph

        # Prepare visualization data
        viz_data = render_network_visualization(
            display_graph,
            st.session_state.query_genes,
            sizing_metric=sizing_metric,
            layout=layout,
            max_intermediates=max_intermediates,
            disease_curie=st.session_state.disease_curie,
            base_node_size=base_node_size,
            use_metric_sizing=use_metric_sizing,
            edge_width=edge_width,
        )

        if viz_data:
            # Build dynamic key that changes when cluster selection or layout changes
            # This forces Cytoscape to re-run the layout algorithm instead of preserving old positions
            cluster_key = selected_cluster.replace(" ", "_").replace("(", "").replace(")", "")
            cytoscape_key = f"biograph_network_{cluster_key}_{layout}"

            # Render component
            streamlit_cytoscape(
                viz_data["elements"],
                layout=viz_data["layout"],
                node_styles=viz_data["node_styles"],
                edge_styles=viz_data["edge_styles"],
                key=cytoscape_key
            )

            st.caption("""
            **:material/lightbulb: How to explore:**
            - **Drag** to pan • **Scroll** to zoom • **Click** node or edge to select and view information
            - **Fullscreen** button in top-right • **Export JSON** for external tools
            """)
        else:
            st.error("Failed to render visualization")
    
    with tab_communities:
        st.header("Community Detection Results")
        
        st.write(f"Detected **{st.session_state.clustering_results.num_communities} communities** using Louvain algorithm")
        st.write(f"Modularity: **{st.session_state.clustering_results.modularity:.3f}**")
        
        # Display each community
        for comm in st.session_state.clustering_results.communities[:10]:  # Limit to top 10
            with st.expander(f"Community {comm.community_id} ({comm.size} nodes, density={comm.density:.3f})"):
                if comm.top_nodes:
                    st.write("**Top nodes by PageRank:**")
                    for node_info in comm.top_nodes:
                        st.write(f"- {node_info['label']} (PageRank: {node_info['pagerank']:.4f})")
                else:
                    st.write(f"Nodes: {', '.join(comm.nodes[:10])}")
                    if len(comm.nodes) > 10:
                        st.caption(f"... and {len(comm.nodes) - 10} more")

# Footer
st.divider()
st.markdown("**BioGraph Explorer** | Built with Streamlit, NetworkX, TCT, and Cytoscape.js")
