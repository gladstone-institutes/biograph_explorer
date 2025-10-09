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

# Page config
st.set_page_config(
    page_title="BioGraph Explorer",
    page_icon=":material/biotech:",
    layout="wide",
    initial_sidebar_state="expanded"
)

# No custom CSS needed - using Material Design theme from .streamlit/config.toml

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
disease_curie = "MONDO:0004975"  # Default: Alzheimer's

if input_method == "Example Dataset":
    dataset_choice = st.sidebar.selectbox(
        "Select Example",
        ["Alzheimer's Disease (15 genes)", "COVID-19 (10 genes)"]
    )
    
    if dataset_choice == "Alzheimer's Disease (15 genes)":
        csv_path = "data/test_genes/alzheimers_genes.csv"
        disease_curie = "MONDO:0004975"
    else:
        csv_path = "data/test_genes/covid19_genes.csv"
        disease_curie = "MONDO:0100096"
    
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
                    kg = builder.build_from_trapi_edges(
                        cached_response.edges,
                        cached_response.input_genes,
                        curie_to_symbol
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

# Disease CURIE input
disease_curie = st.sidebar.text_input(
    "Disease CURIE",
    value=disease_curie,
    help="Optional disease CURIE (e.g., MONDO:0004975 for Alzheimer's)"
)

# Intermediate entity type selector
st.sidebar.markdown("---")
st.sidebar.subheader(":material/tune: Query Configuration")

# Query pattern selector
query_pattern = st.sidebar.radio(
    "Query Pattern",
    options=[
        "1-hop (Neighborhood Discovery)",
        "2-hop (Gene → Intermediate → Disease)"
    ],
    index=1,  # Default to 2-hop
    help="1-hop finds all direct connections to genes. 2-hop finds therapeutic targets between genes and disease."
)

if query_pattern == "2-hop (Gene → Intermediate → Disease)":
    st.sidebar.markdown("**Path:** Gene → **[Intermediate]** → Disease")

    intermediate_types = st.sidebar.multiselect(
        "Intermediate Entity Types",
        options=[
            "Protein",
            "ChemicalEntity (Drugs/Metabolites)",
            "Gene",
            "PhenotypicFeature",
            "Pathway",
            "BiologicalProcess",
            "AnatomicalEntity",
        ],
        default=["Protein", "ChemicalEntity (Drugs/Metabolites)", "Gene"],
        help="Select intermediate entity types for 2-hop query. These are potential therapeutic targets."
    )

    if not disease_curie:
        st.sidebar.warning(":material/warning: Disease CURIE required for 2-hop queries")
else:
    st.sidebar.markdown("**Path:** Gene → **[Any Connection]**")
    intermediate_types = []  # Empty for 1-hop

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
        
        client = TRAPIClient(cache_dir=Path("data/cache"), max_workers=5)
        builder = GraphBuilder()
        engine = ClusteringEngine()

        # Convert intermediate type selections to biolink categories
        intermediate_categories = None
        if query_pattern == "2-hop (Gene → Intermediate → Disease)" and intermediate_types:
            # Map UI labels to biolink categories
            type_mapping = {
                "Protein": "biolink:Protein",
                "ChemicalEntity (Drugs/Metabolites)": "biolink:ChemicalEntity",
                "PhenotypicFeature": "biolink:PhenotypicFeature",
                "Pathway": "biolink:Pathway",
                "BiologicalProcess": "biolink:BiologicalProcess",
                "Gene": "biolink:Gene",
                "AnatomicalEntity": "biolink:AnatomicalEntity",
            }
            intermediate_categories = [type_mapping[t] for t in intermediate_types if t in type_mapping]

            # Validate disease_curie for 2-hop queries
            if not disease_curie:
                raise ValidationError("Disease CURIE is required for 2-hop queries")

        # Step 2: Query TRAPI
        status_text.text(f"Querying Translator APIs for {len(validated_genes)} genes...")
        progress_bar.progress(20)

        def progress_callback(msg):
            status_text.text(msg)

        response = client.query_gene_neighborhood(
            validated_genes,
            disease_curie=disease_curie,
            intermediate_categories=intermediate_categories,
            progress_callback=progress_callback
        )
        
        progress_bar.progress(60)

        # Show query pattern used
        query_pattern_used = response.metadata.get("query_pattern", "unknown")
        if query_pattern_used == "2-hop":
            intermediate_cats = response.metadata.get("intermediate_categories", [])
            intermediate_str = ", ".join([c.replace("biolink:", "") for c in intermediate_cats])
            st.success(f":material/check_circle: 2-hop query: Found {len(response.edges)} edges from {response.apis_succeeded}/{response.apis_queried} APIs")
            st.info(f":material/timeline: Query: Gene → [{intermediate_str}] → {disease_curie}")
        else:
            st.success(f":material/check_circle: 1-hop query: Found {len(response.edges)} edges from {response.apis_succeeded}/{response.apis_queried} APIs")
        
        # Step 3: Build graph
        status_text.text("Building knowledge graph...")
        progress_bar.progress(70)

        # Get curie_to_symbol mapping from response metadata
        curie_to_symbol = response.metadata.get("curie_to_symbol", {})

        kg = builder.build_from_trapi_edges(response.edges, response.input_genes, curie_to_symbol)
        
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
        status_text.text(":material/check_circle: Analysis complete!")
        
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
    tab1, tab2, tab3 = st.tabs([":material/analytics: Overview", ":material/hub: Network", ":material/group_work: Communities"])
    
    with tab1:
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
    
    with tab2:
        st.header("Knowledge Graph Visualization")

        # Visualization controls
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
                default_value = min(200, intermediate_count)  # Default to 200 or total if less
            else:
                max_slider_value = 500
                default_value = 200

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

        # Render st-link-analysis visualization
        from biograph_explorer.ui.network_viz import render_network_visualization
        from st_link_analysis import st_link_analysis

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
                # Create subgraph with only nodes from this cluster
                cluster_nodes = selected_community.nodes
                display_graph = st.session_state.graph.subgraph(cluster_nodes).copy()

                st.info(f":material/info: Showing Cluster {cluster_id}: {display_graph.number_of_nodes()} nodes, {display_graph.number_of_edges()} edges")
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
        )

        if viz_data:
            # Render component
            st_link_analysis(
                viz_data["elements"],
                layout=viz_data["layout"],
                node_styles=viz_data["node_styles"],
                edge_styles=viz_data["edge_styles"],
                key="biograph_network"
            )

            st.caption("""
            **:material/lightbulb: How to explore:**
            - **Drag** to pan • **Scroll** to zoom • **Click** node to select
            - **Fullscreen** button in top-right • **Export JSON** for external tools
            """)
        else:
            st.error("Failed to render visualization")
    
    with tab3:
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
