"""BioGraph Explorer - Streamlit Application

Main Streamlit app for multi-gene TRAPI query integration with NetworkX clustering.
"""

import streamlit as st
from pathlib import Path
import pandas as pd
import logging
import json
from datetime import datetime
from typing import Dict, Any

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

from biograph_explorer.core import TRAPIClient, GraphBuilder
from biograph_explorer.utils import validate_gene_list, validate_disease_curie, ValidationError
from biograph_explorer.utils.biolink_predicates import (
    GRANULARITY_PRESETS,
    get_allowed_predicates_for_display,
)


def notify(message: str, msg_type: str = "info", icon: str = None):
    """Show toast notification and persist to session state for Overview tab.

    Args:
        message: The message text
        msg_type: One of "success", "info", "warning", "error"
        icon: Optional material icon name (without :material/: wrapper)
    """
    default_icons = {
        "success": "check_circle",
        "info": "info",
        "warning": "warning",
        "error": "error",
    }
    icon = icon or default_icons.get(msg_type, "info")

    # Show toast
    st.toast(message, icon=f":material/{icon}:")

    # Persist to session state
    if 'query_messages' not in st.session_state:
        st.session_state.query_messages = []
    st.session_state.query_messages.append({
        "type": msg_type,
        "message": message,
        "timestamp": datetime.now(),
        "icon": icon,
    })


def _build_query_export(response) -> Dict[str, Any]:
    """Build comprehensive query export for download.

    Includes:
    - TRAPI query_graph structure
    - Input parameters (genes, disease, intermediate categories)
    - Predicate filter settings
    - API timings and success rates
    - Category mismatch statistics (if any)
    """
    metadata = response.metadata

    return {
        "query_id": response.query_id,
        "timestamp": response.timestamp.isoformat() if response.timestamp else None,
        "input": {
            "gene_symbols": metadata.get("gene_symbols", []),
            "normalized_genes": metadata.get("normalized_genes", {}),
            "disease_curie": response.target_disease,
            "intermediate_categories": metadata.get("intermediate_categories", []),
            "query_pattern": metadata.get("query_pattern", ""),
        },
        "trapi_query": metadata.get("query_json", {}),
        "predicate_settings": {
            "predicate_min_depth": metadata.get("predicate_min_depth"),
            "exclude_literature": metadata.get("exclude_literature"),
            "exclude_coexpression": metadata.get("exclude_coexpression"),
            "exclude_homology": metadata.get("exclude_homology"),
            "predicates_used_count": metadata.get("predicates_used", 0),
        },
        "results": {
            "total_edges": len(response.edges),
            "edges_before_filter": metadata.get("edges_before_filter", 0),
            "edges_removed": metadata.get("edges_removed", 0),
            "apis_queried": response.apis_queried,
            "apis_succeeded": response.apis_succeeded,
        },
        "api_timings": metadata.get("api_timings", []),
        "category_mismatch_stats": metadata.get("category_mismatch_stats", {}),
        "stage1_metadata": metadata.get("stage1_metadata"),
        "bioprocess_count": metadata.get("bioprocess_count"),
    }


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
if 'query_genes' not in st.session_state:
    st.session_state.query_genes = []
if 'response' not in st.session_state:
    st.session_state.response = None
if 'selected_node' not in st.session_state:
    st.session_state.selected_node = None
if 'disease_curie' not in st.session_state:
    st.session_state.disease_curie = None
if 'annotation_metadata' not in st.session_state:
    st.session_state.annotation_metadata = None
if 'hpa_metadata' not in st.session_state:
    st.session_state.hpa_metadata = None
if 'annotation_filters' not in st.session_state:
    st.session_state.annotation_filters = {}
if 'debug_mode' not in st.session_state:
    st.session_state.debug_mode = False
if 'query_messages' not in st.session_state:
    st.session_state.query_messages = []
if 'category_mismatch_detail' not in st.session_state:
    st.session_state.category_mismatch_detail = None
if 'collapse_counter' not in st.session_state:
    st.session_state.collapse_counter = 0

# Title
st.markdown("### :material/biotech: BioGraph Explorer")
st.caption("Multi-gene TRAPI query integration with NetworkX graph analysis and visualization")

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

                    # Restore node annotations from cache OR fetch from API
                    status_text.text("Loading node annotations...")
                    if cached_response.node_annotations:
                        # Restore annotations from cache
                        for node, features in cached_response.node_annotations.items():
                            if node in kg.graph.nodes:
                                kg.graph.nodes[node]['annotation_features'] = features
                        st.session_state.annotation_metadata = cached_response.annotation_metadata
                        st.session_state.hpa_metadata = cached_response.hpa_metadata
                        logger.info(f"Restored {len(cached_response.node_annotations)} annotations from cache")
                    else:
                        # Annotations not in cache - fetch from API
                        try:
                            from biograph_explorer.core import NodeAnnotator
                            annotator = NodeAnnotator(cache_dir=Path("data/cache/annotations"))
                            kg.graph, annotation_metadata = annotator.annotate_graph(
                                kg.graph,
                                fields="all",
                                progress_callback=lambda msg: status_text.text(f"Annotating: {msg}")
                            )
                            st.session_state.annotation_metadata = annotation_metadata

                            # Save annotations to cache for next time
                            node_annotations = {}
                            for node in kg.graph.nodes():
                                features = kg.graph.nodes[node].get('annotation_features', {})
                                if features:
                                    node_annotations[node] = features
                            cached_response.node_annotations = node_annotations
                            cached_response.annotation_metadata = annotation_metadata
                            client._cache_response(cached_response)
                            logger.info(f"Fetched and cached {len(node_annotations)} annotations")
                        except Exception as e:
                            logger.warning(f"Node annotation failed (non-fatal): {e}")
                            st.session_state.annotation_metadata = None

                    progress_bar.progress(90)

                    # Save to session state
                    st.session_state.graph = kg.graph
                    st.session_state.knowledge_graph = kg  # Store KnowledgeGraph object
                    st.session_state.query_genes = cached_response.input_genes
                    st.session_state.response = cached_response
                    st.session_state.disease_curie = cached_response.target_disease

                    # Validate publication extraction (detect duplicate edges)
                    from biograph_explorer.utils.publication_utils import validate_publication_extraction
                    pub_validation = validate_publication_extraction(cached_response.edges)
                    st.session_state.pub_validation = pub_validation
                    if pub_validation["has_issues"]:
                        notify(
                            f"Duplicate edges detected: {pub_validation['edges_lost_to_collisions']} duplicates found "
                            f"({pub_validation['publications_at_risk']} publications on duplicates). "
                            f"See Query Log for details.",
                            msg_type="info",
                            icon="info"
                        )

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

# API Timeout slider (in minutes for user-friendliness)
api_timeout_minutes = st.sidebar.slider(
    "API Timeout (minutes)",
    min_value=1,
    max_value=15,
    value=10,
    step=1,
    help="Timeout for each TRAPI API request. Increase if queries are timing out."
)
api_timeout = api_timeout_minutes * 60  # Convert to seconds for client

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

# Debug mode toggle
st.sidebar.markdown("---")
st.sidebar.checkbox(
    "Debug mode",
    key="debug_mode",
    help="Show all node/edge attributes including internal styling data"
)

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
    # Clear previous query messages
    st.session_state.query_messages = []
    st.session_state.category_mismatch_detail = None

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

        client = TRAPIClient(cache_dir=Path("data/cache"), timeout=api_timeout)
        builder = GraphBuilder()

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
                        notify(f"Using cached BiologicalProcesses: {len(bp_curies)} BPs from {most_recent['timestamp_str']}", "info", "cached")
                    except Exception as e:
                        notify(f"Could not load cached BPs: {e}. Running fresh discovery...", "warning")
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
                notify(f"Stage 1: Found {bp_count} BiologicalProcesses (filtered {edges_before} -> {edges_after} edges)", "info", "filter_alt")

            edges_removed = response.metadata.get("edges_removed", 0)
            filter_info = f" (filtered {edges_removed} vague relationships)" if edges_removed > 0 else ""
            intermediate_cats = response.metadata.get("intermediate_categories", [])
            intermediate_str = ", ".join([c.replace("biolink:", "") for c in intermediate_cats])

            notify(f"2-hop query: Found {len(response.edges)} edges from {response.apis_succeeded}/{response.apis_queried} APIs{filter_info}", "success")
            notify(f"Query: Gene -> [{intermediate_str}] -> {bp_count} Disease BiologicalProcesses", "info", "timeline")

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
            notify(f"2-hop query: Found {len(response.edges)} edges from {response.apis_succeeded}/{response.apis_queried} APIs{filter_info}", "success")
            notify(f"Query: Gene -> [{intermediate_str}] -> {disease_curie}", "info", "timeline")

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
            notify(f"1-hop query: Found {len(response.edges)} edges from {response.apis_succeeded}/{response.apis_queried} APIs{filter_info}", "success")

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

        progress_bar.progress(75)

        # Step 3.5: Annotate nodes with metadata from Node Annotator API
        status_text.text("Annotating nodes with metadata...")
        try:
            from biograph_explorer.core import NodeAnnotator
            annotator = NodeAnnotator(cache_dir=Path("data/cache/annotations"))
            kg.graph, annotation_metadata = annotator.annotate_graph(
                kg.graph,
                fields="all",
                progress_callback=lambda msg: status_text.text(f"Annotating: {msg}")
            )
            st.session_state.annotation_metadata = annotation_metadata
            annotated_pct = annotation_metadata.get('annotation_rate', 0) * 100
            logger.info(f"Node annotation complete: {annotated_pct:.1f}% annotated")

            # Step 3.6: Add HPA cell type specificity annotations (BEFORE cache save)
            status_text.text("Fetching HPA cell type specificity data...")
            try:
                kg.graph, hpa_metadata = annotator.annotate_with_hpa(
                    kg.graph,
                    progress_callback=lambda msg: status_text.text(f"HPA: {msg}")
                )
                st.session_state.hpa_metadata = hpa_metadata
                hpa_count = hpa_metadata.get('hpa_annotated_count', 0)
                logger.info(f"HPA annotation complete: {hpa_count} genes with cell type data")
            except Exception as e:
                logger.warning(f"HPA annotation failed (non-fatal): {e}")
                st.session_state.hpa_metadata = None

            # Save ALL annotations (including HPA) to response for caching
            node_annotations = {}
            for node in kg.graph.nodes():
                features = kg.graph.nodes[node].get('annotation_features', {})
                if features:
                    node_annotations[node] = features
            response.node_annotations = node_annotations
            response.annotation_metadata = annotation_metadata
            response.hpa_metadata = st.session_state.hpa_metadata

            # Re-save cache with all annotations included
            status_text.text("Saving annotations to cache...")
            client._cache_response(response)
            logger.info(f"Cache updated with {len(node_annotations)} node annotations")

        except Exception as e:
            logger.warning(f"Node annotation failed (non-fatal): {e}")
            st.session_state.annotation_metadata = None

        progress_bar.progress(90)

        # Save to session state
        st.session_state.graph = kg.graph
        st.session_state.knowledge_graph = kg  # Store KnowledgeGraph object
        st.session_state.query_genes = response.input_genes
        st.session_state.response = response
        st.session_state.disease_curie = disease_curie  # Store for visualization sampling
        
        progress_bar.progress(100)
        status_text.markdown(":material/check_circle: Analysis complete!")

        # Validate publication extraction (detect duplicate edges)
        from biograph_explorer.utils.publication_utils import validate_publication_extraction
        pub_validation = validate_publication_extraction(response.edges)
        st.session_state.pub_validation = pub_validation
        if pub_validation["has_issues"]:
            notify(
                f"Duplicate edges detected: {pub_validation['edges_lost_to_collisions']} duplicates found "
                f"({pub_validation['publications_at_risk']} publications on duplicates). "
                f"See Query Log for details.",
                msg_type="info",
                icon="info"
            )

        # Show category mismatch warning if applicable
        category_stats = response.metadata.get("category_mismatch_stats", {})
        if category_stats.get("mismatched_count", 0) > 0:
            mismatch_count = category_stats["mismatched_count"]
            total_count = category_stats.get("total_intermediates", 0)
            actual_counts = category_stats.get("actual_category_counts", {})
            requested = category_stats.get("requested_categories", [])

            # Format category names without biolink: prefix for readability
            actual_str = ", ".join([f"{k.replace('biolink:', '')}: {v}" for k, v in actual_counts.items()])
            requested_str = ", ".join([r.replace("biolink:", "") for r in requested])

            # Short toast + store detail for Overview tab
            notify(f"Category mismatch: {mismatch_count}/{total_count} intermediates don't match requested categories", "warning")
            st.session_state.category_mismatch_detail = {
                "mismatch_count": mismatch_count,
                "total_count": total_count,
                "requested": requested_str,
                "actual": actual_str,
            }

    except ValidationError as e:
        notify(f"Validation error: {e}", "error")
    except Exception as e:
        notify(f"Error during query: {e}", "error")
        import traceback
        st.code(traceback.format_exc())

# Display results
if st.session_state.graph:
    st.divider()
    
    # Tabs for different views
    tab_network, tab_overview = st.tabs([":material/hub: Network", ":material/analytics: Overview"])

    with tab_overview:
        st.markdown("#### Analysis Overview")
        
        # Key metrics
        col1, col2, col3 = st.columns(3)

        with col1:
            st.metric("Nodes", st.session_state.graph.number_of_nodes(), border=True)

        with col2:
            st.metric("Edges", st.session_state.graph.number_of_edges(), border=True)

        with col3:
            import networkx as nx
            st.metric("Density", f"{nx.density(st.session_state.graph):.4f}", border=True)

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

        # Publication Evidence section
        st.subheader("Publication Evidence")
        if st.session_state.graph and st.session_state.graph.number_of_edges() > 0:
            from biograph_explorer.utils.publication_utils import (
                get_publication_frequency,
                get_publication_frequency_by_category,
                format_publication_display,
            )
            from biograph_explorer.ui.network_viz import CATEGORY_COLORS

            # Get frequency by category for stacked bars
            pub_freq = get_publication_frequency(st.session_state.graph)
            pub_freq_by_cat = get_publication_frequency_by_category(st.session_state.graph)

            if pub_freq_by_cat:
                # Get unique categories for filter dropdown
                all_categories_in_pubs = set()
                for cat_counts in pub_freq_by_cat.values():
                    all_categories_in_pubs.update(cat_counts.keys())
                all_categories_in_pubs = sorted(all_categories_in_pubs)

                # Category filter for bar chart
                col_filter, col_spacer = st.columns([2, 4])
                with col_filter:
                    bar_category_filter = st.selectbox(
                        "Filter by Category",
                        options=["All Categories"] + all_categories_in_pubs,
                        index=0,
                        key="pub_bar_category_filter",
                        help="Filter bar chart to show only edges of specific intermediate type"
                    )

                # Calculate totals for sorting (filtered if category selected)
                pub_totals = {}
                for pub_id, cat_counts in pub_freq_by_cat.items():
                    if bar_category_filter == "All Categories":
                        pub_totals[pub_id] = sum(cat_counts.values())
                    elif bar_category_filter in cat_counts:
                        pub_totals[pub_id] = cat_counts[bar_category_filter]

                # Get top 10 publications (only those with counts after filtering)
                top_pubs = sorted(
                    [(k, v) for k, v in pub_totals.items() if v > 0],
                    key=lambda x: x[1],
                    reverse=True
                )[:10]
                top_pub_ids = [p[0] for p in top_pubs]

                if top_pub_ids:
                    # Build DataFrame for stacked bar chart
                    rows = []
                    for pub_id in top_pub_ids:
                        cat_counts = pub_freq_by_cat[pub_id]
                        for category, count in cat_counts.items():
                            if bar_category_filter == "All Categories" or category == bar_category_filter:
                                rows.append({
                                    "Publication": format_publication_display(pub_id),
                                    "Category": category,
                                    "Edge Count": count,
                                })

                    pub_df = pd.DataFrame(rows)

                    st.caption(f"Top {len(top_pub_ids)} Most Cited Publications")
                    import altair as alt
                    chart = alt.Chart(pub_df).mark_bar().encode(
                        x=alt.X('Edge Count:Q', title='Edge Count', axis=alt.Axis(tickMinStep=1)),
                        y=alt.Y('Publication:N', sort='-x', title=None),
                        color=alt.Color(
                            'Category:N',
                            scale=alt.Scale(
                                domain=list(CATEGORY_COLORS.keys()),
                                range=list(CATEGORY_COLORS.values())
                            ),
                            legend=alt.Legend(title="Intermediate Type")
                        ),
                        tooltip=['Publication', 'Category', 'Edge Count'],
                        order=alt.Order('Category:N')
                    ).properties(
                        height=min(300, len(top_pub_ids) * 30)
                    )
                    st.altair_chart(chart, width="stretch")
                else:
                    st.info("No publications found for selected category")

                # Summary metrics
                col_pub1, col_pub2 = st.columns(2)
                with col_pub1:
                    st.metric("Unique Publications", len(pub_freq), border=True)
                with col_pub2:
                    # Handle both DiGraph and MultiDiGraph edge access
                    graph = st.session_state.graph
                    if isinstance(graph, nx.MultiDiGraph):
                        edges_with_pubs = sum(
                            1 for u, v, key, data in graph.edges(keys=True, data=True)
                            if data.get('publications')
                        )
                    else:
                        edges_with_pubs = sum(
                            1 for u, v, data in graph.edges(data=True)
                            if data.get('publications')
                        )
                    st.metric("Edges with Publications", edges_with_pubs, border=True)
            else:
                st.info("No publication data found in edges")
        else:
            st.info("No graph data available for publication analysis")

        # Query Details and Download
        st.subheader("Query Details")
        if st.session_state.response:
            query_export = _build_query_export(st.session_state.response)
            query_json_str = json.dumps(query_export, indent=2, default=str)

            col_download, col_preview = st.columns([1, 3])
            with col_download:
                st.download_button(
                    label=":material/download: Download Query JSON",
                    data=query_json_str,
                    file_name=f"biograph_query_{st.session_state.response.query_id}.json",
                    mime="application/json",
                    help="Download TRAPI query structure, parameters, and API timings"
                )

            with st.expander("Preview Query JSON", expanded=False):
                st.json(query_export)
        else:
            st.info("No query data available")

        # Query Log section - display persisted messages
        st.subheader("Query Log")

        if st.session_state.query_messages:
            # Group messages by type for organized display
            messages_by_type = {"error": [], "warning": [], "success": [], "info": []}
            for msg in st.session_state.query_messages:
                messages_by_type[msg["type"]].append(msg)

            # Display errors first (most important)
            for msg in messages_by_type["error"]:
                st.error(f":material/{msg['icon']}: {msg['message']}")

            # Then warnings
            for msg in messages_by_type["warning"]:
                st.warning(f":material/{msg['icon']}: {msg['message']}")

            # Category mismatch detail (if exists)
            if st.session_state.category_mismatch_detail:
                detail = st.session_state.category_mismatch_detail
                with st.expander("Category Mismatch Details"):
                    st.markdown(f"**Requested:** {detail['requested']}")
                    st.markdown(f"**Actual:** {detail['actual']}")
                    st.caption("This may be due to knowledge providers conflating similar entity types (e.g., Gene/Protein).")

            # Publication extraction debug details (if duplicate edges exist)
            if st.session_state.get('pub_validation') and st.session_state.pub_validation.get('has_issues'):
                pub_val = st.session_state.pub_validation
                with st.expander("Duplicate Edge Details", expanded=False):
                    st.markdown("**ℹ️ Duplicate edges detected in TRAPI data**")
                    st.markdown(
                        "Multiple edges with identical (subject, object, predicate) tuples were found. "
                        "All edges are preserved in the MultiDiGraph, but this may indicate upstream data issues."
                    )

                    # Summary stats
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        st.metric("Raw TRAPI Edges", pub_val.get('total_edges', 0))
                    with col2:
                        st.metric("Unique (S,O,P) Tuples", pub_val.get('unique_node_pairs', 0))
                    with col3:
                        st.metric("Duplicate Edges", pub_val.get('edges_lost_to_collisions', 0))

                    pubs_on_dups = pub_val.get('publications_at_risk', 0)
                    if pubs_on_dups > 0:
                        st.info(f"**{pubs_on_dups} publications** are on duplicate edges (all preserved)")

                    # Show sample duplicates
                    samples = pub_val.get('sample_collisions', [])
                    if samples:
                        st.markdown("---")
                        st.markdown("**Sample duplicate edge groups:**")
                        for i, sample in enumerate(samples, 1):
                            with st.container(border=True):
                                st.markdown(f"**Duplicate {i}:** `{sample['subject']}` → `{sample['object']}`")
                                st.markdown(f"**Predicate:** `{sample.get('predicate', 'unknown')}`")
                                st.markdown(f"**{sample.get('duplicate_count', sample.get('edge_count', 0))} instances** of this edge:")

                                for pub_info in sample.get('pubs_by_edge', []):
                                    st.markdown(f"- {pub_info}")

                    st.caption(
                        "**Note:** All edges are preserved in MultiDiGraph. Duplicates may indicate "
                        "the same relationship reported by multiple knowledge sources."
                    )

            # Then success messages
            for msg in messages_by_type["success"]:
                st.success(f":material/{msg['icon']}: {msg['message']}")

            # Finally info messages in an expander (less critical)
            if messages_by_type["info"]:
                with st.expander(f"Info Messages ({len(messages_by_type['info'])})", expanded=False):
                    for msg in messages_by_type["info"]:
                        st.info(f":material/{msg['icon']}: {msg['message']}")

            # Clear button
            if st.button(":material/delete: Clear Query Log", key="clear_query_log"):
                st.session_state.query_messages = []
                st.session_state.category_mismatch_detail = None
                st.session_state.pub_validation = None
                st.rerun()
        else:
            st.caption("No query messages to display")

    with tab_network:
        st.markdown("#### Knowledge Graph Visualization")

        # Visualization controls - Row 1
        col1, col2, col3, col4, col5 = st.columns([2, 2, 2, 2, 2])

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
            # Category filter dropdown
            if st.session_state.knowledge_graph and st.session_state.knowledge_graph.node_categories:
                all_categories = list(st.session_state.knowledge_graph.node_categories.keys())
                all_categories.sort()
            else:
                all_categories = []

            category_filter = st.selectbox(
                "Intermediate Category",
                options=["All Categories"] + all_categories,
                index=0,
                help="Filter to show only intermediates of a specific type"
            )

        with col5:
            # Publication filter dropdown - filtered by selected category
            from biograph_explorer.utils.publication_utils import (
                get_publication_frequency_by_category,
                format_publication_display,
            )

            if st.session_state.graph:
                # Get publication counts broken down by category
                pub_by_category = get_publication_frequency_by_category(st.session_state.graph)

                if pub_by_category:
                    # Filter publications based on selected category
                    if category_filter and category_filter != "All Categories":
                        # Only show publications that have edges in the selected category
                        filtered_pubs = {
                            pub_id: cat_counts.get(category_filter, 0)
                            for pub_id, cat_counts in pub_by_category.items()
                            if category_filter in cat_counts
                        }
                    else:
                        # Show all publications with total counts
                        filtered_pubs = {
                            pub_id: sum(cat_counts.values())
                            for pub_id, cat_counts in pub_by_category.items()
                        }

                    # Sort by frequency (descending)
                    sorted_pubs = sorted(filtered_pubs.items(), key=lambda x: x[1], reverse=True)
                    # Format: "PMID 12345 (7 edges)"
                    pub_options = [
                        f"{format_publication_display(pub_id)} ({count} edges)"
                        for pub_id, count in sorted_pubs
                    ]
                    # Keep raw IDs for filtering
                    pub_ids = [pub_id for pub_id, _ in sorted_pubs]
                else:
                    pub_options = []
                    pub_ids = []
            else:
                pub_options = []
                pub_ids = []

            if pub_options:
                pub_filter_display = st.selectbox(
                    "Filter by Publication",
                    options=["All Publications"] + pub_options,
                    index=0,
                    help="Show only edges citing this publication"
                )
                # Extract raw pub_id from selection
                if pub_filter_display != "All Publications":
                    pub_filter_idx = pub_options.index(pub_filter_display)
                    pub_filter = pub_ids[pub_filter_idx]
                else:
                    pub_filter = None
            else:
                pub_filter = None
                no_pub_msg = "No publications in this category" if (category_filter and category_filter != "All Categories") else "No publications available"
                st.selectbox(
                    "Filter by Publication",
                    options=[no_pub_msg],
                    disabled=True,
                    help="No publication metadata found in current graph"
                )

        # Visualization controls - Row 2 (Node sizing)
        col6, col7, col8, col9 = st.columns([2, 2, 3, 1])

        with col6:
            base_node_size = st.slider(
                "Base Node Size",
                min_value=10,
                max_value=100,
                value=30,
                step=5,
                help="Base size for all nodes in pixels"
            )

        with col7:
            use_metric_sizing = st.checkbox(
                "Use Metric-Based Sizing",
                value=True,
                help="When enabled, node sizes vary by the selected metric. When disabled, all nodes are the same size."
            )

        with col8:
            edge_width = st.slider(
                "Edge Width",
                min_value=1,
                max_value=10,
                value=2,
                step=1,
                help="Width of edges in pixels"
            )

        with col9:
            # Add collapse button to allow users to re-collapse expanded edges
            if st.button(":material/unfold_less: Collapse", help="Re-collapse all expanded parallel edges"):
                st.session_state.collapse_counter += 1

        # Visualization controls - Row 3 (Annotation Filters)
        # Dynamically discover filterable attributes from annotation metadata
        if st.session_state.annotation_metadata:
            filterable_attrs = st.session_state.annotation_metadata.get('filterable_attributes', {})
            searchable_attrs = st.session_state.annotation_metadata.get('searchable_attributes', {})

            has_filters = filterable_attrs or searchable_attrs
            if has_filters:
                with st.expander(":material/filter_alt: Annotation Filters", expanded=False):
                    # Section 1: GO Term Filters (searchable multiselect)
                    if searchable_attrs:
                        st.markdown("**GO Term Filters**")
                        st.caption("Filter nodes by Gene Ontology terms (nodes with ANY selected term will match)")

                        go_cols = st.columns(min(3, len(searchable_attrs)))
                        go_labels = {
                            'go_bp': 'Biological Process',
                            'go_mf': 'Molecular Function',
                            'go_cc': 'Cellular Component',
                        }

                        for i, (go_key, go_values) in enumerate(searchable_attrs.items()):
                            col_idx = i % len(go_cols)
                            with go_cols[col_idx]:
                                label = go_labels.get(go_key, go_key)
                                selected_go = st.multiselect(
                                    f"{label} ({len(go_values)} terms)",
                                    options=go_values,
                                    default=[],
                                    key=f"go_filter_{go_key}",
                                    help="Type to search. Nodes with ANY selected term will match.",
                                    placeholder="Search GO terms..."
                                )
                                if selected_go:
                                    st.session_state.annotation_filters[go_key] = selected_go
                                elif go_key in st.session_state.annotation_filters:
                                    del st.session_state.annotation_filters[go_key]

                    # Section 2: HPA Cell Type Specificity Filters
                    hpa_meta = getattr(st.session_state, 'hpa_metadata', None)
                    if hpa_meta and st.session_state.graph:
                        # Collect HPA filter options
                        hpa_filter_options = {
                            'hpa_cell_type_specificity': hpa_meta.get('hpa_specificity_categories', []),
                            'hpa_immune_cell_specificity': hpa_meta.get('hpa_immune_specificities', []),
                            'hpa_tissue_specificity': hpa_meta.get('hpa_tissue_specificities', []),
                        }

                        # Collect unique cell types from graph nodes (not in metadata)
                        all_cell_types = set()
                        all_immune_cells = set()
                        all_tissues = set()
                        for node in st.session_state.graph.nodes():
                            features = st.session_state.graph.nodes[node].get('annotation_features', {})
                            for ct, _ in features.get('hpa_top_cell_types', []):
                                all_cell_types.add(ct)
                            for ic, _ in features.get('hpa_top_immune_cells', []):
                                all_immune_cells.add(ic)
                            for t, _ in features.get('hpa_top_tissues', []):
                                all_tissues.add(t)

                        hpa_filter_options['hpa_top_cell_types'] = sorted(all_cell_types)
                        hpa_filter_options['hpa_top_immune_cells'] = sorted(all_immune_cells)
                        hpa_filter_options['hpa_top_tissues'] = sorted(all_tissues)

                        # Only show HPA section if there are options
                        has_hpa_options = any(hpa_filter_options.values())
                        if has_hpa_options:
                            st.markdown("---")
                            st.markdown("**HPA Cell Type Specificity Filters**")

                            # Row 1: Category dropdowns
                            cat_cols = st.columns(3)
                            with cat_cols[0]:
                                cell_spec_options = hpa_filter_options.get('hpa_cell_type_specificity', [])
                                if cell_spec_options:
                                    selected_cell_spec = st.selectbox(
                                        "Cell Type Specificity",
                                        options=["All"] + cell_spec_options,
                                        key="hpa_cell_type_specificity_filter",
                                        help="Filter by HPA cell type specificity category"
                                    )
                                    if selected_cell_spec != "All":
                                        st.session_state.annotation_filters['hpa_cell_type_specificity'] = [selected_cell_spec]
                                    elif 'hpa_cell_type_specificity' in st.session_state.annotation_filters:
                                        del st.session_state.annotation_filters['hpa_cell_type_specificity']

                            with cat_cols[1]:
                                immune_spec_options = hpa_filter_options.get('hpa_immune_cell_specificity', [])
                                if immune_spec_options:
                                    selected_immune_spec = st.selectbox(
                                        "Immune Cell Specificity",
                                        options=["All"] + immune_spec_options,
                                        key="hpa_immune_cell_specificity_filter",
                                        help="Filter by HPA immune/blood cell specificity"
                                    )
                                    if selected_immune_spec != "All":
                                        st.session_state.annotation_filters['hpa_immune_cell_specificity'] = [selected_immune_spec]
                                    elif 'hpa_immune_cell_specificity' in st.session_state.annotation_filters:
                                        del st.session_state.annotation_filters['hpa_immune_cell_specificity']

                            with cat_cols[2]:
                                tissue_spec_options = hpa_filter_options.get('hpa_tissue_specificity', [])
                                if tissue_spec_options:
                                    selected_tissue_spec = st.selectbox(
                                        "Tissue Specificity",
                                        options=["All"] + tissue_spec_options,
                                        key="hpa_tissue_specificity_filter",
                                        help="Filter by HPA tissue specificity category"
                                    )
                                    if selected_tissue_spec != "All":
                                        st.session_state.annotation_filters['hpa_tissue_specificity'] = [selected_tissue_spec]
                                    elif 'hpa_tissue_specificity' in st.session_state.annotation_filters:
                                        del st.session_state.annotation_filters['hpa_tissue_specificity']

                            # Row 2: Searchable cell type multiselects
                            st.caption("Filter by specific cell types (nodes with ANY selected type will match)")
                            search_cols = st.columns(3)

                            with search_cols[0]:
                                cell_type_options = hpa_filter_options.get('hpa_top_cell_types', [])
                                if cell_type_options:
                                    selected_cell_types = st.multiselect(
                                        f"Cell Types ({len(cell_type_options)})",
                                        options=cell_type_options,
                                        default=[],
                                        key="hpa_cell_types_filter",
                                        help="Search for specific cell types (e.g., Hepatocytes, Microglia)",
                                        placeholder="Search cell types..."
                                    )
                                    if selected_cell_types:
                                        st.session_state.annotation_filters['hpa_top_cell_types'] = selected_cell_types
                                    elif 'hpa_top_cell_types' in st.session_state.annotation_filters:
                                        del st.session_state.annotation_filters['hpa_top_cell_types']

                            with search_cols[1]:
                                immune_cell_options = hpa_filter_options.get('hpa_top_immune_cells', [])
                                if immune_cell_options:
                                    selected_immune_cells = st.multiselect(
                                        f"Immune Cells ({len(immune_cell_options)})",
                                        options=immune_cell_options,
                                        default=[],
                                        key="hpa_immune_cells_filter",
                                        help="Search for specific immune cell types",
                                        placeholder="Search immune cells..."
                                    )
                                    if selected_immune_cells:
                                        st.session_state.annotation_filters['hpa_top_immune_cells'] = selected_immune_cells
                                    elif 'hpa_top_immune_cells' in st.session_state.annotation_filters:
                                        del st.session_state.annotation_filters['hpa_top_immune_cells']

                            with search_cols[2]:
                                tissue_options = hpa_filter_options.get('hpa_top_tissues', [])
                                if tissue_options:
                                    selected_tissues = st.multiselect(
                                        f"Tissues ({len(tissue_options)})",
                                        options=tissue_options,
                                        default=[],
                                        key="hpa_tissues_filter",
                                        help="Search for specific tissues",
                                        placeholder="Search tissues..."
                                    )
                                    if selected_tissues:
                                        st.session_state.annotation_filters['hpa_top_tissues'] = selected_tissues
                                    elif 'hpa_top_tissues' in st.session_state.annotation_filters:
                                        del st.session_state.annotation_filters['hpa_top_tissues']

                    # Show active filters summary
                    active_count = sum(1 for v in st.session_state.annotation_filters.values() if v)
                    if active_count > 0:
                        st.caption(f":material/filter_list: {active_count} filter(s) active - connected nodes will also be shown")

        # Render streamlit-cytoscape visualization
        from biograph_explorer.ui.network_viz import render_network_visualization, filter_graph_by_annotations, filter_graph_by_category, filter_graph_by_publication
        from streamlit_cytoscape import streamlit_cytoscape

        display_graph = st.session_state.graph

        # Apply annotation filters if any are active
        if st.session_state.annotation_filters:
            active_filters = {k: v for k, v in st.session_state.annotation_filters.items() if v}
            if active_filters:
                pre_filter_count = display_graph.number_of_nodes()
                display_graph = filter_graph_by_annotations(
                    display_graph,
                    active_filters,
                )
                post_filter_count = display_graph.number_of_nodes()
                if post_filter_count < pre_filter_count:
                    matched_count = sum(
                        1 for n in display_graph.nodes()
                        if display_graph.nodes[n].get('in_filter', False)
                    )
                    st.info(
                        f":material/filter_alt: Annotation filter: {matched_count} matched nodes + "
                        f"{post_filter_count - matched_count} connected nodes shown"
                    )

        # Apply category filter if not "All Categories"
        if category_filter and category_filter != "All Categories":
            pre_cat_count = display_graph.number_of_nodes()
            display_graph = filter_graph_by_category(
                display_graph,
                category_filter,
                st.session_state.query_genes,
                st.session_state.disease_curie
            )
            post_cat_count = display_graph.number_of_nodes()
            if post_cat_count == 0:
                st.warning(f":material/filter_alt_off: No nodes match the category '{category_filter}'.")
            elif post_cat_count < pre_cat_count:
                # Count intermediates and connected query genes separately
                query_gene_set = set(st.session_state.query_genes) if st.session_state.query_genes else set()
                category_intermediate_count = sum(
                    1 for node in display_graph.nodes()
                    if node not in query_gene_set
                    and node != st.session_state.disease_curie
                    and display_graph.nodes[node].get('category') == category_filter
                )
                connected_query_gene_count = sum(
                    1 for node in display_graph.nodes()
                    if node in query_gene_set
                )
                gene_suffix = f" + {connected_query_gene_count} connected query genes" if connected_query_gene_count > 0 else ""
                st.info(f":material/category: Showing {category_intermediate_count} {category_filter} intermediates{gene_suffix}")

        # Apply publication filter if selected (filters to relevant intermediates + highlights edges)
        highlighted_count = 0  # Track for meta-edge coloring
        if pub_filter:
            pre_pub_nodes = display_graph.number_of_nodes()
            display_graph = filter_graph_by_publication(
                display_graph,
                pub_filter,
                st.session_state.query_genes,
                st.session_state.disease_curie
            )
            # Count highlighted edges
            if isinstance(display_graph, nx.MultiDiGraph):
                highlighted_count = sum(
                    1 for _, _, _, d in display_graph.edges(keys=True, data=True)
                    if d.get('_has_filtered_pub')
                )
            else:
                highlighted_count = sum(
                    1 for _, _, d in display_graph.edges(data=True)
                    if d.get('_has_filtered_pub')
                )
            total_edges = display_graph.number_of_edges()
            if display_graph.number_of_nodes() == 0:
                st.warning(
                    f":material/filter_alt_off: No nodes match the selected filters. "
                    f"The publication '{format_publication_display(pub_filter)}' may not exist in the filtered category."
                )
            else:
                st.info(
                    f":material/highlight: {highlighted_count} edges citing "
                    f"{format_publication_display(pub_filter)} (teal), "
                    f"{total_edges - highlighted_count} context edges (gray)"
                )

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
            # Build dynamic key that changes when layout or filters change
            # This forces Cytoscape to re-run the layout algorithm instead of preserving old positions
            # EXCLUDE styling-only parameters (base_node_size, use_metric_sizing, edge_width, collapse_counter)
            # to preserve manual node positions when only visual properties change

            # Include annotation filter state in key so layout re-runs when filters change
            filter_state = st.session_state.annotation_filters
            if filter_state:
                # Create a deterministic hash from active filters
                filter_items = sorted((k, tuple(sorted(v)) if isinstance(v, list) else v)
                                      for k, v in filter_state.items() if v)
                filter_hash = hash(tuple(filter_items)) if filter_items else 0
            else:
                filter_hash = 0

            pub_filter_hash = hash(pub_filter) if pub_filter else 0
            # Only include graph structure and layout parameters in key
            # Styling parameters (node size, edge width, collapse_counter) excluded to preserve positions
            cytoscape_key = (
                f"biograph_network_{layout}_{filter_hash}_{category_filter}_{pub_filter_hash}_"
                f"{max_intermediates}_{sizing_metric}"
                # Note: base_node_size, use_metric_sizing, edge_width, collapse_counter excluded to preserve positions
            )

            # Compute priority predicate for edge collapsing
            # When publication filter is active, prioritize predicates from edges with the publication
            # Otherwise, use the most specific predicate (highest biolink depth) as the meta-edge label
            from biograph_explorer.utils.biolink_predicates import get_predicate_depths
            all_predicates = set()
            filtered_pub_predicates = set()  # Predicates from edges with the filtered publication

            if isinstance(display_graph, nx.MultiDiGraph):
                for _, _, _, data in display_graph.edges(keys=True, data=True):
                    pred = data.get("predicate", "").replace("biolink:", "")
                    if pred:
                        all_predicates.add(pred)
                        # Track predicates from publication-filtered edges
                        if data.get('_has_filtered_pub'):
                            filtered_pub_predicates.add(pred)
            else:
                for _, _, data in display_graph.edges(data=True):
                    pred = data.get("predicate", "").replace("biolink:", "")
                    if pred:
                        all_predicates.add(pred)
                        if data.get('_has_filtered_pub'):
                            filtered_pub_predicates.add(pred)

            # Use publication-filtered predicates when pub_filter is active, else all predicates
            predicates_to_use = filtered_pub_predicates if (pub_filter and filtered_pub_predicates) else all_predicates

            # Sort by depth (most specific first), then alphabetically for ties
            if predicates_to_use:
                depths = get_predicate_depths()
                sorted_predicates = sorted(
                    predicates_to_use,
                    key=lambda p: (-depths.get(p.lower().replace(" ", "_"), 0), p)
                )
                priority_label = sorted_predicates[0]
            else:
                priority_label = None

            # Render component with edge collapsing enabled
            # When debug_mode is True, show all attributes (hide_underscore_attrs=False)
            # Meta-edges inherit styling from EdgeStyle (width, font-size are set statically there)
            # This keeps meta_edge_style stable across slider changes, preserving collapse state

            # Note: Edge highlighting for publication filter is handled via data selector
            # 'edge[_has_filtered_pub = true]' in get_highlighted_edge_style()
            # This approach colors only edges with the specific publication, not all edges
            # with a matching predicate.

            # Meta-edge styling: only set properties that don't change with sliders
            # Width and font-size are inherited from EdgeStyle to keep this dict stable
            meta_edge_style = {
                "text-rotation": "autorotate",  # Align label with edge line
            }

            streamlit_cytoscape(
                viz_data["elements"],
                layout=viz_data["layout"],
                node_styles=viz_data["node_styles"],
                edge_styles=viz_data["edge_styles"],
                key=cytoscape_key,
                hide_underscore_attrs=not st.session_state.debug_mode,
                # Edge collapsing: disable when publication filter is active to show individual edge highlighting
                edge_actions=["collapse", "expand"],
                collapse_parallel_edges=(not pub_filter),
                priority_edge_label=priority_label,
                meta_edge_style=meta_edge_style,
            )

            st.caption("""
            **:material/lightbulb: How to explore:**
            - **Drag** to pan • **Scroll** to zoom • **Click** node or edge to select and view information
            - **Double-click** collapsed edge to expand parallel edges • **Fullscreen** in top-right
            """)
        else:
            st.error("Failed to render visualization")

# Footer
st.divider()
st.markdown("**BioGraph Explorer** | Built with Streamlit, NetworkX, TCT, and Cytoscape.js")
