"""Network visualization using streamlit-cytoscape component.

Features:
- Interactive Cytoscape.js graph rendering
- Multiple layout algorithms (cose, fcose, circle, grid, etc.)
- Dynamic node sizing by centrality or gene frequency
- Node coloring and icons by category
- Edge styling with predicates
- Material Icons for node types
- Built-in fullscreen and JSON export

Phase 2 Status: Implemented with streamlit-cytoscape
"""

from typing import Optional, List, Dict, Any
import networkx as nx
from networkx.readwrite import json_graph
from pathlib import Path
import logging
import base64

try:
    from streamlit_cytoscape import NodeStyle, EdgeStyle
    STREAMLIT_CYTOSCAPE_AVAILABLE = True
except ImportError:
    STREAMLIT_CYTOSCAPE_AVAILABLE = False

logger = logging.getLogger(__name__)

# Color scheme for node categories
# Colors chosen for visual distinction and accessibility
CATEGORY_COLORS = {
    "Gene": "#1C91D4",           # Blue
    "Disease": "#B14380",        # Magenta/Pink
    "Protein": "#00F6B3",        # Cyan/Teal
    "ChemicalEntity": "#D55E00", # Orange
    "BiologicalProcess": "#D5C711",  # Yellow/Gold
    "AnatomicalEntity": "#c9246f",   # Rose/Pink
    "Pathway": "#8B5CF6",        # Purple (distinct from BiologicalProcess)
    "MolecularActivity": "#10B981",  # Emerald Green
    "CellularComponent": "#F97316",  # Amber/Orange-Red
    "PhenotypicFeature": "#EC4899",  # Pink/Fuchsia
    "Cluster": "#005F45",        # Dark Green
    "Other": "#666666",          # Gray
}
# Get the absolute path to the assets directory
_ASSETS_DIR = Path(__file__).parent.parent / "assets"


def _svg_to_data_uri(svg_path: Path) -> str:
    """Convert SVG file to data URI for embedding in browser.

    Args:
        svg_path: Path to SVG file

    Returns:
        Data URI string (e.g., "url('data:image/svg+xml;base64,...')")
    """
    try:
        with open(svg_path, 'rb') as f:
            svg_data = f.read()
        encoded = base64.b64encode(svg_data).decode('utf-8')
        return f"url('data:image/svg+xml;base64,{encoded}')"
    except Exception as e:
        logger.warning(f"Failed to load SVG icon {svg_path}: {e}")
        return None


# SVG icons for node categories (converted to data URIs)
# New categories use existing icons that semantically match their type
CATEGORY_ICONS = {
    "Gene": _svg_to_data_uri(_ASSETS_DIR / "genetics.svg"),
    "Disease": _svg_to_data_uri(_ASSETS_DIR / "sick.svg"),
    "ChemicalEntity": _svg_to_data_uri(_ASSETS_DIR / "experiment.svg"),
    "Protein": _svg_to_data_uri(_ASSETS_DIR / "polymer.svg"),
    "Pathway": _svg_to_data_uri(_ASSETS_DIR / "arrow_split.svg"),
    "BiologicalProcess": _svg_to_data_uri(_ASSETS_DIR / "arrow_split.svg"),
    "AnatomicalEntity": _svg_to_data_uri(_ASSETS_DIR / "immunology.svg"),
    "MolecularActivity": _svg_to_data_uri(_ASSETS_DIR / "experiment.svg"),  # Chemical/activity
    "CellularComponent": _svg_to_data_uri(_ASSETS_DIR / "immunology.svg"),  # Cell-related
    "PhenotypicFeature": _svg_to_data_uri(_ASSETS_DIR / "health_metrics.svg"),  # Phenotype/health
    "Cluster": None,
    "Other": None,
}


def _format_attribute_value_safe(attr: Dict[str, Any]) -> Optional[str]:
    """Safely format any attribute value for display.

    Args:
        attr: Attribute dictionary with 'value' key

    Returns:
        Formatted string or None if value cannot be displayed
    """
    value = attr.get('value')

    if value is None:
        return None

    # Handle primitive types
    if isinstance(value, (str, int, bool)):
        return str(value)

    if isinstance(value, float):
        return f"{value:.4f}"

    # Handle lists
    if isinstance(value, list):
        # List of primitives - join them
        if all(isinstance(v, (str, int, float, bool)) for v in value):
            return ", ".join(str(v) for v in value)
        # List of objects - just show count
        return f"{len(value)} items"

    # Handle dicts - skip them (too complex)
    if isinstance(value, dict):
        return None

    return None


def _flatten_attributes_for_display(attributes: List[Dict[str, Any]]) -> List[str]:
    """Extract displayable attributes, skipping containers and already-extracted fields.

    Recursively processes nested attribute structures (up to 1 level) and filters
    out attributes that are containers or have already been extracted to dedicated fields.

    Args:
        attributes: List of TRAPI attribute dictionaries

    Returns:
        List of formatted attribute strings (e.g., "knowledge_level: not_provided")
    """
    # Attribute types to skip (containers or already extracted)
    SKIP_TYPES = {
        'biolink:publications',
        'biolink:supporting_text',
        'biolink:extraction_confidence_score',
        'biolink:has_supporting_study_result',  # Container
        'biolink:supporting_study',  # Container
        'biolink:primary_knowledge_source',  # Redundant with sources
        'EDAM:data_0951',  # Statistical data container
    }

    # Attribute names to skip (legacy format)
    SKIP_NAMES = {'publications', 'sentences', 'tmkp_confidence_score', 'tmkp_ids'}

    display_attrs = []

    def process_attrs(attrs: List[Dict[str, Any]], depth: int = 0):
        """Recursively process attributes up to max depth of 1."""
        if depth > 1:  # Max nesting depth in data is 1
            return

        for attr in attrs:
            attr_type = attr.get('attribute_type_id', '')
            attr_name = attr.get('original_attribute_name', '')

            # Skip containers and already-extracted fields
            if attr_type in SKIP_TYPES or attr_name in SKIP_NAMES:
                # But recurse into containers to extract nested displayable attributes
                if 'attributes' in attr and attr.get('attributes'):
                    process_attrs(attr.get('attributes', []), depth + 1)
                continue

            # Format this attribute for display
            formatted = _format_attribute_value_safe(attr)
            if formatted:
                # Use original_attribute_name if available, else use attribute_type_id
                label = attr_name if attr_name else attr_type.replace('biolink:', '')
                display_attrs.append(f"{label}: {formatted}")

            # Recurse if nested (shouldn't happen for non-containers, but be safe)
            if 'attributes' in attr and attr.get('attributes'):
                process_attrs(attr.get('attributes', []), depth + 1)

    process_attrs(attributes)
    return display_attrs


def prepare_cytoscape_elements(
    graph: nx.DiGraph,
    query_genes: List[str],
    sizing_metric: str = "gene_frequency",
    base_node_size: int = 30,
    use_metric_sizing: bool = True,
    edge_width: int = 2,
) -> Dict[str, List[Dict]]:
    """Convert NetworkX graph to Cytoscape.js elements format.

    Uses nx.cytoscape_data() as base, then enriches with custom attributes.
    Internal styling attributes are prefixed with '_' and automatically hidden
    by streamlit-cytoscape's infopanel (controlled via hide_underscore_attrs parameter).

    Args:
        graph: NetworkX DiGraph to convert
        query_genes: List of query gene node IDs (for highlighting)
        sizing_metric: Metric for node sizing (gene_frequency, pagerank, betweenness, degree)
        base_node_size: Base size for nodes in pixels (default: 30)
        use_metric_sizing: If True, size nodes by metric; if False, use uniform size
        edge_width: Width of edges in pixels (default: 2)

    Returns:
        Dictionary with "nodes" and "edges" lists in Cytoscape.js format
    """
    # Get base Cytoscape format from NetworkX
    # First, create a copy without non-serializable, complex, or internal attributes
    graph_copy = graph.copy()
    # Exclude attributes that shouldn't be in Cytoscape element data:
    # - translator_node: TranslatorNode objects can't be JSON serialized
    # - node_annotations: legacy attribute with nested dicts (replaced by annotation_features)
    # - annotation_features: we add clean formatted version to element data instead
    # - is_query_gene, is_disease_associated_bp: internal flags (we add prefixed versions)
    # - curie, value: redundant with 'id'
    exclude_attrs = {
        'translator_node', 'node_annotations', 'annotation_features',
        'is_query_gene', 'is_disease_associated_bp', 'curie', 'value'
    }
    for node in graph_copy.nodes():
        for attr in exclude_attrs:
            if attr in graph_copy.nodes[node]:
                del graph_copy.nodes[node][attr]

    cyto_data = json_graph.cytoscape_data(graph_copy)
    elements = cyto_data["elements"]

    # Calculate metric values for sizing
    metric_values = {}
    if sizing_metric in ["pagerank", "betweenness", "degree"]:
        for node in graph.nodes():
            metric_values[node] = graph.nodes[node].get(sizing_metric, 0)
    elif sizing_metric == "gene_frequency":
        for node in graph.nodes():
            metric_values[node] = graph.nodes[node].get("gene_frequency", 0)
    else:
        # Default to degree
        metric_values = dict(graph.degree())

    # Normalize to 0-1 range
    if metric_values:
        max_val = max(metric_values.values()) if max(metric_values.values()) > 0 else 1
        normalized_metrics = {k: v / max_val for k, v in metric_values.items()}
    else:
        normalized_metrics = {node: 0.5 for node in graph.nodes()}

    # Enrich nodes with custom attributes
    # Attributes prefixed with '_' are internal styling attributes (hidden from infopanel by default)
    for node_element in elements["nodes"]:
        node_id = node_element["data"]["id"]
        node_attrs = graph.nodes[node_id]

        # Add category as label (for NodeStyle matching) - internal, hidden from user
        category = node_attrs.get("category", "Other")
        node_element["data"]["label"] = category  # Used by NodeStyle for matching
        # Also add as 'category' for user display
        node_element["data"]["category"] = category

        # Add display name (for caption)
        original_symbol = node_attrs.get("original_symbol", "")
        display_label = node_attrs.get("label", node_id)
        node_element["data"]["name"] = original_symbol if original_symbol else display_label

        # Add query gene flag (internal - for styling)
        is_query_gene = node_attrs.get("is_query_gene", False)
        node_element["data"]["_is_query_gene"] = is_query_gene

        # Check if this is a disease-associated BiologicalProcess (internal - for styling)
        is_disease_bp = node_attrs.get("is_disease_associated_bp", False)
        node_element["data"]["_is_disease_associated_bp"] = is_disease_bp

        # Set shape: triangles for query genes AND disease-associated BPs (internal - for styling)
        use_triangle = is_query_gene or is_disease_bp
        node_element["data"]["_node_shape"] = "triangle" if use_triangle else "ellipse"

        # Adjust icon sizing/position for triangles (internal - for styling)
        if use_triangle:
            node_element["data"]["_bg_width"] = "55%"
            node_element["data"]["_bg_height"] = "55%"
            node_element["data"]["_bg_position_y"] = "80%"
        else:
            node_element["data"]["_bg_width"] = "70%"
            node_element["data"]["_bg_height"] = "70%"
            node_element["data"]["_bg_position_y"] = "50%"

        # Calculate node size
        if use_metric_sizing:
            # Scale around base_node_size: range is base_node_size ± 15px
            min_size = max(10, base_node_size - 15)
            max_size = base_node_size + 15
            size_range = max_size - min_size
            node_size = min_size + (normalized_metrics.get(node_id, 0.5) * size_range)
        else:
            # Uniform sizing - all nodes same size
            node_size = base_node_size

        # Query-derived nodes (triangles) get larger size (10% bigger, minimum base_node_size)
        if use_triangle:
            node_size = max(node_size * 1.1, base_node_size)

        node_element["data"]["_size"] = round(node_size, 1)

        # Calculate font size proportional to node size (internal - for styling)
        font_size = max(8, min(14, round(node_size * 0.25)))
        node_element["data"]["_font_size"] = font_size

        # Add metrics for tooltip/inspection (user-facing)
        node_element["data"]["gene_frequency"] = node_attrs.get("gene_frequency", 0)
        node_element["data"]["pagerank"] = node_attrs.get("pagerank", 0)
        node_element["data"]["betweenness"] = node_attrs.get("betweenness", 0)
        node_element["data"]["degree"] = graph.degree(node_id)

        # Add cluster membership flag (internal - for opacity styling)
        in_cluster = node_attrs.get("in_cluster", True)
        node_element["data"]["_in_cluster"] = in_cluster
        # Reduce opacity for nodes from other clusters (internal - for styling)
        node_element["data"]["_opacity"] = 1.0 if in_cluster else 0.6

        # Add annotation features as top-level fields for tooltip display (user-facing)
        # Using 'ann_' prefix to avoid conflicts with existing fields
        # Each field is a primitive (string/number) so Cytoscape.js can render it directly
        annotation_features = node_attrs.get('annotation_features', {})
        if annotation_features:
            for key, value in annotation_features.items():
                if isinstance(value, list):
                    # Format lists as readable strings
                    count = len(value)
                    if count <= 3:
                        formatted = ', '.join(str(v) for v in value)
                    else:
                        formatted = f"{', '.join(str(v) for v in value[:3])}... ({count} total)"
                    node_element["data"][f"ann_{key}"] = formatted
                else:
                    node_element["data"][f"ann_{key}"] = value

    # Enrich edges with custom attributes
    # Attributes prefixed with '_' are internal styling attributes (hidden from infopanel by default)
    for idx, edge_element in enumerate(elements["edges"]):
        source = edge_element["data"]["source"]
        target = edge_element["data"]["target"]

        # Add unique ID for streamlit-cytoscape (REQUIRED by Cytoscape.js)
        # Note: 'id' must remain unprefixed as it's a structural requirement, not just display
        edge_element["data"]["id"] = f"e{idx}"

        # Get edge data from graph
        if graph.has_edge(source, target):
            edge_attrs = graph[source][target]

            # Add predicate (cleaned) as label - user-facing
            predicate = edge_attrs.get("predicate", "")
            edge_element["data"]["label"] = predicate.replace("biolink:", "")
            # Store full predicate as internal (redundant with label)
            if predicate:
                edge_element["data"]["_predicate"] = predicate

            # Add full sources information (formatted as readable strings)
            sources = edge_attrs.get("sources", [])
            if sources:
                # Format sources as readable text
                sources_text = []
                for src in sources:
                    resource_id = src.get("resource_id", "").replace("infores:", "")
                    resource_role = src.get("resource_role", "")
                    upstream = src.get("upstream_resource_ids", [])
                    if upstream:
                        upstream_text = ", ".join([u.replace("infores:", "") for u in upstream])
                        sources_text.append(f"{resource_id} ({resource_role}, upstream: {upstream_text})")
                    else:
                        sources_text.append(f"{resource_id} ({resource_role})")
                edge_element["data"]["sources"] = "; ".join(sources_text)

            # Add publications (formatted with clickable links)
            publications = edge_attrs.get("publications", [])
            if publications:
                # Separate PMIDs/PMCIDs from image URLs and already-formatted links
                pmid_links = []
                image_links = []
                for pub in publications:
                    if pub.startswith("https://") or pub.startswith("http://"):
                        # Already a full URL - check if it's an image or keep as-is
                        if any(ext in pub.lower() for ext in ['.jpg', '.png', '.gif', '.svg']):
                            image_links.append(f'<a href="{pub}" target="_blank">Figure</a>')
                        else:
                            # It's already a link, just make it clickable
                            pmid_links.append(f'<a href="{pub}" target="_blank">{pub}</a>')
                    else:
                        # Convert PMID/PMCID to clickable PubMed/PMC link
                        if pub.startswith("PMID:"):
                            pmid = pub.replace("PMID:", "")
                            pmid_links.append(f'<a href="https://pubmed.ncbi.nlm.nih.gov/{pmid}" target="_blank">{pub}</a>')
                        elif pub.startswith("PMC:") or pub.startswith("PMCID:"):
                            # Extract the PMC ID, handling cases like:
                            # PMC:8721920 -> PMC8721920
                            # PMCID:8721920 -> PMC8721920
                            # PMCID:PMC8721920 -> PMC8721920 (already has PMC)
                            # PMC:PMC8721920 -> PMC8721920 (already has PMC)

                            # Remove prefix
                            if pub.startswith("PMCID:"):
                                pmcid = pub[6:]  # Remove "PMCID:"
                            else:  # PMC:
                                pmcid = pub[4:]  # Remove "PMC:"

                            # Ensure it starts with PMC (add if missing)
                            if not pmcid.startswith("PMC"):
                                pmcid = "PMC" + pmcid

                            pmid_links.append(f'<a href="https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}" target="_blank">{pub}</a>')
                        else:
                            pmid_links.append(pub)  # Keep as-is if format unknown

                # Format for display
                pub_parts = []
                if pmid_links:
                    pub_parts.append("Publications: " + ", ".join(pmid_links))
                if image_links:
                    pub_parts.append(f"Figures: {', '.join(image_links)}")

                edge_element["data"]["publications"] = " | ".join(pub_parts) if pub_parts else ""

            # Add supporting sentences (formatted as joined text)
            sentences = edge_attrs.get("sentences", [])
            if sentences:
                # Ensure all items are strings and join with line breaks
                sentence_strs = []
                for s in sentences:
                    if isinstance(s, str) and s:
                        sentence_strs.append(s)
                    elif isinstance(s, list):
                        # Handle lists of sentences
                        sentence_strs.extend([str(item) for item in s if item])
                    elif s:
                        sentence_strs.append(str(s))

                if sentence_strs:
                    edge_element["data"]["sentences"] = " | ".join(sentence_strs)

            # Add confidence scores (formatted as key: value pairs)
            confidence_scores = edge_attrs.get("confidence_scores", {})
            if confidence_scores and isinstance(confidence_scores, dict):
                scores_text = []
                for key, value in confidence_scores.items():
                    # Format the key to be more readable
                    readable_key = key.replace("tmkp_", "").replace("_", " ").title()
                    if isinstance(value, float):
                        scores_text.append(f"{readable_key}: {value:.4f}")
                    else:
                        scores_text.append(f"{readable_key}: {value}")

                if scores_text:
                    edge_element["data"]["confidence_scores"] = ", ".join(scores_text)
                else:
                    # Empty dict - don't show the field at all
                    if "confidence_scores" in edge_element["data"]:
                        del edge_element["data"]["confidence_scores"]
            else:
                # Not a dict or empty - remove if exists
                if "confidence_scores" in edge_element["data"]:
                    del edge_element["data"]["confidence_scores"]

            # Add qualifiers (formatted as readable text)
            qualifiers = edge_attrs.get("qualifiers", [])
            if qualifiers:
                qualifiers_text = []
                for q in qualifiers:
                    q_type = q.get("qualifier_type_id", "").replace("biolink:", "")
                    q_value = q.get("qualifier_value", "").replace("biolink:", "")
                    qualifiers_text.append(f"{q_type}: {q_value}")

                if qualifiers_text:
                    edge_element["data"]["qualifiers"] = "; ".join(qualifiers_text)
                else:
                    # Empty qualifiers - remove the field
                    if "qualifiers" in edge_element["data"]:
                        del edge_element["data"]["qualifiers"]
            else:
                # No qualifiers - remove if exists
                if "qualifiers" in edge_element["data"]:
                    del edge_element["data"]["qualifiers"]

            # Format qualifiers as readable text for display (only if qualifiers exist)
            if qualifiers:
                qualified_predicate = ""
                object_direction = ""
                object_aspect = ""
                for qualifier in qualifiers:
                    q_type = qualifier.get("qualifier_type_id", "")
                    q_value = qualifier.get("qualifier_value", "")
                    if q_type == "biolink:qualified_predicate":
                        qualified_predicate = q_value.replace("biolink:", "")
                    elif q_type == "biolink:object_direction_qualifier":
                        object_direction = q_value
                    elif q_type == "biolink:object_aspect_qualifier":
                        object_aspect = q_value

                # Create human-readable qualifier text only if complete
                if qualified_predicate and object_direction and object_aspect:
                    edge_element["data"]["qualified_relationship"] = f"{qualified_predicate} {object_direction} {object_aspect}"

            # Add complete attributes array (formatted as readable text) - internal/debug
            attributes = edge_attrs.get("attributes", [])
            if attributes:
                # Format attributes using helper function
                attributes_text = _flatten_attributes_for_display(attributes)
                if attributes_text:  # Only add if non-empty
                    edge_element["data"]["_attributes"] = "; ".join(attributes_text)

            # Add query_result_id only if present - internal
            query_result_id = edge_attrs.get("query_result_id")
            if query_result_id is not None:
                edge_element["data"]["_query_result_id"] = query_result_id

        # Add edge width and scaled font size - internal styling
        edge_element["data"]["_edge_width"] = edge_width
        # Scale font: 8px at width 1, up to 14px at width 10
        edge_font_size = max(8, min(14, 6 + edge_width))
        edge_element["data"]["_edge_font_size"] = edge_font_size

    # Note: Internal styling attributes are prefixed with '_' and automatically
    # hidden by streamlit-cytoscape's infopanel (hide_underscore_attrs=True by default).
    # The debug_mode parameter can be passed to streamlit_cytoscape() to control this.

    return {"nodes": elements["nodes"], "edges": elements["edges"]}


def create_node_styles(graph: nx.DiGraph) -> List["NodeStyle"]:
    """Create NodeStyle objects for each category in the graph.

    Args:
        graph: NetworkX graph containing nodes with "category" attribute

    Returns:
        List of NodeStyle objects, one per unique category
    """
    if not STREAMLIT_CYTOSCAPE_AVAILABLE:
        raise ImportError("streamlit-cytoscape not installed. Run: pip install streamlit-cytoscape")

    # Extract unique categories
    categories = set()
    for node in graph.nodes():
        category = graph.nodes[node].get("category", "Other")
        categories.add(category)

    # Create NodeStyle for each category
    # custom_styles enables dynamic node sizing from node data
    node_styles = []
    for category in sorted(categories):
        color = CATEGORY_COLORS.get(category, CATEGORY_COLORS["Other"])
        icon = CATEGORY_ICONS.get(category)

        # Use "name" as caption (displays gene symbol or node label)
        # custom_styles reads size from node data for dynamic sizing
        # Label styling: dynamic font size + background box for visibility
        # Note: Internal styling attributes use underscore prefix (e.g., _size, _node_shape)
        label_styles = {
            "font-size": "data(_font_size)",
            "text-background-color": "#ffffff",
            "text-background-opacity": 0.95,
            "text-background-padding": "2px",
            "text-background-shape": "roundrectangle",
        }

        if icon:
            node_styles.append(
                NodeStyle(
                    label=category,
                    color=color,
                    caption="name",
                    icon=icon,
                    custom_styles={
                        "width": "data(_size)",
                        "height": "data(_size)",
                        "shape": "data(_node_shape)",
                        "background-width": "data(_bg_width)",
                        "background-height": "data(_bg_height)",
                        "background-position-y": "data(_bg_position_y)",
                        "opacity": "data(_opacity)",
                        **label_styles,
                    }
                )
            )
        else:
            node_styles.append(
                NodeStyle(
                    label=category,
                    color=color,
                    caption="name",
                    custom_styles={
                        "width": "data(_size)",
                        "height": "data(_size)",
                        "shape": "data(_node_shape)",
                        "opacity": "data(_opacity)",
                        **label_styles,
                    }
                )
            )

    return node_styles


def create_edge_styles(graph: nx.DiGraph) -> List["EdgeStyle"]:
    """Create EdgeStyle objects for each predicate in the graph.

    Args:
        graph: NetworkX graph containing edges with "predicate" attribute

    Returns:
        List of EdgeStyle objects, one per unique predicate
    """
    if not STREAMLIT_CYTOSCAPE_AVAILABLE:
        raise ImportError("streamlit-cytoscape not installed. Run: pip install streamlit-cytoscape")

    # Extract unique predicates
    predicates = set()
    for _, _, edge_attrs in graph.edges(data=True):
        predicate = edge_attrs.get("predicate", "").replace("biolink:", "")
        if predicate:
            predicates.add(predicate)

    # Edge styling: dynamic width and font size
    # Note: Internal styling attributes use underscore prefix
    edge_custom_styles = {
        "width": "data(_edge_width)",
        "font-size": "data(_edge_font_size)",
    }

    # Create EdgeStyle for each predicate
    edge_styles = []
    for predicate in sorted(predicates):
        edge_styles.append(
            EdgeStyle(
                label=predicate,     # Matches edge data["label"]
                caption="label",     # Display the label attribute
                directed=True,       # Show arrows
                custom_styles=edge_custom_styles,
            )
        )

    # If no predicates found, add a default style
    if not edge_styles:
        edge_styles.append(
            EdgeStyle(
                label="default",     # Default label
                caption="label",     # Display the label attribute
                directed=True,       # Show arrows
                custom_styles=edge_custom_styles,
            )
        )

    return edge_styles


def get_layout_config(layout_name: str = "dagre") -> Dict[str, Any]:
    """Get layout configuration for Cytoscape.js.

    Args:
        layout_name: Name of layout algorithm (dagre, fcose, cola, cose-bilkent, breadthfirst, cose, circle, grid, concentric)

    Returns:
        Layout configuration dictionary
    """
    # Default configuration based on demos
    layout = {
        "name": layout_name,
        "animate": "end",
        "nodeDimensionsIncludeLabels": False
    }

    # Add layout-specific options
    if layout_name == "dagre":
        layout.update({
            "rankDir": "TB",  # Top to bottom (genes → intermediates → disease)
            "ranker": "network-simplex",  # Optimal ranking algorithm
            "nodeSep": 50,  # Separation between nodes on same rank
            "edgeSep": 10,  # Separation between edges
            "rankSep": 75,  # Separation between ranks (levels)
            "fit": True,
            "padding": 30
        })
    elif layout_name == "cose":
        layout.update({
            "nodeRepulsion": 400000,
            "idealEdgeLength": 100,
            "edgeElasticity": 100,
            "nestingFactor": 5,
            "gravity": 80,
            "numIter": 1000,
            "initialTemp": 200,
            "coolingFactor": 0.95,
            "minTemp": 1.0
        })
    elif layout_name == "fcose":
        layout.update({
            "quality": "default",
            "randomize": True,
            "animate": "end",
            "fit": True,
            "padding": 30,
            "nodeSeparation": 75,
            "idealEdgeLength": 50,
            "edgeElasticity": 0.45,
            "nestingFactor": 0.1,
            "gravity": 0.25,
            "numIter": 2500,
            "tile": True,
            "tilingPaddingVertical": 10,
            "tilingPaddingHorizontal": 10,
            "gravityRangeCompound": 1.5,
            "gravityCompound": 1.0,
            "gravityRange": 3.8
        })
    elif layout_name == "cola":
        layout.update({
            "animate": True,
            "refresh": 1,
            "maxSimulationTime": 4000,
            "ungrabifyWhileSimulating": False,
            "fit": True,
            "padding": 30,
            "nodeDimensionsIncludeLabels": False,
            "randomize": False,
            "avoidOverlap": True,
            "handleDisconnected": True,
            "convergenceThreshold": 0.01,
            "nodeSpacing": 10,
            "flow": None,
            "alignment": None,
            "gapInequalities": None
        })
    elif layout_name == "cose-bilkent":
        layout.update({
            "quality": "default",
            "nodeDimensionsIncludeLabels": False,
            "refresh": 30,
            "fit": True,
            "padding": 30,
            "randomize": True,
            "nodeRepulsion": 4500,
            "idealEdgeLength": 50,
            "edgeElasticity": 0.45,
            "nestingFactor": 0.1,
            "gravity": 0.25,
            "numIter": 2500,
            "tile": True,
            "tilingPaddingVertical": 10,
            "tilingPaddingHorizontal": 10,
            "gravityRangeCompound": 1.5,
            "gravityCompound": 1.0,
            "gravityRange": 3.8,
            "initialEnergyOnIncremental": 0.5
        })

    return layout


def render_network_visualization(
    graph: nx.DiGraph,
    query_genes: List[str],
    sizing_metric: str = "gene_frequency",
    layout: str = "dagre",
    max_intermediates: int = 200,
    highlight_nodes: Optional[List[str]] = None,
    disease_curie: Optional[str] = None,
    base_node_size: int = 30,
    use_metric_sizing: bool = True,
    edge_width: int = 2,
) -> Dict[str, Any]:
    """Prepare network visualization data for streamlit-cytoscape component.

    Args:
        graph: NetworkX DiGraph to visualize
        query_genes: List of query gene node IDs (highlighted)
        sizing_metric: Node sizing metric (gene_frequency, pagerank, betweenness, degree)
        layout: Layout algorithm name (cose, fcose, circle, grid, breadthfirst, concentric)
        max_intermediates: Maximum intermediate nodes to display (default: 200)
        highlight_nodes: Optional list of node IDs to highlight (for citations)
        disease_curie: Optional disease CURIE (always included if present)
        base_node_size: Base size for nodes in pixels (default: 30)
        use_metric_sizing: If True, size nodes by metric; if False, use uniform size
        edge_width: Width of edges in pixels (default: 2)

    Returns:
        Dictionary with elements, node_styles, edge_styles, and layout config

    Example:
        >>> viz_data = render_network_visualization(graph, ["NCBIGene:1803"], max_intermediates=200)
        >>> st_link_analysis(viz_data["elements"], viz_data["layout"],
        ...                  viz_data["node_styles"], viz_data["edge_styles"])
    """
    if not STREAMLIT_CYTOSCAPE_AVAILABLE:
        logger.error("streamlit-cytoscape not installed")
        return None

    if graph.number_of_nodes() == 0:
        logger.warning("Empty graph - nothing to visualize")
        return None

    try:
        # Sample graph by top intermediates
        graph = sample_graph_for_visualization(
            graph,
            query_genes,
            max_intermediates=max_intermediates,
            disease_curie=disease_curie
        )

        # Prepare Cytoscape elements
        elements = prepare_cytoscape_elements(
            graph, query_genes, sizing_metric, base_node_size, use_metric_sizing, edge_width
        )

        # Create node and edge styles
        node_styles = create_node_styles(graph)
        edge_styles = create_edge_styles(graph)

        # Get layout configuration
        layout_config = get_layout_config(layout)

        return {
            "elements": elements,
            "node_styles": node_styles,
            "edge_styles": edge_styles,
            "layout": layout_config
        }

    except Exception as e:
        logger.error(f"Error creating visualization: {e}")
        import traceback
        traceback.print_exc()
        return None


def sample_graph_for_visualization(
    graph: nx.DiGraph,
    query_genes: List[str],
    max_intermediates: int = 200,
    disease_curie: Optional[str] = None,
) -> nx.DiGraph:
    """Sample graph by filtering top intermediate nodes ranked by query gene connectivity.

    Sampling strategy:
    1. Always include ALL query genes
    2. Identify intermediate nodes (non-query genes)
    3. Score intermediates by number of edges connecting to query genes
    4. Keep top-N intermediates based on max_intermediates parameter
    5. Include disease node if present
    6. Build subgraph with selected nodes + their edges

    Args:
        graph: Full NetworkX graph
        query_genes: Query gene node IDs (always included)
        max_intermediates: Maximum number of intermediate nodes to include (default: 200)
        disease_curie: Optional disease CURIE (always included if present)

    Returns:
        Sampled subgraph with query genes + top intermediates by connectivity
    """
    # Identify query genes present in graph
    genes_in_graph = set(g for g in query_genes if g in graph.nodes())

    if not genes_in_graph:
        logger.warning("No query genes found in graph")
        return graph

    logger.info(f"Found {len(genes_in_graph)} query genes in graph")

    # Identify intermediate nodes (exclude query genes and disease)
    intermediate_nodes = []
    for node in graph.nodes():
        if node not in genes_in_graph and node != disease_curie:
            intermediate_nodes.append(node)

    total_intermediates = len(intermediate_nodes)
    logger.info(f"Found {total_intermediates} intermediate nodes")

    # If we have fewer intermediates than the limit, return full graph
    if total_intermediates <= max_intermediates:
        logger.info(f"Total intermediates ({total_intermediates}) <= max ({max_intermediates}) - returning full graph")
        return graph

    logger.info(f"Sampling intermediates: {total_intermediates} → {max_intermediates}")

    # Score each intermediate by connections to query genes
    intermediate_scores = {}
    for node in intermediate_nodes:
        # Count edges to/from query genes
        score = 0
        for gene in genes_in_graph:
            if graph.has_edge(gene, node):
                score += 1
            if graph.has_edge(node, gene):
                score += 1

        intermediate_scores[node] = score

    # Sort intermediates by score (descending) and take top-N
    top_intermediates = sorted(
        intermediate_scores.items(),
        key=lambda x: x[1],
        reverse=True
    )[:max_intermediates]

    logger.info(
        f"Top intermediate scores - max: {top_intermediates[0][1] if top_intermediates else 0}, "
        f"min: {top_intermediates[-1][1] if top_intermediates else 0}"
    )

    # Build set of nodes to include
    nodes_to_include = set(genes_in_graph)  # Always include query genes
    nodes_to_include.update(node for node, _ in top_intermediates)  # Add top intermediates

    # Always include disease if present
    if disease_curie and disease_curie in graph.nodes():
        nodes_to_include.add(disease_curie)
        logger.info(f"Including disease node: {disease_curie}")

    # Create subgraph with selected nodes
    sampled_graph = graph.subgraph(nodes_to_include).copy()

    logger.info(
        f"Sampled graph: {sampled_graph.number_of_nodes()} nodes, "
        f"{sampled_graph.number_of_edges()} edges "
        f"(query genes: {len(genes_in_graph)}, intermediates: {len(top_intermediates)})"
    )

    return sampled_graph


def get_node_details(node_id: str, graph: nx.DiGraph) -> Dict[str, Any]:
    """Extract detailed information about a node.

    Args:
        node_id: Node identifier (CURIE)
        graph: NetworkX graph containing the node

    Returns:
        Dictionary with node details including edges and sources
    """
    if node_id not in graph.nodes():
        return {"error": f"Node {node_id} not found in graph"}

    attrs = graph.nodes[node_id]

    # Get edges
    in_edges = list(graph.in_edges(node_id, data=True))
    out_edges = list(graph.out_edges(node_id, data=True))

    # Group edges by predicate
    edges_by_predicate = {}
    all_sources = set()

    for src, tgt, data in in_edges + out_edges:
        predicate = data.get("predicate", "unknown").replace("biolink:", "")
        sources = data.get("knowledge_source", [])

        if predicate not in edges_by_predicate:
            edges_by_predicate[predicate] = []

        edge_info = {
            "direction": "incoming" if tgt == node_id else "outgoing",
            "source" if tgt == node_id else "target": src if tgt == node_id else tgt,
            "source_label" if tgt == node_id else "target_label": graph.nodes[src if tgt == node_id else tgt].get(
                "label", src if tgt == node_id else tgt
            ),
            "sources": sources,
        }
        edges_by_predicate[predicate].append(edge_info)

        # Collect all sources
        all_sources.update(sources)

    return {
        "node_id": node_id,
        "label": attrs.get("label", node_id),
        "original_symbol": attrs.get("original_symbol", ""),
        "category": attrs.get("category", "Unknown"),
        "is_query_gene": attrs.get("is_query_gene", False),
        "metrics": {
            "gene_frequency": attrs.get("gene_frequency", 0),
            "pagerank": attrs.get("pagerank", 0),
            "betweenness": attrs.get("betweenness", 0),
            "degree": graph.degree(node_id),
            "in_degree": graph.in_degree(node_id),
            "out_degree": graph.out_degree(node_id),
        },
        "edges_by_predicate": edges_by_predicate,
        "total_edges": len(in_edges) + len(out_edges),
        "knowledge_sources": list(all_sources),
    }


def filter_graph_by_annotations(
    graph: nx.DiGraph,
    annotation_filters: Dict[str, List[Any]],
) -> nx.DiGraph:
    """Filter graph by annotation attributes.

    Filtering strategy:
    1. Find nodes matching ALL selected annotation criteria (AND logic)
    2. Include only matched nodes plus disease nodes (always preserved)
    3. Edges between included nodes are preserved

    Args:
        graph: Full NetworkX graph
        annotation_filters: Dict of filter_key -> list of selected values.
            Example: {"go_bp": ["cell cycle"], "go_mf": ["ATP binding"]}

    Returns:
        Filtered subgraph with 'in_filter' attribute marking matched nodes
    """
    if not annotation_filters:
        return graph

    # Remove empty filters
    active_filters = {k: v for k, v in annotation_filters.items() if v}
    if not active_filters:
        return graph

    logger.info(f"Filtering graph by annotations: {active_filters}")

    # Find nodes matching ALL filters (AND logic)
    matched_nodes = set()
    for node in graph.nodes():
        node_attrs = graph.nodes[node]
        matches_all = True

        for filter_key, selected_values in active_filters.items():
            # Get node's annotation feature value for this filter
            annotation_features = node_attrs.get('annotation_features', {})
            node_value = annotation_features.get(filter_key)

            if node_value is None:
                matches_all = False
                break

            # Handle list-valued attributes (e.g., go_bp, go_mf, go_cc, alias)
            if isinstance(node_value, list):
                # Check if ANY of the node's values match ANY selected value
                node_value_strs = {str(v) for v in node_value}
                if not node_value_strs.intersection(selected_values):
                    matches_all = False
                    break
            else:
                # Scalar value - check for exact match
                if str(node_value) not in selected_values:
                    matches_all = False
                    break

        if matches_all:
            matched_nodes.add(node)

    # Always include disease nodes (they don't have GO annotations but are important)
    disease_nodes = {
        node for node in graph.nodes()
        if graph.nodes[node].get('category') == 'Disease'
    }

    # Identify query genes
    query_genes = {
        node for node in graph.nodes()
        if graph.nodes[node].get('is_query_gene', False)
    }

    if not matched_nodes:
        logger.warning("No nodes match annotation filters - returning empty graph")
        return nx.DiGraph()

    # Build node set:
    # 1. Matched nodes (intermediates that pass the filter)
    # 2. Disease nodes (always included)
    # 3. Query genes connected to matched nodes
    nodes_to_include = set(matched_nodes)
    nodes_to_include.update(disease_nodes)

    # Add query genes that are connected to matched intermediate nodes
    for node in matched_nodes:
        for neighbor in graph.predecessors(node):
            if neighbor in query_genes:
                nodes_to_include.add(neighbor)
        for neighbor in graph.successors(node):
            if neighbor in query_genes:
                nodes_to_include.add(neighbor)

    # Create subgraph with selected nodes
    filtered_graph = graph.subgraph(nodes_to_include).copy()

    # Mark nodes for styling (matched vs connected neighbor)
    for node in filtered_graph.nodes():
        filtered_graph.nodes[node]['in_filter'] = node in matched_nodes

    logger.info(
        f"Annotation filter result: {len(matched_nodes)} matched nodes, "
        f"{filtered_graph.number_of_nodes()} total nodes (with connections), "
        f"{filtered_graph.number_of_edges()} edges"
    )

    return filtered_graph


def filter_graph_by_category(
    graph: nx.DiGraph,
    selected_category: str,
    query_genes: List[str],
    disease_curie: Optional[str] = None,
) -> nx.DiGraph:
    """Filter graph to show only intermediates of the selected category.

    Strategy:
    1. Find all nodes matching the selected category (intermediates)
    2. Include query genes ONLY if they connect to at least one category node
    3. Always include disease node if present and connected
    4. Preserve edges between included nodes

    Args:
        graph: Full NetworkX graph
        selected_category: Category to filter by (e.g., "Protein", "BiologicalProcess")
        query_genes: Query gene CURIEs (included only if connected to category)
        disease_curie: Optional disease CURIE (always included if present)

    Returns:
        Filtered subgraph with only selected category intermediates
    """
    if not selected_category:
        return graph

    query_gene_set = set(query_genes) if query_genes else set()

    # First, find all nodes matching the selected category
    category_nodes = set()
    for node in graph.nodes():
        if graph.nodes[node].get('category') == selected_category:
            category_nodes.add(node)

    if not category_nodes:
        logger.warning(f"No nodes match category '{selected_category}'")
        return graph

    nodes_to_include = set(category_nodes)

    # Include query genes only if they connect to at least one category node
    for query_gene in query_gene_set:
        if query_gene not in graph:
            continue
        # Check if this query gene has any edge to/from a category node
        neighbors = set(graph.predecessors(query_gene)) | set(graph.successors(query_gene))
        if neighbors & category_nodes:
            nodes_to_include.add(query_gene)

    # Always include disease node if present
    if disease_curie and disease_curie in graph:
        nodes_to_include.add(disease_curie)

    filtered_graph = graph.subgraph(nodes_to_include).copy()

    logger.info(
        f"Category filter '{selected_category}': {filtered_graph.number_of_nodes()} nodes, "
        f"{filtered_graph.number_of_edges()} edges"
    )

    return filtered_graph
