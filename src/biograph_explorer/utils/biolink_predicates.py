"""Biolink predicate hierarchy utilities for granularity-based filtering.

This module provides functions to:
- Load and cache the biolink model from GitHub (version-aware)
- Parse predicate hierarchy from is_a relationships
- Filter predicates by depth (granularity level)
- Get human-readable predicate lists for UI display

The biolink model defines predicates in a tree structure rooted at 'related_to'.
Higher depth = more specific predicates (e.g., 'causes' is more specific than 'affects').

Granularity Presets:
- All Relationships (Level 1): Include everything except root 'related_to'
- Standard (Level 2): Exclude very vague predicates (87% MetaKG coverage)
- Specific Only (Level 3): Only mechanistic relationships (52% MetaKG coverage)
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Set

import requests
import yaml

logger = logging.getLogger(__name__)

# Configuration
DATA_DIR = Path(__file__).parent.parent.parent.parent.parent / "data"
BIOLINK_CACHE_PATH = DATA_DIR / "biolink-model.yaml"
BIOLINK_VERSION_PATH = DATA_DIR / "biolink-model-version.txt"
GITHUB_REPO = "biolink/biolink-model"

# Granularity presets mapping UI labels to depth levels
GRANULARITY_PRESETS = {
    "All Relationships": {
        "min_depth": 1,
        "description": "Include all predicates except root 'related_to'",
    },
    "Standard": {
        "min_depth": 2,
        "description": "Exclude very vague predicates (87% coverage)",
    },
    "Specific Only": {
        "min_depth": 3,
        "description": "Only mechanistic relationships (52% coverage)",
    },
}

# ============================================================================
# Disease → BiologicalProcess Predicate Filtering
# ============================================================================
# Based on analysis in notebooks/disease_bioprocess_exploration.ipynb
# These predicates provide mechanistic/causal relationships vs text mining noise

# Informative predicates to KEEP for Disease → BiologicalProcess queries
DISEASE_BP_INFORMATIVE_PREDICATES = [
    "biolink:affects",
    "biolink:causes",
    "biolink:disrupts",
    "biolink:disease_has_basis_in",
    "biolink:related_to",
    "biolink:correlated_with",
    "biolink:contributes_to",
    "biolink:positively_correlated_with",
    "biolink:negatively_correlated_with",
    "biolink:positively_associated_with",
    "biolink:negatively_associated_with",
]

# Noise predicates to FILTER OUT for Disease → BiologicalProcess queries
# These are mostly text mining artifacts or overly broad relationships
DISEASE_BP_NOISE_PREDICATES = [
    "biolink:occurs_together_in_literature_with",  # 60%+ of results - text mining
    "biolink:manifestation_of",  # Too broad
    "biolink:coexists_with",  # Non-specific
    "biolink:actively_involves",  # Often duplicates
]

# Default intermediate categories for Gene → [Intermediate] → BiologicalProcess queries
# Based on MetaKG exploration in notebooks/intermediate_bioprocess_exploration.ipynb
# All have high coverage (5+ APIs) on both Gene→Intermediate and Intermediate→BP hops
DEFAULT_BP_INTERMEDIATE_CATEGORIES = [
    "biolink:ChemicalEntity",
    "biolink:Protein",
    "biolink:Gene",
    "biolink:Pathway",
    "biolink:MolecularActivity",
    "biolink:CellularComponent",
    "biolink:AnatomicalEntity",
    "biolink:PhenotypicFeature",
]

# Module-level cache to avoid reloading
_predicate_depths_cache: Optional[Dict[str, int]] = None
_predicates_cache: Optional[Dict[str, Dict]] = None


def _normalize_predicate_name(name: str) -> str:
    """Normalize predicate name: spaces to underscores, lowercase.

    Biolink model uses spaces in slot names (e.g., "related to"),
    but the API uses underscores (e.g., "related_to").
    """
    return name.replace(" ", "_").lower().strip()


def get_latest_biolink_version() -> Optional[str]:
    """Fetch latest release tag from GitHub API."""
    try:
        resp = requests.get(
            f"https://api.github.com/repos/{GITHUB_REPO}/releases/latest",
            timeout=10
        )
        if resp.status_code == 200:
            return resp.json().get("tag_name")
        else:
            logger.warning(f"GitHub API returned status {resp.status_code}")
    except requests.RequestException as e:
        logger.warning(f"Could not check GitHub for updates: {e}")
    return None


def get_local_version() -> Optional[str]:
    """Read cached version if exists."""
    if BIOLINK_VERSION_PATH.exists():
        return BIOLINK_VERSION_PATH.read_text().strip()
    return None


def _fetch_biolink_model(version: str) -> str:
    """Fetch biolink-model.yaml for specific version/tag."""
    url = f"https://raw.githubusercontent.com/{GITHUB_REPO}/{version}/biolink-model.yaml"
    logger.info(f"Fetching biolink model from: {url}")
    resp = requests.get(url, timeout=60)
    resp.raise_for_status()
    return resp.text


def load_biolink_model() -> dict:
    """Load biolink model, updating cache if newer version available.

    Returns:
        Parsed biolink model as dict

    Raises:
        RuntimeError: If no local cache and cannot download
    """
    DATA_DIR.mkdir(exist_ok=True)

    local_version = get_local_version()
    latest_version = get_latest_biolink_version()

    need_update = False

    if latest_version:
        if local_version is None:
            logger.info(f"No local cache found. Downloading {latest_version}...")
            need_update = True
        elif latest_version != local_version:
            logger.info(f"Update available: {local_version} -> {latest_version}")
            need_update = True
        else:
            logger.debug(f"Local cache is current: {local_version}")
    else:
        logger.warning("Could not check for updates. Using local cache if available.")

    if need_update and latest_version:
        try:
            yaml_content = _fetch_biolink_model(latest_version)
            BIOLINK_CACHE_PATH.write_text(yaml_content)
            BIOLINK_VERSION_PATH.write_text(latest_version)
            logger.info(f"Successfully cached biolink model version {latest_version}")
        except Exception as e:
            logger.warning(f"Failed to download update: {e}")
            if not BIOLINK_CACHE_PATH.exists():
                raise RuntimeError("No local cache and cannot download biolink model")
            logger.info("Falling back to existing cache.")

    if not BIOLINK_CACHE_PATH.exists():
        raise RuntimeError(
            f"No biolink-model.yaml found at {BIOLINK_CACHE_PATH}. "
            "Check network connection."
        )

    return yaml.safe_load(BIOLINK_CACHE_PATH.read_text())


def _extract_predicates(model: dict) -> Dict[str, Dict]:
    """Extract predicates (slots) and their hierarchy relationships.

    Filters to only include slots that are part of the predicate hierarchy
    (those that eventually trace back to 'related_to' via is_a).
    """
    slots = model.get("slots", {})
    predicates = {}

    # First pass: collect all slots with is_a relationships
    for name, definition in slots.items():
        if definition is None:
            continue

        normalized_name = _normalize_predicate_name(name)
        is_a_raw = definition.get("is_a")
        is_a = _normalize_predicate_name(is_a_raw) if is_a_raw else None

        # Always include 'related_to' as root, plus any slot with is_a
        if normalized_name == "related_to" or is_a:
            predicates[normalized_name] = {
                "is_a": is_a,
                "description": definition.get("description", ""),
                "inverse": _normalize_predicate_name(definition.get("inverse", "")) if definition.get("inverse") else None,
                "symmetric": definition.get("symmetric", False),
                "original_name": name,
            }

    # Second pass: filter to only predicates in the related_to hierarchy
    def traces_to_related_to(name: str, visited: Optional[Set[str]] = None) -> bool:
        """Check if predicate eventually inherits from related_to."""
        if visited is None:
            visited = set()
        if name in visited:
            return False  # Cycle detection
        visited.add(name)

        if name == "related_to":
            return True
        if name not in predicates:
            return False
        parent = predicates[name].get("is_a")
        if parent:
            return traces_to_related_to(parent, visited)
        return False

    # Filter to related_to hierarchy only
    related_to_predicates = {}
    for name, info in predicates.items():
        if name == "related_to" or traces_to_related_to(name):
            related_to_predicates[name] = info

    logger.debug(f"Extracted {len(related_to_predicates)} predicates in related_to hierarchy")
    return related_to_predicates


def _calculate_depths(predicates: Dict[str, Dict]) -> Dict[str, int]:
    """Calculate depth of each predicate from root (related_to)."""
    depths = {}

    def calculate_depth(name: str, visited: Optional[Set[str]] = None) -> int:
        if visited is None:
            visited = set()
        if name in visited:
            return 0  # Cycle - shouldn't happen
        visited.add(name)

        if name in depths:
            return depths[name]

        if name == "related_to":
            depths[name] = 0
            return 0

        if name not in predicates:
            return 0

        parent = predicates[name].get("is_a")
        if parent:
            parent_depth = calculate_depth(parent, visited)
            depths[name] = parent_depth + 1
        else:
            depths[name] = 0

        return depths[name]

    for name in predicates:
        calculate_depth(name)

    return depths


def get_predicate_depths() -> Dict[str, int]:
    """Get mapping of predicate names to their depths in the hierarchy.

    Cached at module level to avoid recomputing.

    Returns:
        Dict mapping predicate name (underscored) to depth (0 = root)
    """
    global _predicate_depths_cache, _predicates_cache

    if _predicate_depths_cache is None:
        model = load_biolink_model()
        _predicates_cache = _extract_predicates(model)
        _predicate_depths_cache = _calculate_depths(_predicates_cache)
        logger.info(f"Loaded {len(_predicate_depths_cache)} predicate depths")

    return _predicate_depths_cache


def filter_predicates_by_granularity(
    predicates: List[str],
    min_depth: int,
    exclude_literature: bool = False,
    exclude_coexpression: bool = False,
    exclude_homology: bool = False,
) -> List[str]:
    """Filter predicate list by granularity level.

    Args:
        predicates: List of predicates (may include biolink: prefix)
        min_depth: Minimum depth required (0 = all, 1 = exclude root, etc.)
        exclude_literature: If True, exclude predicates with 'literature' in name
        exclude_coexpression: If True, exclude 'coexpressed_with' predicate
        exclude_homology: If True, exclude 'homologous_to' and related predicates

    Returns:
        Filtered list of predicates (preserving original format)
    """
    depths = get_predicate_depths()
    result = []

    for pred in predicates:
        # Normalize: remove biolink: prefix if present for lookup
        pred_name = pred.replace("biolink:", "").strip()
        pred_name = _normalize_predicate_name(pred_name)

        # Check depth (predicates not in our tree pass through - don't filter unknown)
        if pred_name in depths:
            if depths[pred_name] < min_depth:
                continue

        # Check literature exclusion
        if exclude_literature and "literature" in pred_name.lower():
            continue

        # Check coexpression exclusion
        if exclude_coexpression and "coexpressed" in pred_name.lower():
            continue

        # Check homology exclusion (covers homologous_to, orthologous_to, paralogous_to, xenologous_to)
        if exclude_homology and ("homolog" in pred_name.lower() or pred_name.lower().endswith("logous_to")):
            continue

        result.append(pred)

    return result


def get_allowed_predicates_for_display(
    min_depth: int,
    exclude_literature: bool = False,
    exclude_coexpression: bool = False,
    exclude_homology: bool = False,
) -> List[str]:
    """Get list of predicates allowed at given granularity level.

    For UI display purposes - shows what predicates will be included.

    Args:
        min_depth: Minimum depth required
        exclude_literature: If True, exclude literature co-occurrence
        exclude_coexpression: If True, exclude 'coexpressed_with' predicate
        exclude_homology: If True, exclude 'homologous_to' and related predicates

    Returns:
        Sorted list of allowed predicate names (human-readable format)
    """
    depths = get_predicate_depths()
    allowed = []

    for pred, depth in depths.items():
        if depth >= min_depth:
            if exclude_literature and "literature" in pred.lower():
                continue
            if exclude_coexpression and "coexpressed" in pred.lower():
                continue
            if exclude_homology and ("homolog" in pred.lower() or pred.lower().endswith("logous_to")):
                continue
            allowed.append(pred)

    return sorted(allowed)


def get_excluded_predicates_for_display(min_depth: int) -> List[str]:
    """Get list of predicates excluded at given granularity level.

    For UI display purposes - shows what predicates will be filtered out.

    Args:
        min_depth: Minimum depth required

    Returns:
        Sorted list of excluded predicate names
    """
    depths = get_predicate_depths()
    excluded = []

    for pred, depth in depths.items():
        if depth < min_depth:
            excluded.append(pred)

    return sorted(excluded)


def get_predicate_info(predicate_name: str) -> Optional[Dict]:
    """Get info about a specific predicate.

    Args:
        predicate_name: Predicate name (with or without biolink: prefix)

    Returns:
        Dict with is_a, description, depth or None if not found
    """
    global _predicates_cache

    # Ensure cache is populated
    get_predicate_depths()

    if _predicates_cache is None:
        return None

    # Normalize name
    name = predicate_name.replace("biolink:", "").strip()
    name = _normalize_predicate_name(name)

    if name in _predicates_cache:
        info = _predicates_cache[name].copy()
        info["depth"] = _predicate_depths_cache.get(name, 0)
        return info

    return None


def get_most_specific_predicate(predicates: List[str]) -> str:
    """Select the most specific (deepest) predicate from a list.

    Used for meta-edge labels when collapsing parallel edges in visualization.
    Higher depth = more specific predicate. Ties broken alphabetically.

    Args:
        predicates: List of predicate strings (with or without biolink: prefix)

    Returns:
        Most specific predicate (original format preserved)
    """
    if not predicates:
        return ""
    if len(predicates) == 1:
        return predicates[0]

    depths = get_predicate_depths()

    # Score each predicate by depth
    scored = []
    for pred in predicates:
        # Normalize for lookup
        pred_name = pred.replace("biolink:", "").strip()
        pred_name = _normalize_predicate_name(pred_name)
        depth = depths.get(pred_name, 0)  # Unknown predicates get depth 0
        scored.append((depth, pred))

    # Sort by depth descending, then alphabetically for ties
    scored.sort(key=lambda x: (-x[0], x[1]))

    return scored[0][1]


def filter_disease_bp_predicates(
    edges: List[Dict],
    use_informative_only: bool = True,
) -> List[Dict]:
    """Filter Disease → BiologicalProcess edges to informative predicates.

    Args:
        edges: List of TRAPI edge dictionaries with 'predicate' key
        use_informative_only: If True, keep only informative predicates.
                              If False, return all edges unfiltered.

    Returns:
        Filtered list of edges
    """
    if not use_informative_only:
        return edges

    # Normalize informative predicates for comparison
    informative_set = set()
    for pred in DISEASE_BP_INFORMATIVE_PREDICATES:
        # Store both with and without prefix for matching
        normalized = pred.replace("biolink:", "").lower()
        informative_set.add(normalized)
        informative_set.add(f"biolink:{normalized}")

    filtered = []
    for edge in edges:
        predicate = edge.get("predicate", "")
        # Normalize for comparison
        pred_normalized = predicate.replace("biolink:", "").lower()

        if pred_normalized in informative_set or predicate.lower() in informative_set:
            filtered.append(edge)

    logger.debug(
        f"Disease-BP predicate filter: {len(edges)} -> {len(filtered)} edges "
        f"({len(edges) - len(filtered)} removed)"
    )
    return filtered
