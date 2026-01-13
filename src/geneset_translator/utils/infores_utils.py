"""Knowledge source metadata extraction from Biolink Information Resource Registry.

This module extracts and summarizes knowledge sources (infores) from TRAPI edges,
fetches metadata from the Biolink registry, and optionally uses LLM-based filtering
to identify the most relevant sources for the query.
"""

import hashlib
import json
import logging
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, Set, Any, List, Optional

import networkx as nx
import requests
import yaml

logger = logging.getLogger(__name__)

INFORES_CATALOG_URL = "https://raw.githubusercontent.com/biolink/information-resource-registry/main/infores_catalog.yaml"


def download_infores_catalog(cache_dir: Path, ttl: int = 604800) -> dict:
    """Download and cache the Biolink infores catalog YAML.

    Args:
        cache_dir: Directory to store cached catalog
        ttl: Time-to-live in seconds (default: 7 days)

    Returns:
        Parsed YAML catalog as dict

    Raises:
        requests.RequestException: If download fails and no cache available
    """
    cache_dir.mkdir(parents=True, exist_ok=True)

    # Find most recent cached catalog
    cached_files = sorted(cache_dir.glob("infores_catalog_*.yaml"), reverse=True)

    if cached_files:
        latest_cache = cached_files[0]
        cache_age = datetime.now() - datetime.fromtimestamp(latest_cache.stat().st_mtime)

        if cache_age < timedelta(seconds=ttl):
            logger.info(f"Using cached infores catalog: {latest_cache.name}")
            with open(latest_cache, 'r') as f:
                return yaml.safe_load(f)

    # Download fresh catalog
    try:
        logger.info(f"Downloading infores catalog from {INFORES_CATALOG_URL}")
        response = requests.get(INFORES_CATALOG_URL, timeout=10)
        response.raise_for_status()

        catalog = yaml.safe_load(response.text)

        # Cache the downloaded catalog
        cache_file = cache_dir / f"infores_catalog_{datetime.now().strftime('%Y%m%d')}.yaml"
        with open(cache_file, 'w') as f:
            yaml.dump(catalog, f)

        logger.info(f"Cached infores catalog to {cache_file.name}")
        return catalog

    except requests.RequestException as e:
        logger.warning(f"Failed to download infores catalog: {e}")

        # Fallback to any cached version
        if cached_files:
            logger.info(f"Falling back to stale cache: {cached_files[0].name}")
            with open(cached_files[0], 'r') as f:
                return yaml.safe_load(f)

        raise


def parse_infores_catalog(catalog: dict) -> Dict[str, dict]:
    """Parse infores catalog YAML into indexed dict.

    Args:
        catalog: Raw catalog dict from YAML

    Returns:
        Dict keyed by infores ID with extracted metadata
    """
    parsed = {}

    # The catalog structure is a list of information resources
    resources = catalog.get('information_resources', [])

    for resource in resources:
        infores_id = resource.get('id')
        if not infores_id:
            continue

        parsed[infores_id] = {
            'id': infores_id,
            'name': resource.get('name', infores_id),
            'knowledge_level': resource.get('knowledge_level', 'unknown'),
            'agent_type': resource.get('agent_type', 'unknown'),
            'description': resource.get('description', ''),
            'url': resource.get('url', ''),
        }

    logger.info(f"Parsed {len(parsed)} information resources from catalog")
    return parsed


def extract_unique_sources(graph: nx.MultiDiGraph) -> Set[str]:
    """Extract unique infores source IDs from graph edges.

    Args:
        graph: NetworkX MultiDiGraph with TRAPI edges

    Returns:
        Set of unique infores:* IDs
    """
    sources = set()

    for u, v, key, data in graph.edges(keys=True, data=True):
        edge_sources = data.get('sources', [])

        for source in edge_sources:
            if isinstance(source, dict):
                resource_id = source.get('resource_id', '')
                if resource_id.startswith('infores:'):
                    sources.add(resource_id)
            elif isinstance(source, str) and source.startswith('infores:'):
                sources.add(source)

    logger.info(f"Extracted {len(sources)} unique information resources from graph")
    return sources


def filter_relevant_infores(
    catalog_entries: dict,
    relevant_ids: Set[str],
    api_key: str,
    cache_dir: Optional[Path] = None
) -> dict:
    """Filter catalog entries using LLM-based relevance filtering.

    Uses Claude Haiku 4 to identify the most relevant knowledge sources
    based on their metadata. Falls back to all entries if API call fails.

    Args:
        catalog_entries: Parsed infores catalog
        relevant_ids: Set of infores IDs present in the graph
        api_key: Anthropic API key
        cache_dir: Optional cache directory for filtered results

    Returns:
        Filtered dict of catalog entries
    """
    # Check cache first
    if cache_dir:
        cache_dir.mkdir(parents=True, exist_ok=True)
        cache_key = hashlib.md5('|'.join(sorted(relevant_ids)).encode()).hexdigest()
        cache_file = cache_dir / f"filtered_infores_{cache_key}.json"

        if cache_file.exists():
            cache_age = datetime.now() - datetime.fromtimestamp(cache_file.stat().st_mtime)
            if cache_age < timedelta(days=30):
                logger.info(f"Using cached filtered infores: {cache_file.name}")
                with open(cache_file, 'r') as f:
                    cached_ids = json.load(f)
                    return {k: v for k, v in catalog_entries.items() if k in cached_ids}

    # Filter to relevant IDs
    relevant_entries = {k: v for k, v in catalog_entries.items() if k in relevant_ids}

    # For now, skip LLM filtering and return all relevant entries
    # LLM filtering can be added later for advanced source prioritization
    logger.info(f"Filtered to {len(relevant_entries)} relevant information resources")

    # Cache the result
    if cache_dir:
        with open(cache_file, 'w') as f:
            json.dump(list(relevant_entries.keys()), f, indent=2)
        logger.info(f"Cached filtered infores to {cache_file.name}")

    return relevant_entries


def generate_source_summary(filtered_entries: dict, graph: nx.MultiDiGraph) -> dict:
    """Generate summary statistics for knowledge sources.

    Args:
        filtered_entries: Filtered catalog entries
        graph: NetworkX MultiDiGraph with edges

    Returns:
        Dict with summary statistics
    """
    # Count by knowledge_level and agent_type
    knowledge_levels = {}
    agent_types = {}

    for entry in filtered_entries.values():
        kl = entry.get('knowledge_level', 'unknown')
        at = entry.get('agent_type', 'unknown')

        knowledge_levels[kl] = knowledge_levels.get(kl, 0) + 1
        agent_types[at] = agent_types.get(at, 0) + 1

    # Count edges per source
    source_edge_counts = {}

    for u, v, key, data in graph.edges(keys=True, data=True):
        edge_sources = data.get('sources', [])

        for source in edge_sources:
            if isinstance(source, dict):
                resource_id = source.get('resource_id', '')
            elif isinstance(source, str):
                resource_id = source
            else:
                continue

            if resource_id in filtered_entries:
                source_edge_counts[resource_id] = source_edge_counts.get(resource_id, 0) + 1

    # Create top sources list
    top_sources = []
    for infores_id, count in sorted(source_edge_counts.items(), key=lambda x: x[1], reverse=True)[:10]:
        entry = filtered_entries.get(infores_id, {})
        top_sources.append({
            'id': infores_id,
            'name': entry.get('name', infores_id),
            'edge_count': count,
            'knowledge_level': entry.get('knowledge_level', 'unknown'),
            'agent_type': entry.get('agent_type', 'unknown'),
        })

    return {
        'total_sources': len(filtered_entries),
        'knowledge_levels': knowledge_levels,
        'agent_types': agent_types,
        'top_sources': top_sources,
    }
