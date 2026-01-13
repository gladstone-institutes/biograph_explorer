"""Human Protein Atlas client for cell type specificity annotations.

This module provides a client for fetching cell type specificity data from
the Human Protein Atlas (HPA) via their JSON API. It also handles the
conversion of NCBI Gene IDs to Ensembl IDs using MyGene.info.

HPA provides single-cell RNA expression data including:
- Cell type specificity category (e.g., "Cell type enhanced", "Low specificity")
- Cell type distribution (e.g., "Detected in many")
- nCPM values per cell type
"""

import json
import logging
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

import networkx as nx
import requests

logger = logging.getLogger(__name__)

# API endpoints
MYGENE_API_URL = "https://mygene.info/v3/query"
HPA_API_URL = "https://www.proteinatlas.org"


@dataclass
class HPAAnnotation:
    """Human Protein Atlas annotation for a gene.

    Attributes:
        ensembl_id: Ensembl gene ID (ENSG...)
        gene_symbol: Gene symbol from HPA
        cell_type_specificity: Specificity category (e.g., "Cell type enhanced")
        cell_type_distribution: Distribution category (e.g., "Detected in many")
        cell_type_ncpm: Dict mapping cell type names to nCPM expression values
        top_cell_types: List of (cell_type, ncpm) tuples sorted by expression (top 5)
        disease_involvement: List of disease associations (e.g., ["Alzheimer disease"])
        protein_class: List of protein functional categories
        tissue_specificity: Tissue specificity category (e.g., "Tissue enhanced")
        tissue_ntpm: Dict mapping tissue names to nTPM expression values
        top_tissues: List of (tissue, ntpm) tuples sorted by expression (top 3)
        immune_cell_specificity: Immune/blood cell specificity (e.g., "Immune cell enhanced")
        immune_cell_distribution: Distribution in immune cells (e.g., "Detected in many")
        immune_cell_ntpm: Dict mapping immune cell types to nTPM expression values
        top_immune_cells: List of (cell_type, ntpm) tuples sorted by expression (top 5)
        immune_lineage_specificity: Blood lineage specificity (e.g., "Group enriched")
        top_immune_lineages: List of (lineage, ntpm) tuples sorted by expression (top 3)
    """
    ensembl_id: str
    gene_symbol: Optional[str] = None
    # Single cell specificity
    cell_type_specificity: Optional[str] = None
    cell_type_distribution: Optional[str] = None
    cell_type_ncpm: Dict[str, float] = field(default_factory=dict)
    top_cell_types: List[Tuple[str, float]] = field(default_factory=list)
    # Disease and protein classification
    disease_involvement: List[str] = field(default_factory=list)
    protein_class: List[str] = field(default_factory=list)
    # Tissue specificity
    tissue_specificity: Optional[str] = None
    tissue_ntpm: Dict[str, float] = field(default_factory=dict)
    top_tissues: List[Tuple[str, float]] = field(default_factory=list)
    # Immune/blood cell specificity
    immune_cell_specificity: Optional[str] = None
    immune_cell_distribution: Optional[str] = None
    immune_cell_ntpm: Dict[str, float] = field(default_factory=dict)
    top_immune_cells: List[Tuple[str, float]] = field(default_factory=list)
    immune_lineage_specificity: Optional[str] = None
    top_immune_lineages: List[Tuple[str, float]] = field(default_factory=list)

    @classmethod
    def from_hpa_response(cls, ensembl_id: str, data: dict) -> "HPAAnnotation":
        """Create HPAAnnotation from HPA API response.

        Args:
            ensembl_id: The Ensembl ID used to query HPA.
            data: Raw JSON response from HPA API.

        Returns:
            HPAAnnotation instance.
        """
        # Single cell data
        cell_type_ncpm = data.get("RNA single cell type specific nCPM", {}) or {}
        top_cell_types = sorted(
            cell_type_ncpm.items(),
            key=lambda x: x[1],
            reverse=True
        )[:5]  # Top 5

        # Tissue data
        tissue_ntpm = data.get("RNA tissue specific nTPM", {}) or {}
        top_tissues = sorted(
            tissue_ntpm.items(),
            key=lambda x: x[1],
            reverse=True
        )[:3]  # Top 3

        # Immune/blood cell data
        immune_cell_ntpm = data.get("RNA blood cell specific nTPM", {}) or {}
        top_immune_cells = sorted(
            immune_cell_ntpm.items(),
            key=lambda x: x[1],
            reverse=True
        )[:5]  # Top 5

        # Immune lineage data
        immune_lineage_ntpm = data.get("RNA blood lineage specific nTPM", {}) or {}
        top_immune_lineages = sorted(
            immune_lineage_ntpm.items(),
            key=lambda x: x[1],
            reverse=True
        )[:3]  # Top 3

        # Disease and protein class (ensure they are lists)
        disease_involvement = data.get("Disease involvement") or []
        if isinstance(disease_involvement, str):
            disease_involvement = [disease_involvement]

        protein_class = data.get("Protein class") or []
        if isinstance(protein_class, str):
            protein_class = [protein_class]

        return cls(
            ensembl_id=ensembl_id,
            gene_symbol=data.get("Gene"),
            # Single cell
            cell_type_specificity=data.get("RNA single cell type specificity"),
            cell_type_distribution=data.get("RNA single cell type distribution"),
            cell_type_ncpm=cell_type_ncpm,
            top_cell_types=top_cell_types,
            # Disease and protein
            disease_involvement=disease_involvement,
            protein_class=protein_class,
            # Tissue
            tissue_specificity=data.get("RNA tissue specificity"),
            tissue_ntpm=tissue_ntpm,
            top_tissues=top_tissues,
            # Immune/blood cells
            immune_cell_specificity=data.get("RNA blood cell specificity"),
            immune_cell_distribution=data.get("RNA blood cell distribution"),
            immune_cell_ntpm=immune_cell_ntpm,
            top_immune_cells=top_immune_cells,
            immune_lineage_specificity=data.get("RNA blood lineage specificity"),
            top_immune_lineages=top_immune_lineages,
        )


@dataclass
class HPATiming:
    """Timing information for HPA batch operations."""
    operation: str
    item_count: int
    duration_seconds: float
    success_count: int = 0
    error: Optional[str] = None


class HPAClient:
    """Client for Human Protein Atlas single-cell RNA data.

    Fetches cell type specificity annotations for genes using:
    1. MyGene.info for NCBI Gene ID -> Ensembl ID conversion
    2. HPA JSON API for cell type specificity data

    Example:
        >>> client = HPAClient()
        >>> ensembl_ids = client.get_ensembl_ids(["348", "920"])
        >>> hpa_data = client.fetch_cell_type_specificity(list(ensembl_ids.values()))
        >>> print(hpa_data["ENSG00000130203"].cell_type_specificity)
        'Cell type enhanced'
    """

    def __init__(
        self,
        cache_dir: Path = Path("data/cache/hpa"),
        timeout: int = 30,
        max_workers: int = 4,
    ):
        """Initialize HPA client.

        Args:
            cache_dir: Directory for caching HPA responses.
            timeout: API timeout in seconds per request.
            max_workers: Maximum parallel requests for HPA API.
        """
        self.cache_dir = cache_dir
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.timeout = timeout
        self.max_workers = max_workers

    def get_ensembl_ids(
        self,
        ncbi_gene_ids: List[str],
        progress_callback: Optional[Callable[[str], None]] = None,
    ) -> Dict[str, List[str]]:
        """Convert NCBI Gene IDs to Ensembl IDs using MyGene.info.

        Args:
            ncbi_gene_ids: List of NCBI Gene IDs (numeric strings, e.g., ["348", "920"]).
            progress_callback: Optional callback for progress updates.

        Returns:
            Dict mapping NCBI Gene ID -> List of Ensembl Gene IDs (ENSG...).
            Returns ALL Ensembl IDs for each gene (some genes have multiple).
            Only includes genes where mapping was found.
        """
        if not ncbi_gene_ids:
            return {}

        logger.info(f"Looking up Ensembl IDs for {len(ncbi_gene_ids)} genes via MyGene.info")
        if progress_callback:
            progress_callback(f"Looking up Ensembl IDs for {len(ncbi_gene_ids)} genes")

        start_time = time.time()
        results: Dict[str, List[str]] = {}

        try:
            # MyGene.info batch query
            response = requests.post(
                MYGENE_API_URL,
                data={
                    "q": ",".join(ncbi_gene_ids),
                    "scopes": "entrezgene",
                    "fields": "ensembl.gene,symbol",
                    "species": "human",
                },
                timeout=self.timeout,
            )
            response.raise_for_status()

            data = response.json()

            for hit in data:
                if "notfound" in hit and hit["notfound"]:
                    continue

                query_id = str(hit.get("query", ""))
                ensembl = hit.get("ensembl", {})

                # Collect ALL Ensembl gene IDs (some genes have multiple)
                ensembl_genes = []
                if isinstance(ensembl, list):
                    for e in ensembl:
                        if isinstance(e, dict) and e.get("gene"):
                            ensembl_genes.append(e["gene"])
                elif isinstance(ensembl, dict) and ensembl.get("gene"):
                    ensembl_genes.append(ensembl["gene"])

                if ensembl_genes and query_id:
                    results[query_id] = ensembl_genes

            duration = time.time() - start_time
            total_ensembl = sum(len(v) for v in results.values())
            logger.info(
                f"MyGene.info lookup complete: {len(results)}/{len(ncbi_gene_ids)} "
                f"genes mapped ({total_ensembl} total Ensembl IDs) in {duration:.2f}s"
            )

        except requests.Timeout:
            logger.warning("MyGene.info lookup timed out")
        except requests.RequestException as e:
            logger.warning(f"MyGene.info lookup failed: {e}")
        except Exception as e:
            logger.warning(f"MyGene.info lookup error: {e}")

        return results

    def get_ensembl_ids_from_uniprot(
        self,
        uniprot_ids: List[str],
        progress_callback: Optional[Callable[[str], None]] = None,
    ) -> Dict[str, str]:
        """Convert UniProtKB IDs to Ensembl IDs using MyGene.info.

        Args:
            uniprot_ids: List of UniProtKB accession IDs (e.g., ["P05067", "P02489"]).
            progress_callback: Optional callback for progress updates.

        Returns:
            Dict mapping UniProtKB ID -> Ensembl Gene ID (ENSG...).
            Only includes proteins where mapping was found.
            For proteins with multiple gene mappings, uses the first (highest score).
        """
        if not uniprot_ids:
            return {}

        logger.info(f"Looking up Ensembl IDs for {len(uniprot_ids)} proteins via MyGene.info")
        if progress_callback:
            progress_callback(f"Looking up Ensembl IDs for {len(uniprot_ids)} proteins")

        start_time = time.time()
        results = {}

        try:
            # MyGene.info batch query with uniprot scope
            response = requests.post(
                MYGENE_API_URL,
                data={
                    "q": ",".join(uniprot_ids),
                    "scopes": "uniprot",
                    "fields": "ensembl.gene,symbol,entrezgene",
                    "species": "human",
                },
                timeout=self.timeout,
            )
            response.raise_for_status()

            data = response.json()

            # Track which queries we've already processed (some proteins map to multiple genes)
            processed_queries = set()

            for hit in data:
                if "notfound" in hit and hit["notfound"]:
                    continue

                query_id = str(hit.get("query", ""))

                # Skip if we already have a result for this query (use first/highest score)
                if query_id in processed_queries:
                    continue

                ensembl = hit.get("ensembl", {})

                # Handle case where ensembl is a list (multiple transcripts)
                if isinstance(ensembl, list):
                    ensembl = ensembl[0] if ensembl else {}

                ensembl_gene = ensembl.get("gene") if isinstance(ensembl, dict) else None

                if ensembl_gene and query_id:
                    results[query_id] = ensembl_gene
                    processed_queries.add(query_id)

            duration = time.time() - start_time
            logger.info(
                f"MyGene.info UniProtKB lookup complete: {len(results)}/{len(uniprot_ids)} "
                f"proteins mapped in {duration:.2f}s"
            )

        except requests.Timeout:
            logger.warning("MyGene.info UniProtKB lookup timed out")
        except requests.RequestException as e:
            logger.warning(f"MyGene.info UniProtKB lookup failed: {e}")
        except Exception as e:
            logger.warning(f"MyGene.info UniProtKB lookup error: {e}")

        return results

    def fetch_cell_type_specificity(
        self,
        ensembl_ids: List[str],
        progress_callback: Optional[Callable[[str], None]] = None,
    ) -> Dict[str, HPAAnnotation]:
        """Fetch cell type specificity data from HPA for multiple genes.

        Args:
            ensembl_ids: List of Ensembl Gene IDs (e.g., ["ENSG00000130203"]).
            progress_callback: Optional callback for progress updates.

        Returns:
            Dict mapping Ensembl ID -> HPAAnnotation.
            Only includes genes where HPA data was found.
        """
        if not ensembl_ids:
            return {}

        logger.info(f"Fetching HPA data for {len(ensembl_ids)} genes")
        if progress_callback:
            progress_callback(f"Fetching HPA cell type data for {len(ensembl_ids)} genes")

        results = {}
        completed = 0
        start_time = time.time()

        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            future_to_id = {
                executor.submit(self._fetch_single_gene, ensembl_id): ensembl_id
                for ensembl_id in ensembl_ids
            }

            for future in as_completed(future_to_id):
                ensembl_id = future_to_id[future]
                completed += 1

                try:
                    annotation = future.result()
                    if annotation:
                        results[ensembl_id] = annotation
                except Exception as e:
                    logger.debug(f"Failed to fetch HPA data for {ensembl_id}: {e}")

                if progress_callback and completed % 10 == 0:
                    progress_callback(
                        f"Fetched HPA data: {completed}/{len(ensembl_ids)} "
                        f"({len(results)} successful)"
                    )

        duration = time.time() - start_time
        logger.info(
            f"HPA fetch complete: {len(results)}/{len(ensembl_ids)} genes "
            f"in {duration:.2f}s"
        )

        return results

    def _fetch_single_gene(self, ensembl_id: str) -> Optional[HPAAnnotation]:
        """Fetch HPA data for a single gene.

        Args:
            ensembl_id: Ensembl Gene ID.

        Returns:
            HPAAnnotation if found, None otherwise.
        """
        try:
            url = f"{HPA_API_URL}/{ensembl_id}.json"
            response = requests.get(url, timeout=self.timeout)

            if response.status_code == 404:
                return None

            response.raise_for_status()
            data = response.json()

            return HPAAnnotation.from_hpa_response(ensembl_id, data)

        except requests.RequestException as e:
            logger.debug(f"HPA request failed for {ensembl_id}: {e}")
            return None
        except json.JSONDecodeError:
            logger.debug(f"Invalid JSON from HPA for {ensembl_id}")
            return None

    def _apply_hpa_annotations(
        self,
        graph: nx.DiGraph,
        hpa_data: Dict[str, HPAAnnotation],
        ensembl_to_curie: Dict[str, str],
        stats: Dict[str, Any],
    ) -> int:
        """Apply HPA annotations to graph nodes.

        Args:
            graph: NetworkX DiGraph to annotate.
            hpa_data: Dict of Ensembl ID -> HPAAnnotation.
            ensembl_to_curie: Dict mapping Ensembl ID -> node CURIE.
            stats: Dict to accumulate statistics (modified in place).

        Returns:
            Number of nodes annotated.
        """
        annotated_count = 0

        for ensembl_id, annotation in hpa_data.items():
            curie = ensembl_to_curie.get(ensembl_id)
            if not curie or curie not in graph.nodes:
                continue

            # Initialize annotation_features if needed
            if "annotation_features" not in graph.nodes[curie]:
                graph.nodes[curie]["annotation_features"] = {}

            features = graph.nodes[curie]["annotation_features"]

            # Core identifier
            features["hpa_ensembl_id"] = annotation.ensembl_id

            # Single cell specificity
            if annotation.cell_type_specificity:
                features["hpa_cell_type_specificity"] = annotation.cell_type_specificity
                stats["all_specificity_categories"].add(annotation.cell_type_specificity)
            if annotation.cell_type_distribution:
                features["hpa_cell_type_distribution"] = annotation.cell_type_distribution

            # Top 5 cell types with nCPM values
            if annotation.top_cell_types:
                features["hpa_top_cell_types"] = annotation.top_cell_types
                stats["all_cell_types"].update(ct for ct, _ in annotation.top_cell_types)

            # Disease involvement
            if annotation.disease_involvement:
                features["hpa_disease_involvement"] = annotation.disease_involvement
                stats["all_diseases"].update(annotation.disease_involvement)

            # Protein classification
            if annotation.protein_class:
                features["hpa_protein_class"] = annotation.protein_class
                stats["all_protein_classes"].update(annotation.protein_class)

            # Tissue specificity
            if annotation.tissue_specificity:
                features["hpa_tissue_specificity"] = annotation.tissue_specificity
                stats["all_tissue_specificities"].add(annotation.tissue_specificity)

            # Top 3 tissues with nTPM values
            if annotation.top_tissues:
                features["hpa_top_tissues"] = annotation.top_tissues

            # Immune/blood cell specificity
            if annotation.immune_cell_specificity:
                features["hpa_immune_cell_specificity"] = annotation.immune_cell_specificity
                stats["all_immune_specificities"].add(annotation.immune_cell_specificity)
            if annotation.immune_cell_distribution:
                features["hpa_immune_cell_distribution"] = annotation.immune_cell_distribution

            # Top 5 immune cells with nTPM values
            if annotation.top_immune_cells:
                features["hpa_top_immune_cells"] = annotation.top_immune_cells
                stats["all_immune_cells"].update(ct for ct, _ in annotation.top_immune_cells)

            # Immune lineage specificity
            if annotation.immune_lineage_specificity:
                features["hpa_immune_lineage_specificity"] = annotation.immune_lineage_specificity

            # Top 3 immune lineages with nTPM values
            if annotation.top_immune_lineages:
                features["hpa_top_immune_lineages"] = annotation.top_immune_lineages

            # Store full data separately (for detailed views, not serialized to Cytoscape)
            graph.nodes[curie]["hpa_cell_type_ncpm"] = annotation.cell_type_ncpm
            graph.nodes[curie]["hpa_tissue_ntpm"] = annotation.tissue_ntpm
            graph.nodes[curie]["hpa_immune_cell_ntpm"] = annotation.immune_cell_ntpm
            graph.nodes[curie]["hpa_annotation"] = annotation

            annotated_count += 1

        return annotated_count

    def annotate_graph(
        self,
        graph: nx.DiGraph,
        progress_callback: Optional[Callable[[str], None]] = None,
    ) -> Tuple[nx.DiGraph, Dict[str, Any]]:
        """Add HPA cell type specificity annotations to gene and protein nodes.

        Annotates all gene nodes (NCBIGene: or HGNC: prefix) and protein nodes
        (UniProtKB: prefix) with HPA data. Protein nodes are mapped to their
        encoding genes via MyGene.info.

        Adds the following node attributes:
        - annotation_features['hpa_*']: Flattened features for filtering
        - hpa_cell_type_ncpm: Full dict of cell type -> nCPM values
        - hpa_tissue_ntpm: Full dict of tissue -> nTPM values

        Args:
            graph: NetworkX DiGraph with gene/protein nodes.
            progress_callback: Optional callback for progress updates.

        Returns:
            Tuple of (annotated graph, metadata dict).
        """
        # Find all gene and protein nodes
        gene_nodes = []
        protein_nodes = []
        for node in graph.nodes():
            prefix = node.split(":")[0] if ":" in node else ""
            if prefix in ["NCBIGene", "HGNC"]:
                gene_nodes.append(node)
            elif prefix == "UniProtKB":
                protein_nodes.append(node)

        total_nodes = len(gene_nodes) + len(protein_nodes)
        logger.info(
            f"Found {len(gene_nodes)} gene nodes and {len(protein_nodes)} protein nodes "
            f"for HPA annotation"
        )

        if total_nodes == 0:
            return graph, {
                "hpa_annotated_count": 0,
                "hpa_total_genes": 0,
                "hpa_total_proteins": 0,
            }

        # Initialize statistics tracking
        stats = {
            "all_specificity_categories": set(),
            "all_cell_types": set(),
            "all_diseases": set(),
            "all_protein_classes": set(),
            "all_tissue_specificities": set(),
            "all_immune_specificities": set(),
            "all_immune_cells": set(),
        }

        # Collect all Ensembl IDs and build mappings
        # For genes with multiple Ensembl IDs, we'll try all of them
        all_ensembl_ids = set()

        # Process gene nodes (NCBIGene -> Ensembl)
        # ncbi_to_ensembl now returns Dict[str, List[str]] (multiple Ensembl IDs per gene)
        ncbi_ids = []
        curie_to_ncbi = {}
        for node in gene_nodes:
            if node.startswith("NCBIGene:"):
                ncbi_id = node.replace("NCBIGene:", "")
                ncbi_ids.append(ncbi_id)
                curie_to_ncbi[node] = ncbi_id

        ncbi_to_ensembl: Dict[str, List[str]] = {}
        if ncbi_ids:
            if progress_callback:
                progress_callback(f"Converting {len(ncbi_ids)} gene IDs to Ensembl format")
            ncbi_to_ensembl = self.get_ensembl_ids(ncbi_ids, progress_callback)

            # Collect ALL Ensembl IDs for all genes
            for ncbi_id, ensembl_list in ncbi_to_ensembl.items():
                all_ensembl_ids.update(ensembl_list)

        # Process protein nodes (UniProtKB -> Ensembl)
        uniprot_ids = []
        curie_to_uniprot = {}
        for node in protein_nodes:
            if node.startswith("UniProtKB:"):
                uniprot_id = node.replace("UniProtKB:", "")
                uniprot_ids.append(uniprot_id)
                curie_to_uniprot[node] = uniprot_id

        uniprot_to_ensembl = {}
        if uniprot_ids:
            if progress_callback:
                progress_callback(f"Converting {len(uniprot_ids)} protein IDs to Ensembl format")
            uniprot_to_ensembl = self.get_ensembl_ids_from_uniprot(uniprot_ids, progress_callback)

            # Build protein mapping (proteins still use single Ensembl ID)
            for uniprot_id, ensembl_id in uniprot_to_ensembl.items():
                all_ensembl_ids.add(ensembl_id)

        total_mapped = len(ncbi_to_ensembl) + len(uniprot_to_ensembl)
        if not all_ensembl_ids:
            logger.warning("No Ensembl IDs found, skipping HPA annotation")
            return graph, {
                "hpa_annotated_count": 0,
                "hpa_total_genes": len(gene_nodes),
                "hpa_total_proteins": len(protein_nodes),
                "hpa_ensembl_mapped": 0,
            }

        # Fetch HPA data for all unique Ensembl IDs
        # This will try all Ensembl IDs and return data for the ones HPA has
        hpa_data = self.fetch_cell_type_specificity(list(all_ensembl_ids), progress_callback)

        # Build mapping from Ensembl ID to CURIE, prioritizing Ensembl IDs that have HPA data
        # For genes with multiple Ensembl IDs, use the first one that has HPA data
        gene_ensembl_to_curie = {}
        for curie, ncbi_id in curie_to_ncbi.items():
            if ncbi_id in ncbi_to_ensembl:
                ensembl_list = ncbi_to_ensembl[ncbi_id]
                # Find the first Ensembl ID that has HPA data
                for ensembl_id in ensembl_list:
                    if ensembl_id in hpa_data:
                        gene_ensembl_to_curie[ensembl_id] = curie
                        break
                else:
                    # None had HPA data, use the first one anyway (for consistency)
                    if ensembl_list:
                        gene_ensembl_to_curie[ensembl_list[0]] = curie

        # Apply annotations to gene nodes
        gene_annotated = 0
        if gene_nodes:
            gene_annotated = self._apply_hpa_annotations(
                graph, hpa_data, gene_ensembl_to_curie, stats
            )

        # Apply annotations to protein nodes
        protein_annotated = 0
        if protein_nodes:
            protein_ensembl_to_curie = {}
            for curie, uniprot_id in curie_to_uniprot.items():
                if uniprot_id in uniprot_to_ensembl:
                    ensembl_id = uniprot_to_ensembl[uniprot_id]
                    protein_ensembl_to_curie[ensembl_id] = curie
            protein_annotated = self._apply_hpa_annotations(
                graph, hpa_data, protein_ensembl_to_curie, stats
            )

        total_annotated = gene_annotated + protein_annotated

        metadata = {
            "hpa_annotated_count": total_annotated,
            "hpa_total_genes": len(gene_nodes),
            "hpa_total_proteins": len(protein_nodes),
            "hpa_genes_annotated": gene_annotated,
            "hpa_proteins_annotated": protein_annotated,
            "hpa_ensembl_mapped": total_mapped,
            "hpa_specificity_categories": sorted(stats["all_specificity_categories"]),
            "hpa_tissue_specificities": sorted(stats["all_tissue_specificities"]),
            "hpa_cell_types_found": len(stats["all_cell_types"]),
            "hpa_diseases_found": sorted(stats["all_diseases"]),
            "hpa_protein_classes_found": sorted(stats["all_protein_classes"]),
            "hpa_immune_specificities": sorted(stats["all_immune_specificities"]),
            "hpa_immune_cells_found": len(stats["all_immune_cells"]),
            "timestamp": datetime.now().isoformat(),
        }

        logger.info(
            f"HPA annotation complete: {gene_annotated}/{len(gene_nodes)} genes, "
            f"{protein_annotated}/{len(protein_nodes)} proteins annotated "
            f"({len(stats['all_specificity_categories'])} specificity categories)"
        )

        if progress_callback:
            progress_callback(
                f"HPA annotation complete: {gene_annotated} genes, {protein_annotated} proteins"
            )

        return graph, metadata

    def save_cache(
        self,
        hpa_data: Dict[str, HPAAnnotation],
        filename: Optional[str] = None,
    ) -> Path:
        """Save HPA annotations to cache file.

        Args:
            hpa_data: Dict of Ensembl ID -> HPAAnnotation.
            filename: Optional filename (auto-generated if not provided).

        Returns:
            Path to saved cache file.
        """
        if filename is None:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            filename = f"hpa_annotations_{timestamp}.json"

        cache_path = self.cache_dir / filename

        # Convert to serializable format
        cache_data = {}
        for ensembl_id, annotation in hpa_data.items():
            cache_data[ensembl_id] = {
                "ensembl_id": annotation.ensembl_id,
                "gene_symbol": annotation.gene_symbol,
                "cell_type_specificity": annotation.cell_type_specificity,
                "cell_type_distribution": annotation.cell_type_distribution,
                "cell_type_ncpm": annotation.cell_type_ncpm,
                "top_cell_types": annotation.top_cell_types,
                "disease_involvement": annotation.disease_involvement,
                "protein_class": annotation.protein_class,
                "tissue_specificity": annotation.tissue_specificity,
                "tissue_ntpm": annotation.tissue_ntpm,
                "top_tissues": annotation.top_tissues,
                "immune_cell_specificity": annotation.immune_cell_specificity,
                "immune_cell_distribution": annotation.immune_cell_distribution,
                "immune_cell_ntpm": annotation.immune_cell_ntpm,
                "top_immune_cells": annotation.top_immune_cells,
                "immune_lineage_specificity": annotation.immune_lineage_specificity,
                "top_immune_lineages": annotation.top_immune_lineages,
            }

        with open(cache_path, "w") as f:
            json.dump(cache_data, f, indent=2)

        logger.info(f"Saved HPA cache to {cache_path}")
        return cache_path

    def load_cache(self, filename: str) -> Dict[str, HPAAnnotation]:
        """Load HPA annotations from cache file.

        Args:
            filename: Cache filename to load.

        Returns:
            Dict of Ensembl ID -> HPAAnnotation.
        """
        cache_path = self.cache_dir / filename

        with open(cache_path) as f:
            cache_data = json.load(f)

        annotations = {}
        for ensembl_id, data in cache_data.items():
            annotations[ensembl_id] = HPAAnnotation(
                ensembl_id=data["ensembl_id"],
                gene_symbol=data.get("gene_symbol"),
                cell_type_specificity=data.get("cell_type_specificity"),
                cell_type_distribution=data.get("cell_type_distribution"),
                cell_type_ncpm=data.get("cell_type_ncpm", {}),
                top_cell_types=[tuple(x) for x in data.get("top_cell_types", [])],
                disease_involvement=data.get("disease_involvement", []),
                protein_class=data.get("protein_class", []),
                tissue_specificity=data.get("tissue_specificity"),
                tissue_ntpm=data.get("tissue_ntpm", {}),
                top_tissues=[tuple(x) for x in data.get("top_tissues", [])],
                immune_cell_specificity=data.get("immune_cell_specificity"),
                immune_cell_distribution=data.get("immune_cell_distribution"),
                immune_cell_ntpm=data.get("immune_cell_ntpm", {}),
                top_immune_cells=[tuple(x) for x in data.get("top_immune_cells", [])],
                immune_lineage_specificity=data.get("immune_lineage_specificity"),
                top_immune_lineages=[tuple(x) for x in data.get("top_immune_lineages", [])],
            )

        logger.info(f"Loaded {len(annotations)} HPA annotations from {cache_path}")
        return annotations
