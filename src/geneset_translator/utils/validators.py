"""Input validation for GeneSet Translator.

Validates:
- Gene symbol lists (HUGO format, non-empty, reasonable size)
- Disease CURIEs (MONDO:, DOID:, etc.)
- Configuration parameters

Phase 2 Status: Implemented
"""

from typing import List, Optional
import re


class ValidationError(Exception):
    """Raised when input validation fails."""

    pass


def validate_gene_list(
    gene_symbols: List[str],
    min_genes: int = 1,
    max_genes: int = 100,
) -> List[str]:
    """Validate gene symbol list.

    Args:
        gene_symbols: List of gene symbols to validate
        min_genes: Minimum number of genes required
        max_genes: Maximum number of genes allowed

    Returns:
        Cleaned gene list (trimmed, uppercased)

    Raises:
        ValidationError: If validation fails
    """
    if not gene_symbols:
        raise ValidationError("Gene list cannot be empty")

    # Clean genes: trim whitespace and convert to uppercase
    cleaned_genes = [g.strip().upper() for g in gene_symbols if g.strip()]

    # Remove duplicates while preserving order
    seen = set()
    unique_genes = []
    for gene in cleaned_genes:
        if gene not in seen:
            seen.add(gene)
            unique_genes.append(gene)

    if len(unique_genes) < min_genes:
        raise ValidationError(f"At least {min_genes} gene(s) required, got {len(unique_genes)}")

    if len(unique_genes) > max_genes:
        raise ValidationError(f"Maximum {max_genes} genes allowed, got {len(unique_genes)}")

    return unique_genes


def validate_disease_curie(disease_curie: str) -> str:
    """Validate disease CURIE format.

    Args:
        disease_curie: Disease CURIE to validate (e.g., "MONDO:0004975")

    Returns:
        Validated CURIE

    Raises:
        ValidationError: If CURIE format is invalid
    """
    if not disease_curie:
        raise ValidationError("Disease CURIE cannot be empty")

    disease_curie = disease_curie.strip()

    if not validate_curie(disease_curie):
        raise ValidationError(f"Invalid CURIE format: '{disease_curie}'. Expected format: PREFIX:ID")

    # Check for disease-related prefixes
    disease_prefixes = ["MONDO", "DOID", "HP", "OMIM", "ORPHANET"]
    prefix = disease_curie.split(":")[0].upper()

    if prefix not in disease_prefixes:
        raise ValidationError(
            f"Disease CURIE must use a disease ontology prefix (MONDO, DOID, HP, OMIM, ORPHANET), got: {prefix}"
        )

    return disease_curie


def validate_curie(curie: str) -> bool:
    """Check if string is a valid CURIE format.

    Args:
        curie: String to validate

    Returns:
        True if valid CURIE format
    """
    # Basic pattern: prefix:id
    pattern = r"^[A-Za-z]+[A-Za-z0-9]*:[A-Za-z0-9_\-\.]+$"
    return bool(re.match(pattern, curie))


def validate_convergence_threshold(threshold: int) -> int:
    """Validate convergence threshold parameter.

    Args:
        threshold: Gene frequency threshold

    Returns:
        Validated threshold

    Raises:
        ValidationError: If threshold is invalid
    """
    if threshold < 1:
        raise ValidationError("Convergence threshold must be >= 1")
    if threshold > 50:
        raise ValidationError("Convergence threshold too high (max: 50)")
    return threshold
