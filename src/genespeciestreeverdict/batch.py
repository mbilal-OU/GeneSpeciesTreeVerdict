from __future__ import annotations

from pathlib import Path

from .analysis import analyze_locus
from .io import discover_gene_trees, read_mapping_table
from .models import LocusResult
from .reporting import write_batch_outputs


def analyze_directory(
    species_tree_path: str | Path,
    gene_tree_dir: str | Path,
    *,
    mapping_path: str | Path | None = None,
    delimiter: str | None = None,
    species_field: int = 0,
    min_coverage: float = 0.95,
    rooted_comparison: bool = False,
    reconcile: bool = True,
    discordance_review_threshold: float | None = None,
    recursive: bool = False,
    outdir: str | Path = "gstv_batch",
) -> tuple[list[LocusResult], list[dict[str, str]]]:
    files = discover_gene_trees(gene_tree_dir, recursive=recursive)
    if not files:
        raise ValueError(f"No Newick-like gene-tree files found in {gene_tree_dir}")
    records = read_mapping_table(mapping_path) if mapping_path else None

    results: list[LocusResult] = []
    errors: list[dict[str, str]] = []
    for path in files:
        try:
            result = analyze_locus(
                species_tree_path,
                path,
                mapping_records=records,
                delimiter=delimiter,
                species_field=species_field,
                min_coverage=min_coverage,
                rooted_comparison=rooted_comparison,
                reconcile=reconcile,
                discordance_review_threshold=discordance_review_threshold,
            )
            results.append(result)
        except Exception as exc:
            errors.append({"family": path.stem, "path": str(path), "error": str(exc)})

    write_batch_outputs(results, errors, outdir)
    return results, errors
