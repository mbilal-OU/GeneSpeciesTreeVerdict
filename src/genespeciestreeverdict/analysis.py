from __future__ import annotations

from collections import Counter
from pathlib import Path

from .decision import make_verdict
from .io import (
    MappingRecord,
    family_name_from_path,
    mapping_for_family,
    read_mapping_table,
    read_tree,
    terminal_names,
)
from .models import LocusResult
from .reconcile import lca_duplication_loss_reconciliation, species_overlap_events
from .splits import compare_topologies


def analyze_locus(
    species_tree_path: str | Path,
    gene_tree_path: str | Path,
    *,
    mapping_path: str | Path | None = None,
    mapping_records: list[MappingRecord] | None = None,
    family: str | None = None,
    delimiter: str | None = None,
    species_field: int = 0,
    min_coverage: float = 0.95,
    rooted_comparison: bool = False,
    reconcile: bool = True,
    discordance_review_threshold: float | None = None,
) -> LocusResult:
    if not 0 < min_coverage <= 1:
        raise ValueError("min_coverage must be in the interval (0, 1].")
    if discordance_review_threshold is not None and not 0 <= discordance_review_threshold <= 1:
        raise ValueError("discordance_review_threshold must be between 0 and 1.")

    species_tree = read_tree(species_tree_path)
    gene_tree = read_tree(gene_tree_path)
    family = family or family_name_from_path(gene_tree_path)

    if mapping_records is None and mapping_path is not None:
        mapping_records = read_mapping_table(mapping_path)

    mapping = mapping_for_family(
        gene_tree,
        species_tree,
        records=mapping_records,
        family=family,
        delimiter=delimiter,
        species_field=species_field,
    )

    counts = Counter(mapping.values())
    species_total = len(terminal_names(species_tree))
    unique_species = len(counts)
    coverage = unique_species / species_total if species_total else 0.0
    max_copies = max(counts.values(), default=0)
    multicopy_species = sorted(species for species, count in counts.items() if count > 1)

    overlap_events = species_overlap_events(gene_tree, mapping)
    topology = compare_topologies(
        species_tree,
        gene_tree,
        mapping,
        rooted=rooted_comparison,
    )

    if reconcile:
        reconciliation = lca_duplication_loss_reconciliation(species_tree, gene_tree, mapping)
    else:
        from .models import ReconciliationSummary

        reconciliation = ReconciliationSummary(
            attempted=False, note="Reconciliation disabled by user."
        )

    verdict = make_verdict(
        coverage=coverage,
        min_coverage=min_coverage,
        max_copies=max_copies,
        overlap_duplications=len(overlap_events),
        reconciliation=reconciliation,
        topology=topology,
        discordance_review_threshold=discordance_review_threshold,
    )

    warnings: list[str] = []
    if reconcile:
        warnings.append(
            "DL reconciliation assumes the supplied species-tree root and gene-tree root are meaningful."
        )
    if max_copies == 1 and topology.comparable and topology.rf_distance:
        warnings.append(
            "Single-copy discordance can arise from ILS, HGT/recombination, tree error, model error, "
            "or an incorrect reference species tree; do not label it paralogy from RF distance alone."
        )

    return LocusResult(
        family=family,
        species_total=species_total,
        gene_leaves=len(terminal_names(gene_tree)),
        mapped_genes=len(mapping),
        unique_species=unique_species,
        coverage=coverage,
        copy_counts=dict(sorted(counts.items())),
        max_copies=max_copies,
        multicopy_species=multicopy_species,
        species_overlap_duplications=len(overlap_events),
        species_overlap_events=overlap_events,
        topology=topology,
        reconciliation=reconciliation,
        verdict=verdict,
        warnings=warnings,
    )
