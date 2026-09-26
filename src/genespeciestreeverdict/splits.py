from __future__ import annotations

from dataclasses import dataclass

from Bio.Phylo.BaseTree import Clade, Tree

from .models import TopologyMetrics

Split = frozenset[str]


@dataclass(slots=True)
class SplitCollection:
    universe: frozenset[str]
    splits: set[Split]
    supports: dict[Split, float | None]


def _canonical_unrooted(side: set[str], universe: set[str]) -> Split:
    other = universe - side
    if len(side) < len(other):
        return frozenset(side)
    if len(other) < len(side):
        return frozenset(other)
    left = tuple(sorted(side))
    right = tuple(sorted(other))
    return frozenset(side if left <= right else other)


def split_to_text(split: Split, universe: set[str]) -> str:
    left = ",".join(sorted(split))
    right = ",".join(sorted(universe - set(split)))
    return f"{left} | {right}"


def descendant_taxa(clade: Clade, leaf_to_taxon: dict[str, str], restrict_to: set[str]) -> set[str]:
    taxa: set[str] = set()
    for terminal in clade.get_terminals():
        name = str(terminal.name)
        taxon = leaf_to_taxon.get(name)
        if taxon in restrict_to:
            taxa.add(taxon)
    return taxa


def extract_splits(
    tree: Tree,
    leaf_to_taxon: dict[str, str],
    universe: set[str],
    rooted: bool = False,
) -> SplitCollection:
    splits: set[Split] = set()
    supports: dict[Split, float | None] = {}
    root = tree.root

    for clade in tree.find_clades(order="preorder"):
        if clade is root or clade.is_terminal():
            continue
        side = descendant_taxa(clade, leaf_to_taxon, universe)
        if rooted:
            if not (1 < len(side) < len(universe)):
                continue
            key = frozenset(side)
        else:
            if len(side) < 2 or len(universe - side) < 2:
                continue
            key = _canonical_unrooted(side, universe)
        splits.add(key)
        value = float(clade.confidence) if clade.confidence is not None else None
        if key not in supports or (
            value is not None and (supports[key] is None or value > supports[key])
        ):
            supports[key] = value

    return SplitCollection(frozenset(universe), splits, supports)


def compare_topologies(
    species_tree: Tree,
    gene_tree: Tree,
    gene_to_species: dict[str, str],
    rooted: bool = False,
) -> TopologyMetrics:
    species_taxa = {str(node.name) for node in species_tree.get_terminals()}
    mapped_taxa = set(gene_to_species.values())
    shared_taxa = sorted(species_taxa & mapped_taxa)
    common = set(shared_taxa)

    min_taxa = 3 if rooted else 4
    if len(common) < min_taxa:
        return TopologyMetrics(
            comparable=False,
            rooted=rooted,
            shared_taxa=shared_taxa,
            note=f"At least {min_taxa} shared taxa are required for this comparison mode.",
        )

    if len(mapped_taxa) != len(gene_to_species):
        return TopologyMetrics(
            comparable=False,
            rooted=rooted,
            shared_taxa=shared_taxa,
            note="RF/split comparison is disabled for multi-copy gene trees.",
        )

    species_map = {str(node.name): str(node.name) for node in species_tree.get_terminals()}
    species = extract_splits(species_tree, species_map, common, rooted=rooted)
    gene = extract_splits(gene_tree, gene_to_species, common, rooted=rooted)

    shared = species.splits & gene.splits
    species_only = species.splits - gene.splits
    gene_only = gene.splits - species.splits
    rf = len(species_only) + len(gene_only)
    denominator = len(species.splits) + len(gene.splits)
    nrf = rf / denominator if denominator else 0.0

    species_recovery = len(shared) / len(species.splits) if species.splits else 1.0
    gene_agreement = len(shared) / len(gene.splits) if gene.splits else 1.0
    universe = set(shared_taxa)

    species_text = sorted(split_to_text(split, universe) for split in species.splits)
    gene_text = sorted(split_to_text(split, universe) for split in gene.splits)
    species_only_text = sorted(split_to_text(split, universe) for split in species_only)
    gene_only_text = sorted(split_to_text(split, universe) for split in gene_only)
    support_by_text = {
        split_to_text(split, universe): gene.supports.get(split) for split in gene_only
    }

    return TopologyMetrics(
        comparable=True,
        rooted=rooted,
        shared_taxa=shared_taxa,
        rf_distance=rf,
        normalized_rf=nrf,
        species_split_count=len(species.splits),
        gene_split_count=len(gene.splits),
        shared_split_count=len(shared),
        species_split_recovery=species_recovery,
        gene_split_agreement=gene_agreement,
        species_splits=species_text,
        gene_splits=gene_text,
        species_only_splits=species_only_text,
        gene_only_splits=gene_only_text,
        gene_only_split_supports=support_by_text,
    )
