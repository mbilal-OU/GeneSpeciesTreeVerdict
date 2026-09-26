from __future__ import annotations

from itertools import combinations
from typing import Any

from Bio.Phylo.BaseTree import Clade, Tree

from .models import ReconciliationSummary


def _node_ids(tree: Tree) -> dict[int, str]:
    labels: dict[int, str] = {}
    for index, node in enumerate(tree.find_clades(order="preorder"), start=1):
        labels[id(node)] = f"N{index}"
    return labels


def _descendant_species(node: Clade, mapping: dict[str, str]) -> set[str]:
    return {mapping[str(terminal.name)] for terminal in node.get_terminals()}


def species_overlap_events(gene_tree: Tree, mapping: dict[str, str]) -> list[dict[str, Any]]:
    """Detect candidate duplications by overlap between child species sets.

    This is a conservative, topology-local screen: a node is marked when two child
    lineages contain at least one of the same species. It is strong evidence that the
    family is not single-copy at that point, but it is not a full evolutionary model.
    """

    ids = _node_ids(gene_tree)
    events: list[dict[str, Any]] = []
    for node in gene_tree.find_clades(order="preorder"):
        if node.is_terminal() or len(node.clades) < 2:
            continue
        child_sets = [_descendant_species(child, mapping) for child in node.clades]
        overlaps: set[str] = set()
        for left, right in combinations(child_sets, 2):
            overlaps.update(left & right)
        if overlaps:
            events.append(
                {
                    "node": ids[id(node)],
                    "overlapping_species": sorted(overlaps),
                    "descendant_species": sorted(_descendant_species(node, mapping)),
                    "descendant_genes": sorted(str(x.name) for x in node.get_terminals()),
                    "support": float(node.confidence) if node.confidence is not None else None,
                }
            )
    return events


def _species_descendant_sets(species_tree: Tree) -> dict[int, frozenset[str]]:
    return {
        id(node): frozenset(str(terminal.name) for terminal in node.get_terminals())
        for node in species_tree.find_clades(order="preorder")
    }


def _species_parent_map(species_tree: Tree) -> dict[int, Clade | None]:
    parents: dict[int, Clade | None] = {id(species_tree.root): None}
    for parent in species_tree.find_clades(order="preorder"):
        for child in parent.clades:
            parents[id(child)] = parent
    return parents


def _smallest_species_clade(
    species_tree: Tree,
    target: set[str],
    descendants: dict[int, frozenset[str]],
) -> Clade:
    candidates: list[tuple[int, Clade]] = []
    for node in species_tree.find_clades(order="preorder"):
        node_set = descendants[id(node)]
        if target.issubset(node_set):
            candidates.append((len(node_set), node))
    if not candidates:
        raise ValueError(f"No species-tree LCA contains: {sorted(target)}")
    candidates.sort(key=lambda item: item[0])
    return candidates[0][1]


def _distance_down(parent: Clade, child: Clade, species_parents: dict[int, Clade | None]) -> int:
    if parent is child:
        return 0
    distance = 0
    current: Clade | None = child
    while current is not None and current is not parent:
        current = species_parents[id(current)]
        distance += 1
    if current is None:
        raise ValueError("Gene-node mapping is not nested in the rooted species tree.")
    return distance


def lca_duplication_loss_reconciliation(
    species_tree: Tree,
    gene_tree: Tree,
    mapping: dict[str, str],
) -> ReconciliationSummary:
    """Perform rooted LCA mapping under a duplication-loss (DL) model.

    Important: the inferred duplication/loss counts are conditional on this DL model.
    Discordance caused by ILS, HGT, gene-tree error, or an incorrect reference species
    tree can also be represented as duplications/losses. Therefore these counts are
    evidence for review, not proof of the biological cause.
    """

    species_names = {str(node.name) for node in species_tree.get_terminals()}
    unknown = sorted(set(mapping.values()) - species_names)
    if unknown:
        return ReconciliationSummary(
            attempted=False,
            note="Reconciliation skipped because mapped species are absent from the species tree: "
            + ", ".join(unknown),
        )

    species_desc = _species_descendant_sets(species_tree)
    species_parents = _species_parent_map(species_tree)
    node_ids = _node_ids(gene_tree)

    mapping_by_gene_node: dict[int, Clade] = {}
    for gene_node in gene_tree.find_clades(order="postorder"):
        target = _descendant_species(gene_node, mapping)
        mapping_by_gene_node[id(gene_node)] = _smallest_species_clade(
            species_tree, target, species_desc
        )

    duplication_nodes: set[int] = set()
    events: list[dict[str, Any]] = []
    speciations = 0

    for gene_node in gene_tree.find_clades(order="preorder"):
        if gene_node.is_terminal():
            continue
        parent_map = mapping_by_gene_node[id(gene_node)]
        child_maps = [mapping_by_gene_node[id(child)] for child in gene_node.clades]
        is_duplication = any(child_map is parent_map for child_map in child_maps)
        if is_duplication:
            duplication_nodes.add(id(gene_node))
            events.append(
                {
                    "node": node_ids[id(gene_node)],
                    "mapped_species_clade": sorted(species_desc[id(parent_map)]),
                    "descendant_species": sorted(_descendant_species(gene_node, mapping)),
                    "descendant_genes": sorted(
                        str(terminal.name) for terminal in gene_node.get_terminals()
                    ),
                    "support": (
                        float(gene_node.confidence)
                        if gene_node.confidence is not None
                        else None
                    ),
                }
            )
        else:
            speciations += 1

    inferred_losses = 0
    for gene_parent in gene_tree.find_clades(order="preorder"):
        if gene_parent.is_terminal():
            continue
        parent_map = mapping_by_gene_node[id(gene_parent)]
        parent_is_duplication = id(gene_parent) in duplication_nodes
        for gene_child in gene_parent.clades:
            child_map = mapping_by_gene_node[id(gene_child)]
            distance = _distance_down(parent_map, child_map, species_parents)
            inferred_losses += distance if parent_is_duplication else max(distance - 1, 0)

    return ReconciliationSummary(
        attempted=True,
        lca_duplications=len(duplication_nodes),
        inferred_losses=inferred_losses,
        speciations=speciations,
        duplication_events=events,
        note=(
            "Counts are conditional on a rooted duplication-loss LCA model; they do not by "
            "themselves distinguish paralogy from ILS, HGT, gene-tree error, or species-tree error."
        ),
    )
