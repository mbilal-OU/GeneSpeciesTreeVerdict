# Core concepts

## Species tree
A species tree represents the evolutionary relationships among sampled species or genomes under a stated inference framework.

## Gene tree
A gene tree represents the history inferred for homologous copies of one gene family. It can differ from the species tree for biological and analytical reasons.

## Orthologs
Orthologous copies trace their relevant relationship through speciation. A conventional concatenated species-tree marker set usually aims to compare the same orthologous lineage across taxa.

## Paralogs
Paralogous copies trace their relevant relationship through gene duplication. A gene family can contain orthologous subgroups and paralogous relationships at the same time.

## Core is not the same as single-copy ortholog
A family can be present in every genome and still have two or more copies in some taxa. Therefore:

```text
core presence + single-copy structure + defensible orthology
                    ↓
      conventional phylogenomic marker candidate
```

## Discordance
Gene-tree/species-tree disagreement does not identify its own cause. The same conflict may be compatible with duplication/loss, ILS, HGT, recombination, error or an incorrect reference tree.

## Reconciliation
Reconciliation maps a gene tree into a species tree under an explicit evolutionary-event model. GeneSpeciesTreeVerdict v0.1 includes a simple rooted duplication-loss LCA diagnostic and labels its output as model-conditional.
