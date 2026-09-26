# Simulated 10-genome tutorial

This deterministic tree-level tutorial is designed to teach what the software measures before you analyze real data.

- `geneA_ortholog_concordant`: single-copy ortholog; topology matches the reference.
- `geneD_singlecopy_discordant`: one copy/species, but the gene tree conflicts with the species tree. Discordance alone is **not** paralogy.
- `geneB_recent_duplication`: a recent extra copy in G3.
- `geneC_ancient_duplication_loss`: an ancient duplicated family with differential copy retention.
- `geneP_mixed_pseudoortholog`: one sequence/species but inconsistent ancient copies; illustrates why one tip/species is not sufficient proof of orthology.

These are topology-level simulations, not sequence simulations.
