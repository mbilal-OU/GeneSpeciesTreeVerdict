# Simulated 10-genome tutorial

The bundled tutorial uses a fixed 10-genome reference tree and five gene-family scenarios. It is deliberately small enough to inspect by eye before using real datasets.

## 1. Reference species tree

Ten genomes, G1–G10, define the reference history. Treat this as the current species-tree hypothesis rather than unquestionable truth.

## 2. Concordant single-copy ortholog

`geneA_ortholog_concordant` has one mapped copy per genome and recovers the same unrooted splits as the species tree.

Expected result: `PASS`, topology `CONCORDANT`.

## 3. Single-copy but discordant

`geneD_singlecopy_discordant` still has exactly one copy per species, but G1/G3 and G2/G4 are grouped differently from the reference species tree.

This scenario teaches the most important rule in the repository:

> One copy per species does not guarantee topological concordance, and discordance does not by itself prove paralogy.

Expected result: `REVIEW`, not an automatic paralogy call.

## 4. Recent duplication

`geneB_recent_duplication` gives G3 two copies while the other genomes retain one. Species overlap at an internal gene-tree node provides direct multi-copy evidence.

Expected result: `RESOLVE`.

## 5. Ancient duplication with differential retention

`geneC_ancient_duplication_loss` contains two deep copy lineages for part of the taxon set and one retained lineage in the others.

Expected result: `RESOLVE`.

## 6. Mixed pseudo-ortholog example

`geneP_mixed_pseudoortholog` has exactly one sampled sequence per genome, but the selected tips come from two ancient copy lineages. Copy number alone therefore looks clean while topology is maximally discordant.

Expected result: `REVIEW`. In real work, sequence context, gene-family reconstruction, synteny, annotation, reconciliation and broader sampling may be needed to establish the cause.

## Run it

```bash
gstv tutorial --outdir tutorial_demo

gstv batch \
  --species-tree tutorial_demo/species_tree.nwk \
  --gene-tree-dir tutorial_demo/gene_trees \
  --mapping tutorial_demo/mapping.tsv \
  --outdir tutorial_demo/results
```

Inspect `locus_summary.tsv`, `report.md`, `reference_split_support.tsv` and the per-locus reasons.
