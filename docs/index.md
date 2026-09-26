# GeneSpeciesTreeVerdict

**Gene tree vs. species tree: diagnose the conflict before deciding how to use the locus.**

GeneSpeciesTreeVerdict is a diagnostic layer for phylogenomics. It compares gene trees with a reference species tree, separates topology discordance from direct copy-number evidence, performs a conservative rooted duplication-loss diagnostic, and returns a transparent locus-handling assessment.

## Core principle

**Discordance is a question, not a diagnosis.** A different gene-tree topology can reflect duplication/loss, HGT, incomplete lineage sorting, recombination, tree-estimation error, model misspecification, sampling, or an incorrect species-tree hypothesis.

## Start here

1. Run the [simulated tutorial](tutorial.md).
2. Read the [core concepts](concepts.md).
3. Follow the [real-data workflow](real-data-workflow.md).
4. Use the [decision framework](decision-framework.md) to interpret `PASS`, `REVIEW`, and `RESOLVE`.
