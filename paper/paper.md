---
title: 'GeneSpeciesTreeVerdict: transparent diagnosis of gene-tree/species-tree discordance for phylogenomic locus assessment'
tags:
  - phylogenetics
  - phylogenomics
  - gene tree
  - species tree
  - reconciliation
authors:
  - name: Muhammad Bilal
    affiliation: 1
affiliations:
  - name: To be finalized before submission
    index: 1
date: Draft only
bibliography: paper.bib
---

# Summary

GeneSpeciesTreeVerdict is an open-source diagnostic toolkit for comparing gene trees with a reference species tree while separating topology discordance from copy-number and duplication evidence. The software is designed around a common phylogenomic decision problem: whether a locus should be treated as a conventional single-copy marker, reviewed because of ambiguous discordance, or resolved because the family is multi-copy. The package emphasizes transparent evidence and explicitly avoids treating tree disagreement as proof of paralogy.

# Statement of need

A mature version of this section will document the gap between full reconciliation/inference tools and the practical need for an interpretable front-end QC layer for locus-by-locus phylogenomic decisions.

# Functionality

The initial implementation includes taxon coverage, copy-number screening, split/RF comparison, species-overlap duplication screening, rooted duplication-loss LCA diagnostics, batch reference-branch support and reproducible simulated tutorials.

# Validation

This section will be completed after simulation benchmarking, real-data case studies and external-user validation.

# Acknowledgements

To be completed.

# References
