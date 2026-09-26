# Changelog

All notable changes to this project will be documented here.

## [Unreleased]

### Planned
- quantitative simulation benchmark across broader evolutionary regimes;
- external reconciliation adapters;
- real-data adapters for common pangenome/orthology outputs;
- stronger uncertainty handling;
- stable packaging and archival release after validation.

## [0.2.0] - 2026-09-25

### Added
- evidence-integration layer separating observation, interpretation and limitation;
- machine-readable `evidence_matrix.tsv` for single-locus and batch analyses;
- explicit supported vs. unsupported conclusions in every verdict;
- machine-readable `primary_concern` categories;
- question-driven next-analysis routing with example specialist tool classes;
- `next_analyses.tsv` and batch `recommendations.tsv` outputs;
- `gstv explain` knowledge topics for discordance, paralogy, single-copy status, RF distance, reconciliation, ILS, HGT and marker selection;
- documentation defining the niche relative to GeneRax, AleRax, ASTRAL-Pro3, DiscoVista and other specialist tools;
- regression tests for evidence interpretation and causal-overclaim safeguards.

### Changed
- project positioning from a tree-comparison toolkit to an evidence-integration and phylogenomic decision-support framework;
- locus reports now state what the evidence supports, what it does not establish, and which scientific question should be tested next;
- batch summaries now expose primary concerns and top routed analyses.

## [0.1.0] - 2026-09-25

### Added
- single-locus and batch CLI;
- explicit gene-to-species mapping validation;
- taxon coverage and copy-number summaries;
- RF and split comparison for single-copy trees;
- species-overlap duplication screening;
- rooted LCA duplication-loss diagnostic;
- PASS / REVIEW / RESOLVE assessment framework;
- reference species-tree split support across full-taxon single-copy loci;
- simulated 10-genome tutorial with five scenarios;
- JSON, TSV, Markdown and tree-plot outputs;
- automated tests, CI, documentation and maintainer guidance.
