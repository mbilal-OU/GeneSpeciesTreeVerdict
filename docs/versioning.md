# Versioning and scientific stability

GeneSpeciesTreeVerdict follows semantic versioning, with additional care for changes that alter scientific interpretation.

- **Patch** releases fix defects without intentionally changing the public API or scientific decision policy.
- **Minor** releases add backwards-compatible functionality or new diagnostics.
- **Major** releases may change CLI/API contracts, output schemas, or core interpretation rules.

During the `0.x` phase, interfaces may still evolve. Any change to a threshold, decision rule, reconciliation assumption, or default that can change locus verdicts must be highlighted in the changelog and release notes, even if the software version remains pre-1.0.

For reproducible research, record the exact GeneSpeciesTreeVerdict version, command, parameters, reference species tree, and mapping used in an analysis.
