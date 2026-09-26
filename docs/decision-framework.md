# Decision framework

GeneSpeciesTreeVerdict is intentionally not a binary "gene tree right / species tree right" machine. It returns a workflow assessment relative to the supplied species-tree hypothesis.

## PASS
Conditions in v0.1:
- one mapped copy per species among observed taxa;
- configured taxon-coverage threshold met;
- no direct species-overlap duplication candidate;
- no rooted DL reconciliation signal that triggers review.

A `PASS` does not require RF = 0. Gene-tree discordance is reported separately.

## REVIEW
Typical reasons:
- taxon coverage below the configured threshold;
- single-copy discordance that requires duplication(s) under the simple DL model;
- user-selected normalized-RF review threshold exceeded.

`REVIEW` means "investigate before inclusion," not "bad gene" and not "paralog."

## RESOLVE
Triggered by direct structural evidence such as:
- more than one mapped gene copy in a species;
- child-lineage species overlap at an internal gene-tree node.

Recommended action: resolve orthology/paralogy or use a method that explicitly accepts multi-copy families before conventional concatenation.

## Why there is no automatic EXCLUDE in v0.1
Automatic exclusion is deliberately avoided in the first release because exclusion depends on study design and downstream method. A multicopy family may be unusable for ordinary single-copy concatenation yet valuable for duplication/loss analysis or multi-copy species-tree methods.

## Thresholds
There is no default normalized-RF cutoff for rejection. RF depends on taxon sampling and does not identify the biological cause of discordance. A study-specific review threshold can be supplied explicitly when justified in advance.
