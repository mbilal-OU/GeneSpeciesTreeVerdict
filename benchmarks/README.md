# GeneSpeciesTreeVerdict benchmarks

This directory contains reproducible validation configurations for the decision-support layer.

## Two validation levels

### 1. Native known-truth structural stress tests

The built-in benchmark generator creates seeded species trees and controlled gene-tree perturbations where copy number, taxon coverage, and topology manipulation are known exactly. It is designed to test whether GeneSpeciesTreeVerdict makes safe workflow decisions.

Run the fast configuration:

```bash
gstv benchmark --config benchmarks/smoke.yaml --outdir benchmark_results/smoke
```

Run the larger validation grid:

```bash
gstv benchmark --config benchmarks/full.yaml --outdir benchmark_results/full
```

The full configuration contains 4 taxon counts × 7 scenarios × 100 replicates = 2,800 locus-level benchmark cases.

### 2. External process-realistic simulations

The native generator is deliberately not described as a full biological simulator of ILS, HGT, recombination, or sequence evolution. For process-level validation, simulate data using a specialist simulator and provide a truth manifest:

```bash
gstv benchmark --truth-manifest path/to/truth_manifest.tsv --outdir benchmark_results/external
```

Good external targets include:

- **SimPhy** — hierarchical species/locus/gene-tree simulation with incomplete lineage sorting, gene duplication/loss, HGT and gene conversion.
- **Zombi** — species-tree, genome and sequence simulation with duplication, loss, transfer, originations and genome rearrangements.

The SpeciesRax validation study is an example of using SimPhy to benchmark species-tree methods across DTL rates, species counts and gene-family properties.

## Native scenarios

| Scenario | Controlled truth | Expected workflow state |
|---|---|---|
| `concordant_singlecopy` | single copy, matching topology | PASS |
| `branch_length_only` | matching topology, changed branch lengths | PASS |
| `singlecopy_discordance` | single copy, topology conflict, no simulated duplication | REVIEW |
| `multicopy_duplication` | explicit sampled duplication | RESOLVE |
| `pseudoorthology` | one sampled copy/species from deep hidden lineages | REVIEW |
| `missing_taxa` | coverage below threshold | REVIEW |
| `gene_tree_error` | single-copy topology perturbation, no simulated duplication | REVIEW |

The expected labels are **workflow states**, not biological truth labels.

## Outputs

Each benchmark writes:

- `truth_manifest.tsv` — exact seeds, truth metadata and input paths;
- `locus_results.tsv` — observed GSTV decisions and metrics;
- `performance_summary.tsv` — headline safety/performance metrics;
- `scenario_metrics.tsv` — accuracy by scenario;
- `class_metrics.tsv` — precision, recall and F1 for PASS/REVIEW/RESOLVE;
- `confusion_matrix.tsv` — expected-vs-observed workflow states;
- `failure_cases.tsv` — cases that violated the expected workflow state;
- `scenario_accuracy.png` — quick visual summary;
- `benchmark_report.md` — human-readable report.

## Safety metrics

Two metrics are especially important for this project:

1. **RESOLVE sensitivity for true sampled multicopy duplication** — does the framework catch families that should not be treated as conventional single-copy markers?
2. **False-RESOLVE rate when no duplication was simulated** — does topology discordance accidentally become a paralogy-like exclusion decision?

The second metric directly tests the core scientific safeguard of the project.
