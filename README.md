# GeneSpeciesTreeVerdict

[![CI](https://github.com/mbilal-OU/GeneSpeciesTreeVerdict/actions/workflows/ci.yml/badge.svg)](https://github.com/mbilal-OU/GeneSpeciesTreeVerdict/actions/workflows/ci.yml)
![Python](https://img.shields.io/badge/Python-3.10--3.13-blue)
![License](https://img.shields.io/badge/License-MIT-green)
![Status](https://img.shields.io/badge/status-alpha-orange)

**From gene-tree/species-tree discordance to a transparent, testable phylogenomic decision.**

GeneSpeciesTreeVerdict (`gstv`) is an open-source **evidence-integration, validation, and decision-support framework** for researchers who already have gene trees and a reference species tree and need to answer a practical question:

> **My gene tree disagrees with my species tree. What evidence do I actually have, what can I legitimately conclude, how should I handle this locus, and what analysis should I run next?**

The software does **not** assume the supplied species tree is automatically correct. It does **not** call a locus paralogous simply because its topology differs from the species tree. It does **not** replace specialist reconciliation, orthology, coalescent, HGT, or recombination methods.

Instead, it connects those analyses through an auditable workflow and now includes a benchmark layer that can test the workflow against known truth.

---

## Why this is different

Many excellent tools already solve specialist phylogenetic problems:

| Tool / method | Main task |
|---|---|
| [GeneRax](https://github.com/BenoitMorel/GeneRax) | species-tree-aware gene-family tree inference and DTL reconciliation |
| [AleRax](https://github.com/BenoitMorel/AleRax) | probabilistic gene/species-tree reconciliation using gene-tree distributions |
| [ASTRAL / ASTRAL-Pro3](https://github.com/chaoszhang/ASTER) | species-tree inference from gene trees, including multi-copy families in ASTRAL-Pro3 |
| [Notung](https://www.cs.cmu.edu/~durand/Notung/) / RANGER-DTL-class methods | reconciliation under explicit evolutionary-event models |
| [DiscoVista](https://github.com/esayyari/DiscoVista) | visualization of phylogenetic discordance |
| [SimPhy](https://academic.oup.com/sysbio/article/65/2/334/2427219) | process-realistic phylogenomic simulation including ILS, duplication/loss and HGT |
| [Zombi](https://academic.oup.com/bioinformatics/article/36/4/1286/5578480) | species-tree, genome and sequence simulation with duplication/loss/transfer |

GeneSpeciesTreeVerdict asks a different question:

> **Given all the evidence I have for this locus, what does it support, what does it not prove, and what should I investigate next?**

The core design principle is:

```text
measurement != mechanism != decision
```

A normalized RF distance measures topological difference. It does not identify its cause.

A duplication-loss reconciliation can show that a DL history explains the supplied trees under a model. It does not prove that duplication was the true biological cause.

One sequence per species is useful. It does not guarantee orthology.

That separation is the reason this project exists.

Read the full rationale: **[Why GeneSpeciesTreeVerdict exists](docs/why-this-tool.md)**.

---

## What v0.3 does

For every locus, GeneSpeciesTreeVerdict integrates:

```text
                           LOCUS
                             |
          +------------------+------------------+
          |                  |                  |
       Sampling          Copy number         Topology
          |                  |                  |
       coverage        repeated copies      RF / splits
          |                  |                  |
          +------------------+------------------+
                             |
                   duplication screens
                             |
                   rooted DL reconciliation
                             |
                     evidence matrix
                             |
          +------------------+------------------+
          |                  |                  |
   what is supported   what is NOT proven   next question
          |                  |                  |
          +------------------+------------------+
                             |
                    PASS / REVIEW / RESOLVE
                             |
                 specialist method if needed
                             |
                       benchmark layer
                             |
                 known truth vs observed
```

### Per-locus evidence

- taxon coverage;
- copy number by species;
- multi-copy species;
- species-overlap duplication candidates;
- RF and normalized RF distance for comparable single-copy trees;
- shared, reference-only, and gene-only splits;
- reference-split recovery;
- available support values on conflicting gene-tree splits;
- rooted duplication-loss LCA reconciliation;
- explicit interpretation limitations.

### Decision support

Every locus receives:

- a workflow state: `PASS`, `REVIEW`, or `RESOLVE`;
- a machine-readable `primary_concern`;
- statements the evidence **supports**;
- statements the evidence **does not establish**;
- a recommended action;
- prioritized scientific questions to test next;
- example specialist tools/workflows appropriate to those questions.

### Reverse-testing the species tree

Batch mode does not merely ask whether genes agree with the reference tree. For full-taxon single-copy loci, it also asks:

> **How consistently do independent gene trees support each branch of the supplied species tree?**

This helps distinguish “one problematic locus” from “a reference branch that is repeatedly disputed across loci.”

### Validation layer

v0.3 adds a reproducible benchmark system with two deliberately separate levels:

1. **Native known-truth structural stress tests** — seeded tests where copy number, coverage and topology manipulations are known exactly.
2. **External truth-manifest mode** — lets independently simulated datasets from tools such as SimPhy or Zombi be scored without pretending GSTV's lightweight generator is a complete evolutionary simulator.

This separation reduces circular validation.

---

## Installation

The project is currently an alpha research package and has not yet been declared a stable PyPI/Bioconda release.

Clone and install from GitHub:

```bash
git clone https://github.com/mbilal-OU/GeneSpeciesTreeVerdict.git
cd GeneSpeciesTreeVerdict
python -m pip install -e .
```

For development:

```bash
python -m pip install -e ".[dev]"
```

Check the installation:

```bash
gstv version
```

---

## Quick start

### 1. Analyze one gene tree

```bash
gstv analyze \
  --species-tree species_tree.nwk \
  --gene-tree gene001.nwk \
  --mapping mapping.tsv \
  --outdir gene001_result
```

Main outputs:

```text
gene001_result/
├── result.json
├── summary.tsv
├── evidence_matrix.tsv
├── next_analyses.tsv
├── report.md
└── tree_comparison.png
```

### 2. Analyze many loci

```bash
gstv batch \
  --species-tree species_tree.nwk \
  --gene-tree-dir gene_trees/ \
  --mapping mapping.tsv \
  --outdir gstv_batch
```

Important batch outputs:

```text
gstv_batch/
├── locus_summary.tsv
├── evidence_matrix.tsv
├── recommendations.tsv
├── reference_split_support.tsv
├── alternative_splits.tsv
├── results.json
└── report.md
```

### 3. Explain a concept

```bash
gstv explain list
gstv explain paralogy
gstv explain single-copy
gstv explain discordance
gstv explain rf
gstv explain reconciliation
gstv explain ils
gstv explain hgt
gstv explain marker-selection
```

The explanations follow the same rule as the analysis engine:

> observation first, causal claim only when justified.

### 4. Generate the teaching dataset

```bash
gstv tutorial --outdir gstv_tutorial
```

The deterministic 10-genome tutorial includes:

- a concordant single-copy ortholog;
- a discordant single-copy locus;
- a recent duplication;
- an ancient duplication with differential retention;
- a mixed pseudo-ortholog example demonstrating why one sequence per genome does not guarantee orthology.

Then run:

```bash
gstv batch \
  --species-tree gstv_tutorial/species_tree.nwk \
  --gene-tree-dir gstv_tutorial/gene_trees \
  --mapping gstv_tutorial/mapping.tsv \
  --outdir gstv_tutorial/results
```

### 5. Benchmark the decision engine

Run the fast reproducible benchmark:

```bash
gstv benchmark \
  --config benchmarks/smoke.yaml \
  --outdir benchmark_results/smoke
```

Run the larger validation grid:

```bash
gstv benchmark \
  --config benchmarks/full.yaml \
  --outdir benchmark_results/full
```

The full configuration currently contains:

```text
4 taxon counts × 7 scenarios × 100 replicates = 2,800 benchmark cases
```

Use externally simulated truth:

```bash
gstv benchmark \
  --truth-manifest external_simulation/truth_manifest.tsv \
  --outdir benchmark_results/external
```

Benchmark outputs include:

```text
benchmark_results/
├── truth_manifest.tsv
├── locus_results.tsv
├── performance_summary.tsv
├── scenario_metrics.tsv
├── class_metrics.tsv
├── confusion_matrix.tsv
├── failure_cases.tsv
├── scenario_accuracy.png
└── benchmark_report.md
```

Key safety metrics include:

- `RESOLVE` sensitivity for explicit sampled multicopy duplication;
- false-`RESOLVE` rate when no duplication was simulated;
- false-`RESOLVE` rate for single-copy/no-duplication truth.

The last two directly test the project's central safeguard: **discordance must not silently become a paralogy decision**.

See **[Benchmarking and scientific validation](docs/benchmarking.md)**.

---

## Understanding the workflow states

### `PASS`

The locus is structurally compatible with the configured conventional single-copy marker workflow.

`PASS` does **not** mean:

- proven orthology;
- biologically perfect locus;
- no alignment/model bias;
- “true gene tree.”

### `REVIEW`

The locus contains evidence that requires interpretation before a marker decision—for example low coverage, supported topology conflict, or a model-dependent reconciliation signal.

`REVIEW` is intentionally broad because the program refuses to turn generic discordance into a unique causal label.

### `RESOLVE`

The family contains direct multi-copy/duplication structure that should be resolved before conventional single-copy concatenation.

Depending on the scientific goal, the next step may instead be to use a method designed for multi-copy families.

---

## Evidence matrix

The key evidence output is `evidence_matrix.tsv`.

Each row contains:

```text
family
evidence_key
domain
signal
observation
interpretation
limitation
```

This makes the decision auditable. A reviewer can see *why* a locus was flagged and *which assumption* limits that inference.

See **[Evidence matrix](docs/evidence-matrix.md)**.

---

## Question-driven tool routing

GeneSpeciesTreeVerdict does not say “run GeneRax” simply because a tree differs.

It asks what remains unresolved.

| Evidence pattern | Next question | Example method class |
|---|---|---|
| multi-copy family | Which copies are orthologous? | orthology / DTL reconciliation |
| single-copy conflict | Is the conflict strongly supported? | gene-tree support / sensitivity analysis |
| recurrent conflict across loci | Is the reference branch uncertain? | split/quartet / species-tree analysis |
| DL screening signal | Does a specialist DTL model support this history? | GeneRax / AleRax-class reconciliation |
| plausible ILS across loci | Does coalescent heterogeneity explain the pattern? | ASTRAL-family methods |
| low coverage | Why are taxa missing? | upstream genome/annotation/orthogroup QC |

See **[Tool routing](docs/tool-routing.md)**.

---

## Input mapping

Gene-tree leaves must map unambiguously to species-tree taxa.

A mapping table can use either:

```text
gene_id    species_id
```

or:

```text
family    gene_id    species_id
```

The three-column form is recommended for batch analyses with reused sequence identifiers.

GeneSpeciesTreeVerdict fails loudly on ambiguous or missing mappings rather than guessing.

---

## Scientific boundaries

GeneSpeciesTreeVerdict currently does **not** claim to perform:

- full sequence-based orthology inference;
- probabilistic DTL inference comparable to GeneRax/AleRax;
- formal ILS-vs-HGT model selection;
- recombination detection;
- species-tree inference;
- universal branch-support calibration;
- automatic biological truth assignment.

The native DL reconciliation is deliberately a transparent screening model.

The native benchmark is deliberately a **structural stress test**, not a complete generative model of ILS, HGT, recombination or sequence evolution.

For that reason, the project should be described as an **evidence-integration, validation, and decision-support framework**, not as a replacement for specialist phylogenetic inference methods.

---

## Validation status

The software is currently **alpha**.

v0.3 provides reproducible native known-truth benchmarking, exact simulation seeds, failure-case reporting and an external simulator interface. This is an important step beyond software unit testing, but it is not the end of scientific validation.

Before a stable methodological claim or software paper, the project should additionally include:

1. replicated SimPhy experiments across ILS, duplication/loss/transfer rates, taxa and gene-tree error;
2. independent microbial DTL simulation with a simulator such as Zombi;
3. sequence simulation followed by gene-tree re-estimation;
4. empirical benchmark datasets with independently curated orthology/paralogy or reconciliation evidence;
5. explicit reporting of failure regimes and ambiguous cases;
6. external-user validation.

A stable release, archival DOI and software-paper claims should follow those validation steps rather than precede them.

---

## Documentation

- [Why this tool exists](docs/why-this-tool.md)
- [Quick start](docs/quickstart.md)
- [Evidence matrix](docs/evidence-matrix.md)
- [Tool routing](docs/tool-routing.md)
- [Benchmarking and validation](docs/benchmarking.md)
- [Core concepts](docs/concepts.md)
- [Simulated tutorial](docs/tutorial.md)
- [Real-data workflow](docs/real-data-workflow.md)
- [Decision framework](docs/decision-framework.md)
- [Interpreting discordance](docs/interpretation.md)
- [Scientific background](docs/scientific-background.md)
- [Limitations](docs/limitations.md)
- [FAQ](docs/faq.md)
- [Maintainer guide](docs/maintainer-guide.md)

---

## Development and reproducibility

The repository uses:

- Python 3.10–3.13 CI;
- `pytest` and coverage;
- Ruff linting/formatting;
- a benchmark smoke run inside CI;
- strict MkDocs builds;
- package build validation;
- Dependabot;
- issue and PR templates;
- scientific-correctness issue template;
- `CITATION.cff`;
- changelog and roadmap;
- contributor, security, and maintainer guidance.

Run the local quality suite with:

```bash
make test
make lint
make docs
make build
```

---

## Contributing

Scientific criticism is particularly welcome.

If you believe a rule overinterprets an evolutionary signal, use the **scientific correctness** issue template and provide a minimal example, relevant assumptions, and supporting references where possible.

See [CONTRIBUTING.md](CONTRIBUTING.md).

---

## Citation

Until an archival release/DOI is created, cite the exact software version or Git commit used in an analysis. See [`CITATION.cff`](CITATION.cff).

---

## License

MIT License. See [LICENSE](LICENSE).

---

## One-sentence scope

> **GeneSpeciesTreeVerdict integrates gene-tree/species-tree evidence, states what that evidence can and cannot establish, benchmarks its workflow decisions against known truth, and routes each locus toward a transparent phylogenomic next step.**
