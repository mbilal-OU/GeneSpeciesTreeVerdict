# GeneSpeciesTreeVerdict

[![CI](https://github.com/mbilal-OU/GeneSpeciesTreeVerdict/actions/workflows/ci.yml/badge.svg)](https://github.com/mbilal-OU/GeneSpeciesTreeVerdict/actions/workflows/ci.yml)
![Python](https://img.shields.io/badge/Python-3.10--3.13-blue)
![License](https://img.shields.io/badge/License-MIT-green)
![Status](https://img.shields.io/badge/status-alpha-orange)

**From gene-tree/species-tree discordance to a transparent phylogenomic decision.**

GeneSpeciesTreeVerdict (`gstv`) is an open-source **evidence-integration and decision-support framework** for researchers who already have gene trees and a reference species tree and need to answer a practical question:

> **My gene tree disagrees with my species tree. What evidence do I actually have, what can I legitimately conclude, how should I handle this locus, and what analysis should I run next?**

The software does **not** assume the species tree is automatically correct. It does **not** call a locus paralogous simply because its topology differs from the species tree. And it does **not** attempt to replace specialist reconciliation, orthology, coalescent, HGT, or recombination methods.

Instead, it connects those analyses through a transparent workflow.

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

## What v0.2 does

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

---

## Example interpretation

Imagine a locus with one sequence per species but strong topological conflict:

```text
Verdict: REVIEW
Primary concern: MODEL_DEPENDENT_DUPLICATION_SIGNAL

Evidence
Taxon coverage             SUPPORT
Copy number                SUPPORT
Species overlap            NEUTRAL
Topology                    FLAG
DL reconciliation           FLAG

What is supported
- the locus is topologically discordant with the supplied species tree;
- a duplication-loss history can reconcile the supplied rooted trees under the native model.

What is not established
- the locus is definitely paralogous;
- duplication is definitely the biological cause;
- the supplied species tree is necessarily correct.

Top next questions
1. Is the conflicting gene-tree signal well supported?
2. Is the reference branch broadly supported across independent loci?
3. Does a specialist DTL analysis support the same explanation?
4. Are ILS, HGT/recombination, or other processes plausible in this dataset?
```

That is the intended product: **a defensible research workflow, not a one-number filter**.

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

### 3. Learn a concept from the CLI

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

The explanations use the same rule as the analysis engine:

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

The key v0.2 output is `evidence_matrix.tsv`.

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

Examples:

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

For that reason, the project should be described as an **evidence-integration and decision-support framework**, not as a replacement for specialist phylogenetic inference methods.

---

## Validation status

The software is currently **alpha**.

The test suite covers software behavior and the deterministic teaching scenarios. Broader simulation benchmarking across duplication/loss, gene-tree error, missing taxa, and other discordance regimes remains an active research objective.

A stable release, archival DOI, and software-paper claims should follow quantitative benchmarking, at least one empirical case study, and external-user validation.

See the public roadmap and issues before interpreting this repository as a fully validated phylogenetic method.

---

## Documentation

- [Why this tool exists](docs/why-this-tool.md)
- [Quick start](docs/quickstart.md)
- [Evidence matrix](docs/evidence-matrix.md)
- [Tool routing](docs/tool-routing.md)
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

> **GeneSpeciesTreeVerdict integrates gene-tree/species-tree evidence, states what that evidence can and cannot establish, and routes each locus toward a transparent phylogenomic workflow decision.**
