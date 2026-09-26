# GeneSpeciesTreeVerdict

[![CI](https://github.com/mbilal-OU/GeneSpeciesTreeVerdict/actions/workflows/ci.yml/badge.svg)](https://github.com/mbilal-OU/GeneSpeciesTreeVerdict/actions/workflows/ci.yml)
![Python](https://img.shields.io/badge/Python-3.10--3.13-blue)
![License](https://img.shields.io/badge/License-MIT-green)
![Status](https://img.shields.io/badge/status-alpha-orange)

**Gene tree vs. species tree: diagnose the conflict before deciding how to use the locus.**

GeneSpeciesTreeVerdict is an open-source phylogenomic diagnostic toolkit for comparing gene trees with a reference species tree, screening copy-number and duplication signals, and producing a transparent **locus-handling assessment** for conventional phylogenomic workflows.

It is built around a question that appears constantly in phylogenomics:

> **My gene tree disagrees with my species tree. What does that mean, and should I use this locus?**

The software does **not** assume that the species tree is automatically correct, and it does **not** call a gene paralogous merely because its topology differs from the species tree.

---

## Why this tool exists

Gene-tree/species-tree discordance can arise from several distinct processes, including gene duplication and loss, horizontal gene transfer, incomplete lineage sorting, recombination, gene-tree uncertainty, model misspecification, taxon sampling, or an incorrect reference species tree. A simple "same tree / different tree" test is therefore not enough.

GeneSpeciesTreeVerdict separates the evidence into layers:

```text
                    Reference species tree
                             +
                         Gene tree(s)
                             |
         +-------------------+-------------------+
         |                   |                   |
     Taxon coverage      Copy number         Topology
         |                   |                   |
         |            species overlap       RF / splits
         |                   |                   |
         +-------------------+-------------------+
                             |
                    Rooted DL reconciliation
                             |
                             v
                  Transparent locus assessment

                  PASS | REVIEW | RESOLVE
```

The assessment answers a **workflow question**: how should this locus be handled in a conventional single-copy species-tree analysis? It is not a declaration of biological truth.

---

## What v0.1 does

For each gene family, GeneSpeciesTreeVerdict reports:

- taxon coverage relative to the reference species tree;
- copy number per species and multi-copy taxa;
- candidate duplication nodes from child-lineage species overlap;
- unrooted or rooted split comparison;
- Robinson-Foulds (RF) distance and normalized RF distance;
- reference splits recovered by the gene tree;
- conflicting gene-tree splits and their support values when available;
- rooted least-common-ancestor duplication-loss reconciliation;
- model-conditional duplication and inferred loss counts;
- a `PASS`, `REVIEW`, or `RESOLVE` locus-handling assessment;
- an explicit explanation of **why** that assessment was produced.

Batch mode also asks the reverse question:

> **How consistently do full-taxon single-copy gene trees support each branch of the supplied species tree?**

This produces `reference_split_support.tsv` and recurring alternative splits, so the reference tree itself can be scrutinized rather than treated as unquestionable truth.

---

## The three verdicts

### `PASS`

The locus is single-copy across mapped taxa and meets the configured coverage threshold. It is structurally suitable as a candidate conventional phylogenomic marker.

**Important:** a `PASS` locus may still have a discordant gene tree. Discordance is reported rather than automatically filtered.

### `REVIEW`

The locus needs interpretation before inclusion. Examples include:

- lower-than-target taxon coverage;
- single-copy gene-tree/species-tree discordance that requires duplication events under a DL reconciliation model;
- user-defined high-discordance review thresholds.

A `REVIEW` result does **not** mean "paralog."

### `RESOLVE`

The family is multi-copy and/or has direct species-overlap duplication evidence. The complete family should not be treated as an ordinary single-copy marker until orthology/paralogy is resolved, or a paralogy-aware method is used.

---

## Scientific rule at the center of the project

> **Gene-tree/species-tree discordance alone is not evidence of paralogy.**

A rooted duplication-loss reconciliation can explain discordance by invoking duplications and losses, but that explanation is conditional on the DL model. The same topological conflict may instead reflect ILS, HGT, recombination, tree-estimation error, or a problem with the reference species tree.

GeneSpeciesTreeVerdict therefore reports these signals separately rather than collapsing them into one binary label.

---

## Installation

### Development install

```bash
git clone https://github.com/mbilal-OU/GeneSpeciesTreeVerdict.git
cd GeneSpeciesTreeVerdict
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
python -m pip install -U pip
python -m pip install -e ".[dev]"
```

The runtime requires Python 3.10+ and uses Biopython for Newick parsing and tree handling.

---

## Start with the simulated tutorial

Generate the built-in 10-genome tutorial:

```bash
gstv tutorial --outdir tutorial_demo
```

Then analyze every simulated family:

```bash
gstv batch \
  --species-tree tutorial_demo/species_tree.nwk \
  --gene-tree-dir tutorial_demo/gene_trees \
  --mapping tutorial_demo/mapping.tsv \
  --outdir tutorial_demo/results
```

The tutorial contains five deliberately different evolutionary situations:

| Scenario | Biological point | Expected handling |
|---|---|---|
| `geneA_ortholog_concordant` | single-copy ortholog, matching topology | `PASS` |
| `geneD_singlecopy_discordant` | one copy/species, but discordant | `REVIEW`, not automatic paralogy |
| `geneB_recent_duplication` | recent extra copy in one species | `RESOLVE` |
| `geneC_ancient_duplication_loss` | duplicated family + differential retention | `RESOLVE` |
| `geneP_mixed_pseudoortholog` | one tip/species but inconsistent ancient copies | `REVIEW` |

See [the tutorial explanation](docs/tutorial.md) for the full visual logic.

---

## Analyze one real gene tree

If gene-tree tip labels are identical to species-tree tip labels:

```bash
gstv analyze \
  --species-tree species_tree.nwk \
  --gene-tree gene001.treefile \
  --outdir gene001_verdict
```

If gene tips are sequence IDs, provide an explicit mapping:

```text
gene_id<TAB>species_id
seq_0001<TAB>Genome_001
seq_0002<TAB>Genome_002
```

Then:

```bash
gstv analyze \
  --species-tree species_tree.nwk \
  --gene-tree gene001.treefile \
  --mapping gene001_mapping.tsv \
  --outdir gene001_verdict
```

Outputs include:

```text
gene001_verdict/
├── result.json
├── summary.tsv
├── report.md
└── tree_comparison.png
```

---

## Analyze hundreds or thousands of gene trees

```bash
gstv batch \
  --species-tree species_tree.nwk \
  --gene-tree-dir gene_trees/ \
  --mapping mappings.tsv \
  --min-coverage 0.95 \
  --outdir verdict_results/
```

A combined mapping may use three columns:

```text
family<TAB>gene_id<TAB>species_id
OG000001<TAB>seq1<TAB>Genome_001
OG000001<TAB>seq2<TAB>Genome_002
OG000002<TAB>seq3<TAB>Genome_001
```

Batch outputs:

```text
verdict_results/
├── locus_summary.tsv
├── results.json
├── report.md
├── reference_split_support.tsv
├── alternative_splits.tsv
└── errors.tsv                     # only if some loci fail validation
```

`reference_split_support.tsv` is intentionally limited to **full-taxon, single-copy** loci so that partially sampled trees are not treated as if they tested branches they cannot actually resolve.

---

## RF distance is evidence, not a verdict

GeneSpeciesTreeVerdict computes split-based RF distance, but **does not use RF disagreement alone to reject a locus**.

By default, discordance is reported only. If a project has a pre-defined analysis policy, a user can explicitly request a review threshold:

```bash
gstv analyze \
  --species-tree species.nwk \
  --gene-tree gene.nwk \
  --discordance-review-threshold 0.50
```

Threshold choice is intentionally left to the study design because the interpretation of discordance depends on taxon number, gene-tree uncertainty, evolutionary process, and the downstream method.

---

## Rooting matters

RF comparison is unrooted by default:

```bash
gstv analyze --species-tree species.nwk --gene-tree gene.nwk
```

Rooted split comparison can be requested with:

```bash
--rooted-comparison
```

The duplication-loss LCA reconciliation interprets the input species-tree and gene-tree roots as biologically meaningful. If the roots are uncertain, use the reconciliation counts as hypotheses to investigate, not definitive event calls.

---

## Where this fits among established tools

GeneSpeciesTreeVerdict is designed as a **diagnostic and triage layer**, not as a replacement for full probabilistic or parsimony reconciliation/inference packages.

| Tool / method | Primary role | Relationship to GeneSpeciesTreeVerdict |
|---|---|---|
| **GeneRax** | species-tree-aware gene-family inference under duplication, transfer and loss | follow-up when a full DTL likelihood model is needed |
| **AleRax** | gene/species-tree co-estimation and reconciliation under a probabilistic DTL model | follow-up for probabilistic reconciliation/co-estimation |
| **Treerecs** | gene-tree rooting/correction and duplication-loss reconciliation | follow-up for DL-aware correction/reconciliation |
| **Notung** | reconciliation and gene-tree rearrangement under evolutionary-event models | follow-up for event-aware reconciliation |
| **ASTRAL-Pro / ASTRAL-Pro3** | species-tree inference from multi-copy gene-family trees | alternative when paralogous families should contribute to species-tree inference |
| **ETE / TreeKO** | tree analysis and duplication-aware tree comparison | conceptually related to duplication-aware comparison |
| **GeneSpeciesTreeVerdict** | transparent QC, discordance diagnosis, per-locus handling, and reference-branch evidence | front-end diagnostic layer |

See [Scientific background](docs/scientific-background.md) for citations and methodological boundaries.

---

## What this software does **not** claim

GeneSpeciesTreeVerdict v0.1 does not:

- prove that the supplied species tree is the true species history;
- infer HGT from topology alone;
- infer ILS from topology alone;
- replace multispecies-coalescent methods;
- replace GeneRax/AleRax/Notung/Treerecs/ASTRAL-Pro;
- determine orthology solely from RF distance;
- turn a multi-copy family into a single-copy marker automatically;
- simulate nucleotide or amino-acid sequence evolution yet.

These boundaries are features, not omissions: the software is designed to report what the input evidence can actually support.

---

## Documentation

- [Quick start](docs/quickstart.md)
- [Core concepts](docs/concepts.md)
- [Simulated tutorial](docs/tutorial.md)
- [Real-data workflow](docs/real-data-workflow.md)
- [Input formats](docs/input-formats.md)
- [Decision framework](docs/decision-framework.md)
- [How to interpret discordance](docs/interpretation.md)
- [Scientific background and related tools](docs/scientific-background.md)
- [Limitations](docs/limitations.md)
- [FAQ](docs/faq.md)
- [Maintainer guide](docs/maintainer-guide.md)
- [Roadmap](ROADMAP.md)

---

## Development

```bash
make install
make test
make lint
make docs
```

The repository uses:

- `pytest` for automated tests;
- `ruff` for linting and formatting;
- GitHub Actions for Python-version matrix CI;
- MkDocs Material for documentation;
- Dependabot for dependency maintenance;
- `CITATION.cff` for software citation metadata;
- issue and pull-request templates for reproducible support and contribution workflows.

See [CONTRIBUTING.md](CONTRIBUTING.md) and the [maintainer guide](docs/maintainer-guide.md).

---

## Citation

This project is under active early development. If you use an unreleased version, cite the repository and record the exact commit SHA. A `CITATION.cff` file is included so GitHub can generate citation metadata.

A versioned release and archival DOI are planned after the initial API and validation benchmark are stabilized.

---

## License

MIT License. See [LICENSE](LICENSE).

---

## Project philosophy

> **Presence makes a gene core; copy number and evolutionary history determine whether it is an appropriate conventional phylogenomic marker.**

> **A gene tree that disagrees with a species tree is a biological or analytical question to investigate—not a paralog by definition.**
