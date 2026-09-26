# From a flag to the next analysis

GeneSpeciesTreeVerdict is designed to route questions to specialist methods rather than duplicate their algorithms.

The recommendations in `next_analyses.tsv` and `recommendations.tsv` are **question-driven**. Tool names are examples, not endorsements, and the researcher remains responsible for checking assumptions, versions, and suitability.

## Routing map

| Evidence pattern | Main unresolved question | Typical next-analysis class | Example tools/workflows |
|---|---|---|---|
| Multi-copy family / species-overlap duplication | Which copies are orthologous? | orthology + gene-family reconciliation | OrthoFinder; GeneRax; AleRax |
| Multi-copy gene families used to infer a species tree | Can the species tree be inferred without forcing single-copy loci? | paralogy-aware species-tree inference | ASTRAL-Pro3 |
| Single-copy topological discordance | Is the conflict well supported? | gene-tree support / sensitivity analysis | IQ-TREE/RAxML-class support analyses |
| Recurrent discordance across loci | Is the reference branch itself uncertain? | cross-locus split/quartet or species-tree analysis | GeneSpeciesTreeVerdict batch summaries; DiscoVista; ASTRAL-family methods |
| DL reconciliation requires events | Is a DTL history supported under a specialist model? | DTL reconciliation | GeneRax; AleRax; RANGER-DTL/Notung-class workflows |
| ILS is biologically plausible across many loci | Does coalescent heterogeneity explain the pattern? | multispecies-coalescent species-tree analysis | ASTRAL-family methods |
| HGT/recombination is biologically plausible | Is the locus history mosaic or transfer-like? | transfer/recombination-aware analysis | DTL-aware reconciliation and dataset-appropriate recombination/HGT methods |
| Low coverage | Why are taxa missing? | upstream genome/annotation/orthogroup QC | pangenome/orthology tables; completeness QC |
| Concordant single-copy locus | Are there other sources of phylogenetic bias? | routine sequence-level QC | alignment, trimming, model, rate/composition workflows |

## Why routing is conditional

No decision rule can choose a specialist method correctly from topology alone.

For example, a discordant single-copy microbial gene could reflect HGT, but a similar pattern in a rapid radiation could reflect ILS. Both could also be caused by weak gene-tree signal. The same RF distance therefore leads to different biological questions depending on the dataset.

GeneSpeciesTreeVerdict deliberately outputs **questions** rather than a causal label.

## Tool-specific context

### GeneRax

GeneRax performs species-tree-aware maximum-likelihood gene-family tree inference under duplication, transfer, and loss. It is appropriate when a DTL model and the sequence-level gene-family data are central to the question.

<https://github.com/BenoitMorel/GeneRax>

### AleRax

AleRax performs probabilistic gene/species-tree reconciliation from gene-tree distributions and can account for gene-tree uncertainty in a reconciliation framework.

<https://github.com/BenoitMorel/AleRax>

### ASTRAL-Pro3

ASTRAL-Pro3 is designed for species-tree inference from multi-copy gene-family trees and does not require users to force each family into a single-copy marker before species-tree inference.

<https://github.com/chaoszhang/ASTER>

### DiscoVista

DiscoVista focuses on interpretable visualization of phylogenetic discordance. It is complementary when the key problem is understanding how conflicts are distributed across branches, genes, or datasets.

<https://github.com/esayyari/DiscoVista>

### OrthoFinder

OrthoFinder is useful when the analysis begins with gene/protein sequences and orthogroup/orthology inference is still required. GeneSpeciesTreeVerdict itself does not perform full sequence-based orthology inference.

## Example routed result

```text
Locus: OG000243
Verdict: RESOLVE
Primary concern: MULTICOPY_OR_DUPLICATION_SIGNAL

Top next question:
Which homologous copies are orthologous across the sampled species?

Suggested analysis classes:
- DTL-aware gene/species-tree reconciliation
- sequence-based orthology inference if needed
- paralogy-aware species-tree inference if the goal is the species tree

Example tools:
- GeneRax / AleRax
- OrthoFinder
- ASTRAL-Pro3
```

The goal is not to say that all three should be run. The goal is to explain **which scientific question each tool would answer**.
