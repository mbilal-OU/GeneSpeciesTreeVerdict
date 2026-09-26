# Why GeneSpeciesTreeVerdict exists

GeneSpeciesTreeVerdict is **not another gene-tree/species-tree reconciliation algorithm**.

It is an evidence-integration and decision-support layer for a practical question that remains after many specialist tools have done their job:

> **My gene tree and species tree disagree. What evidence do I actually have, what can I legitimately conclude, how should I handle this locus, and what analysis should I run next?**

## The gap it is designed to fill

Specialist phylogenetic programs usually solve a well-defined inference problem. Examples include:

| Tool / method class | Primary task |
|---|---|
| GeneRax | species-tree-aware gene-family tree inference and DTL reconciliation |
| AleRax | probabilistic gene/species-tree reconciliation using gene-tree distributions |
| ASTRAL / ASTRAL-Pro3 | species-tree inference from gene trees, including multi-copy families in ASTRAL-Pro3 |
| Notung / RANGER-DTL-class tools | gene/species-tree reconciliation under explicit event models |
| DiscoVista | visualization and interpretation of phylogenetic discordance |

These are complementary to GeneSpeciesTreeVerdict.

GeneSpeciesTreeVerdict instead asks whether the evidence available for a locus supports a workflow decision such as:

- continue as a candidate conventional single-copy marker;
- review a supported or model-dependent conflict;
- resolve a multi-copy family before concatenation;
- reconsider a poorly supported reference species-tree branch;
- route the locus to a specialist reconciliation, coalescent, orthology, or recombination-aware analysis.

## What makes the workflow different

The software keeps four layers separate:

1. **Observation** — what was measured or seen.
2. **Interpretation** — what that observation can reasonably support.
3. **Limitation** — what the observation does *not* establish.
4. **Next question** — what analysis would discriminate among remaining explanations.

For example, a high normalized RF distance is recorded as topological discordance. It is **not** translated into "paralogy". A duplication-loss reconciliation that requires a duplication is recorded as a model-dependent reconciliation signal. It is **not** translated into "a biological duplication has been proven".

## The central design principle

```text
measurement != mechanism != decision
```

A gene tree can disagree with a species tree because of duplication/loss, incomplete lineage sorting, horizontal transfer, recombination, gene-tree error, model misspecification, taxon sampling, or error/uncertainty in the reference species tree.

GeneSpeciesTreeVerdict therefore does not use a rule such as:

```text
RF > 0 -> paralog -> remove
```

Instead, it integrates copy number, prevalence, topology, branch-support context, reconciliation evidence, and cross-locus reference-tree support.

## Example

A single-copy locus might produce:

```text
Verdict: REVIEW
Primary concern: MODEL_DEPENDENT_DUPLICATION_SIGNAL

Evidence
- Coverage: 100%                          SUPPORT
- Copy number: one sequence/species       SUPPORT
- Species-overlap duplication: none       NEUTRAL
- Topology: nRF = 0.57                    FLAG
- DL-LCA reconciliation: 1 duplication    FLAG

Supported conclusion
- The locus is topologically discordant.
- A duplication-loss history can reconcile the supplied rooted trees under the native DL model.

Not established
- The locus is definitely paralogous.
- The conflict is definitely caused by duplication.
- The supplied species tree is necessarily correct.

Next questions
1. Is the conflicting gene-tree signal strongly supported?
2. Is the reference branch supported across independent loci?
3. Does a specialist DTL reconciliation support the same explanation?
4. Are ILS, HGT/recombination, or other processes plausible for this dataset?
```

That is the intended contribution: **turning heterogeneous phylogenetic evidence into a transparent research workflow without overclaiming causal diagnosis**.

## What GeneSpeciesTreeVerdict does not replace

It does not replace:

- sequence alignment and alignment QC;
- substitution-model testing;
- full orthology inference from sequences;
- probabilistic DTL reconciliation;
- multispecies-coalescent species-tree inference;
- recombination/HGT detection;
- biological interpretation by the researcher.

When one of those analyses is needed, the program should make that need explicit and route the user toward the appropriate method class.

## Current maturity

Version 0.2 is an **alpha decision-support framework**. Its rules are intentionally transparent and conservative. The built-in tutorial provides known structural scenarios, while broader simulation benchmarking and empirical validation remain active development goals.

Do not interpret `PASS`, `REVIEW`, or `RESOLVE` as statements of evolutionary truth. They are workflow states relative to the supplied data, reference tree, and configured analysis.

## Specialist tools referenced in the documentation

- GeneRax: <https://github.com/BenoitMorel/GeneRax>
- AleRax: <https://github.com/BenoitMorel/AleRax>
- ASTER / ASTRAL-Pro3: <https://github.com/chaoszhang/ASTER>
- DiscoVista: <https://github.com/esayyari/DiscoVista>
- Notung: <https://www.cs.cmu.edu/~durand/Notung/>

These projects have different goals and assumptions. Their inclusion here is for workflow routing and comparison, not endorsement or equivalence.
