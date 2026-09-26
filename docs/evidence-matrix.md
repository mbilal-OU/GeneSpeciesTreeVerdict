# Evidence matrix

Every locus is decomposed into independent evidence domains. This is the main safeguard against turning one metric into an unjustified biological diagnosis.

## Output columns

`evidence_matrix.tsv` contains:

| Column | Meaning |
|---|---|
| `family` | locus / gene-family identifier |
| `evidence_key` | machine-readable signal name |
| `domain` | evidence domain |
| `signal` | `SUPPORT`, `NEUTRAL`, `INFO`, `FLAG`, or `STRONG_FLAG` |
| `observation` | what was measured |
| `interpretation` | what the observation supports |
| `limitation` | what the observation does not establish |

## Evidence domains

### Taxon coverage

Coverage asks whether the locus meets the configured prevalence threshold.

**It can support:** whether the locus satisfies a chosen marker-definition threshold.

**It cannot establish:** orthology, absence of duplication, or phylogenetic informativeness.

### Copy number

A multi-copy family is a direct structural warning for a conventional single-copy concatenation workflow.

A one-copy-per-species family is useful, but single-copy status alone does not prove orthology. Ancient duplication followed by differential loss can create pseudo-orthologous sampling.

### Species-overlap duplication screen

If child lineages below an internal gene-tree node contain the same species, the family has a direct repeated-species signal and should be resolved before conventional single-copy use.

This is a local screen, not a complete historical reconstruction.

### Topological comparison

For comparable single-copy trees, the software reports RF distance, normalized RF distance, shared splits, reference-only splits, gene-only splits, and reference-split recovery.

Topological discordance is an observation. It does not identify duplication, ILS, HGT, recombination, gene-tree error, or species-tree error by itself.

### Branch-support context

For conflicting gene-tree splits, GeneSpeciesTreeVerdict records whether support values are available. It does not impose one universal "strong support" cutoff because support scales and interpretations vary among methods.

### Duplication-loss reconciliation

The native rooted LCA reconciliation asks how many duplication/loss events are required under a simple DL model.

The result is conditional on:

- the supplied gene-tree topology and root;
- the supplied species-tree topology and root;
- the gene-to-species mapping;
- the duplication-loss model.

It is therefore a **screening signal**, not proof that the inferred events were the actual biological cause.

## Why the matrix matters

Consider two loci with the same nRF value of 0.60.

**Locus A** may be single-copy, have weakly supported conflicting branches, and require no duplication under reconciliation.

**Locus B** may be multi-copy in three species and contain species-overlap duplication nodes.

A single RF threshold would treat them similarly. The evidence matrix does not.

## Machine-readable use

The long-form table is intended for:

- filtering and plotting across hundreds or thousands of loci;
- reproducible review of why a locus was flagged;
- comparison of decision-rule versions;
- future calibration against simulation benchmarks;
- auditability in a methods supplement.
