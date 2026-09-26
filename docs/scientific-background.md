# Scientific background and related methods

## Why gene trees and species trees differ

Gene-tree/species-tree discordance is expected under multiple evolutionary processes. Maddison's classic formulation emphasized lineage sorting, duplication/extinction and horizontal processes as reasons gene histories need not reproduce the species history exactly.

- Maddison WP. 1997. **Gene Trees in Species Trees.** *Systematic Biology* 46:523–536. DOI: [10.1093/sysbio/46.3.523](https://doi.org/10.1093/sysbio/46.3.523)

Reconciliation formalizes the mapping of gene histories into a species history under an explicit event model. Duplication-transfer-loss models explain discordance with gene duplication, horizontal transfer and loss, while extensions can incorporate incomplete lineage sorting.

- Bansal MS, Alm EJ, Kellis M. 2012. **Efficient algorithms for the reconciliation problem with gene duplication, horizontal transfer and loss.** *Bioinformatics* 28:i283–i291. DOI: [10.1093/bioinformatics/bts225](https://doi.org/10.1093/bioinformatics/bts225)
- Stolzer M et al. 2012. **Inferring duplications, losses, transfers and incomplete lineage sorting with nonbinary species trees.** *Bioinformatics* 28:i409–i415. DOI: [10.1093/bioinformatics/bts386](https://doi.org/10.1093/bioinformatics/bts386)

## Established software ecosystem

### GeneRax
GeneRax performs species-tree-aware maximum-likelihood gene-family tree inference under duplication, transfer and loss. It is appropriate when the study requires a full DTL-aware gene-tree inference framework rather than a lightweight diagnostic screen.

- Morel B, Kozlov AM, Stamatakis A, Szöllősi GJ. 2020. **GeneRax: A Tool for Species-Tree-Aware Maximum Likelihood-Based Gene Family Tree Inference under Gene Duplication, Transfer, and Loss.** *Molecular Biology and Evolution* 37:2763–2774. DOI: [10.1093/molbev/msaa141](https://doi.org/10.1093/molbev/msaa141)

### AleRax
AleRax jointly handles reconciled gene trees and a shared species tree under a probabilistic duplication-transfer-loss model and explicitly addresses gene-tree uncertainty.

- Morel B, Williams TA, Stamatakis A, Szöllősi GJ. 2024. **AleRax: a tool for gene and species tree co-estimation and reconciliation under a probabilistic model of gene duplication, transfer, and loss.** *Bioinformatics* 40:btae162. DOI: [10.1093/bioinformatics/btae162](https://doi.org/10.1093/bioinformatics/btae162)

### Treerecs
Treerecs integrates duplication-loss reconciliation with gene-tree rooting, branch contraction/resolution and correction.

- Comte N et al. 2020. **Treerecs: an integrated phylogenetic tool, from sequences to reconciliations.** *Bioinformatics* 36:4822–4824. DOI: [10.1093/bioinformatics/btaa615](https://doi.org/10.1093/bioinformatics/btaa615)

### Notung
Notung is an established reconciliation package used for duplication/loss and extended event-aware gene-tree/species-tree analysis. GeneSpeciesTreeVerdict does not attempt to replace its event-reconstruction capabilities.

### ASTRAL-Pro / ASTER
ASTRAL-Pro-family methods allow species-tree inference from multi-copy gene-family trees, which is fundamentally different from forcing all families into a single-copy concatenation framework.

- Zhang C, Mirarab S. 2022. **ASTRAL-Pro 2: ultrafast species tree reconstruction from multi-copy gene family trees.** *Bioinformatics* 38:4949–4950. DOI: [10.1093/bioinformatics/btac620](https://doi.org/10.1093/bioinformatics/btac620)

## Role of GeneSpeciesTreeVerdict

GeneSpeciesTreeVerdict is positioned **before** or **between** these specialized methods:

1. validate mappings and copy structure;
2. quantify topological concordance/discordance;
3. identify direct multi-copy evidence;
4. run a transparent, deliberately simple DL diagnostic;
5. explain which follow-up question is appropriate.

The project is therefore a diagnostic and teaching layer rather than a new claim that one reconciliation model solves all causes of discordance.

## Interpretation standard

A method that omits a process can redistribute that discordance into the events it does allow. For example, a duplication-loss-only model can explain discordance using duplication/loss even when the true cause is ILS, HGT, estimation error or an incorrect species tree. GeneSpeciesTreeVerdict therefore labels its reconciliation output as model-conditional and avoids treating event counts as proof.
