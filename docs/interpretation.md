# Interpreting gene-tree/species-tree discordance

## First question: is the family multi-copy?
If yes, direct paralogy/duplication structure is plausible and the complete family should not be treated as an ordinary one-sequence-per-species marker without resolution.

## If it is single-copy, do not jump to paralogy
A one-copy gene tree can disagree with a species tree because of:
- incomplete lineage sorting;
- horizontal gene transfer or homologous replacement;
- recombination;
- limited phylogenetic signal;
- alignment problems;
- substitution-model misspecification;
- long-branch attraction or other estimation artifacts;
- incorrect rooting;
- incorrect or weakly supported species-tree branches.

## Use branch support
A highly discordant topology with weak internal support is not equivalent to a strongly supported alternative history. v0.1 records support on conflicting gene splits where available; more complete support-aware decision logic is on the roadmap.

## Treat the species tree as a hypothesis
Batch mode reports how often full-taxon single-copy gene trees recover each reference split. A branch with weak gene-tree support across loci deserves scrutiny even if the original species-tree analysis had high concatenation support.

## Ask a method-matched question
- Want duplication/loss history? Use reconciliation methods.
- Expect HGT in microbes? Use a DTL-aware framework.
- Expect ILS? Use coalescent-aware methods.
- Want multi-copy gene families to contribute to species-tree inference? Consider methods such as ASTRAL-Pro/ASTER.
- Want conventional concatenation? Resolve to defensible orthologous marker sets first.
