# Limitations

GeneSpeciesTreeVerdict v0.1 is an alpha diagnostic toolkit.

## Model limitations
- The built-in reconciliation is a rooted least-common-ancestor duplication-loss diagnostic, not a full probabilistic DTL/ILS model.
- Duplication/loss counts are conditional on the supplied trees and their roots.
- HGT, ILS and recombination are not inferred directly.
- RF distance cannot identify the cause of discordance.
- Branch-length information is not used in the decision engine.

## Data limitations
- Tree quality depends on upstream alignment and phylogenetic inference.
- Wrong gene→species mappings invalidate interpretation.
- Unsampled/extinct lineages can alter reconciliation interpretations.
- Pseudo-orthology can remain difficult when ancient duplicates have been differentially lost.

## Workflow limitations
- No automatic ortholog-subtree extraction in v0.1.
- No sequence simulation yet.
- No direct Roary/PIRATE/Panaroo/OrthoFinder adapters yet.
- No automated HGT/ILS model selection.

## Interpretation rule
Never report a `REVIEW` or DL-duplication count as proof of paralogy without independent support appropriate to the biological question.
