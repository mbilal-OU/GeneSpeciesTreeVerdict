# Limitations

GeneSpeciesTreeVerdict v0.3 remains an **alpha evidence-integration and decision-support framework**. The benchmark layer improves transparency and testability, but it does not by itself establish biological validity across all evolutionary regimes.

## Model limitations
- The built-in reconciliation is a rooted least-common-ancestor duplication-loss diagnostic, not a full probabilistic DTL/ILS model.
- Duplication/loss counts are conditional on the supplied trees and their roots.
- HGT, ILS and recombination are not inferred directly.
- RF distance cannot identify the cause of discordance.
- Branch-length information is not used in the current decision engine.
- A `REVIEW` state is intentionally non-causal: it means the locus needs interpretation, not that a specific evolutionary mechanism has been identified.

## Data limitations
- Tree quality depends on upstream alignment and phylogenetic inference.
- Wrong gene→species mappings invalidate interpretation.
- Unsampled/extinct lineages can alter reconciliation interpretations.
- Pseudo-orthology can remain difficult when ancient duplicates have been differentially lost.
- One sampled sequence per species does not prove orthology.

## Benchmark limitations
- The native v0.3 benchmark is a **structural stress test**, not a complete evolutionary simulator.
- Its topology perturbations do not mechanistically simulate ILS, HGT, recombination or sequence evolution.
- The native benchmark can validate software behavior against known copy-number, coverage and topology truth, but it cannot provide process-specific sensitivity/specificity claims.
- Process-level validation should use the external truth-manifest interface with independent simulators such as SimPhy or Zombi.
- Sequence-level validation still requires simulation followed by alignment and gene-tree re-estimation so that inference error is introduced realistically.
- Empirical validation on curated real datasets remains necessary before a stable methodological claim.

## Workflow limitations
- No automatic ortholog-subtree extraction yet.
- No direct Roary/PIRATE/Panaroo/OrthoFinder adapters yet.
- No automated HGT/ILS model selection.
- Specialist reconciliation and coalescent tools are routed as follow-up methods rather than embedded replacements.

## Interpretation rule
Never report a `REVIEW` state, high RF distance, or native DL-duplication count as proof of paralogy without independent support appropriate to the biological question.

Similarly, do not treat a high native benchmark score as proof that GSTV can identify real ILS, HGT, or duplication histories. Benchmark claims must remain matched to the truth model used to generate the validation data.
