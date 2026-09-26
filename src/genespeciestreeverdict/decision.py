from __future__ import annotations

from .models import ReconciliationSummary, TopologyMetrics, Verdict


def make_verdict(
    *,
    coverage: float,
    min_coverage: float,
    max_copies: int,
    overlap_duplications: int,
    reconciliation: ReconciliationSummary,
    topology: TopologyMetrics,
    discordance_review_threshold: float | None = None,
) -> Verdict:
    reasons: list[str] = []

    if topology.comparable:
        if topology.rf_distance == 0:
            topology_signal = "CONCORDANT"
        else:
            topology_signal = "DISCORDANT"
    else:
        topology_signal = "NOT_COMPARABLE"

    if max_copies > 1 or overlap_duplications > 0:
        status = "RESOLVE"
        reasons.append(
            "The family is multi-copy and/or contains species-overlap duplication candidates."
        )
        recommended = (
            "Do not concatenate the full family as a conventional single-copy marker. "
            "Resolve orthology/paralogy first, or analyze the family with a paralogy-aware method."
        )
    elif reconciliation.attempted and reconciliation.lca_duplications > 0:
        status = "REVIEW"
        reasons.append(
            "Rooted DL reconciliation requires one or more duplication events, but this is "
            "model-dependent and can also reflect other sources of discordance."
        )
        recommended = (
            "Review this locus before concatenation. Check gene-tree support, rooting, HGT/ILS or "
            "recombination where relevant, and the reference species-tree hypothesis."
        )
    elif coverage < min_coverage:
        status = "REVIEW"
        reasons.append(
            f"Taxon coverage ({coverage:.1%}) is below the configured threshold ({min_coverage:.1%})."
        )
        recommended = (
            "Review missing taxa and your core-gene threshold before including this locus in the "
            "target marker set."
        )
    else:
        status = "PASS"
        reasons.append("The locus is single-copy across mapped taxa and meets the coverage threshold.")
        recommended = (
            "This locus is structurally suitable as a candidate conventional phylogenomic marker. "
            "Topology discordance should still be interpreted rather than automatically filtered."
        )

    if (
        status == "PASS"
        and discordance_review_threshold is not None
        and topology.comparable
        and topology.normalized_rf is not None
        and topology.normalized_rf >= discordance_review_threshold
    ):
        status = "REVIEW"
        reasons.append(
            "Normalized RF distance exceeds the user-selected discordance review threshold."
        )
        recommended = (
            "Review branch support and biological alternatives before deciding whether to retain the locus."
        )

    if topology_signal == "DISCORDANT":
        reasons.append(
            "Gene-tree/species-tree discordance is present; discordance alone is not evidence of paralogy."
        )

    caution = (
        "GeneSpeciesTreeVerdict produces a locus-handling assessment relative to the supplied species "
        "tree. It does not declare either tree to be the biological truth."
    )
    return Verdict(
        status=status,
        topology_signal=topology_signal,
        reasons=reasons,
        recommended_action=recommended,
        caution=caution,
    )
