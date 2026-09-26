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
    supported: list[str] = []
    unsupported: list[str] = [
        "The supplied species tree is the unquestionable biological truth.",
        "Topology discordance by itself identifies paralogy, ILS, HGT, or any other unique cause.",
    ]

    if topology.comparable:
        if topology.rf_distance == 0:
            topology_signal = "CONCORDANT"
            supported.append("The compared gene-tree and species-tree bipartitions are concordant.")
        else:
            topology_signal = "DISCORDANT"
            supported.append(
                "The gene tree contains topological discordance relative to the reference tree."
            )
    else:
        topology_signal = "NOT_COMPARABLE"

    if max_copies > 1 or overlap_duplications > 0:
        status = "RESOLVE"
        primary_concern = "MULTICOPY_OR_DUPLICATION_SIGNAL"
        reasons.append(
            "The family is multi-copy and/or contains species-overlap duplication candidates."
        )
        supported.append(
            "The complete sampled family is not a conventional single-copy marker in its current form."
        )
        unsupported.append(
            "The current evidence identifies which individual copy is the correct ortholog."
        )
        recommended = (
            "Do not concatenate the full family as a conventional single-copy marker. Resolve "
            "orthology/paralogy first, or use a method designed for multi-copy gene families."
        )
    elif reconciliation.attempted and reconciliation.lca_duplications > 0:
        status = "REVIEW"
        primary_concern = "MODEL_DEPENDENT_DUPLICATION_SIGNAL"
        reasons.append(
            "Rooted DL reconciliation requires one or more duplication events, but this is "
            "model-dependent and can also reflect other sources of discordance."
        )
        supported.append(
            "A duplication-loss history can reconcile the supplied rooted trees under the native DL model."
        )
        unsupported.append("A biological duplication event has been uniquely demonstrated.")
        recommended = (
            "Review this locus before concatenation. Check gene-tree support, rooting, alternative "
            "biological processes, and the reference species-tree hypothesis; use a specialist "
            "reconciliation method when a DTL explanation is a key hypothesis."
        )
    elif coverage < min_coverage:
        status = "REVIEW"
        primary_concern = "LOW_COVERAGE"
        reasons.append(
            f"Taxon coverage ({coverage:.1%}) is below the configured threshold ({min_coverage:.1%})."
        )
        supported.append("The locus does not meet the configured prevalence threshold.")
        unsupported.append(
            "Missingness is necessarily biological rather than technical or annotation-related."
        )
        recommended = (
            "Review missing taxa, annotation/orthogroup assignment, and your core-gene threshold before "
            "including this locus in the target marker set."
        )
    else:
        status = "PASS"
        primary_concern = "NO_STRUCTURAL_RED_FLAG"
        reasons.append(
            "The locus is single-copy across mapped taxa and meets the coverage threshold."
        )
        supported.append(
            "The locus is structurally suitable as a candidate conventional single-copy marker."
        )
        unsupported.append(
            "Single-copy status proves orthology or guarantees that the locus is free of phylogenetic bias."
        )
        recommended = (
            "This locus is structurally suitable as a candidate conventional phylogenomic marker. "
            "Continue normal alignment/model/data-quality checks; topology discordance should be "
            "interpreted rather than automatically filtered."
        )

    if (
        status == "PASS"
        and discordance_review_threshold is not None
        and topology.comparable
        and topology.normalized_rf is not None
        and topology.normalized_rf >= discordance_review_threshold
    ):
        status = "REVIEW"
        primary_concern = "USER_FLAGGED_TOPOLOGICAL_DISCORDANCE"
        reasons.append(
            "Normalized RF distance exceeds the user-selected discordance review threshold."
        )
        recommended = "Review branch support and biological alternatives before deciding whether to retain the locus."

    if topology_signal == "DISCORDANT":
        reasons.append(
            "Gene-tree/species-tree discordance is present; discordance alone is not evidence of paralogy."
        )

    caution = (
        "GeneSpeciesTreeVerdict produces an evidence-based locus-handling assessment relative to the "
        "supplied data and species-tree hypothesis. It does not declare either tree or any inferred "
        "evolutionary cause to be biological truth."
    )
    return Verdict(
        status=status,
        topology_signal=topology_signal,
        primary_concern=primary_concern,
        reasons=reasons,
        supported_conclusions=supported,
        unsupported_conclusions=unsupported,
        recommended_action=recommended,
        caution=caution,
    )
