from __future__ import annotations

from .models import EvidenceItem, ReconciliationSummary, TopologyMetrics


def build_evidence(
    *,
    coverage: float,
    min_coverage: float,
    max_copies: int,
    multicopy_species: list[str],
    overlap_duplications: int,
    topology: TopologyMetrics,
    reconciliation: ReconciliationSummary,
) -> list[EvidenceItem]:
    """Build a transparent evidence matrix for one locus.

    Each row separates the observed signal from its interpretation and limitation.
    This prevents a single metric (for example RF distance) from being treated as a
    cause-specific diagnosis.
    """

    evidence: list[EvidenceItem] = []

    if coverage >= min_coverage:
        coverage_signal = "SUPPORT"
        coverage_interpretation = "Taxon representation meets the configured marker threshold."
    else:
        coverage_signal = "FLAG"
        coverage_interpretation = "Missing taxa may affect marker suitability and tree comparison."
    evidence.append(
        EvidenceItem(
            key="taxon_coverage",
            domain="sampling",
            signal=coverage_signal,
            observation=f"Coverage is {coverage:.1%}; configured threshold is {min_coverage:.1%}.",
            interpretation=coverage_interpretation,
            limitation="Coverage describes prevalence only; it does not establish orthology.",
        )
    )

    if max_copies > 1:
        species_text = ", ".join(multicopy_species) if multicopy_species else "one or more species"
        evidence.append(
            EvidenceItem(
                key="copy_number",
                domain="copy_number",
                signal="STRONG_FLAG",
                observation=(
                    f"Maximum copy number is {max_copies}; multicopy species: {species_text}."
                ),
                interpretation=(
                    "The sampled family is not single-copy and should not be treated as a conventional "
                    "single-copy marker without resolving homolog relationships."
                ),
                limitation=(
                    "Copy number alone does not identify which copies are orthologs, nor the historical "
                    "event that produced them."
                ),
            )
        )
    else:
        evidence.append(
            EvidenceItem(
                key="copy_number",
                domain="copy_number",
                signal="SUPPORT",
                observation="Exactly one mapped sequence is present per sampled species.",
                interpretation="The observed family is single-copy in the sampled genomes.",
                limitation=(
                    "One sequence per species does not guarantee orthology; differential loss can create "
                    "pseudo-orthologous sampling."
                ),
            )
        )

    if overlap_duplications > 0:
        evidence.append(
            EvidenceItem(
                key="species_overlap",
                domain="duplication_screen",
                signal="STRONG_FLAG",
                observation=f"Detected {overlap_duplications} species-overlap duplication candidate(s).",
                interpretation=(
                    "At least one internal gene-tree node has descendant child lineages containing the "
                    "same species, which is direct evidence that the family contains repeated species copies."
                ),
                limitation=(
                    "This topology-local screen is not a complete evolutionary reconciliation and should "
                    "not be used to infer a precise duplication history by itself."
                ),
            )
        )
    else:
        evidence.append(
            EvidenceItem(
                key="species_overlap",
                domain="duplication_screen",
                signal="NEUTRAL",
                observation="No species-overlap duplication candidate was detected.",
                interpretation="No direct repeated-species signal is visible in this screen.",
                limitation="Absence of species overlap does not rule out ancient duplication plus loss.",
            )
        )

    if not topology.comparable:
        evidence.append(
            EvidenceItem(
                key="topology",
                domain="topology",
                signal="INFO",
                observation=topology.note or "Topology comparison was not performed.",
                interpretation="No RF/split-based conclusion is available for this locus.",
                limitation="A non-comparable tree is not evidence for either concordance or discordance.",
            )
        )
    elif topology.rf_distance == 0:
        evidence.append(
            EvidenceItem(
                key="topology",
                domain="topology",
                signal="SUPPORT",
                observation="Gene and reference species trees have RF distance 0 on shared taxa.",
                interpretation="The compared bipartitions are topologically concordant.",
                limitation=(
                    "Topological concordance does not prove orthology or exclude shared systematic error."
                ),
            )
        )
    else:
        nrf = topology.normalized_rf if topology.normalized_rf is not None else 0.0
        recovery = (
            f"{topology.species_split_recovery:.1%}"
            if topology.species_split_recovery is not None
            else "NA"
        )
        evidence.append(
            EvidenceItem(
                key="topology",
                domain="topology",
                signal="FLAG",
                observation=(
                    f"Trees are discordant: RF={topology.rf_distance}, nRF={nrf:.3f}, "
                    f"reference-split recovery={recovery}."
                ),
                interpretation="At least one bipartition differs between the gene and reference trees.",
                limitation=(
                    "Discordance is not cause-specific: duplication/loss, ILS, HGT/recombination, tree "
                    "error, model misspecification, taxon sampling, or species-tree error can contribute."
                ),
            )
        )

    supported_conflicts = sum(
        support is not None for support in topology.gene_only_split_supports.values()
    )
    if topology.comparable and topology.gene_only_splits:
        evidence.append(
            EvidenceItem(
                key="conflict_support_context",
                domain="branch_support",
                signal="INFO",
                observation=(
                    f"{len(topology.gene_only_splits)} conflicting gene-tree split(s) were found; "
                    f"support values are available for {supported_conflicts}."
                ),
                interpretation="Branch support should be inspected before interpreting a topological conflict.",
                limitation=(
                    "Support scales differ among inference methods, so GeneSpeciesTreeVerdict does not impose "
                    "a universal strong-support cutoff."
                ),
            )
        )

    if reconciliation.attempted:
        if reconciliation.lca_duplications > 0:
            evidence.append(
                EvidenceItem(
                    key="dl_reconciliation",
                    domain="reconciliation",
                    signal="FLAG",
                    observation=(
                        f"Rooted DL-LCA reconciliation requires {reconciliation.lca_duplications} "
                        f"duplication(s) and {reconciliation.inferred_losses} inferred loss(es)."
                    ),
                    interpretation=(
                        "A duplication-loss history can explain the supplied rooted trees under this model."
                    ),
                    limitation=(
                        "The count is conditional on the supplied roots, trees, and DL model; ILS, HGT, "
                        "gene-tree error, or species-tree error can be represented as apparent DL events."
                    ),
                )
            )
        else:
            evidence.append(
                EvidenceItem(
                    key="dl_reconciliation",
                    domain="reconciliation",
                    signal="SUPPORT",
                    observation="Rooted DL-LCA reconciliation requires no duplication event.",
                    interpretation="The supplied rooted trees do not require duplication under this DL model.",
                    limitation=(
                        "A zero-duplication DL reconciliation does not prove orthology or exclude processes "
                        "outside the model."
                    ),
                )
            )
    else:
        evidence.append(
            EvidenceItem(
                key="dl_reconciliation",
                domain="reconciliation",
                signal="INFO",
                observation="Duplication-loss reconciliation was not run.",
                interpretation="No model-based DL evidence is available for this locus.",
                limitation="Copy number and topology evidence remain interpretable independently.",
            )
        )

    return evidence
