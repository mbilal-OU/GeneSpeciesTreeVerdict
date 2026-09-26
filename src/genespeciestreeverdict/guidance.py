from __future__ import annotations

from .models import NextAnalysis, ReconciliationSummary, TopologyMetrics


def recommend_next_analyses(
    *,
    coverage: float,
    min_coverage: float,
    max_copies: int,
    overlap_duplications: int,
    topology: TopologyMetrics,
    reconciliation: ReconciliationSummary,
) -> list[NextAnalysis]:
    """Route a locus to question-driven follow-up analyses.

    Tool names are examples of specialist methods that address the stated question.
    They are not endorsements and are never treated as substitutes for biological
    interpretation or method-specific assumptions.
    """

    recommendations: list[NextAnalysis] = []

    if max_copies > 1 or overlap_duplications > 0:
        recommendations.append(
            NextAnalysis(
                priority="HIGH",
                question="Which homologous copies are orthologous across the sampled species?",
                action=(
                    "Resolve the multi-copy family before conventional single-copy concatenation. "
                    "Inspect the gene tree, mapping, copy distribution, and duplication history."
                ),
                tools=[
                    "GeneRax or AleRax for DTL-aware gene/species-tree reconciliation",
                    "OrthoFinder when sequence-based orthogroup/orthology inference is needed",
                    "ASTRAL-Pro3 when the goal is species-tree inference from multi-copy gene families",
                ],
                rationale=(
                    "Observed multi-copy structure is direct evidence that the complete family cannot be "
                    "treated as one conventional single-copy marker without resolving relationships."
                ),
            )
        )

    if coverage < min_coverage:
        recommendations.append(
            NextAnalysis(
                priority="HIGH",
                question="Is missingness biological, technical, or caused by the marker-definition threshold?",
                action=(
                    "Audit missing taxa, annotation completeness, orthogroup assignment, and the chosen "
                    "core/prevalence threshold before making a marker decision."
                ),
                tools=[
                    "Upstream pangenome/orthology tables",
                    "Genome/annotation completeness QC",
                ],
                rationale="Low taxon coverage can change both marker eligibility and topology comparison.",
            )
        )

    if topology.comparable and topology.rf_distance not in (None, 0):
        recommendations.append(
            NextAnalysis(
                priority="HIGH" if max_copies == 1 else "MEDIUM",
                question="Is the conflicting gene-tree signal itself well supported?",
                action=(
                    "Inspect support on the conflicting splits and repeat gene-tree inference or alignment/model "
                    "QC when support is weak or sensitive to analysis choices."
                ),
                tools=[
                    "IQ-TREE/RAxML-class support analyses",
                    "DiscoVista or quartet/split summaries for discordance visualization",
                ],
                rationale=(
                    "A topological difference should not be interpreted biologically until uncertainty in the "
                    "gene tree has been considered."
                ),
            )
        )
        recommendations.append(
            NextAnalysis(
                priority="MEDIUM",
                question="Is the reference species-tree branch broadly supported across independent loci?",
                action=(
                    "Use the batch reference-split summary and, where appropriate, a species-tree method that "
                    "models gene-tree heterogeneity rather than treating the supplied reference as truth."
                ),
                tools=[
                    "ASTRAL-family methods for coalescent-aware species-tree analysis",
                    "AleRax/SpeciesRax-class approaches when gene-family evolution is central",
                ],
                rationale=(
                    "Repeated conflict across many loci can indicate uncertainty in the reference branch rather "
                    "than a collection of individually 'bad' genes."
                ),
            )
        )

    if reconciliation.attempted and reconciliation.lca_duplications > 0 and max_copies == 1:
        recommendations.append(
            NextAnalysis(
                priority="MEDIUM",
                question="Does a duplication/loss model remain plausible under a specialist reconciliation?",
                action=(
                    "Treat the native DL-LCA result as a screening signal and test the family with a specialist "
                    "reconciliation method if duplication/loss is a biologically relevant hypothesis."
                ),
                tools=[
                    "GeneRax for sequence-aware DTL gene-tree inference/reconciliation",
                    "AleRax for probabilistic gene/species-tree reconciliation from gene-tree distributions",
                    "RANGER-DTL/Notung-class reconciliation workflows where their assumptions fit",
                ],
                rationale=(
                    "A simple DL-LCA reconciliation is intentionally conservative and cannot uniquely identify "
                    "the process responsible for discordance."
                ),
            )
        )

    if topology.comparable and topology.rf_distance not in (None, 0) and max_copies == 1:
        recommendations.append(
            NextAnalysis(
                priority="MEDIUM",
                question="Could the discordance reflect ILS, HGT/recombination, or another process?",
                action=(
                    "Choose process-specific analyses based on the organism, sampling design, and biology; do "
                    "not infer a cause from RF distance alone."
                ),
                tools=[
                    "ASTRAL-family analyses when ILS across many loci is plausible",
                    "DTL-aware reconciliation when horizontal transfer/duplication-loss is plausible",
                    "Recombination-aware locus screening when within-locus mosaic history is plausible",
                ],
                rationale="Single-copy discordance has multiple legitimate biological and technical causes.",
            )
        )

    if (
        coverage >= min_coverage
        and max_copies == 1
        and topology.comparable
        and topology.rf_distance == 0
        and reconciliation.lca_duplications == 0
    ):
        recommendations.append(
            NextAnalysis(
                priority="ROUTINE",
                question="Is the locus otherwise suitable for the intended phylogenomic analysis?",
                action=(
                    "Proceed to normal alignment, trimming, substitution-model, compositional, rate, and other "
                    "dataset-specific QC. GeneSpeciesTreeVerdict does not replace sequence-level QC."
                ),
                tools=["Your standard alignment and phylogenetic inference workflow"],
                rationale=(
                    "The locus is structurally single-copy and topologically concordant, but those properties "
                    "do not test every source of phylogenetic bias."
                ),
            )
        )

    if not recommendations:
        recommendations.append(
            NextAnalysis(
                priority="ROUTINE",
                question="What additional evidence is needed for this locus?",
                action="Review the evidence matrix and resolve any INFO/FLAG rows before making a final choice.",
                tools=["No specialist tool is automatically required"],
                rationale="The available signals do not justify a more specific routing recommendation.",
            )
        )

    order = {"HIGH": 0, "MEDIUM": 1, "ROUTINE": 2}
    return sorted(recommendations, key=lambda item: order.get(item.priority, 99))
