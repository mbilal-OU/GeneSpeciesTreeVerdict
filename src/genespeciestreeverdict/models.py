from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any


@dataclass(slots=True)
class TopologyMetrics:
    """Topology-comparison statistics for a single-copy gene tree."""

    comparable: bool
    rooted: bool
    shared_taxa: list[str] = field(default_factory=list)
    rf_distance: int | None = None
    normalized_rf: float | None = None
    species_split_count: int = 0
    gene_split_count: int = 0
    shared_split_count: int = 0
    species_split_recovery: float | None = None
    gene_split_agreement: float | None = None
    species_splits: list[str] = field(default_factory=list)
    gene_splits: list[str] = field(default_factory=list)
    species_only_splits: list[str] = field(default_factory=list)
    gene_only_splits: list[str] = field(default_factory=list)
    gene_only_split_supports: dict[str, float | None] = field(default_factory=dict)
    note: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(slots=True)
class ReconciliationSummary:
    """A conservative rooted duplication-loss reconciliation summary."""

    attempted: bool
    lca_duplications: int = 0
    inferred_losses: int = 0
    speciations: int = 0
    duplication_events: list[dict[str, Any]] = field(default_factory=list)
    note: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(slots=True)
class EvidenceItem:
    """One transparent line of evidence used in a locus assessment."""

    key: str
    domain: str
    signal: str
    observation: str
    interpretation: str
    limitation: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(slots=True)
class NextAnalysis:
    """A question-driven follow-up analysis suggestion."""

    priority: str
    question: str
    action: str
    tools: list[str] = field(default_factory=list)
    rationale: str = ""

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(slots=True)
class Verdict:
    """Workflow recommendation, not a claim of biological truth."""

    status: str
    topology_signal: str
    primary_concern: str
    reasons: list[str]
    supported_conclusions: list[str]
    unsupported_conclusions: list[str]
    recommended_action: str
    caution: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(slots=True)
class LocusResult:
    family: str
    species_total: int
    gene_leaves: int
    mapped_genes: int
    unique_species: int
    coverage: float
    copy_counts: dict[str, int]
    max_copies: int
    multicopy_species: list[str]
    species_overlap_duplications: int
    species_overlap_events: list[dict[str, Any]]
    topology: TopologyMetrics
    reconciliation: ReconciliationSummary
    evidence: list[EvidenceItem]
    next_analyses: list[NextAnalysis]
    verdict: Verdict
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)
