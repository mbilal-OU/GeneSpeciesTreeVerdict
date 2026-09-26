"""GeneSpeciesTreeVerdict: evidence-based gene-tree/species-tree decision support."""

from .analysis import analyze_locus
from .models import EvidenceItem, LocusResult, NextAnalysis, TopologyMetrics, Verdict

__all__ = [
    "EvidenceItem",
    "LocusResult",
    "NextAnalysis",
    "TopologyMetrics",
    "Verdict",
    "analyze_locus",
]
__version__ = "0.2.0"
