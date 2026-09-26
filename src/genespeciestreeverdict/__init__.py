"""GeneSpeciesTreeVerdict: diagnostic gene-tree/species-tree comparison."""

from .analysis import analyze_locus
from .models import LocusResult, TopologyMetrics, Verdict

__all__ = ["LocusResult", "TopologyMetrics", "Verdict", "analyze_locus"]
__version__ = "0.1.0"
