__version__ = "0.2.0.dev0"
"""Auditable ProbeHarbour candidate records and scoring contracts."""

from .records import audit_candidate
from .rules import RuleSet, ScoreExplanation
from .service import AuditService
from .models import Candidate, Characterization, GenomicCoordinates, Provenance
from .thermo import ThermoParameters, characterize
from .api import load_rules, select
from .facade import ProbeStore, annotate_transcripts, build_store_from_bed, design_gene, design_junction, load_panel, map_probes, read_bed
from .blast_backend import BlastEvidence, BlastHit, LocalBlastBackend, check_specificity
from .panels import Panel
from .panel_compatibility import check_internal_compatibility, check_panel_compatibility, nearest_probe_distances, reverse_complement
from . import plotting as pl
from .guides import design_guides, guide_summary, top_guides
from . import guide_plotting as gpl

__all__ = ["AuditService", "BlastEvidence", "BlastHit", "Candidate", "Characterization", "GenomicCoordinates", "LocalBlastBackend", "Panel", "ProbeStore", "Provenance", "RuleSet", "ScoreExplanation", "ThermoParameters", "annotate_transcripts", "audit_candidate", "build_store_from_bed", "characterize", "check_internal_compatibility", "check_panel_compatibility", "check_specificity", "design_gene", "design_junction", "load_panel", "load_rules", "map_probes", "nearest_probe_distances", "pl", "read_bed", "reverse_complement", "select"]
__all__ += ["design_guides", "guide_summary", "top_guides", "gpl"]
