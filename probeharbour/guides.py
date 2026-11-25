"""Design detection probes across a guide and its flanking target sequence.

This module is maintained directly in Python (it is not an nbdev export).
It keeps ProbeHarbour's target-sense ``rbd5 + rbd3`` convention. The explicitly
named reverse-complement columns describe binding sequences, not complete
vendor oligos, barcodes, scaffolds, or chemical modifications.
"""

from __future__ import annotations

from dataclasses import asdict
from itertools import groupby
import math
from numbers import Integral, Real
from pathlib import Path
from typing import Iterable

import pandas as pd

from .facade import ProbeStore
from .models import Candidate, Provenance
from .rules import RuleSet
from .thermo import ThermoParameters, characterize

__all__ = ["design_guides", "guide_summary", "top_guides"]

_REQUIRED = ("guide_id", "guide_sequence", "scaffold_sequence")
_OPTIONAL_SEQUENCES = ("prefix_sequence", "transcribed_leader")
_RC = str.maketrans("ACGT", "TGCA")
_SCOPE = (
    "Exhaustive ungapped, full-length, target-sense comparisons to the other "
    "supplied targets only; mismatch distance is not a hybridization probability. "
    "Partial/gapped binding, genome/transcriptome off-targets, RNA accessibility "
    "and Cas-protein occupancy are not evaluated."
)
_LIMITS = [
    "Scores reuse inherited ProbeHarbour chemistry heuristics and Primer3 DNA "
    "thermodynamics; they are not validated detection probabilities for sgRNA.",
    "Eligibility means the requested design filters passed; it is not assay validation.",
    "Reverse-complement binding sequences are not order-ready vendor probe constructs.",
    "The supplied prefix is excluded from RNA targets unless explicitly supplied "
    "as a transcribed_leader; no extra initiation G is inferred.",
]

# A stable empty schema allows reporting and exports even when no window fits.
_COLUMNS = [
    "candidate_id", "guide_id", "guide_sequence", "scaffold_sequence",
    "prefix_sequence", "transcribed_leader", "target_sequence", "target_length",
    "guide_start", "guide_end", "scaffold_start", "scaffold_end", "boundary",
    "boundary_position", "start", "end", "ligation_site", "shift_nt",
    "guide_overlap_nt", "boundary_spanning", "on_boundary", "on_junction",
    "junction_spanning", "rbd5", "rbd3", "binder", "arm_length", "lig_dinuc",
    "reverse_complement_binder", "probe_arm_upstream_5to3",
    "probe_arm_downstream_5to3", "tm_5p", "tm_3p", "tm_full", "gc_percent",
    "p3_hairpin_dg", "p3_homodimer_dg", "p3_heterodimer_dg",
    "has_simple_repeat", "has_5g_repeat", "parameters_id", "base_score",
    "modifiers", "recalculated_score", "ruleset_id", "chemistry_qc_pass",
    "nearest_other_guide_id", "nearest_other_mismatches", "nearest_other_start",
    "nearest_other_mismatch_positions", "nearest_other_upstream_mismatches",
    "nearest_other_downstream_mismatches", "other_guides_screened",
    "other_guides_total", "exact_other_match_count", "exact_other_guide_count",
    "upstream_arm_shared_count", "downstream_arm_shared_count",
    "panel_screen_status", "panel_specificity_pass", "reference_screen",
    "target_max_homopolymer", "target_max_g_run", "binder_max_homopolymer",
    "eligible", "exclusion_reasons", "warnings", "rank", "eligible_rank",
]


def _normalize_sequence(value: object, field: str, guide_id: str, *, optional: bool = False) -> str:
    if optional and (value is None or (not isinstance(value, str) and pd.isna(value))):
        return ""
    if not isinstance(value, str):
        raise ValueError(f"{guide_id}: {field} must be a sequence string")
    sequence = "".join(value.split()).upper().replace("U", "T")
    if not sequence and optional:
        return ""
    if not sequence or any(base not in "ACGT" for base in sequence):
        raise ValueError(f"{guide_id}: {field} must contain only A/C/G/T (or U), and cannot be empty")
    return sequence


def _read_guides(guides: pd.DataFrame | str | Path, target_mode: str) -> tuple[list[dict], list[dict], str | None]:
    source_path = None
    if isinstance(guides, pd.DataFrame):
        frame = guides.copy()
    elif isinstance(guides, (str, Path)):
        source_path = str(Path(guides).expanduser())
        # Empty optional CSV cells remain empty strings; numeric IDs/sequences
        # are still rejected by biological validation after loading as strings.
        frame = pd.read_csv(source_path, dtype=str, keep_default_na=False)
    else:
        raise TypeError("guides must be a pandas DataFrame or a CSV path")
    if frame.columns.duplicated().any():
        raise ValueError("Guide input has duplicate column names")
    missing = set(_REQUIRED) - set(frame.columns)
    if missing:
        raise ValueError(f"Guide input is missing required columns: {sorted(missing)}")
    if frame.empty:
        raise ValueError("Guide input must contain at least one guide")
    records, raw_records, seen = [], [], set()
    for index, row in enumerate(frame.to_dict("records"), start=1):
        raw_id = row["guide_id"]
        if not isinstance(raw_id, str) or not raw_id.strip():
            raise ValueError(f"Row {index}: guide_id must be a non-empty string")
        guide_id = raw_id.strip()
        if any(char in guide_id for char in "\r\n\t"):
            raise ValueError(f"Row {index}: guide_id cannot contain tabs or line breaks")
        if guide_id in seen:
            raise ValueError(f"Duplicate guide_id: {guide_id!r}")
        seen.add(guide_id)
        record = {"guide_id": guide_id}
        raw = {"guide_id": raw_id}
        for field in ("guide_sequence", "scaffold_sequence", *_OPTIONAL_SEQUENCES):
            value = row.get(field, "")
            record[field] = _normalize_sequence(value, field, guide_id, optional=field in _OPTIONAL_SEQUENCES)
            raw[field] = value if isinstance(value, str) else None
        if target_mode == "full_construct" and record["transcribed_leader"]:
            raise ValueError(
                f"{guide_id}: full_construct uses prefix_sequence + guide_sequence + scaffold_sequence; "
                "omit transcribed_leader or use target_mode='rna'"
            )
        upstream = record["transcribed_leader"] if target_mode == "rna" else record["prefix_sequence"]
        guide_end = len(upstream) + len(record["guide_sequence"])
        record.update({
            "target_sequence": upstream + record["guide_sequence"] + record["scaffold_sequence"],
            "guide_start": len(upstream), "guide_end": guide_end,
            "scaffold_start": guide_end,
            "scaffold_end": guide_end + len(record["scaffold_sequence"]),
        })
        record["target_length"] = len(record["target_sequence"])
        record["guide_length"] = len(record["guide_sequence"])
        records.append(record)
        raw_records.append(raw)
    return records, raw_records, source_path


def _integer(value: object, name: str, *, minimum: int | None = None, maximum: int | None = None) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise ValueError(f"{name} must be an integer")
    result = int(value)
    if minimum is not None and result < minimum:
        raise ValueError(f"{name} must be >= {minimum}")
    if maximum is not None and result > maximum:
        raise ValueError(f"{name} must be <= {maximum}")
    return result


def _max_run(sequence: str, *, base: str | None = None) -> int:
    return max((sum(1 for _ in group) for char, group in groupby(sequence)
                if base is None or char == base), default=0)


def _validate_parameters(parameters: ThermoParameters) -> None:
    if not isinstance(parameters, ThermoParameters):
        raise TypeError("parameters must be a ThermoParameters instance")
    if not isinstance(parameters.parameters_id, str) or not parameters.parameters_id.strip():
        raise ValueError("ThermoParameters.parameters_id must be a non-empty string")
    for field in ("mv_conc", "dv_conc", "dntp_conc", "dna_conc"):
        value = getattr(parameters, field)
        if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value) or value < 0:
            raise ValueError(f"ThermoParameters.{field} must be finite and nonnegative")
        if field == "dna_conc" and value == 0:
            raise ValueError("ThermoParameters.dna_conc must be greater than zero")
    temperature = parameters.structure_temp_c
    if isinstance(temperature, bool) or not isinstance(temperature, Real) or not math.isfinite(temperature) or not 0 <= temperature <= 100:
        raise ValueError("ThermoParameters.structure_temp_c must be finite and between 0 and 100 Celsius")
    if parameters.salt_corrections_method not in ("schildkraut", "santalucia", "owczarzy"):
        raise ValueError("ThermoParameters.salt_corrections_method must be schildkraut, santalucia, or owczarzy")


def _characterization(candidate: Candidate, parameters: ThermoParameters) -> dict:
    values = characterize(candidate, parameters).as_dict()
    for field in ("tm_5p", "tm_3p", "tm_full", "gc_percent", "p3_hairpin_dg", "p3_homodimer_dg", "p3_heterodimer_dg"):
        value = values.get(field)
        if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value):
            raise RuntimeError(f"Primer3 returned missing or non-finite {field} for {candidate.candidate_id}")
    return values


def _panel_evidence(binder: str, arm_length: int, other_targets: list[dict]) -> dict:
    """Compare sense targets: a complementary probe recognizes their similarity."""
    nearest = None
    exact_sites = 0
    exact_guides = set()
    screened = 0
    upstream_shared = downstream_shared = 0
    for target in other_targets:
        sequence, guide_id = target["target_sequence"], target["guide_id"]
        upstream_shared += int(binder[:arm_length] in sequence)
        downstream_shared += int(binder[arm_length:] in sequence)
        if len(sequence) < len(binder):
            continue
        screened += 1
        for start in range(len(sequence) - len(binder) + 1):
            other = sequence[start:start + len(binder)]
            distance = sum(a != b for a, b in zip(binder, other))
            key = (distance, guide_id, start)
            if nearest is None or key < nearest[:3]:
                nearest = (*key, other)
            if distance == 0:
                exact_sites += 1
                exact_guides.add(guide_id)
    positions = ([index for index, (a, b) in enumerate(zip(binder, nearest[3])) if a != b]
                 if nearest is not None else [])
    return {
        "nearest_other_guide_id": nearest[1] if nearest else None,
        "nearest_other_mismatches": nearest[0] if nearest else None,
        "nearest_other_start": nearest[2] if nearest else None,
        "nearest_other_mismatch_positions": positions,
        "nearest_other_upstream_mismatches": sum(index < arm_length for index in positions) if nearest else None,
        "nearest_other_downstream_mismatches": sum(index >= arm_length for index in positions) if nearest else None,
        "other_guides_screened": screened,
        "other_guides_total": len(other_targets),
        "exact_other_match_count": exact_sites,
        "exact_other_guide_count": len(exact_guides),
        "upstream_arm_shared_count": upstream_shared,
        "downstream_arm_shared_count": downstream_shared,
        "panel_screen_status": ("compared" if screened else
                                "single_guide_only" if not other_targets else "no_full_length_comparisons"),
    }


def design_guides(
    guides: pd.DataFrame | str | Path,
    *,
    chemistry: str = "prime_5k",
    arm_length: int | None = None,
    slides: Iterable[int] = range(-5, 6),
    boundary: str = "guide_scaffold",
    target_mode: str = "rna",
    min_guide_overlap: int = 15,
    min_panel_mismatches: int = 3,
    min_score: float = 0.0,
    require_qc_pass: bool = False,
    parameters: ThermoParameters = ThermoParameters(),
    name: str = "guide_probes",
) -> ProbeStore:
    """Generate, characterize, screen and rank guide-detection probe candidates.

    The default target is ``transcribed_leader + guide_sequence + scaffold_sequence``.
    ``prefix_sequence`` remains input metadata. With ``target_mode='full_construct'``
    the target is exactly ``prefix_sequence + guide_sequence + scaffold_sequence``;
    a nonempty transcribed leader is rejected to avoid ambiguous duplication.

    A candidate consists of two equal, adjacent arms with a cut shifted around
    the named boundary. Only in-bounds windows strictly spanning that boundary
    and covering ``min_guide_overlap`` guide bases are generated. There is no
    padding or silent guide-length correction. All coordinates are 0-based and
    half-open; ``ligation_site`` is a cut between bases.

    Inherited chemistry modifiers reduce ``recalculated_score``. They become
    hard exclusions only with ``require_qc_pass=True``. Exact cross-guide
    matches always fail, including when ``min_panel_mismatches=0``. Mismatch
    thresholds apply only to supplied targets that admit a full-length window.
    A single-guide run explicitly records that no cross-guide comparison exists.

    ``rank`` orders all generated candidates per guide (eligible first);
    ``eligible_rank`` is present only for eligible candidates. Neither field
    implies that a probe has passed experimental or reference-level validation.
    """
    if target_mode not in {"rna", "full_construct"}:
        raise ValueError("target_mode must be 'rna' or 'full_construct'")
    if boundary not in {"guide_scaffold", "prefix_guide", "both"}:
        raise ValueError("boundary must be 'guide_scaffold', 'prefix_guide', or 'both'")
    if not isinstance(name, str) or not name.strip():
        raise ValueError("name must be a non-empty string")
    _validate_parameters(parameters)
    if not isinstance(require_qc_pass, bool):
        raise ValueError("require_qc_pass must be a boolean")
    ruleset = RuleSet.from_packaged(chemistry)
    arm_length = _integer(arm_length if arm_length is not None else
                          (16 if chemistry == "prime_5k" else 20), "arm_length", minimum=2, maximum=30)
    min_guide_overlap = _integer(min_guide_overlap, "min_guide_overlap", minimum=1, maximum=2 * arm_length)
    min_panel_mismatches = _integer(min_panel_mismatches, "min_panel_mismatches", minimum=0, maximum=2 * arm_length)
    if isinstance(min_score, bool) or not isinstance(min_score, Real) or not math.isfinite(min_score) or not 0 <= min_score <= 1:
        raise ValueError("min_score must be a finite number between 0 and 1")
    try:
        offsets = [_integer(value, "slide") for value in slides]
    except TypeError as exc:
        raise ValueError("slides must be a non-empty iterable of integers") from exc
    if not offsets or len(offsets) != len(set(offsets)):
        raise ValueError("slides must be non-empty and contain no duplicate offsets")
    normalized, raw_inputs, source_path = _read_guides(guides, target_mode)
    if boundary in {"prefix_guide", "both"}:
        missing = [record["guide_id"] for record in normalized if not record["guide_start"]]
        if missing:
            field = "transcribed_leader" if target_mode == "rna" else "prefix_sequence"
            raise ValueError(f"prefix_guide boundary requires a non-empty {field}: {missing}")
    rows, generation_counts = [], []
    for guide in normalized:
        target = guide["target_sequence"]
        other_targets = [record for record in normalized if record["guide_id"] != guide["guide_id"]]
        boundaries = (["prefix_guide", "guide_scaffold"] if boundary == "both" else [boundary])
        counts = {"guide_id": guide["guide_id"], "attempted_windows": 0,
                  "skipped_out_of_bounds": 0, "skipped_not_spanning": 0,
                  "skipped_guide_overlap": 0, "generated_candidates": 0}
        for anchor in boundaries:
            position = guide["guide_end"] if anchor == "guide_scaffold" else guide["guide_start"]
            for shift in offsets:
                counts["attempted_windows"] += 1
                cut = position + shift
                start, end = cut - arm_length, cut + arm_length
                if start < 0 or end > len(target):
                    counts["skipped_out_of_bounds"] += 1
                    continue
                if not start < position < end:
                    counts["skipped_not_spanning"] += 1
                    continue
                overlap = max(0, min(end, guide["guide_end"]) - max(start, guide["guide_start"]))
                if overlap < min_guide_overlap:
                    counts["skipped_guide_overlap"] += 1
                    continue
                candidate = Candidate(
                    candidate_id=f"{guide['guide_id']}:{anchor}:slide:{shift:+d}",
                    rbd5=target[start:cut], rbd3=target[cut:end],
                    provenance=Provenance("guide-detection", name, ruleset_id=ruleset.ruleset_id),
                    ligation_site=cut, shift_nt=shift, on_junction=False, junction_spanning=False,
                )
                row = {key: value for key, value in guide.items() if key != "guide_length"}
                row.update({
                    "candidate_id": candidate.candidate_id, "boundary": anchor,
                    "boundary_position": position, "start": start, "end": end,
                    "ligation_site": cut, "shift_nt": shift, "guide_overlap_nt": overlap,
                    "boundary_spanning": True, "on_boundary": cut == position,
                    "on_junction": False, "junction_spanning": False,
                    "rbd5": candidate.rbd5, "rbd3": candidate.rbd3,
                    "binder": candidate.binder, "arm_length": arm_length,
                    "lig_dinuc": candidate.lig_dinuc,
                    "reverse_complement_binder": candidate.binder.translate(_RC)[::-1],
                    "probe_arm_upstream_5to3": candidate.rbd5.translate(_RC)[::-1],
                    "probe_arm_downstream_5to3": candidate.rbd3.translate(_RC)[::-1],
                })
                row.update(_characterization(candidate, parameters))
                explanation = ruleset.explain(row)
                row.update({"base_score": explanation.base_score, "modifiers": list(explanation.modifiers),
                            "recalculated_score": explanation.recalculated_score,
                            "ruleset_id": explanation.ruleset_id,
                            "chemistry_qc_pass": not explanation.modifiers})
                row.update(_panel_evidence(candidate.binder, arm_length, other_targets))
                distance = row["nearest_other_mismatches"]
                reasons = []
                if row["exact_other_match_count"]:
                    reasons.append("exact_match_to_other_guide")
                elif distance is not None and distance < min_panel_mismatches:
                    reasons.append("insufficient_panel_mismatches")
                row["panel_specificity_pass"] = (None if distance is None else not reasons)
                if row["recalculated_score"] < min_score:
                    reasons.append("below_min_score")
                if require_qc_pass and explanation.modifiers:
                    reasons.append("chemistry_qc_modifiers")
                warnings = []
                if row["panel_screen_status"] != "compared":
                    warnings.append(row["panel_screen_status"])
                if row["upstream_arm_shared_count"]:
                    warnings.append("upstream_arm_shared_with_other_guides")
                if row["downstream_arm_shared_count"]:
                    warnings.append("downstream_arm_shared_with_other_guides")
                if guide["guide_length"] != 20:
                    warnings.append("guide_length_is_not_20_nt_preserved_as_supplied")
                target_max_g = _max_run(target, base="G")
                if target_max_g >= 5:
                    warnings.append("target_contains_g_run_of_at_least_5_nt")
                row.update({
                    "reference_screen": "not_run",
                    "target_max_homopolymer": _max_run(target), "target_max_g_run": target_max_g,
                    "binder_max_homopolymer": _max_run(candidate.binder),
                    "eligible": not reasons, "exclusion_reasons": reasons, "warnings": warnings,
                })
                rows.append(row)
                counts["generated_candidates"] += 1
        generation_counts.append(counts)
    frame = pd.DataFrame(rows)
    if not frame.empty:
        # Stable deterministic tie-breaking: chemistry, panel distance, guide
        # coverage, then absolute shift and ID. Missing distance is never infinity.
        frame["_guide_order"] = frame["guide_id"].map({row["guide_id"]: index for index, row in enumerate(normalized)})
        frame["_distance"] = pd.to_numeric(frame["nearest_other_mismatches"], errors="raise").fillna(-1)
        frame["_absolute_shift"] = frame["shift_nt"].abs()
        frame = frame.sort_values(
            ["_guide_order", "eligible", "recalculated_score", "_distance", "guide_overlap_nt", "_absolute_shift", "candidate_id"],
            ascending=[True, False, False, False, False, True, True], kind="stable",
        ).drop(columns=["_guide_order", "_distance", "_absolute_shift"]).reset_index(drop=True)
        frame["rank"] = frame.groupby("guide_id", sort=False).cumcount() + 1
        frame["eligible_rank"] = pd.Series(pd.NA, index=frame.index, dtype="Int64")
        eligible_index = frame.index[frame["eligible"]]
        frame.loc[eligible_index, "eligible_rank"] = frame.loc[eligible_index].groupby("guide_id", sort=False).cumcount() + 1
        frame = frame.reindex(columns=_COLUMNS)
    else:
        frame = pd.DataFrame(columns=_COLUMNS)
    metadata = {
        "guides": normalized,
        "guide_design": {
            "schema_version": "1", "source_path": source_path, "raw_inputs": raw_inputs,
            "config": {"chemistry": chemistry, "arm_length": arm_length, "slides": offsets,
                       "boundary": boundary, "target_mode": target_mode,
                       "min_guide_overlap": min_guide_overlap, "min_panel_mismatches": min_panel_mismatches,
                       "min_score": float(min_score), "require_qc_pass": require_qc_pass},
            "parameters": asdict(parameters), "ruleset_id": ruleset.ruleset_id,
            "generation_counts": generation_counts, "specificity_scope": _SCOPE,
            "reference_screen": "not_run", "rna_accessibility": "not_evaluated",
            "sequence_convention": "rbd5/rbd3/binder are target-sense 5-prime to 3-prime, with U written as T",
            "coordinate_convention": "zero-based half-open; ligation_site is an interbase cut",
            "limits": _LIMITS.copy(),
        },
    }
    return ProbeStore(frame, name=name, metadata=metadata, panel_role="guide_detection")


def top_guides(store: ProbeStore, n: int = 3, *, eligible_only: bool = True) -> pd.DataFrame:
    """Return at most ``n`` distinct binding sequences per guide, preserving rank.

    The default never fills an empty eligible set with failed candidates.
    ``eligible_only=False`` deliberately includes failures with their flags.
    """
    n = _integer(n, "n", minimum=1)
    if not isinstance(eligible_only, bool):
        raise ValueError("eligible_only must be a boolean")
    if "guide_design" not in store.metadata:
        raise ValueError("Expected a ProbeStore returned by design_guides")
    if store.obs.empty:
        return store.obs.copy()
    rows = store.obs.loc[store.obs["eligible"].eq(True)] if eligible_only else store.obs
    rows = rows.sort_values("rank", kind="stable").drop_duplicates(["guide_id", "binder"])
    return rows.groupby("guide_id", sort=False).head(n).copy().reset_index(drop=True)


def guide_summary(store: ProbeStore) -> pd.DataFrame:
    """Summarize every input guide, including guides with no feasible candidates."""
    if "guide_design" not in store.metadata or "guides" not in store.metadata:
        raise ValueError("Expected a ProbeStore returned by design_guides")
    counts = {row["guide_id"]: row for row in store.metadata["guide_design"]["generation_counts"]}
    result = []
    for guide in store.metadata["guides"]:
        rows = store.obs.loc[store.obs["guide_id"].eq(guide["guide_id"])]
        eligible = rows.loc[rows["eligible"].eq(True)].sort_values("rank", kind="stable")
        best = eligible.iloc[0] if len(eligible) else None
        generated = counts.get(guide["guide_id"], {})
        result.append({
            "guide_id": guide["guide_id"], "guide_length": guide["guide_length"],
            "target_length": guide["target_length"], "n_candidates": len(rows), "n_eligible": len(eligible),
            "best_score": float(best["recalculated_score"]) if best is not None else None,
            "best_candidate_id": str(best["candidate_id"]) if best is not None else None,
            "status": "eligible" if len(eligible) else "no_eligible_candidates" if len(rows) else "no_candidates",
            "skipped_out_of_bounds": generated.get("skipped_out_of_bounds", 0),
            "skipped_not_spanning": generated.get("skipped_not_spanning", 0),
            "skipped_guide_overlap": generated.get("skipped_guide_overlap", 0),
        })
    return pd.DataFrame(result)
