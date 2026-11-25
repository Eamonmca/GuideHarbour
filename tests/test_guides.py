"""Behavioral contracts for guide-target design; no biological efficacy claims."""

import json
from dataclasses import replace

import pandas as pd
import pytest

import probeharbour.guides as guide_module
from probeharbour.guides import design_guides, guide_summary, top_guides
from probeharbour.models import Characterization
from probeharbour.thermo import ThermoParameters


SCAFFOLD = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAGTGGCACCGAGTCGGTGC"
PREFIX = "TGGAAAGGACGAAACACCG"
GUIDE = "GCTGGCCGAAGCGACGAAGA"
REAL_CHARACTERIZE = guide_module.characterize


def _frame(*guides, prefix=PREFIX):
    return pd.DataFrame([{"guide_id": f"guide-{i + 1}", "guide_sequence": sequence,
                          "scaffold_sequence": SCAFFOLD, "prefix_sequence": prefix}
                         for i, sequence in enumerate(guides or (GUIDE,))])


@pytest.fixture(autouse=True)
def fake_thermo(monkeypatch):
    def fake(candidate, parameters):
        return Characterization(candidate.candidate_id, 60.0, 60.0, 74.0, 50.0,
                                0.0, 0.0, 0.0, False, False, parameters.parameters_id)
    monkeypatch.setattr(guide_module, "characterize", fake)


def test_rna_target_excludes_prefix_and_keeps_unmodified_input():
    source = _frame()
    original = source.copy(deep=True)
    store = design_guides(source)
    pd.testing.assert_frame_equal(source, original)
    assert store.metadata["guides"][0]["target_sequence"] == GUIDE + SCAFFOLD
    assert store.metadata["guides"][0]["prefix_sequence"] == PREFIX
    assert set(store.obs["guide_start"]) == {0}
    assert store.metadata["guide_design"]["reference_screen"] == "not_run"
    assert store.metadata["guide_design"]["rna_accessibility"] == "not_evaluated"
    json.dumps(store.metadata)  # Persistable provenance, including raw inputs.


def test_sequence_and_physical_binding_arm_orientations_are_explicit():
    row = design_guides(_frame(), slides=[0]).obs.iloc[0]
    reverse_complement = lambda value: value.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    assert row["start"] == 4 and row["ligation_site"] == 20 and row["end"] == 36
    assert row["rbd5"] == GUIDE[-16:] and row["rbd3"] == SCAFFOLD[:16]
    assert row["binder"] == (GUIDE + SCAFFOLD)[4:36]
    assert row["lig_dinuc"] == row["rbd5"][-1] + row["rbd3"][0]
    assert row["probe_arm_upstream_5to3"] == reverse_complement(row["rbd5"])
    assert row["probe_arm_downstream_5to3"] == reverse_complement(row["rbd3"])
    assert row["reverse_complement_binder"] == reverse_complement(row["binder"])
    assert row["reverse_complement_binder"] == row["probe_arm_downstream_5to3"] + row["probe_arm_upstream_5to3"]
    assert row["on_boundary"] and not row["on_junction"] and not row["junction_spanning"]


def test_default_overlap_rejects_scaffold_dominated_windows_without_padding():
    store = design_guides(_frame())
    assert set(store.obs["shift_nt"]) == {-4, -3, -2, -1, 0, 1}
    assert store.obs["guide_overlap_nt"].between(15, 20).all()
    assert (store.obs["start"] >= 0).all() and (store.obs["end"] <= store.obs["target_length"]).all()
    counts = store.metadata["guide_design"]["generation_counts"][0]
    assert counts["attempted_windows"] == 11
    assert counts["skipped_out_of_bounds"] == 1 and counts["skipped_guide_overlap"] == 4
    assert counts["generated_candidates"] == 6


def test_twenty_one_base_guides_are_preserved_and_normalization_is_recorded():
    sequence = "GCAGCTAGTAAAGTCGACACC"
    frame = _frame(sequence.lower().replace("t", "u") + " \n")
    store = design_guides(frame)
    assert store.metadata["guides"][0]["guide_length"] == 21
    assert set(store.obs["guide_sequence"]) == {sequence}
    assert set(store.obs["shift_nt"]) == {-5, -4, -3, -2, -1, 0, 1}
    assert all("guide_length_is_not_20_nt_preserved_as_supplied" in flags for flags in store.obs["warnings"])
    assert store.metadata["guide_design"]["raw_inputs"][0]["guide_sequence"] == frame.iloc[0]["guide_sequence"]


def test_full_construct_and_explicit_rna_leader_are_distinct_modes():
    full = design_guides(_frame(), target_mode="full_construct", boundary="both")
    assert set(full.obs["target_sequence"]) == {PREFIX + GUIDE + SCAFFOLD}
    assert set(full.obs["boundary"]) == {"prefix_guide", "guide_scaffold"}
    assert set(full.obs["guide_start"]) == {19}
    with pytest.raises(ValueError, match="transcribed_leader"):
        design_guides(_frame(), boundary="prefix_guide")
    leader_frame = _frame().assign(transcribed_leader="AACCGGTTAACCGGTT")
    rna = design_guides(leader_frame, boundary="both")
    assert set(rna.obs["target_sequence"]) == {"AACCGGTTAACCGGTT" + GUIDE + SCAFFOLD}
    with pytest.raises(ValueError, match="omit transcribed_leader"):
        design_guides(leader_frame, target_mode="full_construct")


def test_window_must_strictly_span_requested_boundary():
    store = design_guides(_frame(), slides=[16], min_guide_overlap=1)
    assert store.obs.empty
    assert store.metadata["guide_design"]["generation_counts"][0]["skipped_not_spanning"] == 1
    assert guide_summary(store).iloc[0]["status"] == "no_candidates"
    assert top_guides(store).empty
    assert "reverse_complement_binder" in store.obs


def test_v1_arm_length_defaults_to_twenty_and_cannot_extend_before_rna_start():
    store = design_guides(_frame(), chemistry="xenium_v1")
    assert set(store.obs["arm_length"]) == {20}
    assert set(store.obs["shift_nt"]) == {0, 1, 2, 3, 4, 5}
    assert store.metadata["guide_design"]["generation_counts"][0]["skipped_out_of_bounds"] == 5


def test_identical_guides_fail_exact_screen_even_when_mismatch_threshold_is_zero():
    store = design_guides(_frame(GUIDE, GUIDE), min_panel_mismatches=0)
    assert not store.obs["eligible"].any()
    assert (store.obs["nearest_other_mismatches"] == 0).all()
    assert all("exact_match_to_other_guide" in reasons for reasons in store.obs["exclusion_reasons"])
    assert top_guides(store).empty
    assert len(top_guides(store, n=1, eligible_only=False)) == 2
    summary = guide_summary(store)
    assert set(summary["status"]) == {"no_eligible_candidates"}
    assert summary["best_score"].isna().all()
    assert store.obs["eligible_rank"].isna().all()


def test_near_matches_have_per_arm_evidence_and_threshold_excludes_them():
    other = GUIDE[:10] + ("A" if GUIDE[10] != "A" else "C") + GUIDE[11:]
    store = design_guides(_frame(GUIDE, other), slides=[0])
    assert (store.obs["nearest_other_mismatches"] == 1).all()
    assert (store.obs["nearest_other_upstream_mismatches"] == 1).all()
    assert (store.obs["nearest_other_downstream_mismatches"] == 0).all()
    assert all(value == [6] for value in store.obs["nearest_other_mismatch_positions"])
    assert not store.obs["eligible"].any()
    relaxed = design_guides(_frame(GUIDE, other), slides=[0], min_panel_mismatches=1)
    assert relaxed.obs["eligible"].all()
    assert (relaxed.obs["downstream_arm_shared_count"] == 1).all()


def test_panel_comparison_searches_all_offsets_not_only_guide_coordinate():
    frame = _frame(GUIDE, "TACACACACACACACACACA")
    binder = (GUIDE + SCAFFOLD)[4:36]
    frame.loc[1, "scaffold_sequence"] = "ACACAC" + binder + SCAFFOLD
    store = design_guides(frame, slides=[0])
    first = store.obs.loc[store.obs["guide_id"].eq("guide-1")].iloc[0]
    assert first["nearest_other_mismatches"] == 0
    assert first["nearest_other_start"] == 26
    assert first["nearest_other_guide_id"] == "guide-2"
    assert not first["eligible"]


def test_single_guide_is_not_reported_as_reference_validated_or_cross_compared():
    row = design_guides(_frame(), slides=[0]).obs.iloc[0]
    assert row["other_guides_screened"] == 0
    assert row["panel_screen_status"] == "single_guide_only"
    assert row["nearest_other_mismatches"] is None
    assert row["panel_specificity_pass"] is None
    assert row["reference_screen"] == "not_run"
    assert "single_guide_only" in row["warnings"]


def test_chemistry_modifiers_are_penalties_unless_requested_as_hard_filter(monkeypatch):
    def failed_tm(candidate, parameters):
        values = {
            "candidate_id": candidate.candidate_id, "tm_5p": 55, "tm_3p": 55,
            "tm_full": 60, "gc_percent": 50, "p3_hairpin_dg": 0,
            "p3_homodimer_dg": 0, "p3_heterodimer_dg": 0, "has_simple_repeat": False,
            "has_5g_repeat": False, "parameters_id": parameters.parameters_id,
        }
        return Characterization(**values)
    monkeypatch.setattr(guide_module, "characterize", failed_tm)
    soft = design_guides(_frame(), slides=[0])
    assert soft.obs.iloc[0]["eligible"] and not soft.obs.iloc[0]["chemistry_qc_pass"]
    assert "does_not_meet_tm_constraints" in soft.obs.iloc[0]["modifiers"]
    hard = design_guides(_frame(), slides=[0], require_qc_pass=True)
    assert not hard.obs.iloc[0]["eligible"]
    assert "chemistry_qc_modifiers" in hard.obs.iloc[0]["exclusion_reasons"]


def test_low_score_has_no_eligible_fallback_and_summary_keeps_short_guides():
    store = design_guides(_frame(GUIDE, "ACT"), min_score=1.0, slides=[0])
    assert len(store.obs) == 1 and not store.obs.iloc[0]["eligible"]
    assert top_guides(store).empty
    summary = guide_summary(store)
    assert list(summary["status"]) == ["no_eligible_candidates", "no_candidates"]


def test_g_run_advisory_does_not_silently_change_inherited_score():
    store = design_guides(_frame("TGCGATGGGGGGGTGGGTAGC"), slides=[0])
    row = store.obs.iloc[0]
    assert row["target_max_g_run"] == 7 and row["target_max_homopolymer"] == 7
    assert "target_contains_g_run_of_at_least_5_nt" in row["warnings"]
    assert not row["has_5g_repeat"]  # Stubbed inherited metric remains independent.


@pytest.mark.parametrize("column,value", [
    ("guide_id", ""), ("guide_id", None), ("guide_id", 12),
    ("guide_sequence", None), ("guide_sequence", float("nan")),
    ("guide_sequence", "ACTN"), ("guide_sequence", 1234),
    ("scaffold_sequence", ""), ("prefix_sequence", 123),
])
def test_invalid_input_is_rejected(column, value):
    frame = _frame().astype(object)
    frame.loc[0, column] = value
    with pytest.raises(ValueError):
        design_guides(frame)


def test_missing_duplicate_and_empty_inputs_are_rejected():
    with pytest.raises(ValueError, match="missing required"):
        design_guides(_frame().drop(columns="scaffold_sequence"))
    with pytest.raises(ValueError, match="at least one"):
        design_guides(_frame().iloc[:0])
    with pytest.raises(ValueError, match="Duplicate guide_id"):
        design_guides(pd.concat([_frame(), _frame()], ignore_index=True))
    with pytest.raises(ValueError, match="Duplicate guide_id"):
        design_guides(pd.concat([_frame(), _frame().assign(guide_id=" guide-1 ")], ignore_index=True))


@pytest.mark.parametrize("kwargs", [
    {"slides": []}, {"slides": [0, 0]}, {"slides": [1.2]}, {"slides": [True]},
    {"arm_length": 0}, {"arm_length": 31}, {"arm_length": 16.5},
    {"min_guide_overlap": 0}, {"min_guide_overlap": 33},
    {"min_panel_mismatches": -1}, {"min_panel_mismatches": 33},
    {"min_score": float("nan")}, {"min_score": 2},
    {"boundary": "exon"}, {"target_mode": "dna"}, {"require_qc_pass": "yes"},
])
def test_invalid_configuration_is_rejected(kwargs):
    with pytest.raises(ValueError):
        design_guides(_frame(), **kwargs)


def test_csv_input_allows_blank_optional_cells_and_preserves_id(tmp_path):
    path = tmp_path / "guides.csv"
    _frame().assign(guide_id="LINC00152/ CYTOR-sg3", transcribed_leader="").to_csv(path, index=False)
    store = design_guides(path)
    assert store.metadata["guides"][0]["guide_id"] == "LINC00152/ CYTOR-sg3"
    assert store.metadata["guide_design"]["source_path"] == str(path)


@pytest.mark.parametrize("change", [
    {"mv_conc": -1}, {"dv_conc": float("nan")}, {"dntp_conc": float("inf")},
    {"dna_conc": 0}, {"mv_conc": "80"}, {"dv_conc": True},
    {"structure_temp_c": -1}, {"structure_temp_c": 101}, {"structure_temp_c": float("nan")},
    {"parameters_id": " "}, {"salt_corrections_method": "unknown"},
])
def test_invalid_thermodynamic_parameters_are_rejected(change):
    with pytest.raises(ValueError, match="ThermoParameters"):
        design_guides(_frame(), parameters=replace(ThermoParameters(), **change))


def test_nonfinite_characterization_never_yields_a_score(monkeypatch):
    def invalid(candidate, parameters):
        return Characterization(candidate.candidate_id, 60.0, 60.0, float("nan"), 50.0,
                                0.0, 0.0, 0.0, False, False, parameters.parameters_id)
    monkeypatch.setattr(guide_module, "characterize", invalid)
    with pytest.raises(RuntimeError, match="non-finite tm_full"):
        design_guides(_frame())


def test_real_primer3_properties_and_probestore_roundtrip(monkeypatch, tmp_path):
    pytest.importorskip("primer3")
    pytest.importorskip("pyarrow")
    monkeypatch.setattr(guide_module, "characterize", REAL_CHARACTERIZE)
    store = design_guides(_frame(), slides=[0])
    row = store.obs.iloc[0]
    assert isinstance(row["tm_full"], float) and row["tm_full"] > row["tm_5p"]
    assert row["parameters_id"] == "primer3-xenium-schildkraut-v1"
    directory = store.save(tmp_path / "store")
    loaded = guide_module.ProbeStore.load(directory)
    assert loaded.metadata["guides"] == store.metadata["guides"]
    assert guide_summary(loaded).iloc[0]["n_candidates"] == 1
    assert len(top_guides(loaded)) == 1
    assert loaded.to_candidates()[0].binder == row["binder"]


def test_persisted_panel_keeps_all_guide_summaries_and_screening_evidence(tmp_path):
    pytest.importorskip("pyarrow")
    from probeharbour.facade import load_panel
    store = design_guides(_frame(GUIDE, "GCAGCTAGTAAAGTCGACACC", "ACT"))
    loaded = load_panel(store.save(tmp_path / "panel"))
    pd.testing.assert_frame_equal(guide_summary(store), guide_summary(loaded))
    assert list(top_guides(loaded)["candidate_id"]) == list(top_guides(store)["candidate_id"])
    first = loaded.obs.iloc[0]
    assert list(first["nearest_other_mismatch_positions"]) == store.obs.iloc[0]["nearest_other_mismatch_positions"]
    assert loaded.metadata["guide_design"]["parameters"] == store.metadata["guide_design"]["parameters"]
