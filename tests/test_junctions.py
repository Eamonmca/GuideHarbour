import pytest

from probeharbour.junctions import infer_junction_evidence, spacing_evidence


def test_on_junction_and_spanning_are_separate_explicit_facts():
    exact = infer_junction_evidence(binder_start=84, binder_end=116, ligation_site=100, annotated_positions=[100, 120])
    shifted = infer_junction_evidence(binder_start=84, binder_end=116, ligation_site=98, annotated_positions=[100])
    assert exact.on_junction and exact.junction_spanning
    assert not shifted.on_junction and shifted.junction_spanning
    assert shifted.spanned_positions == (100,)


def test_half_open_binder_endpoint_is_not_spanned():
    evidence = infer_junction_evidence(binder_start=84, binder_end=116, ligation_site=100, annotated_positions=[116])
    assert not evidence.junction_spanning


def test_spacing_always_keeps_the_nearest_probe_evidence():
    result = spacing_evidence([
        {"candidate_id": "a", "ligation_site": 100},
        {"candidate_id": "b", "ligation_site": 106},
        {"candidate_id": "c", "ligation_site": 220},
    ], threshold_nt=10)
    assert result["a"].violates_spacing
    assert result["a"].nearest_candidate_id == "b"
    assert result["a"].nearest_distance_nt == 6
    assert not result["c"].violates_spacing


def test_junction_coordinates_must_be_coherent():
    with pytest.raises(ValueError, match="ligation_site"):
        infer_junction_evidence(binder_start=20, binder_end=40, ligation_site=40, annotated_positions=[])
