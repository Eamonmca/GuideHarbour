import pytest

from probeharbour.bed import candidates_to_bed
from probeharbour.models import Candidate, GenomicCoordinates, Provenance
from probeharbour.sequence import gc_percent, has_5g_repeat, has_simple_repeat, ligation_dinucleotide


def make_candidate(**overrides):
    values = {
        "candidate_id": "candidate-1",
        "rbd5": "ACTGACTGACTGACTG",
        "rbd3": "TCAGTCAGTCAGTCAG",
        "provenance": Provenance(source_id="fixture", dataset_id="test-v1"),
        "stored_raw_design_score": 0.5,
    }
    values.update(overrides)
    return Candidate(**values)


def test_candidate_sequence_facts_are_derived_not_repeated_fields():
    candidate = make_candidate()
    assert candidate.lig_dinuc == "GT"
    assert candidate.arm_length == 16
    assert candidate.binder == "ACTGACTGACTGACTGTCAGTCAGTCAGTCAG"


def test_sequence_qc_is_deterministic():
    assert ligation_dinucleotide("ACT", "GTA") == "TG"
    assert gc_percent("AACC") == 50.0
    assert has_simple_repeat("ACTACTACT")
    assert has_5g_repeat("ACCCCCC", "ACTG")
    assert not has_5g_repeat("ACGGGGG", "ACTG")


def test_candidate_requires_explicit_valid_sequences():
    with pytest.raises(ValueError, match="A/C/G/T"):
        make_candidate(rbd5="ACNG")


def test_bed12_export_keeps_coordinate_convention_and_score():
    candidate = make_candidate(coordinates=GenomicCoordinates(
        chrom="chr1", start=100, end=132, strand="+",
        rbd5_start=100, rbd5_end=116, rbd3_start=116, rbd3_end=132,
    ))
    assert candidates_to_bed([candidate]) == (
        "chr1\t100\t132\tcandidate-1\t500\t+\t100\t132\t85,168,104\t2\t16,16\t0,16\n"
    )
