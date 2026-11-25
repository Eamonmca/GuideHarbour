import pytest

from probeharbour.models import Candidate, Provenance
from probeharbour.thermo import characterize


def test_primer3_characterization_has_parameter_provenance():
    pytest.importorskip("primer3")
    candidate = Candidate(
        candidate_id="thermo-fixture",
        rbd5="ACTGACTGACTGACTG",
        rbd3="TCAGTCAGTCAGTCAG",
        provenance=Provenance(source_id="fixture", dataset_id="test-v1"),
    )
    result = characterize(candidate)
    assert result.parameters_id == "primer3-xenium-schildkraut-v1"
    assert result.tm_full == 73.9
    assert result.tm_full > result.tm_5p
    assert result.candidate_id == candidate.candidate_id
