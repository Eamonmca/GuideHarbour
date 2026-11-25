import json
from pathlib import Path

from probeharbour.records import audit_candidate


def test_spacing_claim_without_neighbour_evidence_is_not_accepted():
    fixtures = json.loads((Path(__file__).parent / "fixtures" / "regression_candidates.json").read_text())
    result = audit_candidate(fixtures["synthetic_incomplete_spacing_evidence"])
    assert not result["ok"]
    assert result["issues"][0]["code"] == "spacing_evidence_incomplete"
