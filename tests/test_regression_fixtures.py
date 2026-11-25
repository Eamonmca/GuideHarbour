import json
from pathlib import Path

from probeharbour.records import audit_candidate


FIXTURE = Path(__file__).parent / "fixtures" / "regression_candidates.json"


def test_three_panel_regression_fixtures_are_present_and_explicit():
    fixtures = json.loads(FIXTURE.read_text())
    sources = {row["source_panel"] for key, row in fixtures.items() if not key.startswith("synthetic_")}
    assert sources == {"Prime 5K mouse", "Xenium V1 mouse", "96-probe validation panel"}
    for key, row in fixtures.items():
        assert row["raw_design_score"] is not None, key
        assert "on_junction" in row, key
        assert "junction_spanning" in row, key


def test_regression_fixture_audits_are_clean():
    fixtures = json.loads(FIXTURE.read_text())
    for key, row in fixtures.items():
        if key.startswith("synthetic_"):
            continue
        assert audit_candidate(row)["ok"], key
