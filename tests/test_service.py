from pathlib import Path

from probeharbour.rules import RuleSet
from probeharbour.service import AuditService


def test_mcp_facing_service_returns_json_safe_score_explanation():
    import json
    rules = RuleSet.from_packaged("xenium_v1")
    candidates = json.loads((Path(__file__).parent / "fixtures" / "regression_candidates.json").read_text())
    service = AuditService(rules, candidates)
    result = service.explain_score("v1_mouse_panel_Aldh1a2_cw379_ps1")
    assert result == {
        "ruleset_id": "xenium-v1-current-v1",
        "ligation_dinucleotide": "GA",
        "base_score": 1.0,
        "modifiers": [],
        "recalculated_score": 1.0,
        "stored_raw_design_score": 0.5,
    }
