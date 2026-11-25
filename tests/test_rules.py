import pytest

from probeharbour.rules import RuleSet, RuleValidationError


@pytest.mark.parametrize("chemistry", ["prime_5k", "xenium_v1"])
def test_packaged_rule_files_have_one_complete_ligbase_table(chemistry):
    rules = RuleSet.from_packaged(chemistry)
    assert len(rules.dimer_scores) == 16


def test_v1_current_rules_recalculate_and_keep_historical_score_separate():
    rules = RuleSet.from_packaged("xenium_v1")
    result = rules.explain({"lig_dinuc": "GA", "tm_full": 76.6, "p3_homodimer_dg": -0.977,
                            "has_simple_repeat": False, "has_5g_repeat": False})
    assert result.base_score == 1.0
    assert result.modifiers == ()
    assert result.recalculated_score == 1.0


def test_duplicate_ligbase_tables_are_rejected():
    payload = {"ligbase_scores": {"AA": 1.0}, "_ligbase": {"AA": 1.0}}
    with pytest.raises(RuleValidationError, match="exactly one"):
        RuleSet.from_mapping(payload, ruleset_id="bad")


def test_user_tier_calibration_is_versioned_and_keeps_qc_modifiers():
    current = RuleSet.from_packaged("xenium_v1")
    calibration = {
        "calibration_id": "partner-cohort-2026-07",
        "tiers": {
            "strong": {"score": 1.0, "dinucleotides": ["AA", "AC", "AG", "AT"]},
            "weak": {"score": 0.2, "dinucleotides": ["CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"]},
        },
    }
    calibrated = current.with_calibration(calibration)
    assert calibrated.ruleset_id == "xenium-v1-current-v1+partner-cohort-2026-07"
    assert calibrated.dimer_scores["AA"] == 1.0
    assert calibrated.dimer_scores["GG"] == 0.2
    assert calibrated.flags == current.flags


def test_user_calibration_must_cover_each_dinucleotide_once():
    current = RuleSet.from_packaged("prime_5k")
    with pytest.raises(RuleValidationError, match="exactly 16"):
        current.with_calibration({"calibration_id": "incomplete", "ligbase_scores": {"AA": 1.0}})


def test_composable_predicates_and_alternative_scoring_strategy():
    all_dimers = [a + b for a in "ACGT" for b in "ACGT"]
    rules = RuleSet.from_mapping({
        "ruleset_id": "custom",
        "ligbase_scores": {dimer: 1.0 for dimer in all_dimers},
        "scoring": {"strategy": "multiplicative_penalty", "factor_per_modifier": 0.8},
        "modifiers": {"problem": {"any": [
            {"type": "comparison", "field": "specificity", "operator": "lt", "value": 90},
            {"all": [{"type": "regex", "field": "sequence", "pattern": "GGGGG"}, {"not": {"type": "boolean", "field": "allowed_g_run"}}]}
        ]}},
    })
    assert rules.explain({"lig_dinuc": "AA", "specificity": 80, "sequence": "ACTG"}).recalculated_score == 0.8
    assert rules.explain({"lig_dinuc": "AA", "specificity": 100, "sequence": "GGGGG", "allowed_g_run": True}).recalculated_score == 1.0
