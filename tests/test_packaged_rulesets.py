import json
from importlib.resources import files

from probeharbour.rules import RuleSet


def test_both_canonical_rulesets_are_packaged_and_valid():
    root = files("probeharbour.rulesets")
    names = {entry.name for entry in root.iterdir() if entry.name.endswith(".json")}
    assert names == {"prime_5k.json", "xenium_v1.json"}
    for name in names:
        rules = RuleSet.from_mapping(json.loads(root.joinpath(name).read_text()))
        assert len(rules.dimer_scores) == 16
