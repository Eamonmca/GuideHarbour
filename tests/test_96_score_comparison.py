from pathlib import Path
import os

import pandas as pd
import pytest

from probeharbour.rules import RuleSet

STORE = Path(os.environ.get("PROBEHARBOUR_TEST_STORES", "internal_tests/stores")) / \
    "96_isoform_panel_scorefix_20260710/probes.parquet"


def test_honest_prospective_score_reproduces_clean_legacy_workflow():
    if not STORE.is_file():
        pytest.skip("set PROBEHARBOUR_TEST_STORES to run private score calibration")
    data = pd.read_parquet(STORE)
    rules = RuleSet.from_packaged("prime_5k")
    data["prospective_score"] = data.apply(
        lambda row: rules.explain(row.to_dict()).recalculated_score, axis=1
    )
    data = data.dropna(subset=["score_10x_csv", "prospective_score"])
    difference = data["prospective_score"] - data["score_10x_csv"]

    assert len(data) == 95
    assert data[["score_10x_csv", "prospective_score"]].corr().iloc[0, 1] == pytest.approx(0.9966350624)
    assert difference.abs().mean() == pytest.approx(0.0047873684)
    assert (difference.abs() <= 1e-4).sum() == 38
    assert (difference.abs() > 0.01).sum() == 2
    assert data.loc[difference.abs().idxmax(), "probe_name"] == "Mbnl2_ex6"
    assert difference.abs().max() == pytest.approx(0.2083)
