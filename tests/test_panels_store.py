from pathlib import Path
import os

import pandas as pd
import pytest

from probeharbour.panels import load_panel
from probeharbour.store import RunManifest, RunStore


STORES = Path(os.environ.get("PROBEHARBOUR_TEST_STORES", "internal_tests/stores"))


@pytest.mark.parametrize(
    ("store_name", "expected_rows"),
    [
        ("xenium_prime_5k_mouse", 211),
        ("xenium_v1_mouse_combined", 794),
        ("96_isoform_panel_scorefix_20260710", 96),
    ],
)
def test_load_real_panel_store_without_changing_historic_scores(store_name, expected_rows):
    parquet_path = STORES / store_name / "probes.parquet"
    if not parquet_path.is_file():
        pytest.skip("set PROBEHARBOUR_TEST_STORES to run private panel regressions")
    stored_scores = pd.read_parquet(parquet_path, columns=["raw_design_score"])["raw_design_score"]

    panel = load_panel(parquet_path.parent)

    assert panel.name == store_name
    assert panel.source_path == parquet_path
    assert len(panel.candidates) == expected_rows
    pd.testing.assert_series_equal(panel.candidates["raw_design_score"], stored_scores)
    assert panel.summary().loc[0, "score_column"] == "raw_design_score"


def test_load_panel_reports_missing_parquet(tmp_path):
    with pytest.raises(FileNotFoundError, match="probes.parquet"):
        load_panel(tmp_path)


def test_run_store_round_trip_and_candidate_order(tmp_path):
    store = RunStore(tmp_path / "runs.sqlite")
    store.create_run(RunManifest("run-1", "dataset-1", "prime_5k@1", "calibration-1"))
    store.add_records("run-1", [
        {"candidate_id": "candidate-b", "historic_score": 0.4, "recalculated_score": 0.5},
        {"candidate_id": "candidate-a", "historic_score": 0.7, "recalculated_score": 0.6},
    ])

    records = store.get_records("run-1")

    assert [record["candidate_id"] for record in records] == ["candidate-a", "candidate-b"]
    assert records[1]["historic_score"] == 0.4
    assert records[1]["recalculated_score"] == 0.5
