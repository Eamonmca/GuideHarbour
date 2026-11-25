"""Behavioral checks for guide report failures, orientation and CLI exports."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import pytest

from probeharbour.guide_cli import main
from probeharbour.guide_plotting import guide_overview, plot_guide_candidates, write_guide_report
from probeharbour.guides import design_guides


SCAFFOLD = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAGTGGCACCGAGTCGGTGC"


def _inputs():
    return pd.DataFrame([
        {"guide_id": "one / two", "guide_sequence": "GCTGGCCGAAGCGACGAAGA", "scaffold_sequence": SCAFFOLD, "prefix_sequence": "TGGAAAGGACGAAACACCG"},
        {"guide_id": "one___two", "guide_sequence": "ACTTCTCAACCGTCCCTGCA", "scaffold_sequence": SCAFFOLD, "prefix_sequence": "TGGAAAGGACGAAACACCG"},
    ])


def test_report_displays_sequences_and_escapes_guide_names(tmp_path):
    pytest.importorskip("primer3")
    inputs = _inputs()
    unsafe_id = '<img src=x onerror="alert(1)">'
    inputs.loc[0, "guide_id"] = unsafe_id
    store = design_guides(inputs)
    report = write_guide_report(store, tmp_path / "report.html").read_text()
    assert unsafe_id not in report
    assert "&lt;img src=x onerror=&quot;alert(1)&quot;&gt;" in report
    assert "Reverse complement of whole target window" in report
    row = store.obs.loc[store.obs["eligible"]].iloc[0]
    assert row["reverse_complement_binder"] in report
    assert "Panel mismatches" in report
    assert "transcriptome" in report
    assert "src=\"http" not in report and "href=\"http" not in report
    assert '<option value="guide-0">' in report


def test_report_keeps_guides_when_all_candidates_are_filtered(tmp_path):
    pytest.importorskip("primer3")
    inputs = _inputs()
    inputs.loc[1, "guide_sequence"] = inputs.loc[0, "guide_sequence"]
    store = design_guides(inputs)
    assert len(store.obs) > 0
    assert not store.obs["eligible"].any()
    html = write_guide_report(store, tmp_path / "report.html").read_text()
    assert html.count("No eligible candidate under these settings.") == 2
    assert "excluded candidate(s)" in html
    assert "one / two" in html and "one___two" in html
    figure = guide_overview(store)
    assert any("no eligible candidate" in text.get_text() for text in figure.axes[0].texts)
    plt.close(figure)


def test_empty_candidate_figures_keep_target_geometry(tmp_path):
    store = design_guides(_inputs(), slides=range(100, 103))
    assert store.obs.empty
    html = write_guide_report(store, tmp_path / "empty.html").read_text()
    assert "0 scored" in html
    overview = guide_overview(store)
    detail = plot_guide_candidates(store, "one / two")
    assert len(overview.axes[0].patches) >= 2
    assert any("No eligible candidate" in text.get_text() for text in detail.axes[0].texts)
    with pytest.raises(KeyError, match="Unknown guide_id"):
        plot_guide_candidates(store, "missing")
    plt.close(overview)
    plt.close(detail)


def test_cli_exports_distinct_fasta_ids_and_preserves_inputs(tmp_path):
    pytest.importorskip("primer3")
    inputs = _inputs()
    source = tmp_path / "guides.csv"
    destination = tmp_path / "results"
    inputs.to_csv(source, index=False)
    assert main([str(source), "--output", str(destination), "--top", "2"]) == 0
    expected = {"candidates.csv", "top_candidates.csv", "guide_summary.csv", "targets.csv", "targets.fasta", "report.html", "overview.png", "overview.svg", "manifest.json"}
    assert expected <= {path.name for path in destination.iterdir()}
    targets = pd.read_csv(destination / "targets.csv")
    assert targets["target_id"].nunique() == len(inputs)
    assert targets["guide_id"].tolist() == inputs["guide_id"].tolist()
    fasta_ids = [line[1:] for line in (destination / "targets.fasta").read_text().splitlines() if line.startswith(">")]
    assert fasta_ids == targets["target_id"].tolist()
    assert all(" " not in identifier and "/" not in identifier for identifier in fasta_ids)
    assert all(sequence == guide + SCAFFOLD for sequence, guide in zip(targets["target_sequence"], inputs["guide_sequence"]))
    manifest = json.loads((destination / "manifest.json").read_text())
    assert manifest["n_guides"] == 2
    assert manifest["metadata"]["guide_design"]["config"]["target_mode"] == "rna"
    assert "transcriptome" in manifest["scope"]
    before = (destination / "manifest.json").read_bytes()
    with pytest.raises(SystemExit) as error:
        main([str(source), "--output", str(destination)])
    assert error.value.code == 2
    assert (destination / "manifest.json").read_bytes() == before


def test_cli_rejects_invalid_source_choice_and_slide_bounds(tmp_path):
    with pytest.raises(SystemExit) as error:
        main([])
    assert error.value.code == 2
    with pytest.raises(SystemExit):
        main(["guides.csv", "--example"])
    with pytest.raises(SystemExit):
        main(["--example", "--slide-min", "5", "--slide-max", "-5"])


def test_plot_limits_are_explicit(tmp_path):
    store = design_guides(_inputs(), slides=range(100, 101))
    with pytest.raises(ValueError, match="positive integer"):
        guide_overview(store, n=0)
    with pytest.raises(ValueError, match="positive integer"):
        write_guide_report(store, tmp_path / "report.html", n=-1)
