from io import BytesIO, StringIO
import zipfile

import pandas as pd

from probeharbour import ProbeStore
from probeharbour import _dashboard_demo_launcher, _dashboard_launcher
from probeharbour.dashboard import (
    add_to_basket,
    candidate_review_table,
    export_panel_bundle,
    junction_summary,
    load_cart_csv,
    reannotate_junction_summary,
    remove_from_basket,
    selected_panel,
)


def candidate_frame():
    return pd.DataFrame({
        "candidate_id": ["c1", "c2", "c3"],
        "probe_name": ["p1", "p2", "p3"], "gene_name": ["GeneA"] * 3,
        "transcript_id": ["tx1"] * 3, "junction_id": ["J1", "J1", "J2"],
        "rbd5": ["AACCGTAACCGTAACC", "AACCGTAACCGTAACG", "TTAACCGGTTAACCGG"],
        "rbd3": ["CCGTTGCACCGTTGCA", "ACGTTGCACCGTTGCA", "GGCCTTAAGGCCTTAA"],
        "shift_nt": [0, 1, 0], "lig_dinuc": ["CC", "GA", "GG"],
        "probeharbour_score": [.4, .3, .8],
        "probeharbour_score_modifiers": [[], ["has_simple_repeat"], []],
        "junction_specificity_class": ["specific", "specific", "partially specific"],
        "junction_supporting_transcripts": [["tx1"], ["tx1"], ["tx1"]],
        "junction_transcript_count": [1, 1, 2], "gene_transcript_count": [3, 3, 3],
        "target_coverage": [33.3, 33.3, 66.7],
        "coordinates": [
            {"chrom": "1", "start": 100, "end": 308, "strand": "+",
             "rbd5_start": 100, "rbd5_end": 116, "rbd3_start": 292, "rbd3_end": 308},
            {"chrom": "1", "start": 101, "end": 309, "strand": "+",
             "rbd5_start": 101, "rbd5_end": 117, "rbd3_start": 293, "rbd3_end": 309},
            {"chrom": "1", "start": 500, "end": 708, "strand": "+",
             "rbd5_start": 500, "rbd5_end": 516, "rbd3_start": 692, "rbd3_end": 708},
        ],
    })


def test_installed_dashboard_commands_launch_the_current_workflow():
    assert _dashboard_launcher.APP_FILENAME == "dashboard_demo.py"
    assert _dashboard_demo_launcher.APP_FILENAME == "dashboard_demo.py"


def test_dashboard_summary_and_review_defaults_are_deterministic():
    store = ProbeStore(candidate_frame())

    junctions = junction_summary(store)
    review = candidate_review_table(store)

    assert junctions.set_index("junction_id").loc["J1", "candidate_count"] == 2
    assert junctions["junction_label"].tolist() == ["J1", "J2"]
    assert junctions.set_index("junction_id").loc["J1", "supporting_transcripts"] == ["tx1"]
    assert review.loc[review["candidate_id"].eq("c1"), "include"].item()
    assert not review.loc[review["candidate_id"].eq("c2"), "include"].item()
    assert review.loc[review["candidate_id"].eq("c3"), "include"].item()


def test_junction_summary_can_be_reannotated_against_visible_transcripts():
    summary = pd.DataFrame({
        "junction_label": ["J1", "J2", "J3"],
        "junction_id": ["j1", "j2", "j3"],
        "junction_class": ["partially specific", "specific", "specific"],
        "supporting_transcripts": [["tx1", "tx2"], ["tx2"], ["tx3"]],
        "transcripts_with_junction": [2, 1, 1],
        "annotated_transcripts": [3, 3, 3],
        "transcript_coverage_percent": [66.7, 33.3, 33.3],
    })

    recalculated = reannotate_junction_summary(summary, ["tx1", "tx2"])
    indexed = recalculated.set_index("junction_id")

    assert list(recalculated["junction_id"]) == ["j1", "j2"]
    assert indexed.loc["j1", "junction_class"] == "constitutive"
    assert indexed.loc["j2", "junction_class"] == "specific"
    assert indexed.loc["j1", "transcripts_with_junction"] == 2
    assert indexed.loc["j1", "annotated_transcripts"] == 2
    assert indexed.loc["j1", "transcript_coverage_percent"] == 100


def test_confirmed_dashboard_panel_exports_portable_bundle():
    panel = selected_panel(ProbeStore(candidate_frame()), ["c1", "c3"], name="GeneA panel")

    payload = export_panel_bundle(panel)

    with zipfile.ZipFile(BytesIO(payload)) as archive:
        names = set(archive.namelist())
    assert "GeneA_panel/GeneA_panel.csv" in names
    assert "GeneA_panel/GeneA_panel.bed" in names
    assert "GeneA_panel/GeneA_panel.fasta" in names
    assert "GeneA_panel/probestore/probes.parquet" in names
    assert "GeneA_panel/manifest.json" in names


def test_dashboard_basket_accumulates_deduplicates_and_removes_candidates():
    store = ProbeStore(candidate_frame())
    first = selected_panel(store, ["c1"], name="working panel")
    second = selected_panel(store, ["c1", "c3"], name="working panel")
    first.obs["panel_record_id"] = "local:0"
    second.obs["panel_record_id"] = ["local:0", "local:1"]

    basket = add_to_basket(None, first, name="working panel")
    basket = add_to_basket(basket, second, name="working panel")

    assert basket.obs["candidate_id"].tolist() == ["c1", "c3"]
    assert basket.obs["panel_role"].eq("add_on").all()
    assert "panel_record_id" not in basket.obs
    reduced = remove_from_basket(basket, ["c1"])
    assert reduced.obs["candidate_id"].tolist() == ["c3"]


def test_dashboard_cart_csv_roundtrip_restores_structured_evidence():
    source = candidate_frame().iloc[:2].copy()
    source["panel_compatibility_reasons"] = [["duplicate_arm"], []]
    csv_buffer = StringIO()
    source.to_csv(csv_buffer, index=False)
    csv_buffer.seek(0)

    cart = load_cart_csv(csv_buffer, name="uploaded cart")

    assert cart.name == "uploaded cart"
    assert isinstance(cart.obs.loc[0, "coordinates"], dict)
    assert cart.obs.loc[0, "panel_compatibility_reasons"] == ["duplicate_arm"]
    assert cart.obs["panel_role"].eq("add_on").all()
    assert "panel_record_id" not in cart.obs


def test_dashboard_loads_sequence_only_existing_addon_csv_schema():
    csv_buffer = StringIO(
        "name,rbdm5,rbdm3\n"
        "Snap25_ex5b,AACCGTAACCGTAACC,CCGTTGCACCGTTGCA\n"
    )

    cart = load_cart_csv(csv_buffer, name="existing add-on")

    assert cart.obs.loc[0, "probe_name"] == "Snap25_ex5b"
    assert cart.obs.loc[0, "candidate_id"] == "Snap25_ex5b"
    assert cart.obs.loc[0, "gene_name"] == "Snap25"
    assert cart.obs.loc[0, "rbd5"] == "AACCGTAACCGTAACC"
    assert cart.obs.loc[0, "rbd3"] == "CCGTTGCACCGTTGCA"
