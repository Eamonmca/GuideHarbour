import matplotlib
matplotlib.use("Agg")

import pandas as pd
import pytest

import probeharbour as ps


class FakeTranscriptData:
    def get_transcripts_by_gene_name(self, gene):
        return ["tx1"] if gene == "GeneA" else []

    def get_transcripts_by_gene_id(self, gene):
        return []

    def get_exon_coords_and_strand(self, transcript_id):
        return [[108, 100], [308, 300]], 1

    def get_transcript_info(self, transcript_id):
        return {"transcript_name": "GeneA-201", "chromosome": "1"}


class TwoTranscriptData(FakeTranscriptData):
    def get_transcripts_by_gene_name(self, gene):
        return ["tx1", "tx2"] if gene == "GeneA" else []

    def get_transcript_info(self, transcript_id):
        return {"transcript_name": transcript_id, "chromosome": "1"}


class PartlyMalformedTranscriptData(FakeTranscriptData):
    def get_exon_coords_and_strand(self, transcript_id):
        return [[float("nan"), 150], [108, 100], [float("inf"), 250], [308, 300]], 1


def test_scanpy_like_plotting_namespace_returns_figures():
    frame = pd.DataFrame({
        "probe_name": ["p1", "p2"], "gene_name": ["GeneA", "GeneA"],
        "rbd5": ["AACCGTAA", "GGCCAATT"], "rbd3": ["CCGTTGCA", "TTAACCGG"],
        "junction_spanning": [True, False],
        "tm_5p": [55.0, 56.0], "tm_3p": [56.0, 55.5], "tm_full": [74.0, 75.0],
        "p3_hairpin_dg": [-1.0, -2.0], "p3_homodimer_dg": [-2.0, -3.0],
        "has_simple_repeat": [False, True], "has_5g_repeat": [False, False],
        "external_score": [.5, .7], "probeharbour_score": [.52, .69],
        "rbd5_chrom_start": [100, 200], "rbd5_chrom_end": [108, 208],
        "rbd3_chrom_start": [300, 208], "rbd3_chrom_end": [308, 216],
    })
    store = ps.ProbeStore(frame, name="demo")

    figures = [
        ps.pl.panel_overview(store),
        ps.pl.thermodynamics(store),
        ps.pl.probe_map(store, gene="GeneA"),
        ps.pl.gene_map("GeneA", FakeTranscriptData(), store),
        ps.pl.score_comparison(store, "external_score", "probeharbour_score"),
    ]

    assert all(figure.axes for figure in figures)


def test_plotting_handles_empty_data():
    figure = ps.pl.panel_overview(pd.DataFrame())
    assert len(figure.axes) == 1


def test_gene_map_focuses_on_mapped_probes_and_supporting_transcripts():
    probes = pd.DataFrame({
        "probe_name": ["event_a"], "gene_name": ["GeneA"], "transcript_id": ["tx1"],
        "rbd5_chrom_start": [100], "rbd5_chrom_end": [108],
        "rbd3_chrom_start": [300], "rbd3_chrom_end": [308],
    })

    figure = ps.pl.gene_map("GeneA", TwoTranscriptData(), probes, show_overview=False)
    ax = figure.axes[-1]

    assert len(ax.texts) == 2  # one transcript label and one probe annotation
    assert ax.get_xlim()[1] - ax.get_xlim()[0] <= 1_200
    assert "e" not in ax.get_xticklabels()[0].get_text().lower()


def test_gene_map_shows_bed_and_nested_designed_coordinates_together():
    probes = pd.DataFrame([
        {
            "probe_name": "existing_1", "gene_name": "GeneA", "panel_role": "existing",
            "chrom_start": 120, "chrom_end": 152,
            "bed_block_starts": "0", "bed_block_sizes": "32",
        },
        {
            "probe_name": "existing_junction", "gene_name": "GeneA",
            "panel_role": "existing", "junction_spanning": True,
            "chrom_start": 100, "chrom_end": 308,
            "bed_block_starts": "0,200", "bed_block_sizes": "8,8",
        },
        {
            # Designed candidates and BED panels may use different gene fields
            # after merging; both must remain visible in the combined map.
            "candidate_id": "designed_1", "gene": "GeneA", "panel_role": "add_on",
            "transcript_id": "tx1",
            "coordinates": {
                "chrom": "1", "start": 100, "end": 308, "strand": "+",
                "rbd5_start": 100, "rbd5_end": 108,
                "rbd3_start": 300, "rbd3_end": 308,
            },
        },
        {
            "probe_name": "malformed_missing_geometry", "gene_name": "GeneA",
            "panel_role": "existing", "chrom_start": float("inf"), "chrom_end": float("inf"),
        },
    ])

    figure = ps.pl.gene_map("GeneA", FakeTranscriptData(), probes, show_overview=False)
    detail = figure.axes[-1]
    labels = [text.get_text() for text in detail.texts]

    assert "10x · existing_1" in labels
    assert "10x · existing_junction" in labels
    assert "designed_1 · add-on" in labels
    coloured_segments = {
        line.get_color(): tuple(line.get_xdata())
        for line in detail.lines
        if line.get_color() in {"#2d5fa6", "#c05820", "#7c3aed"}
    }
    assert coloured_segments["#2d5fa6"] == (100, 108)
    assert coloured_segments["#c05820"] == (300, 308)
    purple_segments = [
        tuple(line.get_xdata()) for line in detail.lines
        if line.get_color() == "#7c3aed"
    ]
    assert (120.0, 152.0) in purple_segments
    assert (100.0, 108.0) in purple_segments
    assert (300.0, 308.0) in purple_segments
    from matplotlib.patches import PathPatch
    assert any(
        isinstance(patch, PathPatch)
        and patch.get_edgecolor()[:3] == matplotlib.colors.to_rgb("#7c3aed")
        for patch in detail.patches
    )
    legend_labels = [
        text.get_text() for text in figure.axes[0].get_legend().get_texts()
    ]
    assert "existing 10x exon probe" in legend_labels
    assert "existing 10x junction probe" in legend_labels
    import math
    assert all(math.isfinite(value) for value in detail.get_xlim())

    full_figure = ps.pl.gene_map(
        "GeneA", FakeTranscriptData(), probes, label_mode="full", show_overview=False
    )
    full_labels = [text.get_text() for text in full_figure.axes[-1].texts]
    assert "existing_1 [existing]" in full_labels
    assert "designed_1 [add_on]" in full_labels


def test_gene_map_ignores_nonfinite_exon_geometry():
    probes = pd.DataFrame({
        "probe_name": ["event_a"], "gene_name": ["GeneA"], "transcript_id": ["tx1"],
        "rbd5_chrom_start": [100], "rbd5_chrom_end": [108],
        "rbd3_chrom_start": [300], "rbd3_chrom_end": [308],
    })

    figure = ps.pl.gene_map(
        "GeneA", PartlyMalformedTranscriptData(), probes, show_overview=False
    )

    import math
    assert all(math.isfinite(value) for value in figure.axes[-1].get_xlim())


def test_gene_map_links_full_locus_overview_to_detail():
    probes = pd.DataFrame({
        "probe_name": ["event_a"], "gene_name": ["GeneA"], "transcript_id": ["tx1"],
        "rbd5_chrom_start": [100], "rbd5_chrom_end": [108],
        "rbd3_chrom_start": [300], "rbd3_chrom_end": [308],
    })

    figure = ps.pl.gene_map("GeneA", TwoTranscriptData(), probes)

    assert len(figure.axes) == 3
    overview, legend, detail = figure.axes
    assert not legend.axison
    legend_labels = [text.get_text() for text in legend.get_legend().get_texts()]
    assert "custom/designed junction arc" in legend_labels
    assert {"RBD5", "RBD3"} <= set(legend_labels)
    assert overview.get_xlim() == (100.0, 308.0)
    assert detail.get_xlim()[0] <= 100
    assert detail.get_xlim()[1] >= 308
    assert "1 probe-supporting transcript" in overview.get_title(loc="left")
    assert "same 1 probe-supporting transcript" in detail.get_title(loc="left")
    assert len(overview.texts) == 1
    assert "5′ → 3′ · positive strand" in overview.texts[0].get_text()
    from matplotlib.patches import PathPatch
    assert any(isinstance(patch, PathPatch) for patch in overview.patches)

    compressed = ps.pl.gene_map(
        "GeneA", TwoTranscriptData(), probes, overview_scale="compressed"
    )
    compressed_overview = compressed.axes[0]
    assert "introns compressed" in compressed_overview.get_title(loc="left")
    assert not compressed_overview.get_xticks().size
    assert compressed_overview.get_xlabel() == ""


def test_gene_map_can_write_vector_svg(tmp_path):
    probes = pd.DataFrame({
        "probe_name": ["event_a"], "gene_name": ["GeneA"], "transcript_id": ["tx1"],
        "rbd5_chrom_start": [100], "rbd5_chrom_end": [108],
        "rbd3_chrom_start": [300], "rbd3_chrom_end": [308],
    })
    destination = tmp_path / "figures" / "gene_map.svg"

    figure = ps.pl.gene_map(
        "GeneA", FakeTranscriptData(), probes,
        save=destination, dpi=600, transparent=True,
    )

    assert figure.axes
    assert destination.is_file()
    assert "<svg" in destination.read_text()[:500]


def test_panel_conflict_alignment_reports_exact_complementary_segment(tmp_path):
    query = pd.Series({
        "probe_name": "Clta_ex6",
        "rbd5": "GGCAAACAACAGCACA",
        "rbd3": "AACATAAACCATCCTT",
    })
    reference = pd.Series({
        "probe_name": "panel_Dnajc5_cw13378_ps1",
        "rbd5": "CTGTTTGTGCTGTTGC",
        "rbd3": "TTTAACTGCTGCTGTG",
    })
    destination = tmp_path / "conflict.svg"

    figure = ps.pl.panel_conflict_alignment(query, reference, save=destination)

    evidence = figure.probeharbour_alignment
    assert evidence["longest_complementary_run"] == 14
    assert evidence["match_sequence"] == "CAACAGCACAAACA"
    assert evidence["query_start"] == 6
    assert evidence["reference_rc_start"] == 17
    assert evidence["threshold"] == 14
    assert evidence["query_ligation_dinucleotide"] == "AA"
    assert evidence["reference_ligation_dinucleotide"] == "CT"
    assert evidence["query_ligation_boundary"] == 16
    assert evidence["reference_rc_ligation_boundary"] == 16
    assert not evidence["pass"]
    assert destination.is_file()


def test_junction_map_links_short_labels_classes_and_selection():
    junctions = pd.DataFrame({
        "junction_label": ["J1"], "junction_id": ["1:108-300:+"],
        "junction_class": ["specific"],
    })

    figure = ps.pl.junction_map(
        "GeneA", FakeTranscriptData(), junctions,
        selected_junctions=["1:108-300:+"], overview_scale="genomic",
    )

    ax = figure.axes[0]
    assert any(text.get_text() == "J1" for text in ax.texts)
    from matplotlib.patches import PathPatch
    arcs = [patch for patch in ax.patches if isinstance(patch, PathPatch)]
    assert len(arcs) == 1
    assert arcs[0].get_edgecolor()[:3] == matplotlib.colors.to_rgb("#c026d3")
    assert len(ax.collections) == 0
    selected_label = next(text for text in ax.texts if text.get_text() == "J1")
    assert selected_label.get_color() == "#c026d3"
    assert "selected" in [
        text.get_text() for text in ax.get_legend().get_texts()
    ]


def test_junction_map_can_filter_transcripts_and_renumber_display_labels():
    junctions = pd.DataFrame({
        "junction_label": ["J7"], "junction_id": ["1:108-300:+"],
        "junction_class": ["specific"],
    })

    figure = ps.pl.junction_map(
        "GeneA", TwoTranscriptData(), junctions,
        transcript_ids=["tx2"], renumber_junctions=True,
        overview_scale="genomic",
    )

    labels = [text.get_text() for text in figure.axes[0].texts]
    assert any(value.startswith("tx2") for value in labels)
    assert not any(value.startswith("tx1") for value in labels)
    assert "J1" in labels


def test_design_comparison_shows_paired_scores_and_exports_svg(tmp_path):
    comparison = pd.DataFrame({
        "probe_name": ["event_a", "event_b"],
        "control_score": [.5, .8], "selected_score": [.7, .8],
        "same_sequence": [False, True],
    })
    destination = tmp_path / "comparison.svg"

    figure = ps.pl.design_comparison(
        comparison, "control_score", "selected_score",
        same_sequence="same_sequence", save=destination,
    )

    assert destination.is_file()
    assert len(figure.axes[0].collections) == 2
    annotations = [text.get_text() for text in figure.axes[0].texts]
    assert any("alternative RBD" in value for value in annotations)
    assert any("same RBD" in value for value in annotations)


def test_candidate_record_cards_are_public_and_export_svg(tmp_path):
    candidates = pd.DataFrame({
        "candidate_id": ["c1", "c2"], "probe_name": ["external", "selected"],
        "gene_id": ["Pkm", "Pkm"], "transcript_id": ["tx1", "tx1"],
        "junction_id": ["J8", "J8"], "panel_role": ["control", "designed"],
        "rbd5": ["AACCGTAACCGTAACC", "AACCGTAACCGTAACC"],
        "rbd3": ["CCGTTGCACCGTTGCA", "TCGTTGCACCGTTGCA"],
        "shift_nt": [0, 1], "tm_full": [74.1, 74.5],
        "probeharbour_score": [.5, .6], "pass_tm_full": [True, True],
        "pass_repeat": [True, True], "has_5g_repeat": [False, False],
        "coordinates": [
            {"chrom": "5", "start": 100, "end": 316, "strand": "+",
             "rbd5_start": 100, "rbd5_end": 116, "rbd3_start": 300, "rbd3_end": 316},
            {"chrom": "5", "start": 101, "end": 317, "strand": "+",
             "rbd5_start": 101, "rbd5_end": 117, "rbd3_start": 301, "rbd3_end": 317},
        ],
    })
    single = tmp_path / "candidate.svg"
    comparison = tmp_path / "candidate_comparison.svg"

    detailed = ps.pl.candidate_card(candidates, "c2", save=single)
    compact = ps.pl.candidate_cards(candidates, save=comparison)

    assert single.is_file() and comparison.is_file()
    assert any("not calculated" in text.get_text() for text in detailed.axes[0].texts)
    assert not any("unknown gene" in text.get_text() for text in detailed.axes[0].texts)
    assert not any("vnot calculated" in text.get_text() for text in detailed.axes[0].texts)
    assert any(text.get_text() == "genomic 5′ RBD block" for text in detailed.axes[0].texts)
    assert any(text.get_text() == "genomic 3′ RBD block" for text in detailed.axes[0].texts)
    assert not any(text.get_text() == "genomic RBD blocks" for text in detailed.axes[0].texts)
    assert len(compact.axes) == 2
