from pathlib import Path
import os

import pandas as pd
import pytest

import probeharbour as ps


STORES = Path(os.environ.get("PROBEHARBOUR_TEST_STORES", "internal_tests/stores"))
PRIVATE_STORE = STORES / "96_isoform_panel_scorefix_20260710" / "probes.parquet"
requires_private_store = pytest.mark.skipif(
    not PRIVATE_STORE.is_file(),
    reason="set PROBEHARBOUR_TEST_STORES to run private panel regressions",
)


def _write_indexed_fasta(path, records):
    """Write a tiny fixed-width FASTA and matching .fai for importer tests."""
    offset = 0
    index_lines = []
    with path.open("wb") as handle:
        for name, sequence in records.items():
            header = f">{name}\n".encode()
            handle.write(header)
            offset += len(header)
            handle.write(f"{sequence}\n".encode())
            index_lines.append(f"{name}\t{len(sequence)}\t{offset}\t{len(sequence)}\t{len(sequence) + 1}")
            offset += len(sequence) + 1
    Path(f"{path}.fai").write_text("\n".join(index_lines) + "\n")


def test_read_bed_reconstructs_stranded_and_blocked_sequences(tmp_path):
    fasta = tmp_path / "reference.fa"
    _write_indexed_fasta(fasta, {"1": "AAAACCCCGGTTATAAACCCGGGG"})
    bed = tmp_path / "panel.bed"
    bed.write_text(
        "track name=demo\n"
        "chr1\t4\t12\tGeneA (ENSG1) | 7 | 2\t0\t+\t4\t12\t0\t2\t4,4\t0,4\n"
        "chr1\t12\t20\tGeneB (ENSG2) | 8 | 1\t0\t-\n"
    )

    store = ps.read_bed(bed, reference_fasta=fasta, name="demo", panel_role="existing")

    assert store.obs.loc[0, "rbd5"] + store.obs.loc[0, "rbd3"] == "AACCGGGG"
    assert store.obs.loc[1, "rbd5"] + store.obs.loc[1, "rbd3"] == "ATAAACCC"
    assert store.obs.loc[0, "junction_spanning"]
    assert store.obs.loc[0, "chrom_start"] == 4
    assert store.obs.loc[0, "gene_name"] == "GeneA"
    assert store.obs["panel_role"].eq("existing").all()
    assert store.metadata["coordinate_convention"] == "BED: zero-based, half-open"


def test_bed_score_is_not_mistaken_for_a_design_score(tmp_path):
    fasta = tmp_path / "reference.fa"
    _write_indexed_fasta(fasta, {"1": "AAAACCCCGGGGTTTT"})
    bed = tmp_path / "panel.bed"
    bed.write_text("chr1\t0\t8\tGeneA (ENSG1) | 7 | 1\t0\t-\n")

    summary = ps.read_bed(bed, reference_fasta=fasta).summary()

    assert pd.isna(summary.loc[0, "median_score"])
    assert summary.loc[0, "score_column"] is None


class FakePanelTranscriptData:
    def get_transcripts_by_gene_id(self, gene_id):
        return ["tx1", "tx2"] if gene_id == "ENSG1" else []

    def get_transcripts_by_gene_name(self, gene_name):
        return []

    def get_exon_coords_and_strand(self, transcript_id):
        exons = [[8, 0], [24, 12]] if transcript_id == "tx1" else [[8, 0], [24, 16]]
        return exons, -1

    def get_transcript_info(self, transcript_id):
        return {"chromosome": "1"}

    def get_chromosome(self, transcript_id):
        return "1"


def test_build_store_from_bed_annotates_transcripts_and_keeps_failures(tmp_path):
    fasta = tmp_path / "reference.fa"
    _write_indexed_fasta(fasta, {"1": "AAAACCCCGGGGTTTTAAAACCCC"})
    bed = tmp_path / "panel.bed"
    bed.write_text(
        "chr1\t4\t16\tGeneA (ENSG1) | 7 | 2\t0\t+\t4\t16\t0\t2\t4,4\t0,8\n"
        "chr1\t0\t8\tMissing (NO_GENE) | 8 | 1\t0\t+\n"
    )

    store = ps.build_store_from_bed(
        bed, reference_fasta=fasta, transcript_data=FakePanelTranscriptData(),
        characterize=False, name="annotated",
    )

    assert store.obs.loc[0, "mapping_status"] == "mapped"
    assert store.obs.loc[0, "supporting_transcripts"] == ["tx1"]
    assert store.obs.loc[0, "on_junction"]
    assert store.obs.loc[0, "isoform_specific_junction"]
    assert store.obs.loc[1, "mapping_reason"] == "gene_not_found"
    assert len(store.obs) == 2
    assert store.metadata["transcript_annotation"]["unmapped"] == 1
    assert store.metadata["build_pipeline"] == ["read_bed", "annotate_transcripts"]


@requires_private_store
def test_dataframe_first_panel_api_preserves_stored_scores():
    store = ps.load_panel(STORES / "96_isoform_panel_scorefix_20260710")
    before = store.obs["raw_design_score"].copy()

    summary = store.summary()
    explained = store.explain("prime_5k")

    pd.testing.assert_series_equal(store.obs["raw_design_score"], before)
    assert repr(store).startswith("ProbeStore(")
    assert summary.loc[0, "probes"] == 96
    assert summary.loc[0, "score_column"] == "raw_design_score"
    assert len(explained) == 96
    assert "stored_raw_design_score" in explained
    assert "recalculated_score" in explained


def test_design_junction_characterize_select_and_fasta(tmp_path):
    store = ps.design_junction("ACGT" * 30, 60, slides=range(-1, 2), candidate_prefix="demo")

    store.characterize()
    explanation = store.explain()
    store.obs["recalculated_score"] = explanation["recalculated_score"].values
    decisions = store.select(n_best=1)
    fasta = store.write_fasta(tmp_path / "demo.fasta")

    assert len(store.obs) == 3
    assert {"binder", "arm_length", "lig_dinuc", "tm_5p", "tm_3p", "tm_full", "gc_percent"} <= set(store.obs)
    store.annotate_scores("prime_5k", score_column="test_score")
    assert {"test_score", "test_score_modifiers", "test_score_base_score", "test_score_ruleset_id"} <= set(store.obs)
    assert decisions["selected"].sum() == 1
    assert fasta.read_text().count(">") == 3


def test_designed_store_with_empty_extra_can_be_saved(tmp_path):
    store = ps.design_junction("ACGT" * 30, 60, slides=range(-1, 2), candidate_prefix="demo")

    saved = store.save(tmp_path / "designed")
    loaded = ps.load_panel(saved)

    assert len(loaded.obs) == 3
    assert loaded.obs["extra"].isna().all()


def test_summary_uses_populated_custom_score_column():
    store = ps.design_junction("ACGT" * 30, 60, slides=range(-1, 2), candidate_prefix="demo")
    store.obs["probeharbour_score"] = [0.2, 0.8, 0.5]

    summary = store.summary()

    assert summary.loc[0, "score_column"] == "probeharbour_score"
    assert summary.loc[0, "median_score"] == 0.5


def test_characterize_keeps_duplicate_display_names_but_uses_unique_ids():
    store = ps.ProbeStore(pd.DataFrame({
        "probe_name": ["same", "same"],
        "rbd5": ["AACCGTAA", "AACCGTAA"],
        "rbd3": ["CCGTTGCA", "CCGTTGCA"],
    }))

    store.characterize()

    assert store.obs["probe_name"].tolist() == ["same", "same"]
    assert store.obs["candidate_id"].is_unique
    assert store.obs["tm_full"].notna().all()


class FakeTranscriptData:
    def get_transcripts_by_gene_id(self, gene):
        return []

    def get_transcripts_by_gene_name(self, gene):
        return ["tx1"] if gene == "Demo" else []

    def get_transcript_sequence(self, transcript_id):
        return "ACGT" * 30

    def get_exon_coords_and_strand(self, transcript_id):
        return ([[1040, 1000], [2080, 2000]], 1)

    def get_transcript_info(self, transcript_id):
        return {"chromosome": "1"}


def test_design_gene_uses_public_annotation_adapter_and_writes_bed(tmp_path):
    store = ps.design_gene("Demo", FakeTranscriptData(), slides=[0])

    bed = store.write_bed(tmp_path / "demo.bed")

    assert len(store.obs) == 1
    assert store.obs.loc[0, "gene_id"] == "Demo"
    assert store.obs.loc[0, "on_junction"]
    assert store.obs.loc[0, "junction_transcript_count"] == 1
    assert store.obs.loc[0, "target_coverage"] == 100
    assert store.obs.loc[0, "junction_specificity_class"] == "specific"
    coordinates = store.obs.loc[0, "coordinates"]
    assert coordinates["rbd5_start"] is not None
    assert coordinates["rbd3_start"] is not None
    assert bed.read_text().startswith("1\t")


def test_sequence_only_probe_maps_to_actual_transcript_and_arm_coordinates():
    store = ps.ProbeStore(pd.DataFrame({
        "probe_name": ["external"], "gene_name": ["Demo"],
        "rbd5": ["ACGT" * 4], "rbd3": ["ACGT" * 4],
    }))

    store.map_probes(FakeTranscriptData())

    row = store.obs.iloc[0]
    assert row["mapping_status"] == "mapped"
    assert row["transcript_id"] == "tx1"
    assert row["on_junction"]
    assert row["junction_shift"] == 0
    assert row["rbd5_chrom_start"] == 1024
    assert row["rbd5_chrom_end"] == 1040
    assert row["rbd3_chrom_start"] == 2000
    assert row["rbd3_chrom_end"] == 2016
    assert row["junction_transcript_count"] == 1
    assert row["target_coverage"] == 100
    assert row["junction_specificity_class"] == "specific"


@requires_private_store
def test_score_annotation_never_overwrites_historic_score():
    store = ps.load_panel(STORES / "96_isoform_panel_scorefix_20260710")
    historic = store.obs["raw_design_score"].copy()

    store.annotate_scores("prime_5k", score_column="prospective_score")

    pd.testing.assert_series_equal(store.obs["raw_design_score"], historic)
    assert store.obs["prospective_score"].notna().all()
    with pytest.raises(ValueError, match="already exists"):
        store.annotate_scores("prime_5k", score_column="prospective_score")


@requires_private_store
def test_store_query_save_load_merge_and_csv(tmp_path):
    store = ps.load_panel(STORES / "96_isoform_panel_scorefix_20260710")
    subset = store.query("gene_name == 'Ptbp2'")
    saved = subset.save(tmp_path / "ptbp2")
    loaded = ps.ProbeStore.load(saved)
    merged = ps.ProbeStore.merge([loaded, loaded], name="twice")
    csv_path = merged.write_csv(tmp_path / "twice.csv")

    assert len(subset.obs) == 2
    assert loaded.name == subset.name
    assert len(merged.obs) == 4
    assert merged.metadata["merged_from"] == [loaded.name, loaded.name]
    assert csv_path.is_file()


def test_load_panel_directory_and_merge_preserve_source_provenance(tmp_path):
    source = ps.ProbeStore(
        pd.DataFrame({"rbd5": ["AAAAAAAA"], "rbd3": ["CCCCCCCC"]}),
        name="bed_source",
        source_path=tmp_path / "source.bed",
        metadata={"reference_fasta": "/reference/genome.fa", "input_format": "bed"},
        panel_role="existing",
    )
    saved = source.save(tmp_path / "saved")

    reloaded = ps.load_panel(saved)
    merged = ps.ProbeStore.merge([reloaded], name="combined")

    assert reloaded.metadata["reference_fasta"] == "/reference/genome.fa"
    assert merged.metadata["source_stores"][0]["metadata"]["reference_fasta"] == "/reference/genome.fa"


def test_load_panel_accepts_csv(tmp_path):
    path = tmp_path / "addon.csv"
    pd.DataFrame({
        "probe_name": ["p1"], "gene": ["GeneA"],
        "rbd5": ["AACCGTAA"], "rbd3": ["CCGTTGCA"],
    }).to_csv(path, index=False)

    store = ps.load_panel(path, name="add_on")

    assert store.name == "add_on"
    assert store.source_path == path
    assert store.obs.loc[0, "probe_name"] == "p1"


def test_load_panel_rejects_unknown_format(tmp_path):
    path = tmp_path / "panel.txt"
    path.write_text("not a panel")
    with pytest.raises(ValueError, match="expected CSV or Parquet"):
        ps.load_panel(path)


def test_load_panel_normalizes_raw_schema_and_carries_role(tmp_path):
    path = tmp_path / "raw.csv"
    pd.DataFrame({
        "name": ["p1"], "left_arm": ["AACCGTAA"], "right_arm": ["CCGTTGCA"],
    }).to_csv(path, index=False)
    store = ps.load_panel(
        path, name="addon",
        column_map={"name": "probe_name", "left_arm": "rbd5", "right_arm": "rbd3"},
        panel_role="add_on",
    )
    assert store.obs.loc[0, "probe_name"] == "p1"
    assert store.obs.loc[0, "panel_role"] == "add_on"
    assert store.metadata["column_map"]["left_arm"] == "rbd5"


def test_merge_preserves_store_role_when_roles_are_not_repeated():
    existing = ps.ProbeStore(
        pd.DataFrame({"rbd5": ["AAAAAAAA"], "rbd3": ["CCCCCCCC"]}),
        name="source-a", panel_role="existing",
    )
    addon = ps.ProbeStore(
        pd.DataFrame({"rbd5": ["GGGGGGGG"], "rbd3": ["TTTTTTTT"]}),
        name="source-b", panel_role="add_on",
    )
    merged = ps.ProbeStore.merge([existing, addon])
    assert merged.obs["panel_role"].tolist() == ["existing", "add_on"]
