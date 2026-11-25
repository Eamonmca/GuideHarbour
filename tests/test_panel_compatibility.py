import pandas as pd

from probeharbour import ProbeStore, check_panel_compatibility, reverse_complement


def _panel(rows, name):
    return ProbeStore(pd.DataFrame(rows), name=name)


def test_reverse_complement():
    assert reverse_complement("AACCGT") == "ACGGTT"


def test_exact_duplicate_and_target_overlap_are_reported_separately():
    reference = _panel([
        {"probe_name": "old-1", "gene": "GeneA", "rbd5": "AACCGTAA", "rbd3": "CCGTTGCA"},
    ], "original")
    add_on = _panel([
        {"probe_name": "new-1", "gene": "GeneA", "rbd5": "AACCGTAA", "rbd3": "CCGTTGCA"},
    ], "add-on")

    add_on.check_panel_compatibility(reference, calculate_heterodimer=False)
    row = add_on.obs.iloc[0]
    assert bool(row["duplicate_binder"])
    assert bool(row["duplicate_arm"])
    assert bool(row["target_gene_already_present"])
    assert not bool(row["panel_compatibility_pass"])
    assert "duplicate_binder" in row["panel_compatibility_reasons"]
    assert "target_already_present" not in row["panel_compatibility_reasons"]


def test_reverse_complement_conflict_and_clean_probe():
    old_binder = "AACCGTAACCGTTGCA"
    reference = _panel([
        {"probe_name": "old-1", "gene": "GeneA", "rbd5": old_binder[:8], "rbd3": old_binder[8:]},
    ], "original")
    complement = reverse_complement(old_binder)
    query = _panel([
        {"probe_name": "conflict", "gene": "GeneB", "rbd5": complement[:8], "rbd3": complement[8:]},
        {"probe_name": "clean", "gene": "GeneC", "rbd5": "AAAAAAAA", "rbd3": "CCCCCCCC"},
    ], "add-on")

    result = check_panel_compatibility(query, reference, calculate_heterodimer=False)
    conflict = result.set_index("candidate_id").loc["conflict"]
    clean = result.set_index("candidate_id").loc["clean"]
    assert bool(conflict["reverse_complement_binder"])
    assert conflict["longest_panel_complement"] == 16
    assert not bool(conflict["panel_compatibility_pass"])
    assert bool(clean["panel_compatibility_pass"])


def test_default_complementary_run_threshold_is_fourteen():
    reference_binder = "ACGTCAGTCCGATGCA"
    reference = _panel([{
        "probe_name": "old-1",
        "rbd5": reference_binder[:8],
        "rbd3": reference_binder[8:],
    }], "original")
    reference_rc = reverse_complement(reference_binder)
    # The first 13 bases match exactly; the suffix is selected to break the run.
    query_binder = reference_rc[:13] + (
        "AAA" if reference_rc[13] != "A" else "CCC"
    )
    query = _panel([{
        "probe_name": "new-1",
        "rbd5": query_binder[:8],
        "rbd3": query_binder[8:],
    }], "add-on")

    result = check_panel_compatibility(
        query, reference, calculate_heterodimer=False,
    ).iloc[0]

    assert result["longest_panel_complement"] == 13
    assert "long_complementary_run" not in result["panel_compatibility_reasons"]


def test_primer3_heterodimer_columns_are_populated():
    reference = _panel([
        {"probe_name": "old-1", "gene": "GeneA", "rbd5": "AACCGTAA", "rbd3": "CCGTTGCA"},
    ], "original")
    query = _panel([
        {"probe_name": "new-1", "gene": "GeneB", "rbd5": "TGCAACGG", "rbd3": "TTACGGTT"},
    ], "add-on")
    result = check_panel_compatibility(query, reference)
    assert pd.notna(result.loc[0, "worst_panel_heterodimer_dg"])
    assert result.loc[0, "worst_panel_heterodimer_probe"] == "old-1"


def test_query_to_existing_panel_spacing_is_evaluated_without_full_panel_self_audit():
    reference = _panel([{
        "probe_name": "old-1", "rbd5": "AACCGTAA", "rbd3": "CCGTTGCA",
        "chrom": "1", "chrom_start": 100, "chrom_end": 132,
    }], "original")
    query = _panel([{
        "probe_name": "new-1", "rbd5": "ATATCGCG", "rbd3": "GCGCATAT",
        "coordinates": {"chrom": "1", "start": 160, "end": 192, "strand": "+"},
    }], "add-on")

    result = check_panel_compatibility(
        query, reference, calculate_heterodimer=False, min_distance=75,
    ).iloc[0]

    assert result["nearest_existing_panel_probe_name"] == "old-1"
    assert result["nearest_existing_panel_probe_distance"] == 28
    assert bool(result["existing_panel_spacing_violation"])
    assert "existing_panel_spacing" in result["panel_compatibility_reasons"]

    query.obs["panel_record_id"] = ["add-on:new-1:0"]
    query.check_panel_compatibility(reference, calculate_heterodimer=False, min_distance=75)
    assert query.obs.loc[0, "nearest_existing_panel_probe_distance"] == 28
    assert not bool(query.obs.loc[0, "existing_panel_compatibility_pass"])


def test_junction_probe_distance_uses_bound_arms_not_spanned_intron():
    reference = _panel([{
        "probe_name": "intronic-old", "rbd5": "AACCGTAA", "rbd3": "CCGTTGCA",
        "chrom": "1", "chrom_start": 500, "chrom_end": 532,
    }], "original")
    query = _panel([{
        "probe_name": "junction-new", "rbd5": "ATATCGCG", "rbd3": "GCGCATAT",
        "coordinates": {
            "chrom": "1", "start": 100, "end": 1032, "strand": "+",
            "rbd5_start": 100, "rbd5_end": 116,
            "rbd3_start": 1016, "rbd3_end": 1032,
        },
    }], "add-on")

    result = check_panel_compatibility(
        query, reference, calculate_heterodimer=False, min_distance=75,
    ).iloc[0]

    assert result["nearest_existing_panel_probe_distance"] == 384
    assert not bool(result["existing_panel_spacing_violation"])


def test_role_aware_merge_and_internal_spacing_exclude_self():
    existing = _panel([
        {"probe_name": "old-1", "gene": "GeneA", "rbd5": "AACCGTAA", "rbd3": "CCGTTGCA",
         "chrom": "1", "chrom_start": 100, "chrom_end": 132, "strand": "+"},
        {"probe_name": "old-2", "gene": "GeneB", "rbd5": "GGGGAATT", "rbd3": "TTAACCCC",
         "chrom": "1", "chrom_start": 400, "chrom_end": 432, "strand": "+"},
    ], "original")
    add_on = _panel([
        {"probe_name": "new-1", "gene": "GeneC", "rbd5": "ATATCGCG", "rbd3": "GCGCATAT",
         "chrom": "1", "chrom_start": 160, "chrom_end": 192, "strand": "+"},
    ], "addon")

    combined = ProbeStore.merge([existing, add_on], roles=["existing", "add_on"], name="combined")
    combined.check_internal_compatibility(calculate_heterodimer=False, min_distance=75)
    new = combined.obs.set_index("probe_name").loc["new-1"]
    assert new["panel_role"] == "add_on"
    assert new["nearest_probe_name"] == "old-1"
    assert new["nearest_probe_distance"] == 28
    assert new["nearest_existing_probe_distance"] == 28
    assert pd.isna(new["nearest_add_on_probe_distance"])
    assert bool(new["panel_spacing_violation"])
    assert "panel_spacing_violation" in new["panel_compatibility_reasons"]
    assert not bool(new["panel_compatibility_pass"])


def test_internal_sequence_check_finds_addon_addon_conflict_not_self():
    binder = "AACCGTAACCGTTGCA"
    add_on = _panel([
        {"probe_name": "new-1", "rbd5": binder[:8], "rbd3": binder[8:]},
        {"probe_name": "new-2", "rbd5": binder[:8], "rbd3": binder[8:]},
    ], "addon")
    combined = ProbeStore.merge([add_on], roles=["add_on"])
    combined.check_internal_compatibility(calculate_heterodimer=False)
    assert combined.obs["duplicate_binder"].all()
    assert combined.obs["duplicate_panel_probe"].notna().all()


def test_internal_audit_does_not_overwrite_existing_panel_evidence():
    external_binder = "AACCGTAACCGTTGCA"
    existing = _panel([{
        "probe_name": "existing-10x",
        "rbd5": external_binder[:8],
        "rbd3": external_binder[8:],
    }], "existing")
    external_complement = reverse_complement(external_binder)
    add_on = _panel([
        {
            "probe_name": "new-1",
            "rbd5": external_complement[:8],
            "rbd3": external_complement[8:],
        },
        {
            "probe_name": "new-2",
            "rbd5": "AAAAAAAA",
            "rbd3": "CCCCCCCC",
        },
    ], "add-on")

    add_on.check_panel_compatibility(existing, calculate_heterodimer=False)
    assert add_on.obs.loc[0, "existing_panel_longest_complement"] == 16
    assert (
        add_on.obs.loc[0, "existing_panel_closest_probe"]
        == "existing-10x"
    )

    add_on.check_internal_compatibility(calculate_heterodimer=False)

    # The generic columns now describe add-on/add-on evidence, while the
    # external namespace remains internally consistent.
    assert add_on.obs.loc[0, "longest_panel_complement"] < 12
    assert add_on.obs.loc[0, "existing_panel_longest_complement"] == 16
    assert "long_complementary_run" in add_on.obs.loc[
        0, "existing_panel_compatibility_reasons"
    ]
    assert (
        add_on.obs.loc[0, "existing_panel_closest_probe"]
        == "existing-10x"
    )


def test_internal_compatibility_repairs_duplicate_panel_record_ids():
    add_on = _panel([
        {"candidate_id": "gene-a", "probe_name": "gene-a", "rbd5": "AACCGTAA",
         "rbd3": "CCGTTGCA", "panel_record_id": "add_on:0"},
        {"candidate_id": "gene-b", "probe_name": "gene-b", "rbd5": "ATATCGCG",
         "rbd3": "GCGCATAT", "panel_record_id": "add_on:0"},
    ], "addon")

    add_on.check_internal_compatibility(calculate_heterodimer=False)

    assert add_on.obs["panel_record_id"].is_unique
    assert add_on.obs["panel_compatibility_pass"].notna().all()


def test_merge_requires_one_role_per_store():
    panel = _panel([{"rbd5": "AAAAAAAA", "rbd3": "CCCCCCCC"}], "one")
    try:
        ProbeStore.merge([panel], roles=["existing", "add_on"])
    except ValueError as error:
        assert "one value per store" in str(error)
    else:
        raise AssertionError("Expected mismatched roles to fail")
