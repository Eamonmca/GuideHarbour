# %% [markdown]
# # 4. Design a literature-backed isoform panel
#
# This tutorial designs probes for four well-characterized alternative-splicing
# events and checks them against probes from a real custom panel:
#
# - **Nrxn1 SS4 inclusion/skipping**, which changes neurexin ligand interactions;
# - **Gria2 flip/flop**, which changes AMPA-receptor kinetics;
# - **Snap25 exon 5a/5b**, a developmental neuronal isoform switch; and
# - **Pkm exon 9/10**, the mutually exclusive PKM1/PKM2 event.
#
# The published biology motivates the targets; the annotation and reference
# sequence determine what ProbeHarbour can actually design. The external probes
# are used only as positive controls after candidate generation.

# %%
from pathlib import Path
import os
import sys
import pandas as pd
from allos.transcript_data import TranscriptData
import probeharbour as ps

# %% [markdown]
# ## References and input files
#
# These are primary studies describing the selected events:
#
# - [Nrxn1 SS4 regulation](https://pubmed.ncbi.nlm.nih.gov/22196734/)
# - [Gria2 flip/flop receptor kinetics](https://pubmed.ncbi.nlm.nih.gov/10704491/)
# - [Snap25 exon 5 developmental switch](https://pubmed.ncbi.nlm.nih.gov/15470145/)
# - [Pkm mutually exclusive exon regulation](https://pubmed.ncbi.nlm.nih.gov/22307054/)
#
# All genomic inputs below use GRCm38. Mixing assemblies would produce plausible
# but incorrect coordinates, so the paths are checked before any design work.

# %%
PROJECT_ROOT = Path(ps.__file__).resolve().parents[1]
REFERENCE_ROOT = Path(os.environ.get("PROBEHARBOUR_REFERENCE_ROOT", "reference/mouse_grcm38"))
GTF = REFERENCE_ROOT / "Mus_musculus.GRCm38.102.gtf.gz"
FASTA = REFERENCE_ROOT / "Mus_musculus.GRCm38.dna.primary_assembly.fa"
CONTROL_CSV = PROJECT_ROOT / "tutorials/data/96ISOS.csv"
MOUSE_5K_BED = Path(os.environ.get(
    "PROBEHARBOUR_EXISTING_BED",
    "reference/panels/5K/xenium_mouse_5K_gene_expression_panel_probe_locations.bed",
))
OUTPUT = Path("/tmp/probeharbour_tutorials/literature_isoform_panel")
FIGURE_DIR = PROJECT_ROOT / "figures"
MOUSE_CDNA_BLAST_DB = Path(os.environ.get(
    "PROBEHARBOUR_BLAST_DB", "reference/blast/GRCm38_ens102_cdna"
))
BLASTN = Path(sys.executable).with_name("blastn")

assert GTF.is_file() and FASTA.is_file() and CONTROL_CSV.is_file() and MOUSE_5K_BED.is_file()
assert Path(f"{FASTA}.fai").is_file(), "The reference FASTA requires a .fai index"
assert MOUSE_CDNA_BLAST_DB.with_suffix(".nhr").is_file() and BLASTN.is_file()

# %% [markdown]
# ## Load the transcript reference
#
# `TranscriptData` provides exon models and spliced transcript sequences.
# ProbeHarbour uses its public methods: there is no prepared coordinate table.

# %%
td = TranscriptData(
    gtf_file=GTF,
    reference_fasta=FASTA,
    cache_dir="/tmp/probeharbour-allos-cache",
)

# %% [markdown]
# ## Load the experimental positive controls
#
# The CSV contains only a name and the two 16-nt RBD arms. We normalize those
# visible columns, retain the eight probes of interest, and map each exact
# `rbd5 + rbd3` sequence back to the transcript annotation.

# %%
genes = ["Nrxn1", "Gria2", "Snap25", "Pkm"]
controls = ps.load_panel(
    CONTROL_CSV,
    name="literature_event_controls",
    column_map={"name": "probe_name", "rbdm5": "rbd5", "rbdm3": "rbd3"},
)
controls.obs["gene_name"] = controls.obs["probe_name"].str.split("_", n=1).str[0]
controls = ps.ProbeStore(
    controls.obs.loc[controls.obs["gene_name"].isin(genes)].copy(),
    name="literature_event_controls",
)
controls.map_probes(td)

controls.obs[[
    "probe_name", "gene_name", "mapping_status", "transcript_id",
    "junction_id", "junction_shift", "junction_spanning",
]]

# %%
assert len(controls.obs) == 8
assert controls.obs["mapping_status"].eq("mapped").all()

# %% [markdown]
# ## See the isoforms before designing
#
# Grey blocks are annotated exons. Blue and orange marks show where the two
# arms of each mapped control bind. Each alternative form maps to a different
# transcript junction rather than merely carrying a different label.

# %%
for gene in genes:
    ps.pl.gene_map(
        gene,
        td,
        controls,
        show_overview=False,
        title=f"{gene}: annotated isoforms and external positive controls",
    )

# %% [markdown]
# This is the compact event-region view. Set `show_overview=True` (the default)
# to add the same supporting transcripts across their full locus. Long introns
# are compressed automatically so exons remain readable; the overview says
# explicitly when it is schematic, while the detailed panel always retains
# exact genomic scaling. Use `overview_scale="genomic"` to disable compression.

# %% [markdown]
# ## Generate candidates independently
#
# `design_gene` enumerates every annotated exon-exon junction for every
# transcript of a gene. Thirteen placements are considered at each junction:
# shifts −6 through +6. The external sequences are not passed to the designer.
# This explicit range matters because the Nrxn1 SS4 control uses shift +6.

# %%
designed_stores = []
for gene in genes:
    gene_design = ps.design_gene(
        gene,
        td,
        chemistry="prime_5k",
        slides=range(-6, 7),
        dataset_id="tutorial_04_literature_events",
    )
    gene_design.obs["gene_name"] = gene
    designed_stores.append(gene_design)

designed = ps.ProbeStore.merge(designed_stores, name="literature_isoform_candidates")
designed.summary()

# %% [markdown]
# ## Test exact reproduction
#
# An exact match requires both RBD arms to agree. Matching a gene or approximate
# coordinate is not counted. We also compare the transcript junction inferred
# from the external sequence with the junction assigned during design.

# %%
control_keys = controls.obs[[
    "probe_name", "gene_name", "rbd5", "rbd3", "junction_id", "junction_shift",
    "chrom", "strand", "rbd5_chrom_start", "rbd5_chrom_end",
    "rbd3_chrom_start", "rbd3_chrom_end",
]].rename(columns={
    "junction_id": "control_junction_id",
    "junction_shift": "control_junction_shift",
})

design_keys = designed.obs[[
    "candidate_id", "gene_name", "rbd5", "rbd3", "junction_id", "shift_nt",
    "transcript_id", "coordinates",
]].rename(columns={"junction_id": "designed_junction_id"})
for arm in ("rbd5", "rbd3"):
    for edge in ("start", "end"):
        design_keys[f"designed_{arm}_chrom_{edge}"] = design_keys["coordinates"].map(
            lambda value: value[f"{arm}_{edge}"]
        )
design_keys["designed_chrom"] = design_keys["coordinates"].map(lambda value: value["chrom"])
design_keys["designed_strand"] = design_keys["coordinates"].map(lambda value: value["strand"])

reproduction = control_keys.merge(
    design_keys,
    on=["gene_name", "rbd5", "rbd3"],
    how="left",
)
reproduction["exact_rbd_match"] = reproduction["candidate_id"].notna()
reproduction["same_junction"] = (
    reproduction["control_junction_id"] == reproduction["designed_junction_id"]
)
reproduction["exact_coordinates"] = (
    reproduction["chrom"].astype(str).eq(reproduction["designed_chrom"].astype(str))
    & reproduction["strand"].eq(reproduction["designed_strand"])
    & reproduction["rbd5_chrom_start"].eq(reproduction["designed_rbd5_chrom_start"])
    & reproduction["rbd5_chrom_end"].eq(reproduction["designed_rbd5_chrom_end"])
    & reproduction["rbd3_chrom_start"].eq(reproduction["designed_rbd3_chrom_start"])
    & reproduction["rbd3_chrom_end"].eq(reproduction["designed_rbd3_chrom_end"])
)
reproduction[[
    "probe_name", "gene_name", "exact_rbd_match", "same_junction", "exact_coordinates",
    "control_junction_shift", "shift_nt", "transcript_id",
]]

# %%
per_probe = reproduction.groupby("probe_name", as_index=False).agg(
    gene_name=("gene_name", "first"),
    exact_rbd_match=("exact_rbd_match", "any"),
    same_junction=("same_junction", "any"),
    exact_coordinates=("exact_coordinates", "any"),
    matching_design_rows=("candidate_id", "count"),
)
per_probe

# %%
assert per_probe["exact_rbd_match"].all()
assert per_probe["same_junction"].all()
assert per_probe["exact_coordinates"].all()

# %% [markdown]
# All eight external probes are regenerated with exact RBDs and arm coordinates.
# More than one design row
# can carry the same sequence when annotated transcripts share a junction; that
# is transcript redundancy, not eight different experimental probes.
#
# Exact reproduction is a useful completeness check, but it does not prove that
# an external probe is optimal. We therefore characterize and score the complete
# independently generated candidate set next.

# %% [markdown]
# ## Characterize and score candidates
#
# Characterization calculates sequence composition and thermodynamics.
# `annotate_scores` applies the packaged Prime 5K rules to a new column. The
# explanation column records the modifiers responsible for each score.

# %%
designed.characterize().annotate_scores(
    "prime_5k",
    score_column="probeharbour_score",
)

designed.obs[[
    "candidate_id", "gene_name", "junction_id", "shift_nt", "tm_full",
    "gc_percent", "probeharbour_score", "probeharbour_score_modifiers",
]].head()

# %% [markdown]
# Restrict the final mini-panel to the eight biologically chosen junctions, then
# retain the highest-scoring unique sequence for each junction. This separates
# biological target selection from sequence-level probe ranking.

# %%
target_junctions = set(controls.obs["junction_id"])
target_candidates = designed.obs.loc[
    designed.obs["junction_id"].isin(target_junctions)
].copy()
target_candidates = target_candidates.sort_values(
    ["junction_id", "probeharbour_score", "candidate_id"],
    ascending=[True, False, True],
).drop_duplicates(["junction_id", "rbd5", "rbd3"])

selected = target_candidates.groupby("junction_id", as_index=False).head(1).copy()
event_names = (
    controls.obs[["junction_id", "probe_name"]]
    .drop_duplicates("junction_id")
    .set_index("junction_id")["probe_name"]
)
selected["probe_name"] = selected["junction_id"].map(event_names) + "_design"
selected["selection_status"] = "selected: highest score for target junction"

# Evaluate the proposed probes in the context in which they would be ordered:
# first against every existing 10x binder, then within the combined physical
# panel for nearest-neighbour spacing and pairwise sequence conflicts.
existing = ps.read_bed(
    MOUSE_5K_BED,
    reference_fasta=FASTA,
    name="xenium_prime_5k_mouse",
    panel_role="existing",
)
selected_store = ps.ProbeStore(
    selected, name="literature_isoform_panel", panel_role="add_on"
).check_panel_compatibility(existing)

# Full-binder transcriptome BLAST supplies an independent specificity record,
# including the exact command, database version and retained raw hits.
blast = ps.LocalBlastBackend(
    MOUSE_CDNA_BLAST_DB,
    database_id="GRCm38_Ensembl_cDNA",
    database_version="102",
    executable=BLASTN,
)
ps.check_specificity(selected_store, blast)
selected_store.check_internal_compatibility()
selected = selected_store.obs.copy()
selected[[
    "probe_name", "gene_name", "junction_id", "candidate_id", "shift_nt",
    "rbd5", "rbd3", "probeharbour_score", "panel_compatibility_pass",
    "nearest_existing_panel_probe_distance", "existing_panel_spacing_violation",
    "nearest_probe_distance", "panel_spacing_violation",
]].sort_values(["gene_name", "junction_id"])

# %% [markdown]
# ## Did selection improve on the exact external match?
#
# Reproducing a control proves that the candidate generator is complete, but it
# does not prove that the control is the highest-scoring placement. We now score
# the eight external controls prospectively with exactly the same rules used for
# the generated pool, then compare each one with the independently selected
# candidate at the same junction. The external score is still never used during
# candidate generation or selection.

# %%
scored_controls = controls.copy().characterize().annotate_scores(
    "prime_5k",
    score_column="control_probeharbour_score",
)
scored_controls.check_panel_compatibility(existing)
ps.check_specificity(scored_controls, blast)
scored_controls.check_internal_compatibility()

control_comparison = scored_controls.obs[[
    "probe_name", "junction_id", "rbd5", "rbd3", "junction_shift",
    "control_probeharbour_score", "control_probeharbour_score_modifiers",
]].merge(
    selected[[
        "junction_id", "probe_name", "rbd5", "rbd3", "shift_nt",
        "probeharbour_score", "probeharbour_score_modifiers",
    ]].rename(columns={
        "probe_name": "selected_probe_name",
        "rbd5": "selected_rbd5",
        "rbd3": "selected_rbd3",
        "shift_nt": "selected_shift",
        "probeharbour_score": "selected_probeharbour_score",
        "probeharbour_score_modifiers": "selected_score_modifiers",
    }),
    on="junction_id",
    how="inner",
)
control_comparison["same_sequence"] = (
    control_comparison["rbd5"].eq(control_comparison["selected_rbd5"])
    & control_comparison["rbd3"].eq(control_comparison["selected_rbd3"])
)
control_comparison["score_gain"] = (
    control_comparison["selected_probeharbour_score"]
    - control_comparison["control_probeharbour_score"]
)

control_comparison[[
    "probe_name", "junction_shift", "selected_shift", "same_sequence",
    "control_probeharbour_score", "selected_probeharbour_score", "score_gain",
]].sort_values("score_gain", ascending=False)

# %%
comparison_figure = ps.pl.design_comparison(
    control_comparison,
    "control_probeharbour_score",
    "selected_probeharbour_score",
    label="probe_name",
    same_sequence="same_sequence",
    reference_label="exact external match",
    designed_label="selected candidate",
    title="Exact external match versus selected candidate",
    save=FIGURE_DIR / "external_match_vs_selected_design.svg",
)

# %% [markdown]
# ## Inspect the candidates as records, not just points
#
# A score difference is only useful when its underlying sequence and evidence
# can be inspected. `ps.pl.candidate_cards` puts the exact external probe beside
# the independently selected candidate for the event with the largest score
# gain. The compact cards show sequence, placement, thermodynamics, QC,
# specificity and panel-spacing evidence. An analysis that has not yet been run
# is labelled **not calculated**; absence of evidence is never displayed as a
# pass.

# %%
comparison_event = control_comparison.sort_values(
    "score_gain", ascending=False
).iloc[0]
comparison_junction = comparison_event["junction_id"]

external_card = scored_controls.obs.loc[
    scored_controls.obs["junction_id"].eq(comparison_junction)
].iloc[[0]].copy()
external_card["card_label"] = "Exact external match"
external_card["panel_role"] = "external control"
external_card["selection_status"] = "imported custom/external probe"
external_card["probeharbour_score"] = external_card["control_probeharbour_score"]
external_card["probeharbour_score_base_score"] = external_card[
    "control_probeharbour_score_base_score"
]
external_card["probeharbour_score_modifiers"] = external_card[
    "control_probeharbour_score_modifiers"
]
external_card["probeharbour_score_ruleset_id"] = external_card[
    "control_probeharbour_score_ruleset_id"
]

selected_card = selected.loc[
    selected["junction_id"].eq(comparison_junction)
].iloc[[0]].copy()
selected_card["card_label"] = "ProbeHarbour-selected candidate"
selected_card["panel_role"] = "designed"

candidate_pair = pd.concat([external_card, selected_card], ignore_index=True)
ps.pl.candidate_cards(
    candidate_pair,
    title=f"{comparison_event['probe_name']}: exact match versus selected design",
    card_title="card_label",
    save=FIGURE_DIR / "external_match_vs_selected_candidate_cards.svg",
)

# The detailed card is useful for reviewing one candidate or attaching an
# auditable design record to a report. SVG output remains fully vector based;
# use `dpi=600` when requesting a PNG instead.
ps.pl.candidate_card(
    selected_card,
    title="Designed isoform-probe candidate",
    subtitle="A self-contained view of the sequence, target, score and available evidence",
    save=FIGURE_DIR / "selected_candidate_record_card.svg",
)
ps.pl.candidate_card(
    external_card,
    title="Custom/external isoform-probe candidate",
    subtitle="The same evidence contract applied to an imported probe sequence",
    save=FIGURE_DIR / "custom_candidate_record_card.svg",
)

# %% [markdown]
# The result is deliberately not summarized as “all selected probes are better.”
# Five selected candidates are the exact external sequence, two are alternative
# sequences with the same prospective score, and one (`Pkm_ex9-202`) improves by
# 0.0212. Equal scores reflect the resolution of the current ligation/flag model;
# transcriptome specificity or panel compatibility may still separate tied
# candidates and should be evaluated before ordering.

# %%
assert control_comparison["same_sequence"].sum() == 5
assert control_comparison["score_gain"].ge(0).all()
assert control_comparison["score_gain"].gt(0).sum() == 1

# %% [markdown]
# ## View the proposed add-on beside the existing 10x panel
#
# A probe design is ultimately added to a physical panel. We therefore load the
# existing mouse 5K BED and merge it with the selected design before plotting.
# Purple tracks are exact BED blocks from the existing panel. Blue/orange tracks
# are the two arms of mapped or designed probes. Role labels make their source
# explicit, and BED probes remain visible even without transcript assignment.
# The upper tracks show the same supporting transcripts across their complete
# locus and mark every supplied probe; their shaded interval is the exact
# genomic detail expanded below.

# %%
panel = ps.ProbeStore(selected, name="literature_isoform_panel", panel_role="add_on")
combined_view = ps.ProbeStore.merge(
    [existing, panel],
    name="mouse_5k_plus_literature_isoforms",
    roles=["existing", "add_on"],
)

for gene in genes:
    ps.pl.gene_map(
        gene,
        td,
        combined_view,
        title=f"{gene}: existing 10x probes and proposed add-on probes",
        save=FIGURE_DIR / f"{gene}_existing_and_add_on.svg",
    )

# %% [markdown]
# Each call above also writes a publication-ready SVG to the project-local,
# Git-ignored `figures/` directory. SVG is vector output and therefore
# resolution-independent. For a high-resolution raster figure, use
# `save="figure.png", dpi=600`; `transparent=True` is available for either.

# %% [markdown]
# ## Save a portable result
#
# The selected rows retain genomic coordinates, provenance, characterization,
# and score explanations. Saving a `ProbeStore` writes `probes.parquet` and
# `metadata.json`; BED and FASTA exports are suitable for downstream inspection.

# %%
saved = panel.save(OUTPUT)
panel.write_bed(OUTPUT / "literature_isoform_panel.bed")
panel.write_fasta(OUTPUT / "literature_isoform_panel.fasta")

print(f"Saved ProbeStore: {saved}")
panel.summary()

# %% [markdown]
# This example establishes three separate claims:
#
# 1. the events have an explicit biological rationale;
# 2. all eight positive-control RBDs and arm coordinates can be regenerated;
# 3. final selection uses public ProbeHarbour calculations, not the controls.
#
# This tutorial performs the transcriptome BLAST and both existing-panel and
# add-on/add-on compatibility checks. A production order should exclude or
# deliberately resolve every reported conflict; genomic BLAST can be added when
# intergenic or intronic hybridization is also in scope.
