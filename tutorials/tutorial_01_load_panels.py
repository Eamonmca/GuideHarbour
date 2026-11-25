# %% [markdown]
# # 1. Load an existing 10x panel and a custom panel
#
# This tutorial starts with the files users actually have:
#
# - a 10x Xenium probe-locations BED;
# - the matching reference FASTA and GTF; and
# - a custom CSV containing probe names and two RBD arms.
#
# We create one `ProbeStore` for each panel, map the sequence-only custom probes
# back to transcripts and genomic coordinates, merge the panels, and save the
# result. Probe annotations live in `store.obs`, like observations in Scanpy.

# %%
from pathlib import Path
import os
import pandas as pd
from allos.transcript_data import TranscriptData
import probeharbour as ps

# %% [markdown]
# ## Input files
#
# BED, GTF, and FASTA must use the same genome assembly. The mouse 5K BED and
# custom panel in this example were designed against GRCm38/mm10.

# %%
PROJECT_ROOT = Path(ps.__file__).resolve().parents[1]
PANEL_ROOT = Path(os.environ.get("PROBEHARBOUR_PANEL_ROOT", "reference/panels"))
REFERENCE_ROOT = Path(os.environ.get("PROBEHARBOUR_REFERENCE_ROOT", "reference/mouse_grcm38"))

MOUSE_5K_BED = PANEL_ROOT / "5K/xenium_mouse_5K_gene_expression_panel_probe_locations.bed"
MOUSE_FASTA = REFERENCE_ROOT / "Mus_musculus.GRCm38.dna.primary_assembly.fa"
MOUSE_GTF = REFERENCE_ROOT / "Mus_musculus.GRCm38.102.gtf.gz"
CUSTOM_CSV = PROJECT_ROOT / "tutorials/data/96ISOS.csv"
OUTPUT = Path("/tmp/probeharbour_tutorials/mouse_5k_plus_96")

paths = pd.DataFrame({
    "input": ["10x BED", "reference FASTA", "annotation GTF", "custom CSV"],
    "path": [MOUSE_5K_BED, MOUSE_FASTA, MOUSE_GTF, CUSTOM_CSV],
})
paths["exists"] = paths["path"].map(Path.is_file)
paths

# %%
assert paths["exists"].all(), paths.loc[~paths["exists"]]
assert Path(f"{MOUSE_FASTA}.fai").is_file(), "The FASTA requires a .fai index"

# %% [markdown]
# ## Load the existing 10x panel from BED
#
# The 10x BED contains genomic locations but not RBD sequences. `read_bed`
# retrieves the targeted sequence from the FASTA, respects BED12 exon blocks,
# and records the source files and coordinate convention in metadata.
#
# `junction_spanning=True` means the BED record contains multiple genomic
# blocks. It is direct geometry from the 10x file—not an inferred score.

# %%
existing = ps.read_bed(
    MOUSE_5K_BED,
    reference_fasta=MOUSE_FASTA,
    name="xenium_prime_5k_mouse",
    panel_role="existing",
)
existing

# %%
existing.obs[[
    "probe_name", "gene_name", "rbd5", "rbd3", "chrom",
    "chrom_start", "chrom_end", "strand", "bed_block_count",
    "junction_spanning", "panel_role",
]].head()

# %% [markdown]
# ## Load the sequence-only custom panel
#
# The raw CSV uses `name`, `rbdm5`, and `rbdm3`. We state the normalization
# explicitly; no prepared replacement CSV is hidden behind the tutorial.

# %%
custom = ps.load_panel(
    CUSTOM_CSV,
    name="96_custom_probes",
    column_map={"name": "probe_name", "rbdm5": "rbd5", "rbdm3": "rbd3"},
    panel_role="add_on",
)
custom.obs.head()

# %% [markdown]
# The names begin with the gene symbol (`Ptbp2_ex10`, `Abat_ex13skip`, …), so
# we extract that visible prefix to constrain transcript searching. If a user
# already has a `gene_name` column, this line is unnecessary.

# %%
custom.obs["gene_name"] = custom.obs["probe_name"].str.split("_", n=1).str[0]

# %% [markdown]
# ## Map custom sequences to transcripts and coordinates
#
# `map_probes` searches each exact `rbd5 + rbd3` sequence in the gene's spliced
# transcripts. It calculates the actual arm coordinates and junction shift from
# the sequence match. Shifted probes are never forced back onto the junction.
# Unmapped rows remain in `.obs` with a reason.

# %%
td = TranscriptData(
    gtf_file=MOUSE_GTF,
    reference_fasta=MOUSE_FASTA,
    cache_dir="/tmp/probeharbour-allos-cache",
)
custom.map_probes(td)
custom.obs["mapping_status"].value_counts(dropna=False)

# %%
custom.obs[[
    "probe_name", "gene_name", "mapping_status", "transcript_id",
    "junction_id", "junction_shift", "junction_spanning", "on_junction",
    "chrom", "chrom_start", "chrom_end",
]].head(12)

# %% [markdown]
# A coordinate map makes the sequence→transcript result tangible. Blue and
# orange segments are the two RBD arms; the dotted gap is the spliced intron.

# %%
ps.pl.gene_map("Ptbp2", td, custom, title="Ptbp2 isoforms and mapped custom probes")

# %% [markdown]
# For a compact probe-only view, use `ps.pl.probe_map(custom, gene="Ptbp2")`.

# %% [markdown]
# 93 probes map exactly in this Ensembl GRCm38.102 annotation. The three
# retained unmapped rows (`Rbfox1_ex19`, `Nlgn1_A1A2`, and `Nlgn1_A1`) are
# visible annotation/reference mismatches, not silently discarded records.

# %%
assert custom.obs["mapping_status"].eq("mapped").sum() == 93
custom.obs.loc[
    custom.obs["mapping_status"].ne("mapped"),
    ["probe_name", "gene_name", "mapping_reason"],
]

# %% [markdown]
# ## Merge, save, and reopen
#
# Merging preserves `existing` and `add_on` roles plus source-store identity.
# Saving writes `probes.parquet` and `metadata.json`; loading does not recompute
# or overwrite the annotations.

# %%
combined = ps.ProbeStore.merge([existing, custom], name="mouse_5k_plus_96")
combined.obs.groupby(["panel_role", "source_store"]).size().rename("probes").reset_index()

# %%
saved = combined.save(OUTPUT)
reloaded = ps.load_panel(saved)
assert len(reloaded.obs) == len(existing.obs) + len(custom.obs)
print(f"Saved ProbeStore: {saved}")
reloaded.summary()
