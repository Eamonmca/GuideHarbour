# %% [markdown]
# # 2. Examine an existing 10x panel
#
# This tutorial treats the mouse 5K BED-derived ProbeStore as a dataset. We
# inspect panel size, gene coverage, junction geometry, sequence composition,
# and thermodynamics without changing the original input.

# %%
from pathlib import Path
import os
import matplotlib.pyplot as plt
import pandas as pd
import probeharbour as ps

# %%
PANEL_ROOT = Path(os.environ.get("PROBEHARBOUR_PANEL_ROOT", "reference/panels"))
REFERENCE_ROOT = Path(os.environ.get("PROBEHARBOUR_REFERENCE_ROOT", "reference/mouse_grcm38"))
BED = PANEL_ROOT / "5K/xenium_mouse_5K_gene_expression_panel_probe_locations.bed"
FASTA = REFERENCE_ROOT / "Mus_musculus.GRCm38.dna.primary_assembly.fa"

panel = ps.read_bed(BED, reference_fasta=FASTA, name="xenium_prime_5k_mouse")
panel.summary()

# %% [markdown]
# ## Panel composition
#
# `junction_spanning` comes directly from BED12 block geometry. A single-block
# probe is continuous in the genome; a multi-block probe crosses at least one
# genomic gap.

# %%
junction_summary = (
    panel.obs.groupby(["bed_block_count", "junction_spanning"], as_index=False)
    .agg(probes=("probe_name", "size"), genes=("gene_name", "nunique"))
)
junction_summary["percent"] = 100 * junction_summary["probes"] / len(panel.obs)
junction_summary

# %%
probes_per_gene = panel.obs.groupby("gene_name").size().sort_values(ascending=False)
probes_per_gene.describe(percentiles=[.1, .25, .5, .75, .9, .95, .99])

# %%
probes_per_gene.head(15).rename("probes").reset_index()

# %% [markdown]
# ## Sequence statistics
#
# The reconstructed target region is `rbd5 + rbd3`. These simple summaries do
# not depend on a scoring model.

# %%
panel.obs["binder"] = panel.obs["rbd5"] + panel.obs["rbd3"]
panel.obs["gc_percent_simple"] = panel.obs["binder"].map(
    lambda sequence: 100 * (sequence.count("G") + sequence.count("C")) / len(sequence)
)
panel.obs["lig_dinuc_simple"] = panel.obs["rbd5"].str[-1] + panel.obs["rbd3"].str[0]

panel.obs[["binder", "gc_percent_simple"]].describe()

# %%
panel.obs["lig_dinuc_simple"].value_counts().sort_index().rename("probes").reset_index()

# %% [markdown]
# `ps.pl` is the Scanpy-like plotting namespace. `panel_overview` accepts the
# ProbeStore directly and summarizes the same dataframe columns shown above.

# %%
ps.pl.panel_overview(panel, title="10x mouse 5K panel overview")

# %% [markdown]
# ## Thermodynamic characterization
#
# `characterize` adds calculated fields in place. It does not invent or
# overwrite a historic 10x score. The parameter set used is recorded alongside
# the results.

# %%
panel.characterize()
panel.obs[[
    "probe_name", "gene_name", "tm_5p", "tm_3p", "tm_full",
    "gc_percent", "p3_hairpin_dg", "p3_homodimer_dg",
    "has_simple_repeat", "has_5g_repeat",
]].head()

# %%
thermodynamic_summary = panel.obs[[
    "tm_5p", "tm_3p", "tm_full", "gc_percent",
    "p3_hairpin_dg", "p3_homodimer_dg",
]].describe().T
thermodynamic_summary

# %% [markdown]
# ## Thermodynamic visual overview

# `thermodynamics` expects characterization columns and tells the user to run
# `characterize()` if they are absent.

# %%
ps.pl.thermodynamics(panel, title="10x mouse 5K thermodynamics")

# %% [markdown]
# At this point `panel.obs` is an ordinary pandas table, so users can group,
# plot, filter, or export any statistic without learning a second query system.
