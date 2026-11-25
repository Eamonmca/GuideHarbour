# %% [markdown]
# # 3. Compare ProbeHarbour scores with an external CSV
#
# This tutorial loads the real 10x all-candidates export as an external table,
# recalculates sequence features and a prospective ProbeHarbour score, and
# compares that score with `raw_design_score` from the CSV.
#
# The vendor score is an outcome only. It is never supplied to characterization
# or scoring, and it is never overwritten.

# %%
from pathlib import Path
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import probeharbour as ps

# %%
CSV = Path(os.environ.get("PROBEHARBOUR_ALL_CANDIDATES_CSV", "reference/all_candidates.csv"))
assert CSV.is_file(), (
    "Set PROBEHARBOUR_ALL_CANDIDATES_CSV to the vendor all-candidates CSV. "
    "The calibration file is not distributed with ProbeHarbour."
)

raw = pd.read_csv(CSV)
raw.head()

# %% [markdown]
# The required RBD columns already use canonical names. The export has no
# unique probe-name column because each row is one slide around a design target,
# so we create a transparent row identifier from the target name and row number.

# %%
store = ps.load_panel(CSV, name="3J6VC_all_candidates")
store.obs["probe_name"] = [
    f"{gene}:row:{index}" for index, gene in enumerate(store.obs["gene"], start=1)
]
store

# %% [markdown]
# Keep an immutable copy of the external score before doing any calculations.

# %%
historic_score = store.obs["raw_design_score"].copy()

# %% [markdown]
# `characterize` recalculates Tm, GC, repeats, hairpins, and dimers from sequence.
# `annotate_scores` then applies the packaged Prime 5K rules and writes to a new
# column named `probeharbour_score`.

# %%
store.characterize().annotate_scores(
    "prime_5k",
    score_column="probeharbour_score",
)
pd.testing.assert_series_equal(store.obs["raw_design_score"], historic_score)

# %%
comparison = store.obs.dropna(subset=["raw_design_score", "probeharbour_score"]).copy()
comparison["difference"] = comparison["probeharbour_score"] - comparison["raw_design_score"]
comparison["absolute_difference"] = comparison["difference"].abs()

metrics = pd.DataFrame([{
    "probes": len(comparison),
    "pearson_r": comparison[["raw_design_score", "probeharbour_score"]].corr().iloc[0, 1],
    "mae": comparison["absolute_difference"].mean(),
    "median_absolute_error": comparison["absolute_difference"].median(),
    "within_0.01": comparison["absolute_difference"].le(.01).mean(),
    "within_0.05": comparison["absolute_difference"].le(.05).mean(),
    "exact_at_4dp": comparison["raw_design_score"].round(4).eq(
        comparison["probeharbour_score"].round(4)
    ).mean(),
    "maximum_absolute_error": comparison["absolute_difference"].max(),
}])
metrics

# %% [markdown]
# ## Visual comparison
#
# The dashed line is exact agreement. Points away from it are useful diagnostic
# cases; they are not silently removed to improve the correlation.

# %%
ps.pl.score_comparison(
    comparison,
    "raw_design_score",
    "probeharbour_score",
    label_x="Score in external CSV",
    label_y="ProbeHarbour score",
    title="External versus prospective score",
)

# %%
fig, ax = plt.subplots(figsize=(5, 3.5))
comparison["difference"].hist(ax=ax, bins=40, color="#16a34a")
ax.axvline(0, color="#64748b", linestyle="--", linewidth=1)
ax.set(title="Score residuals", xlabel="ProbeHarbour − external", ylabel="probes")
fig.tight_layout()

# %% [markdown]
# ## Inspect the discrepancies
#
# The largest differences should be examined rather than explained away. Common
# causes include different Tm implementations near a threshold, vendor-only
# specificity evidence, or flags whose exact internal definition is unavailable.

# %%
largest_differences = comparison.nlargest(20, "absolute_difference")[[
    "probe_name", "gene", "rbd5", "rbd3", "junction",
    "raw_design_score", "probeharbour_score", "difference",
    "probe_flags", "combined_tm", "tm_full",
    "probeharbour_score_modifiers",
]]
largest_differences

# %%
comparison.assign(
    external_flag_present=comparison["probe_flags"].notna()
).groupby("external_flag_present").agg(
    probes=("probe_name", "size"),
    mae=("absolute_difference", "mean"),
    maximum_error=("absolute_difference", "max"),
)

# %% [markdown]
# High correlation demonstrates that the public rules reproduce the main score
# behavior. It does not imply that ProbeHarbour has access to undisclosed vendor
# calculations, so close-but-not-identical is the scientifically honest target.
