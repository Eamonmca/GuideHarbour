# %% [markdown]
# # GuideHarbour: design probes for the supplied guide constructs
# This worked example uses the packaged 15-guide table. A probe detects the
# guide-bearing construct, not the endogenous gene named in its label.

# %%
from importlib.resources import files
from pathlib import Path
import json

import matplotlib.pyplot as plt
import pandas as pd
import probeharbour as ph

guides = pd.read_csv(files("probeharbour").joinpath("data/example_guides.csv"), keep_default_na=False)
print(guides.assign(guide_length=guides.guide_sequence.str.len())[["guide_id", "guide_length"]].to_string(index=False))

# %% [markdown]
# ## Design across the guide/scaffold boundary
# RNA mode excludes the promoter/ligation prefix, preserves the supplied guide
# lengths, and uses any explicit transcribed_leader. The defaults retain at least
# 15 guide bases and require at least 3 differences from another full target
# window in this guide set. These are configurable project heuristics.

# %%
store = ph.design_guides(guides, chemistry="prime_5k", min_guide_overlap=15, min_panel_mismatches=3)
print(ph.guide_summary(store).to_string(index=False))

# %%
top = ph.top_guides(store, n=3)
columns = ["guide_id", "eligible_rank", "recalculated_score", "guide_overlap_nt", "tm_full",
           "nearest_other_mismatches", "rbd5", "rbd3", "modifiers"]
print(top[columns].to_string(index=False))

# %% [markdown]
# ## Visualise the candidates
# Top candidates are overlapping alternatives. rbd5/rbd3 are sense target
# segments; the explicit probe_arm_* columns hold reverse-complement binding
# pieces. They are not complete synthetic padlock probes.

# %%
overview = ph.gpl.guide_overview(store, n=3)
overview

# %%
detail = ph.gpl.plot_guide_candidates(store, "HIF1A-sg1", n=3)
detail

# %% [markdown]
# ## Compare the inherited v1 chemistry rules
# The v1 default uses 20 nt per arm, so fewer windows fit near the 5-prime end.
# Scores from different rulesets are not calibrated detection probabilities.

# %%
v1 = ph.design_guides(guides, chemistry="xenium_v1")
comparison = ph.guide_summary(store)[["guide_id", "n_candidates", "best_score"]].merge(
    ph.guide_summary(v1)[["guide_id", "n_candidates", "best_score"]],
    on="guide_id", suffixes=("_5k", "_v1"))
print(comparison.to_string(index=False))

# %% [markdown]
# ## Export and record the run
# No transcriptome screening or RNA accessibility analysis has been performed.
# The HTML report is self-contained and can be opened locally.

# %%
output = Path("guide_results_notebook")
output.mkdir(exist_ok=True)
ph.gpl.write_guide_report(store, output / "report.html", n=3)
store.write_csv(output / "candidates.csv")
top.to_csv(output / "top_candidates.csv", index=False)
ph.guide_summary(store).to_csv(output / "guide_summary.csv", index=False)
(output / "metadata.json").write_text(json.dumps(store.metadata, indent=2), encoding="utf-8")
overview.savefig(output / "overview.png", dpi=180, bbox_inches="tight")
detail.savefig(output / "HIF1A-sg1.png", dpi=180, bbox_inches="tight")
plt.close(overview)
plt.close(detail)
print(f"Results: {output.resolve()}")

# %% [markdown]
# ## Optional context and QC settings
# Only use full_construct when prefix + guide + scaffold is truly the intended
# target: ph.design_guides(guides, target_mode="full_construct", boundary="both").
# For an RNA with a known upstream leader, populate transcribed_leader and keep
# target_mode="rna". Use require_qc_pass=True to exclude all inherited chemistry
# modifiers, or min_score=0.5 for a score threshold. Always inspect empty results.
