# GuideHarbour

**Design, score and visualise detection probes for an existing guide RNA and its adjacent sequence.**

This is a personal adaptation of [ProbeHarbour](https://github.com/cobioda/ProbeHarbour). It keeps the original `probeharbour` Python API and `ProbeStore`, adds a sequence-first guide workflow, and requires no GTF or genome reference for local design. It detects the supplied guide constructs; it does not design or score CRISPR perturbations.

The original code and Apache-2.0 license are retained. See [upstream provenance](UPSTREAM.md) and the [original ProbeHarbour README](README_ProbeHarbour.md).

## Start with the included guides

Clone your separate repository, then install in a fresh virtual environment:

```bash
git clone https://github.com/Eamonmca/GuideHarbour.git
cd GuideHarbour
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[guide-ui]'
guideharbour-dashboard
```

On Windows PowerShell, activate the environment with `.venv\Scripts\Activate.ps1` instead of the `source` command.

The dashboard starts with the supplied 15-guide table. Edit sequences or upload another CSV, select the chemistry and window settings, then choose **Design and rank probes**. Each guide has its own shortlist, sequence plot and exclusion reasons. Downloads include all candidates, the shortlist and an offline HTML report.

For command-line use, installation with `.[guides]` is sufficient:

```bash
guideharbour-design --example --output guide_results --top 3
# Or supply a CSV:
guideharbour-design examples/guides.csv --output my_results --chemistry prime_5k
```

If shell entry points are unavailable:

```bash
python -m probeharbour.guide_cli --example --output guide_results
python -m probeharbour._guide_dashboard_launcher
```

The CSV in [examples/guides.csv](examples/guides.csv) preserves all supplied sequences: thirteen guides are 20 nt, **SNHG12-sg3 and Neg sg2 are 21 nt**, and all have the supplied 75 nt scaffold. No guide is trimmed and no initiation G is inserted.

## Python and notebooks

```python
import pandas as pd
import probeharbour as ph

guides = pd.read_csv("examples/guides.csv", keep_default_na=False)
store = ph.design_guides(guides, chemistry="prime_5k")
summary = ph.guide_summary(store)
top = ph.top_guides(store, n=3)
top[["guide_id", "recalculated_score", "guide_overlap_nt", "rbd5", "rbd3"]]

ph.gpl.guide_overview(store, n=3)
ph.gpl.plot_guide_candidates(store, "HIF1A-sg1", n=3)
ph.gpl.write_guide_report(store, "guide_report.html", n=3)
```

Use the worked [notebook](examples/guide_design.ipynb) or its [Python equivalent](examples/guide_design.py) to inspect settings, run both chemistries and export results.

## Target definitions

Input columns are:

| Column | Meaning |
| --- | --- |
| `guide_id` | Unique label for this guide, preserved in results |
| `guide_sequence` | Supplied guide/protospacer, sense 5′→3′ |
| `scaffold_sequence` | Adjacent downstream sequence, sense 5′→3′ |
| `prefix_sequence` | Optional upstream construct/promoter/cloning sequence |
| `transcribed_leader` | Optional upstream sequence explicitly known to be in the RNA |

All sequences accept lowercase, whitespace and RNA U (normalised to DNA T). Ambiguous bases, missing required sequences and duplicate guide IDs raise errors. Distinct IDs with identical sequences remain separate targets and are detected as specificity conflicts.

**Default RNA mode** designs against `transcribed_leader + guide_sequence + scaffold_sequence`. The supplied promoter/ligation prefix is retained as metadata and excluded from the target. To target a confirmed transcript containing the upstream cassette, or the full construct, choose that explicitly:

```python
store = ph.design_guides(guides, target_mode="full_construct", boundary="both")
```

Full-construct mode uses exactly `prefix_sequence + guide_sequence + scaffold_sequence`. Do not combine that mode with a `transcribed_leader`. Available boundaries are `guide_scaffold`, `prefix_guide` and `both`; an RNA upstream boundary requires an explicit leader.

## What is scored and selected?

1. Slide the ligation cut around the requested boundary and retain only complete, in-bounds equal-arm windows that actually span it.
2. Require a configurable number of guide bases in the binder (default **15 nt**).
3. Calculate the inherited Primer3 thermodynamic and sequence properties and apply the selected ProbeHarbour ruleset.
4. Compare each full binder against every equal-length window in every other supplied target. Exclude exact collisions and candidates closer than the configured mismatch threshold (default **3** differences).
5. Rank eligible candidates separately for each guide, while retaining excluded candidates, score modifiers and rejection reasons.

The default arm lengths are **16 + 16 nt for Prime 5K** and **20 + 20 nt for Xenium v1**. Slides run from −5 to +5 nt. Out-of-bounds windows are skipped; missing RNA sequence is never padded with promoter sequence. Changing the arm length is supported, but does not recalibrate the chemistry score.

Guide overlap and mismatch thresholds are project heuristics, not experimental specificity guarantees. Chemistry modifiers reduce the score; they are not all hard exclusion filters. The default minimum score is 0; use `min_score` or `--min-score` to impose a score cutoff. A guide with no eligible candidates stays visible in the summary and report.

## Output files

| File | Contents |
| --- | --- |
| `report.html` | Offline interactive report with guide selector, positions, sequences, scores and explanations |
| `overview.png`, `overview.svg` | Overview of the top alternatives for all guides |
| `candidates.csv` | Every generated/scored candidate and its local eligibility evidence |
| `top_candidates.csv` | Up to the requested number of eligible alternatives per guide |
| `guide_summary.csv` | One row per input guide, including guides with no candidates |
| `targets.fasta` | Complete target sequences in sense orientation |
| `targets.csv` | Mapping of sanitised FASTA IDs back to original guide IDs |
| `manifest.json` | Input targets, configuration, thermodynamic parameters, package versions and evidence scope |

Existing outputs are protected unless `--overwrite` is supplied. CSV is for analysis/review; the exported sequences are not complete manufactured padlock probes. Current Xenium advanced guide-target requests use sense target FASTA; see the linked [10x instructions](https://www.10xgenomics.com/support/software/xenium-panel-designer/latest/tutorials/advanced-custom-designs).

## Reading the sequences and evidence

- `rbd5`, `rbd3` and `binder` follow ProbeHarbour's **target-sense** convention. They are not the physical antisense oligo sequence.
- `probe_arm_upstream_5to3` and `probe_arm_downstream_5to3` are separately named reverse complements of the target segments. `reverse_complement_binder` is the reverse complement of the full target window. None includes a vendor backbone or terminal chemistry.
- A scaffold arm may be shared by all guides; specificity is assessed for the adjacent pair, and shared individual arms are reported.
- “Top 3” means **overlapping alternatives**, not three independently validated or simultaneously compatible probes.
- The panel comparison is an exhaustive **ungapped sequence screen within the supplied guide targets**, not a transcriptome/genome screen, alignment with gaps, RNA structure prediction or assay validation. Negative control guides are still positive targets for detecting their constructs.
- Primer3 properties and inherited chemistry scores are useful prioritisation evidence. They do not establish accessibility of a structured or Cas-bound guide RNA.

See [design decisions and limitations](docs/guide_design.md) for sources, scoring details and target definitions.

## Development

```bash
python -m pip install -e '.[test,guide-ui]' build twine
python -m pytest -q
python -m build
python -m twine check dist/*
```

New guide modules are maintained directly as Python files. Existing `probeharbour` modules marked as generated retain their notebook sources in `nbs/`; regenerate those if changing the inherited engine. Use an isolated environment for this fork: it intentionally retains the `probeharbour` distribution/import name and replaces an upstream installation in the same environment.
