# %% [markdown]
# # 5. Rebuild and characterize published Xenium panels
#
# This tutorial shows the complete, auditable route from published 10x
# probe-location BED files to `ProbeStore` objects. It:
#
# - reconstructs each target sequence from the matching genome FASTA;
# - retains the source BED path and BED12 junction geometry;
# - optionally checks the coordinates against a matching GTF annotation;
# - calculates thermodynamic properties with Primer3; and
# - makes overview figures with the public `ps.pl` plotting API.
#
# BED, FASTA, and GTF files must describe the same genome assembly. ProbeHarbour
# does not download or silently substitute reference files.

# %%
from pathlib import Path
import os
import sys

import pandas as pd

# Prefer the checkout containing this tutorial over an older installed release.
# Set PROBEHARBOUR_REPO when launching Jupyter from outside the repository.
checkout_candidates = [
    Path(os.environ["PROBEHARBOUR_REPO"]).expanduser()
    if "PROBEHARBOUR_REPO" in os.environ else None,
    Path.cwd(),
    Path.cwd().parent,
    Path("/data/analysis/data_mcandrew/ProbeHarbour_Cobioda"),
]
CHECKOUT_ROOT = next(
    (
        path.resolve()
        for path in checkout_candidates
        if path is not None and (path / "probeharbour/facade.py").is_file()
    ),
    None,
)
if CHECKOUT_ROOT is not None:
    sys.path.insert(0, str(CHECKOUT_ROOT))

import probeharbour as ps

required_api = {"read_bed", "load_panel", "ProbeStore", "pl"}
missing_api = sorted(required_api - set(dir(ps)))
if missing_api:
    raise ImportError(
        f"Imported ProbeHarbour from {ps.__file__}, but it lacks {missing_api}. "
        "Restart the kernel after setting PROBEHARBOUR_REPO to the Cobioda checkout."
    )

print("Using ProbeHarbour:", Path(ps.__file__).resolve())

# %% [markdown]
# ## Configure the published panels and reference genomes
#
# Set the three environment variables below, or replace the relative defaults.
# V1 is distributed as several tissue BED files, so its registry entries are
# directories; all BEDs in each directory are combined into one store.

# %%
LOCAL_PANEL_ROOT = Path("/data/analysis/data_mcandrew/00_isoprobe/10x_panel")
LOCAL_MOUSE_REFERENCE = Path("/data/analysis/data_mcandrew/00_allos_manuscript/allos/data")

PANEL_ROOT = Path(os.environ.get(
    "PROBEHARBOUR_PANEL_ROOT",
    LOCAL_PANEL_ROOT if LOCAL_PANEL_ROOT.is_dir() else "reference/panels",
))
HUMAN_REFERENCE = Path(os.environ.get(
    "PROBEHARBOUR_HUMAN_REFERENCE",
    PANEL_ROOT / "Annotation" if (PANEL_ROOT / "Annotation").is_dir()
    else "reference/human_grch38",
))
MOUSE_REFERENCE = Path(os.environ.get(
    "PROBEHARBOUR_MOUSE_REFERENCE",
    LOCAL_MOUSE_REFERENCE if LOCAL_MOUSE_REFERENCE.is_dir()
    else "reference/mouse_grcm38",
))

MOUSE_FASTA = Path(os.environ.get(
    "PROBEHARBOUR_MOUSE_FASTA",
    MOUSE_REFERENCE / "Mus_musculus.GRCm38.dna.primary_assembly.fa",
))
HUMAN_FASTA = Path(os.environ.get(
    "PROBEHARBOUR_HUMAN_FASTA",
    HUMAN_REFERENCE / "Homo_sapiens.GRCh38.dna.primary_assembly.fa",
))


def one_bed(pattern):
    """Resolve one BED while tolerating browser-added suffixes such as ` (1)`."""
    matches = sorted(PANEL_ROOT.glob(pattern))
    return matches[:1]

PANEL_REGISTRY = {
    "prime_5k_mouse": {
        "species": "mouse",
        "chemistry": "prime_5k",
        "fasta": MOUSE_FASTA,
        "gtf": MOUSE_REFERENCE / "Mus_musculus.GRCm38.102.gtf.gz",
        "beds": one_bed("5K/xenium_mouse_5K_gene_expression_panel_probe_locations*.bed"),
    },
    "prime_5k_human": {
        "species": "human",
        "chemistry": "prime_5k",
        "fasta": HUMAN_FASTA,
        "gtf": HUMAN_REFERENCE / "Homo_sapiens.GRCh38.110.gtf.gz",
        "beds": one_bed("5K/xenium_human_5K_gene_expression_panel_probe_locations*.bed"),
    },
    "xenium_v1_mouse": {
        "species": "mouse",
        "chemistry": "xenium_v1",
        "fasta": MOUSE_FASTA,
        "gtf": MOUSE_REFERENCE / "Mus_musculus.GRCm38.102.gtf.gz",
        "beds": sorted((PANEL_ROOT / "V1/mouse").glob("*.bed")),
    },
    "xenium_v1_human": {
        "species": "human",
        "chemistry": "xenium_v1",
        "fasta": HUMAN_FASTA,
        "gtf": HUMAN_REFERENCE / "Homo_sapiens.GRCh38.98.gtf.gz",
        "beds": sorted((PANEL_ROOT / "V1/human").glob("*.bed")),
    },
}

# %%
input_manifest = pd.DataFrame([
    {
        "panel": name,
        "species": spec["species"],
        "chemistry": spec["chemistry"],
        "bed_files": len(spec["beds"]),
        "all_beds_exist": bool(spec["beds"]) and all(path.is_file() for path in spec["beds"]),
        "fasta_exists": spec["fasta"].is_file(),
        "gtf_exists": spec["gtf"].is_file(),
    }
    for name, spec in PANEL_REGISTRY.items()
])
input_manifest["ready"] = (
    input_manifest["all_beds_exist"] & input_manifest["fasta_exists"]
)
input_manifest

# %% [markdown]
# Inspect this table before continuing. Missing V1 BEDs usually mean the panel
# root has a different directory layout. A missing GTF only disables the
# optional transcript audit; BED loading and thermodynamics require BED+FASTA.

# %%
AVAILABLE_PANEL_NAMES = input_manifest.loc[input_manifest["ready"], "panel"].tolist()
if not AVAILABLE_PANEL_NAMES:
    unresolved = input_manifest.loc[
        ~input_manifest["ready"],
        ["panel", "all_beds_exist", "fasta_exists"],
    ]
    raise FileNotFoundError(
        "No complete BED+FASTA panel inputs were found. "
        "Set PROBEHARBOUR_PANEL_ROOT and the species FASTA variables.\n"
        f"{unresolved.to_string(index=False)}"
    )

print("Panels ready to load:", ", ".join(AVAILABLE_PANEL_NAMES))

# %% [markdown]
# ## One reusable loading call
#
# `ps.read_bed` reads every record and preserves the originating BED in
# `source_bed`. Passing a list is the supported way to merge the V1 tissue
# panels. The FASTA index (`.fai`) is created or reused by the underlying reader.

# %%
def load_existing_panel(panel_name):
    spec = PANEL_REGISTRY[panel_name]
    store = ps.read_bed(
        spec["beds"],
        reference_fasta=spec["fasta"],
        name=panel_name,
        panel_role="existing",
    )
    store.metadata.update({
        "species": spec["species"],
        "chemistry": spec["chemistry"],
        "annotation_gtf": str(spec["gtf"]),
    })
    return store


panels = {name: load_existing_panel(name) for name in AVAILABLE_PANEL_NAMES}
pd.concat([store.summary() for store in panels.values()], ignore_index=True)

# %% [markdown]
# The input audit below is intentionally simple. It makes duplicate records
# across merged V1 BEDs, chemistry-specific binder lengths, and BED-derived
# junction counts visible before annotation or thermodynamic calculation.

# %%
def panel_audit(store):
    obs = store.obs
    binder = obs["rbd5"].astype(str) + obs["rbd3"].astype(str)
    expected_binder_length = {
        "prime_5k": 32,  # two 16-nt arms
        "xenium_v1": 40,  # two 20-nt arms
    }[store.metadata["chemistry"]]
    return {
        "panel": store.name,
        "probes": len(obs),
        "genes": obs["gene_name"].nunique(dropna=True),
        "source_beds": obs["source_bed"].nunique(),
        "duplicate_records": int(
            obs.duplicated(["source_bed", "chrom", "chrom_start", "chrom_end", "probe_name"]).sum()
        ),
        "expected_binder_nt": expected_binder_length,
        "unexpected_binder_length": int(
            binder.str.len().ne(expected_binder_length).sum()
        ),
        "junction_spanning": int(obs["junction_spanning"].fillna(False).sum()),
    }


audit = pd.DataFrame(panel_audit(store) for store in panels.values())
audit

# %% [markdown]
# ## Optional transcript-annotation audit
#
# BED geometry does not require a GTF: a BED12 record with multiple blocks is
# directly junction-spanning. A matching GTF is useful for asking whether an
# annotated transcript supports those blocks. This step retains every input
# probe and gives unsupported rows an explicit `mapping_reason`.
#
# Run it for the species you need. Loading both transcriptomes can use
# substantial memory.

# %%
# from allos.transcript_data import TranscriptData
#
# mouse_td = TranscriptData(
#     gtf_file=PANEL_REGISTRY["prime_5k_mouse"]["gtf"],
#     reference_fasta=PANEL_REGISTRY["prime_5k_mouse"]["fasta"],
#     cache_dir="/tmp/probeharbour-allos-mouse",
# )
# panels["prime_5k_mouse"].annotate_transcripts(mouse_td)
# panels["prime_5k_mouse"].obs["mapping_status"].value_counts(dropna=False)

# %% [markdown]
# ## Thermodynamics
#
# `characterize()` adds Primer3 Tm, secondary-structure, GC, repeat, and
# parameter-ID columns in place. It never replaces a vendor or historic score.

# %%
for store in panels.values():
    store.characterize()

thermodynamic_summary = pd.concat(
    [
        store.obs[["tm_5p", "tm_3p", "tm_full"]]
        .describe()
        .loc[["count", "mean", "50%", "std", "min", "max"]]
        .assign(panel=name)
        .reset_index(names="statistic")
        for name, store in panels.items()
    ],
    ignore_index=True,
)
thermodynamic_summary

# %% [markdown]
# ## Plot with basic ProbeHarbour API calls
#
# The plotting functions accept a `ProbeStore` directly and return a Matplotlib
# figure. Keeping the figure object makes saving explicit and reproducible.

# %%
overview_figures = {
    name: ps.pl.panel_overview(store, title=name.replace("_", " ").title())
    for name, store in panels.items()
}

# %%
thermo_figures = {
    name: ps.pl.thermodynamics(store, title=f"{name.replace('_', ' ').title()} thermodynamics")
    for name, store in panels.items()
}

# %% [markdown]
# The same calls can save publication-independent PDFs. This tutorial writes
# one overview and one thermodynamics PDF per panel.

# %%
FIGURE_DIR = Path(os.environ.get("PROBEHARBOUR_FIGURE_DIR", "/tmp/probeharbour_existing_panels"))
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

for name, figure in overview_figures.items():
    figure.savefig(FIGURE_DIR / f"{name}_overview.pdf", bbox_inches="tight")
for name, figure in thermo_figures.items():
    figure.savefig(FIGURE_DIR / f"{name}_thermodynamics.pdf", bbox_inches="tight")

sorted(FIGURE_DIR.glob("*.pdf"))

# %% [markdown]
# ## Save reusable stores
#
# Saving after characterization preserves the calculated columns and provenance
# metadata, so later notebooks can start with `ps.load_panel(path)` without
# repeating sequence reconstruction or Primer3 calculations.

# %%
STORE_DIR = Path(os.environ.get("PROBEHARBOUR_STORE_DIR", "/tmp/probeharbour_existing_stores"))
saved_stores = {
    name: store.save(STORE_DIR / name)
    for name, store in panels.items()
}
saved_stores

# %%
example_name = next(iter(saved_stores))
reloaded = ps.load_panel(saved_stores[example_name])
assert len(reloaded.obs) == len(panels[example_name].obs)
assert {"tm_5p", "tm_3p", "tm_full", "junction_spanning"} <= set(reloaded.obs)
reloaded.summary()
