"""ProbeHarbour Streamlit dashboard: gene to reviewed, exportable panel."""
from __future__ import annotations

import gzip
import importlib
import inspect
import ast
import os
from pathlib import Path
import re
import sys

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

try:
    from allos.transcript_data import TranscriptData
except ImportError:
    TranscriptData = None

import probeharbour as ps

# Streamlit reruns this script while retaining already imported local modules.
# During development, refresh the plotting namespace only when the live process
# predates an API parameter used by this dashboard.
if "transcript_ids" not in inspect.signature(ps.pl.junction_map).parameters:
    ps.pl = importlib.reload(ps.pl)

import probeharbour.dashboard as dashboard_api
if getattr(dashboard_api, "CART_CSV_SCHEMA_VERSION", 0) < 2:
    dashboard_api = importlib.reload(dashboard_api)

from probeharbour.dashboard import (
    add_to_basket,
    candidate_review_table,
    export_panel_bundle,
    junction_summary,
    load_cart_csv,
    reannotate_junction_summary,
    remove_from_basket,
    selected_panel,
)

_logo_candidates = (
    Path(__file__).resolve().parent / "assets" / "ProbeHarbour_dashboard_logo.png",
    Path(__file__).resolve().parents[1] / "assets" / "ProbeHarbour_dashboard_logo.png",
)
DASHBOARD_LOGO = next((path for path in _logo_candidates if path.is_file()), None)

st.set_page_config(page_title="ProbeHarbour · workflow demo", page_icon="⚓", layout="wide")

# This is an intentionally separate UI experiment. The production dashboard
# remains in dashboard/app.py and is not altered by these styles.
st.markdown("""
<style>
  .stApp { background: #ffffff; color: #08224c; }
  [data-testid="stSidebar"] {
    background: #ffffff;
    border-right: 1px solid #d7e0e5;
  }
  [data-testid="stSidebar"] [data-testid="stMetric"] {
    background: #ffffff;
    border: 1px solid #d7e0e5;
    border-radius: .45rem;
    padding: .55rem .65rem;
  }
  div[data-testid="stVerticalBlockBorderWrapper"] {
    border-color: #d7e0e5;
    border-radius: .5rem;
  }
  .ph-kicker {
    color: #008f91;
    font-size: .74rem;
    font-weight: 700;
    letter-spacing: .09em;
    text-transform: uppercase;
    margin-bottom: -.25rem;
  }
  .ph-workflow {
    position: sticky;
    top: 2.8rem;
    z-index: 900;
    display: grid;
    grid-template-columns: repeat(5, minmax(0, 1fr));
    gap: .45rem;
    padding: .65rem;
    margin: .75rem 0 1.25rem 0;
    background: rgba(255,255,255,.96);
    border: 1px solid #d7e0e5;
    border-radius: .6rem;
    box-shadow: 0 5px 18px rgba(31,41,55,.06);
    backdrop-filter: blur(8px);
  }
  .ph-step {
    display: block;
    padding: .58rem .7rem;
    border: 1px solid #d7e0e5;
    border-radius: .4rem;
    color: #64748b !important;
    background: #f8fafc;
    font-size: .84rem;
    font-weight: 650;
    text-align: center;
    text-decoration: none !important;
  }
  .ph-step.done { color: #007f82 !important; background: #f0fbfb; border-color: #8bd1d1; }
  .ph-step.current { color: #ffffff !important; background: #008f91; border-color: #008f91; }
  .ph-step.future { color: #94a3b8 !important; cursor: default; }
  .ph-cart-title { color: #008f91; font-size: .75rem; font-weight: 750;
                   letter-spacing: .08em; text-transform: uppercase; }
  .ph-cart-actions { display: grid; grid-template-columns: 1fr 1fr; gap: .4rem; margin: .6rem 0 1rem; }
  .ph-cart-actions a { display: block; padding: .48rem; text-align: center; text-decoration: none !important;
                     color: #007f82 !important; background: white; border: 1px solid #8bc7c8;
                     border-radius: .4rem; font-size: .82rem; font-weight: 650; }
  .ph-cart-actions a:last-child { background: #008f91; color: white !important; }
</style>
""", unsafe_allow_html=True)

if DASHBOARD_LOGO is not None:
    st.image(str(DASHBOARD_LOGO), width=460)
else:
    st.title("ProbeHarbour")
st.caption("Design, evaluate and export isoform-aware Xenium probe panels.")

_basket_for_stage = st.session_state.get("probe_basket")
_has_basket = _basket_for_stage is not None and not _basket_for_stage.obs.empty
_current_stage = (
    5 if st.session_state.get("prepared_export") is not None else
    4 if _has_basket else
    3 if st.session_state.get("scored_pool") is not None else
    2 if st.session_state.get("design_pool") is not None else 1
)
_steps = ((1, "Gene", "gene"), (2, "Junction", "junction"),
          (3, "Candidate", "candidate"), (4, "Cart", "cart"),
          (5, "Export", "export"))
_step_html = []
for _number, _label, _anchor in _steps:
    _state = "current" if _number == _current_stage else ("done" if _number < _current_stage else "future")
    _text = f"✓ {_label}" if _number < _current_stage else f"{_number} {_label}"
    if _number <= _current_stage:
        _step_html.append(f'<a class="ph-step {_state}" href="#{_anchor}">{_text}</a>')
    else:
        _step_html.append(f'<span class="ph-step {_state}">{_text}</span>')
st.markdown('<div class="ph-workflow">' + ''.join(_step_html) + '</div>', unsafe_allow_html=True)

with st.expander("New to ProbeHarbour? Start here", expanded=True):
    st.markdown(
        """
        Build an add-on panel one gene at a time:

        1. **Choose a gene** to load its annotated transcripts.
        2. **Choose junctions** that distinguish the isoforms you want to detect.
        3. **Review candidates** scored under the selected Xenium chemistry.
        4. **Add confirmed probes to the cart** and repeat for more genes.
        5. **Export once** when the complete panel has been reviewed.

        The Snap25 example is preloaded so you can explore the workflow immediately.
        Nothing is added to the panel until you press **Add selected probes to design cart**.
        """
    )


DEFAULT_GTF = Path(os.environ.get("PROBEHARBOUR_GTF", ""))
DEFAULT_FASTA = Path(os.environ.get("PROBEHARBOUR_FASTA", ""))
DEFAULT_BED = Path(os.environ.get("PROBEHARBOUR_EXISTING_BED", ""))
DEFAULT_BLAST_DB = Path(os.environ.get("PROBEHARBOUR_BLAST_DB", ""))
_bundled_csv = Path(__file__).resolve().parent / "data" / "96ISOS.csv"
if not _bundled_csv.is_file():
    _bundled_csv = Path(__file__).resolve().parents[1] / "tutorials" / "data" / "96ISOS.csv"
DEFAULT_ADDON_CSV = Path(os.environ.get("PROBEHARBOUR_ADDON_CSV", str(_bundled_csv)))


@st.cache_resource(show_spinner="Loading transcript annotation…")
def load_transcripts(gtf: str, fasta: str):
    if TranscriptData is None:
        raise ImportError(
            "Transcript-aware design requires the Cobioda Allos implementation "
            "that provides allos.transcript_data.TranscriptData."
        )
    return TranscriptData(
        gtf_file=Path(gtf), reference_fasta=Path(fasta),
        cache_dir="/tmp/probestore-dashboard-allos-cache",
    )


@st.cache_data(show_spinner="Indexing gene names…")
def gene_names(gtf: str):
    opener = gzip.open if str(gtf).endswith(".gz") else open
    names = set()
    pattern = re.compile(r'gene_name "([^"]+)"')
    with opener(gtf, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            match = pattern.search(line)
            if match:
                names.add(match.group(1))
    return sorted(names, key=str.casefold)


@st.cache_resource(show_spinner="Loading existing panel…")
def load_existing_panel(bed: str, fasta: str):
    return ps.read_bed(
        Path(bed), reference_fasta=Path(fasta), name="existing_10x_panel",
        panel_role="existing",
    )


def evaluate_addon_csv(source, *, name: str, chemistry: str, gtf: str, fasta: str,
                       existing_bed: str | None = None):
    """Normalize, map and evaluate an existing add-on CSV for cart use."""
    imported = load_cart_csv(source, name=name)
    has_nested_mapping = (
        "coordinates" in imported.obs
        and imported.obs["coordinates"].map(lambda value: isinstance(value, dict)).any()
    )
    has_flat_mapping = {
        "rbd5_chrom_start", "rbd5_chrom_end", "rbd3_chrom_start", "rbd3_chrom_end",
    } <= set(imported.obs)
    if not has_nested_mapping and not has_flat_mapping:
        if not Path(gtf).is_file() or not Path(fasta).is_file():
            raise FileNotFoundError("Reference GTF and FASTA are required to map sequence-only probes")
        imported.map_probes(load_transcripts(gtf, fasta))
    if "tm_full" not in imported.obs:
        imported.characterize()
    if "probeharbour_score" not in imported.obs:
        imported.annotate_scores(chemistry, score_column="probeharbour_score")
    if existing_bed and Path(existing_bed).is_file():
        imported.check_panel_compatibility(
            load_existing_panel(existing_bed, fasta), calculate_heterodimer=False,
        )
    return imported


def existing_panel_compatibility(frame: pd.DataFrame) -> pd.Series:
    """Return compatibility evidence specifically against the existing 10x panel."""
    if "existing_panel_compatibility_pass" in frame:
        return frame["existing_panel_compatibility_pass"]
    return pd.Series(pd.NA, index=frame.index, dtype="boolean")


def existing_panel_conflicts(frame: pd.DataFrame) -> pd.DataFrame:
    """Return one readable row per probe that fails the existing-panel audit."""
    compatibility = existing_panel_compatibility(frame)
    failed = frame.loc[compatibility.eq(False)].copy()
    if failed.empty:
        return pd.DataFrame()

    reason_labels = {
        "duplicate_binder": "duplicate full binder",
        "reverse_complement_binder": "reverse-complement full binder",
        "duplicate_arm": "duplicate RBD arm",
        "long_complementary_run": "complementary run ≥ 14 nt",
        "strong_panel_heterodimer": "heterodimer ΔG < −10.5 kcal/mol",
        "target_already_present": "target gene already represented",
        "existing_panel_spacing": "existing 10x probe < 75 bp away",
    }

    def reasons(value) -> str:
        if isinstance(value, str):
            try:
                parsed = ast.literal_eval(value)
                value = parsed if isinstance(parsed, (list, tuple, set)) else [value]
            except (SyntaxError, ValueError):
                value = [part.strip() for part in value.split(",") if part.strip()]
        if not isinstance(value, (list, tuple, set)):
            value = [] if pd.isna(value) else [value]
        return "; ".join(reason_labels.get(str(item), str(item)) for item in value)

    def first(*columns, default=pd.NA):
        for column in columns:
            if column in failed:
                return failed[column]
        return pd.Series(default, index=failed.index)

    def coalesce(*columns):
        result = pd.Series(pd.NA, index=failed.index, dtype="object")
        for column in columns:
            if column in failed:
                result = result.fillna(failed[column])
        return result

    result = pd.DataFrame(index=failed.index)
    result["Probe"] = first("candidate_id", "probe_name", "panel_record_id")
    result["Gene"] = first("gene_name", "gene", "gene_id")
    result["Triggered check(s)"] = first(
        "existing_panel_compatibility_reasons",
        "panel_compatibility_reasons",
    ).map(reasons)
    result["Conflicting existing 10x probe"] = coalesce(
        "existing_panel_duplicate_probe",
        "existing_panel_duplicate_arm_probe",
        "existing_panel_worst_heterodimer_probe",
        "existing_panel_closest_probe",
    )
    result["Nearest existing 10x probe"] = first(
        "nearest_existing_panel_probe_name"
    )
    result["Distance (bp)"] = pd.to_numeric(
        first("nearest_existing_panel_probe_distance"), errors="coerce"
    )
    result["Longest complementary run (nt)"] = pd.to_numeric(
        first("existing_panel_longest_complement"),
        errors="coerce",
    )
    result["Worst heterodimer ΔG (kcal/mol)"] = pd.to_numeric(
        first("existing_panel_worst_heterodimer_dg"),
        errors="coerce",
    )
    return result.reset_index(drop=True)


with st.sidebar:
    st.markdown('<div class="ph-cart-title">Working ProbeSet</div>', unsafe_allow_html=True)
    demo_basket = st.session_state.get("probe_basket")
    if demo_basket is None or demo_basket.obs.empty:
        st.markdown("#### Design cart")
        st.caption("Your cart is empty. Confirm a candidate to begin building the panel.")
    else:
        demo_gene_values = pd.Series(pd.NA, index=demo_basket.obs.index, dtype="object")
        for demo_gene_column in ("gene_name", "gene", "gene_id"):
            if demo_gene_column in demo_basket.obs:
                demo_gene_values = demo_gene_values.fillna(demo_basket.obs[demo_gene_column])
        demo_compatibility = existing_panel_compatibility(demo_basket.obs)
        demo_metrics = st.columns(3)
        demo_metrics[0].metric("Probes", len(demo_basket.obs))
        demo_metrics[1].metric("Genes", demo_gene_values.dropna().astype(str).nunique())
        demo_metrics[2].metric(
            "10x conflicts", int(demo_compatibility.eq(False).sum())
        )
        demo_active_gene = st.session_state.get("dashboard_gene")
        if demo_active_gene:
            demo_gene_probes = int(demo_gene_values.astype(str).eq(str(demo_active_gene)).sum())
            st.caption(f"Current gene · {demo_active_gene}: {demo_gene_probes} probe(s) in cart")
        st.markdown(
            '<div class="ph-cart-actions"><a href="#cart">Review cart</a>'
            '<a href="#export">Checkout</a></div>', unsafe_allow_html=True,
        )
    st.divider()
    st.header("Design configuration")
    chemistry = st.selectbox("Xenium chemistry", ["prime_5k", "xenium_v1"],
                             format_func=lambda value: {"prime_5k": "Prime 5K", "xenium_v1": "Xenium V1"}[value])
    arm_length = st.number_input("RBD arm length (nt)", 8, 40,
                                 16 if chemistry == "prime_5k" else 20)
    slide = st.number_input("Maximum slide (±nt)", 0, 20, 6)
    st.divider()
    with st.expander("Reference files", expanded=False):
        gtf_path = st.text_input("GENCODE/Ensembl GTF", str(DEFAULT_GTF))
        fasta_path = st.text_input("Reference FASTA", str(DEFAULT_FASTA))
        existing_bed = st.text_input("Existing 10x BED (optional)", str(DEFAULT_BED))
        blast_database = st.text_input("cDNA BLAST database (optional)", str(DEFAULT_BLAST_DB))
    run_blast = st.checkbox("Run transcriptome BLAST", value=False,
                            help="Recommended before final export; requires a local BLAST database.")
    st.divider()
    st.header("Design cart")
    panel_name = st.text_input("ProbeSet name", "xenium_add_on_panel", key="basket_name")
    st.caption("Start from the bundled 96-probe isoform add-on panel, or build a new cart below.")
    with st.expander("What is the 96-probe example?", expanded=False):
        st.markdown(
            """
            This is an existing mouse isoform add-on panel supplied as RBD5/RBD3
            sequences. ProbeHarbour imports all 96 records, maps their complete
            binders to the configured transcript annotation, calculates the
            selected chemistry score, and checks them against the existing 10x
            panel and one another.

            **To inspect it:** load the panel, open **Review design cart**, choose
            a gene under **View probes for**, and enable **Show this gene on the
            transcript map**. The map displays the add-on probes together with
            the supporting transcripts and, when an existing BED is configured,
            the existing 10x probes.
            """
        )
    if st.button(
        "🛒 Load the 96-probe add-on panel",
        disabled=not DEFAULT_ADDON_CSV.is_file(),
        help=str(DEFAULT_ADDON_CSV),
    ):
        try:
            with st.spinner("Loading, mapping and evaluating the 96-probe add-on panel…"):
                imported_cart = evaluate_addon_csv(
                    DEFAULT_ADDON_CSV, name=panel_name, chemistry=chemistry,
                    gtf=gtf_path, fasta=fasta_path, existing_bed=existing_bed,
                )
                current_cart = st.session_state.get("probe_basket")
                loaded_cart = add_to_basket(current_cart, imported_cart, name=panel_name)
                loaded_cart.check_internal_compatibility()
            st.session_state.probe_basket = loaded_cart
            st.session_state.pop("prepared_export", None)
            st.session_state.open_cart_next = True
            st.session_state.cart_notice = (
                f"Loaded {len(imported_cart.obs)} probes from the existing add-on panel; "
                f"the cart now contains {len(loaded_cart.obs)} probes."
            )
            st.rerun()
        except Exception as error:
            st.error(f"Could not load the 96-probe add-on panel: {error}")
    with st.expander("Load a different add-on CSV", expanded=False):
        uploaded_cart = st.file_uploader(
            "Add-on panel CSV", type=["csv"], key="cart_csv_upload",
            help=("Accepts ProbeHarbour exports and sequence-only add-on files using "
                  "name, rbdm5 and rbdm3 columns."),
        )
        cart_import_mode = st.radio(
            "Import action", ["Merge into cart", "Replace cart"], horizontal=True,
            disabled=uploaded_cart is None,
        )
        if st.button("Load uploaded CSV", disabled=uploaded_cart is None):
            try:
                with st.spinner("Loading, mapping and evaluating the uploaded add-on panel…"):
                    uploaded_cart.seek(0)
                    imported_cart = evaluate_addon_csv(
                        uploaded_cart, name=panel_name, chemistry=chemistry,
                        gtf=gtf_path, fasta=fasta_path, existing_bed=existing_bed,
                    )
                    current_cart = (st.session_state.get("probe_basket")
                                    if cart_import_mode == "Merge into cart" else None)
                    loaded_cart = add_to_basket(current_cart, imported_cart, name=panel_name)
                    loaded_cart.check_internal_compatibility()
                st.session_state.probe_basket = loaded_cart
                st.session_state.pop("prepared_export", None)
                st.session_state.open_cart_next = True
                st.session_state.cart_notice = (
                    f"Loaded {len(imported_cart.obs)} probes from the uploaded CSV; "
                    f"the cart now contains {len(loaded_cart.obs)} probes."
                )
                st.rerun()
            except Exception as error:
                st.error(f"Could not load the uploaded add-on CSV: {error}")
    sidebar_basket = st.session_state.get("probe_basket")
    if sidebar_basket is None or sidebar_basket.obs.empty:
        st.caption("Empty · add reviewed candidates below to begin building your ProbeSet.")
    else:
        st.metric("Probes in cart", len(sidebar_basket.obs))
    if "show_design_cart" not in st.session_state:
        st.session_state.show_design_cart = False
    if st.session_state.pop("open_cart_next", False):
        st.session_state.show_design_cart = True
    show_basket = st.toggle(
        "Review design cart", key="show_design_cart",
        disabled=sidebar_basket is None or sidebar_basket.obs.empty,
        help="Open the persistent working cart above the design workflow.",
    )


cart_notice = st.session_state.pop("cart_notice", None)
if cart_notice:
    st.toast(cart_notice, icon="🛒")
if sidebar_basket is not None and not sidebar_basket.obs.empty:
    with st.expander(
        f"Design cart · {len(sidebar_basket.obs)} probes",
        expanded=show_basket,
    ):
        basket_gene_values = pd.Series(pd.NA, index=sidebar_basket.obs.index, dtype="object")
        for basket_gene_column in ("gene_name", "gene", "gene_id"):
            if basket_gene_column in sidebar_basket.obs:
                basket_gene_values = basket_gene_values.fillna(
                    sidebar_basket.obs[basket_gene_column]
                )
        basket_gene_count = basket_gene_values.dropna().astype(str).nunique()
        compatibility = existing_panel_compatibility(sidebar_basket.obs)
        flagged_count = int(compatibility.eq(False).sum())
        basket_metrics = st.columns(4)
        basket_metrics[0].metric("Probes", len(sidebar_basket.obs))
        basket_metrics[1].metric("Genes", basket_gene_count)
        basket_metrics[2].metric(
            "Junctions",
            sidebar_basket.obs["junction_id"].nunique()
            if "junction_id" in sidebar_basket.obs else 0,
        )
        basket_metrics[3].metric("Existing 10x conflicts", flagged_count)
        basket_preview = sidebar_basket.obs.copy()
        basket_preview["gene"] = basket_gene_values
        st.caption(
            "Choose a gene below to inspect its probe records. Turn on the "
            "transcript map to verify where the add-on probes bind and how they "
            "sit relative to the existing 10x panel."
        )
        if flagged_count:
            st.warning(
                f"{flagged_count} probe(s) conflict with the configured existing "
                "10x panel. Review those records before final export."
            )
            conflict_details = existing_panel_conflicts(basket_preview)
            with st.expander(
                f"Show exact existing-panel conflicts ({len(conflict_details)})",
                expanded=True,
            ):
                st.caption(
                    "A probe fails when any listed check crosses its threshold. "
                    "Distance is the gap between physically bound RBD blocks; "
                    "overlapping blocks therefore have a distance of 0 bp."
                )
                st.dataframe(
                    conflict_details,
                    hide_index=True,
                    use_container_width=True,
                    column_config={
                        "Distance (bp)": st.column_config.NumberColumn(format="%d"),
                        "Longest complementary run (nt)": st.column_config.NumberColumn(
                            format="%d"
                        ),
                        "Worst heterodimer ΔG (kcal/mol)": st.column_config.NumberColumn(
                            format="%.2f"
                        ),
                    },
                )
                sequence_conflicts = conflict_details.loc[
                    conflict_details["Triggered check(s)"].astype(str).str.contains(
                        "binder|RBD arm|complementary|heterodimer",
                        case=False, regex=True,
                    )
                    & conflict_details["Conflicting existing 10x probe"].notna()
                ].reset_index(drop=True)
                if not sequence_conflicts.empty and existing_bed and Path(existing_bed).is_file():
                    show_alignment = st.toggle(
                        "Inspect a sequence conflict",
                        value=False,
                        help=(
                            "Align the selected add-on binder to the reverse complement "
                            "of the existing 10x probe that triggered the sequence screen."
                        ),
                    )
                    if show_alignment:
                        alignment_index = st.selectbox(
                            "Probe pair",
                            range(len(sequence_conflicts)),
                            format_func=lambda index: (
                                f"{sequence_conflicts.iloc[index]['Probe']} ↔ "
                                f"{sequence_conflicts.iloc[index]['Conflicting existing 10x probe']}"
                            ),
                        )
                        alignment_record = sequence_conflicts.iloc[alignment_index]
                        query_identifier = str(alignment_record["Probe"])
                        reference_identifier = str(
                            alignment_record["Conflicting existing 10x probe"]
                        )
                        reference_panel = load_existing_panel(existing_bed, fasta_path)

                        def dashboard_probe_row(frame, identifier):
                            matches = pd.Series(False, index=frame.index)
                            for identifier_column in (
                                "candidate_id", "probe_name", "name", "panel_record_id",
                            ):
                                if identifier_column in frame:
                                    matches |= (
                                        frame[identifier_column].notna()
                                        & frame[identifier_column].astype(str).eq(identifier)
                                    )
                            return frame.loc[matches].iloc[0] if matches.any() else None

                        query_row = dashboard_probe_row(
                            sidebar_basket.obs, query_identifier
                        )
                        reference_row = dashboard_probe_row(
                            reference_panel.obs, reference_identifier
                        )
                        if query_row is None or reference_row is None:
                            st.warning(
                                "The selected probe pair could not be recovered from "
                                "the current cart and existing-panel records."
                            )
                        else:
                            alignment_figure = ps.pl.panel_conflict_alignment(
                                query_row,
                                reference_row,
                                threshold=14,
                                title=(
                                    f"{query_identifier}: existing-panel "
                                    "complementarity audit"
                                ),
                            )
                            st.pyplot(alignment_figure, use_container_width=True)
                            plt.close(alignment_figure)
        elif not compatibility.notna().any():
            st.info(
                "Compatibility against an existing 10x panel has not been evaluated."
            )
        else:
            st.success("No conflicts with the configured existing 10x panel were detected.")

        gene_summary = basket_preview.groupby("gene", dropna=False).agg(
            probes=("candidate_id", "size"),
            junctions=("junction_id", "nunique"),
            median_score=("probeharbour_score", "median"),
        )
        basket_preview["existing_10x_compatible"] = compatibility
        if compatibility.notna().any():
            gene_summary["existing_10x_conflicts"] = basket_preview.groupby(
                "gene", dropna=False
            )["existing_10x_compatible"].apply(
                lambda values: int(values.eq(False).sum())
            )
        st.markdown("**ProbeSet summary by gene**")
        st.dataframe(
            gene_summary.reset_index(), hide_index=True, use_container_width=True,
            column_config={"median_score": st.column_config.NumberColumn(format="%.4f")},
        )

        cart_gene_options = ["All genes"] + sorted(
            basket_preview["gene"].dropna().astype(str).unique(), key=str.casefold
        )
        cart_gene = st.selectbox(
            "View probes for", cart_gene_options,
            key="cart_gene_filter_" + str(abs(hash(tuple(cart_gene_options)))),
            help="Focus the cart on one gene without changing or removing other probes.",
        )
        focused_cart = (basket_preview if cart_gene == "All genes" else
                        basket_preview.loc[basket_preview["gene"].astype(str).eq(cart_gene)])
        show_cart_map = st.checkbox(
            "Show this gene on the transcript map",
            value=False,
            disabled=cart_gene == "All genes",
            key="show_cart_gene_map_" + str(abs(hash(tuple(cart_gene_options)))),
        )
        if show_cart_map and cart_gene != "All genes":
            if Path(gtf_path).is_file() and Path(fasta_path).is_file():
                cart_td = load_transcripts(gtf_path, fasta_path)
                cart_display = sidebar_basket
                if existing_bed and Path(existing_bed).is_file():
                    cart_display = ps.ProbeStore.merge(
                        [load_existing_panel(existing_bed, fasta_path), sidebar_basket],
                        name=f"existing_plus_{sidebar_basket.name}",
                        roles=["existing", "add_on"],
                    )
                cart_figure = ps.pl.gene_map(
                    cart_gene, cart_td, cart_display,
                    title=f"{cart_gene}: existing 10x and add-on probes in cart",
                )
                st.pyplot(cart_figure, use_container_width=True)
                st.caption(
                    "Existing 10x probes retain their BED geometry: purple junction "
                    "arcs span exon blocks, while purple diamonds mark exon-contained probes."
                )
                plt.close(cart_figure)
            else:
                st.warning("A valid GTF and FASTA are required for the transcript map.")

        cart_columns = [
            "gene", "candidate_id", "junction_id", "shift_nt", "probeharbour_score",
            "existing_10x_compatible",
        ]
        cart_items = focused_cart[
            [column for column in cart_columns if column in focused_cart]
        ].copy()
        cart_items.insert(0, "remove", False)
        cart_key = "cart_editor_" + str(abs(hash((cart_gene, *tuple(
            sidebar_basket.obs["candidate_id"].astype(str)
        )))))
        st.markdown("**Probes in this cart**")
        edited_cart = st.data_editor(
            cart_items, hide_index=True, use_container_width=True, key=cart_key,
            disabled=[column for column in cart_items if column != "remove"],
            column_config={
                "remove": st.column_config.CheckboxColumn("Remove"),
                "probeharbour_score": st.column_config.NumberColumn("Score", format="%.4f"),
                "existing_10x_compatible": st.column_config.CheckboxColumn(
                    "Compatible with existing 10x"
                ),
            },
        )
        cart_removals = edited_cart.loc[edited_cart["remove"], "candidate_id"].astype(str).tolist()
        if st.button("Remove selected probes from cart", disabled=not cart_removals):
            reduced = remove_from_basket(sidebar_basket, cart_removals, name=panel_name)
            if reduced.obs.empty:
                st.session_state.pop("probe_basket", None)
                remaining_count = 0
            else:
                reduced.check_internal_compatibility()
                st.session_state.probe_basket = reduced
                remaining_count = len(reduced.obs)
            st.session_state.pop("prepared_export", None)
            st.session_state.cart_notice = (
                f"Removed {len(cart_removals)} probe(s); {remaining_count} remain in the cart."
            )
            st.rerun()
        st.caption(
            "The cart is a working design space. Add, inspect and remove probes here; "
            "export is a separate finalization step."
        )


missing_reference = [path for path in (gtf_path, fasta_path) if not Path(path).is_file()]
if missing_reference:
    st.error("Reference files not found: " + ", ".join(missing_reference))
    st.stop()

genes = gene_names(gtf_path)

def reset_gene_workspace():
    """Clear gene-specific views while preserving the multi-gene design cart."""
    selected_gene = st.session_state.get("gene_picker")
    if selected_gene == st.session_state.get("dashboard_gene"):
        return
    for state_key in (
        "design_pool", "scored_pool", "confirmed_panel", "prepared_export",
        "open_cart_next",
    ):
        st.session_state.pop(state_key, None)
    gene_widget_prefixes = (
        "visible_transcripts_", "renumber_junctions_", "junction_class_filter_",
        "junction_selection_", "junction_inspector_", "junction_annotation_basis_",
        "use_visible_annotation_", "candidate_editor_", "scored_junctions_",
    )
    for state_key in list(st.session_state):
        if state_key.startswith(gene_widget_prefixes):
            st.session_state.pop(state_key, None)
    st.session_state.dashboard_gene = selected_gene
    st.session_state.gene_reset_notice = (
        f"Ready to begin {selected_gene}. The previous gene workspace was cleared; "
        "probes already added to the design cart were preserved."
    )


st.markdown('<div id="gene"></div>', unsafe_allow_html=True)
st.subheader("1 · Choose a gene")
st.caption(
    "Start with one gene. ProbeHarbour reads its annotated transcript models and "
    "enumerates every exon–exon junction that can be targeted."
)
demo_gene = "Snap25" if "Snap25" in genes else None
demo_index = genes.index(demo_gene) if demo_gene is not None else None
gene = st.selectbox(
    "Gene", genes, index=demo_index, placeholder="Type to search the annotation…",
    help="Snap25 is loaded initially as an isoform-junction design demonstration.",
    key="gene_picker", on_change=reset_gene_workspace,
)
if gene_reset_notice := st.session_state.pop("gene_reset_notice", None):
    st.info(gene_reset_notice)
enumerate_clicked = st.button("Enumerate annotated junctions", type="primary", disabled=gene is None)
load_initial_demo = "design_pool" not in st.session_state and gene == demo_gene
if enumerate_clicked or load_initial_demo:
    with st.spinner(f"Enumerating {gene} junction candidates…"):
        td = load_transcripts(gtf_path, fasta_path)
        st.session_state.design_pool = ps.design_gene(
            gene, td, chemistry=chemistry, arm_length=int(arm_length),
            slides=range(-int(slide), int(slide) + 1), dataset_id="streamlit-dashboard",
        )
        st.session_state.dashboard_gene = gene
        st.session_state.pop("scored_pool", None)
        st.session_state.pop("confirmed_panel", None)

design_pool = st.session_state.get("design_pool")
if design_pool is None:
    st.info("Select a gene to begin. ProbeHarbour will show its annotated exon–exon junctions before designing a panel.")
    st.stop()

td = load_transcripts(gtf_path, fasta_path)
active_gene = st.session_state.get("dashboard_gene", gene)
summary = junction_summary(design_pool)
st.markdown('<div id="junction"></div>', unsafe_allow_html=True)
st.subheader("2 · Select biological junctions")
st.caption(
    "First decide which isoform event matters biologically. The plot shows each "
    "transcript as exons joined by coloured junction arcs; the controls below only "
    "change what you see and select."
)

st.caption(
    "Arc colours: green = specific, orange = partially specific, grey = "
    "constitutive. J-labels are short display names; stable genomic IDs are "
    "retained for design and export."
)

annotated_transcript_ids = list(td.get_transcripts_by_gene_name(active_gene))
transcript_names = {
    str(transcript_id): str(
        td.get_transcript_info(transcript_id).get("transcript_name") or transcript_id
    )
    for transcript_id in annotated_transcript_ids
}
annotation_basis_key = f"junction_annotation_basis_{active_gene}"
visible_transcript_ids = st.multiselect(
    "Transcripts to compare",
    [str(value) for value in annotated_transcript_ids],
    default=[str(value) for value in annotated_transcript_ids],
    format_func=lambda value: transcript_names.get(str(value), str(value)),
    key=f"visible_transcripts_{active_gene}",
    help="Remove transcripts that are not relevant to the isoform event.",
)
reannotation_enabled = st.toggle(
    "Reclassify junctions using only these transcripts",
    value=False,
    key=f"use_visible_annotation_{active_gene}",
    disabled=not visible_transcript_ids,
    help=(
        "Off: transcript removal simplifies only the plot. On: support counts "
        "and junction classes are recalculated from the transcripts above."
    ),
)
previous_annotation_basis = st.session_state.get(annotation_basis_key)
annotation_basis = (
    list(visible_transcript_ids)
    if reannotation_enabled and visible_transcript_ids else None
)
if annotation_basis is None:
    st.session_state.pop(annotation_basis_key, None)
else:
    st.session_state[annotation_basis_key] = annotation_basis
if previous_annotation_basis != annotation_basis:
    for state_key in ("scored_pool", "confirmed_panel", "prepared_export"):
        st.session_state.pop(state_key, None)
    st.session_state.pop(f"scored_junctions_{active_gene}", None)
st.caption(
    f"Showing {len(visible_transcript_ids)} of {len(annotated_transcript_ids)} transcripts · "
    + (
        "junction classes recalculated from this set"
        if annotation_basis is not None else "junction classes use the full annotation"
    )
)
working_summary = (
    reannotate_junction_summary(summary, annotation_basis)
    if annotation_basis is not None else summary.copy()
)
renumber_visible = False
visible_set = set(visible_transcript_ids)

def supporting_transcript_ids(value):
    if isinstance(value, str):
        try:
            import ast
            value = ast.literal_eval(value)
        except (SyntaxError, ValueError):
            value = [value]
    if not isinstance(value, (list, tuple, set)):
        return []
    return [str(item) for item in value]

def is_visible_junction(value):
    supporting = supporting_transcript_ids(value)
    return True if not supporting else bool(visible_set.intersection(supporting))

if "supporting_transcripts" in working_summary:
    visible_summary = working_summary.loc[
        working_summary["supporting_transcripts"].map(is_visible_junction)
    ].copy().reset_index(drop=True)
else:
    visible_summary = working_summary.copy()
if renumber_visible:
    visible_summary["junction_label"] = [f"J{index + 1}" for index in range(len(visible_summary))]

junction_lookup = visible_summary.set_index("junction_id")
class_filter = st.segmented_control(
    "Show junction classes",
    ["Isoform-discriminating", "Specific", "Partially specific", "All"],
    default="Isoform-discriminating",
    key=f"junction_class_filter_{active_gene}",
    help=("Filters the junction choices. Existing selections "
          "remain selected when the filter changes."),
)

def junction_matches_class(junction_id):
    junction_class = str(junction_lookup.loc[junction_id].get(
        "junction_class", "unclassified"
    )).lower().replace("_", " ")
    if class_filter == "All":
        return True
    if class_filter == "Specific":
        return junction_class == "specific"
    if class_filter == "Partially specific":
        return junction_class == "partially specific"
    return junction_class in {"specific", "partially specific"}

visible_junction_options = visible_summary["junction_id"].astype(str).tolist()
filtered_junction_options = [
    value for value in visible_junction_options if junction_matches_class(value)
]
junction_key = f"junction_selection_{active_gene}"
if junction_key in st.session_state:
    st.session_state[junction_key] = [
        value for value in st.session_state[junction_key] if value in visible_junction_options
    ]

if active_gene == "Snap25":
    suggested_snap25 = visible_summary.loc[
        visible_summary["junction_label"].astype(str).isin(["J4", "J5"]),
        "junction_id",
    ].astype(str).tolist()
    if st.button(
        "Use suggested Snap25 targets", key=f"suggest_snap25_{active_gene}",
        disabled=not suggested_snap25,
        help="Selects the alternative exon-5 junctions used in this demonstration.",
    ):
        st.session_state[junction_key] = suggested_snap25
        st.rerun()

retained_selections = [
    value for value in st.session_state.get(junction_key, [])
    if value in visible_junction_options
]
junction_options = [
    value for value in visible_junction_options
    if value in filtered_junction_options or value in retained_selections
]

def junction_option_label(value):
    row = junction_lookup.loc[value]
    label = row.get("junction_label", "junction")
    junction_class = str(row.get("junction_class", "unclassified")).replace("_", " ")
    supporting = row.get("transcripts_with_junction")
    annotated = row.get("annotated_transcripts")
    coverage = (f" · {int(supporting)}/{int(annotated)} transcripts"
                if pd.notna(supporting) and pd.notna(annotated) else "")
    return f"{label} · {junction_class} · {value}{coverage}"

junction_widget_options = {"key": junction_key}
if junction_key not in st.session_state:
    junction_widget_options["default"] = []
selected_junctions = st.multiselect(
    "Junctions to design probes for", junction_options, **junction_widget_options,
    format_func=junction_option_label,
    placeholder="Search by J-number or stable genomic junction ID…",
    help=("Type a junction ID such as J4 or a stable genomic ID such as "
          "2:136763621-136770056:+. Class and transcript details remain visible "
          "in each result."),
)

scored_junctions_key = f"scored_junctions_{active_gene}"
selected_junction_signature = tuple(sorted(str(value) for value in selected_junctions))
previously_scored_junctions = tuple(
    st.session_state.get(scored_junctions_key, ())
)
selection_changed_after_design = (
    "scored_pool" in st.session_state
    and previously_scored_junctions != selected_junction_signature
)
if selection_changed_after_design:
    for state_key in ("scored_pool", "confirmed_panel", "prepared_export"):
        st.session_state.pop(state_key, None)
    st.session_state.pop(scored_junctions_key, None)

if selected_junctions:
    selected_labels = [str(junction_lookup.loc[value].get("junction_label", value))
                       for value in selected_junctions]
    st.success(
        f"{len(selected_junctions)} target junction(s) selected: "
        + ", ".join(selected_labels)
    )
    inspected_junction = st.selectbox(
        "Inspect selected junction",
        selected_junctions,
        format_func=junction_option_label,
        key=f"junction_inspector_{active_gene}",
    )
    inspected = junction_lookup.loc[inspected_junction]
    inspected_support = supporting_transcript_ids(
        inspected.get("supporting_transcripts", [])
    )
    inspected_names = [transcript_names.get(value, value) for value in inspected_support]
    inspector = st.container(border=True)
    inspector_columns = inspector.columns(3)
    inspector_columns[0].metric(
        "Junction", str(inspected.get("junction_label", "junction"))
    )
    inspector_columns[1].metric(
        "Class", str(inspected.get("junction_class", "unclassified")).replace("_", " ")
    )
    support_count = inspected.get("transcripts_with_junction", len(inspected_support))
    annotation_count = inspected.get("annotated_transcripts", len(annotated_transcript_ids))
    inspector_columns[2].metric(
        "Transcript support", f"{int(support_count)}/{int(annotation_count)}"
        if pd.notna(support_count) and pd.notna(annotation_count) else "not available",
    )
    inspector.caption("Supporting transcripts · " + (
        ", ".join(inspected_names) if inspected_names else "not available"
    ))
    inspector.caption(f"Stable genomic junction · {inspected_junction}")
else:
    st.info(
        "Choose one or more junctions above. A good starting point is a green "
        "specific junction supported by the transcript you want to distinguish."
    )

junction_figure = ps.pl.junction_map(
    active_gene, td, visible_summary, selected_junctions=selected_junctions,
    transcript_ids=visible_transcript_ids,
    title=f"{active_gene}: choose annotated isoform junctions",
)
st.pyplot(junction_figure, use_container_width=True)
plt.close(junction_figure)

with st.expander("View complete junction table", expanded=False):
    st.dataframe(
        visible_summary.drop(columns=["supporting_transcripts"], errors="ignore"),
        hide_index=True, use_container_width=True,
    )

design_button_label = (
    f"Design probes for {len(selected_junctions)} selected junction"
    f"{'s' if len(selected_junctions) != 1 else ''} →"
    if selected_junctions else "Select junctions to continue"
)
if selection_changed_after_design:
    st.info(
        "The target-junction selection changed, so the previous candidate table "
        "was cleared. Design the currently selected junctions to continue."
    )
if st.button(design_button_label, type="primary", disabled=not selected_junctions):
    with st.spinner("Calculating thermodynamics and applying chemistry rules…"):
        subset = design_pool.query("junction_id in @selected", local_dict={"selected": selected_junctions})
        if annotation_basis is not None:
            reannotation = visible_summary.set_index("junction_id")
            reannotation_columns = {
                "junction_class": "junction_specificity_class",
                "transcripts_with_junction": "junction_transcript_count",
                "annotated_transcripts": "gene_transcript_count",
                "transcript_coverage_percent": "target_coverage",
                "supporting_transcripts": "junction_supporting_transcripts",
            }
            for summary_column, candidate_column in reannotation_columns.items():
                if summary_column in reannotation:
                    subset.obs[candidate_column] = subset.obs["junction_id"].map(
                        reannotation[summary_column]
                    )
        subset.characterize().annotate_scores(chemistry, score_column="probeharbour_score")
        if existing_bed and Path(existing_bed).is_file():
            existing = load_existing_panel(existing_bed, fasta_path)
            subset.check_panel_compatibility(existing)
        if run_blast:
            if not Path(f"{blast_database}.nhr").is_file():
                st.error("BLAST database index (.nhr) was not found.")
                st.stop()
            backend = ps.LocalBlastBackend(
                blast_database, database_id=Path(blast_database).name,
                database_version="user-configured", executable=Path(sys.executable).with_name("blastn"),
            )
            ps.check_specificity(subset, backend)
        st.session_state.scored_pool = subset
        st.session_state[scored_junctions_key] = selected_junction_signature
        st.session_state.pop("confirmed_panel", None)

scored_pool = st.session_state.get("scored_pool")
if scored_pool is None:
    st.stop()

junction_label_map = visible_summary.set_index("junction_id")["junction_label"]
scored_pool.obs["display_junction_label"] = scored_pool.obs["junction_id"].map(
    junction_label_map
).fillna("unlabelled")

st.markdown('<div id="candidate"></div>', unsafe_allow_html=True)
st.subheader("3 · Review and confirm candidates")
st.caption(
    "ProbeHarbour proposes one preselected candidate per target junction. Review "
    "the score and QC evidence, then keep or clear each Include checkbox."
)
with st.expander("What should I check before adding probes?", expanded=False):
    st.markdown(
        """
        - Confirm that the **junction** and transcript event are the ones you intended.
        - Prefer stronger **ProbeHarbour scores**, but inspect the underlying QC flags
          rather than treating the score as experimental validation.
        - Check **panel compatibility** for conflicts with the existing 10x panel and
          with probes already in your cart.
        - Open the candidate record later in the cart for the complete sequence,
          thermodynamic and specificity evidence available for that probe.
        """
    )
review = candidate_review_table(scored_pool)
review.insert(
    review.columns.get_loc("junction_rank"),
    "target_junction",
    review["junction_id"].map(junction_label_map).fillna("unlabelled"),
)
edited = st.data_editor(
    review, hide_index=True, use_container_width=True, key=f"candidate_editor_{active_gene}",
    disabled=[column for column in review.columns if column != "include"],
    column_config={
        "include": st.column_config.CheckboxColumn("Include in panel"),
        "target_junction": st.column_config.TextColumn(
            "Target junction",
            help="The gene-wide J-label shown on the junction-selection plot.",
        ),
        "candidate_id": None,
        "junction_id": st.column_config.TextColumn("Stable genomic junction"),
    },
)

if st.button("🛒 Add selected probes to design cart", type="primary"):
    ids = edited.loc[edited["include"], "candidate_id"].astype(str).tolist()
    additions = selected_panel(scored_pool, ids, name=f"{active_gene}_confirmed")
    basket = add_to_basket(
        st.session_state.get("probe_basket"), additions, name=panel_name,
    )
    basket.check_internal_compatibility()
    st.session_state.confirmed_panel = additions
    st.session_state.probe_basket = basket
    st.session_state.pop("prepared_export", None)
    st.session_state.open_cart_next = True
    st.session_state.cart_notice = (
        f"Added {len(additions.obs)} probe(s); the cart now contains {len(basket.obs)} probes."
    )
    st.rerun()

basket = st.session_state.get("probe_basket")
if basket is None or basket.obs.empty:
    st.stop()

st.markdown('<div id="cart"></div>', unsafe_allow_html=True)
st.subheader("4 · Detailed cart inspection")
st.caption(
    "This is your working add-on panel. Inspect probes by gene, remove unwanted "
    "items, then return to Gene to add targets from another locus."
)
gene_values = pd.Series(pd.NA, index=basket.obs.index, dtype="object")
for gene_column in ("gene_name", "gene", "gene_id"):
    if gene_column in basket.obs:
        gene_values = gene_values.fillna(basket.obs[gene_column])
basket_genes = sorted(gene_values.dropna().astype(str).unique(), key=str.casefold)
basket_display = basket.obs.copy()
basket_display["display_gene"] = gene_values
basket_display["target_junction"] = (
    basket_display["display_junction_label"]
    if "display_junction_label" in basket_display else pd.Series(
        pd.NA, index=basket_display.index, dtype="object"
    )
)
active_gene_rows = gene_values.astype(str).eq(str(active_gene))
basket_display.loc[active_gene_rows, "target_junction"] = (
    basket_display.loc[active_gene_rows, "junction_id"].map(junction_label_map)
    .fillna(basket_display.loc[active_gene_rows, "target_junction"])
)

def cart_candidate_label(candidate_id):
    row = basket_display.loc[
        basket_display["candidate_id"].astype(str).eq(str(candidate_id))
    ].iloc[0]
    gene_label = row.get("display_gene")
    target_label = row.get("target_junction")
    shift_value = pd.to_numeric(row.get("shift_nt"), errors="coerce")
    shift_label = (
        f" · slide {int(shift_value):+d} nt" if pd.notna(shift_value) else ""
    )
    return f"{gene_label} · {target_label}{shift_label}"

metrics = st.columns(4)
metrics[0].metric("Basket probes", len(basket.obs))
metrics[1].metric("Genes", len(basket_genes))
metrics[2].metric("Target junctions", basket.obs["junction_id"].nunique())
existing_compatibility = existing_panel_compatibility(basket.obs)
existing_compatibility_rate = (
    f"{existing_compatibility.dropna().astype(bool).mean():.0%}"
    if existing_compatibility.notna().any() else "Not evaluated"
)
metrics[3].metric("Existing 10x compatible", existing_compatibility_rate)

plot_tab, cards_tab, table_tab = st.tabs(["Gene map", "Candidate record", "Complete table"])
with plot_tab:
    map_index = basket_genes.index(active_gene) if active_gene in basket_genes else 0
    map_gene = st.selectbox("Basket gene to display", basket_genes, index=map_index)
    display_store = basket
    if existing_bed and Path(existing_bed).is_file():
        display_store = ps.ProbeStore.merge(
            [load_existing_panel(existing_bed, fasta_path), basket],
            name=f"existing_plus_{basket.name}", roles=["existing", "add_on"],
        )
    figure = ps.pl.gene_map(
        map_gene, td, display_store,
        title=f"{map_gene}: existing 10x and confirmed add-on probes",
    )
    st.pyplot(figure, use_container_width=True)
    st.caption(
        "Existing 10x probes retain their BED geometry: purple junction arcs "
        "span exon blocks, while purple diamonds mark exon-contained probes."
    )
    plt.close(figure)
with cards_tab:
    card_id = st.selectbox(
        "Candidate",
        basket.obs["candidate_id"].astype(str).tolist(),
        format_func=cart_candidate_label,
    )
    card_row = basket_display.loc[
        basket_display["candidate_id"].astype(str).eq(str(card_id))
    ].iloc[0].copy()
    card_row["probe_name"] = cart_candidate_label(card_id)
    figure = ps.pl.candidate_card(
        basket, card_row, title="Confirmed ProbeHarbour candidate"
    )
    st.pyplot(figure, use_container_width=True)
    plt.close(figure)
with table_tab:
    st.caption(
        "Target junction uses the same gene-wide J-label as the selection plot. "
        "Stable genomic coordinates remain available for unambiguous export."
    )
    preferred_cart_columns = [
        "display_gene", "target_junction", "junction_id", "shift_nt", "rbd5", "rbd3",
        "lig_dinuc", "tm_full", "probeharbour_score", "specificity_score",
        "existing_panel_compatibility_pass", "selection_status",
    ]
    st.dataframe(
        basket_display[[
            column for column in preferred_cart_columns if column in basket_display
        ]],
        hide_index=True,
        use_container_width=True,
        column_config={
            "display_gene": "Gene",
            "target_junction": "Target junction",
            "junction_id": "Stable genomic junction",
            "shift_nt": "Slide (nt)",
            "probeharbour_score": st.column_config.NumberColumn(
                "ProbeHarbour score", format="%.4f"
            ),
        },
    )
    with st.expander("Advanced audit fields", expanded=False):
        st.caption(
            "Internal transcript-local candidate IDs and full provenance are "
            "retained here for reproducibility; they are not gene-wide J-labels."
        )
        st.dataframe(basket.obs, hide_index=True, use_container_width=True)

st.markdown('<div id="export"></div>', unsafe_allow_html=True)
st.subheader("5 · Finalize ProbeSet")
st.caption(
    "Export only after reviewing the complete multi-gene cart. The ZIP contains the "
    "panel records and assay export files; changing the cart invalidates the prepared package."
)
export_flag_count = int(existing_compatibility.eq(False).sum())
if export_flag_count:
    st.warning(
        f"The cart contains {export_flag_count} probe(s) that conflict with the "
        "configured existing 10x panel."
    )
if st.button("Prepare export package"):
    st.session_state.prepared_export = {
        "name": panel_name,
        "payload": export_panel_bundle(basket, panel_name=panel_name),
        "probe_count": len(basket.obs),
    }
prepared_export = st.session_state.get("prepared_export")
if prepared_export is not None:
    st.success(
        f"Export prepared for {prepared_export['probe_count']} probes. "
        "It will be invalidated if the cart changes."
    )
    st.download_button(
        "Download finalized ProbeSet (.zip)", prepared_export["payload"],
        file_name=f"{prepared_export['name']}.zip", mime="application/zip",
    )
