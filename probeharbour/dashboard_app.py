"""ProbeHarbour Streamlit dashboard: gene to reviewed, exportable panel."""
from __future__ import annotations

import gzip
import importlib
import inspect
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
    remove_from_basket,
    selected_panel,
)


st.set_page_config(page_title="ProbeHarbour", page_icon="⚓", layout="wide")
st.title("ProbeHarbour")
st.caption("Design, evaluate and export isoform-aware Xenium probe panels.")


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


with st.sidebar:
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
        compatibility = (sidebar_basket.obs["panel_compatibility_pass"]
                         if "panel_compatibility_pass" in sidebar_basket.obs
                         else pd.Series(pd.NA, index=sidebar_basket.obs.index))
        flagged_count = int(compatibility.eq(False).sum())
        basket_metrics = st.columns(4)
        basket_metrics[0].metric("Probes", len(sidebar_basket.obs))
        basket_metrics[1].metric("Genes", basket_gene_count)
        basket_metrics[2].metric(
            "Junctions",
            sidebar_basket.obs["junction_id"].nunique()
            if "junction_id" in sidebar_basket.obs else 0,
        )
        basket_metrics[3].metric("Compatibility flags", flagged_count)
        basket_preview = sidebar_basket.obs.copy()
        basket_preview["gene"] = basket_gene_values
        if flagged_count:
            st.warning(
                f"{flagged_count} probe(s) have a panel compatibility flag. "
                "Review or replace them before final export."
            )
        else:
            st.success("No panel compatibility conflicts are currently detected.")

        gene_summary = basket_preview.groupby("gene", dropna=False).agg(
            probes=("candidate_id", "size"),
            junctions=("junction_id", "nunique"),
            median_score=("probeharbour_score", "median"),
        )
        if "panel_compatibility_pass" in basket_preview:
            gene_summary["compatibility_flags"] = basket_preview.groupby(
                "gene", dropna=False
            )["panel_compatibility_pass"].apply(lambda values: int(values.eq(False).sum()))
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
                plt.close(cart_figure)
            else:
                st.warning("A valid GTF and FASTA are required for the transcript map.")

        cart_columns = [
            "gene", "candidate_id", "junction_id", "shift_nt", "probeharbour_score",
            "panel_compatibility_pass",
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
                "panel_compatibility_pass": st.column_config.CheckboxColumn("Compatible"),
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

st.subheader("1 · Choose a gene")
demo_gene = "Snap25" if "Snap25" in genes else None
demo_index = genes.index(demo_gene) if demo_gene is not None else None
gene = st.selectbox(
    "Gene", genes, index=demo_index, placeholder="Type to search the annotation…",
    help="Snap25 is loaded initially as an isoform-junction design demonstration.",
)
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
st.subheader("2 · Select biological junctions")
st.caption("Use the transcript models and J-labels to identify the exact annotated junction to target.")

annotated_transcript_ids = list(td.get_transcripts_by_gene_name(active_gene))
transcript_names = {
    str(transcript_id): str(
        td.get_transcript_info(transcript_id).get("transcript_name") or transcript_id
    )
    for transcript_id in annotated_transcript_ids
}
visible_transcript_ids = st.multiselect(
    "Transcripts shown",
    [str(value) for value in annotated_transcript_ids],
    default=[str(value) for value in annotated_transcript_ids],
    format_func=lambda value: transcript_names.get(str(value), str(value)),
    key=f"visible_transcripts_{active_gene}",
    help="Hide transcript models that are not relevant to the isoform event you are designing.",
)
renumber_visible = st.checkbox(
    "Renumber the visible junctions from J1",
    value=False,
    key=f"renumber_junctions_{active_gene}",
    help="Changes only the short display labels. Stable genomic junction IDs are retained for design and export.",
)
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

if "supporting_transcripts" in summary:
    visible_summary = summary.loc[
        summary["supporting_transcripts"].map(is_visible_junction)
    ].copy().reset_index(drop=True)
else:
    visible_summary = summary.copy()
if renumber_visible:
    visible_summary["junction_label"] = [f"J{index + 1}" for index in range(len(visible_summary))]

all_junction_options = visible_summary["junction_id"].astype(str).tolist()
junction_key = f"junction_selection_{active_gene}"
if junction_key in st.session_state:
    st.session_state[junction_key] = [
        value for value in st.session_state[junction_key] if value in all_junction_options
    ]
junction_lookup = visible_summary.set_index("junction_id")
junction_search = st.text_input(
    "Find a junction",
    key=f"junction_search_{active_gene}",
    placeholder="Search J4, specific, genomic coordinates, Pkm-202, ENSMUST…",
    help=("Searches the visible J-label, junction class, stable genomic ID, "
          "and every supporting transcript ID or name."),
).strip().lower()
search_tokens = junction_search.split()

def junction_matches_search(junction_id):
    row = junction_lookup.loc[junction_id]
    label = str(row.get("junction_label", "")).lower()
    supporting = supporting_transcript_ids(row.get("supporting_transcripts", []))
    supporting_names = [transcript_names.get(value, value) for value in supporting]
    searchable = " ".join([
        label, str(junction_id), str(row.get("junction_class", "")),
        str(row.get("transcripts_with_junction", "")),
        str(row.get("annotated_transcripts", "")),
        *supporting, *supporting_names,
    ]).lower()
    for token in search_tokens:
        if re.fullmatch(r"j\d+", token):
            if token != label:
                return False
        elif token not in searchable:
            return False
    return True

search_matches = [value for value in all_junction_options if junction_matches_search(value)]
retained_selections = [
    value for value in st.session_state.get(junction_key, [])
    if value in all_junction_options
]
junction_options = [
    value for value in all_junction_options
    if value in search_matches or value in retained_selections
]
if junction_search:
    st.caption(
        f"{len(search_matches)} matching junction(s) · "
        f"{len(retained_selections)} current selection(s) remain available"
    )

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
    "Junctions to target", junction_options, **junction_widget_options,
    format_func=junction_option_label,
    placeholder=("Choose from matching junctions" if junction_search
                 else "Choose annotated junctions"),
)
junction_figure = ps.pl.junction_map(
    active_gene, td, visible_summary, selected_junctions=selected_junctions,
    transcript_ids=visible_transcript_ids,
    title=f"{active_gene}: choose annotated isoform junctions",
)
st.pyplot(junction_figure, use_container_width=True)
plt.close(junction_figure)

st.dataframe(
    visible_summary.drop(columns=["supporting_transcripts"], errors="ignore"),
    hide_index=True, use_container_width=True,
)

if st.button("Design and score candidates", type="primary", disabled=not selected_junctions):
    with st.spinner("Calculating thermodynamics and applying chemistry rules…"):
        subset = design_pool.query("junction_id in @selected", local_dict={"selected": selected_junctions})
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
        st.session_state.pop("confirmed_panel", None)

scored_pool = st.session_state.get("scored_pool")
if scored_pool is None:
    st.stop()

st.subheader("3 · Review and confirm candidates")
st.caption("The highest-scoring candidate per junction is preselected. Change the checkboxes before confirmation.")
review = candidate_review_table(scored_pool)
edited = st.data_editor(
    review, hide_index=True, use_container_width=True, key=f"candidate_editor_{active_gene}",
    disabled=[column for column in review.columns if column != "include"],
    column_config={"include": st.column_config.CheckboxColumn("Include in panel")},
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

st.subheader("4 · Detailed cart inspection")
gene_values = pd.Series(pd.NA, index=basket.obs.index, dtype="object")
for gene_column in ("gene_name", "gene", "gene_id"):
    if gene_column in basket.obs:
        gene_values = gene_values.fillna(basket.obs[gene_column])
basket_genes = sorted(gene_values.dropna().astype(str).unique(), key=str.casefold)
metrics = st.columns(4)
metrics[0].metric("Basket probes", len(basket.obs))
metrics[1].metric("Genes", len(basket_genes))
metrics[2].metric("Target junctions", basket.obs["junction_id"].nunique())
metrics[3].metric("Panel-compatible", f"{basket.obs['panel_compatibility_pass'].mean():.0%}")

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
    plt.close(figure)
with cards_tab:
    card_id = st.selectbox("Candidate", basket.obs["candidate_id"].astype(str).tolist())
    figure = ps.pl.candidate_card(basket, card_id, title="Confirmed ProbeHarbour candidate")
    st.pyplot(figure, use_container_width=True)
    plt.close(figure)
with table_tab:
    st.dataframe(basket.obs, hide_index=True, use_container_width=True)

st.subheader("5 · Finalize ProbeSet")
st.caption(
    "Preparing an export is an explicit checkout step. Adding or removing a probe invalidates "
    "the previously prepared package."
)
export_flag_count = int(basket.obs["panel_compatibility_pass"].eq(False).sum())
if export_flag_count:
    st.warning(f"The cart currently contains {export_flag_count} compatibility-flagged probe(s).")
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
