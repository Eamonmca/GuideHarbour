"""Sequence-aware plots and an offline report for guide detection probes.

This module is maintained as Python, independently of the upstream nbdev
notebooks. All coordinates refer to the sense target (zero-based, half-open).
"""
from __future__ import annotations

from html import escape
from pathlib import Path
from typing import Any

import pandas as pd

__all__ = ["guide_overview", "plot_guide_candidates", "write_guide_report"]

GUIDE = "#087f8c"
SCAFFOLD = "#dce4eb"
LEADER = "#c4b5fd"
UPSTREAM = "#28558c"
DOWNSTREAM = "#d97732"
INK = "#17324d"


def _records(store) -> list[dict]:
    """Read normalized target records, including targets with zero candidates."""
    records = store.metadata.get("guides")
    if records is None:
        raise ValueError("Expected a guide design ProbeStore with metadata['guides']")
    return list(records)


def _top(store, n: int, eligible_only: bool) -> pd.DataFrame:
    if isinstance(n, bool) or not isinstance(n, int) or n < 1:
        raise ValueError("n must be a positive integer")
    from .guides import top_guides
    return top_guides(store, n=n, eligible_only=eligible_only)


def _value(value: Any, digits: int | None = None) -> str:
    if value is None or (not isinstance(value, (list, dict, tuple)) and pd.isna(value)):
        return "not assessed"
    if digits is not None:
        return f"{float(value):.{digits}f}"
    if isinstance(value, (list, tuple)):
        return "; ".join(str(part) for part in value) or "none"
    return str(value)


def _target_segments(record: dict) -> list[tuple[int, int, str]]:
    start, end = int(record["guide_start"]), int(record["guide_end"])
    length = len(record["target_sequence"])
    return [(0, start, LEADER), (start, end, GUIDE), (end, length, SCAFFOLD)]


def _plot_target(ax, record: dict, y: float, *, height: float = 0.20) -> None:
    for start, end, color in _target_segments(record):
        if end > start:
            ax.barh(y, end - start, left=start, height=height, color=color, edgecolor="none")


def guide_overview(store, *, n: int = 3, eligible_only: bool = True):
    """Return a figure with the best target window and top scores per guide.

    A missing shortlist stays visibly missing. The score is inherited chemistry
    evidence; a panel mismatch filter does not establish transcriptome specificity.
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    records = _records(store)
    top = _top(store, n, eligible_only)
    if not records:
        fig, ax = plt.subplots(figsize=(10, 3))
        ax.text(.5, .5, "No guide targets supplied", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
        return fig
    fig, (ax, score_ax) = plt.subplots(
        1, 2, figsize=(13.8, max(4.5, 1.7 + .43 * len(records))),
        gridspec_kw={"width_ratios": [3.6, 1.2], "wspace": .13}, sharey=True,
    )
    # Keep the complete target visible and use the best candidate only on the
    # sequence track; score dots show the remaining shortlist candidates.
    max_length = max(len(record["target_sequence"]) for record in records)
    for index, record in enumerate(records):
        _plot_target(ax, record, index, height=.23)
        rows = top.loc[top["guide_id"].eq(record["guide_id"])] if not top.empty else top
        if rows.empty:
            ax.text(max_length + 2, index, "no eligible candidate" if eligible_only else "no candidate",
                    va="center", fontsize=8, color="#a54538")
            continue
        row = rows.iloc[0]
        start, cut, end = int(row["start"]), int(row["ligation_site"]), int(row["end"])
        ax.barh(index + .24, cut - start, left=start, height=.12, color=UPSTREAM)
        ax.barh(index + .24, end - cut, left=cut, height=.12, color=DOWNSTREAM)
        ax.plot([cut, cut], [index + .12, index + .36], color=INK, lw=1)
        for rank, (_, candidate) in enumerate(rows.iterrows()):
            score = float(candidate["recalculated_score"])
            score_ax.scatter(score, index + .13 * rank - .06, s=46 if rank == 0 else 17,
                             color=GUIDE if rank == 0 else "#97b5be", zorder=3,
                             marker="o" if bool(candidate["eligible"]) else "x")
        score_ax.text(1.06, index, f"{float(row['recalculated_score']):.3f}", va="center", fontsize=8, color=INK)
    ax.set_yticks(range(len(records)), [record["guide_id"] for record in records], fontsize=9)
    ax.set_ylim(len(records) - .45, -.75)
    ax.set_xlim(-1, max_length + 24)
    ax.set_xlabel("Sense target position (0-based; DNA alphabet)", fontsize=9)
    ax.set_title("Best candidate on each target", loc="left", color=INK, fontsize=12, pad=14)
    score_ax.set_xlim(-.04, 1.26)
    score_ax.set_xticks([0, .25, .5, .75, 1.0])
    score_ax.set_xlabel("Chemistry score", fontsize=9)
    score_ax.set_title(f"Top {n} eligible" if eligible_only else f"Top {n} (includes failures)",
                       loc="left", color=INK, fontsize=12, pad=14)
    score_ax.tick_params(axis="y", left=False, labelleft=False)
    score_ax.grid(axis="x", color="#edf0f4", zorder=0)
    for axis in (ax, score_ax):
        axis.spines[["top", "right", "left"]].set_visible(False)
        axis.spines["bottom"].set_color("#cbd5df")
        axis.tick_params(colors=INK, length=0)
    fig.suptitle("GuideHarbour · guide detection candidates", x=.06, ha="left", color=INK, fontsize=17, fontweight="bold")
    handles = [Patch(color=GUIDE, label="Guide"), Patch(color=SCAFFOLD, label="Scaffold"),
               Patch(color=LEADER, label="Explicit leader / construct prefix"),
               Patch(color=UPSTREAM, label="Upstream target arm"),
               Patch(color=DOWNSTREAM, label="Downstream target arm")]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=8, bbox_to_anchor=(.5, .005))
    fig.subplots_adjust(left=.18, right=.96, top=.88, bottom=.17)
    return fig


def plot_guide_candidates(store, guide_id: str, *, n: int = 3, eligible_only: bool = True):
    """Return a target-coordinate figure of the top candidates for one guide."""
    import matplotlib.pyplot as plt

    records = {record["guide_id"]: record for record in _records(store)}
    if guide_id not in records:
        raise KeyError(f"Unknown guide_id: {guide_id!r}")
    record = records[guide_id]
    top = _top(store, n, eligible_only)
    rows = top.loc[top["guide_id"].eq(guide_id)] if not top.empty else top
    fig, ax = plt.subplots(figsize=(13, max(3.4, 2.4 + .65 * len(rows))))
    _plot_target(ax, record, 0, height=.28)
    labels = ["Target 5′ → 3′"]
    for index, (_, row) in enumerate(rows.iterrows(), start=1):
        start, cut, end = int(row["start"]), int(row["ligation_site"]), int(row["end"])
        ax.barh(index, cut - start, left=start, height=.27, color=UPSTREAM)
        ax.barh(index, end - cut, left=cut, height=.27, color=DOWNSTREAM)
        ax.plot([cut, cut], [index - .2, index + .2], color=INK, lw=1)
        ax.text(end + 1, index, f"score {row['recalculated_score']:.4f} · {int(row['guide_overlap_nt'])} guide nt · "
                f"{_value(row['nearest_other_mismatches'])} nearest mismatches", va="center", fontsize=9)
        labels.append(f"#{index} · shift {int(row['shift_nt']):+d}")
    if rows.empty:
        ax.text(0, 1, "No eligible candidate for these settings." if eligible_only else "No candidates generated.", color="#a54538")
    ax.set_yticks(range(len(labels)), labels)
    ax.set_ylim(max(1.5, len(rows) + .7), -.7)
    ax.set_xlim(-1, len(record["target_sequence"]) + 30)
    ax.set_xlabel("Sense target position (0-based, half-open)")
    ax.set_title(f"{guide_id} · top {n} guide detection candidates", loc="left", fontsize=15, color=INK, pad=15)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    fig.tight_layout()
    return fig


def _sequence_html(record: dict, row: dict | None = None) -> str:
    parts = []
    start, end = int(record["guide_start"]), int(record["guide_end"])
    for index, base in enumerate(record["target_sequence"]):
        kind = "guide" if start <= index < end else "leader" if index < start else "scaffold"
        classes = [kind]
        if row is not None and int(row["start"]) <= index < int(row["end"]):
            classes.append("arm-up" if index < int(row["ligation_site"]) else "arm-down")
        if row is not None and index == int(row["ligation_site"]):
            classes.append("cut")
        parts.append(f'<span class="{" ".join(classes)}" title="Position {index}">{escape(base)}</span>')
    return '<div class="sequence" aria-label="Sense target sequence">' + "".join(parts) + "</div>"


def _metrics_html(row: dict) -> str:
    metrics = [
        ("Chemistry score", _value(row["recalculated_score"], 4)),
        ("Guide overlap", f"{int(row['guide_overlap_nt'])} nt"),
        ("Full-window Tm proxy", f"{_value(row.get('tm_full'), 1)} °C"),
        ("Nearest other guide", _value(row.get("nearest_other_guide_id"))),
        ("Panel mismatches", _value(row.get("nearest_other_mismatches"))),
        ("Ligation dinucleotide", row.get("lig_dinuc", "")),
    ]
    return '<dl class="metrics">' + "".join(
        f"<div><dt>{escape(label)}</dt><dd>{escape(str(value))}</dd></div>" for label, value in metrics
    ) + "</dl>"


def _candidate_html(record: dict, row: dict, index: int) -> str:
    status = "Eligible" if bool(row["eligible"]) else "Excluded"
    badge_class = "badge" if bool(row["eligible"]) else "badge excluded"
    details = [
        ("Upstream target arm, 5′→3′ (rbd5)", row["rbd5"]),
        ("Downstream target arm, 5′→3′ (rbd3)", row["rbd3"]),
        ("Oligo piece binding upstream arm, 5′→3′", row["probe_arm_upstream_5to3"]),
        ("Oligo piece binding downstream arm, 5′→3′", row["probe_arm_downstream_5to3"]),
        ("Reverse complement of whole target window, 5′→3′", row["reverse_complement_binder"]),
    ]
    seq_details = "".join(f"<dt>{escape(label)}</dt><dd><code>{escape(str(value))}</code></dd>" for label, value in details)
    warnings = _value(row.get("warnings", []))
    return f'''<article class="candidate">
      <div class="candidate-heading"><h3>Candidate {index} <span class="{badge_class}">{status}</span></h3>
      <span>shift {int(row['shift_nt']):+d} nt · target [{int(row['start'])}, {int(row['end'])}) · cut {int(row['ligation_site'])}</span></div>
      <div class="candidate-id">{escape(str(row['candidate_id']))}</div>
      {_sequence_html(record, row)}{_metrics_html(row)}
      <p class="small">Active chemistry modifiers: {escape(_value(row.get('modifiers', [])))}.</p>
      <details><summary>Arm sequences and evidence</summary><dl class="seq-details">{seq_details}</dl>
      <p class="small">The reverse-complement pieces show hybridizing sequence orientation. These are not complete order-ready vendor probes: backbone, barcodes, linkers and modifications are not added.</p>
      <p class="small">Shared arm matches in other guide targets: upstream {escape(_value(row.get('upstream_arm_shared_count')))}; downstream {escape(_value(row.get('downstream_arm_shared_count')))}.</p>
      <p class="small">Warnings: {escape(warnings)}.</p>
      <p class="small">Exclusion reasons: {escape(_value(row.get('exclusion_reasons', [])))}.</p></details>
      </article>'''


def write_guide_report(store, path: str | Path, *, n: int = 3, eligible_only: bool = True) -> Path:
    """Write a self-contained offline HTML report; no data leaves the browser.

    Guide IDs and sequences are escaped, never interpolated into JavaScript.
    Per-guide failures remain visible even when no candidate can be shortlisted.
    """
    from .guides import guide_summary

    records = _records(store)
    top = _top(store, n, eligible_only)
    summary = guide_summary(store)
    by_id = {row["guide_id"]: row for row in summary.to_dict("records")}
    design = store.metadata.get("guide_design", {})
    config = design.get("config", {})
    qc_policy = ("Chemistry modifiers are hard filters in this run; eligible candidates have none."
                 if config.get("require_qc_pass", False) else
                 "Eligible means the guide-panel and minimum-score filters pass; chemistry modifiers may still apply as score penalties.")
    sections, options = [], []
    for index, record in enumerate(records):
        guide_id = record["guide_id"]
        status = by_id.get(guide_id, {})
        rows = top.loc[top["guide_id"].eq(guide_id)] if not top.empty else top
        options.append(f'<option value="guide-{index}">{escape(guide_id)}</option>')
        candidates = "".join(_candidate_html(record, row, rank) for rank, row in enumerate(rows.to_dict("records"), start=1))
        if not candidates:
            candidates = '<div class="empty">No eligible candidate under these settings.</div>' if eligible_only else '<div class="empty">No candidates generated under these settings.</div>'
        all_rows = store.obs.loc[store.obs["guide_id"].eq(guide_id)] if not store.obs.empty else store.obs
        excluded = all_rows.loc[~all_rows["eligible"].astype(bool)] if not all_rows.empty else all_rows
        exclusions = "".join(
            f"<tr><td>{escape(str(row['candidate_id']))}</td><td>{escape(_value(row.get('exclusion_reasons', [])))}</td></tr>"
            for row in excluded.to_dict("records")
        )
        extra = f'<details class="exclusions"><summary>{len(excluded)} excluded candidate(s)</summary><table><thead><tr><th>Candidate</th><th>Reason</th></tr></thead><tbody>{exclusions}</tbody></table></details>' if len(excluded) else ""
        generated = status.get("n_candidates", len(all_rows))
        eligible = status.get("n_eligible", int(all_rows["eligible"].sum()) if not all_rows.empty else 0)
        skipped = "; ".join(f"{label}: {int(status.get(field, 0))}" for field, label in [
            ("skipped_out_of_bounds", "out of bounds"), ("skipped_not_spanning", "not spanning boundary"),
            ("skipped_guide_overlap", "insufficient guide overlap")])
        sections.append(f'''<section id="guide-{index}" class="guide-panel"{' hidden' if index else ''}>
            <h2>{escape(guide_id)}</h2><p class="meta">{len(record['guide_sequence'])} nt guide · {len(record['target_sequence'])} nt target · {generated} scored · {eligible} eligible</p>
            <p class="small">Generation skips — {escape(skipped)}.</p>
            <div class="target-card"><strong>Sense target, 5′ → 3′</strong>{_sequence_html(record)}
            <p class="small">Guide shown in teal; scaffold in pale blue; explicit leader or construct prefix in purple. Sequences use the DNA alphabet (T represents RNA U).</p></div>
            {candidates}{extra}</section>''')
    n_eligible = int(store.obs["eligible"].sum()) if not store.obs.empty else 0
    shortlist_guides = top["guide_id"].nunique() if not top.empty else 0
    no_targets = '<div class="empty">No guide targets supplied.</div>' if not records else ""
    page = f'''<!doctype html><html lang="en"><head><meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1"><title>GuideHarbour · probe report</title>
    <style>
    :root{{--ink:#17324d;--muted:#5c6d7d;--teal:#087f8c;--line:#dbe4ec;--bg:#f4f7fa}}
    *{{box-sizing:border-box}}body{{margin:0;background:var(--bg);color:var(--ink);font:15px/1.55 Arial,sans-serif}}
    header{{background:#17324d;color:white;padding:40px max(24px,calc((100vw - 1170px)/2)) 30px}}
    header .eyebrow{{font-size:12px;text-transform:uppercase;letter-spacing:.15em;color:#9ee0db}}
    h1{{margin:8px 0;font-size:34px;line-height:1.2}}header p{{max-width:850px;color:#c6d7e7;margin:8px 0}}
    main{{max-width:1220px;margin:auto;padding:24px}}h2{{font-size:25px;margin:25px 0 4px}}h3{{font-size:17px;margin:0}}
    .stats{{display:flex;gap:14px;flex-wrap:wrap;margin:0 0 22px}}.stat{{background:white;border:1px solid var(--line);border-radius:12px;padding:16px 22px;min-width:160px}}
    .stat strong{{display:block;font-size:26px;line-height:1.3}}.stat span,.meta,.small{{color:var(--muted)}}.small{{font-size:12px;margin:8px 0}}
    .notice{{padding:16px 20px;border-left:4px solid #e1ab47;background:#fff8e9;border-radius:6px;font-size:13px;margin:16px 0 22px}}
    .controls{{display:flex;align-items:center;gap:15px;flex-wrap:wrap;background:white;padding:15px 20px;border:1px solid var(--line);border-radius:10px;position:sticky;top:0;z-index:3}}
    select{{font:inherit;padding:8px 12px;border:1px solid #b8c8d5;border-radius:6px;background:white;min-width:240px;color:var(--ink)}}
    .legend{{font-size:12px;margin-left:auto;color:var(--muted)}}.swatch{{display:inline-block;width:10px;height:10px;margin:0 5px 0 12px;border-radius:2px}}
    .candidate,.target-card{{background:white;border:1px solid var(--line);border-radius:12px;padding:20px;margin:16px 0}}
    .candidate-heading{{display:flex;gap:12px;align-items:center;justify-content:space-between;flex-wrap:wrap}}.candidate-heading>span,.candidate-id{{color:var(--muted);font-size:12px}}.candidate-id{{word-break:break-all;margin:4px 0}}
    .badge{{display:inline-block;background:#e5f4f0;color:#116254;font-size:11px;font-weight:normal;padding:2px 8px;border-radius:20px;vertical-align:middle;margin-left:8px}}.badge.excluded{{background:#fff0e8;color:#9a452a}}
    .sequence{{white-space:nowrap;overflow-x:auto;max-width:100%;padding:13px 0 9px;font:13px/2.2 Consolas,monospace;letter-spacing:0}}
    .sequence>span{{display:inline-block;min-width:10px;text-align:center;border-bottom:3px solid transparent;position:relative}}
    .guide{{background:#daf0ec;color:#00656d}}.scaffold{{background:#edf2f6;color:#66798b}}.leader{{background:#f0eaff;color:#7255a5}}
    .sequence .arm-up{{border-bottom-color:#28558c;color:#17324d;font-weight:bold}}.sequence .arm-down{{border-bottom-color:#d97732;color:#9b4e12;font-weight:bold}}
    .sequence .cut::before{{content:"";position:absolute;left:-1px;top:0;bottom:-5px;border-left:2px solid #17324d}}
    .metrics{{display:grid;grid-template-columns:repeat(6,minmax(0,1fr));gap:14px;margin:15px 0}}dt{{font-size:11px;color:var(--muted)}}dd{{margin:3px 0 0;font-size:16px;font-weight:bold;overflow-wrap:anywhere}}
    details{{margin-top:12px}}summary{{cursor:pointer;color:var(--teal);font-size:13px;font-weight:bold}}.seq-details{{display:grid;grid-template-columns:minmax(220px,1fr) 2fr;gap:8px 16px}}
    .seq-details dt{{font-size:12px}}.seq-details dd{{margin:0;font-size:12px;font-weight:normal;overflow-wrap:anywhere}}code{{font-family:Consolas,monospace}}
    .empty{{padding:25px;background:#fff2ee;border:1px solid #f1d6ca;border-radius:10px;color:#914637;margin-top:16px}}
    table{{border-collapse:collapse;width:100%;font-size:12px;margin:14px 0}}td,th{{border-bottom:1px solid var(--line);padding:8px;text-align:left;word-break:break-word}}
    footer{{font-size:12px;color:var(--muted);padding:25px 0}}[hidden]{{display:none!important}}
    @media(max-width:800px){{.metrics{{grid-template-columns:repeat(2,minmax(0,1fr))}}.seq-details{{display:block}}.seq-details dd{{margin-bottom:10px}}.legend{{margin-left:0}}}}
    @media print{{.controls{{position:static}}.guide-panel[hidden]{{display:block!important}}.candidate{{break-inside:avoid}}body{{background:white}}header{{padding:20px}}}}
    </style></head><body><header><div class="eyebrow">GuideHarbour / sequence-first design</div>
    <h1>Guide detection probe shortlist</h1><p>Inspect the guide-specific sequence in each candidate, its chemistry score and its comparison against the supplied guide panel.</p></header>
    <main><div class="stats"><div class="stat"><strong>{len(records)}</strong><span>guide targets</span></div>
    <div class="stat"><strong>{len(store.obs)}</strong><span>scored candidates</span></div>
    <div class="stat"><strong>{n_eligible}</strong><span>eligible candidates</span></div>
    <div class="stat"><strong>{shortlist_guides} / {len(records)}</strong><span>guides in shortlist</span></div></div>
    <div class="notice"><strong>Research shortlist.</strong> Scores reuse ProbeHarbour's chemistry rules and Primer3 DNA thermodynamic proxies. The specificity screen compares full target windows against the supplied guide targets only; it does not search the transcriptome or model RNA accessibility, Cas9 binding, ligation yield or experimental detection. A shared scaffold arm is reported. These sequences are not complete order-ready probes. Top candidates for one guide are overlapping alternatives, not a validated simultaneous probe set. {escape(qc_policy)}</div>
    <p class="small">Chemistry: {escape(str(config.get('chemistry', design.get('ruleset_id', 'unspecified'))))} · Target mode: {escape(str(config.get('target_mode', 'unspecified')))} · Minimum guide overlap: {escape(str(config.get('min_guide_overlap', 'unspecified')))} nt · Panel mismatch threshold: {escape(str(config.get('min_panel_mismatches', 'unspecified')))}.</p>
    <div class="controls"><label for="guide-picker"><strong>Guide</strong></label><select id="guide-picker"{' disabled' if not records else ''}>{''.join(options)}</select>
    <span class="legend"><span class="swatch" style="background:#28558c"></span>upstream target arm<span class="swatch" style="background:#d97732"></span>downstream target arm</span></div>
    {no_targets}{''.join(sections)}
    <footer>All coordinates are zero-based and half-open on the sense target. Hover over a base to inspect its position. The supplied prefix is excluded in RNA mode; only an explicit transcribed leader is included. Full-construct mode includes the prefix and requires the experimental target to be confirmed. This file works offline and loads no external scripts or assets.</footer></main>
    <script>document.getElementById('guide-picker').addEventListener('change', function(){{document.querySelectorAll('.guide-panel').forEach(function(panel){{panel.hidden = panel.id !== document.getElementById('guide-picker').value;}});}});</script>
    </body></html>'''
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(page, encoding="utf-8")
    return destination
