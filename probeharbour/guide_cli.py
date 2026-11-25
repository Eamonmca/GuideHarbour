"""Run guide probe design and write a reproducible, offline result directory."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from importlib import metadata, resources
import json
from pathlib import Path
import re

import pandas as pd


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Design guide + adjacent-sequence detection probes with ProbeHarbour.")
    parser.add_argument("guides", nargs="?", type=Path, help="CSV with guide_id, guide_sequence, scaffold_sequence; optional prefix_sequence and transcribed_leader")
    parser.add_argument("--example", action="store_true", help="Use the packaged 15-guide example")
    parser.add_argument("--output", type=Path, default=Path("guide_results"), help="Result directory (default: guide_results)")
    parser.add_argument("--chemistry", choices=["prime_5k", "xenium_v1"], default="prime_5k")
    parser.add_argument("--arm-length", type=int, default=None, help="Equal arm length; default 16 for Prime 5K, 20 for Xenium V1")
    parser.add_argument("--slide-min", type=int, default=-5)
    parser.add_argument("--slide-max", type=int, default=5)
    parser.add_argument("--boundary", choices=["guide_scaffold", "prefix_guide", "both"], default="guide_scaffold")
    parser.add_argument("--target-mode", choices=["rna", "full_construct"], default="rna", help="RNA excludes supplied promoter/cloning prefix; full_construct includes it")
    parser.add_argument("--min-guide-overlap", type=int, default=15)
    parser.add_argument("--min-panel-mismatches", type=int, default=3)
    parser.add_argument("--min-score", type=float, default=0.0)
    parser.add_argument("--require-qc-pass", action="store_true", help="Exclude candidates with any chemistry score modifier; default treats modifiers as penalties")
    parser.add_argument("--top", type=int, default=3, help="Eligible candidates per guide in shortlist and report")
    parser.add_argument("--overwrite", action="store_true", help="Replace known output filenames in an existing directory")
    return parser


def _write_csv(frame: pd.DataFrame, path: Path) -> None:
    result = frame.copy()
    for column in result:
        if result[column].map(lambda value: isinstance(value, (list, dict, tuple))).any():
            result[column] = result[column].map(lambda value: json.dumps(value, ensure_ascii=False) if isinstance(value, (list, dict, tuple)) else value)
    result.to_csv(path, index=False)


def _json_safe(value):
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if hasattr(value, "item"):
        value = value.item()
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    if isinstance(value, Path):
        return str(value)
    return value


def main(argv: list[str] | None = None) -> int:
    parser = _parser()
    args = parser.parse_args(argv)
    if bool(args.guides) == bool(args.example):
        parser.error("Supply a guide CSV or --example, but not both")
    if args.top < 1:
        parser.error("--top must be positive")
    if args.slide_min > args.slide_max:
        parser.error("--slide-min must be <= --slide-max")
    output_names = ["candidates.csv", "top_candidates.csv", "guide_summary.csv", "targets.csv", "targets.fasta", "report.html", "overview.png", "overview.svg", "manifest.json"]
    collisions = [name for name in output_names if (args.output / name).exists()]
    if collisions and not args.overwrite:
        parser.error("Output files already exist; choose another --output or pass --overwrite: " + ", ".join(collisions))
    try:
        from .guides import design_guides, guide_summary, top_guides
        from .guide_plotting import guide_overview, write_guide_report
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        if args.example:
            example = resources.files("probeharbour").joinpath("data/example_guides.csv")
            with example.open("r", encoding="utf-8") as handle:
                inputs = pd.read_csv(handle, keep_default_na=False)
            input_description = "packaged example_guides.csv"
        else:
            inputs = pd.read_csv(args.guides, keep_default_na=False)
            input_description = str(args.guides.resolve())
        store = design_guides(
            inputs, chemistry=args.chemistry, arm_length=args.arm_length,
            slides=range(args.slide_min, args.slide_max + 1), boundary=args.boundary,
            target_mode=args.target_mode, min_guide_overlap=args.min_guide_overlap,
            min_panel_mismatches=args.min_panel_mismatches, min_score=args.min_score,
            require_qc_pass=args.require_qc_pass,
        )
        shortlist = top_guides(store, n=args.top, eligible_only=True)
        summary = guide_summary(store)
        args.output.mkdir(parents=True, exist_ok=True)
        _write_csv(store.obs, args.output / "candidates.csv")
        _write_csv(shortlist, args.output / "top_candidates.csv")
        _write_csv(summary, args.output / "guide_summary.csv")
        targets = []
        fasta = []
        for index, record in enumerate(store.metadata["guides"], start=1):
            safe_id = re.sub(r"[^A-Za-z0-9_.-]", "_", record["guide_id"])
            target_id = f"guide_{index:03d}_{safe_id}"
            sequence = record["target_sequence"]
            targets.append({"target_id": target_id, "guide_id": record["guide_id"], "target_sequence": sequence})
            fasta.append(f">{target_id}\n" + "\n".join(sequence[start:start + 80] for start in range(0, len(sequence), 80)))
        pd.DataFrame(targets, columns=["target_id", "guide_id", "target_sequence"]).to_csv(args.output / "targets.csv", index=False)
        (args.output / "targets.fasta").write_text("\n".join(fasta) + ("\n" if fasta else ""), encoding="utf-8")
        write_guide_report(store, args.output / "report.html", n=args.top)
        figure = guide_overview(store, n=args.top)
        figure.savefig(args.output / "overview.png", dpi=180, bbox_inches="tight", facecolor="white")
        figure.savefig(args.output / "overview.svg", bbox_inches="tight", facecolor="white")
        plt.close(figure)
        versions = {}
        for package in ("probeharbour", "primer3-py", "pandas", "matplotlib"):
            try:
                versions[package] = metadata.version(package)
            except metadata.PackageNotFoundError:
                versions[package] = "not installed as distribution"
        manifest = {
            "created_at": datetime.now(timezone.utc).isoformat(), "input": input_description,
            "top_per_guide": args.top, "outputs": output_names, "versions": versions,
            "n_guides": len(store.metadata["guides"]), "n_candidates": len(store.obs),
            "n_eligible": int(store.obs["eligible"].sum()) if not store.obs.empty else 0,
            "metadata": store.metadata,
            "scope": "Research shortlist. Ungapped supplied-guide-panel sequence screen only; transcriptome specificity and RNA accessibility not assessed. No complete vendor oligos constructed.",
        }
        (args.output / "manifest.json").write_text(json.dumps(_json_safe(manifest), indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    except (ValueError, RuntimeError, OSError, ImportError) as exc:
        parser.exit(2, f"GuideHarbour: {exc}\n")
    print(f"Wrote {len(store.obs)} candidates and {len(shortlist)} shortlisted probes for {len(summary)} guides to {args.output.resolve()}")
    print(f"Open {args.output.resolve() / 'report.html'} to inspect the sequence windows and scores.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
