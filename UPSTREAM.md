# Upstream provenance

- Source repository: https://github.com/cobioda/ProbeHarbour
- Source branch: `main`
- Source commit: `a4bee8077e2ced2135ee7b068a4c6715573bfc6e`
- Source retrieved: 24 November 2026
- Original license: Apache License 2.0; retained in `LICENSE`.

GuideHarbour is a separate personal adaptation for guide-RNA detection probe design. The upstream implementation and its history were cloned before adding the guide-specific workflow. The original README is preserved as `README_ProbeHarbour.md`.

New code adds guide/context target validation, candidate enumeration, intra-panel comparisons, per-guide ranking, visual reports, a command-line interface and a local dashboard. The upstream `Candidate`, `ProbeStore`, Primer3 characterization and packaged ligation/chemistry rules are reused. Existing exon-junction workflows remain available.

The guide input example is the 15-row table supplied for this adaptation. The guide/scaffold sequences, 19 nt construct prefix and original guide labels are preserved. Public distribution of project-specific guide inputs should be chosen by the repository owner.

The import and distribution name stays `probeharbour` for compatibility with the inherited API; the fork version is `0.2.0.dev0`. Use a separate virtual environment for upstream and this adaptation.

The new guide modules are normal Python source, not nbdev exports. No new notebook is a source of truth for these modules; the example notebook calls their public API.
