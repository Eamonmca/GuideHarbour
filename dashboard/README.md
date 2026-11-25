> This is the inherited exon-junction dashboard. For guide-RNA design, use `guideharbour-dashboard` and the [main README](../README.md).

# ProbeHarbour dashboard

Run the dashboard from the project root with the approved environment:

From a source checkout:

```bash
python -m streamlit run dashboard/app.py
```

From an installed package:

```bash
python -m pip install 'probeharbour[dashboard] @ git+https://github.com/Eamonmca/GuideHarbour.git@main'
probeharbour-dashboard
```

An alternative workflow-oriented UI demo is included without replacing the
current dashboard:

```bash
probeharbour-dashboard-demo
```

The demo adds a persistent five-stage workflow rail, a compact cart summary at
the top of the sidebar, and a restrained white/teal professional theme. Both
commands use the same ProbeHarbour API and session-state model.

Reference files can be selected in the sidebar or configured with the
`PROBEHARBOUR_GTF`, `PROBEHARBOUR_FASTA`, `PROBEHARBOUR_EXISTING_BED` and
`PROBEHARBOUR_BLAST_DB` environment variables.

The workflow is deliberately staged:

1. choose a chemistry, reference and gene;
2. inspect transcript models and select colour-coded annotated junctions;
3. design, characterize and score every candidate placement;
4. review and explicitly confirm candidate probes;
5. inspect library-generated gene maps and record cards;
6. export a ZIP containing CSV, BED, FASTA and a portable ProbeStore.

Transcriptome BLAST and comparison with an existing 10x BED are optional but
recommended before final export.

The junction explorer uses stable J1/J2 labels shared by the plot and table.
Optional CSV/TSV evidence can be joined with a required `junction_id` column;
additional columns such as `event`, `source`, `PMID`, `URL`, and `notes` are
displayed without changing the underlying annotation.
