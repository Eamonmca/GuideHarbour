# Guide detection: design choices and interpretation

GuideHarbour adapts ProbeHarbour to rank **detection probes for supplied guide RNAs**. It does not design CRISPR spacers, predict editing or repression efficiency, or assess genomic editing off-targets. A label such as `HIF1A-sg1` identifies the guide-bearing target; it does not mean the probe measures endogenous HIF1A RNA.

The workflow starts from an explicitly supplied sequence, enumerates adjacent two-arm binding windows, evaluates their inherited ProbeHarbour chemistry score, and screens their sequence similarity against the other supplied guide targets. `ph.design_guides` returns a `ProbeStore`; `ph.top_guides` returns ranked alternatives; `ph.guide_summary` provides a per-guide overview. The report shows where candidates bind and why they rank as they do. See the main [README](../README.md) for installation and executable examples.

## 1. What sequence is being detected?

The example input preserves the 15 guide sequences, common prefix, and common scaffold supplied for this project. Thirteen guides are 20 bases long. `SNHG12-sg3` and `Neg sg2` are 21 bases long. The common prefix is 19 bases, and the scaffold is 75 bases. These are counts of the supplied strings, not a correction or validation of the original constructs.

| Mode | Target sequence used for design | Intended interpretation |
| --- | --- | --- |
| `rna` (default) | `transcribed_leader + guide_sequence + scaffold_sequence` | A known guide-containing RNA region. The supplied promoter/cloning prefix is retained as metadata and excluded from this target. |
| `full_construct` | `prefix_sequence + guide_sequence + scaffold_sequence` | The exact supplied construct region. Use only when this complete region is the intended assay target. A prefix supplied in a cloning table is not evidence that it occurs in the RNA being measured. |

No extra initiating G is inferred, and no 21-base guide is trimmed to 20 bases. If a known transcribed leader is relevant, supply it explicitly in `transcribed_leader`. `full_construct` rejects a nonempty leader to avoid ambiguous duplication; it uses exactly the prefix, guide, and scaffold. Input sequences are normalized to uppercase DNA letters, with U represented as T and whitespace removed; the original inputs remain recorded in the run metadata. A target region need not describe the entire mature RNA: windows inside the guide and scaffold remain valid substrings if the actual RNA also contains a verified upstream leader.

### Why the prefix is excluded by default

The supplied prefix `TGGAAAGGACGAAACACCG` is also used as a cloning homology arm in published guide-cassette workflows. Its appearance beside a guide does not establish the 5′ boundary of the transcribed guide RNA. A U6-driven guide and a longer Pol-II transcript carrying the cassette can have different boundaries. The default therefore targets the guide and scaffold, while retaining an explicit option for a verified longer construct. See the [AAV-Perturb-seq methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC10567566/) and [GiRAFR/CROP-seq methods](https://www.nature.com/articles/s42003-023-05351-7) for examples using this exact prefix.

The two `Neg` guides are negative controls for the biological perturbation. Their detection probes should still detect the corresponding guide-containing target when it is present. They are not automatically negative controls for probe hybridization.

## 2. From a guide boundary to candidate probes

The same local-window principle used for exon junctions is applied to a guide/context boundary. Each candidate contains a target-oriented upstream arm (`rbd5`) and an immediately adjacent downstream arm (`rbd3`). Their concatenation, `binder`, is a continuous substring of the chosen target.

| Boundary setting | Boundary searched |
| --- | --- |
| `guide_scaffold` | End of the guide and start of the scaffold; the default. |
| `prefix_guide` | End of the available upstream context and start of the guide. In RNA mode, available upstream context is the explicit transcribed leader, not the excluded cloning prefix. |
| `both` | Both boundaries. Each input must contain a nonempty upstream context in the chosen target mode. |

The `prefix_guide` and `both` settings require a nonempty modeled upstream context for every input guide; they raise an error when that context is absent. The ligation position slides around the selected boundary. Windows extending beyond the supplied target are not padded with invented bases. This matters near the 5′ end: a 20-base guide with no leader cannot support a 20-base upstream arm shifted before the guide/scaffold boundary. A 16-base upstream arm allows more upstream shifts, but it does not create additional target sequence. Arm lengths default to 16 bases with the packaged `prime_5k` rules and 20 bases with `xenium_v1`; an explicit `arm_length` overrides this project setting.

Candidates must overlap the guide by at least **15 bases by default**. This is a configurable project heuristic intended to keep the limited guide-specific sequence represented. It is not an experimentally calibrated detection threshold or a manufacturer specification.

### Why guide overlap needs its own rule

All 15 examples share the same scaffold. A candidate shifted far enough into that scaffold can acquire favorable ligation properties while losing guide-specific bases. For example, with 16-base arms, shifting five bases downstream of the guide/scaffold boundary leaves only 11 guide bases within the 32-base binder. A chemistry score alone cannot establish whether this is a good guide discriminator.

Guide coverage, boundary position, chemistry measurements, and competing guide matches should therefore be inspected together. A high score does not compensate for targeting only shared sequence.

## 3. Sequence orientation and physical probe arms

`rbd5`, `rbd3`, and `binder` describe the **target sequence in 5′→3′ transcription orientation**, written with DNA letters. This preserves the upstream ProbeHarbour convention. The ligation dinucleotide is the last base of target `rbd5` followed by the first base of target `rbd3`.

The physical DNA binding arms that hybridize to this RNA are reverse complements of their corresponding target segments. These sequences must be distinguished from the sense target strings in any export or visualization. In particular:

- `rbd5` and `rbd3` are target-oriented fields; they are not complete synthetic oligos.
- The reverse complement of `rbd5 + rbd3` is `reverse_complement(rbd3) + reverse_complement(rbd5)`. Complementation reverses the order.
- A padlock probe contains its platform-specific backbone and terminal chemistry as well as its binding arms. A reverse-complement binder is not, by itself, an order-ready padlock design.
- The inherited chemistry calculations on target-oriented strings and any separately reported physical-arm calculations must retain distinct labels. Neither models a complete assembled probe unless the complete probe sequence has actually been supplied.

Historical 10x Xenium v1 documentation defines its target-arm input in transcription orientation and describes reverse complementation during padlock construction. See [CG000683 Rev D, orientation section](https://cdn.10xgenomics.com/image/upload/v1716400584/support-documents/CG000683_TechNote_Xenium_Custom_Panel_Design_RevD.pdf). GuideHarbour follows the target-oriented convention for compatibility; this does not establish vendor acceptance or experimental performance.

## 4. What the chemistry score means

The guide workflow inherits ProbeHarbour's versioned chemistry rules and characterization machinery. Its upstream baseline is commit [`a4bee8077e2ced2135ee7b068a4c6715573bfc6e`](https://github.com/cobioda/ProbeHarbour/tree/a4bee8077e2ced2135ee7b068a4c6715573bfc6e).

The inherited default `ThermoParameters` profile uses:

| Parameter | Inherited value | Primer3 unit/convention |
| --- | --- | --- |
| `mv_conc` | 80 | mM monovalent cations |
| `dv_conc` | 5 | mM divalent cations |
| `dntp_conc` | 0 | mM dNTPs |
| `dna_conc` | 1500 | nM DNA concentration |
| `salt_corrections_method` | `schildkraut` | Salt-correction method |
| `structure_temp_c` | 50 | °C for structure calculations |

These are **inherited software defaults, not an official or calibrated guide-detection assay profile**. Historical Xenium v1 guidance used a different Primer3 concentration profile. Scores or temperatures calculated with different profiles should not be compared as if they were obtained under the same conditions. Record both the chemistry/ruleset identity and thermodynamic parameters when comparing runs.

The packaged rules start from the ligation dinucleotide's base score and apply named modifiers. With the inherited default scoring strategy:

```text
score = ligation_base_score / (1 + number_of_active_modifiers)
```

The active modifiers depend on the selected ruleset. They may include temperature, repeat, and predicted structure criteria. This is an explainable ranking heuristic. It is **not a probability of detection, a measured ligation efficiency, or a certificate that every chemistry constraint passed**. A candidate can retain a nonzero score despite an active modifier. Inspect modifier names and measurements alongside the ranking. The API's `require_qc_pass=True` option makes active chemistry modifiers a reason for exclusion; its default is `False`. The `chemistry_qc_pass` field means no packaged modifier was active, not that the assay has been validated.

Primer3 characterization also does not directly establish the accessibility of the RNA target. In particular, folding of the guide scaffold, protein occupancy, fixation, and assay conditions can affect access to the target. [10x's guide-capture guidance for Flex](https://kb.10xgenomics.com/s/article/How-do-I-design-custom-probes-for-CRISPR-Guide-Capture-for-use-with-the-GEM-X-Flex-v2-Assay) describes scaffold-structure effects on binding. That is useful evidence for this limitation; its probe architecture and optimization choices are not automatically transferable to Xenium.

## 5. Specificity against the supplied guide panel

GuideHarbour compares the complete binder with possible contiguous windows in the other supplied guide targets. This is a local **intra-panel sequence screen**. It is not a transcriptome-wide search or an experimental measurement of cross-binding.

The default minimum separation is **three mismatches across the full binder**. This is a configurable project heuristic. It is not a manufacturer-approved or calibrated specificity threshold. Mismatch count alone ignores binding thermodynamics, mismatch identities, gaps and bulges, and the positional importance of mismatches near the ligation junction. A threshold of three must not be described as satisfying historical guidance that recommended four mismatches for a different screening context.

The intended target is identified by its guide record. A different guide ID with an identical target sequence remains a competing target and must appear as an exact-match collision. Exact matches are excluded even if the configured mismatch threshold is zero. Guide IDs also keep two guides for the same gene separate. Reported mismatch positions are zero-based offsets within the complete binder; separate upstream and downstream mismatch counts identify which arm contains the differences. Only other targets long enough to contain the complete binder contribute full-length comparisons; the number screened and the screen status make missing comparisons visible.

Interpret panel evidence in the context of the entire supplied pool:

- A perfect competing match is a direct ambiguity between supplied targets.
- A near match is a reason to inspect the competing guide, alignment position, and which arm carries the differences.
- A small mismatch count across a binder dominated by common scaffold deserves particular scrutiny.
- A single-guide input cannot demonstrate discrimination from other guides that were never supplied.
- Adding another guide can change the specificity evidence and eligibility of existing candidates.

Endogenous RNA, vector-derived transcripts, other guide libraries, and contaminating nucleic acids are outside this screen unless explicitly included in a suitable reference analysis. A guide's genomic off-target editing predictions cannot substitute for a detection-probe transcriptome screen. No transcriptome specificity is established by this guide-only workflow.

## 6. Reading the top-candidate output

`ph.top_guides` selects ranked alternatives per guide from the computed candidate set and collapses duplicate binding sequences for the same guide. Its eligibility decisions and the chemistry score answer different questions. A candidate may meet the selected guide-overlap and intra-panel-distance settings while still having chemistry modifiers or untested transcriptome specificity. The default `eligible_only=True` does not fill an empty shortlist with failed candidates. Setting it to `False` deliberately exposes excluded candidates with their reasons.

The top alternatives around one guide usually overlap extensively. They should not be interpreted as a validated collection of probes that can all bind simultaneously. If the experimental plan needs multiple probes per target, evaluate spacing, competition, and the actual platform's panel constraints separately. Reporting fewer than the requested number of alternatives, or no eligible alternative for one guide, is preferable to silently relaxing a setting.

For each shortlisted candidate, inspect its placement, retained guide bases, target-arm sequences, temperature measurements, ligation dinucleotide, active modifiers, and nearest competing target. `Neg sg2` contains a long G run in its supplied guide; candidates covering that run deserve attention to repeat flags. A ranking can help prioritize a pilot, but cannot replace verification of the construct sequence and experimental detection.

## 7. Exports and current vendor workflow

Analysis CSVs and visual reports preserve candidate sequences and ranking evidence. They are useful for review and reproducibility, and are not claims of manufacturer acceptance or ordering readiness.

As accessed on 24 November 2026, 10x's [advanced custom sequence guidance](https://www.10xgenomics.com/support/software/xenium-panel-designer/latest/tutorials/advanced-custom-designs) requests a sense-orientation FASTA for CRISPR guide targets, including the protospacer or another uniquely targetable region and at least 80 bases of sequence. The example guide+scaffold target regions here contain 95 or 96 bases. That length comparison only addresses one input characteristic; it does not establish acceptance, sensitivity, or specificity.

The current [panel-design options page](https://www.10xgenomics.com/support/software/xenium-panel-designer/latest/tutorials/choose-design) describes direct guide detection and indirect detection through co-expressed barcodes. The former direct RBD submission instructions are deprecated according to the [technical-note landing page](https://www.10xgenomics.com/support/instruments/xenium-analyzer/species-standalone-custom-and-advanced-custom-panel-design-for-xenium-in-situ). Use current vendor instructions for submission rather than treating a historical three-column RBD CSV as a guaranteed current interface.

## 8. Scope of this implementation

This version provides reproducible sequence-window enumeration, inherited chemistry characterization and scoring, discrimination checks within the supplied guide panel, and visual comparison of candidates. It does not establish transcriptome-wide specificity, RNA accessibility, in situ retention of guide RNAs, sensitivity, or wet-lab validation. It also does not establish the complete mature transcript sequence from the table's column headings.

The guide-specific source modules and this methods document are new files. They are maintained directly, separately from the inherited nbdev-generated modules and their source notebooks. Preserve that distinction when extending the fork so that exporting upstream notebooks does not silently erase guide-specific changes.
