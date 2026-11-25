#!/usr/bin/env bash
set -euo pipefail

# Canonical ProbeHarbour dashboard launcher for the shared analysis workstation.
# Every path may be overridden by exporting the corresponding variable first.

PROBEHARBOUR_LAUNCH_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROBEHARBOUR_REPO_DIR="$(cd "${PROBEHARBOUR_LAUNCH_DIR}/.." && pwd)"
PROBEHARBOUR_PYTHON="${PROBEHARBOUR_PYTHON:-/data/analysis/data_mcandrew/Allos_new/allos_env/bin/python}"

export PROBEHARBOUR_GTF="${PROBEHARBOUR_GTF:-/data/analysis/data_mcandrew/00_allos_manuscript/allos/data/Mus_musculus.GRCm38.102.gtf.gz}"
export PROBEHARBOUR_FASTA="${PROBEHARBOUR_FASTA:-/data/analysis/data_mcandrew/00_allos_manuscript/allos/data/Mus_musculus.GRCm38.dna.primary_assembly.fa}"
export PROBEHARBOUR_EXISTING_BED="${PROBEHARBOUR_EXISTING_BED:-/data/analysis/data_mcandrew/00_isoprobe/10x_panel/5K/xenium_mouse_5K_gene_expression_panel_probe_locations.bed}"
export PROBEHARBOUR_BLAST_DB="${PROBEHARBOUR_BLAST_DB:-/data/analysis/data_mcandrew/ProbeHarbour/data/blast_dbs/GRCm38_ens102_cdna}"
export PYTHONPATH="${PROBEHARBOUR_REPO_DIR}${PYTHONPATH:+:${PYTHONPATH}}"
export PYTHONNOUSERSITE=1
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-/tmp/probeharbour-cache}"

for PROBEHARBOUR_REQUIRED_PATH in \
    "${PROBEHARBOUR_PYTHON}" \
    "${PROBEHARBOUR_GTF}" \
    "${PROBEHARBOUR_FASTA}" \
    "${PROBEHARBOUR_EXISTING_BED}"; do
    if [[ ! -e "${PROBEHARBOUR_REQUIRED_PATH}" ]]; then
        echo "ProbeHarbour launch error: required path not found: ${PROBEHARBOUR_REQUIRED_PATH}" >&2
        exit 1
    fi
done

exec "${PROBEHARBOUR_PYTHON}" -m streamlit run \
    "${PROBEHARBOUR_LAUNCH_DIR}/app_demo.py" \
    --server.address 0.0.0.0 \
    --server.port "${PROBEHARBOUR_PORT:-8501}" \
    --server.headless true \
    --theme.base light \
    --theme.primaryColor "#008F91" \
    --theme.backgroundColor "#FFFFFF" \
    --theme.secondaryBackgroundColor "#F5F9FA" \
    --theme.textColor "#08224C"
