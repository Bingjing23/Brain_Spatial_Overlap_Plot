#!/bin/bash
# ==============================================================================
# Execution script for render_manuscript_spatial.R
# Uses Bingjing's EXACT manuscript spatial rendering pipeline with simulated data.
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

COUNTS="${SCRIPT_DIR}/example_data/simulated_probe_domain_counts.tsv.gz"
OUTPUT="${SCRIPT_DIR}/output_figures/simulated_run"
PALETTE="${SCRIPT_DIR}/config/trait_palette.tsv"

# brainMapR atlas directory containing coordinate text files
# Default QIMR HPC path or user-specified ATLAS_DIR
ATLAS="${ATLAS_DIR:-/mnt/lustre/working/lab_puyag/bingjinZ/Side_Project_AD/software/Rlibs/brainMapR-076069ce/brainMapR/extdata/atlas}"

mkdir -p "${OUTPUT}"

echo "Starting manuscript spatial rendering pipeline..."
echo "  Counts input : ${COUNTS}"
echo "  Atlas dir    : ${ATLAS}"
echo "  Output dir   : ${OUTPUT}"

# Set option for binary domain alternative if desired:
# export BINARY_DOMAIN_ALTERNATIVE=true

CMD="Rscript render_manuscript_spatial.R \
  --counts \"${COUNTS}\" \
  --outcome-id \"Simulated_Trait_Overlap\" \
  --atlas-dir \"${ATLAS}\" \
  --palette \"${PALETTE}\" \
  --output-dir \"${OUTPUT}\""

if command -v xvfb-run >/dev/null 2>&1; then
  echo "Executing with virtual X11 server (xvfb-run)..."
  xvfb-run -a ${CMD}
else
  echo "Executing Rscript directly..."
  ${CMD}
fi

echo "Done! Outputs generated under: ${OUTPUT}"
