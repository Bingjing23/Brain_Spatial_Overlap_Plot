#!/bin/bash
# Render manuscript spatial figures from an atlas-matched counts table.
# Inputs: COUNTS, PALETTE, ATLAS_DIR, OUTCOME_ID environment variables.
# Outputs: composite, panel, and trait figures under OUTPUT_DIR.
# Run: ATLAS_DIR=/path/to/atlas bash run_render.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
COUNTS="${COUNTS:-${SCRIPT_DIR}/example_data/simulated_probe_domain_counts.tsv.gz}"
OUTPUT="${OUTPUT_DIR:-${SCRIPT_DIR}/output_figures/simulated_run}"
PALETTE="${PALETTE:-${SCRIPT_DIR}/config/trait_palette.tsv}"
OUTCOME_ID="${OUTCOME_ID:-Simulated_Trait_Overlap}"
ATLAS="${ATLAS_DIR:-/working/lab_puyag/bingjinZ/Side_Project_AD/software/Rlibs/brainMapR-076069ce/brainMapR/extdata/atlas}"

mkdir -p "${OUTPUT}"
echo "Counts input: ${COUNTS}"
echo "Atlas dir: ${ATLAS}"
echo "Output dir: ${OUTPUT}"

render_command=(
  Rscript "${SCRIPT_DIR}/render_manuscript_spatial.R"
  --counts "${COUNTS}"
  --outcome-id "${OUTCOME_ID}"
  --atlas-dir "${ATLAS}"
  --palette "${PALETTE}"
  --output-dir "${OUTPUT}"
)

if command -v xvfb-run >/dev/null 2>&1; then
  xvfb-run -a "${render_command[@]}"
else
  "${render_command[@]}"
fi

echo "Figures generated under: ${OUTPUT}"
