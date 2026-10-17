#!/bin/bash
# Render the bundled atlas-matched synthetic overlap table with the main pipeline.
# Input: example_data/simulated_probe_domain_counts.tsv.gz and ATLAS_DIR.
# Output: output_figures/demo (or OUTPUT_DIR).
# Run: ATLAS_DIR=/path/to/atlas bash example_data/run_demo.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export COUNTS="${COUNTS:-${SCRIPT_DIR}/example_data/simulated_probe_domain_counts.tsv.gz}"
export OUTPUT_DIR="${OUTPUT_DIR:-${SCRIPT_DIR}/output_figures/demo}"
export OUTCOME_ID="${OUTCOME_ID:-Simulated_Trait_Overlap}"

bash "${SCRIPT_DIR}/run_render.sh"
