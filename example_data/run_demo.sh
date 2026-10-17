#!/bin/bash
# ==============================================================================
# Example execution script for Brain_Spatial_Overlap_Plot
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${SCRIPT_DIR}"

# 1. Generate example dataset if not present
if [ ! -f "example_data/example_overlap_data.tsv.gz" ]; then
  echo "Generating sample dataset..."
  (cd example_data && Rscript generate_demo_dataset.R)
fi

# 2. Run plot_brain_spatial.R
# NOTE for HPC/Cluster environments (without physical display):
# Prefix command with xvfb-run to provide a virtual X11 server for rgl:
#   xvfb-run -a Rscript plot_brain_spatial.R ...

echo "Running brain spatial overlap visualization (Count mode)..."
Rscript plot_brain_spatial.R \
  --input example_data/example_overlap_data.tsv.gz \
  --output-dir output_figures \
  --prefix demo_trait_overlap_counts \
  --value-col overlap_count \
  --mode count \
  --surface inflated_pre \
  --modalities all \
  --title "Cross-Trait Vertex Overlap (Simulated Data)"

echo "Running brain spatial visualization (Continuous z-score / rGM mode)..."
Rscript plot_brain_spatial.R \
  --input example_data/example_overlap_data.tsv.gz \
  --output-dir output_figures \
  --prefix demo_trait_continuous_z \
  --value-col z_score \
  --mode continuous \
  --surface inflated_pre \
  --modalities all \
  --title "Trait Morphometric Association (Simulated Data)"

echo "Done! Outputs are saved in: ${SCRIPT_DIR}/output_figures"
