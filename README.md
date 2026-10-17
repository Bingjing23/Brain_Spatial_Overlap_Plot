# Brain Spatial Overlap Plot

Render manuscript-style cortical and subcortical maps of cross-trait spatial overlap. The main figure shows overlap count and trait identity across eight modality–hemisphere columns, each with outside and inside views. Cortical coordinates use the `inflated_pre` template; subcortical coordinates use flattened projections adapted from [brainMapR](https://github.com/couvydub/brainMapR).

## Quick start

The bundled `example_data/simulated_probe_domain_counts.tsv.gz` is an atlas-matched synthetic input with 653,994 vertices. The atlas coordinate files are required separately. On QIMR, save the following as `render.pbs` in the repository and submit it with `qsub render.pbs`; adjust the resource request to your cluster:

```bash
#!/bin/bash
#PBS -l select=1:ncpus=1:mem=24gb
#PBS -l walltime=08:00:00
set -euo pipefail
module load R/4.3.1 ImageMagick/7.1.1-11
cd "$PBS_O_WORKDIR"
export ATLAS_DIR=/path/to/brainMapR/extdata/atlas
bash run_render.sh
```

A full render can take several hours. `run_render.sh` uses `xvfb-run` when available for headless `rgl` rendering. Required R packages are `data.table`, `ggplot2`, `gridExtra`, `png`, `RColorBrewer`, and `rgl`. The `brainMapR` R package itself is not required.

For a quick input and layout check, set `RENDER_POINTS_PER_GROUP=8000 RENDER_SUPPLEMENTS=false`. This sparse preview is **not** visually comparable with the manuscript. The default `RENDER_POINTS_PER_GROUP=0` draws every matched vertex; the grey background sphere radius remains `0.95`.

## Input format

Pass your own gzip-compressed TSV through `COUNTS=/path/to/counts.tsv.gz`. The file must contain exactly 653,994 unique atlas `Probe` IDs and these columns:

| Column | Meaning |
| --- | --- |
| `Probe`, `ROIlabel`, `hemi`, `moda` | Atlas identifiers; `hemi` is `lh` or `rh`, and `moda` is `thickness`, `area`, `thick`, or `LogJacs`. |
| `outcome_fdr_mask` | `TRUE`/`FALSE` indicator for outcome-level FDR significance. |
| One column per trait | Binary concordant-overlap indicator (`0`/`1`), named either `<trait_id>` or `<trait_id>__same`. |
| `all_same` (optional) | Number of overlapping traits at each vertex. Derived from positive trait indicators if omitted. |

Define matching `trait_id`, `trait_label`, `hex`, and `order` values in `config/trait_palette.tsv`, then set `PALETTE=/path/to/trait_palette.tsv`. The concordance and significance tests that produce the input indicators are performed upstream of this plotting repository. If trait columns hold counts rather than binary indicators, provide `all_same` explicitly as the sum of those counts.

Set `OUTCOME_ID=Your_Outcome` to label the output. By default, `run_render.sh` uses the bundled synthetic input and palette. The separate `example_data/example_overlap_data.tsv.gz` contains placeholder Probe IDs and is not compatible with the manuscript renderer.

## Figures and settings

Outputs are written under `OUTPUT_DIR`:

- `Main_<outcome>_count_and_trait.png` and `.pdf`: two-panel count and trait-identity figure.
- `Panel_<outcome>_trait_count.*` and `Panel_<outcome>_trait_identity.*`: individual panels.
- `Supplement_<outcome>_<trait>_overlap.*`: trait-specific maps. Set `RENDER_SUPPLEMENTS=false` to omit them.
- `RENDERING_PARAMETERS.tsv`: recorded plotting settings.

PNG files are exported at 500 DPI. The PDFs embed rasterized `rgl` brain views. Colours and point radii are defined near the top of `render_manuscript_spatial.R`; the default settings match the manuscript. The files in `output_figures/reference_Figure_*.png` are manuscript reference images, while the bundled input table is synthetic.
