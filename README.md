# Brain Spatial Overlap Plotting Suite

A standalone, publication-grade neuroimaging visualization toolkit for mapping cross-trait spatial convergence, multi-trait vertex overlap, and genetic/morphometric statistics across cortical and subcortical brain structures.

This pipeline implements the spatial rendering framework and visual specifications developed by **Bingjing Zhou** for the SP-BWAS study, built on **`brainMapR`** ([Couvy-Duchesne et al.](https://github.com/couvydub/brainMapR)).

---

## 1. Overview and Visual Specifications

The toolkit generates standardized multi-panel composite brain maps adhering to publication specifications:

- **Eight-Column Anatomical Layout**:
  1. `Left CT`: Cortical Thickness, Left Hemisphere
  2. `Right CT`: Cortical Thickness, Right Hemisphere
  3. `Left SA`: Surface Area, Left Hemisphere
  4. `Right SA`: Surface Area, Right Hemisphere
  5. `Left Subcortical Thickness`: Planar projection, Left Hemisphere
  6. `Right Subcortical Thickness`: Planar projection, Right Hemisphere
  7. `Left Subcortical Log-Jacobian`: Planar projection, Left Hemisphere
  8. `Right Subcortical Log-Jacobian`: Planar projection, Right Hemisphere
- **Dual Perspectives**:
  - Upper row: Lateral / Outside views
  - Lower row: Medial / Inside views
- **Point-Sphere Visual Encodings**:
  - **Background vertices**: Neutral light grey (`#D0D0D0`, sphere radius = `0.95`, alpha = `1.0`), outlining full anatomical geometry across all non-significant vertices.
  - **Outcome FDR baseline**: Neutral beige (`#D7BE88`, sphere radius = `1.55`), denoting primary outcome FDR-significant locations with zero trait overlap.
  - **Single-trait overlap**: Designated trait-specific hue (`radius = 1.55`).
  - **Multi-trait overlap**: Saturated dark grey (`#333333`, `radius = 1.55`), denoting locations where $\ge 2$ traits overlap simultaneously.
  - **Trait count gradient**: Stepped warm palette (1, 2, 3, 4+ overlapping traits).

---

## 2. Input Data Specifications

Input tables must be formatted as gzip-compressed tab-separated values (`.tsv.gz`). Because the pipeline operates directly on the standard `brainMapR` template, rows must match the authoritative 653,994-vertex probe inventory.

### Required Columns

| Column | Type | Values | Description |
| :--- | :--- | :--- | :--- |
| `Probe` | Character | e.g. `lha_000`, `rht_102` | Vertex identifier matching standard brainMapR probe inventory (653,994 rows). |
| `ROIlabel` | Character | e.g. `precentral`, `Putamen` | Anatomical region label. |
| `hemi` | Character | `lh`, `rh` | Hemisphere identifier. |
| `moda` | Character | `thickness`, `area`, `thick`, `LogJacs` | Imaging modality. |
| `outcome_fdr_mask` | Logical | `TRUE`, `FALSE` | Outcome-level FDR significance indicator. |
| `<trait_id>` | Integer / Logical | `0`, `1` (or `TRUE`, `FALSE`) | Binary indicator of concordant significance for each trait. |
| `all_same` *(optional)* | Integer | $\ge 0$ | Total number of overlapping traits at the vertex. Computed automatically if omitted. |

The bundled domain example uses nonnegative counts within each domain, so its `all_same` can exceed the number of palette rows. For count-valued columns, provide `all_same` explicitly as the sum of those counts. Automatic calculation counts positive columns and is suitable for binary trait indicators.

### Data Preparation Snippet (R / `data.table`)

```R
library(data.table)

# Load existing brainMapR results table (653,994 rows)
df <- fread("brainmap_vertex_stats.tsv")

# Specify primary outcome FDR mask (e.g. AD or PD)
df[, outcome_fdr_mask := outcome_fdr < 0.05]

# Add trait indicators (1 for concordant significance, 0 otherwise)
# e.g., columns "mdd", "scz", "bip", "anx"
traits <- c("mdd", "scz", "bip", "anx")

# Optional: compute total overlap count
df[, all_same := rowSums(.SD > 0L), .SDcols = traits]

# Export compressed table
fwrite(df, "counts_trait_overlap.tsv.gz", sep = "\t", compress = "gzip")
```

---

## 3. Trait and Palette Configuration

Trait identifiers, labels, and display colors are specified in `config/trait_palette.tsv`. The pipeline dynamically adapts to any number of traits ($N \ge 1$):

```tsv
trait_id	trait_label	hex	order
mdd	Major Depression	#A23B72	1
scz	Schizophrenia	#3B6FA1	2
bip	Bipolar Disorder	#C9891B	3
anx	Anxiety	#23847B	4
```

---

## 4. Execution

The main renderer needs the `data.table`, `ggplot2`, `gridExtra`, `png`, `RColorBrewer`, and `rgl` R packages, plus the brainMapR atlas files. The `brainMapR` R package is not required.
The bundled `simulated_probe_domain_counts.tsv.gz` has 653,994 atlas-matched synthetic rows. The separate `example_overlap_data.tsv.gz` uses placeholder Probe IDs and cannot be joined to the atlas.
On QIMR, submit full rendering through PBS. The tested `R/4.3.1` module includes the rendering packages.

### Shell Runner

```bash
cd Brain_Spatial_Overlap_Plot
bash run_render.sh
```

### Direct Rscript Invocation

```bash
Rscript render_manuscript_spatial.R \
  --counts "example_data/simulated_probe_domain_counts.tsv.gz" \
  --outcome-id "Trait_Convergence" \
  --atlas-dir "/path/to/brainMapR/extdata/atlas" \
  --palette "config/trait_palette.tsv" \
  --output-dir "output_figures/run_01"
```

### Manuscript-density synthetic render

The default `RENDER_POINTS_PER_GROUP=0` draws every atlas-matched vertex. The bundled example has 149,953 left and 149,924 right vertices for each cortical modality, and 13,560 vertices per subcortical modality and hemisphere. The grey background sphere radius is 0.95, matching the manuscript script. Submit the full render to PBS on QIMR.

```bash
RENDER_POINTS_PER_GROUP=0 RENDER_SUPPLEMENTS=false \
  OUTCOME_ID=Simulated_Preview \
  OUTPUT_DIR=output_figures/preview \
  ATLAS_DIR=/path/to/brainMapR/extdata/atlas \
  bash run_render.sh
```

`RENDER_SUPPLEMENTS=false` limits this example to the main count and trait-identity panels. Leave it unset for all supplementary trait maps. For a quick input and layout check, `RENDER_POINTS_PER_GROUP=8000` is available, but its sparse cortex is not visually comparable with the manuscript.

### HPC Cluster Job Submission (PBS / QIMR `hpcapp01`)

```bash
#!/bin/bash
#PBS -N brain_spatial_render
#PBS -l select=1:ncpus=1:mem=24gb
#PBS -l walltime=04:00:00
#PBS -j oe

set -euo pipefail
module purge
module load R/4.3.1 ImageMagick/7.1.1-11

cd "${PBS_O_WORKDIR:-.}"

export ATLAS_DIR="/mnt/lustre/working/lab_puyag/bingjinZ/Side_Project_AD/software/Rlibs/brainMapR-076069ce/brainMapR/extdata/atlas"

bash run_render.sh
```

---

## 5. Output Figures

The rendering pipeline exports publication-ready composites and modular sub-panels in both vector PDF (`cairo_pdf`) and 500 DPI PNG formats:

- `Main_<outcome>_count_and_trait.png` / `.pdf`: Complete two-panel composite figure.
  - **Panel a**: Trait overlap count ladder (all observed positive counts).
  - **Panel b**: Trait identity map (trait-specific colors + multi-trait dark grey).
- `Panel_<outcome>_trait_count.png` / `.pdf`: Standalone Panel a.
- `Panel_<outcome>_trait_identity.png` / `.pdf`: Standalone Panel b.
- `Supplement_<outcome>_<trait>_overlap.png` / `.pdf`: Individual spatial convergence maps for each trait.
- `RENDERING_PARAMETERS.tsv`: Machine-readable summary of rendering parameters and sample statistics.

---

## 6. Reference Benchmarks

Authoritative benchmark figures matching the SP-BWAS manuscript are available under `output_figures/`:
- `reference_Figure_domain_identity_and_count.png`: Trait identity and multi-trait overlap layout.
- `reference_Figure_trait_overlap_count.png`: Stepped overlap count layout.
- `reference_Figure_section_breadth.png`: Multi-trait breadth layout.
