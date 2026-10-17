#!/usr/bin/env Rscript

# ==============================================================================
# Standalone Cortical & Subcortical Vertex Overlap / Trait Visualization Suite
#
# Generates publication-grade cortical (inflated) and subcortical (flat projection)
# surface point maps for trait vertex overlap and BWAS / GWAS summary statistics.
# Standardized to match the Bingjing Zhou SP-BWAS & brainMapR visualization styles.
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(png)
  library(rgl)
})

# Determine script root directory
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, args)
  if (length(match) > 0) {
    return(dirname(normalizePath(sub(needle, "", args[match]))))
  }
  return(getwd())
}

script_dir <- get_script_dir()
source(file.path(script_dir, "lib/figure_style.R"))
source(file.path(script_dir, "lib/brain_render_engine.R"))

# Main R Function Interface
render_brain_spatial <- function(
  input_data,
  output_dir = "./output_figures",
  prefix = "brain_spatial_overlap",
  value_col = "overlap_count",
  mask_col = NULL,
  mode = c("count", "continuous", "binary"),
  atlas_dir = NULL,
  surface = "inflated_pre",
  modalities = c("all", "cortical", "subcortical"),
  title = "Cortical and Subcortical Vertex Overlap",
  bg_radius = DEFAULT_BG_RADIUS,
  fg_radius = DEFAULT_FG_RADIUS,
  fg_alpha  = DEFAULT_FG_ALPHA,
  pdf_width = 14.4,
  pdf_height = 4.2
) {
  mode <- match.arg(mode)
  modalities <- match.arg(modalities)
  atlas_dir <- find_atlas_dir(atlas_dir)
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  temp_comp_dir <- file.path(output_dir, paste0(".components_", prefix))
  dir.create(temp_comp_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Load or convert input table
  dt <- if (is.data.frame(input_data)) as.data.table(input_data) else fread(input_data)
  
  req_cols <- c("Probe", "hemi", "moda")
  for (col in req_cols) {
    if (!col %in% names(dt)) stop("Input data is missing required column: '", col, "'")
  }
  if (!value_col %in% names(dt)) {
    stop("Input data is missing specified value column: '", value_col, "'")
  }
  
  # Determine mask / foreground vertices
  if (!is.null(mask_col) && mask_col %in% names(dt)) {
    dt[, .fg_mask := as.logical(get(mask_col))]
  } else {
    if (mode == "count") {
      dt[, .fg_mask := get(value_col) > 0L]
    } else if (mode == "binary") {
      dt[, .fg_mask := as.logical(get(value_col))]
    } else {
      dt[, .fg_mask := !is.na(get(value_col)) & get(value_col) != 0]
    }
  }
  
  # Column specifications
  all_column_specs <- data.table(
    hemisphere = c("lh", "rh", "lh", "rh", "lh", "rh", "lh", "rh"),
    modality   = c("thickness", "thickness", "area", "area", "thick", "thick", "LogJacs", "LogJacs"),
    moda_group = c("cortical", "cortical", "cortical", "cortical", "subcortical", "subcortical", "subcortical", "subcortical"),
    label      = c("Left\nCT", "Right\nCT", "Left\nSA", "Right\nSA",
                   "Left\nSubc. Thick", "Right\nSubc. Thick",
                   "Left\nSubc. LogJac", "Right\nSubc. LogJac")
  )
  
  active_specs <- if (modalities == "all") {
    all_column_specs
  } else {
    all_column_specs[moda_group == modalities]
  }
  
  # Determine color mapping
  max_count <- 1L
  continuous_limit <- 1.0
  if (mode == "count") {
    max_count <- max(1L, as.integer(max(dt[[value_col]][dt$.fg_mask], na.rm = TRUE)))
    count_palette <- publication_count_palette(max_count)
  } else if (mode == "continuous") {
    continuous_limit <- as.numeric(quantile(abs(dt[[value_col]][dt$.fg_mask]), 0.99, na.rm = TRUE))
    if (!is.finite(continuous_limit) || continuous_limit <= 0) continuous_limit <- 1.0
  }
  
  message("--> Rendering view components using surface '", surface, "'...")
  
  # Subcortical atlas coordinates cache
  subcortical_atlas_file <- file.path(atlas_dir, "Atlas_coordinates_ROI_subortical.txt")
  subcortical_atlas <- if (file.exists(subcortical_atlas_file)) {
    sub_raw <- fread(subcortical_atlas_file, showProgress = FALSE)
    setnames(sub_raw, "ProbeID", "Probe")
    sub_raw
  } else NULL
  
  # Render all panels
  for (i in seq_len(nrow(active_specs))) {
    hemi_curr <- active_specs$hemisphere[i]
    moda_curr <- active_specs$modality[i]
    group_curr <- active_specs$moda_group[i]
    
    sub_dt <- dt[hemi == hemi_curr & moda == moda_curr]
    if (nrow(sub_dt) == 0L) {
      warning("No data rows found for modality: ", moda_curr, ", hemi: ", hemi_curr)
      next
    }
    
    # Merge with coordinates
    if (group_curr == "cortical") {
      atlas_file <- file.path(atlas_dir, paste0("Atlas_allCoordinates_corticalVertices_", hemi_curr, "_", moda_curr, ".txt"))
      if (!file.exists(atlas_file)) stop("Missing cortical atlas file: ", atlas_file)
      atlas_c <- fread(atlas_file, showProgress = FALSE)
      setnames(atlas_c, "vertexID", "Probe")
      coords_merged <- merge(sub_dt, atlas_c, by = "Probe")
      
      cols_xyz <- paste(c("Y", "X", "Z"), surface, sep = "_")
      if (!all(cols_xyz %in% names(coords_merged))) {
        stop("Surface columns '", paste(cols_xyz, collapse = ", "), "' not found in atlas.")
      }
      
      for (view in c("outside", "inside")) {
        xyz <- as.matrix(coords_merged[, ..cols_xyz])
        if (view == "outside") xyz[, 1:2] <- -xyz[, 1:2] # Lateral flip
        
        # Determine colors
        cols <- rep(DEFAULT_BG_COLOUR, nrow(coords_merged))
        fg_idx <- which(coords_merged$.fg_mask)
        if (length(fg_idx) > 0L) {
          if (mode == "count") {
            v <- pmax(1L, pmin(max_count, as.integer(coords_merged[[value_col]][fg_idx])))
            cols[fg_idx] <- count_palette[v]
          } else if (mode == "binary") {
            cols[fg_idx] <- "#E65100"
          } else {
            # Continuous diverging
            val_norm <- pmax(-continuous_limit, pmin(continuous_limit, coords_merged[[value_col]][fg_idx]))
            ratio <- (val_norm + continuous_limit) / (2 * continuous_limit)
            col_idx <- ceiling(ratio * 255L) + 1L
            cols[fg_idx] <- publication_diverging_palette(256L)[col_idx]
          }
        }
        
        comp_out <- file.path(temp_comp_dir, paste(moda_curr, hemi_curr, view, "png", sep = "_"))
        render_single_view(xyz, coords_merged$.fg_mask, cols, comp_out,
                           bg_radius = bg_radius, fg_radius = fg_radius,
                           fg_alpha = fg_alpha)
      }
    } else {
      # Subcortical
      if (is.null(subcortical_atlas)) stop("Subcortical atlas file not found.")
      coords_merged <- merge(sub_dt, subcortical_atlas, by = "Probe")
      coords_merged[, `:=`(X = -X, Y = -Y)]
      
      if (requireNamespace("brainMapR", quietly = TRUE)) {
        coords_merged <- as.data.table(brainMapR::addFlatCoordinatesSubcortical(annotBwas = coords_merged, hemi = hemi_curr))
      } else {
        coords_merged[, `:=`(Zf = Z, Xf = X, Yf = Y, Zf2 = Z, Xf2 = X, Yf2 = Y)]
      }
      
      for (view in c("inside", "outside")) {
        cols_xyz <- if (view == "inside") c("Zf", "Xf", "Yf") else c("Zf2", "Xf2", "Yf2")
        xyz <- as.matrix(coords_merged[, ..cols_xyz])
        
        cols <- rep(DEFAULT_BG_COLOUR, nrow(coords_merged))
        fg_idx <- which(coords_merged$.fg_mask)
        if (length(fg_idx) > 0L) {
          if (mode == "count") {
            v <- pmax(1L, pmin(max_count, as.integer(coords_merged[[value_col]][fg_idx])))
            cols[fg_idx] <- count_palette[v]
          } else if (mode == "binary") {
            cols[fg_idx] <- "#E65100"
          } else {
            val_norm <- pmax(-continuous_limit, pmin(continuous_limit, coords_merged[[value_col]][fg_idx]))
            ratio <- (val_norm + continuous_limit) / (2 * continuous_limit)
            col_idx <- ceiling(ratio * 255L) + 1L
            cols[fg_idx] <- publication_diverging_palette(256L)[col_idx]
          }
        }
        
        comp_out <- file.path(temp_comp_dir, paste(moda_curr, hemi_curr, view, "png", sep = "_"))
        render_single_view(xyz, coords_merged$.fg_mask, cols, comp_out,
                           bg_radius = bg_radius, fg_radius = fg_radius,
                           fg_alpha = fg_alpha)
      }
    }
  }
  
  message("--> Assembling publication grid figure...")
  
  # Assemble multi-panel view
  spaced_row <- function(left_grob, panel_grobs) {
    grobs <- list(left_grob)
    widths <- 0.60
    for (idx in seq_along(panel_grobs)) {
      grobs <- c(grobs, list(panel_grobs[[idx]]))
      widths <- c(widths, 1)
      if (idx < length(panel_grobs)) {
        gap <- if (idx %% 2L == 0L) 0.075 else 0.018
        grobs <- c(grobs, list(nullGrob()))
        widths <- c(widths, gap)
      }
    }
    arrangeGrob(grobs = grobs, nrow = 1L, widths = widths, padding = unit(1, "mm"))
  }
  
  header_grob <- spaced_row(
    nullGrob(),
    lapply(active_specs$label, function(lbl) {
      textGrob(lbl, y = unit(0.04, "npc"), just = "bottom",
               gp = gpar(fontsize = FIG_FONT$axis_text, fontface = "bold",
                         col = FIG_TEXT_COLOUR, fontfamily = FIG_FONT_FAMILY, lineheight = 0.9))
    })
  )
  
  row_outside <- spaced_row(
    textGrob("Lateral / Outside", rot = 90,
             gp = gpar(fontsize = FIG_FONT$axis_text, fontface = "bold",
                       col = FIG_TEXT_COLOUR, fontfamily = FIG_FONT_FAMILY)),
    lapply(seq_len(nrow(active_specs)), function(idx) {
      img_file <- file.path(temp_comp_dir, paste(active_specs$modality[idx], active_specs$hemisphere[idx], "outside_png", sep = "_"))
      image_file_to_grob(img_file)
    })
  )
  
  row_inside <- spaced_row(
    textGrob("Medial / Inside", rot = 90,
             gp = gpar(fontsize = FIG_FONT$axis_text, fontface = "bold",
                       col = FIG_TEXT_COLOUR, fontfamily = FIG_FONT_FAMILY)),
    lapply(seq_len(nrow(active_specs)), function(idx) {
      img_file <- file.path(temp_comp_dir, paste(active_specs$modality[idx], active_specs$hemisphere[idx], "inside_png", sep = "_"))
      image_file_to_grob(img_file)
    })
  )
  
  # Legend grob
  legend_grob <- if (mode == "count") {
    make_count_legend_grob(max_count, count_palette)
  } else if (mode == "continuous") {
    make_continuous_legend_grob(continuous_limit, title = value_col)
  } else {
    nullGrob()
  }
  
  # Full composite figure
  title_grob <- textGrob(title, x = unit(0.01, "npc"), just = "left",
                         gp = gpar(fontsize = FIG_FONT$figure_title, fontface = "bold", fontfamily = FIG_FONT_FAMILY))
  
  full_plot <- arrangeGrob(
    title_grob,
    header_grob,
    row_outside,
    row_inside,
    legend_grob,
    ncol = 1L,
    heights = c(0.12, 0.20, 1.0, 1.0, 0.25),
    padding = unit(1, "mm")
  )
  
  pdf_path <- file.path(output_dir, paste0(prefix, ".pdf"))
  png_path <- file.path(output_dir, paste0(prefix, ".png"))
  
  cairo_pdf(pdf_path, width = pdf_width, height = pdf_height, family = FIG_FONT_FAMILY)
  grid.draw(full_plot)
  dev.off()
  
  png(png_path, width = round(pdf_width * 300), height = round(pdf_height * 300), res = 300, bg = "white")
  grid.draw(full_plot)
  dev.off()
  
  message("==> Output written to:")
  message("    PDF: ", pdf_path)
  message("    PNG: ", png_path)
  
  return(invisible(list(pdf = pdf_path, png = png_path)))
}

# ==============================================================================
# CLI Entry Point
# ==============================================================================
if (sys.nframe() == 0L) {
  parse_cli_args <- function(args) {
    keys <- c("input", "output_dir", "prefix", "value_col", "mask_col",
              "mode", "atlas_dir", "surface", "modalities", "title")
    res <- list(
      input = NULL,
      output_dir = "./output_figures",
      prefix = "brain_spatial_overlap",
      value_col = "overlap_count",
      mask_col = NULL,
      mode = "count",
      atlas_dir = NULL,
      surface = "inflated_pre",
      modalities = "all",
      title = "Cortical and Subcortical Vertex Overlap"
    )
    i <- 1L
    while (i <= length(args)) {
      key <- sub("^--", "", gsub("-", "_", args[i]))
      if (key %in% keys) {
        res[[key]] <- args[i + 1L]
        i <- i + 2L
      } else {
        warning("Unknown argument: ", args[i])
        i <- i + 1L
      }
    }
    res
  }
  
  cli_opts <- parse_cli_args(commandArgs(trailingOnly = TRUE))
  if (is.null(cli_opts$input)) {
    cat("\nUsage:\n")
    cat("  Rscript plot_brain_spatial.R --input <data.tsv> [options]\n\n")
    cat("Options:\n")
    cat("  --input        <path>   Input vertex table (.tsv, .csv, or .gz) [Required]\n")
    cat("  --output-dir   <path>   Output directory [default: ./output_figures]\n")
    cat("  --prefix       <string> File prefix for outputs [default: brain_spatial_overlap]\n")
    cat("  --value-col    <string> Column name of overlap count or stat [default: overlap_count]\n")
    cat("  --mask-col     <string> Optional binary/logical mask column [default: value > 0]\n")
    cat("  --mode         <string> 'count', 'continuous', or 'binary' [default: count]\n")
    cat("  --atlas-dir    <path>   brainMapR atlas directory [autodetected if omitted]\n")
    cat("  --surface      <string> 'inflated_pre', 'inflated', 'pial' [default: inflated_pre]\n")
    cat("  --modalities   <string> 'all', 'cortical', or 'subcortical' [default: all]\n")
    cat("  --title        <string> Figure header title\n\n")
    quit(status = 1L)
  }
  
  render_brain_spatial(
    input_data = cli_opts$input,
    output_dir = cli_opts$output_dir,
    prefix     = cli_opts$prefix,
    value_col  = cli_opts$value_col,
    mask_col   = cli_opts$mask_col,
    mode       = cli_opts$mode,
    atlas_dir  = cli_opts$atlas_dir,
    surface    = cli_opts$surface,
    modalities = cli_opts$modalities,
    title      = cli_opts$title
  )
}
