#!/usr/bin/env Rscript

# Render final manuscript spatial figures from authoritative sign-aligned counts.
# Inputs are immutable derived count tables; outputs use inflated_pre cortex and
# the approved opaque point-sphere style for all cortical/subcortical panels.

suppressPackageStartupMessages({
  library(brainMapR)
  library(data.table)
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(png)
  library(RColorBrewer)
  library(rgl)
})
setDTthreads(1L)
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  match <- grep("--file=", args)
  if (length(match) > 0) {
    return(dirname(normalizePath(sub("--file=", "", args[match]))))
  }
  return(getwd())
}
script_dir <- get_script_dir()
source(file.path(script_dir, "lib/manuscript_figure_style.R"))

parse_args <- function(x) {
  required_keys <- c("counts", "outcome_id", "atlas_dir", "output_dir")
  keys <- c(required_keys, "palette")
  if (length(x) %% 2L != 0L) stop("Arguments must be --name value pairs")
  out <- setNames(as.list(rep(NA_character_, length(keys))), keys)
  
  # Default to trait_palette.tsv (or legacy manuscript_domain_palette.tsv)
  default_palette <- file.path(script_dir, "config/trait_palette.tsv")
  if (!file.exists(default_palette)) {
    default_palette <- file.path(script_dir, "config/manuscript_domain_palette.tsv")
  }
  out$palette <- default_palette

  for (i in seq(1L, length(x), by = 2L)) {
    key <- gsub("-", "_", sub("^--", "", x[[i]]))
    if (!key %in% keys) stop("Unknown argument: ", x[[i]])
    out[[key]] <- x[[i + 1L]]
  }
  if (anyNA(unlist(out[required_keys]))) stop("Missing required argument")
  out
}
opt <- parse_args(commandArgs(trailingOnly = TRUE))
display_outcome <- if (opt$outcome_id == "ADvsHC") "AD vs HC" else if (opt$outcome_id == "Conversion5years") "5-year conversion" else opt$outcome_id
if (!file.exists(opt$counts)) stop("Missing counts: ", opt$counts)
if (!dir.exists(opt$atlas_dir)) stop("Missing atlas directory: ", opt$atlas_dir)
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 🎨 1. COLOR & VISUAL STYLING PARAMETERS (Modify here if needed)
# ==============================================================================
# Brain geometry surface style: "inflated_pre" is the standard Freesurfer inflated cortex
style <- "inflated_pre"

# Neutral background vertices (outlines full brain anatomy for all non-significant probes):
background_colour   <- "#D0D0D0" # Neutral light grey
background_radius   <- 0.95      # Point sphere radius (smaller for dense cortical mesh)

# Outcome FDR-significant vertices with zero overlapping traits:
outcome_only_colour <- "#D7BE88" # Neutral beige (distinguishes disease effect from trait overlap)
significant_radius  <- if (opt$outcome_id == "Conversion5years") 1.85 else 1.55 # Sphere radius

# Multi-trait overlap vertices (where 2 or more traits overlap simultaneously):
multi_trait_colour  <- "#333333" # Saturated dark grey

# Sphere opacity (1 = 100% solid/opaque spheres, matching publication contract):
point_alpha         <- 1.0

# ==============================================================================
# 🏷️ 2. TRAIT PALETTE & CONFIGURATION (Single-Level Traits)
# ==============================================================================
palette_path <- normalizePath(opt$palette, mustWork = TRUE)
palette_contract <- fread(palette_path)

# Support both trait_id (new standard) and domain_id (legacy fallback)
if ("domain_id" %in% names(palette_contract) && !"trait_id" %in% names(palette_contract)) {
  setnames(palette_contract, "domain_id", "trait_id")
}
if ("domain_label" %in% names(palette_contract) && !"trait_label" %in% names(palette_contract)) {
  setnames(palette_contract, "domain_label", "trait_label")
}

n_traits <- nrow(palette_contract)
if (n_traits < 1L) stop("Trait palette must contain at least 1 trait entry")
if (!identical(palette_contract$order, seq_len(n_traits)) ||
    uniqueN(palette_contract$trait_id) != n_traits ||
    uniqueN(palette_contract$trait_label) != n_traits ||
    any(!grepl("^#[0-9A-Fa-f]{6}$", palette_contract$hex))) {
  stop("Invalid trait palette: order must be 1..N with unique IDs and valid hex codes")
}

traits <- palette_contract$trait_id
trait_labels <- setNames(palette_contract$trait_label, palette_contract$trait_id)
trait_colours <- setNames(palette_contract$hex, palette_contract$trait_id)
file.copy(palette_path, file.path(opt$output_dir, "TRAIT_PALETTE.tsv"), overwrite = TRUE)

typo_src <- file.path(script_dir, "config/manuscript_figure_typography.tsv")
if (file.exists(typo_src)) {
  file.copy(typo_src, file.path(opt$output_dir, "FIGURE_TYPOGRAPHY.tsv"), overwrite = TRUE)
}

modality_labels <- c(
  thickness = "Cortical thickness",
  area = "Surface area",
  thick = "Local thickness",
  LogJacs = "Log-Jacobian"
)
cortical_modalities <- c("thickness", "area")

# ==============================================================================
# 📥 3. INPUT DATA LOADING (Supports "<trait>" or "<trait>__same")
# ==============================================================================
# Load input counts table
raw_header <- names(fread(opt$counts, nrows = 0L))

base_cols <- c("Probe", "ROIlabel", "hemi", "moda", "outcome_fdr_mask")
if (!all(base_cols %in% raw_header)) {
  missing_base <- setdiff(base_cols, raw_header)
  stop("Missing required base columns in input table: ", paste(missing_base, collapse = ", "))
}

# Identify trait columns: accept either "<trait>__same" or "<trait>"
trait_cols <- character()
for (t in traits) {
  if (paste0(t, "__same") %in% raw_header) {
    trait_cols <- c(trait_cols, paste0(t, "__same"))
  } else if (t %in% raw_header) {
    trait_cols <- c(trait_cols, t)
  } else {
    stop("Input table is missing column for trait '", t, "'. Expected '", t, "' or '", paste0(t, "__same"), "'")
  }
}

cols_to_load <- unique(c(base_cols, trait_cols, intersect(c("all_same", "overlap_count"), raw_header)))
counts <- fread(opt$counts, select = cols_to_load, showProgress = FALSE)

if (nrow(counts) != 653994L || anyDuplicated(counts$Probe)) {
  stop("Input table must have exactly 653,994 unique vertex rows matching standard brainMapR probe inventory")
}
if (anyNA(counts[, ..cols_to_load])) stop("Input table contains NA values")

counts[, outcome_fdr_mask := as.logical(outcome_fdr_mask)]

# Standardize trait columns to "<trait>__same"
for (i in seq_along(traits)) {
  t <- traits[[i]]
  src_col <- trait_cols[[i]]
  dst_col <- paste0(t, "__same")
  if (src_col != dst_col) {
    setnames(counts, src_col, dst_col)
  }
}
standard_trait_cols <- paste0(traits, "__same")
trait_matrix <- as.matrix(counts[, ..standard_trait_cols])

# Derive total overlap count if not already provided
if ("overlap_count" %in% names(counts)) {
  counts[, all_same := as.integer(overlap_count)]
} else if (!"all_same" %in% names(counts)) {
  counts[, all_same := as.integer(rowSums(trait_matrix > 0L))]
}
counts[, all_same := as.integer(all_same)]

# Single trait identification for 1-trait locations
counts[, single_trait := NA_character_]
one <- which(counts$all_same == 1L)
if (length(one)) {
  counts$single_trait[one] <- traits[max.col(trait_matrix[one, , drop = FALSE], ties.method = "first")]
}

prepare_cortical <- function(modality_id, hemisphere_id) {
  atlas_file <- file.path(
    opt$atlas_dir,
    paste0("Atlas_allCoordinates_corticalVertices_", hemisphere_id, "_", modality_id, ".txt")
  )
  atlas <- fread(atlas_file, showProgress = FALSE)
  setnames(atlas, "vertexID", "Probe")
  x <- merge(counts[moda == modality_id & hemi == hemisphere_id], atlas, by = "Probe")
  if (!nrow(x)) stop("Empty cortical annotation")
  setorder(x, order)
  x
}
prepare_subcortical <- function(modality_id, hemisphere_id) {
  atlas <- fread(file.path(opt$atlas_dir, "Atlas_coordinates_ROI_subortical.txt"), showProgress = FALSE)
  setnames(atlas, "ProbeID", "Probe")
  x <- merge(counts[moda == modality_id & hemi == hemisphere_id], atlas, by = "Probe")
  if (!nrow(x)) stop("Empty subcortical annotation")
  x[, `:=`(X = -X, Y = -Y)]
  as.data.table(brainMapR::addFlatCoordinatesSubcortical(annotBwas = x, hemi = hemisphere_id))
}
prepared <- list()
for (modality_id in names(modality_labels)) {
  for (hemisphere_id in c("lh", "rh")) {
    key <- paste(modality_id, hemisphere_id, sep = "__")
    prepared[[key]] <- if (modality_id %in% cortical_modalities) {
      prepare_cortical(modality_id, hemisphere_id)
    } else {
      prepare_subcortical(modality_id, hemisphere_id)
    }
  }
}

# ==============================================================================
# 🗺️ 4. MAP SPECIFICATIONS (Trait Count & Trait Identity)
# ==============================================================================
map_specs <- list(
  count = list(
    label = "Concordant trait overlap count", kind = "numeric", value = "all_same"
  ),
  trait_identity = list(
    label = "Trait overlap identity", kind = "binary_identity", value = NA_character_
  )
)
for (trait_id in traits) {
  map_specs[[paste0("trait_", trait_id)]] <- list(
    label = paste0(trait_labels[[trait_id]], " overlap"),
    kind = "single_trait",
    value = paste0(trait_id, "__same"),
    trait = trait_id
  )
}

mix_with_white <- function(colour, fraction) {
  rgb <- col2rgb(colour)
  mixed <- round(255 * (1 - fraction) + rgb * fraction)
  rgb(mixed[1], mixed[2], mixed[3], maxColorValue = 255)
}
tint_colour <- function(colour, white_fraction = 0.30) {
  mix_with_white(colour, 1 - white_fraction)
}
numeric_palette <- function(maximum, base = NULL) {
  if (maximum < 1L) return(character())
  if (is.null(base)) {
    colorRampPalette(c(
      RColorBrewer::brewer.pal(9, "YlOrRd")[4],
      RColorBrewer::brewer.pal(9, "YlOrRd")[9]
    ))(maximum)
  } else {
    colorRampPalette(c(tint_colour(base), base))(maximum)
  }
}
domain_count_palette <- function(maximum, base) {
  if (maximum < 1L) return(character())
  base_rgb <- col2rgb(base)
  base_hsv <- grDevices::rgb2hsv(
    base_rgb[1L, 1L], base_rgb[2L, 1L], base_rgb[3L, 1L],
    maxColorValue = 255
  )
  if (opt$outcome_id == "ADvsHC") {
    # Keep the canonical hue while increasing separation between count levels.
    shades <- vapply(seq_len(maximum), function(trait_count) {
      if (trait_count == 1L) {
        return(grDevices::hsv(
          h = base_hsv["h", 1L],
          s = min(1, base_hsv["s", 1L] + 0.05),
          v = min(1, base_hsv["v", 1L] + 0.14)
        ))
      }
      darkness <- min(0.20 + 0.30 * (trait_count - 2L), 0.65)
      channel <- round(base_rgb[, 1L] * (1 - darkness))
      grDevices::rgb(channel[1L], channel[2L], channel[3L], maxColorValue = 255)
    }, character(1L))
    return(shades)
  }
  strength <- if (maximum == 1L) 1 else seq(1, 0, length.out = maximum)
  shades <- grDevices::hsv(
    h = rep(base_hsv["h", 1L], maximum),
    s = pmin(1, base_hsv["s", 1L] + 0.05 * strength),
    v = pmin(1, base_hsv["v", 1L] + 0.14 * strength)
  )
  if (maximum > 1L) shades[[maximum]] <- base
  shades
}
breadth_base_palette <- c("#FB6A4A", "#EF3B2C", "#CB181D", "#99000D", "#67000D")
get_breadth_palette <- function(k) {
  if (k <= 0L) return(character())
  if (k <= length(breadth_base_palette)) {
    return(breadth_base_palette[seq_len(k)])
  }
  colorRampPalette(breadth_base_palette)(k)
}

spec_palette <- function(spec) {
  if (spec$kind == "numeric") {
    max_count <- max(1L, max(counts[[spec$value]]))
    return(numeric_palette(max_count))
  }
  if (spec$kind == "single_trait") {
    return(trait_colours[[spec$trait]])
  }
  character()
}

scene_anchors <- function(coordinates) {
  ranges <- apply(coordinates, 2L, range, na.rm = TRUE)
  centre <- colMeans(ranges)
  half <- pmax((ranges[2L, ] - ranges[1L, ]) / 2, 1e-6)
  signs <- as.matrix(expand.grid(c(-1, 1), c(-1, 1), c(-1, 1)))
  sweep(sweep(signs, 2L, half * 1.30, "*"), 2L, centre, "+")
}
draw_layer <- function(coordinates, selected, colours, radius) {
  if (!any(selected)) return(invisible(NULL))
  rgl::spheres3d(
    coordinates[selected, , drop = FALSE], col = colours[selected],
    radius = radius, alpha = point_alpha
  )
}
colours_for_spec <- function(x, spec) {
  # 1. Base grey: all 654k vertices start as neutral background grey
  colours <- rep(background_colour, nrow(x))
  
  # 2. Outcome FDR: vertices passing outcome significance turn beige
  in_mask <- x$outcome_fdr_mask
  colours[in_mask] <- outcome_only_colour
  palette <- spec_palette(spec)
  
  if (spec$kind == "numeric") {
    # Count gradient: 1, 2, 3, 4+ overlapping traits
    values <- x[[spec$value]]
    hit <- in_mask & values > 0L
    if (any(hit)) {
      colours[hit] <- palette[pmin(values[hit], length(palette))]
    }
  } else if (spec$kind == "single_trait") {
    # Single trait map: highlight vertices where this specific trait is active
    hit <- in_mask & x[[spec$value]] > 0L
    if (any(hit)) {
      colours[hit] <- trait_colours[[spec$trait]]
    }
  } else if (spec$kind == "binary_identity") {
    # Trait identity map:
    # A) Single trait: receives that specific trait's unique color
    for (trait_id in traits) {
      hit <- in_mask & x$all_same == 1L & x$single_trait == trait_id
      colours[hit] <- trait_colours[[trait_id]]
    }
    # B) Multi-trait overlap: 2 or more traits overlapping simultaneously receive dark grey
    colours[in_mask & x$all_same > 1L] <- multi_trait_colour
  } else {
    stop("Unsupported map kind: ", spec$kind)
  }
  colours
}
render_component <- function(x, coordinates, spec, output_file) {
  colours <- colours_for_spec(x, spec)
  outside <- !x$outcome_fdr_mask
  inside <- x$outcome_fdr_mask
  rgl::open3d()
  rgl::bg3d(color = "white")
  rgl::par3d(windowRect = c(0, 0, 1200, 1200), zoom = 0.80)
  rgl::spheres3d(scene_anchors(coordinates), col = "#FFFFFF", alpha = 0, radius = 0.01)
  draw_layer(coordinates, outside, colours, background_radius)
  if (opt$outcome_id == "ADvsHC" && spec$kind == "domain_count") {
    positive <- inside & x[[spec$value]] > 0L
    draw_layer(coordinates, inside & !positive, colours, significant_radius)
    positive_radius <- if (spec$domain == "infection") 4.00 else 1.85
    draw_layer(coordinates, positive, colours, positive_radius)
  } else {
    draw_layer(coordinates, inside, colours, significant_radius)
  }
  rgl::par3d(zoom = 0.80)
  rgl::rgl.snapshot(output_file)
  rgl::close3d()
}
component_path <- function(map_id, modality_id, hemisphere_id, view) {
  component_dir <- file.path(opt$output_dir, "components", map_id)
  dir.create(component_dir, recursive = TRUE, showWarnings = FALSE)
  file.path(component_dir, paste(modality_id, hemisphere_id, view, "png", sep = "_"))
}
coordinates_for <- function(x, modality_id, view) {
  if (modality_id %in% cortical_modalities) {
    cols <- paste(c("Y", "X", "Z"), style, sep = "_")
    xyz <- as.matrix(x[, ..cols])
    if (view == "outside") xyz[, 1:2] <- -xyz[, 1:2]
  } else {
    cols <- if (view == "inside") c("Zf", "Xf", "Yf") else c("Zf2", "Xf2", "Yf2")
    xyz <- as.matrix(x[, ..cols])
  }
  if (any(!is.finite(xyz))) stop("Non-finite plotting coordinates")
  xyz
}

map_spec_filter <- Sys.getenv("MAP_SPEC_FILTER", "")
component_only <- tolower(Sys.getenv("COMPONENT_ONLY", "false")) == "true"
render_map_ids <- names(map_specs)
if (nzchar(map_spec_filter)) {
  requested <- strsplit(map_spec_filter, ",", fixed = TRUE)[[1L]]
  if (any(!requested %in% names(map_specs))) {
    stop("Unknown MAP_SPEC_FILTER value: ", map_spec_filter)
  }
  render_map_ids <- requested
}

reuse_components <- tolower(Sys.getenv("REUSE_COMPONENTS", "false")) == "true"
skip_existing_components <-
  tolower(Sys.getenv("SKIP_EXISTING_COMPONENTS", "false")) == "true"
if (reuse_components) {
  expected_components <- unlist(lapply(render_map_ids, function(map_id) {
    unlist(lapply(names(modality_labels), function(modality_id) {
      unlist(lapply(c("lh", "rh"), function(hemisphere_id) {
        vapply(c("inside", "outside"), function(view) {
          component_path(map_id, modality_id, hemisphere_id, view)
        }, character(1L))
      }))
    }))
  }))
  file_sizes <- file.info(expected_components)$size
  expected_count <- length(render_map_ids) * length(modality_labels) * 2L * 2L
  if (length(expected_components) != expected_count || anyNA(file_sizes) || any(file_sizes <= 0L)) {
    stop("REUSE_COMPONENTS requested, but the expected component images are incomplete")
  }
  message("Reusing ", expected_count, " validated component images for ", opt$outcome_id)
} else {
  for (map_id in render_map_ids) {
    spec <- map_specs[[map_id]]
    message("Rendering ", opt$outcome_id, " / ", map_id)
    for (modality_id in names(modality_labels)) {
      for (hemisphere_id in c("lh", "rh")) {
        x <- prepared[[paste(modality_id, hemisphere_id, sep = "__")]]
        for (view in c("inside", "outside")) {
          output_file <- component_path(
            map_id, modality_id, hemisphere_id, view
          )
          if (skip_existing_components && file.exists(output_file) &&
              file.info(output_file)$size > 100000L) {
            message("Reusing validated rendered view: ", output_file)
            next
          }
          render_component(
            x, coordinates_for(x, modality_id, view), spec, output_file
          )
        }
      }
    }
  }
}

if (component_only) {
  marker_id <- gsub("[^A-Za-z0-9_]+", "_", paste(render_map_ids, collapse = "__"))
  writeLines("COMPLETE", file.path(opt$output_dir, paste0("COMPONENTS_", marker_id, ".COMPLETE")))
  message("Completed component-only rendering for ", opt$outcome_id, ": ", paste(render_map_ids, collapse = ", "))
  quit(save = "no", status = 0L)
}

image_grob <- function(path) {
  if (!file.exists(path)) stop("Missing component: ", path)
  image <- readPNG(path)
  nonwhite <- image[, , 1L] < 0.985 |
    image[, , 2L] < 0.985 |
    image[, , 3L] < 0.985
  occupied_rows <- which(rowSums(nonwhite) > 0L)
  occupied_columns <- which(colSums(nonwhite) > 0L)
  if (!length(occupied_rows) || !length(occupied_columns)) {
    stop("Component image has no visible content: ", path)
  }
  cropped <- image[
    min(occupied_rows):max(occupied_rows),
    min(occupied_columns):max(occupied_columns),
    , drop = FALSE
  ]
  content_height <- dim(cropped)[1L]
  content_width <- dim(cropped)[2L]
  channel_count <- dim(cropped)[3L]
  canvas_side <- ceiling(max(content_height, content_width))
  canvas <- array(1, dim = c(canvas_side, canvas_side, channel_count))
  row_start <- floor((canvas_side - content_height) / 2) + 1L
  column_start <- floor((canvas_side - content_width) / 2) + 1L
  canvas[
    row_start:(row_start + content_height - 1L),
    column_start:(column_start + content_width - 1L),
  ] <- cropped
  rasterGrob(
    canvas, width = unit(0.98, "snpc"), height = unit(0.98, "snpc"),
    interpolate = TRUE
  )
}

column_specs <- data.table(
  hemisphere = c("lh", "rh", "lh", "rh", "lh", "rh", "lh", "rh"),
  modality = c(
    "thickness", "thickness", "area", "area",
    "thick", "thick", "LogJacs", "LogJacs"
  ),
  label = c(
    "Left\nCT", "Right\nCT", "Left\nSA", "Right\nSA",
    "Left\nsubcortical\nthickness", "Right\nsubcortical\nthickness",
    "Left\nsubcortical\nLog-Jacobian", "Right\nsubcortical\nLog-Jacobian"
  )
)

spaced_row <- function(left_grob, panel_grobs) {
  grobs <- list(left_grob)
  widths <- 0.60
  for (i in seq_along(panel_grobs)) {
    grobs <- c(grobs, list(panel_grobs[[i]]))
    widths <- c(widths, 1)
    if (i < length(panel_grobs)) {
      gap <- if (i %% 2L == 0L) 0.075 else 0.018
      grobs <- c(grobs, list(nullGrob()))
      widths <- c(widths, gap)
    }
  }
  arrangeGrob(
    grobs = grobs, nrow = 1L, widths = widths,
    padding = unit(1, "mm")
  )
}

make_brain_grid <- function(map_id) {
  header <- spaced_row(
    nullGrob(),
    lapply(column_specs$label, function(label) {
      textGrob(
        label, y = unit(0.04, "npc"), just = "bottom",
        gp = gpar(
          fontsize = FIG_FONT$axis_text, fontface = "bold", col = FIG_TEXT_COLOUR,
          fontfamily = FIG_FONT_FAMILY, lineheight = 0.9
        )
      )
    })
  )
  make_view_row <- function(component_view, label) {
    panels <- lapply(seq_len(nrow(column_specs)), function(i) {
      image_grob(component_path(
        map_id, column_specs$modality[[i]],
        column_specs$hemisphere[[i]], component_view
      ))
    })
    spaced_row(
      textGrob(
        label, rot = 90,
        gp = gpar(
          fontsize = FIG_FONT$axis_text, fontface = "bold", col = FIG_TEXT_COLOUR,
          fontfamily = FIG_FONT_FAMILY
        )
      ),
      panels
    )
  }
  arrangeGrob(
    header,
    make_view_row("outside", "Outside"),
    make_view_row("inside", "Inside"),
    ncol = 1L, heights = c(0.25, 1, 1),
    padding = unit(1, "mm")
  )
}

legend_row_widths <- unit.c(unit(0.50, "null"), unit(0.70, "null"))
legend_label_x <- 0.96

make_count_bar <- function(
  label, values, colours, continuous = FALSE,
  breaks = values, tick_labels = as.character(breaks), slots = NULL
) {
  d <- data.frame(x = values, y = 1, colour = colours)
  plot <- ggplot(d, aes(x, y, fill = colour))
  plot <- if (continuous) {
    plot + geom_raster()
  } else {
    plot + geom_tile(
      width = 1, height = 1, colour = "white",
      linewidth = FIG_LINE$key_border
    )
  }
  plot <- plot +
    scale_fill_identity() +
    scale_x_continuous(
      breaks = breaks, labels = tick_labels, expand = c(0, 0),
      limits = if (is.null(slots)) range(values) + c(-0.5, 0.5) else c(0.5, slots + 0.5)
    ) +
    theme_void(base_family = FIG_FONT_FAMILY, base_size = FIG_FONT$legend_text) +
    theme(
      aspect.ratio = 0.075,
      axis.text.x = element_text(
        size = FIG_FONT$legend_text, face = "bold", colour = FIG_TEXT_COLOUR
      ),
      plot.margin = margin(1, 8, 1, 4), legend.position = "none"
    )
  arrangeGrob(
    textGrob(
      gsub(" / ", "\n", label, fixed = TRUE),
      x = unit(legend_label_x, "npc"), just = "right",
      gp = gpar(
        fontsize = FIG_FONT$legend_text, fontface = "bold",
        col = FIG_TEXT_COLOUR, fontfamily = FIG_FONT_FAMILY,
        lineheight = 0.90
      )
    ),
    ggplotGrob(plot), ncol = 2L, widths = legend_row_widths,
    padding = unit(0.15, "mm")
  )
}

legend_heading <- function(label) {
  arrangeGrob(
    textGrob(
      label, x = unit(legend_label_x, "npc"), just = "right",
      gp = gpar(
        fontsize = FIG_FONT$legend_title, fontface = "bold",
        fontfamily = FIG_FONT_FAMILY, col = FIG_TEXT_COLOUR
      )
    ),
    nullGrob(),
    ncol = 2L, widths = legend_row_widths,
    padding = unit(0.15, "mm")
  )
}

# A one-cell bar gives every discrete key the same geometry as count = 1.
legend_cell_bar <- function(label, colour, maximum, fontface = "plain") {
  plot <- ggplot(
    data.frame(x = 1, y = 1, colour = colour),
    aes(x, y, fill = colour)
  ) +
    geom_tile(
      width = 1, height = 1, colour = "white",
      linewidth = FIG_LINE$key_border
    ) +
    scale_fill_identity() +
    scale_x_continuous(
      breaks = NULL, expand = c(0, 0), limits = c(0.5, maximum + 0.5)
    ) +
    theme_void(base_family = FIG_FONT_FAMILY, base_size = FIG_FONT$legend_text) +
    theme(
      aspect.ratio = 0.075,
      plot.margin = margin(1, 8, 1, 4), legend.position = "none"
    )
  arrangeGrob(
    textGrob(
      label, x = unit(legend_label_x, "npc"), just = "right",
      gp = gpar(
        fontsize = FIG_FONT$legend_text, fontface = fontface,
        fontfamily = FIG_FONT_FAMILY, col = FIG_TEXT_COLOUR
      )
    ),
    ggplotGrob(plot), ncol = 2L, widths = legend_row_widths,
    padding = unit(0.10, "mm")
  )
}

context_bars <- function(maximum, domain_specific = FALSE) {
  overlap_label <- if (domain_specific) {
    "Outcome FDR,\nno overlap in this domain"
  } else {
    "Outcome FDR, no overlap"
  }
  list(
    legend_cell_bar(
      overlap_label, outcome_only_colour, maximum
    ),
    legend_cell_bar(
      "Outside outcome FDR", background_colour, maximum
    )
  )
}

single_count_legend <- function(spec) {
  maximum <- max(counts[[spec$value]])
  if (maximum < 1L) maximum <- 1L
  values <- seq_len(maximum)
  colours <- spec_palette(spec)
  label <- spec$label
  context <- context_bars(maximum)
  arrangeGrob(
    nullGrob(), legend_heading("Concordant traits"),
    make_count_bar(label, values, colours, slots = maximum),
    context[[1L]], context[[2L]], nullGrob(),
    ncol = 1L, heights = c(0.65, 0.16, 0.22, 0.12, 0.12, 0.65),
    padding = unit(0.10, "mm")
  )
}

trait_identity_legend <- function() {
  maximum <- max(1L, max(counts$all_same))
  keys <- lapply(traits, function(trait_id) {
    legend_cell_bar(
      trait_labels[[trait_id]], trait_colours[[trait_id]],
      maximum, fontface = "bold"
    )
  })
  keys[[length(keys) + 1L]] <- legend_cell_bar(
    "Multiple traits", multi_trait_colour, maximum, fontface = "bold"
  )
  context <- context_bars(maximum)
  n_keys <- length(keys)
  key_height <- min(0.12, 0.90 / (n_keys + 2L))
  arrangeGrob(
    grobs = c(
      list(nullGrob(), legend_heading("Trait overlap")),
      keys, context, list(nullGrob())
    ),
    ncol = 1L,
    heights = c(
      0.48, 0.15, rep(key_height, n_keys),
      0.115, 0.115, 0.48
    ),
    padding = unit(0.10, "mm")
  )
}

single_trait_legend <- function(spec) {
  maximum <- max(1L, max(counts$all_same))
  key <- list(legend_cell_bar(
    trait_labels[[spec$trait]], trait_colours[[spec$trait]],
    maximum, fontface = "bold"
  ))
  context <- context_bars(maximum)
  arrangeGrob(
    grobs = c(
      list(nullGrob(), legend_heading("Trait presence")),
      key, context, list(nullGrob())
    ),
    ncol = 1L,
    heights = c(0.55, 0.16, 0.14, 0.12, 0.12, 0.55),
    padding = unit(0.10, "mm")
  )
}

legend_for <- function(spec) {
  if (spec$kind == "numeric") return(single_count_legend(spec))
  if (spec$kind == "binary_identity") return(trait_identity_legend())
  if (spec$kind == "single_trait") return(single_trait_legend(spec))
  stop("Unsupported legend for map kind: ", spec$kind)
}

make_panel <- function(map_id, panel_label = NULL) {
  spec <- map_specs[[map_id]]
  body <- arrangeGrob(
    make_brain_grid(map_id), nullGrob(), legend_for(spec),
    nrow = 1L,
    widths = unit.c(unit(3.55, "null"), unit(3.0, "mm"), unit(0.95, "null")),
    padding = unit(1, "mm")
  )
  panel_title <- paste(c(panel_label, spec$label), collapse = "  ")
  arrangeGrob(
    textGrob(
      panel_title, x = unit(0.01, "npc"), just = "left",
      gp = gpar(
        fontsize = FIG_FONT$panel_title, fontface = "bold", fontfamily = FIG_FONT_FAMILY,
        col = FIG_TEXT_COLOUR
      )
    ),
    body, ncol = 1L, heights = c(0.08, 1)
  )
}

save_grob <- function(grob, stem, width = 14.4, height = 3.8) {
  png_path <- file.path(opt$output_dir, paste0(stem, ".png"))
  pdf_path <- file.path(opt$output_dir, paste0(stem, ".pdf"))
  png(
    png_path, width = round(width * 500), height = round(height * 500),
    res = 500, bg = "white"
  )
  grid.newpage()
  grid.draw(grob)
  dev.off()
  cairo_pdf(
    pdf_path, width = width, height = height,
    family = FIG_FONT_FAMILY, onefile = FALSE
  )
  grid.newpage()
  grid.draw(grob)
  dev.off()
  if (file.info(png_path)$size <= 0 || file.info(pdf_path)$size <= 0) {
    stop("Empty output: ", png_path)
  }
}

# ==============================================================================
# 🖼️ 5. COMPOSE & SAVE PUBLICATION FIGURES
# ==============================================================================
safe_id <- gsub("[^A-Za-z0-9]+", "", opt$outcome_id)

# 1. Main composite publication figure (Panel a: Count, Panel b: Trait Identity)
count_panel <- make_panel("count", "a")
identity_panel <- make_panel("trait_identity", "b")
main_figure <- arrangeGrob(
  textGrob(
    display_outcome, x = unit(0.005, "npc"), just = "left",
    gp = gpar(fontsize = FIG_FONT$figure_title, fontface = "bold", fontfamily = FIG_FONT_FAMILY)
  ),
  count_panel, identity_panel,
  ncol = 1L, heights = c(0.05, 1, 1), padding = unit(1, "mm")
)

save_grob(
  main_figure, paste0("Main_", safe_id, "_count_and_trait"),
  width = 14.4, height = 7.7
)

# 2. Individual sub-panels (for separate slides or manuscripts)
save_grob(
  count_panel, paste0("Panel_", safe_id, "_trait_count"),
  width = 14.4, height = 3.8
)
save_grob(
  identity_panel, paste0("Panel_", safe_id, "_trait_identity"),
  width = 14.4, height = 3.8
)

# 3. Individual trait supplementary figures
for (trait_id in traits) {
  t_panel <- make_panel(paste0("trait_", trait_id), trait_labels[[trait_id]])
  save_grob(
    t_panel, paste0("Supplement_", safe_id, "_", trait_id, "_overlap"),
    width = 14.4, height = 3.8
  )
}

# 4. Detailed captions and parameter summaries
caption <- c(
  paste0(display_outcome, " spatial convergence across ", length(traits), " selected traits."),
  "Panel a gives the integer count of concordant overlapping traits at each outcome-FDR-significant location.",
  "Panel b displays trait identity: single-trait locations receive that trait's distinct color; dark grey indicates locations with 2 or more overlapping traits.",
  "Grey points lie outside the outcome FDR q<0.05 mask; beige points pass outcome FDR but have zero overlapping traits.",
  "Cortical geometry uses inflated_pre; subcortical views use the brainMapR flat projection across 7 bilateral nuclei.",
  sprintf("Background point radius is %.2f; outcome/trait points use radius %.2f.", background_radius, significant_radius),
  "Horizontal layout orders left/right CT, left/right SA, left/right subcortical thickness, and left/right subcortical Log-Jacobian across 8 columns."
)
writeLines(caption, file.path(opt$output_dir, paste0("Main_", safe_id, "_caption.txt")))

fwrite(data.table(
  setting = c("surface", "background_radius", "significant_radius", "alpha",
              "background_colour", "outcome_only_colour", "multi_trait_colour",
              "max_overlapping_traits", "trait_count"),
  value = c(style, background_radius, significant_radius, point_alpha,
            background_colour, outcome_only_colour, multi_trait_colour,
            max(counts$all_same), length(traits))
), file.path(opt$output_dir, "RENDERING_PARAMETERS.tsv"), sep = "\t")

writeLines(capture.output(sessionInfo()), file.path(opt$output_dir, "sessionInfo.txt"))
writeLines("COMPLETE", file.path(opt$output_dir, "COMPLETE.txt"))
message("Successfully completed manuscript spatial rendering for ", opt$outcome_id)
