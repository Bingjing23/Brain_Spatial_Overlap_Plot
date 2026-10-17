# Core rendering engine for cortical and subcortical vertex-level maps
# Compatible with brainMapR atlas definitions and UK Biobank BWAS coordinate standards

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(png)
  library(rgl)
})

# Locate brainMapR atlas files
find_atlas_dir <- function(atlas_dir = NULL) {
  candidates <- c()
  if (!is.null(atlas_dir) && nzchar(atlas_dir)) {
    candidates <- c(candidates, atlas_dir)
  }
  # Check R package system directory if installed
  if (requireNamespace("brainMapR", quietly = TRUE)) {
    pkg_atlas <- system.file("extdata", "atlas", package = "brainMapR")
    if (nzchar(pkg_atlas)) candidates <- c(candidates, pkg_atlas)
  }
  # Common cluster / local paths
  candidates <- c(
    candidates,
    "/mnt/lustre/working/lab_puyag/bingjinZ/Side_Project_AD/software/Rlibs/brainMapR-076069ce/brainMapR/extdata/atlas",
    "/working/lab_puyag/bingjinZ/Side_Project_AD/software/Rlibs/brainMapR-076069ce/brainMapR/extdata/atlas"
  )
  for (path in candidates) {
    if (dir.exists(path)) {
      sub_test <- file.path(path, "Atlas_coordinates_ROI_subortical.txt")
      if (file.exists(sub_test)) return(normalizePath(path))
    }
  }
  stop(
    "Could not locate brainMapR atlas directory.\n",
    "Please provide the path via the '--atlas-dir' argument or ensure package 'brainMapR' is installed."
  )
}

# 3D Scene Anchoring (prevents rgl perspective warping across different views)
scene_anchors <- function(coordinates) {
  coordinate_ranges <- apply(coordinates, 2L, range, na.rm = TRUE)
  scene_center <- colMeans(coordinate_ranges)
  half_extent <- pmax((coordinate_ranges[2L, ] - coordinate_ranges[1L, ]) / 2, 1e-6)
  anchor_signs <- as.matrix(expand.grid(c(-1, 1), c(-1, 1), c(-1, 1)))
  anchors <- sweep(anchor_signs, 2L, half_extent * 1.30, "*")
  sweep(anchors, 2L, scene_center, "+")
}

# Layer drawing with spheres3d
draw_spheres_layer <- function(coordinates, mask, colours, radius, alpha = 1.0) {
  if (!any(mask)) return(invisible(NULL))
  rgl::spheres3d(
    coordinates[mask, , drop = FALSE],
    col = colours[mask],
    radius = radius,
    alpha = alpha
  )
}

# Render a single brain view component to PNG
render_single_view <- function(
  coordinates,
  mask_fg,
  colours_all,
  output_file,
  bg_radius = DEFAULT_BG_RADIUS,
  fg_radius = DEFAULT_FG_RADIUS,
  fg_alpha  = DEFAULT_FG_ALPHA,
  bg_colour = DEFAULT_BG_COLOUR
) {
  mask_bg <- !mask_fg
  rgl::open3d()
  rgl::bg3d(color = "white")
  rgl::par3d(windowRect = c(0, 0, 1200, 1200), zoom = 0.80)
  
  # Invisible anchors to enforce consistent camera bounding box
  rgl::spheres3d(scene_anchors(coordinates), col = "#FFFFFF", alpha = 0, radius = 0.01)
  
  # Draw background non-significant vertices
  draw_spheres_layer(coordinates, mask_bg, rep(bg_colour, length(mask_bg)), bg_radius, alpha = 0.85)
  
  # Draw highlighted foreground vertices
  draw_spheres_layer(coordinates, mask_fg, colours_all, fg_radius, alpha = fg_alpha)
  
  rgl::par3d(zoom = 0.80)
  rgl::rgl.snapshot(output_file)
  rgl::close3d()
}

# Crop whitespace and convert to rasterGrob
image_file_to_grob <- function(path) {
  if (!file.exists(path)) stop("Missing rendered component: ", path)
  img <- readPNG(path)
  nonwhite <- img[, , 1L] < 0.985 | img[, , 2L] < 0.985 | img[, , 3L] < 0.985
  rows <- which(rowSums(nonwhite) > 0L)
  cols <- which(colSums(nonwhite) > 0L)
  if (!length(rows) || !length(cols)) stop("Component image contains only white pixels: ", path)
  
  cropped <- img[min(rows):max(rows), min(cols):max(cols), , drop = FALSE]
  h <- dim(cropped)[1L]
  w <- dim(cropped)[2L]
  channels <- dim(cropped)[3L]
  side <- ceiling(max(h, w))
  
  canvas <- array(1, dim = c(side, side, channels))
  r_start <- floor((side - h) / 2) + 1L
  c_start <- floor((side - w) / 2) + 1L
  canvas[r_start:(r_start + h - 1L), c_start:(c_start + w - 1L), ] <- cropped
  
  rasterGrob(canvas, width = unit(0.98, "snpc"), height = unit(0.98, "snpc"), interpolate = TRUE)
}

# Create a clean horizontal count/scale legend
make_count_legend_grob <- function(max_count, colours) {
  df <- data.frame(count = seq_len(max_count), y = 1, col = colours[seq_len(max_count)])
  p <- ggplot(df, aes(x = count, y = y, fill = col)) +
    geom_tile(width = 1, height = 1, colour = "white", linewidth = FIG_LINE$key_border) +
    scale_fill_identity() +
    scale_x_continuous(breaks = seq_len(max_count), expand = c(0, 0)) +
    labs(x = "Number of concordant / overlapping traits") +
    theme_void(base_family = FIG_FONT_FAMILY, base_size = FIG_FONT$legend_text) +
    theme(
      axis.title.x = element_text(size = FIG_FONT$legend_title, margin = margin(t = 4, b = 2)),
      axis.text.x  = element_text(size = FIG_FONT$legend_text, margin = margin(t = 2)),
      axis.ticks.x = element_blank(),
      plot.margin  = margin(2, 5, 2, 5)
    )
  ggplotGrob(p)
}

make_continuous_legend_grob <- function(limit, title = "Value") {
  vals <- seq(-limit, limit, length.out = 100)
  df <- data.frame(val = vals, y = 1, col = publication_diverging_palette(100))
  p <- ggplot(df, aes(x = val, y = y, fill = col)) +
    geom_raster() +
    scale_fill_identity() +
    scale_x_continuous(
      breaks = c(-limit, 0, limit),
      labels = format(c(-limit, 0, limit), digits = 2),
      expand = c(0, 0)
    ) +
    labs(x = title) +
    theme_void(base_family = FIG_FONT_FAMILY, base_size = FIG_FONT$legend_text) +
    theme(
      axis.title.x = element_text(size = FIG_FONT$legend_title, margin = margin(t = 4, b = 2)),
      axis.text.x  = element_text(size = FIG_FONT$legend_text, margin = margin(t = 2)),
      plot.margin  = margin(2, 5, 2, 5)
    )
  ggplotGrob(p)
}
