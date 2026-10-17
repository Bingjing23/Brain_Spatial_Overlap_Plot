# Shared typography, line weights, and publication palettes
# Standardized for cortical/subcortical spatial figures

FIG_FONT_FAMILY <- "Helvetica"
FIG_FONT <- list(
  figure_title = 10.0,
  panel_title  = 9.0,
  legend_title = 8.0,
  axis_title   = 8.0,
  strip        = 7.0,
  legend_text  = 7.0,
  axis_text    = 7.0,
  annotation   = 6.5
)

FIG_LINE <- list(
  cell_border = 0.22,
  key_border  = 0.35,
  reference   = 0.35
)

FIG_TEXT_COLOUR <- "#202124"
FIG_SECONDARY_TEXT_COLOUR <- "#4A5359"
FIG_KEY_BORDER_COLOUR <- "#6B7378"

# Background point styling
DEFAULT_BG_COLOUR <- "#D0D0D0"
DEFAULT_BG_RADIUS <- 0.95
DEFAULT_FG_RADIUS <- 1.55
DEFAULT_FG_ALPHA  <- 1.0

# Diverging color palette for continuous statistics (e.g. z-scores, rGM, beta)
FIG_DIVERGING_COLOURS <- c(
  low = "#2C7BB6",
  mid = "#F7F7F7",
  high = "#B2182B"
)

publication_diverging_palette <- function(n = 256L) {
  grDevices::colorRampPalette(
    unname(FIG_DIVERGING_COLOURS),
    space = "Lab"
  )(n)
}

# Discrete count color palette for integer overlap (e.g., 1, 2, 3... traits)
publication_count_palette <- function(max_count) {
  if (max_count <= 0L) return(character(0))
  if (max_count == 1L) return("#E65100") # Saturated amber/orange
  if (max_count <= 9L) {
    # Warm YlOrRd ramp
    colours <- RColorBrewer::brewer.pal(max(3L, min(9L, max_count + 1L)), "YlOrRd")
    return(grDevices::colorRampPalette(colours[3:length(colours)])(max_count))
  }
  grDevices::colorRampPalette(c("#FED976", "#FD8D3C", "#E31A1C", "#800026"))(max_count)
}
