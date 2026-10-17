# Shared typography and line-weight contract for manuscript figures.
# Values are in typographic points unless a plotting API documents otherwise.

FIG_FONT_FAMILY <- "Helvetica"
FIG_FONT <- list(
  figure_title = 10.0,
  panel_title = 9.0,
  legend_title = 8.0,
  axis_title = 8.0,
  strip = 7.0,
  legend_text = 7.0,
  axis_text = 7.0,
  annotation = 6.5,
  dense_axis = 5.6
)
FIG_LINE <- list(
  cell_border = 0.22,
  key_border = 0.35,
  reference = 0.35
)
FIG_TEXT_COLOUR <- "#202124"
FIG_SECONDARY_TEXT_COLOUR <- "#4A5359"
FIG_KEY_BORDER_COLOUR <- "#6B7378"

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

draw_publication_diverging_key_base <- function(
  fig, name, limits = c(-1, 1), breaks = c(-1, 0, 1)
) {
  if (length(fig) != 4L || length(limits) != 2L || diff(limits) <= 0) {
    stop("Invalid base colourbar geometry or limits")
  }
  if (any(breaks < limits[1L] | breaks > limits[2L])) {
    stop("Base colourbar breaks must lie within limits")
  }

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(
    fig = fig, new = TRUE, mar = c(0, 0, 0, 0), xpd = NA,
    family = FIG_FONT_FAMILY
  )
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")

  bar_left <- 0.05
  bar_right <- 0.30
  bar_bottom <- 0.18
  bar_top <- 0.73
  colours <- publication_diverging_palette(256L)
  edges <- seq(bar_bottom, bar_top, length.out = length(colours) + 1L)
  graphics::rect(
    xleft = bar_left,
    ybottom = edges[-length(edges)],
    xright = bar_right,
    ytop = edges[-1L],
    col = colours,
    border = NA
  )

  label_y <- bar_bottom +
    (breaks - limits[1L]) / diff(limits) * (bar_top - bar_bottom)
  graphics::text(
    x = 0.40, y = label_y,
    labels = format(breaks, trim = TRUE, scientific = FALSE),
    adj = c(0, 0.5), cex = FIG_FONT$legend_text / FIG_FONT$figure_title,
    col = FIG_TEXT_COLOUR
  )
  graphics::text(
    x = bar_left, y = 0.82, labels = name,
    adj = c(0, 0.5), font = 2,
    cex = FIG_FONT$legend_title / FIG_FONT$figure_title,
    col = FIG_TEXT_COLOUR
  )
}

publication_diverging_fill <- function(
  name, limits, midpoint = 0, bar_width_mm = 4.2, bar_height_mm = 36
) {
  ggplot2::scale_fill_gradient2(
    low = unname(FIG_DIVERGING_COLOURS[["low"]]),
    mid = unname(FIG_DIVERGING_COLOURS[["mid"]]),
    high = unname(FIG_DIVERGING_COLOURS[["high"]]),
    midpoint = midpoint,
    limits = limits,
    name = name,
    guide = ggplot2::guide_colourbar(
      direction = "vertical",
      barwidth = grid::unit(bar_width_mm, "mm"),
      barheight = grid::unit(bar_height_mm, "mm"),
      theme = ggplot2::theme(
        legend.ticks = ggplot2::element_blank(),
        legend.axis.line = ggplot2::element_blank()
      )
    )
  )
}
