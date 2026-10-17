#!/usr/bin/env Rscript

# Generate a synthetic demonstration dataset for testing the brain plotting suite.
# Can read real Probe IDs from the atlas if available, or generate synthetic probes.

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
atlas_dir <- if (length(args) > 0) args[1] else NULL
output_file <- "example_overlap_data.tsv.gz"

cat("Generating mock vertex overlap dataset...\n")

# Try to find atlas directory
if (is.null(atlas_dir) || !dir.exists(atlas_dir)) {
  if (requireNamespace("brainMapR", quietly = TRUE)) {
    atlas_dir <- system.file("extdata", "atlas", package = "brainMapR")
  }
}

if (!is.null(atlas_dir) && dir.exists(atlas_dir)) {
  cat("Reading authentic Probe IDs from atlas at:", atlas_dir, "\n")
  dt_list <- list()
  for (hemi in c("lh", "rh")) {
    for (moda in c("thickness", "area")) {
      f <- file.path(atlas_dir, paste0("Atlas_allCoordinates_corticalVertices_", hemi, "_", moda, ".txt"))
      if (file.exists(f)) {
        x <- fread(f, select = c("vertexID", "ROIlabel"))
        setnames(x, "vertexID", "Probe")
        x[, `:=`(hemi = hemi, moda = moda)]
        dt_list[[paste(hemi, moda)]] <- x
      }
    }
  }
  sub_f <- file.path(atlas_dir, "Atlas_coordinates_ROI_subortical.txt")
  if (file.exists(sub_f)) {
    sub_raw <- fread(sub_f, select = c("ProbeID", "ROIlabel", "hemi", "moda"))
    setnames(sub_raw, "ProbeID", "Probe")
    dt_list[["subcortical"]] <- sub_raw
  }
  dt <- rbindlist(dt_list, use.names = TRUE)
} else {
  cat("Atlas not found; generating synthetic UKB-style probe inventory...\n")
  dt_list <- list()
  for (hemi in c("lh", "rh")) {
    # Cortical CT & SA
    dt_list[[paste0(hemi, "_ct")]] <- data.table(
      Probe = paste0("v_", hemi, "_ct_", seq_len(32768)),
      hemi = hemi, moda = "thickness", ROIlabel = "Cortex"
    )
    dt_list[[paste0(hemi, "_sa")]] <- data.table(
      Probe = paste0("v_", hemi, "_sa_", seq_len(32768)),
      hemi = hemi, moda = "area", ROIlabel = "Cortex"
    )
    # Subcortical
    dt_list[[paste0(hemi, "_thick")]] <- data.table(
      Probe = paste0("sub_", hemi, "_th_", seq_len(10000)),
      hemi = hemi, moda = "thick", ROIlabel = "Subcortex"
    )
    dt_list[[paste0(hemi, "_logj")]] <- data.table(
      Probe = paste0("sub_", hemi, "_lj_", seq_len(10000)),
      hemi = hemi, moda = "LogJacs", ROIlabel = "Subcortex"
    )
  }
  dt <- rbindlist(dt_list)
}

set.seed(42)
n <- nrow(dt)

# Simulate overlap count (0 to 4 traits)
prob_overlap <- 0.05 # 5% vertices show overlap
dt[, overlap_count := 0L]
hits <- sample(n, size = round(n * prob_overlap))
dt[hits, overlap_count := sample(1:4, length(hits), replace = TRUE, prob = c(0.55, 0.25, 0.15, 0.05))]

# Simulate continuous statistic (e.g. z-score or local rGM)
dt[, z_score := rnorm(n, mean = 0, sd = 0.5)]
dt[hits, z_score := rnorm(length(hits), mean = 1.8, sd = 0.8)]

fwrite(dt, output_file, sep = "\t")
cat("Successfully generated:", output_file, "with", n, "probes.\n")
