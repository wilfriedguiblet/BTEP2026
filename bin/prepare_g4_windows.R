#!/usr/bin/env Rscript

# PURPOSE: build the coordinate checkpoint consumed by probe annotation.
# For the paper, pass the original Quadron stable and unstable files explicitly.
# This script does NOT run Quadron or choose a stability-score threshold; class
# labels come solely from the supplied file identity. Retain source hashes and
# the original threshold/version separately if available.
# Coordinate caution: standard BED uses zero-based, half-open intervals. This
# historical implementation preserves the supplied numbers. Downstream
# foverlaps uses closed intervals; do not silently shift coordinates during a
# reproduction. A coordinate-corrected analysis needs separate validation.

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

option_list <- list(
  make_option("--stable", type = "character"),
  make_option("--unstable", type = "character"),
  make_option("--flank", type = "integer", default = 100L),
  make_option("--long-quantile", type = "double", default = 0.99, dest = "long_quantile"),
  make_option("--long-mad-multiplier", type = "double", default = 3, dest = "long_mad_multiplier"),
  make_option("--long-min-bp", type = "integer", default = 0L, dest = "long_min_bp"),
  make_option("--outdir", type = "character", default = ".")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$stable) || is.null(opt$unstable)) stop("--stable and --unstable are required")
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
flank_label <- paste0(opt$flank, "bp")

read_g4 <- function(path, class_label) {
  # The historical input schema includes chromosome/start/end plus motif,
  # score, strand, and g4_score fields. Check the actual input columns before
  # substituting another BED export; the prepared checkpoints have other schemas.
  x <- fread(path, header = FALSE)
  bed_names <- c("chr", "start", "end", "motif", "score", "strand", "g4_score")
  setnames(x, bed_names[seq_len(ncol(x))])
  x <- x[grepl("^chr([0-9]+|X|Y)$", x[["chr"]])]
  x[, `:=`(
    start = as.integer(start),
    end = as.integer(end),
    g4_class = class_label,
    motif_id = sprintf("%s_%07d", class_label, seq_len(.N))
  )]
  x[]
}

core <- rbindlist(list(read_g4(opt$stable, "stable"), read_g4(opt$unstable, "unstable")), fill = TRUE)
setorder(core, chr, start, end)
fwrite(core[, .(chr, start, end, motif_id, g4_class, strand, g4_score, motif)],
       file.path(opt$outdir, "g4_motifs_core.bed"), sep = "\t", col.names = FALSE)

windows <- copy(core)
# G4 padding is 100 bp for the paper. This is NOT the 1000-bp default flank
# used later to measure G-richness around CpGs. Negative starts are clamped to
# zero; chromosome-end clipping is not added by this historical window builder.
windows[, `:=`(
  core_start = start,
  core_end = end,
  start = pmax(0L, start - opt$flank),
  end = end + opt$flank
)]
setorder(windows, chr, start, end)
fwrite(windows[, .(chr, start, end, motif_id, g4_class, strand, g4_score, core_start, core_end)],
  file.path(opt$outdir, paste0("g4_motifs_", flank_label, "_unmerged.bed")), sep = "\t", col.names = FALSE)

reduce_chr <- function(dt) {
  # This is a union of intervals, not a nearest-motif assignment. A chain of
  # overlapping or abutting motifs can make a wide neighborhood. Each resulting
  # neighborhood carries counts from ALL contributing stable/unstable motifs.
  dt <- copy(dt)[order(start, end)]
  groups <- integer(nrow(dt))
  group_id <- 0L
  current_end <- -1L
  for (i in seq_len(nrow(dt))) {
    if (dt$start[i] > current_end) {
      group_id <- group_id + 1L
      current_end <- dt$end[i]
    } else {
      current_end <- max(current_end, dt$end[i])
    }
    groups[i] <- group_id
  }
  dt[, group_id := groups]
  dt[, list(
    start = min(get("start")),
    end = max(get("end")),
    motif_count = length(get("start")),
    stable_motifs = sum(get("g4_class") == "stable"),
    unstable_motifs = sum(get("g4_class") == "unstable")
  ), by = "group_id"]
}

merged <- windows[, reduce_chr(.SD), by = chr]
merged[, width := end - start]
# Long-cluster labels are checkpoint metadata. The paper's two minimal models
# retain these neighborhoods but do not include width/density as covariates.
# mad(..., constant=1) is the raw median absolute deviation, not R's default
# normal-consistency-scaled MAD. Preserve it to reproduce these annotations.
width_median <- median(merged$width)
width_mad <- mad(merged$width, constant = 1)
quantile_cutoff <- as.numeric(quantile(merged$width, probs = opt$long_quantile, names = FALSE))
mad_cutoff <- width_median + opt$long_mad_multiplier * width_mad
long_width_cutoff <- max(quantile_cutoff, mad_cutoff, opt$long_min_bp)
merged[, base_class := fifelse(stable_motifs > 0 & unstable_motifs > 0, "stable_unstable_overlap",
                               fifelse(stable_motifs > 0, "stable_only", "unstable_only"))]
# Keep the base class even when the display class becomes long_g4_cluster:
# later stable_motifs > 0 and unstable_motifs > 0 indicators can BOTH be true.
merged[, is_long_cluster := width >= long_width_cutoff]
merged[, class_label := fifelse(is_long_cluster, "long_g4_cluster", base_class)]
setorder(merged, chr, start, end)
merged[, window_id := sprintf("g4win_%07d", seq_len(.N))]
fwrite(merged[, .(chr, start, end, window_id, class_label, motif_count, stable_motifs, unstable_motifs,
                  width, is_long_cluster, base_class)],
  file.path(opt$outdir, paste0("g4_motifs_", flank_label, "_merged.bed")), sep = "\t", col.names = FALSE)
fwrite(merged[is_long_cluster == TRUE,
              .(chr, start, end, window_id, class_label, motif_count, stable_motifs, unstable_motifs,
                width, is_long_cluster, base_class)],
  file.path(opt$outdir, paste0("g4_motifs_", flank_label, "_long_clusters.bed")), sep = "\t", col.names = FALSE)

summary <- data.table(
  flank_bp = opt$flank,
  long_quantile = opt$long_quantile,
  long_mad_multiplier = opt$long_mad_multiplier,
  long_min_bp = opt$long_min_bp,
  long_width_cutoff = long_width_cutoff,
  median_merged_window_bp = width_median,
  mad_merged_window_bp = width_mad,
  quantile_merged_window_bp = quantile_cutoff,
  core_motifs = nrow(core),
  stable_core_motifs = sum(core$g4_class == "stable"),
  unstable_core_motifs = sum(core$g4_class == "unstable"),
  unmerged_windows = nrow(windows),
  merged_windows = nrow(merged),
  overlapping_or_adjacent_windows_collapsed = nrow(windows) - nrow(merged),
  long_cluster_windows = sum(merged$is_long_cluster),
  long_cluster_total_bp = sum(merged[is_long_cluster == TRUE, width]),
  windows_with_stable_and_unstable = sum(merged$base_class == "stable_unstable_overlap")
)
fwrite(summary, file.path(opt$outdir, "g4_window_overlap_summary.csv"))
