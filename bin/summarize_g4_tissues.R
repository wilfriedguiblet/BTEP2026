#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

option_list <- list(
  make_option("--summaries", type = "character"),
  make_option("--tests", type = "character"),
  make_option("--regressions", type = "character", default = ""),
  make_option("--overlap-summary", type = "character", dest = "overlap_summary"),
  make_option("--outdir", type = "character", default = ".")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

summary_files <- strsplit(opt$summaries, "[, ]+")[[1]]
test_files <- strsplit(opt$tests, "[, ]+")[[1]]
regression_files <- strsplit(opt$regressions, "[, ]+")[[1]]
summary_dt <- rbindlist(lapply(summary_files[nzchar(summary_files)], fread), fill = TRUE)
test_dt <- rbindlist(lapply(test_files[nzchar(test_files)], fread), fill = TRUE)
regression_dt <- rbindlist(lapply(regression_files[nzchar(regression_files)], fread), fill = TRUE)
fwrite(summary_dt, file.path(opt$outdir, "btep_g4_cross_tissue_summary.csv"))
fwrite(test_dt, file.path(opt$outdir, "btep_g4_cross_tissue_tests.csv"))
fwrite(regression_dt, file.path(opt$outdir, "btep_g4_cross_tissue_regression_effects.csv"))

overlap <- if (!is.null(opt$overlap_summary) && file.exists(opt$overlap_summary)) fread(opt$overlap_summary) else data.table()
readme <- c(
  "# BTEP G4 methylation-aging experiment",
  "",
  "This workflow tests whether CpG methylation age-slopes are larger inside annotated G4 motifs plus 1 kb flanks than outside those coordinates.",
  "",
  "## Coordinate-overlap policy",
  "",
  "Primary coordinates of interest are built by expanding every stable and unstable G4Hunter motif by 1 kb on both sides, then merging overlapping or adjacent expanded intervals into non-overlapping windows. Each CpG can therefore enter the primary G4-window analysis only once. The merged BED retains the number of stable and unstable motifs contributing to each window, with class labels `stable_only`, `unstable_only`, or `stable_unstable_overlap`.",
  "",
  "This avoids inflated evidence from double-counting CpGs in dense G4 clusters. For secondary analyses, use the unmerged BED to ask motif-level questions, but assign CpGs by nearest motif or fractional weights when windows overlap.",
  "",
  "## Generated outputs",
  "",
  "- `btep_g4_cross_tissue_summary.csv`: per-tissue slope summaries by G4-window class.",
  "- `btep_g4_cross_tissue_tests.csv`: Wilcoxon tests comparing G4-window CpGs to outside CpGs.",
  "- `btep_g4_cross_tissue_regression_effects.csv`: regression estimates for overlapping G4 predictors and age-slope responses.",
  "- Per-tissue `*_probe_slopes.csv.gz`: probe-level age slopes and G4-window annotations."
)
if (nrow(overlap)) {
  readme <- c(readme, "", "## G4 window counts", "", capture.output(print(overlap)))
}
writeLines(readme, file.path(opt$outdir, "btep_g4_experiment_readme.md"))
