#!/usr/bin/env Rscript

# CURATED ENTRY POINT for the paper's M0/M1 Quadron regressions.
# Run --mode check first. The default never downloads data or fits a model.
# All computations use the existing BTEP scripts with explicit parameters;
# this runner does not contain an alternative implementation of the science.
suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
  library(yaml)
})

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_argument))
btep_root <- normalizePath(file.path(dirname(script_path), ".."))
source(file.path(btep_root, "replication/compare_results.R"))

options_list <- list(
  make_option("--mode", default = "check", help = "check, models, prepared, source, or verify [default %default]"),
  make_option("--config", default = file.path(btep_root, "replication/config.yaml")),
  make_option("--outdir", default = "", help = "Override output_root; relative paths are BTEP-relative"),
  make_option("--allow-downloads", action = "store_true", default = FALSE, dest = "allow_downloads")
)
options <- parse_args(OptionParser(option_list = options_list))
if (!options$mode %in% c("check", "models", "prepared", "source", "verify")) stop("Unknown --mode")
configuration <- read_yaml(options$config)

# canonical_path (from compare_results.R) resolves existing symlinks and rejects
# ambiguous '..' traversal through nonexistent parents before files are created.
resolve_path <- function(path) {
  if (is.null(path) || !nzchar(path)) return("")
  path <- path.expand(path)
  if (!startsWith(path, "/")) path <- file.path(btep_root, path)
  canonical_path(path)
}
required_settings <- c("source_root", "output_root", "expected_tables", "probe_manifest",
                       "g4_flank_bp", "sequence_flank_bp", "age_column", "sample_column",
                       "covariates", "minimum_samples", "absolute_tolerance", "relative_tolerance")
if (!all(required_settings %in% names(configuration))) stop("Incomplete replication configuration")
stopifnot(configuration$g4_flank_bp == 100L,
          configuration$sequence_flank_bp >= 0L,
          configuration$sequence_flank_bp == as.integer(configuration$sequence_flank_bp),
          identical(as.character(configuration$covariates), c("sex", "bmi")),
          configuration$age_column == "age", configuration$sample_column == "sample_id",
          configuration$minimum_samples == 2L)

source_root <- resolve_path(configuration$source_root)
expected_root <- resolve_path(configuration$expected_tables)
output_root <- resolve_path(if (nzchar(options$outdir)) options$outdir else configuration$output_root)
protected_roots <- c(source_root, expected_root, file.path(btep_root, "bin"),
                     file.path(btep_root, "replication"), file.path(btep_root, "G4_methylation_paper_files"))
check_output_root <- function(destination) {
  for (protected_root in protected_roots) {
    if (within_path(destination, protected_root) || within_path(protected_root, destination)) {
      stop("Output overlaps a protected input/script/manuscript path: ", destination)
    }
  }
}
check_output_root(output_root)

cohorts <- data.table(accession = c("GSE61257", "GSE61258", "GSE61259"),
                      tissue = c("adipose", "liver", "muscle"))
cohorts[, prefix := paste(accession, tissue, sep = "_")]
models <- paper_model_spec()
merged_name <- "g4_motifs_100bp_merged.bed"
unmerged_name <- "g4_motifs_100bp_unmerged.bed"
script_names <- c("download_geo_series_matrix.R", "prepare_g4_windows.R",
                  "analyze_g4_methylation_age.R", "compare_g4_effects_adjusted.R")

# Reference results are immutable. Check their recorded hashes before accepting
# them as the comparison target; do not refresh a target from a newly fitted run.
source_manifest <- fread(file.path(btep_root, "G4_methylation_paper_files/source_manifest.csv"))
expected_manifest <- source_manifest[type == "table"]
stopifnot(nrow(expected_manifest) == 6L)
expected_paths <- file.path(expected_root, basename(expected_manifest$bundled))
stopifnot(all(file.exists(expected_paths)),
          all(unname(tools::md5sum(expected_paths)) == expected_manifest$md5))

if (options$mode == "verify") {
  verify_paper_results(output_root, expected_root, file.path(output_root, "verification/numerical_comparison.csv"),
                       configuration$absolute_tolerance, configuration$relative_tolerance, protected_roots)
  quit(save = "no", status = 0L)
}

packages <- c("data.table", "optparse", "yaml", "ggplot2",
              "IlluminaHumanMethylation450kanno.ilmn12.hg19", "BSgenome.Hsapiens.UCSC.hg38",
              "GenomicRanges", "IRanges", "GenomeInfoDb", "Biostrings", "BSgenome")
installed <- vapply(packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
if (!all(installed)) stop("Missing packages: ", paste(packages[!installed], collapse = ", "),
                           ". See replication/install_packages.R")

check_file <- function(path) {
  if (!nzchar(path) || !file.exists(path) || dir.exists(path) || file.info(path)$size == 0) {
    stop("Required non-empty input missing: ", path)
  }
}
check_slopes <- function(directory) {
  for (cohort_index in seq_len(nrow(cohorts))) {
    path <- file.path(directory, paste0(cohorts$prefix[cohort_index], "_probe_slopes.csv.gz"))
    check_file(path)
    # Read only model-relevant columns for preflight. Requiring stored mean_beta
    # avoids the final script's potentially different prepared-matrix fallback.
    columns <- c("probe_id", "chr", "start", "end", "mean_beta", "slope",
                 "stable_motifs", "unstable_motifs", "cohort", "tissue")
    header <- names(fread(path, nrows = 0))
    if (!all(columns %in% header)) stop("Incomplete slope checkpoint: ", path)
    table <- fread(path, select = columns)
    stopifnot(uniqueN(table$probe_id) == nrow(table),
              all(table$tissue == cohorts$tissue[cohort_index]),
              all(table$cohort == cohorts$accession[cohort_index]),
              all(is.finite(table$start)), all(is.finite(table$end)),
              all(table$end >= table$start),
              all(table$stable_motifs >= 0), all(table$unstable_motifs >= 0))
  }
}
check_windows <- function(path, expected_columns) {
  check_file(path)
  if (ncol(fread(path, header = FALSE, nrows = 3L)) != expected_columns) {
    stop("Wrong window-checkpoint schema: ", path)
  }
}

prepared_dir <- file.path(source_root, "prepared_inputs")
windows_dir <- file.path(source_root, "g4_windows")
slopes_dir <- file.path(source_root, "per_tissue")
if (options$mode %in% c("check", "models")) {
  check_slopes(slopes_dir)
  check_windows(file.path(windows_dir, unmerged_name), 9L)
}
if (options$mode == "check") {
  cat("PASS: three Quadron slope checkpoints, unmerged windows, required packages, and six immutable expected tables.\n")
  cat("Candidate sequence flank:", configuration$sequence_flank_bp,
      "bp; historical invocation unknown, numerical replay is required.\n")
  cat("Original stable BED supplied:", nzchar(resolve_path(configuration$quadron_stable_bed)), "\n")
  cat("Original unstable BED supplied:", nzchar(resolve_path(configuration$quadron_unstable_bed)), "\n")
  cat("No data downloaded, models fitted, or files written.\n")
  cat("Next: Rscript BTEP/replication/run.R --mode models\n")
  quit(save = "no", status = 0L)
}

# Refuse reuse of any populated destination: a failed or completed run is
# evidence and must remain intact. For another attempt choose a NEW --outdir.
# 'verify' above is the only mode allowed to revisit existing model outputs.
if (file.exists(output_root) && !dir.exists(output_root)) stop("Output is an existing file")
if (dir.exists(output_root) && length(list.files(output_root, all.files = TRUE, no.. = TRUE))) {
  stop("Output directory is not empty; choose a new --outdir: ", output_root)
}
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
output_root <- normalizePath(output_root, mustWork = TRUE)
check_output_root(output_root)
dir.create(file.path(output_root, "provenance"))
dir.create(file.path(output_root, "logs"))
write_yaml(configuration, file.path(output_root, "provenance/config_requested.yaml"))
write_yaml(list(mode = options$mode, source_root = source_root, expected_root = expected_root,
                output_root = output_root, sequence_flank_bp = configuration$sequence_flank_bp,
                allow_downloads = options$allow_downloads), file.path(output_root, "provenance/run_resolved.yaml"))
writeLines(capture.output(sessionInfo()), file.path(output_root, "provenance/session_info.txt"))
fwrite(data.table(package = packages, version = vapply(packages, function(package) as.character(packageVersion(package)), "")),
       file.path(output_root, "provenance/package_versions.csv"))

input_records <- list()
record_files <- function(paths, role) {
  for (path in unique(paths)) {
    check_file(path)
    input_records[[length(input_records) + 1L]] <<- data.table(
      role = role, path = normalizePath(path), bytes = file.info(path)$size,
      md5 = unname(tools::md5sum(path))
    )
  }
  fwrite(rbindlist(input_records), file.path(output_root, "provenance/input_manifest.csv"))
}
record_files(expected_paths, "paper_expected_table")
code_paths <- c(file.path(btep_root, "bin", script_names), script_path,
                file.path(btep_root, "replication/compare_results.R"), normalizePath(options$config))
record_files(code_paths, "code_or_configuration")
dir.create(file.path(output_root, "provenance/source_used"))
stopifnot(all(file.copy(code_paths, file.path(output_root, "provenance/source_used"), overwrite = FALSE)))

run_step <- function(label, script, arguments) {
  log_path <- file.path(output_root, "logs", paste0(label, ".log"))
  script_path <- file.path(btep_root, "bin", script)
  rscript <- file.path(R.home("bin"), "Rscript")
  # Quote each shell argument separately: spaces in paths and the 'sex|bmi'
  # option must survive intact. Logs capture warnings/errors without hiding them.
  command_text <- paste(c(shQuote(rscript), shQuote(script_path), shQuote(arguments)), collapse = " ")
  cat(command_text, "\n", file = file.path(output_root, "provenance/commands.txt"), append = TRUE)
  message("Running ", label, "; log: ", log_path)
  status <- system2(rscript, shQuote(c(script_path, arguments)), stdout = log_path, stderr = log_path)
  if (status != 0L) {
    cat(paste(tail(readLines(log_path, warn = FALSE), 20L), collapse = "\n"), "\n")
    stop("Step failed: ", label, "; see ", log_path)
  }
}

if (options$mode == "source") {
  # This optional path requires original, user-supplied class BEDs. It does not
  # run Quadron or reconstruct its training data, scoring threshold or version.
  stable_bed <- resolve_path(configuration$quadron_stable_bed)
  unstable_bed <- resolve_path(configuration$quadron_unstable_bed)
  check_file(stable_bed)
  check_file(unstable_bed)
  record_files(c(stable_bed, unstable_bed), "original_quadron_class_files")
  windows_dir <- file.path(output_root, "g4_windows")
  run_step("01_windows", "prepare_g4_windows.R", c(
    "--stable", stable_bed, "--unstable", unstable_bed, "--flank", configuration$g4_flank_bp,
    "--long-quantile", configuration$long_quantile,
    "--long-mad-multiplier", configuration$long_mad_multiplier,
    "--long-min-bp", configuration$long_min_bp, "--outdir", windows_dir
  ))
  prepared_dir <- file.path(output_root, "prepared_inputs")
  for (cohort_index in seq_len(nrow(cohorts))) {
    accession <- cohorts$accession[cohort_index]
    cached_matrix <- resolve_path(configuration$series_matrices[[accession]])
    if (!nzchar(cached_matrix) && !options$allow_downloads) {
      stop("No cached series matrix for ", accession, "; downloads require --allow-downloads")
    }
    if (nzchar(cached_matrix)) record_files(cached_matrix, "geo_processed_series_matrix")
    run_step(paste0("02_import_", accession), "download_geo_series_matrix.R", c(
      "--accession", accession, "--tissue", cohorts$tissue[cohort_index],
      "--prefix", cohorts$prefix[cohort_index], "--outdir", prepared_dir,
      if (nzchar(cached_matrix)) c("--series-matrix", cached_matrix)
    ))
    if (!nzchar(cached_matrix)) {
      record_files(file.path(prepared_dir, paste0(accession, "_series_matrix.txt.gz")), "downloaded_geo_processed_matrix")
    }
  }
}

if (options$mode %in% c("prepared", "source")) {
  # Rebuilding upstream slopes adds assumptions that checkpoint replay avoids.
  # Unmapped probes are legitimate manifest entries and are excluded by the
  # historical analysis script. Validate the mapped subset without filling in
  # or shifting coordinates; also check it against the actual saved checkpoint.
  manifest_path <- resolve_path(configuration$probe_manifest)
  check_file(manifest_path)
  manifest <- fread(manifest_path)
  manifest_columns <- c("probe_id", "chr_hg38", "start_hg38", "end_hg38")
  if (!all(manifest_columns %in% names(manifest))) stop("Manifest must have explicit hg38 columns")
  mapped_manifest <- manifest[grepl("^chr([0-9]+|X|Y)$", chr_hg38)]
  if (uniqueN(manifest$probe_id) != nrow(manifest) ||
      nrow(mapped_manifest) == 0L || anyNA(manifest$probe_id) || any(!nzchar(manifest$probe_id)) ||
      !all(is.finite(mapped_manifest$start_hg38)) || !all(is.finite(mapped_manifest$end_hg38)) ||
      !all(mapped_manifest$start_hg38 >= 0 & mapped_manifest$end_hg38 >= mapped_manifest$start_hg38)) {
    stop("Manifest coordinates/uniqueness fail. Supply the actual validated historical manifest; do not guess a repair.")
  }
  fwrite(data.table(manifest_rows = nrow(manifest), mapped_rows = nrow(mapped_manifest),
                    unmapped_or_nonstandard_rows = nrow(manifest) - nrow(mapped_manifest)),
         file.path(output_root, "provenance/manifest_qc.csv"))
  for (prefix in cohorts$prefix) {
    saved_coordinates <- fread(file.path(source_root, "per_tissue", paste0(prefix, "_probe_slopes.csv.gz")),
                                select = c("probe_id", "chr", "start", "end"))
    coordinates <- merge(saved_coordinates, mapped_manifest[, ..manifest_columns], by = "probe_id", all.x = TRUE)
    if (nrow(coordinates) != nrow(saved_coordinates) ||
        anyNA(coordinates[, ..manifest_columns]) ||
        !all(coordinates$chr == coordinates$chr_hg38 & coordinates$start == coordinates$start_hg38 &
               coordinates$end == coordinates$end_hg38)) {
      stop("Manifest does not match saved probe coordinates for ", prefix)
    }
  }
  record_files(manifest_path, "probe_manifest")
  check_windows(file.path(windows_dir, merged_name), 11L)
  record_files(file.path(windows_dir, merged_name), "merged_quadron_neighborhoods")
  slopes_dir <- file.path(output_root, "per_tissue")
  for (cohort_index in seq_len(nrow(cohorts))) {
    prefix <- cohorts$prefix[cohort_index]
    beta_path <- file.path(prepared_dir, paste0(prefix, "_beta.rds"))
    metadata_path <- file.path(prepared_dir, paste0(prefix, "_metadata.csv"))
    record_files(c(beta_path, metadata_path), "prepared_methylation_inputs")
    metadata <- fread(metadata_path)
    stopifnot(all(c("sample_id", "age", "sex", "bmi") %in% names(metadata)),
              uniqueN(metadata$sample_id) == nrow(metadata))
    eligible <- metadata[complete.cases(metadata[, .(age, sex, bmi)])]
    design <- model.matrix(~ age + sex + bmi, data = eligible)
    if (qr(design)$rank != ncol(design)) stop("Rank-deficient age/sex/BMI design for ", prefix)
    run_step(paste0("03_slopes_", prefix), "analyze_g4_methylation_age.R", c(
      "--beta", beta_path, "--metadata", metadata_path, "--manifest", manifest_path,
      "--g4-windows", file.path(windows_dir, merged_name),
      "--age-col", configuration$age_column, "--sample-col", configuration$sample_column,
      "--covariates", paste(configuration$covariates, collapse = "|"),
      "--min-samples", configuration$minimum_samples, "--cohort", cohorts$accession[cohort_index],
      "--tissue", cohorts$tissue[cohort_index], "--prefix", prefix, "--slopes-only", "--outdir", slopes_dir
    ))
    # Compare the full regenerated checkpoint, including annotation and mean
    # beta, BEFORE attributing a later mismatch to the final regression.
    checkpoint_name <- paste0(prefix, "_probe_slopes.csv.gz")
    comparison <- compare_tables(file.path(slopes_dir, checkpoint_name),
                                  file.path(source_root, "per_tissue", checkpoint_name), "probe_id",
                                  configuration$absolute_tolerance, configuration$relative_tolerance)
    fwrite(comparison, file.path(output_root, "provenance", paste0(prefix, "_slope_parity.csv")))
    if (!all(comparison$passed)) stop("Rebuilt slope checkpoint differs for ", prefix, "; see provenance report")
  }
  check_slopes(slopes_dir)
}

check_windows(file.path(windows_dir, unmerged_name), 9L)
record_files(file.path(windows_dir, unmerged_name), "unmerged_quadron_checkpoint")
staged_dir <- file.path(output_root, "staged_slopes")
dir.create(staged_dir)
slope_paths <- file.path(slopes_dir, paste0(cohorts$prefix, "_probe_slopes.csv.gz"))
record_files(slope_paths, "annotated_age_slope_checkpoint")
stopifnot(all(file.copy(slope_paths, staged_dir, overwrite = FALSE)))

# EXACT final-model entry points. Earlier Nextflow result families are neither
# called nor used as evidence. --minimal-model selects the owning code branch.
# No structure-file path is supplied: formula equality alone would not prevent
# optional upstream cross-reactive filtering in the shared legacy script.
for (model_index in seq_len(nrow(models))) {
  run_step(paste0("04_", models$model[model_index]), "compare_g4_effects_adjusted.R", c(
    "--prepared-dir", prepared_dir, "--slopes-dir", staged_dir,
    "--g4-unmerged", file.path(windows_dir, unmerged_name),
    "--minimal-model", "true", "--include-strand-g-richness", models$g_richness[model_index],
    "--include-chromatin", "false", "--include-g4-architecture", "false",
    "--stratify-g4-density", "false", "--exclude-cross-reactive", "false",
    "--gc-flank-bp", configuration$sequence_flank_bp,
    "--outdir", file.path(output_root, models$directory[model_index]),
    "--prefix", models$prefix[model_index]
  ))
}

# Never rewrite the manuscript or its expected figures as part of numerical
# reproduction. The current plotter may differ cosmetically from the retained
# PNGs (for example, class brackets and fonts); table parity is the science gate.
verify_paper_results(output_root, expected_root, file.path(output_root, "verification/numerical_comparison.csv"),
                     configuration$absolute_tolerance, configuration$relative_tolerance, protected_roots)
for (model_index in seq_len(nrow(models))) {
  for (extension in c("png", "pdf")) {
    check_file(file.path(output_root, models$directory[model_index],
                         paste0(models$prefix[model_index], "_g4_terms.", extension)))
  }
}
message("Completed a numerically verified ", options$mode, " replay in ", output_root)
message("Original manuscript, original results, and original figure files were not overwritten.")