#!/usr/bin/env Rscript

# PACKAGE VERIFICATION EVIDENCE ONLY. Run after the two documented local
# replays have completed. This copies diagnostic records into a non-ignored
# directory for review; it never copies new coefficients over paper results.
# The full commands, code snapshots, and logs remain in the original run dirs.
suppressPackageStartupMessages(library(data.table))
script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_argument))
btep_root <- normalizePath(file.path(dirname(script_path), ".."))
verification_dir <- file.path(btep_root, "replication/verification")

collect_validation <- function() {
  run_names <- c(checkpoint = "quadron_paper_reproduction",
                 prepared = "quadron_paper_reproduction_prepared")
  # Check both runs before publishing any summary. A missing/failed gate is
  # not converted into a positive reproduction statement by this packager.
  for (run_name in run_names) {
    report <- fread(file.path(btep_root, "results", run_name, "verification/numerical_comparison.csv"))
    stopifnot(nrow(report) > 0L, all(report$passed), uniqueN(report$file) == 6L)
  }
  dir.create(verification_dir, showWarnings = FALSE)
  for (route in names(run_names)) {
    run_dir <- file.path(btep_root, "results", run_names[[route]])
    report <- fread(file.path(run_dir, "verification/numerical_comparison.csv"))
    summary <- report[, .(
      columns_checked = .N, all_passed = all(passed),
      maximum_absolute_error = max(max_absolute_error, na.rm = TRUE)
    ), by = .(model, file)]
    fwrite(summary, file.path(verification_dir, paste0(route, "_model_parity.csv")))
    # Paths in these manifests identify the validation machine's files; hashes
    # are the portable identity evidence. They are not public download URLs.
    copied <- file.copy(file.path(run_dir, "provenance/input_manifest.csv"),
                         file.path(verification_dir, paste0(route, "_input_manifest.csv")), overwrite = TRUE)
    stopifnot(copied)
    if (route == "prepared") {
      slope_reports <- list.files(file.path(run_dir, "provenance"), pattern = "_slope_parity[.]csv$",
                                   full.names = TRUE)
      stopifnot(length(slope_reports) == 3L)
      for (slope_report in slope_reports) {
        stopifnot(all(fread(slope_report)$passed))
        stopifnot(file.copy(slope_report, verification_dir, overwrite = TRUE))
      }
    }
    message(route, ": six model tables passed; maximum parsed numeric error = ",
            max(summary$maximum_absolute_error))
  }

  # Parsing is a narrow code-quality gate, not a substitute for the executed
  # model comparisons and fixture tests recorded elsewhere in this directory.
  scripts <- c(list.files(dirname(script_path), pattern = "[.]R$", recursive = TRUE, full.names = TRUE),
               file.path(btep_root, "bin", c("download_geo_series_matrix.R", "prepare_g4_windows.R",
                                             "analyze_g4_methylation_age.R", "compare_g4_effects_adjusted.R")),
               file.path(btep_root, c("render_G4_methylation_paper.R", "validate_G4_methylation_paper.R")))
  for (script in scripts) parse(script)
  stopifnot(yaml::read_yaml(file.path(btep_root, "replication/config.yaml"))$sequence_flank_bp == 1000)
  message("PASS: both replay summaries and three upstream parity reports packaged; all curated R scripts parse.")
}

collect_validation()