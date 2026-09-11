# Offline fixture tests for the three upstream steps. These exercise real
# command-line scripts, not reimplemented approximations of their behavior.
# The values and intervals below are synthetic, not additional study evidence.
script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_argument))
btep_root <- dirname(dirname(dirname(script_path)))
scratch <- tempfile("quadron_inputs_test_")
dir.create(scratch)

run_tests <- function() {
  on.exit(unlink(scratch, recursive = TRUE), add = TRUE)
  run_script <- function(name, arguments) {
    output <- system2(file.path(R.home("bin"), "Rscript"),
                      shQuote(c(file.path(btep_root, "bin", name), arguments)),
                      stdout = TRUE, stderr = TRUE)
    status <- attr(output, "status")
    if (!is.null(status) && status != 0L) stop(paste(output, collapse = "\n"))
  }
  set.seed(41)
  metadata <- data.frame(sample_id = paste0("GSM_TEST", seq_len(8)),
                          age = c(21, 34, 29, 60, 51, 75, 44, 83),
                          sex = rep(c("female", "male"), 4),
                          bmi = c(20, 33, 22, 41, 28, 35, 30, 24))
  beta <- rbind(cg_test1 = 0.3 + 0.001 * metadata$age + rnorm(8, 0, 0.01),
                cg_test2 = 0.7 - 0.002 * metadata$age + rnorm(8, 0, 0.01))
  colnames(beta) <- metadata$sample_id
  beta[2, 3] <- NA_real_

  # Match the minimal GEO text structure consumed by the parser, including
  # quoted sample IDs and numeric characteristics. The explicit local path is
  # what guarantees that this test makes no GEO request.
  geo_row <- function(key, values) paste(c(key, sprintf('"%s"', values)), collapse = "\t")
  matrix_rows <- vapply(seq_len(nrow(beta)), function(probe_index) {
    paste(c(rownames(beta)[probe_index], format(beta[probe_index, ], digits = 17, trim = TRUE)), collapse = "\t")
  }, "")
  lines <- c(geo_row("!Sample_geo_accession", metadata$sample_id),
             geo_row("!Sample_characteristics_ch1", paste("age:", metadata$age)),
             geo_row("!Sample_characteristics_ch1", paste("sex:", metadata$sex)),
             geo_row("!Sample_characteristics_ch1", paste("bmi:", metadata$bmi)),
             "!series_matrix_table_begin", geo_row("ID_REF", metadata$sample_id),
             matrix_rows, "!series_matrix_table_end")
  series_matrix <- file.path(scratch, "series_matrix.txt.gz")
  connection <- gzfile(series_matrix, "wt")
  writeLines(lines, connection)
  close(connection)
  prepared_dir <- file.path(scratch, "prepared")
  run_script("download_geo_series_matrix.R", c("--accession", "GSE61257", "--tissue", "adipose",
                                              "--prefix", "test", "--series-matrix", series_matrix,
                                              "--outdir", prepared_dir))
  imported_beta <- readRDS(file.path(prepared_dir, "test_beta.rds"))
  stopifnot(isTRUE(all.equal(imported_beta, beta)),
            identical(data.table::fread(file.path(prepared_dir, "test_metadata.csv"))$sample_id, metadata$sample_id))

  # The padded intervals [0,210] and [10,220] must become ONE mixed
  # neighborhood. Both class labels are retained; no exclusive priority rule.
  stable <- data.frame(chr = "chr1", start = 100L, end = 110L, motif = "synthetic",
                       score = 1, strand = "+", g4_score = 20)
  unstable <- stable
  unstable$start <- 110L
  unstable$end <- 120L
  unstable$strand <- "-"
  stable_path <- file.path(scratch, "stable.bed")
  unstable_path <- file.path(scratch, "unstable.bed")
  data.table::fwrite(stable, stable_path, sep = "\t", col.names = FALSE)
  data.table::fwrite(unstable, unstable_path, sep = "\t", col.names = FALSE)
  windows_dir <- file.path(scratch, "windows")
  run_script("prepare_g4_windows.R", c("--stable", stable_path, "--unstable", unstable_path,
                                       "--flank", "100", "--outdir", windows_dir))
  merged_path <- file.path(windows_dir, "g4_motifs_100bp_merged.bed")
  merged <- data.table::fread(merged_path, header = FALSE)
  stopifnot(nrow(merged) == 1L, ncol(merged) == 11L, merged$V2 == 0L,
            merged$V3 == 220L, merged$V7 == 1L, merged$V8 == 1L)
  stopifnot(ncol(data.table::fread(file.path(windows_dir, "g4_motifs_100bp_unmerged.bed"), header = FALSE)) == 9L)

  manifest_path <- file.path(scratch, "manifest.csv")
  data.table::fwrite(data.frame(probe_id = rownames(beta), chr_hg38 = "chr1",
                                start_hg38 = c(100L, 500L), end_hg38 = c(101L, 501L)), manifest_path)
  output_dir <- file.path(scratch, "slopes")
  run_script("analyze_g4_methylation_age.R", c(
    "--beta", file.path(prepared_dir, "test_beta.rds"),
    "--metadata", file.path(prepared_dir, "test_metadata.csv"), "--manifest", manifest_path,
    "--g4-windows", merged_path, "--cohort", "GSE61257", "--tissue", "adipose",
    "--prefix", "test", "--covariates", "sex|bmi", "--slopes-only", "--outdir", output_dir
  ))
  slopes <- data.table::fread(file.path(output_dir, "test_probe_slopes.csv.gz"))
  for (probe_index in seq_len(nrow(beta))) {
    expected <- summary(lm(beta[probe_index, ] ~ age + sex + bmi, data = metadata))$coefficients["age", ]
    stopifnot(abs(slopes$slope[probe_index] - expected["Estimate"]) < 1e-12,
              abs(slopes$se[probe_index] - expected["Std. Error"]) < 1e-12,
              abs(slopes$mean_beta[probe_index] - mean(beta[probe_index, ], na.rm = TRUE)) < 1e-12)
  }
  stopifnot(slopes$has_stable_g4[1], slopes$has_unstable_g4[1],
            !slopes$has_stable_g4[2], !slopes$has_unstable_g4[2],
            identical(list.files(output_dir), "test_probe_slopes.csv.gz"))
  message("PASS: offline GEO import, 100-bp mixed-window merge, missing-beta OLS/SE/mean, overlapping class indicators, and slope-only output.")
}

run_tests()