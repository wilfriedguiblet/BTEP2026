# Regression tests for preservation of source data. Use only temporary files;
# the real paper, reference tables and existing run outputs are never touched.
script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_argument))
source(file.path(dirname(dirname(script_path)), "compare_results.R"))

run_tests <- function() {
  scratch <- tempfile("quadron_paths_test_")
  dir.create(scratch)
  on.exit(unlink(scratch, recursive = TRUE), add = TRUE)
  scratch <- normalizePath(scratch)
  reference_dir <- file.path(scratch, "protected")
  output_dir <- file.path(scratch, "output")
  dir.create(reference_dir)
  dir.create(output_dir)
  reference_file <- file.path(reference_dir, "reference.csv")
  writeLines("unchanged_reference", reference_file)
  original_hash <- unname(tools::md5sum(reference_file))
  fails <- function(expression) inherits(tryCatch(force(expression), error = identity), "error")

  # A nonexistent parent followed by '..' must not gain new meaning after
  # dir.create(recursive=TRUE). Existing, unambiguous '..' remains usable for
  # legitimate input paths such as the workspace's shared reference manifest.
  unsafe <- file.path(scratch, "not_created", "..", "protected")
  stopifnot(fails(canonical_path(unsafe)), !dir.exists(file.path(scratch, "not_created")))
  stopifnot(identical(canonical_path(file.path(output_dir, "..", "protected")), reference_dir))
  expected_report <- file.path(output_dir, "verification", "report.csv")
  stopifnot(identical(safe_report_path(expected_report, output_dir, reference_dir), expected_report))

  # File and directory symlinks must not allow a report write to escape the
  # output tree or to replace one of the six reference/model tables.
  dir.create(file.path(output_dir, "verification"))
  report_link <- file.path(output_dir, "verification", "report.csv")
  stopifnot(file.symlink(reference_file, report_link))
  stopifnot(fails(verify_paper_results(output_dir, reference_dir, report_link)))
  unlink(report_link)
  directory_link <- file.path(output_dir, "elsewhere")
  stopifnot(file.symlink(reference_dir, directory_link))
  stopifnot(fails(safe_report_path(file.path(directory_link, "reference.csv"), output_dir, reference_dir)))
  stopifnot(fails(safe_report_path(file.path(output_dir, "minimal_model", "table.csv"), output_dir, reference_dir)))
  dangling <- file.path(output_dir, "dangling")
  stopifnot(file.symlink(file.path(scratch, "absent"), dangling), fails(canonical_path(dangling)))
  stopifnot(identical(unname(tools::md5sum(reference_file)), original_hash))
  message("PASS: missing-parent traversal, report-file symlinks, escaping directory symlinks, model-table destinations, and dangling links rejected without changing reference data.")
}

run_tests()