# Numeric reproduction is different from document validation. This file checks
# re-estimated statistics against the six unrounded tables bundled with the
# paper. It never refits a model, changes an expected result, or judges success
# solely by the direction of an effect or whether q is below 0.05.

# Resolve real ancestors before appending NEW path components. A '..' that
# crosses a nonexistent directory has ambiguous physical meaning: creating the
# directory later could redirect a write into an existing input tree. Reject
# it instead of doing a purely textual cleanup. Existing paths are normalized
# by the filesystem, so symlink semantics are preserved.
canonical_path <- function(path) {
  path <- path.expand(path)
  if (!startsWith(path, "/")) path <- file.path(getwd(), path)
  if (file.exists(path) || dir.exists(path)) return(normalizePath(path, mustWork = TRUE))
  link_target <- Sys.readlink(path)
  if (!is.na(link_target) && nzchar(link_target)) stop("Dangling symlink in path: ", path)
  if (basename(path) %in% c(".", "..")) stop("Unresolved traversal through a missing directory: ", path)
  parent <- dirname(path)
  if (parent == path) stop("Cannot resolve path: ", path)
  file.path(canonical_path(parent), basename(path))
}

within_path <- function(path, parent) identical(path, parent) || startsWith(path, paste0(parent, "/"))

safe_report_path <- function(report_path, actual_root, expected_root, protected_paths = character()) {
  # A verification report may replace a prior REGULAR report, but must never
  # follow a symlink into a reference CSV, script, or regenerated model table.
  # Check both the complete filename and resolved ancestors before any write.
  link_target <- Sys.readlink(report_path)
  if (!is.na(link_target) && nzchar(link_target)) stop("Report file must not be a symlink: ", report_path)
  destination <- canonical_path(report_path)
  actual_root <- canonical_path(actual_root)
  if (!within_path(destination, actual_root) || identical(destination, actual_root)) {
    stop("Report destination escapes the reproduction output directory: ", report_path)
  }
  protected <- c(expected_root, protected_paths,
                 file.path(actual_root, paper_model_spec()$directory))
  for (protected_path in protected) {
    if (within_path(destination, canonical_path(protected_path))) {
      stop("Report destination overlaps protected evidence: ", report_path)
    }
  }
  destination
}

paper_model_spec <- function() {
  data.frame(
    model = c("M0", "M1"),
    directory = c("minimal_model", "minimal_model_strand_g_richness"),
    prefix = c("quadron_minimal_stable_unstable_mean_beta", "quadron_minimal_strand_g_richness"),
    g_richness = c("false", "true"),
    stringsAsFactors = FALSE
  )
}

compare_tables <- function(actual_path, expected_path, keys,
                           absolute_tolerance = 1e-12, relative_tolerance = 1e-7) {
  stopifnot(is.finite(absolute_tolerance), absolute_tolerance >= 0,
            is.finite(relative_tolerance), relative_tolerance >= 0)
  diagnostic <- function(column, mismatches, detail = "", max_absolute_error = NA_real_) {
    data.table::data.table(
      file = basename(expected_path), column = column, passed = mismatches == 0L,
      mismatches = as.integer(mismatches), max_absolute_error = max_absolute_error,
      detail = detail
    )
  }
  if (!file.exists(actual_path) || !file.exists(expected_path)) {
    return(diagnostic("file", 1L, "Actual or expected file is missing"))
  }
  actual <- data.table::fread(actual_path)
  expected <- data.table::fread(expected_path)
  if (!setequal(names(actual), names(expected))) {
    return(diagnostic("columns", 1L, "Column sets differ; no silent column dropping"))
  }
  if (!all(keys %in% names(expected)) ||
      data.table::uniqueN(actual, by = keys) != nrow(actual) ||
      data.table::uniqueN(expected, by = keys) != nrow(expected)) {
    return(diagnostic("keys", 1L, "Missing or duplicate biological/statistical row keys"))
  }

  # Sort by semantic keys, never assume the same row order or compare floating
  # point results by CSV byte identity. Different R versions may serialize the
  # same numbers differently. Extra/missing rows are always a failure.
  data.table::setorderv(actual, keys)
  data.table::setorderv(expected, keys)
  if (nrow(actual) != nrow(expected) ||
      !identical(as.data.frame(actual[, keys, with = FALSE]),
                 as.data.frame(expected[, keys, with = FALSE]))) {
    return(diagnostic("keys", 1L, "Row identities differ"))
  }

  diagnostics <- lapply(names(expected), function(column) {
    observed <- actual[[column]]
    target <- expected[[column]]
    missing_disagrees <- xor(is.na(observed), is.na(target))
    present <- !is.na(observed) & !is.na(target)
    mismatches <- missing_disagrees
    max_error <- NA_real_
    if (is.numeric(target) && is.numeric(observed)) {
      finite <- present & is.finite(observed) & is.finite(target)
      difference <- abs(observed[finite] - target[finite])
      if (length(difference)) max_error <- max(difference)

      # A liberal absolute tolerance would incorrectly equate P=1e-100 with
      # P=1e-15. Use relative error for P/q and exact equality for probe counts.
      if (column %in% c("p_value", "q_value")) {
        allowance <- relative_tolerance * pmax(abs(target[finite]), abs(observed[finite]))
      } else if (column == "n_probes") {
        allowance <- rep(0, sum(finite))
      } else {
        allowance <- absolute_tolerance + relative_tolerance * abs(target[finite])
      }
      mismatches[finite] <- difference > allowance
      nonfinite <- present & !finite
      mismatches[nonfinite] <- observed[nonfinite] != target[nonfinite]
    } else {
      # Formulas, tissue names, predictor labels, and estimability flags must
      # match exactly. R integer/double differences matter only for row keys.
      mismatches[present] <- as.character(observed[present]) != as.character(target[present])
    }
    detail <- if (any(mismatches)) {
      paste("First differing row:", paste(unlist(expected[which(mismatches)[1], keys, with = FALSE]), collapse = "/"))
    } else ""
    diagnostic(column, sum(mismatches), detail, max_error)
  })
  data.table::rbindlist(diagnostics)
}

verify_paper_results <- function(actual_root, expected_root, report_path,
                                 absolute_tolerance = 1e-12, relative_tolerance = 1e-7,
                                 protected_paths = character()) {
  report_path <- safe_report_path(report_path, actual_root, expected_root, protected_paths)
  models <- paper_model_spec()
  table_keys <- list(
    coefficients = c("tissue", "density_group", "term"),
    g4_terms = c("tissue", "density_group", "term"),
    stable_minus_unstable = c("tissue", "density_group", "contrast")
  )
  comparisons <- list()
  for (model_index in seq_len(nrow(models))) {
    for (table_name in names(table_keys)) {
      filename <- paste0(models$prefix[model_index], "_", table_name, ".csv")
      comparison <- compare_tables(
        file.path(actual_root, models$directory[model_index], filename),
        file.path(expected_root, filename), table_keys[[table_name]],
        absolute_tolerance, relative_tolerance
      )
      comparison[, model := models$model[model_index]]
      comparisons[[length(comparisons) + 1L]] <- comparison
    }
  }
  report <- data.table::rbindlist(comparisons)
  dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)
  report_path <- safe_report_path(report_path, actual_root, expected_root, protected_paths)
  # Write then rename in the same directory. This avoids partially replacing an
  # existing report, and a late symlink substitution cannot truncate its target.
  temporary_report <- tempfile("comparison_", tmpdir = dirname(report_path), fileext = ".csv")
  on.exit(unlink(temporary_report), add = TRUE)
  data.table::fwrite(report, temporary_report)
  if (!file.rename(temporary_report, report_path)) stop("Cannot finalize verification report: ", report_path)
  if (!all(report$passed)) {
    print(report[passed == FALSE])
    stop("Numerical reproduction failed. Inspect ", report_path,
         "; the paper's expected results have NOT been changed.")
  }
  message("PASS: all columns of all six regenerated model tables match the paper within the declared tolerances.")
  invisible(report)
}