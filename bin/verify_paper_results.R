# Run only as the final Nextflow task. The six expected CSVs are frozen paper
# assets, never refreshed from a new computation. A successful fit is not a
# successful reproduction unless every row, formula and statistic agrees.
suppressPackageStartupMessages(library(data.table))

compare_tables <- function(actual_path, expected_path, keys,
                           absolute_tolerance = 1e-12, relative_tolerance = 1e-7) {
  diagnostic <- function(column, mismatches, detail = "", max_absolute_error = NA_real_) {
    data.table(file = basename(expected_path), column = column, passed = mismatches == 0L,
               mismatches = as.integer(mismatches), max_absolute_error = max_absolute_error,
               detail = detail)
  }
  if (!file.exists(actual_path) || !file.exists(expected_path)) {
    return(diagnostic("file", 1L, "Actual or expected file is missing"))
  }
  actual <- fread(actual_path)
  expected <- fread(expected_path)
  if (!setequal(names(actual), names(expected))) {
    return(diagnostic("columns", 1L, "Column sets differ"))
  }
  if (!all(keys %in% names(expected)) || uniqueN(actual, by = keys) != nrow(actual) ||
      uniqueN(expected, by = keys) != nrow(expected)) {
    return(diagnostic("keys", 1L, "Missing or duplicate row keys"))
  }

  # Semantic row keys, not CSV ordering, define identity. Extra or missing
  # probes/terms must fail rather than disappear in an inner join.
  setorderv(actual, keys)
  setorderv(expected, keys)
  if (nrow(actual) != nrow(expected) ||
      !identical(as.data.frame(actual[, keys, with = FALSE]),
                 as.data.frame(expected[, keys, with = FALSE]))) {
    return(diagnostic("keys", 1L, "Row identities differ"))
  }
  rbindlist(lapply(names(expected), function(column) {
    observed <- actual[[column]]
    target <- expected[[column]]
    mismatches <- xor(is.na(observed), is.na(target))
    present <- !is.na(observed) & !is.na(target)
    max_error <- NA_real_
    if (is.numeric(target) && is.numeric(observed)) {
      finite <- present & is.finite(observed) & is.finite(target)
      difference <- abs(observed[finite] - target[finite])
      if (length(difference)) max_error <- max(difference)

      # No absolute floor for P/q: 1e-100 must not match 1e-15 merely because
      # both are small. Counts are exact; estimates permit floating-point noise.
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
      mismatches[present] <- as.character(observed[present]) != as.character(target[present])
    }
    detail <- if (any(mismatches)) {
      paste("First differing row:", paste(unlist(expected[which(mismatches)[1], keys, with = FALSE]), collapse = "/"))
    } else ""
    diagnostic(column, sum(mismatches), detail, max_error)
  }))
}

verify_paper_results <- function(actual_root, expected_root, report_path) {
  models <- data.frame(
    model = c("M0", "M1"), directory = c("minimal_model", "minimal_model_strand_g_richness"),
    prefix = c("quadron_minimal_stable_unstable_mean_beta", "quadron_minimal_strand_g_richness")
  )
  table_keys <- list(coefficients = c("tissue", "density_group", "term"),
                     g4_terms = c("tissue", "density_group", "term"),
                     stable_minus_unstable = c("tissue", "density_group", "contrast"))
  comparisons <- list()
  for (model_index in seq_len(nrow(models))) {
    for (table_name in names(table_keys)) {
      filename <- paste0(models$prefix[model_index], "_", table_name, ".csv")
      comparison <- compare_tables(file.path(actual_root, models$directory[model_index], filename),
                                    file.path(expected_root, filename), table_keys[[table_name]])
      comparison[, model := models$model[model_index]]
      comparisons[[length(comparisons) + 1L]] <- comparison
    }
  }
  report <- rbindlist(comparisons)
  # Nextflow provides an isolated task directory. The only write is this new
  # report, never a staged input/reference or a caller-selected arbitrary path.
  link_target <- Sys.readlink(report_path)
  is_link <- !is.na(link_target) && nzchar(link_target)
  if (basename(report_path) != report_path || file.exists(report_path) || is_link) {
    stop("Verification report must be a new filename in the task directory")
  }
  fwrite(report, report_path)
  if (!all(report$passed)) {
    print(report[passed == FALSE])
    stop("Numerical reproduction failed; expected paper results were not changed")
  }
  message("PASS: every column of all six model tables matches the paper within tolerance")
}