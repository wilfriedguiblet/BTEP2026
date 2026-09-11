# Fast, offline tests of the reproduction gate. These deliberately mutate
# synthetic results to confirm that the checker can fail, not just pass.
script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_argument))
source(file.path(dirname(dirname(script_path)), "compare_results.R"))
scratch <- tempfile("quadron_comparison_test_")
dir.create(scratch)

run_tests <- function() {
  on.exit(unlink(scratch, recursive = TRUE), add = TRUE)
  expected <- data.table::data.table(
    tissue = c("adipose", "liver", "muscle"), term = "has_stable_g4TRUE",
    estimate = c(1e-5, 2e-5, 3e-5), std_error = rep(1e-6, 3),
    p_value = c(0.2, 1e-60, 0.004), q_value = c(0.2, 3e-60, 0.006),
    n_probes = rep(100L, 3), model = "slope ~ stable + unstable + mean_beta"
  )
  expected_path <- file.path(scratch, "expected.csv")
  actual_path <- file.path(scratch, "actual.csv")
  data.table::fwrite(expected, expected_path)
  check <- function(candidate) {
    data.table::fwrite(candidate, actual_path)
    all(compare_tables(actual_path, expected_path, c("tissue", "term"))$passed)
  }
  stopifnot(check(expected), check(expected[3:1]))
  perturbed <- data.table::copy(expected)
  perturbed[, estimate := estimate + 1e-14]
  stopifnot(check(perturbed))
  perturbed <- data.table::copy(expected)
  perturbed[2, p_value := 1e-15]
  stopifnot(!check(perturbed))
  perturbed <- data.table::copy(expected)
  perturbed[1, n_probes := 101L]
  stopifnot(!check(perturbed))
  perturbed <- data.table::copy(expected)
  perturbed[1, model := "slope ~ wrong_model"]
  stopifnot(!check(perturbed), !check(expected[-1]), !check(rbind(expected, expected[1])))
  perturbed <- data.table::copy(expected)
  perturbed[1, estimate := NA_real_]
  stopifnot(!check(perturbed))
  message("PASS: row order and harmless rounding accepted; tiny-P changes, formula/count changes, missing/duplicate rows, and missing estimates rejected.")
}

run_tests()