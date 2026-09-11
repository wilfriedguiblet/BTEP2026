#!/usr/bin/env Rscript

# PURPOSE: estimate signed methylation-versus-age slopes for each CpG and
# attach merged-neighborhood annotations. The paper uses these slope tables as
# INPUTS to compare_g4_effects_adjusted.R, not the exploratory regressions at
# the end of this script. The curated runner uses --slopes-only to stop early.
# For exact replay, retain the prepared beta matrix, sample metadata, hg38
# manifest, merged Quadron checkpoint, requested covariates, and their hashes.

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

option_list <- list(
  make_option("--beta", type = "character"),
  make_option("--metadata", type = "character"),
  make_option("--manifest", type = "character"),
  make_option("--g4-windows", type = "character", dest = "g4_windows"),
  make_option("--age-col", type = "character", default = "age", dest = "age_col"),
  make_option("--sample-col", type = "character", default = "sample_id", dest = "sample_col"),
  make_option("--covariates", type = "character", default = ""),
  make_option("--min-samples", type = "integer", default = 2L, dest = "min_samples"),
  make_option("--cohort", type = "character"),
  make_option("--tissue", type = "character"),
  make_option("--prefix", type = "character"),
  make_option("--slopes-only", action = "store_true", default = FALSE, dest = "slopes_only"),
  make_option("--outdir", type = "character", default = ".")
)
opt <- parse_args(OptionParser(option_list = option_list))
required <- c("beta", "metadata", "manifest", "g4_windows", "cohort", "tissue", "prefix")
missing <- required[vapply(required, function(x) is.null(opt[[x]]), logical(1))]
if (length(missing)) stop("Missing required options: ", paste(missing, collapse = ", "))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

read_table_or_rds <- function(path) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) return(readRDS(path))
  x <- fread(path, data.table = FALSE)
  first <- names(x)[1]
  if (!grepl("^cg", names(x)[1]) && any(grepl("^cg", x[[first]]))) {
    rownames(x) <- x[[first]]
    x[[first]] <- NULL
  }
  as.matrix(x)
}

beta <- read_table_or_rds(opt$beta)
storage.mode(beta) <- "numeric"
metadata <- if (grepl("\\.rds$", opt$metadata, ignore.case = TRUE)) as.data.table(readRDS(opt$metadata)) else fread(opt$metadata)
if (!opt$age_col %in% names(metadata)) stop("Age column not found: ", opt$age_col)
if (!opt$sample_col %in% names(metadata)) stop("Sample column not found: ", opt$sample_col)
metadata[, age_for_model := suppressWarnings(as.numeric(get(opt$age_col)))]

covariate_text <- if (is.null(opt$covariates)) "" else opt$covariates
# Legacy behavior silently drops requested covariates absent from metadata.
# The curated paper runner therefore checks age/sex/BMI columns before calling
# this script. Missing values are a different issue and trigger complete cases.
covariates <- unlist(strsplit(covariate_text, "\\|"))
covariates <- covariates[nzchar(covariates) & covariates %in% names(metadata)]
common_samples <- intersect(colnames(beta), metadata[[opt$sample_col]])
# Sample order MUST follow the same ID join in beta and metadata. Subject IDs
# may recur across tissues; this code fits tissues separately, not a donor-
# aware joint model, and does not repair inconsistent subject labels.
metadata <- metadata[match(common_samples, metadata[[opt$sample_col]])]
beta <- beta[, common_samples, drop = FALSE]

model_columns <- c("age_for_model", covariates)
keep_samples <- complete.cases(metadata[, ..model_columns])
metadata <- metadata[keep_samples]
beta <- beta[, keep_samples, drop = FALSE]
if (ncol(beta) < opt$min_samples) stop("Too few complete samples after metadata matching")

manifest <- fread(opt$manifest)
# Prefer explicitly named hg38 columns. The manifest defines the tested probe
# universe; substituting all 450K probes or rebuilding a gene-only manifest can
# change counts. Verify unique probe IDs and the actual chromosome/coordinate
# content upstream, rather than trusting a filename or hg38 column label.
chr_col <- intersect(c("chr_hg38", "chr", "chrom", "chromosome"), names(manifest))[1]
start_col <- intersect(c("start_hg38", "start", "pos", "position"), names(manifest))[1]
end_col <- intersect(c("end_hg38", "end"), names(manifest))[1]
probe_col <- intersect(c("probe_id", "ID", "Name", "IlmnID"), names(manifest))[1]
if (any(is.na(c(chr_col, start_col, end_col, probe_col)))) {
  stop("Manifest must contain probe id and hg38 chr/start/end columns")
}
probe_dt <- manifest[, .(
  probe_id = as.character(get(probe_col)),
  chr = as.character(get(chr_col)),
  start = as.integer(get(start_col)),
  end = as.integer(get(end_col))
)]
probe_dt <- probe_dt[probe_id %in% rownames(beta) & grepl("^chr([0-9]+|X|Y)$", chr)]
beta <- beta[probe_dt$probe_id, , drop = FALSE]
# This mean is taken AFTER sample matching/metadata exclusions, omitting NA
# betas within each probe. It is the same-sample mean, not a young or baseline
# reference. Preserve it in the slope checkpoint for the second-stage models.
probe_dt[, mean_beta := rowMeans(beta, na.rm = TRUE)]

design_dt <- copy(metadata)
if (length(covariates)) {
  for (covar in covariates) {
    if (!is.numeric(design_dt[[covar]])) design_dt[, (covar) := as.factor(get(covar))]
  }
}
formula_text <- paste("~ age_for_model", if (length(covariates)) paste("+", paste(covariates, collapse = " + ")) else "")
design <- model.matrix(as.formula(formula_text), data = design_dt)
age_column <- which(colnames(design) == "age_for_model")

fit_probe <- function(y) {
  # Fit beta on age + the requested covariates using the finite observations
  # for this probe. Requiring one more observation than design columns leaves
  # residual degrees of freedom (normally >=5 specimens for age/sex/BMI).
  # Slopes retain beta/year units. This is cross-sectional association, NOT a
  # within-person aging rate. There is no weighting or empirical-Bayes step.
  ok <- is.finite(y) & complete.cases(design)
  if (sum(ok) < max(opt$min_samples, ncol(design) + 1L)) return(c(NA_real_, NA_real_, NA_real_))
  fit <- lm.fit(design[ok, , drop = FALSE], y[ok])
  rdf <- fit$df.residual
  if (rdf <= 0) return(c(NA_real_, NA_real_, NA_real_))
  # Conventional residual-variance/QR SE. A rank-deficient or pivoted design
  # needs scrutiny: the historical indexing is preserved here. The curated
  # runner checks the specimen-level design rank; probe-level missingness may
  # still change estimability. Later models do not propagate these slope SEs.
  xtx_inv <- tryCatch(chol2inv(fit$qr$qr[seq_len(fit$rank), seq_len(fit$rank), drop = FALSE]), error = function(e) NULL)
  if (is.null(xtx_inv) || age_column > nrow(xtx_inv)) return(c(NA_real_, NA_real_, NA_real_))
  sigma2 <- sum(fit$residuals^2, na.rm = TRUE) / rdf
  se <- sqrt(sigma2 * xtx_inv[age_column, age_column])
  slope <- fit$coefficients[age_column]
  t_stat <- slope / se
  p_value <- 2 * pt(-abs(t_stat), df = rdf)
  c(slope, se, p_value)
}

coef_mat <- t(apply(beta, 1, fit_probe))
probe_dt[, `:=`(
  slope = coef_mat[, 1],
  se = coef_mat[, 2],
  p_value = coef_mat[, 3],
  abs_slope = abs(coef_mat[, 1])
)]

windows <- fread(opt$g4_windows, header = FALSE)
window_names <- c("chr", "start", "end", "window_id", "class_label", "motif_count", "stable_motifs", "unstable_motifs",
                  "width", "is_long_cluster", "base_class")
setnames(windows, window_names[seq_len(ncol(windows))])
if (!"is_long_cluster" %in% names(windows)) windows[, is_long_cluster := FALSE]
if (!"base_class" %in% names(windows)) windows[, base_class := class_label]
probe_intervals <- copy(probe_dt)[, probe_index := .I]
setkey(probe_intervals, chr, start, end)
setkey(windows, chr, start, end)
# foverlaps uses inclusive endpoints. This preserves the historical coordinate
# implementation and is NOT an implicit conversion from standard BED. A probe
# inherits counts from its whole merged neighborhood, not necessarily from a
# direct overlap with each contributing motif. The first overlap is retained.
hits <- foverlaps(probe_intervals, windows, nomatch = 0L)
# in_g4_1kb is a legacy COLUMN NAME, not an instruction to use 1-kb G4 flanks.
# The supplied merged 100-bp Quadron checkpoint determines the paper's geometry.
annotation <- data.table(probe_index = seq_len(nrow(probe_dt)), in_g4_1kb = FALSE,
                         window_id = NA_character_, g4_class = "outside",
                         motif_count = 0L, stable_motifs = 0L, unstable_motifs = 0L,
                         is_long_cluster = FALSE, base_class = "outside")
if (nrow(hits)) {
  hit_map <- hits[, .SD[1], by = probe_index]
  annotation[hit_map$probe_index, `:=`(
    in_g4_1kb = TRUE,
    window_id = hit_map$window_id,
    g4_class = hit_map$class_label,
    motif_count = as.integer(hit_map$motif_count),
    stable_motifs = as.integer(hit_map$stable_motifs),
    unstable_motifs = as.integer(hit_map$unstable_motifs),
    is_long_cluster = as.logical(hit_map$is_long_cluster),
    base_class = hit_map$base_class
  )]
}
probe_dt <- cbind(probe_dt, annotation[, -"probe_index"])
probe_dt[, `:=`(
  cohort = opt$cohort,
  tissue = opt$tissue,
  has_stable_g4 = stable_motifs > 0,
  has_unstable_g4 = unstable_motifs > 0,
  log1p_motif_count = log1p(motif_count)
)]

# Reproduction of the paper needs ONLY the slope/annotation checkpoint. The
# opt-in exit avoids emitting older Wilcoxon and architecture-heavy regressions
# without changing the historical default for other callers of this script.
if (isTRUE(opt$slopes_only)) {
  fwrite(probe_dt, file.path(opt$outdir, sprintf("%s_probe_slopes.csv.gz", opt$prefix)))
  quit(save = "no", status = 0L)
}

# Everything below is an older exploratory output family, NOT manuscript
# evidence. It remains for backwards compatibility but the paper runner skips it.
summary_dt <- probe_dt[is.finite(slope), .(
  n_probes = .N,
  mean_slope = mean(slope),
  median_slope = median(slope),
  mean_abs_slope = mean(abs_slope),
  median_abs_slope = median(abs_slope),
  frac_hypermethylating = mean(slope > 0),
  frac_nominal_age_associated = mean(p_value < 0.05, na.rm = TRUE)
), by = .(cohort, tissue, group = fifelse(in_g4_1kb, g4_class, "outside"))]

inside <- probe_dt[in_g4_1kb == TRUE & is.finite(abs_slope), abs_slope]
outside <- probe_dt[in_g4_1kb == FALSE & is.finite(abs_slope), abs_slope]
tests <- data.table(
  cohort = opt$cohort,
  tissue = opt$tissue,
  comparison = c("g4_1kb_vs_outside_abs_slope", "g4_1kb_vs_outside_slope"),
  statistic = NA_real_,
  p_value = NA_real_,
  alternative = c("greater", "two.sided"),
  n_inside = length(inside),
  n_outside = length(outside),
  formula = formula_text
)
if (length(inside) >= 2 && length(outside) >= 2) {
  wt_abs <- wilcox.test(inside, outside, alternative = "greater")
  wt_signed <- wilcox.test(probe_dt[in_g4_1kb == TRUE, slope], probe_dt[in_g4_1kb == FALSE, slope])
  tests[comparison == "g4_1kb_vs_outside_abs_slope", `:=`(statistic = unname(wt_abs$statistic), p_value = wt_abs$p.value)]
  tests[comparison == "g4_1kb_vs_outside_slope", `:=`(statistic = unname(wt_signed$statistic), p_value = wt_signed$p.value)]
}

fit_annotation_model <- function(response, response_name) {
  model_dt <- probe_dt[is.finite(probe_dt[[response]]), list(
    y = probe_dt[[response]][is.finite(probe_dt[[response]])],
    has_stable_g4 = get("has_stable_g4"),
    has_unstable_g4 = get("has_unstable_g4"),
    is_long_cluster = get("is_long_cluster"),
    log1p_motif_count = get("log1p_motif_count")
  )]
  if (nrow(model_dt) < 10) return(data.table())
  fit <- lm(y ~ has_stable_g4 * has_unstable_g4 + is_long_cluster + log1p_motif_count, data = model_dt)
  coefs <- as.data.table(summary(fit)$coefficients, keep.rownames = "term")
  setnames(coefs, c("term", "estimate", "std_error", "statistic", "p_value"))
  coefs[, `:=`(
    cohort = opt$cohort,
    tissue = opt$tissue,
    response = response_name,
    n_probes = nrow(model_dt),
    model = "response ~ has_stable_g4 * has_unstable_g4 + is_long_cluster + log1p_motif_count"
  )]
  coefs[, c("cohort", "tissue", "response", "term", "estimate", "std_error", "statistic", "p_value", "n_probes", "model"), with = FALSE]
}

regression_effects <- rbindlist(list(
  fit_annotation_model("mean_beta", "baseline_mean_methylation_beta"),
  fit_annotation_model("slope", "age_slope_delta_beta_per_year"),
  fit_annotation_model("abs_slope", "absolute_age_slope_delta_beta_per_year")
), fill = TRUE)

fwrite(probe_dt, file.path(opt$outdir, sprintf("%s_probe_slopes.csv.gz", opt$prefix)))
fwrite(summary_dt, file.path(opt$outdir, sprintf("%s_g4_age_summary.csv", opt$prefix)))
fwrite(tests, file.path(opt$outdir, sprintf("%s_g4_age_tests.csv", opt$prefix)))
fwrite(regression_effects, file.path(opt$outdir, sprintf("%s_g4_regression_effects.csv", opt$prefix)))
