#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(glmnet)
  library(optparse)
})

option_list <- list(
  make_option("--prepared-dir", type = "character", default = "results/g4_methylation_age_100bp/prepared_inputs", dest = "prepared_dir"),
  make_option("--slopes-dir", type = "character", default = "results/g4_methylation_age_100bp/per_tissue", dest = "slopes_dir"),
  make_option("--g4-unmerged", type = "character", default = "results/g4_methylation_age_100bp/g4_windows/g4_motifs_100bp_unmerged.bed", dest = "g4_unmerged"),
  make_option("--chromatin", type = "character", default = "../results/v2_hallmark/chromatin_stratification/v2/v2_chromatin_probe_assignment.csv"),
  make_option("--include-chromatin", type = "character", default = "true", dest = "include_chromatin"),
  make_option("--gc-flank-bp", type = "integer", default = 1000L, dest = "gc_flank_bp"),
  make_option("--outdir", type = "character", default = "results/g4_methylation_age/cofactor_elastic_net"),
  make_option("--prefix", type = "character", default = "g4_cofactor_elastic_net"),
  make_option("--alpha", type = "double", default = 0.5),
  make_option("--nfolds", type = "integer", default = 5L),
  make_option("--seed", type = "integer", default = 42L)
)
opt <- parse_args(OptionParser(option_list = option_list))
include_chromatin <- tolower(opt$include_chromatin) %in% c("true", "t", "1", "yes", "y")
if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE) ||
    !requireNamespace("GenomicRanges", quietly = TRUE) ||
    !requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Missing BSgenome.Hsapiens.UCSC.hg38, GenomicRanges, or Biostrings for GC correction")
}

set.seed(opt$seed)
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(opt$slopes_dir, pattern = "_probe_slopes\\.csv\\.gz$", full.names = TRUE)
if (!length(files)) stop("No *_probe_slopes.csv.gz files found in ", opt$slopes_dir)
slopes <- rbindlist(lapply(files, fread), fill = TRUE)
if (!all(c("probe_id", "tissue", "slope", "stable_motifs", "unstable_motifs", "is_long_cluster", "motif_count") %in% names(slopes))) {
  stop("Probe slope files do not have the required G4/slope columns")
}

if (!"mean_beta" %in% names(slopes)) {
  beta_files <- list.files(opt$prepared_dir, pattern = "_beta\\.rds$", full.names = TRUE)
  if (!length(beta_files)) stop("No *_beta.rds files found in ", opt$prepared_dir)
  mean_beta_rows <- lapply(beta_files, function(path) {
    beta <- readRDS(path)
    tissue <- sub("^GSE[0-9]+_", "", sub("_beta\\.rds$", "", basename(path)))
    data.table(probe_id = rownames(beta), tissue = tissue, mean_beta = rowMeans(beta, na.rm = TRUE))
  })
  mean_beta <- rbindlist(mean_beta_rows)
  slopes <- merge(slopes, mean_beta, by = c("probe_id", "tissue"), all.x = TRUE)
}

if (requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE)) {
  data("Islands.UCSC", package = "IlluminaHumanMethylation450kanno.ilmn12.hg19")
  island_dt <- data.table(
    probe_id = rownames(Islands.UCSC),
    cpg_relation = as.character(Islands.UCSC$Relation_to_Island)
  )
  island_dt[, cpg_context := fifelse(cpg_relation == "Island", "Island",
                              fifelse(cpg_relation %in% c("N_Shore", "S_Shore"), "Shore",
                              fifelse(cpg_relation %in% c("N_Shelf", "S_Shelf"), "Shelf", "OpenSea")))]
  slopes <- merge(slopes, island_dt[, .(probe_id, cpg_context)], by = "probe_id", all.x = TRUE)
  data("Locations", package = "IlluminaHumanMethylation450kanno.ilmn12.hg19")
  locations_df <- as.data.frame(Locations)
  probe_strand <- data.table(probe_id = rownames(locations_df), probe_strand = as.character(locations_df$strand))
  slopes <- merge(slopes, probe_strand, by = "probe_id", all.x = TRUE)
} else {
  slopes[, cpg_context := NA_character_]
  slopes[, probe_strand := NA_character_]
}

if (include_chromatin && file.exists(opt$chromatin)) {
  chrom <- fread(opt$chromatin)
  probe_col <- intersect(c("probe", "probe_id"), names(chrom))[1]
  if (!is.na(probe_col) && "state" %in% names(chrom)) {
    chrom <- unique(chrom[, .(probe_id = as.character(get(probe_col)), chrom_state = as.character(state))], by = "probe_id")
    slopes <- merge(slopes, chrom, by = "probe_id", all.x = TRUE)
  } else {
    slopes[, chrom_state := NA_character_]
  }
} else {
  slopes[, chrom_state := NA_character_]
}

slopes[, `:=`(
  has_stable_g4 = as.integer(stable_motifs > 0),
  has_unstable_g4 = as.integer(unstable_motifs > 0),
  is_long_cluster = as.integer(as.logical(is_long_cluster)),
  log1p_motif_count = log1p(motif_count),
  cpg_context = fifelse(is.na(cpg_context), "Unknown_CpG_context", cpg_context),
  chrom_state = fifelse(is.na(chrom_state), "Unknown_chromatin", chrom_state)
)]

compute_sequence_features <- function(probe_dt, flank_bp) {
  genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  probes <- unique(probe_dt[, c("probe_id", "chr", "start", "end"), with = FALSE], by = "probe_id")
  probes <- probes[probes[["chr"]] %in% names(GenomeInfoDb::seqlengths(genome)) & is.finite(probes[["start"]]) & is.finite(probes[["end"]])]
  probes[["center"]] <- floor((probes[["start"]] + probes[["end"]]) / 2)
  probes[["gc_start"]] <- pmax(1L, probes[["center"]] - flank_bp)
  chr_lengths <- GenomeInfoDb::seqlengths(genome)
  probes[["gc_end"]] <- pmin(chr_lengths[probes[["chr"]]], probes[["center"]] + flank_bp)
  gr <- GenomicRanges::GRanges(probes$chr, IRanges::IRanges(probes$gc_start, probes$gc_end))
  seqs <- BSgenome::getSeq(genome, gr)
  counts <- Biostrings::letterFrequency(seqs, letters = c("G", "C"), as.prob = FALSE)
  probes[["gc_fraction"]] <- rowSums(counts) / Biostrings::width(seqs)
  probes[["g_fraction"]] <- counts[, "G"] / Biostrings::width(seqs)
  probes[, c("probe_id", "gc_fraction", "g_fraction"), with = FALSE]
}
sequence_dt <- compute_sequence_features(slopes, opt$gc_flank_bp)
slopes <- merge(slopes, sequence_dt, by = "probe_id", all.x = TRUE)

compute_g4_strand_features <- function(probe_dt, g4_unmerged) {
  probes <- unique(probe_dt[, c("probe_id", "chr", "start", "end", "probe_strand"), with = FALSE], by = "probe_id")
  probes <- probes[probes[["probe_strand"]] %in% c("+", "-")]
  probes[["probe_index"]] <- seq_len(nrow(probes))
  g4 <- fread(g4_unmerged, header = FALSE)
  setnames(g4, c("chr", "start", "end", "motif_id", "g4_class", "g4_strand", "g4_score", "core_start", "core_end"))
  g4 <- g4[g4[["g4_strand"]] %in% c("+", "-") & g4[["g4_class"]] %in% c("stable", "unstable")]
  setkeyv(probes, c("chr", "start", "end"))
  setkeyv(g4, c("chr", "start", "end"))
  hits <- foverlaps(probes, g4, nomatch = 0L)
  out <- data.table(probe_id = probes$probe_id)
  if (nrow(hits)) {
    hits[["relation"]] <- fifelse(hits[["probe_strand"]] == hits[["g4_strand"]], "same", "opposite")
    counts <- hits[, list(N = .N), by = c("probe_id", "g4_class", "relation")]
    counts[["feature"]] <- paste(counts[["g4_class"]], counts[["relation"]], "strand", sep = "_")
    wide <- dcast(counts, probe_id ~ feature, value.var = "N", fill = 0)
    out <- merge(out, wide, by = "probe_id", all.x = TRUE)
  }
  for (column in c("stable_same_strand", "stable_opposite_strand", "unstable_same_strand", "unstable_opposite_strand")) {
    if (!column %in% names(out)) out[, (column) := 0L]
    out[is.na(get(column)), (column) := 0L]
  }
  out[["log1p_stable_same_strand"]] <- log1p(out[["stable_same_strand"]])
  out[["log1p_stable_opposite_strand"]] <- log1p(out[["stable_opposite_strand"]])
  out[["log1p_unstable_same_strand"]] <- log1p(out[["unstable_same_strand"]])
  out[["log1p_unstable_opposite_strand"]] <- log1p(out[["unstable_opposite_strand"]])
  out
}
strand_dt <- compute_g4_strand_features(slopes, opt$g4_unmerged)
slopes <- merge(slopes, strand_dt, by = "probe_id", all.x = TRUE)
strand_cols <- c("stable_same_strand", "stable_opposite_strand", "unstable_same_strand", "unstable_opposite_strand",
                 "log1p_stable_same_strand", "log1p_stable_opposite_strand", "log1p_unstable_same_strand", "log1p_unstable_opposite_strand")
for (column in strand_cols) slopes[is.na(get(column)), (column) := 0]

cofactor_inventory <- rbindlist(list(
  data.table(cofactor = "has_stable_g4", type = "binary", description = "Probe lies in a merged +/-1 kb G4 window containing at least one stable G4Hunter motif"),
  data.table(cofactor = "has_unstable_g4", type = "binary", description = "Probe lies in a merged +/-1 kb G4 window containing at least one unstable G4Hunter motif"),
  data.table(cofactor = "is_long_cluster", type = "binary", description = "Probe lies in an unusually long merged dense G4 cluster"),
  data.table(cofactor = "log1p_motif_count", type = "numeric", description = "Log-transformed number of G4 motifs contributing to the merged window"),
  data.table(cofactor = "log1p_stable_same_strand", type = "numeric", description = "Log-transformed count of stable +/-1 kb G4 motifs on the same strand as the Illumina probe annotation"),
  data.table(cofactor = "log1p_stable_opposite_strand", type = "numeric", description = "Log-transformed count of stable +/-1 kb G4 motifs on the opposite strand from the Illumina probe annotation"),
  data.table(cofactor = "log1p_unstable_same_strand", type = "numeric", description = "Log-transformed count of unstable +/-1 kb G4 motifs on the same strand as the Illumina probe annotation"),
  data.table(cofactor = "log1p_unstable_opposite_strand", type = "numeric", description = "Log-transformed count of unstable +/-1 kb G4 motifs on the opposite strand from the Illumina probe annotation"),
  data.table(cofactor = "mean_beta", type = "numeric", description = "Mean methylation beta of the CpG across samples in the tissue"),
  data.table(cofactor = "gc_fraction", type = "numeric", description = sprintf("Local hg38 GC fraction in +/- %d bp around the CpG probe", opt$gc_flank_bp)),
  data.table(cofactor = "g_fraction", type = "numeric", description = sprintf("Local hg38 guanine fraction in +/- %d bp around the CpG probe", opt$gc_flank_bp)),
  data.table(cofactor = "cpg_context", type = "categorical", description = "Illumina 450K CpG island relation collapsed to Island/Shore/Shelf/OpenSea")
))
if (include_chromatin) {
  cofactor_inventory <- rbind(
    cofactor_inventory,
    data.table(cofactor = "chrom_state", type = "categorical", description = "Probe-level ChromHMM state from the local CCRRBL-20 chromatin assignment")
  )
}
fwrite(cofactor_inventory, file.path(opt$outdir, sprintf("%s_cofactor_inventory.csv", opt$prefix)))

strand_predictors <- c("log1p_stable_same_strand", "log1p_stable_opposite_strand", "log1p_unstable_same_strand", "log1p_unstable_opposite_strand")
predictor_terms <- c("has_stable_g4", "has_unstable_g4", "is_long_cluster", "log1p_motif_count", strand_predictors, "mean_beta", "gc_fraction", "g_fraction", "cpg_context")
if (include_chromatin) predictor_terms <- c(predictor_terms, "chrom_state")

fit_one_tissue <- function(tissue_name) {
  dt <- slopes[slopes[["tissue"]] == tissue_name & is.finite(slopes[["slope"]]) & is.finite(slopes[["mean_beta"]])]
  complete_cols <- predictor_terms
  dt <- dt[complete.cases(dt[, complete_cols, with = FALSE])]
  if (nrow(dt) < 100) stop("Too few complete probes for tissue ", tissue_name)

  formula <- as.formula(paste("~ (", paste(predictor_terms, collapse = " + "), ")^2"))
  x <- model.matrix(formula, data = dt)[, -1, drop = FALSE]
  keep <- apply(x, 2, function(col) is.finite(sd(col)) && sd(col) > 0)
  x <- x[, keep, drop = FALSE]
  y <- dt$slope

  cvfit <- cv.glmnet(x, y, alpha = opt$alpha, nfolds = opt$nfolds, standardize = TRUE)
  pred_min <- as.numeric(predict(cvfit, newx = x, s = "lambda.min"))
  pred_1se <- as.numeric(predict(cvfit, newx = x, s = "lambda.1se"))
  perf <- data.table(
    tissue = tissue_name,
    n_probes = nrow(dt),
    n_predictors = ncol(x),
    alpha = opt$alpha,
    lambda_min = cvfit$lambda.min,
    lambda_1se = cvfit$lambda.1se,
    train_r2_lambda_min = 1 - sum((y - pred_min)^2) / sum((y - mean(y))^2),
    train_r2_lambda_1se = 1 - sum((y - pred_1se)^2) / sum((y - mean(y))^2),
    cv_mse_min = min(cvfit$cvm),
    null_mse = var(y)
  )

  extract_coef <- function(s_label) {
    cf <- as.matrix(coef(cvfit, s = s_label))
    out <- data.table(term = rownames(cf), coefficient = as.numeric(cf[, 1]))
    out[, `:=`(tissue = tissue_name, lambda = s_label)]
    out[out[["coefficient"]] != 0]
  }
  coefs <- rbindlist(list(extract_coef("lambda.min"), extract_coef("lambda.1se")))
  list(perf = perf, coefs = coefs)
}

fits <- lapply(sort(unique(slopes$tissue)), fit_one_tissue)
performance <- rbindlist(lapply(fits, `[[`, "perf"), fill = TRUE)
coefficients <- rbindlist(lapply(fits, `[[`, "coefs"), fill = TRUE)
coefficients[, abs_coefficient := abs(coefficient)]
coefficients[, term_class := fifelse(grepl(":", term), "interaction", "main_effect")]

fwrite(performance, file.path(opt$outdir, sprintf("%s_model_performance.csv", opt$prefix)))
fwrite(coefficients, file.path(opt$outdir, sprintf("%s_coefficients.csv", opt$prefix)))

plot_dt <- coefficients[lambda == "lambda.1se" & term != "(Intercept)"]
top_terms <- plot_dt[, .(max_abs = max(abs_coefficient)), by = term][order(-max_abs)][seq_len(min(30, .N)), term]
plot_dt <- plot_dt[term %in% top_terms]
plot_dt[, term := factor(term, levels = rev(top_terms))]
plot_dt[, tissue := factor(tissue, levels = c("adipose", "liver", "muscle"))]
max_abs <- max(abs(plot_dt$coefficient), na.rm = TRUE)

heatmap <- ggplot(plot_dt, aes(x = tissue, y = term, fill = coefficient)) +
  geom_tile(color = "white", linewidth = 0.45) +
  scale_fill_gradient2(low = "#2C7BB6", mid = "#F7F7F7", high = "#D7191C",
                       midpoint = 0, limits = c(-max_abs, max_abs), name = "Elastic-net\ncoefficient") +
  labs(
    title = if (include_chromatin) "Wide cofactor elastic-net model of signed CpG methylation aging" else "G4/CpG-context elastic-net model of signed CpG methylation aging",
    subtitle = "Top non-zero lambda.1se coefficients; model includes all pairwise interactions",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), axis.text.y = element_text(size = 8),
        axis.text.x = element_text(face = "bold"), plot.title = element_text(face = "bold"))

ggsave(file.path(opt$outdir, sprintf("%s_top_coefficients_heatmap.png", opt$prefix)), heatmap, width = 10, height = 8, dpi = 300)
ggsave(file.path(opt$outdir, sprintf("%s_top_coefficients_heatmap.pdf", opt$prefix)), heatmap, width = 10, height = 8)

message("Wrote outputs to: ", opt$outdir)
