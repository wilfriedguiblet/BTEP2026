#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

option_list <- list(
  make_option("--slope-file", type = "character", default = "results/g4_methylation_age/per_tissue/GSE61257_adipose_probe_slopes.csv.gz", dest = "slope_file"),
  make_option("--outdir", type = "character", default = "results/g4_methylation_age/context_g4_association"),
  make_option("--prefix", type = "character", default = "cpg_context_g4_association")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
if (!requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE)) {
  stop("Missing IlluminaHumanMethylation450kanno.ilmn12.hg19")
}

slopes <- fread(opt$slope_file)
data("Islands.UCSC", package = "IlluminaHumanMethylation450kanno.ilmn12.hg19")
island <- data.table(
  probe_id = rownames(Islands.UCSC),
  relation_to_island = as.character(Islands.UCSC$Relation_to_Island)
)
island[, cpg_context := fifelse(
  relation_to_island == "Island", "Island",
  fifelse(relation_to_island %in% c("N_Shore", "S_Shore"), "Shore",
          fifelse(relation_to_island %in% c("N_Shelf", "S_Shelf"), "Shelf", "OpenSea"))
)]

d <- merge(slopes, island[, .(probe_id, cpg_context)], by = "probe_id", all.x = TRUE)
d <- d[!is.na(d[["cpg_context"]])]
d[, `:=`(
  has_stable_g4 = stable_motifs > 0,
  has_unstable_g4 = unstable_motifs > 0
)]

cramers_v <- function(x, y) {
  tab <- table(x, y)
  chi <- suppressWarnings(chisq.test(tab, correct = FALSE)$statistic)
  n <- sum(tab)
  k <- min(nrow(tab), ncol(tab))
  as.numeric(sqrt(chi / (n * (k - 1))))
}

summarize_predictor <- function(predictor) {
  overall <- mean(d[[predictor]])
  contexts <- sort(unique(d$cpg_context))
  rows <- lapply(contexts, function(context) {
    in_context <- d$cpg_context == context
    a <- sum(in_context & d[[predictor]])
    b <- sum(in_context & !d[[predictor]])
    c <- sum(!in_context & d[[predictor]])
    e <- sum(!in_context & !d[[predictor]])
    odds_ratio <- ((a + 0.5) * (e + 0.5)) / ((b + 0.5) * (c + 0.5))
    data.table(
      predictor = predictor,
      cpg_context = context,
      n_probes = sum(in_context),
      n_with_g4 = a,
      g4_fraction = a / sum(in_context),
      overall_g4_fraction = overall,
      enrichment_vs_overall = (a / sum(in_context)) / overall,
      odds_ratio_context_vs_rest = odds_ratio
    )
  })
  rbindlist(rows)
}

summary_dt <- rbindlist(list(
  summarize_predictor("has_stable_g4"),
  summarize_predictor("has_unstable_g4")
))
association_dt <- data.table(
  predictor = c("has_stable_g4", "has_unstable_g4"),
  cramers_v = c(
    cramers_v(d$cpg_context, d$has_stable_g4),
    cramers_v(d$cpg_context, d$has_unstable_g4)
  ),
  n_probes = nrow(d)
)
counts_dt <- d[, .N, by = c("cpg_context", "has_stable_g4", "has_unstable_g4")]

fwrite(summary_dt, file.path(opt$outdir, sprintf("%s_by_context.csv", opt$prefix)))
fwrite(association_dt, file.path(opt$outdir, sprintf("%s_cramers_v.csv", opt$prefix)))
fwrite(counts_dt, file.path(opt$outdir, sprintf("%s_counts.csv", opt$prefix)))

cat("Cramer's V\n")
print(association_dt)
cat("\nFractions by CpG context\n")
print(summary_dt[order(predictor, -g4_fraction)])
