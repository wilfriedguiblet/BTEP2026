#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(Biostrings)
  library(optparse)
})

option_list <- list(
  make_option("--slopes", type = "character", default = "results/g4_methylation_age_quadron_100bp/per_tissue/GSE61257_adipose_probe_slopes.csv.gz"),
  make_option("--structure", type = "character", default = "results/g4_methylation_age_quadron_100bp/sensitivity_annotations/probe_structure_annotations.csv.gz"),
  make_option("--chromatin", type = "character", default = "../results/v2_hallmark/chromatin_stratification/v2/v2_chromatin_probe_assignment.csv"),
  make_option("--outdir", type = "character", default = "results/g4_methylation_age_quadron_100bp/confounder_overlap"),
  make_option("--flank", type = "integer", default = 1000L)
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

slopes <- fread(opt$slopes)
structure <- fread(opt$structure)
d <- merge(slopes, structure, by = "probe_id", all.x = TRUE)
if (all(c("chr.x", "start.x", "end.x") %in% names(d))) {
  setnames(d, c("chr.x", "start.x", "end.x"), c("chr", "start", "end"))
}
d[, `:=`(stable_quadron = stable_motifs > 0, unstable_quadron = unstable_motifs > 0)]

data("Islands.UCSC", package = "IlluminaHumanMethylation450kanno.ilmn12.hg19")
islands <- data.table(probe_id = rownames(Islands.UCSC), relation = as.character(Islands.UCSC$Relation_to_Island))
islands[, cpg_context := fifelse(relation == "Island", "Island",
                           fifelse(relation %in% c("N_Shore", "S_Shore"), "Shore",
                           fifelse(relation %in% c("N_Shelf", "S_Shelf"), "Shelf", "OpenSea")))]
d <- merge(d, islands[, .(probe_id, cpg_context)], by = "probe_id", all.x = TRUE)

chrom <- fread(opt$chromatin)
chrom <- unique(chrom[, .(probe_id = as.character(probe), chrom_state = as.character(state))], by = "probe_id")
d <- merge(d, chrom, by = "probe_id", all.x = TRUE)

genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
coords <- unique(d[, c("probe_id", "chr", "start", "end"), with = FALSE], by = "probe_id")
coords <- coords[coords[["chr"]] %in% names(GenomeInfoDb::seqlengths(genome))]
coords[, center := floor((start + end) / 2)]
coords[, seq_start := pmax(1L, center - opt$flank)]
coords[, seq_end := pmin(as.integer(GenomeInfoDb::seqlengths(genome)[chr]), center + opt$flank)]
seqs <- BSgenome::getSeq(genome, GRanges(coords$chr, IRanges(coords$seq_start, coords$seq_end)))
counts <- Biostrings::letterFrequency(seqs, letters = c("G", "C"), as.prob = FALSE)
coords[, `:=`(gc_fraction = rowSums(counts) / Biostrings::width(seqs), g_fraction = counts[, "G"] / Biostrings::width(seqs))]
d <- merge(d, coords[, .(probe_id, gc_fraction, g_fraction)], by = "probe_id", all.x = TRUE)

cramers_v <- function(x, y) {
  tab <- table(x, y, useNA = "no")
  if (min(dim(tab)) < 2) return(NA_real_)
  chi <- suppressWarnings(chisq.test(tab, correct = FALSE)$statistic)
  sqrt(as.numeric(chi) / (sum(tab) * (min(dim(tab)) - 1)))
}

binary_features <- c("rloop_overlap", "segmental_duplication")
categorical_features <- c("repeat_class", "cpg_context", "chrom_state")
numeric_features <- c("mean_beta", "gc_fraction", "g_fraction", "log1p_motif_count")
rows <- list()
index <- 0L
for (target in c("stable_quadron", "unstable_quadron")) {
  for (feature in binary_features) {
    index <- index + 1L
    tab <- table(d[[target]], d[[feature]], useNA = "no")
    rows[[index]] <- data.table(target = target, feature = feature, feature_type = "binary",
                                association = cramers_v(d[[target]], d[[feature]]),
                                positive_fraction_target = mean(d[d[[target]] == TRUE, get(feature)], na.rm = TRUE),
                                positive_fraction_non_target = mean(d[d[[target]] == FALSE, get(feature)], na.rm = TRUE))
  }
  for (feature in categorical_features) {
    index <- index + 1L
    rows[[index]] <- data.table(target = target, feature = feature, feature_type = "categorical",
                                association = cramers_v(d[[target]], d[[feature]]),
                                positive_fraction_target = NA_real_, positive_fraction_non_target = NA_real_)
  }
  for (feature in numeric_features) {
    index <- index + 1L
    in_target <- d[d[[target]] == TRUE, get(feature)]
    outside_target <- d[d[[target]] == FALSE, get(feature)]
    pooled_sd <- sqrt((var(in_target, na.rm = TRUE) + var(outside_target, na.rm = TRUE)) / 2)
    rows[[index]] <- data.table(target = target, feature = feature, feature_type = "numeric",
                                association = (mean(in_target, na.rm = TRUE) - mean(outside_target, na.rm = TRUE)) / pooled_sd,
                                positive_fraction_target = mean(in_target, na.rm = TRUE),
                                positive_fraction_non_target = mean(outside_target, na.rm = TRUE))
  }
}

summary_dt <- rbindlist(rows)
summary_dt[, abs_association := abs(association)]
setorder(summary_dt, target, -abs_association)
fwrite(summary_dt, file.path(opt$outdir, "quadron_confounder_overlap_summary.csv"))

detail <- rbindlist(lapply(c("stable_quadron", "unstable_quadron"), function(target) {
  rbindlist(lapply(c("repeat_class", "cpg_context", "chrom_state"), function(feature) {
    d[, .(n = .N, target_fraction = mean(get(target))), by = feature][, `:=`(target = target, feature_type = feature)]
  }), fill = TRUE)
}), fill = TRUE)
fwrite(detail, file.path(opt$outdir, "quadron_confounder_overlap_detail.csv"))

plot_dt <- summary_dt[feature_type %in% c("binary", "categorical")]
plot_dt[, feature := factor(feature, levels = unique(feature[order(abs(association))]))]
plot <- ggplot(plot_dt, aes(x = association, y = feature, color = target)) +
  geom_point(size = 3, position = position_dodge(width = 0.35)) +
  geom_vline(xintercept = 0, color = "grey55") +
  labs(title = "Quadron annotation overlap with potential confounders",
       subtitle = "Cramer's V for binary/categorical features; higher values indicate stronger association",
       x = "Cramer's V", y = NULL, color = NULL) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))
ggsave(file.path(opt$outdir, "quadron_confounder_overlap.png"), plot, width = 8, height = 4.8, dpi = 300)
