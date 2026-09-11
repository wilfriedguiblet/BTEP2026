#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(optparse)
})

option_list <- list(
  make_option("--slopes-dir", type = "character", default = "results/g4_methylation_age/per_tissue", dest = "slopes_dir"),
  make_option("--chromatin", type = "character", default = "../results/v2_hallmark/chromatin_stratification/v2/v2_chromatin_probe_assignment.csv"),
  make_option("--outdir", type = "character", default = "results/g4_methylation_age/chromatin_capture"),
  make_option("--prefix", type = "character", default = "g4_chromatin_capture")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

slope_files <- list.files(opt$slopes_dir, pattern = "_probe_slopes\\.csv\\.gz$", full.names = TRUE)
if (!length(slope_files)) stop("No slope files found")
slopes <- rbindlist(lapply(slope_files, fread), fill = TRUE)
chrom <- fread(opt$chromatin)
chrom <- unique(chrom[, .(probe_id = as.character(probe), chrom_state = as.character(state))], by = "probe_id")
slopes <- merge(slopes, chrom, by = "probe_id", all.x = TRUE)
slopes[, `:=`(
  has_stable_g4 = stable_motifs > 0,
  has_unstable_g4 = unstable_motifs > 0,
  is_long_cluster = as.logical(is_long_cluster),
  log1p_motif_count = log1p(motif_count),
  chrom_state = fifelse(is.na(chrom_state), "Unknown_chromatin", chrom_state)
)]

extract_terms <- function(fit, tissue_name, model_name) {
  out <- as.data.table(summary(fit)$coefficients, keep.rownames = "term")
  setnames(out, c("term", "estimate", "std_error", "statistic", "p_value"))
  out[, `:=`(tissue = tissue_name, model = model_name)]
  out[]
}

fit_tissue <- function(tissue_name) {
  dt <- slopes[slopes[["tissue"]] == tissue_name & is.finite(slopes[["slope"]])]
  dt <- dt[complete.cases(dt[, c("slope", "has_stable_g4", "has_unstable_g4", "is_long_cluster", "log1p_motif_count", "chrom_state"), with = FALSE])]
  models <- list(
    g4_only = slope ~ has_stable_g4 * has_unstable_g4 + is_long_cluster + log1p_motif_count,
    chromatin_only = slope ~ chrom_state,
    chromatin_plus_g4 = slope ~ chrom_state + has_stable_g4 * has_unstable_g4 + is_long_cluster + log1p_motif_count,
    chromatin_g4_interactions = slope ~ chrom_state * (has_stable_g4 + has_unstable_g4 + is_long_cluster + log1p_motif_count)
  )
  fits <- lapply(models, lm, data = dt)
  perf <- rbindlist(lapply(names(fits), function(model_name) {
    fit <- fits[[model_name]]
    data.table(
      tissue = tissue_name,
      model = model_name,
      n_probes = nrow(dt),
      df = length(coef(fit)),
      r2 = summary(fit)$r.squared,
      adj_r2 = summary(fit)$adj.r.squared,
      aic = AIC(fit)
    )
  }))
  coefs <- rbindlist(lapply(names(fits), function(model_name) extract_terms(fits[[model_name]], tissue_name, model_name)))
  list(perf = perf, coefs = coefs)
}

fits <- lapply(sort(unique(slopes$tissue)), fit_tissue)
performance <- rbindlist(lapply(fits, `[[`, "perf"), fill = TRUE)
coefficients <- rbindlist(lapply(fits, `[[`, "coefs"), fill = TRUE)

wide <- dcast(performance, tissue ~ model, value.var = "r2")
wide[, `:=`(
  g4_increment_after_chromatin = chromatin_plus_g4 - chromatin_only,
  interaction_increment_after_additive = chromatin_g4_interactions - chromatin_plus_g4,
  chromatin_increment_after_g4 = chromatin_plus_g4 - g4_only
)]
fwrite(performance, file.path(opt$outdir, sprintf("%s_model_performance.csv", opt$prefix)))
fwrite(wide, file.path(opt$outdir, sprintf("%s_r2_increments.csv", opt$prefix)))
fwrite(coefficients, file.path(opt$outdir, sprintf("%s_coefficients.csv", opt$prefix)))

plot_dt <- melt(performance, id.vars = c("tissue", "model"), measure.vars = "r2", variable.name = "metric", value.name = "value")
plot_dt[, model := factor(model, levels = c("g4_only", "chromatin_only", "chromatin_plus_g4", "chromatin_g4_interactions"))]
plot <- ggplot(plot_dt, aes(x = model, y = value, fill = model)) +
  geom_col(width = 0.72, color = "grey25", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.3f", value)), vjust = -0.35, size = 3.2) +
  facet_wrap(~ tissue, nrow = 1) +
  scale_fill_manual(values = c(
    g4_only = "#80CDC1",
    chromatin_only = "#DFC27D",
    chromatin_plus_g4 = "#A6611A",
    chromatin_g4_interactions = "#018571"
  )) +
  coord_cartesian(ylim = c(0, max(plot_dt$value) * 1.15)) +
  labs(title = "Does chromatin state capture the G4 age-slope signal?",
       subtitle = "Nested probe-level signed age-slope models; compare R2 increments after adding G4 terms to chromatin",
       x = NULL, y = "R2", fill = NULL) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"), legend.position = "bottom")

ggsave(file.path(opt$outdir, sprintf("%s_model_r2.png", opt$prefix)), plot, width = 10, height = 5, dpi = 300)
ggsave(file.path(opt$outdir, sprintf("%s_model_r2.pdf", opt$prefix)), plot, width = 10, height = 5)

message("Wrote outputs to: ", opt$outdir)
