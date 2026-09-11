#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(optparse)
})

option_list <- list(
  make_option("--slopes-dir", type = "character", default = "results/g4_methylation_age/per_tissue", dest = "slopes_dir"),
  make_option("--outdir", type = "character", default = "results/g4_methylation_age/figures"),
  make_option("--prefix", type = "character", default = "stable_vs_unstable_g4")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
files <- list.files(opt$slopes_dir, pattern = "_probe_slopes\\.csv\\.gz$", full.names = TRUE)
if (!length(files)) stop("No *_probe_slopes.csv.gz files found in ", opt$slopes_dir)

slopes <- rbindlist(lapply(files, fread), fill = TRUE)
plot_dt <- slopes[
  g4_class %in% c("stable_only", "unstable_only") &
    is_long_cluster == FALSE &
    is.finite(abs_slope)
]
if (!nrow(plot_dt)) stop("No stable_only or unstable_only non-long-cluster probes found")

plot_dt[, class := fifelse(g4_class == "stable_only", "Stable G4 (test)", "Unstable G4 (control)")]
plot_dt[, class := factor(class, levels = c("Stable G4 (test)", "Unstable G4 (control)"))]
plot_dt[, tissue := factor(tissue, levels = c("adipose", "liver", "muscle"))]

summary_dt <- plot_dt[, .(
  n_probes = .N,
  mean_abs_slope = mean(abs_slope),
  median_abs_slope = median(abs_slope),
  q25_abs_slope = quantile(abs_slope, 0.25),
  q75_abs_slope = quantile(abs_slope, 0.75),
  mean_signed_slope = mean(slope),
  frac_hypermethylating = mean(slope > 0),
  frac_nominal_age_associated = mean(p_value < 0.05, na.rm = TRUE)
), by = .(cohort, tissue, class)]

tests <- rbindlist(lapply(levels(plot_dt$tissue), function(tissue_name) {
  tissue_dt <- plot_dt[tissue == tissue_name]
  stable <- tissue_dt[class == "Stable G4 (test)", abs_slope]
  unstable <- tissue_dt[class == "Unstable G4 (control)", abs_slope]
  signed_stable <- tissue_dt[class == "Stable G4 (test)", slope]
  signed_unstable <- tissue_dt[class == "Unstable G4 (control)", slope]
  if (length(stable) < 2 || length(unstable) < 2) {
    return(data.table(tissue = tissue_name, comparison = c("abs_slope", "signed_slope"),
                      statistic = NA_real_, p_value = NA_real_, n_stable = length(stable),
                      n_unstable = length(unstable)))
  }
  abs_test <- wilcox.test(stable, unstable, alternative = "greater")
  signed_test <- wilcox.test(signed_stable, signed_unstable, alternative = "two.sided")
  data.table(
    tissue = tissue_name,
    comparison = c("abs_slope_stable_greater_than_unstable", "signed_slope_stable_vs_unstable"),
    statistic = c(unname(abs_test$statistic), unname(signed_test$statistic)),
    p_value = c(abs_test$p.value, signed_test$p.value),
    n_stable = length(stable),
    n_unstable = length(unstable)
  )
}), fill = TRUE)

label_dt <- summary_dt[, .(
  label = paste0("n=", format(n_probes, big.mark = ",")),
  y = q75_abs_slope * 1.18
), by = .(tissue, class)]

pval_dt <- tests[comparison == "abs_slope_stable_greater_than_unstable",
  .(tissue, p_label = paste0("Wilcoxon P=", format.pval(p_value, digits = 2, eps = 1e-300)))]
label_y <- summary_dt[, .(y = max(q75_abs_slope, na.rm = TRUE) * 1.45), by = tissue]
pval_dt <- merge(pval_dt, label_y, by = "tissue", all.x = TRUE)
pval_dt[, x := 1.5]

figure <- ggplot(summary_dt, aes(x = class, y = median_abs_slope, fill = class)) +
  geom_col(width = 0.64, color = "grey20", linewidth = 0.25) +
  geom_errorbar(aes(ymin = q25_abs_slope, ymax = q75_abs_slope), width = 0.18, linewidth = 0.45) +
  geom_text(data = label_dt, aes(y = y, label = label), size = 3.1, vjust = 0, color = "grey15") +
  geom_text(data = pval_dt, aes(x = x, y = y, label = p_label), inherit.aes = FALSE, size = 3.2, color = "grey10") +
  facet_wrap(~ tissue, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = c("Stable G4 (test)" = "#1B9E77", "Unstable G4 (control)" = "#D95F02")) +
  labs(
    title = "Age-related CpG methylation drift near stable versus unstable G4s",
    subtitle = "Non-long-cluster G4 +/- 1 kb windows; bars show median absolute age slope, whiskers show IQR",
    x = NULL,
    y = "Absolute methylation age slope (delta beta per year)",
    fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(face = "bold")
  )

png_file <- file.path(opt$outdir, sprintf("%s_by_tissue.png", opt$prefix))
pdf_file <- file.path(opt$outdir, sprintf("%s_by_tissue.pdf", opt$prefix))
ggsave(png_file, figure, width = 10, height = 4.8, dpi = 300)
ggsave(pdf_file, figure, width = 10, height = 4.8)

fwrite(summary_dt, file.path(opt$outdir, sprintf("%s_plot_summary.csv", opt$prefix)))
fwrite(tests, file.path(opt$outdir, sprintf("%s_tests.csv", opt$prefix)))

message("Wrote: ", png_file)
message("Wrote: ", pdf_file)
