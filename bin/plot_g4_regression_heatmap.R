#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(optparse)
})

option_list <- list(
  make_option("--slopes-dir", type = "character", default = "results/g4_methylation_age/per_tissue", dest = "slopes_dir"),
  make_option("--outdir", type = "character", default = "results/g4_methylation_age/figures"),
  make_option("--prefix", type = "character", default = "g4_regression_heatmap")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
files <- list.files(opt$slopes_dir, pattern = "_probe_slopes\\.csv\\.gz$", full.names = TRUE)
if (!length(files)) stop("No *_probe_slopes.csv.gz files found in ", opt$slopes_dir)

slopes <- rbindlist(lapply(files, fread), fill = TRUE)
required <- c("tissue", "slope", "abs_slope", "stable_motifs", "unstable_motifs", "is_long_cluster", "motif_count")
missing <- setdiff(required, names(slopes))
if (length(missing)) stop("Probe slope files are missing columns: ", paste(missing, collapse = ", "))

slopes[, `:=`(
  has_stable_g4 = stable_motifs > 0,
  has_unstable_g4 = unstable_motifs > 0,
  is_long_cluster = as.logical(is_long_cluster),
  log1p_motif_count = log1p(motif_count)
)]

fit_response <- function(dt, response, response_label) {
  keep <- is.finite(dt[[response]])
  model_dt <- dt[keep, list(
    y = dt[[response]][keep],
    has_stable_g4 = get("has_stable_g4"),
    has_unstable_g4 = get("has_unstable_g4"),
    is_long_cluster = get("is_long_cluster"),
    log1p_motif_count = get("log1p_motif_count")
  )]
  if (nrow(model_dt) < 20) return(data.table())
  fit <- lm(y ~ has_stable_g4 * has_unstable_g4 + is_long_cluster + log1p_motif_count, data = model_dt)
  out <- as.data.table(summary(fit)$coefficients, keep.rownames = "term")
  setnames(out, c("term", "estimate", "std_error", "statistic", "p_value"))
  out[, `:=`(response = response_label, n_probes = nrow(model_dt))]
  out[]
}

effects <- slopes[, rbindlist(list(
  fit_response(.SD, "slope", "Signed age slope")
), fill = TRUE), by = tissue]

effects <- effects[term != "(Intercept)"]
effects[, q_value := p.adjust(p_value, method = "BH"), by = response]
effects[, term_label := fifelse(term == "has_stable_g4TRUE", "Stable G4",
                         fifelse(term == "has_unstable_g4TRUE", "Unstable G4",
                         fifelse(term == "has_stable_g4TRUE:has_unstable_g4TRUE", "Stable x unstable",
                         fifelse(term == "is_long_clusterTRUE", "Long G4 cluster",
                         fifelse(term == "log1p_motif_count", "Motif density", term)))))]
effects[, term_label := factor(term_label, levels = c(
  "Stable G4", "Unstable G4", "Stable x unstable", "Long G4 cluster", "Motif density"
))]
effects[, tissue := factor(tissue, levels = c("adipose", "liver", "muscle"))]
effects[, sig_label := fifelse(q_value < 0.001, "***",
                        fifelse(q_value < 0.01, "**",
                        fifelse(q_value < 0.05, "*", "")))]
effects[, estimate_label := paste0(sprintf("%.2g", estimate), sig_label)]

max_abs <- max(abs(effects$estimate), na.rm = TRUE)
figure <- ggplot(effects, aes(x = tissue, y = term_label, fill = estimate)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = estimate_label), size = 3.4, color = "grey10") +
  scale_fill_gradient2(
    low = "#2C7BB6", mid = "#F7F7F7", high = "#D7191C",
    midpoint = 0, limits = c(-max_abs, max_abs), name = "Regression\nestimate"
  ) +
  labs(
    title = "G4 predictors of signed CpG methylation aging",
    subtitle = "Probe-level regression on signed age slopes; labels show estimates with BH significance (* q<0.05, ** q<0.01, *** q<0.001)",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

png_file <- file.path(opt$outdir, sprintf("%s_signed_only.png", opt$prefix))
pdf_file <- file.path(opt$outdir, sprintf("%s_signed_only.pdf", opt$prefix))
csv_file <- file.path(opt$outdir, sprintf("%s_effects.csv", opt$prefix))

ggsave(png_file, figure, width = 7.2, height = 4.8, dpi = 300)
ggsave(pdf_file, figure, width = 7.2, height = 4.8)
fwrite(effects, csv_file)

message("Wrote: ", png_file)
message("Wrote: ", pdf_file)
message("Wrote: ", csv_file)
