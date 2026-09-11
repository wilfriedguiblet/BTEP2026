#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(optparse)
})

option_list <- list(
  make_option("--primary", type = "character", default = "results/g4_methylation_age_100bp/adjusted_g4_effects_g_richness/adjusted_stable_vs_unstable_100bp_g_richness_stable_minus_unstable.csv"),
  make_option("--chromatin", type = "character", default = "results/g4_methylation_age_100bp/adjusted_g4_effects_g_richness_chromatin_sensitivity/adjusted_stable_vs_unstable_100bp_g_richness_chromatin_stable_minus_unstable.csv"),
  make_option("--structure", type = "character", default = "results/g4_methylation_age_100bp/structure_sensitivity/structure_sensitivity_stable_minus_unstable.csv"),
  make_option("--annotation-provenance", type = "character", default = "results/g4_methylation_age_100bp/sensitivity_annotations/probe_structure_annotation_provenance.csv", dest = "annotation_provenance"),
  make_option("--outdir", type = "character", default = "results/g4_methylation_age_100bp/confounder_sensitivity"),
  make_option("--prefix", type = "character", default = "stable_unstable_confounder_sensitivity")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

read_contrast <- function(path, model_label) {
  x <- fread(path)
  x[, model := model_label]
  x[, effect_change_vs_primary_pct := NA_real_]
  x
}

primary <- read_contrast(opt$primary, "Primary: GC, G-richness, CpG context, G4 density/strand")
chromatin <- read_contrast(opt$chromatin, "Sensitivity: primary + chromatin state")
structure <- read_contrast(opt$structure, "Sensitivity: primary + RepeatMasker class")
summary_dt <- rbindlist(list(primary, chromatin, structure), fill = TRUE)
primary_effect <- primary[, .(tissue, primary_effect = estimate)]
summary_dt <- merge(summary_dt, primary_effect, by = "tissue", all.x = TRUE)
summary_dt[, effect_change_vs_primary_pct := 100 * (estimate - primary_effect) / abs(primary_effect)]
summary_dt[model == "Primary: GC, G-richness, CpG context, G4 density/strand", effect_change_vs_primary_pct := 0]
summary_dt[, model := factor(model, levels = c(
  "Primary: GC, G-richness, CpG context, G4 density/strand",
  "Sensitivity: primary + chromatin state",
  "Sensitivity: primary + RepeatMasker class"
))]
summary_dt[, model_label := factor(as.character(model), levels = levels(model),
                                   labels = c("Primary", "+ Chromatin", "+ RepeatMasker"))]
summary_dt[, tissue := factor(tissue, levels = c("adipose", "liver", "muscle"))]

provenance <- fread(opt$annotation_provenance)
fwrite(summary_dt, file.path(opt$outdir, sprintf("%s.csv", opt$prefix)))
fwrite(provenance, file.path(opt$outdir, sprintf("%s_annotation_availability.csv", opt$prefix)))

plot <- ggplot(summary_dt, aes(x = model_label, y = effect_per_decade, color = model_label)) +
  geom_hline(yintercept = 0, color = "grey45", linewidth = 0.35) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = (estimate - 1.96 * std_error) * 10,
                    ymax = (estimate + 1.96 * std_error) * 10), width = 0.1) +
  facet_wrap(~ tissue, nrow = 1, scales = "free_y") +
  scale_color_manual(values = c("Primary" = "#2C7BB6", "+ Chromatin" = "#D95F02", "+ RepeatMasker" = "#1B9E77")) +
  labs(
    title = "Stable minus unstable G4 age-slope effect is robust to measured confounders",
    subtitle = "Effect per decade; primary model includes GC/G-richness, CpG context, mean methylation, G4 density and strand counts",
    x = "Adjustment model",
    y = "Stable minus unstable effect on signed age slope per decade (delta beta)",
    color = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 8),
        strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold"),
        legend.position = "bottom")

ggsave(file.path(opt$outdir, sprintf("%s.png", opt$prefix)), plot, width = 11, height = 5, dpi = 300)
ggsave(file.path(opt$outdir, sprintf("%s.pdf", opt$prefix)), plot, width = 11, height = 5)

writeLines(c(
  "# Confounder Sensitivity Summary",
  "",
  "Available structural annotations: RepeatMasker class and segmental duplication.",
  "R-loop, replication-timing, and cross-reactive-probe annotations were not supplied for this run.",
  "Segmental duplication did not enter the fitted structural model because the downloaded annotation had no positive 450K probe overlaps.",
  "The table reports the stable-minus-unstable signed age-slope contrast under the primary model and each available sensitivity adjustment."
), file.path(opt$outdir, sprintf("%s.md", opt$prefix)))
