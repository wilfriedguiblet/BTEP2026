#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(optparse)
})

option_list <- list(
  make_option("--prepared-dir", type = "character", default = "results/g4_methylation_age/prepared_inputs", dest = "prepared_dir"),
  make_option("--slopes-dir", type = "character", default = "results/g4_methylation_age/per_tissue", dest = "slopes_dir"),
  make_option("--outdir", type = "character", default = "results/g4_methylation_age/figures"),
  make_option("--prefix", type = "character", default = "methylation_dataset_correlations")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

read_beta <- function(path) {
  beta <- readRDS(path)
  storage.mode(beta) <- "numeric"
  beta
}

beta_files <- list.files(opt$prepared_dir, pattern = "_beta\\.rds$", full.names = TRUE)
metadata_files <- list.files(opt$prepared_dir, pattern = "_metadata\\.csv$", full.names = TRUE)
if (!length(beta_files) || !length(metadata_files)) stop("Prepared beta/metadata files not found")

metadata <- rbindlist(lapply(metadata_files, fread), fill = TRUE)
if (!"subjectid" %in% names(metadata)) stop("Metadata does not contain subjectid")

beta_list <- lapply(beta_files, read_beta)
names(beta_list) <- sub("^GSE[0-9]+_", "", sub("_beta\\.rds$", "", basename(beta_files)))

common_probes <- Reduce(intersect, lapply(beta_list, rownames))
if (length(common_probes) < 10) stop("Too few common CpG probes for correlation")
beta_list <- lapply(beta_list, function(x) x[common_probes, , drop = FALSE])

mean_profiles <- sapply(beta_list, rowMeans, na.rm = TRUE)
mean_cor <- cor(mean_profiles, use = "pairwise.complete.obs", method = "pearson")
mean_cor_dt <- as.data.table(as.table(mean_cor))
setnames(mean_cor_dt, c("tissue_x", "tissue_y", "correlation"))
mean_cor_dt[, metric := "Mean methylation profile"]

slope_files <- list.files(opt$slopes_dir, pattern = "_probe_slopes\\.csv\\.gz$", full.names = TRUE)
if (!length(slope_files)) stop("Probe slope files not found")
slope_dt <- rbindlist(lapply(slope_files, fread), fill = TRUE)
slope_wide <- dcast(slope_dt[probe_id %in% common_probes], probe_id ~ tissue, value.var = "slope")
slope_mat <- as.matrix(slope_wide[, -"probe_id"])
slope_cor <- cor(slope_mat, use = "pairwise.complete.obs", method = "pearson")
slope_cor_dt <- as.data.table(as.table(slope_cor))
setnames(slope_cor_dt, c("tissue_x", "tissue_y", "correlation"))
slope_cor_dt[, metric := "Age-slope profile"]

paired_rows <- list()
pair_index <- 0L
for (tissue_a in names(beta_list)) {
  for (tissue_b in names(beta_list)) {
    if (tissue_a >= tissue_b) next
    meta_a <- metadata[tissue == tissue_a & subjectid != ""]
    meta_b <- metadata[tissue == tissue_b & subjectid != ""]
    shared_subjects <- intersect(meta_a$subjectid, meta_b$subjectid)
    if (!length(shared_subjects)) next
    for (subject in shared_subjects) {
      sample_a <- meta_a[subjectid == subject, sample_id][1]
      sample_b <- meta_b[subjectid == subject, sample_id][1]
      if (!sample_a %in% colnames(beta_list[[tissue_a]]) || !sample_b %in% colnames(beta_list[[tissue_b]])) next
      pair_index <- pair_index + 1L
      paired_rows[[pair_index]] <- data.table(
        tissue_x = tissue_a,
        tissue_y = tissue_b,
        subjectid = subject,
        correlation = cor(beta_list[[tissue_a]][, sample_a], beta_list[[tissue_b]][, sample_b],
                          use = "pairwise.complete.obs", method = "pearson")
      )
    }
  }
}
paired_dt <- if (length(paired_rows)) rbindlist(paired_rows) else data.table()
paired_summary <- if (nrow(paired_dt)) {
  paired_dt[, .(
    n_subjects = .N,
    mean_correlation = mean(correlation, na.rm = TRUE),
    median_correlation = median(correlation, na.rm = TRUE),
    min_correlation = min(correlation, na.rm = TRUE),
    max_correlation = max(correlation, na.rm = TRUE)
  ), by = .(tissue_x, tissue_y)]
} else {
  data.table(tissue_x = character(), tissue_y = character(), n_subjects = integer(),
             mean_correlation = numeric(), median_correlation = numeric(),
             min_correlation = numeric(), max_correlation = numeric())
}

plot_cor <- rbind(mean_cor_dt, slope_cor_dt, fill = TRUE)
plot_cor[, tissue_x := factor(tissue_x, levels = names(beta_list))]
plot_cor[, tissue_y := factor(tissue_y, levels = rev(names(beta_list)))]
plot_cor[, label := sprintf("%.3f", correlation)]

heatmap <- ggplot(plot_cor, aes(x = tissue_x, y = tissue_y, fill = correlation)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = label), size = 4.2, color = "grey10") +
  facet_wrap(~ metric, nrow = 1) +
  scale_fill_gradient2(low = "#2C7BB6", mid = "#F7F7F7", high = "#D7191C", midpoint = 0,
                       limits = c(-1, 1), name = "Pearson r") +
  coord_equal() +
  labs(title = "Cross-tissue methylation dataset correlations", x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"), plot.title = element_text(face = "bold"))

png_file <- file.path(opt$outdir, sprintf("%s_heatmap.png", opt$prefix))
pdf_file <- file.path(opt$outdir, sprintf("%s_heatmap.pdf", opt$prefix))
ggsave(png_file, heatmap, width = 9.5, height = 4.6, dpi = 300)
ggsave(pdf_file, heatmap, width = 9.5, height = 4.6)

if (nrow(paired_dt)) {
  paired_plot <- ggplot(paired_dt, aes(x = interaction(tissue_x, tissue_y, sep = " vs "), y = correlation)) +
    geom_boxplot(width = 0.5, outlier.shape = NA, fill = "#D8B365", color = "grey20") +
    geom_jitter(width = 0.12, height = 0, size = 1.7, alpha = 0.75, color = "grey15") +
    labs(title = "Paired-subject cross-tissue methylation correlations",
         x = NULL, y = "Pearson r across common CpG probes") +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(face = "bold"), plot.title = element_text(face = "bold"))
  ggsave(file.path(opt$outdir, sprintf("%s_paired_subjects.png", opt$prefix)), paired_plot, width = 6, height = 4, dpi = 300)
  ggsave(file.path(opt$outdir, sprintf("%s_paired_subjects.pdf", opt$prefix)), paired_plot, width = 6, height = 4)
}

fwrite(rbind(mean_cor_dt, slope_cor_dt, fill = TRUE), file.path(opt$outdir, sprintf("%s_profile_correlations.csv", opt$prefix)))
fwrite(paired_dt, file.path(opt$outdir, sprintf("%s_paired_subject_correlations.csv", opt$prefix)))
fwrite(paired_summary, file.path(opt$outdir, sprintf("%s_paired_subject_summary.csv", opt$prefix)))

message("Wrote: ", png_file)
message("Wrote: ", pdf_file)
