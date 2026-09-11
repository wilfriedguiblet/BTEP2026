#!/usr/bin/env Rscript

# HUMAN AUDIT: the paper uses only the two minimal Quadron models in this file.
# M0: slope ~ has_stable_g4 + has_unstable_g4 + mean_beta
# M1: M0 + probe_strand_g_fraction + opposite_strand_g_fraction
# Enable --minimal-model true for BOTH; enable --include-strand-g-richness true
# only for M1. Supply the Quadron input paths explicitly: historical defaults
# below belong to exploratory runs and are NOT the paper's input selection.
# The response is an already fitted signed beta-versus-age slope (beta/year),
# not absolute drift, a methylation clock, or a longitudinal change per person.
# Comments identify the implementation that produced the saved results; they
# do not imply that the model assumptions or coordinate provenance are settled.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(optparse)
})

option_list <- list(
  make_option("--prepared-dir", type = "character", default = "results/g4_methylation_age/prepared_inputs", dest = "prepared_dir"),
  make_option("--slopes-dir", type = "character", default = "results/g4_methylation_age/per_tissue", dest = "slopes_dir"),
  make_option("--g4-unmerged", type = "character", default = "results/g4_methylation_age_100bp/g4_windows/g4_motifs_100bp_unmerged.bed", dest = "g4_unmerged"),
  make_option("--chromatin", type = "character", default = "../results/v2_hallmark/chromatin_stratification/v2/v2_chromatin_probe_assignment.csv"),
  make_option("--include-chromatin", type = "character", default = "false", dest = "include_chromatin"),
  make_option("--structure-annotations", type = "character", default = "", dest = "structure_annotations"),
  make_option("--exclude-cross-reactive", type = "character", default = "true", dest = "exclude_cross_reactive"),
  make_option("--include-g4-architecture", type = "character", default = "false", dest = "include_g4_architecture"),
  make_option("--minimal-model", type = "character", default = "false", dest = "minimal_model"),
  make_option("--include-strand-g-richness", type = "character", default = "false", dest = "include_strand_g_richness"),
  make_option("--stratify-g4-density", type = "character", default = "false", dest = "stratify_g4_density"),
  make_option("--gc-flank-bp", type = "integer", default = 1000L, dest = "gc_flank_bp"),
  make_option("--outdir", type = "character", default = "results/g4_methylation_age/adjusted_g4_effects"),
  make_option("--prefix", type = "character", default = "adjusted_stable_vs_unstable")
)
opt <- parse_args(OptionParser(option_list = option_list))
include_chromatin <- tolower(opt$include_chromatin) %in% c("true", "t", "1", "yes", "y")
stratify_g4_density <- tolower(opt$stratify_g4_density) %in% c("true", "t", "1", "yes", "y")
exclude_cross_reactive <- tolower(opt$exclude_cross_reactive) %in% c("true", "t", "1", "yes", "y")
include_g4_architecture <- tolower(opt$include_g4_architecture) %in% c("true", "t", "1", "yes", "y")
minimal_model <- tolower(opt$minimal_model) %in% c("true", "t", "1", "yes", "y")
include_strand_g_richness <- tolower(opt$include_strand_g_richness) %in% c("true", "t", "1", "yes", "y")

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
if (!requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE)) {
  stop("Missing IlluminaHumanMethylation450kanno.ilmn12.hg19")
}
if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE) ||
    !requireNamespace("GenomicRanges", quietly = TRUE) ||
    !requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Missing BSgenome.Hsapiens.UCSC.hg38, GenomicRanges, or Biostrings for GC correction")
}

slope_files <- list.files(opt$slopes_dir, pattern = "_probe_slopes\\.csv\\.gz$", full.names = TRUE)
# Every matching file is read. The curated runner stages exactly the three
# paper cohorts in a fresh directory to prevent accidental extra-cohort input.
if (!length(slope_files)) stop("No probe slope files found")
slopes <- rbindlist(lapply(slope_files, fread), fill = TRUE)

if (!"mean_beta" %in% names(slopes)) {
  # Compatibility fallback only. It averages ALL prepared beta columns, which
  # may differ from means after age/sex/BMI filtering. The paper replay requires
  # mean_beta in the existing slope checkpoint and never relies on this fallback.
  beta_files <- list.files(opt$prepared_dir, pattern = "_beta\\.rds$", full.names = TRUE)
  if (!length(beta_files)) stop("No prepared beta RDS files found")
  mean_beta <- rbindlist(lapply(beta_files, function(path) {
    beta <- readRDS(path)
    tissue <- sub("^GSE[0-9]+_", "", sub("_beta\\.rds$", "", basename(path)))
    data.table(probe_id = rownames(beta), tissue = tissue, mean_beta = rowMeans(beta, na.rm = TRUE))
  }))
  slopes <- merge(slopes, mean_beta, by = c("probe_id", "tissue"), all.x = TRUE)
}

data("Islands.UCSC", package = "IlluminaHumanMethylation450kanno.ilmn12.hg19")
island <- data.table(probe_id = rownames(Islands.UCSC), relation = as.character(Islands.UCSC$Relation_to_Island))
island[, cpg_context := fifelse(relation == "Island", "Island",
                         fifelse(relation %in% c("N_Shore", "S_Shore"), "Shore",
                         fifelse(relation %in% c("N_Shelf", "S_Shelf"), "Shelf", "OpenSea")))]
slopes <- merge(slopes, island[, .(probe_id, cpg_context)], by = "probe_id", all.x = TRUE)

if (include_chromatin && file.exists(opt$chromatin)) {
  chrom <- fread(opt$chromatin)
  probe_col <- intersect(c("probe", "probe_id"), names(chrom))[1]
  if (!is.na(probe_col) && "state" %in% names(chrom)) {
    chrom <- unique(chrom[, .(probe_id = as.character(get(probe_col)), chrom_state = as.character(state))], by = "probe_id")
    slopes <- merge(slopes, chrom, by = "probe_id", all.x = TRUE)
  }
}
if (!"chrom_state" %in% names(slopes)) slopes[, chrom_state := NA_character_]

structure_terms <- character()
if (nzchar(opt$structure_annotations) && file.exists(opt$structure_annotations)) {
  # IMPORTANT: this optional join/filter happens even for the minimal models.
  # The curated runner supplies no structure file, explicitly disables the
  # cross-reactive exclusion flag, and checks numerical results against the
  # manuscript. Formula equality alone cannot certify an identical sample set.
  structure <- fread(opt$structure_annotations)
  structure_cols <- intersect(c("probe_id", "repeat_class", "segmental_duplication", "rloop_overlap", "replication_timing", "cross_reactive"), names(structure))
  slopes <- merge(slopes, structure[, ..structure_cols], by = "probe_id", all.x = TRUE)
  if ("cross_reactive" %in% names(slopes) && exclude_cross_reactive) slopes <- slopes[is.na(cross_reactive) | cross_reactive == FALSE]
  if ("repeat_class" %in% names(slopes)) {
    slopes[, repeat_class := factor(fifelse(is.na(repeat_class), "No_repeat", repeat_class))]
    structure_terms <- c(structure_terms, "repeat_class")
  }
  for (column in c("segmental_duplication", "rloop_overlap")) {
    if (column %in% names(slopes) && any(slopes[[column]] %in% TRUE, na.rm = TRUE)) {
      slopes[, (column) := as.logical(get(column))]
      structure_terms <- c(structure_terms, column)
    }
  }
  if ("replication_timing" %in% names(slopes) && sum(is.finite(slopes[["replication_timing"]])) > 10) structure_terms <- c(structure_terms, "replication_timing")
}

data("Locations", package = "IlluminaHumanMethylation450kanno.ilmn12.hg19")
locations_df <- as.data.frame(Locations)
probe_strand <- data.table(probe_id = rownames(locations_df), probe_strand = as.character(locations_df$strand))
slopes <- merge(slopes, probe_strand, by = "probe_id", all.x = TRUE)

compute_sequence_features <- function(probe_dt, flank_bp) {
  # These sequence covariates are centered on each probe, not on G4 motifs.
  # The implementation uses floor((start+end)/2), then a ONE-BASED INCLUSIVE
  # GRanges interval [center-flank, center+flank], clipped at chromosome ends.
  # At the default flank=1000 this is normally 2001 bases. The source manifest's
  # coordinate convention and hg19 strand-to-hg38 agreement remain audit points;
  # do not insert a +/-1 coordinate change in an exact-reproduction branch.
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
  # Denominator includes the whole extracted sequence length (including any N),
  # rather than counting G/C only among known A/C/G/T bases.
  probes[["gc_fraction"]] <- rowSums(counts) / Biostrings::width(seqs)
  probes[["g_fraction"]] <- counts[, "G"] / Biostrings::width(seqs)
  probes[["c_fraction"]] <- counts[, "C"] / Biostrings::width(seqs)
  probes[, c("probe_id", "gc_fraction", "g_fraction", "c_fraction"), with = FALSE]
}
sequence_dt <- compute_sequence_features(slopes, opt$gc_flank_bp)
slopes <- merge(slopes, sequence_dt, by = "probe_id", all.x = TRUE)
# Probe orientation comes from the 450K hg19 Locations annotation. On '+' the
# two fractions are reference G/C; on '-' they are C/G. Their sum is GC content.
# This is neither transcriptional orientation nor an independent strand assay.
slopes[, probe_strand_g_fraction := fifelse(probe_strand == "+", g_fraction, c_fraction)]
slopes[, opposite_strand_g_fraction := fifelse(probe_strand == "+", c_fraction, g_fraction)]

compute_g4_strand_features <- function(probe_dt, g4_unmerged) {
  # These are MOTIF COUNTS, unlike the base-composition fractions above. They
  # are computed by this shared script even when omitted from minimal formulas.
  # This is why the existing unmerged nine-column BED is still a dependency of
  # a historical-code replay, including M0. Removing it is a refactor, not needed
  # to reproduce the coefficients. No motif counts enter the paper's M0/M1.
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
strand_summary <- unique(slopes[, c("probe_id", "probe_strand", strand_cols), with = FALSE], by = "probe_id")[, lapply(.SD, function(x) {
  if (is.numeric(x) || is.integer(x)) sum(x > 0, na.rm = TRUE) else sum(!is.na(x) & x != "")
}), .SDcols = c("probe_strand", strand_cols)]
strand_summary[, total_probes := uniqueN(slopes$probe_id)]
fwrite(strand_summary, file.path(opt$outdir, sprintf("%s_strand_feature_summary.csv", opt$prefix)))

slopes[, `:=`(
  has_stable_g4 = stable_motifs > 0,
  has_unstable_g4 = unstable_motifs > 0,
  is_long_cluster = as.logical(is_long_cluster),
  log1p_motif_count = log1p(motif_count),
  cpg_context = factor(fifelse(is.na(cpg_context), "Unknown", cpg_context), levels = c("Island", "Shore", "Shelf", "OpenSea", "Unknown")),
  chrom_state = factor(fifelse(is.na(chrom_state), "Unknown_chromatin", chrom_state))
)]

fit_one <- function(tissue_name, density_group = "all") {
  dt <- slopes[slopes[["tissue"]] == tissue_name & is.finite(slopes[["slope"]]) & is.finite(slopes[["mean_beta"]])]
  if (density_group == "standard") dt <- dt[is_long_cluster == FALSE]
  if (density_group == "dense_cluster") dt <- dt[is_long_cluster == TRUE]
  strand_predictors <- c("log1p_stable_same_strand", "log1p_stable_opposite_strand", "log1p_unstable_same_strand", "log1p_unstable_opposite_strand")
  # This is the branch used by both manuscript models. Do not add density,
  # chromatin, CpG-island, repeat, or strand-specific MOTIF COUNT predictors
  # here when reproducing the paper. Such a change defines a new analysis.
  # mean_beta is the within-tissue mean over the same specimens as the slopes,
  # not an independently observed baseline methylation measurement.
  if (minimal_model) {
    cofactor_terms <- "mean_beta"
    if (include_strand_g_richness) {
      cofactor_terms <- c(cofactor_terms, "probe_strand_g_fraction", "opposite_strand_g_fraction")
    }
    architecture_terms <- character()
    density_term <- character()
  } else {
    cofactor_terms <- c("mean_beta", "gc_fraction", "g_fraction", "cpg_context")
    if (include_chromatin) cofactor_terms <- c(cofactor_terms, "chrom_state")
    cofactor_terms <- c(cofactor_terms, structure_terms)
    architecture_terms <- if (include_g4_architecture) c("log1p_motif_count", strand_predictors) else character()
    density_term <- if (include_g4_architecture && density_group == "all") "is_long_cluster" else character()
  }
  complete_cols <- c("slope", "has_stable_g4", "has_unstable_g4", density_term, architecture_terms, cofactor_terms)
  dt <- dt[complete.cases(dt[, complete_cols, with = FALSE])]
  formula_terms <- c("has_stable_g4", "has_unstable_g4", density_term, architecture_terms, cofactor_terms)
  formula_text <- paste("slope ~", paste(formula_terms, collapse = " + "))
  # Each eligible CpG receives equal weight. The first-stage slope SE is NOT
  # propagated, and the residual model does NOT account for genomic correlation
  # or shared-sample estimation error. Preserve this for reproduction; a
  # dependency-aware or weighted fit would be a separately labelled extension.
  fit <- lm(as.formula(formula_text), data = dt)
  coef_dt <- as.data.table(summary(fit)$coefficients, keep.rownames = "term")
  setnames(coef_dt, c("term", "estimate", "std_error", "statistic", "p_value"))
  coef_dt[, `:=`(tissue = tissue_name, density_group = density_group, n_probes = nrow(dt), model = formula_text)]

  vc <- vcov(fit)
  b <- coef(fit)
  stable_term <- "has_stable_g4TRUE"
  unstable_term <- "has_unstable_g4TRUE"
  contrast_estimable <- stable_term %in% names(b) && unstable_term %in% names(b) &&
    is.finite(b[[stable_term]]) && is.finite(b[[unstable_term]])
  # Test the difference DIRECTLY. Comparing a significant stable term with a
  # non-significant unstable term is not a test that the two classes differ.
  # The contrast SE below includes covariance between the fitted coefficients.
  if (contrast_estimable) {
    contrast <- rep(0, length(b))
    names(contrast) <- names(b)
    contrast[stable_term] <- 1
    contrast[unstable_term] <- -1
    estimable <- is.finite(b)
    contrast <- contrast[estimable]
    b <- b[estimable]
    vc <- vc[estimable, estimable, drop = FALSE]
    diff_est <- sum(contrast * b)
    diff_se <- sqrt(as.numeric(t(contrast) %*% vc %*% contrast))
    diff_t <- diff_est / diff_se
    diff_p <- 2 * pt(-abs(diff_t), df = df.residual(fit))
  } else {
    diff_est <- NA_real_
    diff_se <- NA_real_
    diff_t <- NA_real_
    diff_p <- NA_real_
  }
  contrast_dt <- data.table(
    tissue = tissue_name,
    density_group = density_group,
    contrast = "stable_minus_unstable",
    contrast_estimable = contrast_estimable,
    estimate = diff_est,
    std_error = diff_se,
    statistic = diff_t,
    p_value = diff_p,
    effect_per_decade = diff_est * 10,
    n_probes = nrow(dt),
    model = formula_text
  )
  list(coef = coef_dt, contrast = contrast_dt)
}

groups <- if (stratify_g4_density) c("standard", "dense_cluster") else "all"
fits <- unlist(lapply(sort(unique(slopes$tissue)), function(tissue_name) {
  lapply(groups, function(group) fit_one(tissue_name, group))
}), recursive = FALSE)
coef_dt <- rbindlist(lapply(fits, `[[`, "coef"), fill = TRUE)
contrast_dt <- rbindlist(lapply(fits, `[[`, "contrast"), fill = TRUE)
# These are separate BH families. With density_group='all', each family has
# three tissues: one family per coefficient term, plus a distinct family for
# stable-minus-unstable contrasts. M0 and M1 are adjusted separately. The q
# values do not correct across every term, both models, or earlier explorations.
contrast_dt[, q_value := p.adjust(p_value, method = "BH"), by = density_group]
coef_dt[, q_value := p.adjust(p_value, method = "BH"), by = .(term, density_group)]

g4_terms <- coef_dt[term %in% c("has_stable_g4TRUE", "has_unstable_g4TRUE"),
                    .(tissue, density_group, term, estimate, std_error, p_value, q_value, n_probes)]
g4_terms[, predictor := fifelse(term == "has_stable_g4TRUE", "Stable G4", "Unstable G4")]
# Convert a slope DIFFERENCE from beta/year to beta/decade by multiplying by
# ten. Multiply beta/decade by another 100 only when reporting percentage points.
g4_terms[, effect_per_decade := estimate * 10]

fwrite(coef_dt, file.path(opt$outdir, sprintf("%s_coefficients.csv", opt$prefix)))
fwrite(g4_terms, file.path(opt$outdir, sprintf("%s_g4_terms.csv", opt$prefix)))
fwrite(contrast_dt, file.path(opt$outdir, sprintf("%s_stable_minus_unstable.csv", opt$prefix)))

plot_dt <- copy(g4_terms)
plot_dt[, predictor := factor(predictor, levels = c("Stable G4", "Unstable G4"))]
plot_dt[, tissue := factor(tissue, levels = c("adipose", "liver", "muscle"))]
plot_dt[, density_group := factor(density_group, levels = c("standard", "dense_cluster", "all"),
                                  labels = c("Standard G4 neighborhoods", "Dense G4 clusters", "All G4 neighborhoods"))]
plot_dt[, significance := fifelse(q_value < 0.001, "***",
                           fifelse(q_value < 0.01, "**",
                           fifelse(q_value < 0.05, "*", "ns")))]
# Whiskers are ordinary approximate 95% intervals, NOT BH-adjusted intervals.
# Stars are term tests; brackets below are different, direct class comparisons.
# Recreated PNG/PDF bytes may vary with R, ggplot2, fonts, and graphics device.
# Validate numerical tables first; retain the original paper PNGs unchanged.
plot_dt[, label_y := (estimate + 1.96 * std_error) * 10]
plot_dt[, label_y := label_y + ifelse(label_y >= 0, 1, -1) * max(abs(effect_per_decade), na.rm = TRUE) * 0.04]
contrast_display <- copy(contrast_dt)
contrast_display[, density_group := factor(density_group, levels = c("standard", "dense_cluster", "all"),
                                            labels = c("Standard G4 neighborhoods", "Dense G4 clusters", "All G4 neighborhoods"))]
contrast_plot <- merge(
  contrast_display[, .(tissue, density_group, contrast_estimable, contrast_q_value = q_value)],
  plot_dt[, .(panel_top = max((estimate + 1.96 * std_error) * 10, na.rm = TRUE)), by = .(tissue, density_group)],
  by = c("tissue", "density_group"),
  all.x = TRUE
)
contrast_plot[, density_group := factor(density_group, levels = levels(plot_dt$density_group))]
contrast_plot[, tissue := factor(tissue, levels = c("adipose", "liver", "muscle"))]
contrast_plot[, bracket_y := panel_top + pmax(abs(panel_top) * 0.18, 0.00003)]
contrast_plot[, tip_y := bracket_y - pmax(abs(panel_top) * 0.035, 0.000006)]
contrast_plot[, contrast_label := fifelse(!contrast_estimable, "Stable vs unstable: not estimable",
                                   fifelse(contrast_q_value < 0.001, "Stable vs unstable: ***",
                                   fifelse(contrast_q_value < 0.01, "Stable vs unstable: **",
                                   fifelse(contrast_q_value < 0.05, "Stable vs unstable: *", "Stable vs unstable: ns"))))]
plot <- ggplot(plot_dt, aes(x = predictor, y = effect_per_decade, fill = predictor)) +
  geom_hline(yintercept = 0, linewidth = 0.35, color = "grey35") +
  geom_col(width = 0.68, color = "grey20", linewidth = 0.25) +
  geom_errorbar(aes(ymin = (estimate - 1.96 * std_error) * 10,
                    ymax = (estimate + 1.96 * std_error) * 10), width = 0.18) +
  geom_text(aes(y = label_y, label = significance), vjust = ifelse(plot_dt$label_y >= 0, 0, 1), size = 4.2, fontface = "bold") +
  geom_segment(data = contrast_plot, aes(x = 1, xend = 2, y = bracket_y, yend = bracket_y), inherit.aes = FALSE, linewidth = 0.4) +
  geom_segment(data = contrast_plot, aes(x = 1, xend = 1, y = tip_y, yend = bracket_y), inherit.aes = FALSE, linewidth = 0.4) +
  geom_segment(data = contrast_plot, aes(x = 2, xend = 2, y = tip_y, yend = bracket_y), inherit.aes = FALSE, linewidth = 0.4) +
  geom_text(data = contrast_plot, aes(x = 1.5, y = bracket_y, label = contrast_label), inherit.aes = FALSE, vjust = -0.45, size = 3.0, fontface = "bold") +
  facet_grid(density_group ~ tissue, scales = "free_y") +
  scale_fill_manual(values = c("Stable G4" = "#1B9E77", "Unstable G4" = "#D95F02")) +
  labs(title = "Adjusted stable versus unstable G4 effects on signed methylation aging",
      subtitle = if (minimal_model) paste0("Quadron G4 annotation: +/- 100 bp; minimal model: stable G4 + unstable G4 + mean methylation", if (include_strand_g_richness) sprintf(" + probe/oriented G richness (+/- %d bp sequence window)", opt$gc_flank_bp) else "", "; bars and brackets use BH q: * <0.05, ** <0.01, *** <0.001") else sprintf("G4 annotation: +/- 100 bp; adjusted for mean methylation, GC/G +/- %d bp sequence window, CpG context%s%s; bars and brackets use BH q: * <0.05, ** <0.01, *** <0.001", opt$gc_flank_bp, if (include_chromatin) ", chromatin state" else "", if (include_g4_architecture) ", and G4 architecture" else ""),
       x = NULL, y = "Effect on signed age slope per decade (delta beta)", fill = NULL) +
  theme_classic(base_size = 12) +
  theme(strip.text = element_text(face = "bold"), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.title = element_text(face = "bold"), legend.position = "bottom")

plot_height <- if (stratify_g4_density) 7.2 else 4.4
ggsave(file.path(opt$outdir, sprintf("%s_g4_terms.png", opt$prefix)), plot, width = 8.8, height = plot_height, dpi = 300)
ggsave(file.path(opt$outdir, sprintf("%s_g4_terms.pdf", opt$prefix)), plot, width = 8.8, height = plot_height)

cat("Adjusted G4 terms\n")
print(g4_terms)
cat("\nStable minus unstable contrasts\n")
print(contrast_dt)
