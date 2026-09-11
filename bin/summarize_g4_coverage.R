#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(optparse)
})

option_list <- list(
  make_option("--core-bed", type = "character", default = "results/g4_methylation_age_100bp/g4_windows/g4_motifs_core.bed", dest = "core_bed"),
  make_option("--window-bed", type = "character", default = "results/g4_methylation_age_100bp/g4_windows/g4_motifs_100bp_unmerged.bed", dest = "window_bed"),
  make_option("--probe-slopes", type = "character", default = "results/g4_methylation_age_100bp/per_tissue/GSE61257_adipose_probe_slopes.csv.gz", dest = "probe_slopes"),
  make_option("--outdir", type = "character", default = "results/g4_methylation_age_100bp/coverage"),
  make_option("--prefix", type = "character", default = "g4_coverage_100bp")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
  stop("Missing BSgenome.Hsapiens.UCSC.hg38")
}
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
chromosomes <- paste0("chr", c(1:22, "X", "Y"))
chrom_lengths <- GenomeInfoDb::seqlengths(genome)[chromosomes]
reference_bp <- sum(chrom_lengths)

read_bed <- function(path, is_core) {
  x <- fread(path, header = FALSE)
  if (is_core) {
    setnames(x, c("chr", "start", "end", "motif_id", "g4_class", "strand", "score", "motif"))
  } else {
    setnames(x, c("chr", "start", "end", "motif_id", "g4_class", "strand", "score", "core_start", "core_end"))
  }
  x[x[["chr"]] %in% chromosomes & x[["g4_class"]] %in% c("stable", "unstable")]
}

core <- read_bed(opt$core_bed, is_core = TRUE)
windows <- read_bed(opt$window_bed, is_core = FALSE)
probes <- fread(opt$probe_slopes, select = c("probe_id", "chr", "start", "end"))
probes <- unique(probes, by = "probe_id")
probes <- probes[probes[["chr"]] %in% chromosomes]

coverage_for <- function(annotation, scope, dt) {
  subset <- dt[dt[["g4_class"]] == annotation]
  starts <- pmax(1L, subset[["start"]] + 1L)
  ends <- pmin(as.integer(chrom_lengths[subset[["chr"]]]), subset[["end"]])
  valid <- is.finite(starts) & is.finite(ends) & ends >= starts
  ranges <- GRanges(subset[["chr"]][valid], IRanges(starts[valid], ends[valid]))
  reduced <- reduce(ranges)
  query <- GRanges(probes[["chr"]], IRanges(probes[["start"]], probes[["end"]]))
  n_probe_overlap <- length(unique(queryHits(findOverlaps(query, reduced, ignore.strand = TRUE))))
  data.table(
    annotation = annotation,
    coordinate_scope = scope,
    motif_count = nrow(subset),
    invalid_or_out_of_bound_intervals = sum(!valid),
    raw_interval_bp = sum(width(ranges)),
    union_bp = sum(width(reduced)),
    genome_bp = reference_bp,
    genome_coverage_pct = 100 * sum(width(reduced)) / reference_bp,
    assayed_cpg_probes = length(query),
    cpg_probes_covered = n_probe_overlap,
    cpg_probe_coverage_pct = 100 * n_probe_overlap / length(query)
  )
}

summary_dt <- rbindlist(list(
  coverage_for("stable", "motif_core", core),
  coverage_for("unstable", "motif_core", core),
  coverage_for("stable", "motif_plus_100bp", windows),
  coverage_for("unstable", "motif_plus_100bp", windows)
))
summary_dt[, interval_overlap_bp := raw_interval_bp - union_bp]
summary_dt[, interval_overlap_pct := 100 * interval_overlap_bp / raw_interval_bp]
setcolorder(summary_dt, c(
  "annotation", "coordinate_scope", "motif_count", "raw_interval_bp", "union_bp", "interval_overlap_bp",
  "interval_overlap_pct", "genome_bp", "genome_coverage_pct", "assayed_cpg_probes", "cpg_probes_covered", "cpg_probe_coverage_pct"
))

fwrite(summary_dt, file.path(opt$outdir, sprintf("%s.csv", opt$prefix)))
writeLines(c(
  "# G4 Coverage Summary",
  "",
  "Genome coverage is the deduplicated span across hg38 chromosomes 1-22, X, and Y.",
  "CpG probe coverage is the fraction of unique assayed Illumina 450K probes that overlap each annotation.",
  "The `motif_plus_100bp` rows use the current BTEP design: each motif expanded by 100 bp on both sides before deduplicating overlaps."
), file.path(opt$outdir, sprintf("%s.md", opt$prefix)))

print(summary_dt)
