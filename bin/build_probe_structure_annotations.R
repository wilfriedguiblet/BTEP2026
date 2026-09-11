#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

option_list <- list(
  make_option("--manifest", type = "character"),
  make_option("--ucsc-goldenpath-dir", type = "character", default = "/fdb/genomebrowser/goldenPath/hg38", dest = "ucsc_goldenpath_dir"),
  make_option("--rloop-bed", type = "character", default = "", dest = "rloop_bed"),
  make_option("--rloop-url", type = "character", default = "", dest = "rloop_url"),
  make_option("--replication-timing-bed", type = "character", default = "", dest = "replication_timing_bed"),
  make_option("--replication-timing-url", type = "character", default = "", dest = "replication_timing_url"),
  make_option("--cross-reactive-probes", type = "character", default = "", dest = "cross_reactive_probes"),
  make_option("--cross-reactive-probes-url", type = "character", default = "", dest = "cross_reactive_probes_url"),
  make_option("--outdir", type = "character", default = ".")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$manifest)) stop("--manifest is required")
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

manifest <- fread(opt$manifest)
probes <- manifest[, .(probe_id, chr = chr_hg38, start = as.integer(start_hg38), end = as.integer(end_hg38))]
probes <- unique(probes[grepl("^chr([0-9]+|X|Y)$", chr) & is.finite(start)], by = "probe_id")
probes[, probe_index := .I]

overlap_any <- function(intervals) {
  if (!nrow(intervals)) return(rep(FALSE, nrow(probes)))
  query <- copy(probes)
  setkeyv(query, c("chr", "start", "end"))
  setkeyv(intervals, c("chr", "start", "end"))
  hits <- foverlaps(query, intervals, nomatch = 0L)
  out <- rep(FALSE, nrow(probes))
  out[unique(hits$probe_index)] <- TRUE
  out
}

download_if_needed <- function(url, filename) {
  path <- file.path(opt$outdir, filename)
  if (!file.exists(path)) download.file(url, path, mode = "wb", quiet = TRUE)
  path
}

download_first_available <- function(sources) {
  for (source in sources) {
    if (!is.null(source$local_path) && file.exists(source$local_path)) {
      return(c(path = source$local_path, format = source$format, url = paste0("local:", source$local_path)))
    }
    path <- file.path(opt$outdir, source$filename)
    ok <- file.exists(path) || isTRUE(tryCatch({
      download.file(source$url, path, mode = "wb", quiet = TRUE)
      TRUE
    }, error = function(e) FALSE))
    if (ok) return(c(path = path, format = source$format, url = source$url))
  }
  stop("Could not download an hg38 RepeatMasker annotation from UCSC. Provide a local RepeatMasker BED/TSV or restore network access.")
}

resolve_optional_input <- function(path, url, filename) {
  if (nzchar(path) && file.exists(path)) return(path)
  if (nzchar(url)) return(download_if_needed(url, filename))
  ""
}

rmsk_source <- download_first_available(list(
  list(
    local_path = file.path(opt$ucsc_goldenpath_dir, "database", "rmsk.txt.gz"),
    filename = "hg38_rmsk.txt.gz",
    format = "ucsc_table"
  ),
  list(
    url = "https://hgdownload.soeucsc.edu/goldenPath/hg38/database/rmsk.txt.gz",
    filename = "hg38_rmsk.txt.gz",
    format = "ucsc_table"
  ),
  list(
    url = "https://hgdownload.soeucsc.edu/goldenPath/hg38/bigZips/hg38.fa.out.gz",
    filename = "hg38.fa.out.gz",
    format = "repeatmasker_out"
  )
))
if (rmsk_source[["format"]] == "ucsc_table") {
  rmsk <- fread(rmsk_source[["path"]], header = FALSE, select = c(6, 7, 8, 11, 12))
  setnames(rmsk, c("chr", "start", "end", "repeat_class", "repeat_family"))
} else {
  rmsk <- fread(rmsk_source[["path"]], skip = 3, header = FALSE, fill = TRUE, select = c(5, 6, 7, 11))
  setnames(rmsk, c("chr", "start", "end", "repeat_class"))
  rmsk[, repeat_family := repeat_class]
  rmsk[, repeat_class := sub("/.*$", "", repeat_class)]
}
rmsk <- rmsk[grepl("^chr([0-9]+|X|Y)$", chr)]
rmsk[, repeat_class := fifelse(repeat_class %in% c("Simple_repeat", "Low_complexity", "Satellite"), repeat_class,
                         fifelse(repeat_class %in% c("LINE", "SINE", "LTR", "DNA"), repeat_class, "Other_repeat"))]

query <- copy(probes)
setkeyv(query, c("chr", "start", "end"))
setkeyv(rmsk, c("chr", "start", "end"))
repeat_hits <- foverlaps(query, rmsk, nomatch = 0L)
repeat_class <- data.table(probe_index = probes$probe_index, repeat_class = "No_repeat")
if (nrow(repeat_hits)) {
  priority <- c("Satellite", "Simple_repeat", "Low_complexity", "LINE", "SINE", "LTR", "DNA", "Other_repeat")
  repeat_hits[, repeat_rank := match(repeat_class, priority)]
  best_repeat <- repeat_hits[order(probe_index, repeat_rank), .SD[1], by = probe_index]
  repeat_class[best_repeat$probe_index, repeat_class := best_repeat$repeat_class]
}

segdup_source <- download_first_available(list(
  list(
    local_path = file.path(opt$ucsc_goldenpath_dir, "database", "genomicSuperDups.txt.gz"),
    filename = "hg38_genomicSuperDups.txt.gz",
    format = "ucsc_table"
  ),
  list(
    url = "https://hgdownload.soeucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz",
    filename = "hg38_genomicSuperDups.txt.gz",
    format = "ucsc_table"
  )
))
segdup <- fread(segdup_source[["path"]], header = FALSE, select = 2:4)
setnames(segdup, c("chr", "start", "end"))
segdup <- segdup[grepl("^chr([0-9]+|X|Y)$", chr)]

annotations <- cbind(
  probes[, .(probe_id, chr, start, end)],
  repeat_class[, .(repeat_class)],
  segmental_duplication = overlap_any(segdup),
  rloop_overlap = FALSE,
  replication_timing = NA_real_,
  cross_reactive = FALSE
)

rloop_path <- resolve_optional_input(opt$rloop_bed, opt$rloop_url, "optional_rloop.bed")
timing_path <- resolve_optional_input(opt$replication_timing_bed, opt$replication_timing_url, "optional_replication_timing.bed")
cross_path <- resolve_optional_input(opt$cross_reactive_probes, opt$cross_reactive_probes_url, "optional_cross_reactive_probes.txt")

if (nzchar(rloop_path)) {
  rloop <- fread(rloop_path, header = FALSE, select = 1:3)
  setnames(rloop, c("chr", "start", "end"))
  annotations[, rloop_overlap := overlap_any(rloop)]
}
if (nzchar(timing_path)) {
  timing <- fread(timing_path, header = FALSE, select = 1:4)
  setnames(timing, c("chr", "start", "end", "replication_timing"))
  query <- copy(probes)
  setkeyv(query, c("chr", "start", "end"))
  setkeyv(timing, c("chr", "start", "end"))
  timing_hits <- foverlaps(query, timing, nomatch = 0L)
  if (nrow(timing_hits)) {
    timing_mean <- timing_hits[, .(replication_timing = mean(as.numeric(replication_timing), na.rm = TRUE)), by = probe_index]
    annotations[timing_mean$probe_index, replication_timing := timing_mean$replication_timing]
  }
}
if (nzchar(cross_path)) {
  cross_ids <- fread(cross_path, header = FALSE)[[1]]
  annotations[, cross_reactive := probe_id %in% cross_ids]
}

fwrite(annotations, file.path(opt$outdir, "probe_structure_annotations.csv.gz"))
qc <- data.table(
  metric = c("n_probes", "repeat_overlap_fraction", "segmental_duplication_fraction", "rloop_overlap_fraction", "replication_timing_nonmissing_fraction", "cross_reactive_fraction"),
  value = c(
    nrow(annotations),
    mean(annotations$repeat_class != "No_repeat"),
    mean(annotations$segmental_duplication),
    mean(annotations$rloop_overlap),
    mean(is.finite(annotations$replication_timing)),
    mean(annotations$cross_reactive)
  )
)
fwrite(qc, file.path(opt$outdir, "probe_structure_annotation_qc.csv"))
if (!nrow(segdup)) stop("Segmental-duplication input has no canonical hg38 intervals; check source column mapping")
provenance <- data.table(
  annotation = c("RepeatMasker", "segmental_duplication", "rloop_overlap", "replication_timing", "cross_reactive"),
  source = c(rmsk_source[["url"]], segdup_source[["url"]], rloop_path, timing_path, cross_path),
  available = c(TRUE, TRUE, nzchar(rloop_path), nzchar(timing_path), nzchar(cross_path))
)
fwrite(provenance, file.path(opt$outdir, "probe_structure_annotation_provenance.csv"))
