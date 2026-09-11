#!/usr/bin/env Rscript

# PURPOSE: import the deposited, already processed GEO beta matrix and metadata.
# This is NOT raw-IDAT normalization. Exact historical normalization, probe QC,
# and donor exclusions must be checked against the source study separately.
# The paper uses GSE61257/adipose, GSE61258/liver, and GSE61259/muscle only.
# Use --series-matrix to replay a frozen local download; otherwise this script
# downloads the current GEO file. The caller should record its checksum.

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

option_list <- list(
  make_option("--accession", type = "character"),
  make_option("--tissue", type = "character"),
  make_option("--prefix", type = "character"),
  make_option("--series-matrix", type = "character", default = "", dest = "series_matrix"),
  make_option("--outdir", type = "character", default = ".")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$accession) || is.null(opt$prefix)) {
  stop("--accession and --prefix are required")
}
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

series_prefix <- sub("[0-9]{3}$", "nnn", opt$accession)
url <- sprintf(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/%s_series_matrix.txt.gz",
  series_prefix, opt$accession, opt$accession
)
dest <- if (nzchar(opt$series_matrix)) opt$series_matrix else file.path(opt$outdir, basename(url))
if (nzchar(opt$series_matrix)) {
  if (!file.exists(dest)) stop("Local series matrix not found: ", dest)
} else {
  download.file(url, dest, mode = "wb", quiet = TRUE)
}

# A GEO series matrix is not a plain table from its first line: sample-level
# attributes precede the marked beta table. Keep accession order until the
# explicit sample-ID join below. Do not infer ordering from donor labels.
lines <- readLines(gzfile(dest), warn = FALSE)
sample_line <- grep("^!Sample_geo_accession", lines, value = TRUE)[1]
if (is.na(sample_line)) stop("No !Sample_geo_accession line found in ", opt$accession)
sample_ids <- gsub('^"|"$', '', strsplit(sample_line, "\t")[[1]][-1])

metadata <- data.table(sample_id = sample_ids, tissue = opt$tissue)
metadata_lines <- lines[grep("^!Sample_(title|source_name_ch1|characteristics_ch1)", lines)]
# The source series stores age, sex, BMI and subject labels in characteristics.
# Review the generated columns after import: malformed numeric values become NA,
# and later first-stage complete-case filtering can then remove specimens.
# A shared subject label is not automatically proof of an independent donor or
# a replicate; do not silently deduplicate the published metadata here.
for (line in metadata_lines) {
  parts <- strsplit(line, "\t")[[1]]
  key <- sub("^!Sample_", "", parts[1])
  values <- gsub('^"|"$', '', parts[-1])
  if (length(values) != length(sample_ids)) next
  if (key == "title") {
    metadata[, title := values]
  } else if (key == "source_name_ch1") {
    metadata[, source_name := values]
  } else if (key == "characteristics_ch1") {
    split_values <- tstrsplit(values, ":", fixed = TRUE, keep = 1:2)
    field <- tolower(trimws(split_values[[1]][1]))
    field <- gsub("[^a-z0-9_]+", "_", field)
    field <- sub("_+$", "", field)
    if (!nzchar(field)) next
    metadata[, (field) := trimws(sub("^[^:]+:", "", values))]
  }
}

for (column in intersect(c("age", "bmi", "dnamage"), names(metadata))) {
  metadata[, (column) := suppressWarnings(as.numeric(get(column)))]
}
if ("sex" %in% names(metadata)) {
  metadata[, sex := tolower(sex)]
}

begin <- grep("^!series_matrix_table_begin", lines)
end <- grep("^!series_matrix_table_end", lines)
# Preserve the deposited beta values. No cell deconvolution, age correction,
# imputation, clipping, M-value transformation, or new normalization occurs.
if (length(begin) != 1 || length(end) != 1 || end <= begin) {
  stop("Could not find series matrix table boundaries in ", opt$accession)
}
matrix_lines <- lines[(begin + 1):(end - 1)]
beta <- fread(text = paste(matrix_lines, collapse = "\n"), data.table = FALSE)
probe_col <- names(beta)[1]
rownames(beta) <- beta[[probe_col]]
beta[[probe_col]] <- NULL
names(beta) <- gsub('^"|"$', '', names(beta))
beta <- as.matrix(beta)
storage.mode(beta) <- "numeric"

common_samples <- intersect(colnames(beta), metadata$sample_id)
# Align by exact GSM/sample ID, not merely by column position. Both artifacts
# are written in identical sample order so a human can inspect the pairing.
if (length(common_samples) < 2) {
  stop("Fewer than two samples matched between matrix and metadata for ", opt$accession)
}
beta <- beta[, common_samples, drop = FALSE]
metadata <- metadata[match(common_samples, sample_id)]

saveRDS(beta, file.path(opt$outdir, sprintf("%s_beta.rds", opt$prefix)))
fwrite(metadata, file.path(opt$outdir, sprintf("%s_metadata.csv", opt$prefix)))
