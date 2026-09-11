#!/usr/bin/env Rscript

# Run without arguments to CHECK dependencies. Use --install explicitly to
# install missing ones. No existing package is automatically upgraded, and no
# global .Rprofile or repository setting is written. A project-specific R
# library is preferable for long-term reproduction (set R_LIBS_USER beforehand).
# The hg38 BSgenome package is large; the check-only default never downloads it.
arguments <- commandArgs(trailingOnly = TRUE)
if (any(!arguments %in% c("--install", "--record-environment"))) {
  stop("Supported flags: --install, --record-environment")
}
cran_packages <- c("data.table", "ggplot2", "optparse", "yaml", "rmarkdown", "knitr", "xml2")
bioconductor_packages <- c("IlluminaHumanMethylation450kanno.ilmn12.hg19",
                          "BSgenome.Hsapiens.UCSC.hg38", "GenomicRanges", "IRanges",
                          "GenomeInfoDb", "Biostrings", "BSgenome")
installed <- function(package) requireNamespace(package, quietly = TRUE)

if ("--install" %in% arguments) {
  missing_cran <- cran_packages[!vapply(cran_packages, installed, logical(1))]
  if (length(missing_cran)) install.packages(missing_cran, repos = "https://cloud.r-project.org")
  missing_bioc <- bioconductor_packages[!vapply(bioconductor_packages, installed, logical(1))]
  if (length(missing_bioc)) {
    if (!installed("BiocManager")) install.packages("BiocManager", repos = "https://cloud.r-project.org")
    # Let BiocManager choose the release compatible with the current R. This
    # installs a usable environment, NOT a claim to recover the original one.
    # Compare with environment_validated.csv and rerun the numerical gate.
    BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
  }
}

packages <- c(cran_packages, bioconductor_packages)
versions <- vapply(packages, function(package) {
  if (installed(package)) as.character(packageVersion(package)) else NA_character_
}, "")
environment <- data.frame(
  package = c("R", "Pandoc", packages),
  version = c(as.character(getRversion()),
              if (installed("rmarkdown")) as.character(rmarkdown::pandoc_version()) else NA_character_,
              versions),
  platform = R.version$platform,
  stringsAsFactors = FALSE
)
print(environment, row.names = FALSE)
if (anyNA(versions)) stop("Some packages are missing; rerun with --install after reviewing the download requirements")

# This optional manifest records the environment used NOW. It is not a package
# lockfile and does not certify historical package versions. The scientific
# runner writes a second per-run version record beside its outputs.
if ("--record-environment" %in% arguments) {
  script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
  script_path <- normalizePath(sub("^--file=", "", script_argument))
  write.csv(environment, file.path(dirname(script_path), "environment_validated.csv"), row.names = FALSE)
}
message("PASS: analysis and document packages are available; no package installation occurs unless --install was supplied.")