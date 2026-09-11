# DOCUMENT REPLAY ONLY: this script does not refit age slopes or G4 models.
# For numerical re-estimation, use replication/run.R first. The Rmd bundles
# the retained figures/tables and uses their values for prose and tables.
# It intentionally does not replace the original PNGs with refitted plots.
render_g4_paper <- function() {
  document_path <- "G4_methylation_paper.Rmd"
  if (!file.exists(document_path)) document_path <- file.path("BTEP", document_path)
  if (!file.exists(document_path)) stop("Run from BTEP or the workspace root")
  document_path <- normalizePath(document_path)
  previous_directory <- setwd(dirname(document_path))
  on.exit(setwd(previous_directory), add = TRUE)

  required_packages <- c("rmarkdown", "knitr", "data.table", "xml2")
  missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                                quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing_packages)) stop("Missing packages: ", paste(missing_packages, collapse = ", "))
  # Both outputs execute the Rmd's formula, coefficient/contrast, cohort, and
  # BH-family assertions. 'all' means HTML and GitHub Markdown from its YAML.
  rmarkdown::render(basename(document_path), output_format = "all", quiet = TRUE)

  html_path <- "G4_methylation_paper.html"
  document <- xml2::read_html(html_path)
  # R Markdown may exempt KaTeX CDN resources from self_contained embedding.
  # Remove that exemption, then use Pandoc's supported embedding pass. Network
  # access is needed here, but the completed HTML is readable offline.
  math_dependencies <- xml2::xml_find_all(document, "//script[@src] | //link[@rel='stylesheet']")
  xml2::xml_attr(math_dependencies, "data-external") <- NULL
  temporary_html <- tempfile(pattern = "g4_offline_", tmpdir = ".", fileext = ".html")
  on.exit(unlink(temporary_html), add = TRUE)
  xml2::write_html(document, temporary_html)
  rmarkdown::pandoc_self_contained_html(temporary_html, html_path)

  rendered <- xml2::read_html(html_path)
  # Outbound citation links are allowed. Display resources (figures, styles,
  # scripts/fonts) must be embedded rather than fetched when a reader opens it.
  resources <- c(xml2::xml_attr(xml2::xml_find_all(rendered, "//script[@src] | //img[@src]"), "src"),
                 xml2::xml_attr(xml2::xml_find_all(rendered, "//link[@rel='stylesheet']"), "href"))
  stopifnot(all(startsWith(resources, "data:")),
            length(xml2::xml_find_all(rendered, "//img")) == 2L)
  message("Rendered G4_methylation_paper.md and offline G4_methylation_paper.html")
  invisible(html_path)
}

render_g4_paper()