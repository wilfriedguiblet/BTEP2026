# ARTIFACT CHECK, not a numerical-refitting check. A successful result verifies
# that the paper references the approved tables/images and renders consistently.
# Use replication/run.R --mode models/verify for independent numerical parity.
validate_g4_paper <- function() {
  document_path <- "G4_methylation_paper.Rmd"
  if (!file.exists(document_path)) document_path <- file.path("BTEP", document_path)
  if (!file.exists(document_path)) stop("Run from BTEP or the workspace root")
  previous_directory <- setwd(dirname(normalizePath(document_path)))
  on.exit(setwd(previous_directory), add = TRUE)

  files <- paste0("G4_methylation_paper", c(".Rmd", ".md", ".html", ".bib"))
  stopifnot(all(file.exists(files)), all(file.info(files)$size > 0))
  assets <- data.table::fread("G4_methylation_paper_files/source_manifest.csv")
  # Six selected model tables plus two original PNGs form the evidence bundle.
  # MD5 is used as an identity check, not a cryptographic security guarantee.
  stopifnot(nrow(assets) == 8L, sum(assets$type == "figure") == 2L,
            all(assets$md5 == unname(tools::md5sum(assets$bundled))))
  available_sources <- file.exists(assets$source)
  # A portable checkout may retain only bundled assets, not the ignored results
  # tree. Verify original-source identity too whenever those sources are present.
  stopifnot(all(assets$md5[available_sources] ==
                  unname(tools::md5sum(assets$source[available_sources]))))

  markdown <- paste(readLines("G4_methylation_paper.md", warn = FALSE), collapse = "\n")
  # Scope is a requirement: no old result families or unevaluated inline R may
  # leak into the paper. These exclusions do not prohibit literature discussion.
  forbidden_text <- c("g4_methylation_age_100bp/", "G4Hunter", "chromatin_sensitivity",
                      "cofactor_elastic_net", "# Literature synthesis", "`r ", "[@")
  for (forbidden in forbidden_text) stopifnot(!grepl(forbidden, markdown, fixed = TRUE))
  for (image_path in assets$bundled[assets$type == "figure"]) {
    stopifnot(grepl(image_path, markdown, fixed = TRUE))
  }

  document <- xml2::read_html("G4_methylation_paper.html")
  # Check actual HTML resources rather than assuming image syntax after Pandoc
  # conversion. Both figures and all automatic display dependencies are local
  # data URIs; references remain ordinary links for the reader to follow.
  images <- xml2::xml_attr(xml2::xml_find_all(document, "//img"), "src")
  stopifnot(length(images) == 2L, all(startsWith(images, "data:image/png;base64,")))
  resources <- c(xml2::xml_attr(xml2::xml_find_all(document, "//script[@src]"), "src"),
                 xml2::xml_attr(xml2::xml_find_all(document, "//link[@rel='stylesheet']"), "href"))
  stopifnot(all(startsWith(resources, "data:")))
  styles <- paste(xml2::xml_text(xml2::xml_find_all(document, "//style")), collapse = "\n")
  for (external_url in c("url(http", 'url("http', "url('http", "url(//", 'url("//', "url('//")) {
    stopifnot(!grepl(external_url, styles, fixed = TRUE))
  }

  identifiers <- xml2::xml_attr(xml2::xml_find_all(document, "//*[@id]"), "id")
  # Distinguish bibliography/section anchors from relative file attachments and
  # web links. An existing file is checked; external sites are not re-downloaded.
  references <- xml2::xml_find_all(document, "//*[contains(concat(' ', normalize-space(@class), ' '), ' csl-entry ')]")
  stopifnot(length(references) == 8L, length(xml2::xml_find_all(document, "//table")) == 3L)
  targets <- xml2::xml_attr(xml2::xml_find_all(document, "//a[@href]"), "href")
  fragments <- targets[startsWith(targets, "#")]
  stopifnot(all(substring(fragments, 2L) %in% identifiers))
  local_targets <- targets[!grepl("^(#|[A-Za-z][A-Za-z0-9+.-]*:)", targets)]
  stopifnot(all(file.exists(URLdecode(sub("#.*$", "", local_targets)))))

  message("PASS: original figure hashes, six model tables, two displayed figures, three manuscript tables, eight citations, valid links, no excluded results, offline HTML.")
  print(data.frame(file = files, bytes = file.info(files)$size))
  invisible(TRUE)
}

validate_g4_paper()