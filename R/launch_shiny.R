#' Launch the SCP Shiny app
#'
#' Opens the bundled simulation-based power-analysis app in the user's
#' default browser. The app lets the user pick a bundled pilot (or upload
#' their own), set target power and FDR, read off the recommended sample
#' size, explore per-gene rhythmicity, and run pathway enrichment.
#'
#' On first launch the function checks for the app's optional helper packages
#' and offers to install any that are missing (so the gene-symbol mapping,
#' XLSX export, full-app capture, and pathway-enrichment features work out of
#' the box). These are kept out of the package's hard dependencies so a plain
#' \code{install.packages("SCP")} stays lightweight and the API works without
#' the large Bioconductor annotation databases.
#'
#' @param ... Additional arguments passed to \code{shiny::runApp()}.
#' @param install_deps Logical; if \code{TRUE} (default) install any missing
#'   optional app packages before launching. Set \code{FALSE} to skip the
#'   check (missing features then degrade gracefully inside the app).
#' @param species_all Logical; if \code{TRUE} also install the annotation
#'   databases for rat, zebrafish, Drosophila, and Arabidopsis (in addition to
#'   human + mouse) so upload gene-symbol mapping works for every supported
#'   species. Default \code{FALSE} installs only human + mouse (the common case);
#'   the others install on demand or you can pass \code{TRUE} once.
#'
#' @return Invisibly returns \code{NULL}. Called for its side effect of
#'   launching the Shiny app.
#'
#' @examples
#' \dontrun{
#'   launchShiny()
#'   launchShiny(species_all = TRUE)     # also install rat/zebrafish/fly/plant DBs
#'   launchShiny(install_deps = FALSE)   # skip the dependency check
#' }
#' @export
launchShiny <- function(..., install_deps = TRUE, species_all = FALSE) {
  # Optional packages the app uses. CRAN ones install via install.packages();
  # the Bioconductor annotation DBs via BiocManager. Each feature degrades
  # gracefully inside the app if its package is absent.
  cran_pkgs <- c(
    shiny           = "the app interface (required)",
    DT              = "interactive sortable gene table",
    writexl         = "Excel (.xlsx) gene-table export",
    shinyscreenshot = "full-app screenshot capture"
  )
  # Annotation DBs for upload ID->symbol mapping, one per supported species
  # (mirrors the app's .SPECIES_MAP). Each is large; only the ones whose IDs you
  # upload are needed, so by default we install just the two most common
  # (human, mouse) and let `bioc_pkgs` cover the rest when species_all = TRUE.
  bioc_core <- c(
    org.Hs.eg.db  = "human gene IDs -> symbols (uploads)",
    org.Mm.eg.db  = "mouse gene IDs -> symbols (uploads)",
    clusterProfiler = "pathway enrichment (KEGG/GO) with a custom background",
    ReactomePA    = "Reactome pathway enrichment",
    reactome.db   = "Reactome pathway annotations (ReactomePA dep)"
  )
  bioc_extra <- c(
    org.Rn.eg.db    = "rat gene IDs -> symbols (uploads)",
    org.Dr.eg.db    = "zebrafish gene IDs -> symbols (uploads)",
    org.Dm.eg.db    = "Drosophila (fly) gene IDs -> symbols (uploads)",
    org.Ce.eg.db    = "C. elegans gene IDs -> symbols (uploads)",
    org.Sc.sgd.db   = "S. cerevisiae (yeast) gene IDs -> symbols (uploads)",
    org.Pf.plasmo.db = "P. falciparum gene IDs -> symbols (uploads)",
    org.At.tair.db  = "Arabidopsis gene IDs -> symbols (uploads)"
  )
  bioc_pkgs <- if (isTRUE(species_all)) c(bioc_core, bioc_extra) else bioc_core

  if (isTRUE(install_deps)) {
    missing_cran <- names(cran_pkgs)[!vapply(names(cran_pkgs),
      function(p) requireNamespace(p, quietly = TRUE), logical(1))]
    missing_bioc <- names(bioc_pkgs)[!vapply(names(bioc_pkgs),
      function(p) requireNamespace(p, quietly = TRUE), logical(1))]

    if (length(missing_cran) || length(missing_bioc)) {
      message("SCP app: installing missing optional packages so all features work:")
      for (p in missing_cran) message("  - ", p, "  (", cran_pkgs[[p]], ")")
      for (p in missing_bioc) message("  - ", p, "  (", bioc_pkgs[[p]], "; Bioconductor)")

      if (length(missing_cran)) {
        tryCatch(utils::install.packages(missing_cran),
                 error = function(e)
                   warning("Could not install ", paste(missing_cran, collapse = ", "),
                           ": ", conditionMessage(e), call. = FALSE))
      }
      if (length(missing_bioc)) {
        if (!requireNamespace("BiocManager", quietly = TRUE))
          tryCatch(utils::install.packages("BiocManager"), error = function(e) NULL)
        if (requireNamespace("BiocManager", quietly = TRUE)) {
          tryCatch(BiocManager::install(missing_bioc, update = FALSE, ask = FALSE),
                   error = function(e)
                     warning("Could not install ", paste(missing_bioc, collapse = ", "),
                             ": ", conditionMessage(e), call. = FALSE))
        } else {
          warning("BiocManager unavailable; skipping ", paste(missing_bioc, collapse = ", "),
                  ". Gene-symbol mapping for uploads will fall back to original IDs.",
                  call. = FALSE)
        }
      }
    }
  }

  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Install the `shiny` package first: install.packages(\"shiny\"), ",
         "or call launchShiny() with internet access to auto-install it.")
  }
  appdir <- system.file("shiny", package = "SCP")
  if (!nzchar(appdir)) {
    stop("Shiny app directory not found. Reinstall SCP from a source that ",
         "includes `inst/shiny/app.R`.")
  }
  shiny::runApp(appdir, ...)
}
