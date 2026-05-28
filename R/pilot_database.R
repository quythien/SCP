#' SCP Pilot Dataset Database
#'
#' @description
#' The SCP package bundles a curated database of public circadian
#' transcriptomic pilot datasets that users can plug into
#' \code{runSimsSingleCohort()} or \code{runDifferentialPower()} as the
#' parameter source for power simulation. Each bundled pilot is a
#' pre-fit \code{CircadianBioOptions} object stored under
#' \code{inst/extdata/pilots/<species>/}, alongside a CSV manifest that
#' documents provenance (accession, citation, license), study design
#' (active vs passive), and condition arm.
#'
#' Three user-facing functions are exposed:
#' \describe{
#'   \item{\code{\link{scp_pilots}}}{Return the manifest as a data.frame,
#'         optionally filtered by species and ingestion status.}
#'   \item{\code{\link{scp_load_pilot}}}{Load and return one standardized
#'         pilot object given species, dataset, tissue, and condition.}
#'   \item{\code{\link{scp_pilot_search}}}{Convenience search: filter the
#'         manifest by free-text query, species, tissue regex, and
#'         in-vivo vs cell-line sample type.}
#' }
#'
#' @name scp_pilot_database
NULL


# -------------------------------------------------------------------------
# Internal: locate manifest and pilots directory
# -------------------------------------------------------------------------

#' @keywords internal
.scp_pilots_dir <- function() {
  # Prefer installed package extdata; fall back to source-tree inst/extdata.
  p <- system.file("extdata", "pilots", package = "SCP")
  if (nzchar(p) && dir.exists(p)) return(p)
  src <- file.path("inst", "extdata", "pilots")
  if (dir.exists(src)) return(normalizePath(src))
  stop("SCP pilot database not found. Reinstall the SCP package.")
}

#' @keywords internal
.scp_manifest_path <- function() {
  file.path(.scp_pilots_dir(), "manifest.csv")
}


# -------------------------------------------------------------------------
# scp_pilots: load the manifest
# -------------------------------------------------------------------------

#' List Bundled SCP Pilot Datasets
#'
#' @description
#' Reads the bundled pilot manifest (\code{inst/extdata/pilots/manifest.csv})
#' and returns it as a data.frame. Optionally filters by species and
#' ingestion status so a user can quickly enumerate the pilots they
#' actually have available.
#'
#' @param species Optional character. One or more of
#'   \code{"human"}, \code{"baboon"}, \code{"other_primate"},
#'   \code{"mouse"}, \code{"rat"}. \code{NULL} (default) returns all species.
#' @param status Character. Filter on the \code{status} column. Defaults
#'   to \code{"ingested"} (datasets actually present on disk). Use
#'   \code{NULL} or \code{"all"} to return every row, including
#'   \code{"pending"}, \code{"verify_failed"}, \code{"controlled"}, and
#'   \code{"ingest_failed"} entries.
#'
#' @return A data.frame with one row per (dataset, tissue, condition)
#'   pilot. Columns: \code{species}, \code{dataset}, \code{tissue},
#'   \code{condition}, \code{n}, \code{ngenes}, \code{design},
#'   \code{assay}, \code{sample_type}, \code{accession}, \code{citation},
#'   \code{license}, \code{file}, \code{status}.
#'
#' @examples
#' \dontrun{
#'   scp_pilots()
#'   scp_pilots(species = "human")
#'   scp_pilots(status = "all")
#' }
#' @export
scp_pilots <- function(species = NULL, status = "ingested") {
  m <- utils::read.csv(.scp_manifest_path(), stringsAsFactors = FALSE)
  if (!is.null(species)) m <- m[m$species %in% species, , drop = FALSE]
  if (!is.null(status) && !identical(status, "all")) {
    m <- m[m$status %in% status, , drop = FALSE]
  }
  rownames(m) <- NULL
  m
}


# -------------------------------------------------------------------------
# scp_load_pilot
# -------------------------------------------------------------------------

#' Load a Bundled SCP Pilot Dataset
#'
#' @description
#' Loads one standardized pilot object from the bundled pilot database.
#' Each pilot is a pre-fit \code{CircadianBioOptions} S3 object (see
#' \code{\link{estimate_circadian_params}}) and can be passed directly
#' to \code{runSimsSingleCohort()} or \code{runDifferentialPower()}.
#'
#' @param species Character, length 1. One of \code{"human"},
#'   \code{"baboon"}, \code{"other_primate"}, \code{"mouse"},
#'   \code{"rat"}.
#' @param dataset Character, length 1. Manifest \code{dataset} field
#'   (typically GEO accession or short label, e.g. \code{"GTEx"},
#'   \code{"GSE160521"}).
#' @param tissue Character, length 1. Tissue or region label as recorded
#'   in the manifest (e.g. \code{"Liver"}, \code{"NAc"}).
#' @param condition Character, length 1, optional. Condition arm
#'   (e.g. \code{"Control"}, \code{"All"}). If \code{NULL} and the
#'   (species, dataset, tissue) triple identifies more than one row,
#'   the function errors and prints the available conditions.
#'
#' @return A \code{CircadianBioOptions} object with empirical
#'   distributions of mesor, log-sigma, amplitude, phase, and proportion
#'   rhythmic, ready for downstream simulation.
#'
#' @examples
#' \dontrun{
#'   bio <- scp_load_pilot("human", "GTEx", "Liver")
#'   bio <- scp_load_pilot("human", "GSE160521", "NAc", "Control")
#' }
#' @export
scp_load_pilot <- function(species, dataset, tissue, condition = NULL) {
  stopifnot(length(species) == 1, length(dataset) == 1, length(tissue) == 1)
  m <- scp_pilots(species = species, status = "ingested")
  hit <- m$dataset == dataset & m$tissue == tissue
  if (!is.null(condition)) hit <- hit & m$condition == condition
  m <- m[hit, , drop = FALSE]
  if (nrow(m) == 0L) {
    stop(sprintf("No ingested pilot matches species='%s' dataset='%s' tissue='%s'%s",
                 species, dataset, tissue,
                 if (is.null(condition)) "" else sprintf(" condition='%s'", condition)))
  }
  if (nrow(m) > 1L) {
    stop(sprintf("Multiple pilots match; please pass `condition` (available: %s)",
                 paste(sQuote(m$condition), collapse = ", ")))
  }
  path <- file.path(.scp_pilots_dir(), m$file[1])
  if (!file.exists(path)) {
    stop(sprintf("Manifest lists '%s' but the file is missing on disk.", m$file[1]))
  }
  p <- readRDS(path)
  # Harmonize amplitude / phase to match sigma_rhythmic length when the pilot
  # was built with mismatched per-slot gene sets (caught in several agent-
  # ingested GTEx pilots: amplitude held the full-threshold set while
  # sigma_rhythmic was top-K). The simulation pipeline pairs amplitude[g] with
  # sigma_rhythmic[g] per gene, so length mismatch produces NaN p-values.
  if (!is.null(p$amplitude) && !is.null(p$sigma_rhythmic)) {
    n_sigma <- length(p$sigma_rhythmic)
    if (n_sigma > 0L) {
      if (length(p$amplitude) != n_sigma)
        p$amplitude <- if (length(p$amplitude) > n_sigma)
          p$amplitude[seq_len(n_sigma)] else rep_len(p$amplitude, n_sigma)
      if (!is.null(p$phase) && length(p$phase) != n_sigma)
        p$phase <- if (length(p$phase) > n_sigma)
          p$phase[seq_len(n_sigma)] else rep_len(p$phase, n_sigma)
    }
  }
  if (!inherits(p, "CircadianBioOptions"))
    class(p) <- c("CircadianBioOptions", class(p))
  p
}


# -------------------------------------------------------------------------
# scp_pilot_search
# -------------------------------------------------------------------------

#' Search the SCP Pilot Database
#'
#' @description
#' Convenience filter over the pilot manifest. Combines free-text
#' substring matching against dataset, tissue, citation, and accession
#' fields with optional species and tissue-regex filters and an
#' in-vivo vs cell-line restriction. Useful when a user wants to pick
#' a pilot whose design (tissue, sample size, signal strength) matches
#' their planned study.
#'
#' @param query Character, optional. Case-insensitive substring matched
#'   against \code{dataset}, \code{tissue}, \code{accession},
#'   \code{citation}. \code{NULL} matches everything.
#' @param species Character, optional. Restrict to these species.
#' @param tissue_pattern Character, optional. Regular expression matched
#'   against \code{tissue} (case-insensitive).
#' @param sample_type Character. Default \code{"in_vivo"}. Use
#'   \code{NULL} to include cell-line pilots as well.
#'
#' @return A data.frame: rows of the manifest passing every active filter.
#'
#' @examples
#' \dontrun{
#'   scp_pilot_search("liver")
#'   scp_pilot_search(species = "mouse", tissue_pattern = "brain|cortex")
#'   scp_pilot_search(sample_type = NULL)   # include cell lines
#' }
#' @export
scp_pilot_search <- function(query = NULL,
                             species = NULL,
                             tissue_pattern = NULL,
                             sample_type = "in_vivo") {
  m <- scp_pilots(species = species, status = "ingested")
  if (!is.null(sample_type)) m <- m[m$sample_type %in% sample_type, , drop = FALSE]
  if (!is.null(tissue_pattern)) {
    hit <- grepl(tissue_pattern, m$tissue, ignore.case = TRUE)
    if ("tissue_canonical" %in% names(m))
      hit <- hit | grepl(tissue_pattern, m$tissue_canonical, ignore.case = TRUE)
    m <- m[hit, , drop = FALSE]
  }
  if (!is.null(query)) {
    hay <- paste(m$dataset, m$tissue, m$accession, m$citation, sep = " | ")
    if ("tissue_canonical" %in% names(m))
      hay <- paste(hay, m$tissue_canonical, sep = " | ")
    m <- m[grepl(query, hay, ignore.case = TRUE), , drop = FALSE]
  }
  rownames(m) <- NULL
  m
}


# -------------------------------------------------------------------------
# scp_pilot_tod
# -------------------------------------------------------------------------

#' Get the time-of-day (TOD) vector for a bundled SCP pilot
#'
#' Returns the per-sample sampling times for a pilot, in hours mod 24.
#' Prefers the raw \code{p$times} field when the pilot rds preserved it;
#' otherwise parses the manifest \code{tod_phases} column and returns the
#' unique sampling phases.
#'
#' @param species,dataset,tissue,condition manifest tuple identifying the
#'   pilot.
#' @return Numeric vector of times in hours (mod 24); \code{NULL} if no
#'   TOD information is available for this pilot.
#' @examples
#' \dontrun{
#'   tod <- scp_pilot_tod("baboon", "GSE98965", "SCN", "All")
#'   hist(tod, breaks = seq(0, 24, by = 1))
#' }
#' @export
scp_pilot_tod <- function(species, dataset, tissue, condition = NULL) {
  p <- tryCatch(scp_load_pilot(species, dataset, tissue, condition),
                error = function(e) NULL)
  if (!is.null(p)) {
    tt <- p$times %||% p$raw$times %||% NULL
    if (!is.null(tt) && length(tt) > 0L) return(as.numeric(tt))
  }
  m <- scp_pilots(species = species, status = "ingested")
  hit <- m$dataset == dataset & m$tissue == tissue
  if (!is.null(condition)) hit <- hit & m$condition == condition
  m <- m[hit, , drop = FALSE]
  if (nrow(m) == 0L || !("tod_phases" %in% names(m))) return(NULL)
  ph <- m$tod_phases[1]
  if (is.na(ph) || !nzchar(ph) || grepl("phases\\)$", ph)) return(NULL)
  suppressWarnings(as.numeric(strsplit(ph, ",")[[1]]))
}
