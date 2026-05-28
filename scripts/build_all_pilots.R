#' Canonical pilot build pipeline for SCP
#'
#' Documents EXACTLY how every bundled pilot rds in inst/extdata/pilots/
#' should be (re)built so the package's `scp_load_pilot()` results match
#' what the SCP manuscript uses. All three species follow the same
#' two-step pattern:
#'
#'   1. Load a tissue-specific normalized expression matrix + per-sample
#'      time-of-day vector from the species-specific source file.
#'   2. Call `estCircadianParam(expr, times, period = 24)` — the package
#'      function defined in R/estimation.R:351. The resulting object is a
#'      `CircadianBioOptions` list with both summary scalars AND the full
#'      per-gene raw fit table at `p$raw` (M, A, phi, sigma, pvalue,
#'      is_rhythmic, in_estim_set). The sim-facing amplitude/sigma_rhythmic
#'      /phase vectors are populated automatically by estCircadianParam at
#'      its internal threshold (default min_rhythm_pval = 0.01) — there is
#'      no extra filtering step in this script.
#'
#' Source files used (all live on the lab server):
#'
#' GTEx
#'   /home/qtp1/Projects/Collaborative/GTEXdata/TOD Xiangning/data/CPM/CPM.all.norm.RData
#'   Object: CPM.all.norm = list of 54 normalized CPM matrices keyed by
#'           full human-readable tissue name (e.g. "Adrenal Gland", "Liver").
#'   TOD encoding: column names are
#'           GTEX.<donor>.<HHMM>.<sample-meta>, e.g.
#'           "GTEX.11DXY.0526.SM.5EGGQ" -> 05:26 -> 5.43 h.
#'   This is the source that produced the manuscript pilots:
#'           data/gtex_Liver_single_pilot.rds, data/gtex_AdrenalGland_single_pilot.rds
#'   DO NOT use GTEX_Filtered.RDS for this purpose — that file is the
#'   pre-computed per-gene cosinor fit table, NOT the raw expression.
#'
#' Baboon
#'   data/CAMO_PRC_hmb.RData
#'   Object: baboon_withTOD$baboon[[tissue]] (4938-gene CPM-like matrix,
#'           columns named <tissue>_ZT<HH>), baboon_withTOD$tod[[tissue]]
#'           = per-column TOD vector.
#'
#' Mouse GSE54651
#'   data/mice_GSE54651_CPM.RData
#'   Object: list with $count_clean[[tissue]] (log2 CPM matrix) and
#'           $tod[[tissue]] (per-sample TOD in hours).
#'
#' Usage:
#'   Rscript scripts/build_all_pilots.R [species]
#'   species in {gtex, baboon, mouse, all}; default = all.

suppressPackageStartupMessages({
  old_wd <- setwd("code"); source("setup.R"); setwd(old_wd)
})

args <- commandArgs(trailingOnly = TRUE)
which_sp <- if (length(args) == 0) "all" else tolower(args[1])

to_target_name <- function(nm) {
  out <- gsub(" - ", "_", nm, fixed = TRUE)
  out <- gsub("\\(|\\)", "", out)
  out <- gsub(" ", "_", out)
  out <- gsub("__+", "_", out)
  out <- gsub("_+$", "", out)
  out
}

build_one <- function(expr, times, out_path) {
  if (is.null(expr) || ncol(expr) < 10) return(FALSE)
  set.seed(2025)   # GLOBAL_SEED used by every manuscript figure script
  p <- tryCatch(estCircadianParam(expr, times, period = 24, verbose = FALSE),
                error = function(e) { cat("    ERR:", conditionMessage(e), "\n"); NULL })
  if (is.null(p)) return(FALSE)
  p$times <- times
  if (!inherits(p, "CircadianBioOptions"))
    class(p) <- c("CircadianBioOptions", class(p))
  saveRDS(p, out_path)
  TRUE
}

# ============================================================================
# GTEx
# ============================================================================
build_gtex <- function(manifest) {
  src <- "/home/qtp1/Projects/Collaborative/GTEXdata/TOD Xiangning/data/CPM/CPM.all.norm.RData"
  cat(sprintf("\n[GTEx] Loading %s\n", src))
  e <- new.env(); load(src, envir = e)
  cpm <- get(ls(e)[1], e)

  extract <- function(df) {
    if (is.null(df) || nrow(df) == 0L || ncol(df) == 0L) return(NULL)
    ids  <- as.character(colnames(df))
    hhmm <- sapply(strsplit(ids, "\\."),
                   function(x) if (length(x) >= 3) x[3] else NA)
    hrs  <- as.numeric(substr(hhmm, 1, 2)) +
            as.numeric(substr(hhmm, 3, 4)) / 60
    ok   <- !is.na(hrs)
    if (sum(ok) < 10L) return(NULL)
    list(expr = as.matrix(df[, ok]), times = hrs[ok])
  }

  fixed <- 0L
  for (tn in names(cpm)) {
    ts <- extract(cpm[[tn]])
    if (is.null(ts)) { cat(sprintf("  skip %s (no usable samples)\n", tn)); next }
    target <- to_target_name(tn)
    fname  <- if (target == "Adrenal_Gland") "GTEx_Adrenal_All.rds" else
              if (target == "Liver")        "GTEx_Liver_All.rds"   else
              sprintf("GTEx_%s_All.rds", target)
    out_path <- file.path("inst/extdata/pilots/human", fname)
    if (build_one(ts$expr, ts$times, out_path)) {
      base <- if (target == "Adrenal_Gland") "Adrenal" else target
      rows <- which(manifest$species == "human" &
                    manifest$dataset == "GTEx" &
                    manifest$tissue %in% c(target, base))
      if (length(rows) > 0L) {
        manifest$status[rows] <- "ingested"; manifest$note[rows] <- NA
        manifest$n[rows]      <- length(ts$times)
        p <- readRDS(out_path)
        manifest$ngenes[rows] <- p$ngenes
      }
      fixed <- fixed + 1L
    }
  }
  cat(sprintf("[GTEx] Rebuilt %d tissues.\n", fixed))
  manifest
}

# ============================================================================
# Baboon GSE98965 (Mure 2018)
# ============================================================================
build_baboon <- function(manifest) {
  cat("\n[baboon] Loading data/CAMO_PRC_hmb.RData\n")
  e <- new.env(); load("data/CAMO_PRC_hmb.RData", envir = e)
  bb <- get("baboon_withTOD", envir = e)
  fixed <- 0L
  for (tn in names(bb$baboon)) {
    expr  <- as.matrix(bb$baboon[[tn]]); mode(expr) <- "numeric"
    times <- as.numeric(bb$tod[[tn]])
    if (length(times) == 0L || ncol(expr) != length(times)) next
    out_path <- file.path("inst/extdata/pilots/baboon",
                          sprintf("GSE98965_%s_All.rds", tn))
    if (build_one(expr, times, out_path)) {
      rows <- which(manifest$species == "baboon" &
                    manifest$dataset == "GSE98965" &
                    manifest$tissue == tn)
      if (length(rows) > 0L) {
        manifest$status[rows] <- "ingested"; manifest$note[rows] <- NA
        manifest$n[rows]      <- length(times)
        p <- readRDS(out_path)
        manifest$ngenes[rows] <- p$ngenes
      }
      fixed <- fixed + 1L
    }
  }
  cat(sprintf("[baboon] Rebuilt %d tissues.\n", fixed))
  manifest
}

# ============================================================================
# Mouse GSE54651 (Zhang 2014 RNA-seq atlas)
# ============================================================================
build_mouse <- function(manifest) {
  cat("\n[mouse GSE54651] Loading data/mice_GSE54651_CPM.RData\n")
  d <- readRDS("data/mice_GSE54651_CPM.RData")
  fixed <- 0L
  for (tn in names(d$count_clean)) {
    times <- d$tod[[tn]]
    if (length(times) == 0L) next
    expr  <- as.matrix(d$count_clean[[tn]]); mode(expr) <- "numeric"
    if (ncol(expr) != length(times)) next
    out_path <- file.path("inst/extdata/pilots/mouse",
                          sprintf("GSE54651_%s_All.rds", tn))
    if (build_one(expr, as.numeric(times), out_path)) {
      rows <- which(manifest$species == "mouse" &
                    manifest$dataset == "GSE54651" &
                    manifest$tissue == tn)
      if (length(rows) > 0L) {
        manifest$status[rows] <- "ingested"; manifest$note[rows] <- NA
        manifest$n[rows]      <- length(times)
        p <- readRDS(out_path)
        manifest$ngenes[rows] <- p$ngenes
      }
      fixed <- fixed + 1L
    }
  }
  cat(sprintf("[mouse GSE54651] Rebuilt %d tissues.\n", fixed))
  manifest
}

# ============================================================================
# Run
# ============================================================================
mp <- "inst/extdata/pilots/manifest.csv"
m  <- read.csv(mp, stringsAsFactors = FALSE)
if (which_sp %in% c("gtex",   "all")) m <- build_gtex(m)
if (which_sp %in% c("baboon", "all")) m <- build_baboon(m)
if (which_sp %in% c("mouse",  "all")) m <- build_mouse(m)
write.csv(m, mp, row.names = FALSE)
cat(sprintf("\nFinal ingested rows in manifest: %d\n",
            sum(m$status == "ingested")))
