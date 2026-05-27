# =============================================================================
# _pilot_2h_liver.R
#
# Build the K=2 two-harmonic pilots used by Figs 4, 5, 6 and SuppS1.
#
#   pilot_2h_GTExLiver.rds      - Passive human pilot for Figs 4, 5, SuppS1
#                                  (GTEx v10 Liver, N=262 after TOD filter,
#                                  Ensembl rows mapped to HGNC symbols.)
#
#   pilot_2h_Hughes2009.rds     - Active mouse pilot for Fig 6
#                                  (Hughes 2009 mouse liver, GSE11923,
#                                  N=48, B=24 hourly active sampling
#                                  CT18..CT65 mod 24.)
#
# Both pilots are fitted via estCircadianParam2H() with the package defaults
# (top_k=500, min_rhythm_pval=0.01), so they share the same calibration knobs
# and the downstream simulator behavior is consistent.
# =============================================================================

setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")

sink(stderr())
source("code/utils.R")
source("code/options.R")
source("code/simulation.R")
source("code/estimation.R")
source("code/bootstrap_sim.R")
source("code/detection.R")
sink()

RES_DIR <- "output/two_harmonic/results"
if (!dir.exists(RES_DIR)) dir.create(RES_DIR, recursive = TRUE)

# ---- GTEx Liver (passive) ----------------------------------------------------
build_gtex_liver <- function() {
  cat("=== GTEx Liver (passive) ===\n")
  env <- new.env()
  load("/home/qtp1/Projects/Collaborative/GTEXdata/data/v10/cpm/v10_CPM_full.RData",
       envir = env)
  td   <- env$CPM.all[["Liver"]]
  lcpm <- td$counts
  ann  <- td$gene_annotation
  tod  <- as.numeric(td$tod) %% 24
  ok   <- !is.na(tod)
  lcpm <- lcpm[, ok, drop = FALSE]; tod <- tod[ok]
  keep <- rowMeans(lcpm) > 3
  lcpm <- lcpm[keep, , drop = FALSE]

  if (!is.null(ann) && "symbol" %in% colnames(ann)) {
    rownames(lcpm) <- ann$symbol[match(rownames(lcpm), ann$gene_id)]
  }
  first <- !duplicated(rownames(lcpm)) & nzchar(rownames(lcpm)) &
           !is.na(rownames(lcpm))
  lcpm <- lcpm[first, , drop = FALSE]

  cat(sprintf("  Expression matrix: %d genes x %d samples\n",
              nrow(lcpm), ncol(lcpm)))
  cat(sprintf("  Unique TOD bins:   %d\n", length(unique(tod))))

  psi <- estCircadianParam2H(data = lcpm, times = tod, period = 24,
                              min_rhythm_pval = 0.01, top_k = 500,
                              verbose = TRUE)
  attr(psi, "pilot_label") <- "GTEx Liver (passive)"
  attr(psi, "n_pilot")     <- ncol(lcpm)
  out <- file.path(RES_DIR, "pilot_2h_GTExLiver.rds")
  saveRDS(psi, out)
  cat(sprintf("  Saved to %s\n\n", out))
  invisible(psi)
}

# ---- Hughes 2009 mouse liver (active, GSE11923) -----------------------------
build_hughes_2009 <- function() {
  cat("=== Hughes 2009 mouse liver (active, GSE11923) ===\n")
  h <- get(load("base/learnRhythmanalysis-master/data/hughes_2009_liver.rda"))
  gene_col <- which(colnames(h) %in% c("geneName", "gene", "Gene", "symbol"))
  expr <- as.matrix(h[, -gene_col])
  storage.mode(expr) <- "numeric"
  rownames(expr) <- h[[gene_col]]
  ct <- as.numeric(sub("^CT", "", colnames(expr)))
  tod <- ct %% 24

  # Filter low-expressed rows
  keep <- rowMeans(expr, na.rm = TRUE) > 5
  expr <- expr[keep, , drop = FALSE]

  cat(sprintf("  Expression matrix: %d genes x %d samples\n",
              nrow(expr), ncol(expr)))
  cat(sprintf("  Unique TOD bins:   %d (each with %d replicates)\n",
              length(unique(tod)), unique(table(tod))))

  psi <- estCircadianParam2H(data = expr, times = tod, period = 24,
                              min_rhythm_pval = 0.01, top_k = 500,
                              verbose = TRUE)
  attr(psi, "pilot_label") <- "Hughes 2009 mouse liver (active B=24)"
  attr(psi, "n_pilot")     <- ncol(expr)
  out <- file.path(RES_DIR, "pilot_2h_Hughes2009.rds")
  saveRDS(psi, out)
  cat(sprintf("  Saved to %s\n\n", out))
  invisible(psi)
}

# ---- Build both --------------------------------------------------------------
psi_liver  <- build_gtex_liver()
psi_hughes <- build_hughes_2009()

cat("Pilot summaries:\n")
cat(sprintf("  GTEx Liver:    n=%d, top_k_used=%d, median r1=%.3f, A2/A1 med=%.3f\n",
            attr(psi_liver, "n_pilot"),
            psi_liver$diagnostics$top_k_used,
            psi_liver$diagnostics$A1_median /
              psi_liver$diagnostics$sigma_median,
            psi_liver$diagnostics$A2_over_A1_med))
cat(sprintf("  Hughes 2009:   n=%d, top_k_used=%d, median r1=%.3f, A2/A1 med=%.3f\n",
            attr(psi_hughes, "n_pilot"),
            psi_hughes$diagnostics$top_k_used,
            psi_hughes$diagnostics$A1_median /
              psi_hughes$diagnostics$sigma_median,
            psi_hughes$diagnostics$A2_over_A1_med))
