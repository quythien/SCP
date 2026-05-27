# =============================================================================
# Build the K=2 two-harmonic pilot for GTEx Adrenal Gland (passive).
#
#   pilot_2h_GTExAdrenal.rds - Passive human pilot for the Fig 5 row 2
#                              (GTEx v10 Adrenal Gland, after TOD filter).
# =============================================================================

setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
sink(stderr())
source("code/utils.R"); source("code/options.R"); source("code/simulation.R")
source("code/estimation.R"); source("code/bootstrap_sim.R"); source("code/detection.R")
sink()

RES_DIR <- "output/two_harmonic/results"
if (!dir.exists(RES_DIR)) dir.create(RES_DIR, recursive = TRUE)

cat("=== GTEx Adrenal Gland (passive) ===\n")
env <- new.env()
load("/home/qtp1/Projects/Collaborative/GTEXdata/data/v10/cpm/v10_CPM_full.RData",
     envir = env)
td   <- env$CPM.all[["Adrenal_Gland"]]
if (is.null(td)) td <- env$CPM.all[["Adrenal Gland"]]
if (is.null(td)) {
  cat("Available tissue names (first 30):\n")
  print(head(names(env$CPM.all), 30))
  stop("Adrenal not found under expected name")
}
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
                            min_rhythm_pval = 0.01, top_k = 300,
                            verbose = TRUE)
attr(psi, "pilot_label") <- "GTEx Adrenal Gland (passive)"
attr(psi, "n_pilot")     <- ncol(lcpm)
out <- file.path(RES_DIR, "pilot_2h_GTExAdrenal_topK300.rds")
saveRDS(psi, out)
cat(sprintf("  Saved to %s\n", out))
cat(sprintf("  A2/A1 median: %.2f\n", psi$diagnostics$A2_over_A1_med))
