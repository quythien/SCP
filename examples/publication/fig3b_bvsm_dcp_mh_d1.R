#' =======================================================================
#' 15d_bvsm_dcp_mh_d1.R — B vs m: DCP + MH for Mouse D1, N up to 96
#' =======================================================================
#'
#' Re-runs DCP and MH for the Mouse D1 dataset at N = 12–96 (matching the
#' N range used in the JTK D1 rerun). The existing results_D1.rds only
#' covers N ≤ 48. This produces an updated D1 results file for merging
#' into results_all.rds.
#'
#' Run AFTER jtk_d1_n96 screen completes to avoid OOM.
#'
#' USAGE:
#'   Rscript examples/publication/15d_bvsm_dcp_mh_d1.R
#'   SMOKE_TEST=true Rscript examples/publication/15d_bvsm_dcp_mh_d1.R

SMOKE_TEST  <- identical(Sys.getenv("SMOKE_TEST"), "true")
GLOBAL_SEED <- 2025L
set.seed(GLOBAL_SEED)

old_wd <- setwd("code"); source("setup.R"); setwd(old_wd)

B_VALS     <- c(3L, 4L, 6L, 8L, 12L)
N_GRID     <- if (SMOKE_TEST) c(12L, 24L, 48L) else seq(12L, 96L, by = 12L)
NSIMS      <- if (SMOKE_TEST) 3L   else 30L
NGENES     <- if (SMOKE_TEST) 200L else 5000L
FDR_THRESH <- 0.05
N_CORES    <- as.integer(Sys.getenv("MC_CORES", unset = "50"))

cat(sprintf("Mode     : %s\n", if (SMOKE_TEST) "SMOKE" else "PRODUCTION"))
cat(sprintf("Methods  : DCP, MH\n"))
cat(sprintf("Dataset  : Mouse D1 (D1D2 striatum)\n"))
cat(sprintf("B vals   : %s\n", paste(B_VALS, collapse = ", ")))
cat(sprintf("N grid   : %s\n", paste(N_GRID,  collapse = ", ")))
cat(sprintf("nsims    : %d\n", NSIMS))
cat(sprintf("ngenes   : %d\n", NGENES))
cat(sprintf("mc.cores : %d\n\n", N_CORES))

out_dir <- "output/bvsm_method_comparison"
dir.create(file.path(out_dir, "results"), recursive = TRUE, showWarnings = FALSE)

# =====================================================================
# Load Mouse D1 pilot
# =====================================================================
cat("--- Loading Mouse D1 ---\n")
pheno   <- read.csv("data/mouse_clinicalinfo_03082021_rmOutliers.csv", row.names = 1)
prep_d1 <- prepCircadianData("data/mouse_D1D2_logCPMfiltered_counts.csv",
                              times      = "time",
                              input_type = "counts",
                              pheno      = pheno,
                              sample_col = "sample")
d1_samp <- pheno$sample[pheno$cell == "D1"]
mat_d1  <- prep_d1$data[, colnames(prep_d1$data) %in% d1_samp, drop = FALSE]
tod_d1  <- pheno$time[match(colnames(mat_d1), pheno$sample)]
mat_d1  <- mat_d1[rowSums(mat_d1 > 1) >= floor(ncol(mat_d1) / 2), , drop = FALSE]
rm(pheno, prep_d1)
cat(sprintf("Mouse D1: %d genes x %d samples\n", nrow(mat_d1), ncol(mat_d1)))

# =====================================================================
# Run DCP + MH
# =====================================================================
bio <- estCircadianParam(mat_d1, times = tod_d1, period = 24, verbose = TRUE)
bio$ngenes <- NGENES

design <- CircadianDesignOptions(
  sample_sizes = N_GRID,
  nsims        = NSIMS,
  design       = "active",
  cts          = seq(0, 24 * (1 - 1/B_VALS[1]), length.out = B_VALS[1]),
  B_values     = B_VALS
)
analysis <- CircadianAnalysisOptions(
  alpha           = FDR_THRESH,
  p.adjust.method = "BH",
  fdr_thresholds  = FDR_THRESH
)

set.seed(GLOBAL_SEED)
res <- runSingleCohortPower(bio, design, analysis,
                             methods  = c("DCP", "MH"),
                             alpha2   = 0,
                             mc.cores = N_CORES,
                             plot     = FALSE,
                             verbose  = TRUE)

saveRDS(res$power_df,
        file.path(out_dir, "results", "results_D1_dcp_mh_n96.rds"))
cat("Saved: results_D1_dcp_mh_n96.rds\n")

# =====================================================================
# Merge with JTK D1 results and rebuild results_D1.rds
# =====================================================================
jtk_file <- file.path(out_dir, "results", "results_D1_jtk_n96.rds")
if (file.exists(jtk_file)) {
  jtk_df <- readRDS(jtk_file)
  combined <- rbind(res$power_df, jtk_df)
  saveRDS(combined, file.path(out_dir, "results", "results_D1.rds"))
  cat("Rebuilt: results_D1.rds (DCP + MH + JTK)\n")
} else {
  cat("Note: JTK results not found yet — results_D1.rds not rebuilt\n")
  cat(sprintf("       Expected at: %s\n", jtk_file))
}

cat("\n=== Done ===\n")
