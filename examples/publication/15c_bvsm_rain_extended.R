#' =======================================================================
#' 15c_bvsm_rain_extended.R — RAIN B vs m: Extended N + B + Waveform Violation
#' =======================================================================
#'
#' Extends 15b_bvsm_rain.R with:
#'   - N up to 96 to cross the 80% power threshold for all datasets
#'   - Extended B: c(3, 4, 6, 8, 12, 18, 24) — full range of active designs
#'   - alpha2 in {0, 0.5, 1.0} to show waveform violation effect
#'
#' Datasets (single group):
#'   A. Mouse LIV  (GSE54651)  r~2.88  strong
#'   B. Baboon LUN (CAMO)      r~1.72  moderate
#'   C. Mouse D1   (D1D2)      r~0.65  weak / brain  ← primary for Figure 4
#'
#' Provides results for:
#'   - Figure 3 (alpha2=0, all datasets, full N range to 80% crossing)
#'   - Figure 4 (RAIN + MH, alpha2 variation, D1 primary)
#'
#' Note: only valid (N, B) cells where N %% B == 0 are run. B=18 requires
#' N in {18,36,54,72,90}; B=24 requires N in {24,48,72,96}.
#'
#' USAGE:
#'   Rscript examples/publication/15c_bvsm_rain_extended.R
#'   SMOKE_TEST=true Rscript examples/publication/15c_bvsm_rain_extended.R

# =====================================================================
# 0. Setup
# =====================================================================
SMOKE_TEST  <- identical(Sys.getenv("SMOKE_TEST"), "true")
GLOBAL_SEED <- 2025L
set.seed(GLOBAL_SEED)

old_wd <- setwd("code"); source("setup.R"); setwd(old_wd)

B_VALS      <- c(3L, 4L, 6L, 8L, 12L, 18L, 24L)
N_GRID      <- if (SMOKE_TEST) c(12L, 24L, 48L) else seq(12L, 96L, by = 12L)
ALPHA2_VALS <- if (SMOKE_TEST) c(0, 1.0) else c(0, 0.5, 0.75, 1.0)
NSIMS       <- if (SMOKE_TEST) 3L   else 30L
NGENES      <- if (SMOKE_TEST) 200L else 5000L
FDR_THRESH  <- 0.05
N_CORES     <- as.integer(Sys.getenv("MC_CORES", unset = "60"))

cat(sprintf("Mode     : %s\n", if (SMOKE_TEST) "SMOKE" else "PRODUCTION"))
cat(sprintf("Method   : RAIN\n"))
cat(sprintf("B vals   : %s\n", paste(B_VALS,      collapse = ", ")))
cat(sprintf("N grid   : %s\n", paste(N_GRID,       collapse = ", ")))
cat(sprintf("alpha2   : %s\n", paste(ALPHA2_VALS,  collapse = ", ")))
cat(sprintf("nsims    : %d\n", NSIMS))
cat(sprintf("ngenes   : %d\n", NGENES))
cat(sprintf("mc.cores : %d\n\n", N_CORES))

# Valid (N, B) cells
valid_grid <- expand.grid(N = N_GRID, B = B_VALS)
valid_grid <- valid_grid[valid_grid$N %% valid_grid$B == 0, ]
cat(sprintf("Valid (N,B) cells: %d\n\n", nrow(valid_grid)))

out_dir <- "output/bvsm_method_comparison"
dir.create(file.path(out_dir, "results"), recursive = TRUE, showWarnings = FALSE)

# =====================================================================
# 1. Load pilot datasets
# =====================================================================
cat("--- Loading datasets ---\n")

## A: Mouse LIV (GSE54651)
dat_mouse <- readRDS("data/mice_GSE54651_CPM.RData")
prep_liv  <- prepCircadianData(dat_mouse$count_clean[["LIV"]],
                               times      = dat_mouse$tod[["LIV"]],
                               input_type = "log2")
mat_liv  <- prep_liv$data[rowSums(prep_liv$data > 0) >= 4, , drop = FALSE]
tod_liv  <- prep_liv$times
rm(dat_mouse, prep_liv)

## B: Baboon LUN (CAMO)
load("data/CAMO_PRC_hmb.RData")
prep_lun <- prepCircadianData(baboon_withTOD$baboon[["LUN"]],
                               times      = baboon_withTOD$tod[["LUN"]],
                               input_type = "cpm")
mat_lun  <- prep_lun$data[rowSums(prep_lun$data > 0) >= 6, , drop = FALSE]
tod_lun  <- prep_lun$times
rm(baboon_withTOD, gtex, mice, prep_lun)

## C: Mouse D1 (D1D2 striatum)
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

datasets <- list(
  LIV = list(mat = mat_liv, tod = tod_liv, label = "Mouse LIV (r~2.88)"),
  LUN = list(mat = mat_lun, tod = tod_lun, label = "Baboon LUN (r~1.72)"),
  D1  = list(mat = mat_d1,  tod = tod_d1,  label = "Mouse D1 (r~0.65)")
)
rm(mat_liv, mat_lun, mat_d1, tod_liv, tod_lun, tod_d1)

# =====================================================================
# 2. Run RAIN power with extended N, B, and alpha2 sweep
# =====================================================================
cat("\n--- Running RAIN extended (N<=96, B<=24, alpha2 sweep) ---\n")

all_results <- list()

for (ds_name in names(datasets)) {
  ds <- datasets[[ds_name]]
  cat(sprintf("\n==============================\n"))
  cat(sprintf("Dataset: %s\n", ds$label))
  cat(sprintf("==============================\n"))

  bio <- estCircadianParam(ds$mat, times = ds$tod, period = 24, verbose = TRUE)
  bio$ngenes <- NGENES   # use full matrix for estimation, cap simulation gene count

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
                               methods     = "RAIN",
                               alpha2      = ALPHA2_VALS,
                               mc.cores    = N_CORES,
                               plot        = FALSE,
                               verbose     = TRUE)

  all_results[[ds_name]] <- res$power_df

  # Save per alpha2 slice
  for (a2 in ALPHA2_VALS) {
    sub <- res$power_df[res$power_df$alpha2 == a2, ]
    fname <- sprintf("results_RAIN_ext_%s_a2_%.1f.rds", ds_name, a2)
    saveRDS(sub, file.path(out_dir, "results", fname))
    cat(sprintf("  Saved: %s\n", fname))

    # Quick power summary for alpha2=0
    if (a2 == 0) {
      cat(sprintf("  Power summary (alpha2=0):\n"))
      for (n in sort(unique(sub$N))) {
        row <- sub[sub$N == n, ]
        row <- row[order(row$B), ]
        vals <- sprintf("B%d=%.1f%%", row$B, 100 * row$power)
        cat(sprintf("    N=%3d: %s\n", n, paste(vals, collapse = "  ")))
      }
    }
  }

  # Combined across alpha2 for this dataset
  saveRDS(res$power_df,
          file.path(out_dir, "results", sprintf("results_RAIN_ext_%s_all.rds", ds_name)))
}

# Full combined table
rain_ext_df <- do.call(rbind, all_results)
saveRDS(rain_ext_df, file.path(out_dir, "results", "results_RAIN_extended_all.rds"))
cat(sprintf("\nSaved: results_RAIN_extended_all.rds  (%d rows)\n", nrow(rain_ext_df)))
cat("\n=== Done ===\n")
