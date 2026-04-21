#' =======================================================================
#' 15b_bvsm_rain.R — B vs m Tradeoff: RAIN only
#' =======================================================================
#'
#' RAIN-only companion to 15_bvsm_method_comparison.R.
#' N capped at 48 to avoid permutation explosion at large m.
#' Results merge with 15_bvsm_method_comparison output for Figure 1.
#'
#' Datasets (single group):
#'   A. Mouse LIV  (GSE54651)  r~2.88  strong
#'   B. Baboon LUN (CAMO)      r~1.72  moderate
#'   C. Mouse D1   (D1D2)      r~0.65  weak / brain
#'
#' USAGE:
#'   Rscript examples/publication/15b_bvsm_rain.R
#'   SMOKE_TEST=true Rscript examples/publication/15b_bvsm_rain.R

# =====================================================================
# 0. Setup
# =====================================================================
SMOKE_TEST  <- identical(Sys.getenv("SMOKE_TEST"), "true")
GLOBAL_SEED <- 2025L
set.seed(GLOBAL_SEED)

old_wd <- setwd("code"); source("setup.R"); setwd(old_wd)

B_VALS     <- c(3L, 4L, 6L, 8L, 12L)
N_GRID     <- if (SMOKE_TEST) c(12L, 24L) else seq(12L, 48L, by = 12L)
ALPHA2_VALS <- 0   # cosinor truth only — sufficient for Figure 1
NSIMS      <- if (SMOKE_TEST) 3L   else 30L
NGENES     <- if (SMOKE_TEST) 100L else 1000L
FDR_THRESH <- 0.05
N_CORES    <- as.integer(Sys.getenv("MC_CORES", unset = "60"))
omega      <- 2 * pi / 24

cat(sprintf("Mode     : %s\n", if (SMOKE_TEST) "SMOKE" else "PRODUCTION"))
cat(sprintf("Method   : RAIN (nr.series)\n"))
cat(sprintf("B vals   : %s\n", paste(B_VALS, collapse = ", ")))
cat(sprintf("N grid   : %s\n", paste(N_GRID,  collapse = ", ")))
cat(sprintf("alpha2   : %s\n", paste(ALPHA2_VALS, collapse = ", ")))
cat(sprintf("nsims    : %d\n", NSIMS))
cat(sprintf("ngenes   : %d\n", NGENES))
cat(sprintf("mc.cores : %d\n\n", N_CORES))

out_dir <- "output/bvsm_method_comparison"
dir.create(file.path(out_dir, "results"), recursive = TRUE, showWarnings = FALSE)

# =====================================================================
# 1. Load pilot datasets — single group from each pair
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
# 2. Estimate pilot parameters per dataset
# =====================================================================
cat("\n--- Estimating pilot parameters ---\n")

pilot_params <- lapply(names(datasets), function(ds_name) {
  ds   <- datasets[[ds_name]]
  fits <- fitCosinorAll(ds$mat, times = ds$tod, period = 24)

  rhy_idx  <- !is.na(fits$pvalue) & fits$pvalue < 0.05
  A_all    <- as.numeric(fits$A)
  sig_all  <- as.numeric(fits$sigma)

  rhy_pool <- which(rhy_idx  & is.finite(A_all)   & A_all > 0 & is.finite(sig_all))
  arr_pool <- which(!rhy_idx & is.finite(sig_all) & sig_all > 0)

  pi_R  <- mean(rhy_idx, na.rm = TRUE)
  n_rhy <- max(10L, round(NGENES * pi_R))
  n_arr <- NGENES - n_rhy

  set.seed(7L)
  g_rhy <- sample(rhy_pool, min(n_rhy, length(rhy_pool)))
  g_arr <- sample(arr_pool, min(n_arr, length(arr_pool)))

  r_med <- median(A_all[g_rhy] / sig_all[g_rhy], na.rm = TRUE)
  cat(sprintf("  %-5s  pi_R=%.0f%%  rhy=%d  arr=%d  median_r=%.2f\n",
              ds_name, 100*pi_R, length(g_rhy), length(g_arr), r_med))

  list(
    A_rhy   = A_all[g_rhy],
    sig_rhy = sig_all[g_rhy],
    sig_arr = sig_all[g_arr],
    n_rhy   = length(g_rhy),
    n_arr   = length(g_arr)
  )
})
names(pilot_params) <- names(datasets)

# =====================================================================
# 3. RAIN power cell
# =====================================================================
power_cell_rain <- function(A_rhy, sig_rhy, sig_arr, B, N, alpha2 = 0) {
  if (N %% B != 0) return(NA_real_)
  m     <- N / B
  cts   <- seq(0, 24 * (1 - 1/B), length.out = B)
  cts_n <- rep(cts, each = m)
  ngr   <- length(A_rhy); nga <- length(sig_arr); ng <- ngr + nga
  gnames <- c(paste0("R", seq_len(ngr)), paste0("A", seq_len(nga)))

  sim_power <- replicate(NSIMS, {
    expr <- matrix(NA_real_, ng, N)
    for (g in seq_len(ngr)) {
      phi <- runif(1, 0, 2*pi)
      mu  <- A_rhy[g] * (cos(omega*cts_n - phi) +
                          alpha2 * cos(2*omega*cts_n - phi))
      expr[g, ] <- rnorm(N, mu, sig_rhy[g])
    }
    for (g in seq_len(nga))
      expr[ngr+g, ] <- rnorm(N, 0, sig_arr[g])
    rownames(expr) <- gnames

    pvals <- tryCatch(
      detect_RAIN(expr, cts_n, gene_names = gnames, period = 24),
      error = function(e) rep(1, ng)
    )
    adj <- p.adjust(pvals, method = "BH")
    sum(adj[seq_len(ngr)] <= FDR_THRESH, na.rm = TRUE) / ngr
  })
  mean(sim_power, na.rm = TRUE)
}

# =====================================================================
# 4. Run grid
# =====================================================================
all_results <- list()

for (ds_name in names(datasets)) {
  pp <- pilot_params[[ds_name]]
  cat(sprintf("\n==============================\n"))
  cat(sprintf("Dataset: %s\n", datasets[[ds_name]]$label))
  cat(sprintf("==============================\n"))

  grid <- expand.grid(N = N_GRID, B = B_VALS)
  grid <- grid[grid$N %% grid$B == 0, ]
  cat(sprintf("  %d valid (N,B) cells\n", nrow(grid)))

  t0 <- proc.time()[["elapsed"]]
  pwr_list <- parallel::mclapply(
    seq_len(nrow(grid)),
    function(i) {
      set.seed(GLOBAL_SEED + i)
      power_cell_rain(pp$A_rhy, pp$sig_rhy, pp$sig_arr,
                      B = grid$B[i], N = grid$N[i], alpha2 = 0)
    },
    mc.cores = N_CORES
  )

  grid$power   <- unlist(pwr_list)
  grid$alpha2  <- 0
  grid$method  <- "RAIN"
  grid$dataset <- ds_name

  elapsed <- proc.time()[["elapsed"]] - t0
  cat(sprintf("  Done in %.0fs\n", elapsed))

  # Print summary
  for (n in N_GRID) {
    row <- grid[grid$N == n, ]
    vals <- sprintf("B%d=%.1f%%", row$B, 100*row$power)
    cat(sprintf("  N=%3d: %s\n", n, paste(vals, collapse="  ")))
  }

  all_results[[ds_name]] <- grid
  saveRDS(grid,
          file.path(out_dir, "results", sprintf("results_RAIN_%s.rds", ds_name)))
  cat(sprintf("Saved: results_RAIN_%s.rds\n", ds_name))
}

# Combined
rain_df <- do.call(rbind, all_results)
saveRDS(rain_df, file.path(out_dir, "results", "results_RAIN_all.rds"))
cat("\nSaved: results_RAIN_all.rds\n")
cat("\n=== Done ===\n")
