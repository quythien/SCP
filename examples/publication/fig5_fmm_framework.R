#' =======================================================================
#' fig5_fmm_framework.R — Full FMM framework on baboon LUN (3-panel data)
#' =======================================================================
#'
#' Generates the simulation data for the 3-panel Fig 5:
#'   Panel A — marginal power vs N, lines per BH-FDR alpha
#'   Panel B — stratified power by r-tilde (r-strata)
#'   Panel C — eta sweep (omega_g ~ Beta(1, eta))
#'
#' Pipeline (Section 2.5 of paper):
#'   FMM-estimated pilot params (estCircadianParamFMM)
#'   -> FMM-generated data
#'   -> K = 2 harmonic F-test detection (method = "FMM", K = 2)
#'   -> raw pvalues retained for offline BH re-thresholding.
#'
#' Panels A and B share a single simulation (varying N at eta = eta_hat);
#' Panel C uses a separate sim with eta varied.
#'
#' OUTPUT:
#'   output/fmm_framework/results/fig5_fmm_framework_<ts>.rds
#'
#' USAGE:
#'   SMOKE_TEST=true Rscript examples/publication/fig5_fmm_framework.R
#'   Rscript examples/publication/fig5_fmm_framework.R

SMOKE_TEST  <- identical(Sys.getenv("SMOKE_TEST"), "true")
NGENES      <- if (SMOKE_TEST) 500L  else 2000L
NSIMS       <- if (SMOKE_TEST) 5L    else 20L
N_CORES     <- as.integer(Sys.getenv("MC_CORES", unset = "20"))
TOP_K_FMM   <- if (SMOKE_TEST) 50L   else 300L
RHYTHM_PVAL <- 0.01
GLOBAL_SEED <- 2025L

# Sample-size grid (Panels A and B share this)
N_GRID    <- if (SMOKE_TEST) c(24L, 48L) else
             c(12L, 24L, 36L, 48L, 60L, 72L, 96L, 120L, 144L)
# FDR thresholds for Panel A (one curve per alpha)
FDR_GRID  <- c(0.01, 0.05, 0.10, 0.20)
# eta grid for Panel C (one curve per eta)
ETA_GRID  <- if (SMOKE_TEST) c(0, 1, 5)  else c(0, 0.5, 1, 2, 5, 20)

cat(sprintf("Mode      : %s\n", if (SMOKE_TEST) "SMOKE" else "PRODUCTION"))
cat(sprintf("NGENES    : %d\n", NGENES))
cat(sprintf("NSIMS     : %d\n", NSIMS))
cat(sprintf("N_GRID    : %s\n", paste(N_GRID,   collapse = ", ")))
cat(sprintf("ETA_GRID  : %s\n", paste(ETA_GRID, collapse = ", ")))
cat(sprintf("FDR_GRID  : %s\n", paste(FDR_GRID, collapse = ", ")))
cat(sprintf("MC_CORES  : %d\n", N_CORES))

old_wd <- setwd("code"); source("setup.R"); setwd(old_wd)

out_dir_res <- "output/fmm_framework/results"
dir.create(out_dir_res, recursive = TRUE, showWarnings = FALSE)

# ====================================================================
# 1. Baboon LUN pilot + cached FMM fit (anchor for all three panels)
# ====================================================================
cat("\n=== Step 1: load Baboon LUN + FMM fit ===\n")
load("data/CAMO_PRC_hmb.RData")
prep_lun <- prepCircadianData(baboon_withTOD$baboon[["LUN"]],
                               times = baboon_withTOD$tod[["LUN"]] %% 24,
                               input_type = "cpm")
mat_lun  <- prep_lun$data[rowSums(prep_lun$data > 0) >= 6, , drop = FALSE]
tod_lun  <- prep_lun$times
rm(baboon_withTOD, gtex, mice, prep_lun)

set.seed(GLOBAL_SEED)
g_idx   <- sample(nrow(mat_lun), min(NGENES, nrow(mat_lun)))
mat_sub <- mat_lun[g_idx, ]

bio_rds <- sprintf("output/sensitivity/results/bio_lun_fmm_NG%d_K%d.rds",
                   NGENES, TOP_K_FMM)
if (file.exists(bio_rds) && !identical(Sys.getenv("REFIT"), "true")) {
  cat(sprintf("Loading cached FMM fit: %s\n", bio_rds))
  bio_lun <- readRDS(bio_rds)
} else {
  cat("Fitting FMM per top-K rhythmic gene...\n")
  bio_lun <- estCircadianParamFMM(mat_sub, tod_lun,
                                   min_rhythm_pval = RHYTHM_PVAL,
                                   top_k    = TOP_K_FMM,
                                   mc.cores = N_CORES,
                                   verbose  = TRUE)
  bio_lun$ngenes <- NGENES
  saveRDS(bio_lun, bio_rds)
}

eta_hat         <- bio_lun$diagnostics$beta_hat
sigma_alpha_hat <- bio_lun$diagnostics$sigma_alpha_hat
cat(sprintf("LUN anchor: eta_hat=%.3f  sigma_alpha_hat=%.3fh\n",
            eta_hat, sigma_alpha_hat))

# Active design B = 12
B_VAL <- 12L
cts_active <- seq(0, 24 * (1 - 1/B_VAL), length.out = B_VAL)

# Anchor bio.opts: FMM truth using empirical eta_hat
set_fmm_truth <- function(bio_base, eta, sd_hours) {
  bio_k <- bio_base
  bio_k$omega_rhythmic <- NULL
  bio_k$alpha_rhythmic <- NULL
  if (eta == 0) {
    bio_k$omega_dist <- list(family = "fixed", value = 1.0)
  } else {
    bio_k$omega_dist <- list(family = "beta", a = 1, b = eta)
  }
  bio_k$alpha_dist <- list(family = "normal", mean = 0, sd_hours = sd_hours)
  bio_k
}

bio_anchor <- set_fmm_truth(bio_lun, eta = eta_hat, sd_hours = sigma_alpha_hat)

# r-strata for Panel B (adaptive, anchored to the pilot)
r_strata    <- makeAdaptiveRStrata(bio_lun, bin_width = 0.25)
strata_lbls <- paste0("(",
                      round(head(r_strata, -1), 2), ", ",
                      round(tail(r_strata, -1), 2), "]")
cat(sprintf("r_strata: %s\n",
            paste(round(r_strata[is.finite(r_strata)], 2), collapse = ", ")))

# ====================================================================
# 2. Panels A and B: shared single-pilot sim at eta = eta_hat
#    Run once at FDR = 0.05; pvalues retained for offline re-thresholding.
# ====================================================================
cat("\n=== Step 2: N sweep at eta = eta_hat (full FMM, K = 2) ===\n")
analysis_AB <- CircadianAnalysisOptions(alpha            = 0.05,
                                         p.adjust.method = "BH",
                                         fdr_thresholds  = FDR_GRID,
                                         r_strata        = r_strata,
                                         strata_labels   = strata_lbls)

# We assemble the sim cell-by-cell so each N gets its own cts_design.
# The runner returns pvalues per (N, gene, sim) and r_values per replicate.
G <- bio_anchor$ngenes
n_N <- length(N_GRID)
n_sims <- NSIMS
n_strata <- length(strata_lbls)
n_fdr <- length(FDR_GRID)

# Storage for re-thresholded power
power_marg <- array(NA_real_, dim = c(n_N, n_fdr, n_sims))
power_strat <- array(NA_real_, dim = c(n_N, n_strata, n_fdr, n_sims))

for (j in seq_along(N_GRID)) {
  N <- N_GRID[j]
  if (N %% B_VAL != 0L) next
  m <- N %/% B_VAL
  cts_design <- rep(cts_active, each = m)
  design <- CircadianDesignOptions(sample_sizes = N, nsims = n_sims,
                                    design = "active",
                                    cts = cts_design, B_values = B_VAL)
  set.seed(GLOBAL_SEED)
  res <- runSimsSingleCohort(bio_anchor, design, analysis_AB,
                              method   = "FMM", K = 2L,
                              mc.cores = N_CORES, verbose = FALSE)
  # res$pvalues is [1 x G x nsims] for this single-N call
  pvals <- res$pvalues[1, , ]                     # G x nsims
  rvals_list <- res$r_values_list[[1]]            # list of length nsims

  for (s in seq_len(n_sims)) {
    p_s <- pvals[, s]
    r_s <- rvals_list[[s]]
    is_rh <- r_s > 0
    # r-stratum assignment for rhythmic genes
    r_for_strat <- ifelse(is_rh, r_s, 0)
    xgr <- cut(r_for_strat, breaks = r_strata, include.lowest = TRUE, labels = FALSE)
    xgr[!is_rh] <- NA

    for (a_i in seq_len(n_fdr)) {
      fdr_g <- p.adjust(p_s, method = "BH")
      disc  <- fdr_g <= FDR_GRID[a_i]
      n_rh  <- sum(is_rh)
      TD_total <- sum(disc & is_rh, na.rm = TRUE)
      power_marg[j, a_i, s] <- if (n_rh > 0) TD_total / n_rh else NA_real_

      for (k in seq_len(n_strata)) {
        in_k <- !is.na(xgr) & xgr == k
        n_tgt <- sum(is_rh & in_k)
        TD_k  <- sum(disc & is_rh & in_k, na.rm = TRUE)
        power_strat[j, k, a_i, s] <- if (n_tgt > 0) TD_k / n_tgt else NA_real_
      }
    }
  }
  cat(sprintf("  N=%3d: power @0.05 = %.3f\n", N,
              mean(power_marg[j, which(FDR_GRID == 0.05), ], na.rm = TRUE)))
}

# ====================================================================
# 3. Panel C: eta sweep at FDR = 0.05
# ====================================================================
cat("\n=== Step 3: eta sweep at FDR = 0.05 (K = 2) ===\n")
analysis_C <- CircadianAnalysisOptions(alpha           = 0.05,
                                        p.adjust.method = "BH",
                                        fdr_thresholds  = 0.05)
power_eta <- matrix(NA_real_, nrow = length(N_GRID), ncol = length(ETA_GRID),
                    dimnames = list(paste0("N=", N_GRID),
                                    paste0("eta=", ETA_GRID)))
for (k in seq_along(ETA_GRID)) {
  bio_k <- set_fmm_truth(bio_lun, eta = ETA_GRID[k], sd_hours = sigma_alpha_hat)
  pwr_curve <- vapply(N_GRID, function(N) {
    if (N %% B_VAL != 0L) return(NA_real_)
    m <- N %/% B_VAL
    cts_design <- rep(cts_active, each = m)
    design <- CircadianDesignOptions(sample_sizes = N, nsims = NSIMS,
                                      design = "active",
                                      cts = cts_design, B_values = B_VAL)
    set.seed(GLOBAL_SEED)
    res <- runSimsSingleCohort(bio_k, design, analysis_C,
                                method   = "FMM", K = 2L,
                                mc.cores = N_CORES, verbose = FALSE)
    mean(res$marginal_power, na.rm = TRUE)
  }, numeric(1))
  power_eta[, k] <- pwr_curve
  cat(sprintf("  eta=%5.2f: %s\n", ETA_GRID[k],
              paste(sprintf("%.2f", pwr_curve), collapse = " ")))
}

# ====================================================================
# 4. Save
# ====================================================================
out <- list(
  pilot           = "baboon LUN (active, GSE98965)",
  detector        = "FMM (K = 2 harmonic F-test)",
  framework       = "full FMM-based SCP (Section 2.5)",
  GLOBAL_SEED     = GLOBAL_SEED,
  NGENES          = NGENES,
  NSIMS           = NSIMS,
  N_GRID          = N_GRID,
  ETA_GRID        = ETA_GRID,
  FDR_GRID        = FDR_GRID,
  r_strata        = r_strata,
  strata_labels   = strata_lbls,
  eta_hat         = eta_hat,
  sigma_alpha_hat = sigma_alpha_hat,
  # Panel A: marginal power [N x FDR x sim]; mean over sims for plotting.
  power_marg      = power_marg,
  # Panel B: stratified power [N x stratum x FDR x sim]
  power_strat     = power_strat,
  # Panel C: marginal power [N x eta] at FDR = 0.05
  power_eta       = power_eta
)
ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
rds <- file.path(out_dir_res, sprintf("fig5_fmm_framework_%s.rds", ts))
saveRDS(out, rds)
cat(sprintf("\nSaved: %s\n", rds))
