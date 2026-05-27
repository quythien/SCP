#' =======================================================================
#' fig6_phase_mse.R — Active design phase-MSE under cosinor and K-harmonic
#' =======================================================================
#'
#' Driver script for Fig 6 Panel B (phase MSE) and supporting analysis
#' for Panels A and C if cached results are unavailable.
#'
#' KEY QUESTION:
#' Does increasing sampling density B (number of distinct time points)
#' reduce per-gene phase estimation error MSE(\hat\phi_g), at fixed
#' total sample size N? This separates the "B improves detection" claim
#' (which is approximately FALSE above the identifiability threshold)
#' from the "B improves phase estimation" claim (which is TRUE).
#'
#' DESIGN:
#'   - Active sampling grid: B in {4, 6, 8, 12, 24}; total N in {48, 96}.
#'   - For each (B, N), simulate FMM-truth data on baboon KIM pilot,
#'     run detect_DCP for cosinor phase estimate and detect_FMM(K=2) for
#'     K-harmonic phase estimate.
#'   - Extract \hat\phi_g per rhythmic gene, compute circular MSE against
#'     FMM truth \alpha_g (converted from radians to hours).
#'   - Aggregate over NSIMS replicate simulations.
#'
#' OUTPUT:
#'   output/phase_mse/results/fig6_phase_mse_<ts>.rds
#'   output/phase_mse/figures/fig6_phase_mse.pdf  (built by replot script)
#'
#' USAGE:
#'   SMOKE_TEST=true  Rscript examples/publication/fig6_phase_mse.R
#'   Rscript examples/publication/fig6_phase_mse.R
#'
#' REPRODUCIBILITY:
#'   - GLOBAL_SEED fixed at 2025L.
#'   - Each (B, N) cell reseeds before its inner sim loop.
#'   - All parameters (NGENES, NSIMS, B_GRID, N_GRID, pilot) are surfaced
#'     at the top of this file. No silent magic.

SMOKE_TEST  <- identical(Sys.getenv("SMOKE_TEST"), "true")
NGENES      <- if (SMOKE_TEST) 500L  else 2000L
NSIMS       <- if (SMOKE_TEST) 3L    else 100L
N_CORES     <- as.integer(Sys.getenv("MC_CORES", unset = "20"))
TOP_K_FMM   <- if (SMOKE_TEST) 50L   else 300L
RHYTHM_PVAL <- 0.01
GLOBAL_SEED <- 2025L
PERIOD      <- 24

# Sampling-grid and sample-size sweeps
B_GRID <- if (SMOKE_TEST) c(4L, 12L)       else c(4L, 6L, 8L, 12L, 24L)
N_GRID <- if (SMOKE_TEST) c(24L)           else c(24L, 48L, 72L, 96L, 120L, 144L)

cat(sprintf("Mode      : %s\n", if (SMOKE_TEST) "SMOKE" else "PRODUCTION"))
cat(sprintf("NGENES    : %d\n", NGENES))
cat(sprintf("NSIMS     : %d\n", NSIMS))
cat(sprintf("MC_CORES  : %d\n", N_CORES))
cat(sprintf("B_GRID    : %s\n", paste(B_GRID, collapse = ", ")))
cat(sprintf("N_GRID    : %s\n", paste(N_GRID, collapse = ", ")))

old_wd <- setwd("code"); source("setup.R"); setwd(old_wd)

out_dir_res <- "output/phase_mse/results"
dir.create(out_dir_res, recursive = TRUE, showWarnings = FALSE)

# ====================================================================
# 1. Load baboon KIM pilot + fit FMM
# ====================================================================
cat("\n=== Step 1: load + FMM-fit baboon KIM pilot ===\n")
load("data/CAMO_PRC_hmb.RData")
prep <- prepCircadianData(baboon_withTOD$baboon[["KIM"]],
                          times = baboon_withTOD$tod[["KIM"]] %% 24,
                          input_type = "cpm")
mat_kim <- prep$data[rowSums(prep$data > 0) >= 6, , drop = FALSE]
tod_kim <- prep$times
rm(baboon_withTOD, gtex, mice, prep)

set.seed(GLOBAL_SEED)
g_idx   <- sample(nrow(mat_kim), min(NGENES, nrow(mat_kim)))
mat_sub <- mat_kim[g_idx, ]

bio_rds <- sprintf("output/active_vs_passive/results/bio_baboon_KIM_NG%d_K%d.rds",
                   NGENES, TOP_K_FMM)
if (file.exists(bio_rds) && !identical(Sys.getenv("REFIT"), "true")) {
  cat(sprintf("Loading cached FMM fit: %s\n", bio_rds))
  bio_kim <- readRDS(bio_rds)
} else {
  cat("Fitting FMM per top-K rhythmic gene...\n")
  bio_kim <- estCircadianParamFMM(mat_sub, tod_kim,
                                   min_rhythm_pval = RHYTHM_PVAL,
                                   top_k    = TOP_K_FMM,
                                   mc.cores = N_CORES,
                                   verbose  = TRUE)
  bio_kim$ngenes <- NGENES
  dir.create(dirname(bio_rds), recursive = TRUE, showWarnings = FALSE)
  saveRDS(bio_kim, bio_rds)
}

# ====================================================================
# 2. Phase-MSE helper
# ====================================================================
# Circular distance in hours on the 24h cycle:
#   d(a, b) = min(|a - b|, period - |a - b|)
# MSE on circular metric: mean(d^2) over rhythmic genes.
circular_diff <- function(phi_hat, phi_truth, period = 24) {
  d <- abs(phi_hat - phi_truth) %% period
  pmin(d, period - d)
}

# Per-(B, N) phase-MSE evaluation under cosinor (DCP) and K-harmonic (K=2).
# Reports both mean and median MSE; also restricts to high-SNR rhythmic genes
# (true r_g >= R_FLOOR) since low-SNR genes give essentially random phi_hat
# regardless of sampling density and dominate the mean MSE.
R_FLOOR <- 1.0

run_one_cell <- function(B, N) {
  if (N %% B != 0L) return(NULL)
  m <- N %/% B
  cts_active <- rep(seq(0, PERIOD * (1 - 1/B), length.out = B), each = m)

  mse_dcp_all   <- numeric(NSIMS); mse_dcp_hi   <- numeric(NSIMS)
  mse_kfmm_all  <- numeric(NSIMS); mse_kfmm_hi  <- numeric(NSIMS)
  med_dcp_all   <- numeric(NSIMS); med_dcp_hi   <- numeric(NSIMS)
  med_kfmm_all  <- numeric(NSIMS); med_kfmm_hi  <- numeric(NSIMS)

  for (s in seq_len(NSIMS)) {
    set.seed(GLOBAL_SEED + 1000L * B + 7L * N + s)
    sim <- simCircadianFMM(bio_kim, cts_active, omega = 1.0, beta = pi)

    expr   <- sim$expr
    is_rh  <- as.logical(sim$is_rhythmic)
    # With simulator default beta = pi, FMM "alpha" is the trough time;
    # the actual peak (which cosinor phi_hat estimates) is at alpha + pi
    # radians, i.e., alpha_hours + 12 hours.
    alpha_truth_rad <- sim$alpha_g                      # radians, length n_rhythmic
    phi_truth_hr    <- (alpha_truth_rad * (PERIOD / (2 * pi)) + PERIOD / 2) %% PERIOD
    r_truth         <- sim$r_values[is_rh]               # length n_rhythmic
    hi_snr          <- !is.na(r_truth) & r_truth >= R_FLOOR

    # --- cosinor detection (DCP) ---
    cos_fit <- fitCosinorAll_fast(expr, cts_active, period = PERIOD)
    phi_hat_cos <- cos_fit$phi[is_rh]
    d_cos       <- circular_diff(phi_hat_cos, phi_truth_hr, PERIOD)
    mse_dcp_all[s] <- mean(d_cos^2, na.rm = TRUE)
    med_dcp_all[s] <- median(d_cos^2, na.rm = TRUE)
    mse_dcp_hi[s]  <- mean(d_cos[hi_snr]^2, na.rm = TRUE)
    med_dcp_hi[s]  <- median(d_cos[hi_snr]^2, na.rm = TRUE)

    # --- K-harmonic (K=2): direct OLS for first-harmonic phase ---
    omega0 <- 2 * pi / PERIOD
    X_k2 <- cbind(1,
                  cos(omega0 * cts_active),    sin(omega0 * cts_active),
                  cos(2 * omega0 * cts_active), sin(2 * omega0 * cts_active))
    XtX_inv <- tryCatch(solve(crossprod(X_k2)),
                        error = function(e) NULL)
    if (is.null(XtX_inv)) {
      mse_kfmm_all[s] <- NA_real_; med_kfmm_all[s] <- NA_real_
      mse_kfmm_hi[s]  <- NA_real_; med_kfmm_hi[s]  <- NA_real_
    } else {
      coef_k2 <- tcrossprod(expr %*% X_k2, XtX_inv)  # G x 5
      phi_hat_k2 <- (atan2(coef_k2[, 3], coef_k2[, 2]) / omega0) %% PERIOD
      d_k2 <- circular_diff(phi_hat_k2[is_rh], phi_truth_hr, PERIOD)
      mse_kfmm_all[s] <- mean(d_k2^2, na.rm = TRUE)
      med_kfmm_all[s] <- median(d_k2^2, na.rm = TRUE)
      mse_kfmm_hi[s]  <- mean(d_k2[hi_snr]^2, na.rm = TRUE)
      med_kfmm_hi[s]  <- median(d_k2[hi_snr]^2, na.rm = TRUE)
    }
  }

  # Per-replicate vectors retained so the replot script can compute
  # 95% Monte Carlo bands (2.5% / 97.5% quantiles across replicates).
  list(B = B, N = N, R_FLOOR = R_FLOOR,
       # Per-replicate raw vectors (length NSIMS)
       mse_dcp_all_raw   = mse_dcp_all,
       mse_dcp_hi_raw    = mse_dcp_hi,
       mse_kfmm_all_raw  = mse_kfmm_all,
       mse_kfmm_hi_raw   = mse_kfmm_hi,
       med_dcp_all_raw   = med_dcp_all,
       med_dcp_hi_raw    = med_dcp_hi,
       med_kfmm_all_raw  = med_kfmm_all,
       med_kfmm_hi_raw   = med_kfmm_hi,
       # Summaries used by older plots
       mse_dcp_all_mean  = mean(mse_dcp_all,  na.rm = TRUE),
       mse_dcp_hi_mean   = mean(mse_dcp_hi,   na.rm = TRUE),
       mse_kfmm_all_mean = mean(mse_kfmm_all, na.rm = TRUE),
       mse_kfmm_hi_mean  = mean(mse_kfmm_hi,  na.rm = TRUE),
       mse_dcp_all_se    = sd(mse_dcp_all,    na.rm = TRUE) / sqrt(NSIMS),
       mse_dcp_hi_se     = sd(mse_dcp_hi,     na.rm = TRUE) / sqrt(NSIMS),
       mse_kfmm_all_se   = sd(mse_kfmm_all,   na.rm = TRUE) / sqrt(NSIMS),
       mse_kfmm_hi_se    = sd(mse_kfmm_hi,    na.rm = TRUE) / sqrt(NSIMS),
       med_dcp_all_mean  = mean(med_dcp_all,  na.rm = TRUE),
       med_dcp_hi_mean   = mean(med_dcp_hi,   na.rm = TRUE),
       med_kfmm_all_mean = mean(med_kfmm_all, na.rm = TRUE),
       med_kfmm_hi_mean  = mean(med_kfmm_hi,  na.rm = TRUE))
}

# ====================================================================
# 3. Sweep (B, N) grid
# ====================================================================
cat("\n=== Step 2: phase-MSE sweep (active KIM, B x N) ===\n")
results <- list()
for (N in N_GRID) {
  for (B in B_GRID) {
    if (N %% B != 0L) {
      cat(sprintf("  Skip B=%2d N=%3d (N not divisible by B)\n", B, N))
      next
    }
    cat(sprintf("  Cell B=%2d N=%3d\n", B, N))
    t0 <- Sys.time()
    cell <- run_one_cell(B, N)
    dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (!is.null(cell)) {
      cat(sprintf("    DCP  all: mean=%.2f (SE %.2f)  med=%.2f  | hi-SNR: mean=%.2f  med=%.2f\n",
                  cell$mse_dcp_all_mean,  cell$mse_dcp_all_se, cell$med_dcp_all_mean,
                  cell$mse_dcp_hi_mean,   cell$med_dcp_hi_mean))
      cat(sprintf("    K=2  all: mean=%.2f (SE %.2f)  med=%.2f  | hi-SNR: mean=%.2f  med=%.2f  | %.1fs\n",
                  cell$mse_kfmm_all_mean, cell$mse_kfmm_all_se, cell$med_kfmm_all_mean,
                  cell$mse_kfmm_hi_mean,  cell$med_kfmm_hi_mean, dt))
      results[[sprintf("B%d_N%d", B, N)]] <- cell
    }
  }
}

# ====================================================================
# 4. Save
# ====================================================================
out <- list(
  GLOBAL_SEED = GLOBAL_SEED,
  NGENES      = NGENES,
  NSIMS       = NSIMS,
  TOP_K_FMM   = TOP_K_FMM,
  PERIOD      = PERIOD,
  pilot       = "baboon KIM (active, GSE98965)",
  B_GRID      = B_GRID,
  N_GRID      = N_GRID,
  results     = results
)
ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
rds <- file.path(out_dir_res, sprintf("fig6_phase_mse_%s.rds", ts))
saveRDS(out, rds)
cat(sprintf("\nSaved: %s\n", rds))
