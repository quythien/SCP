#' =======================================================================
#' fig6_mh_phase_mse.R — Active design phase MSE under MH truth
#' =======================================================================
#'
#' Mirror of fig6_phase_mse.R, but the data-generating model is the
#' multi-harmonic (MH) cosine series instead of FMM. Phase truth is the
#' simulator's phase_g (in hours, same as cosinor convention) — no
#' beta-shift correction needed (unlike FMM with beta = pi).
#'
#' Active KIM pilot (B in {4, 6, 8, 12, 24}) at fixed total N:
#'   - Panels A/C of Fig 6 (DCP power, K = 2 power): not regenerated
#'     here; use the existing FMM cache or regenerate offline.
#'   - Panel B (phase MSE) is the focus.
#'
#' For each (B, N) cell:
#'   (i)   simulate MH data with harmonics = c(ALPHA2_ANCHOR, 0)
#'   (ii)  estimate per-gene phi via cosinor OLS (DCP first-harmonic)
#'   (iii) also estimate via direct K = 2 OLS (same data, K = 2 model)
#'   (iv)  compare against truth phi_g (no beta shift); compute circular MSE
#'   (v)   retain per-replicate median squared errors for downstream CIs
#'
#' OUTPUT:
#'   output/mh_alternative/results/fig6_mh_phase_mse_<ts>.rds
#'
#' USAGE:
#'   SMOKE_TEST=true Rscript examples/publication/mh_alternative/fig6_mh_phase_mse.R
#'   Rscript examples/publication/mh_alternative/fig6_mh_phase_mse.R

SMOKE_TEST   <- identical(Sys.getenv("SMOKE_TEST"), "true")
NGENES       <- if (SMOKE_TEST) 500L  else 2000L
NSIMS        <- if (SMOKE_TEST) 3L    else 100L
N_CORES      <- as.integer(Sys.getenv("MC_CORES", unset = "20"))
RHYTHM_PVAL  <- 0.01
GLOBAL_SEED  <- 2025L
PERIOD       <- 24
ALPHA2_ANCHOR <- 0.3
R_FLOOR      <- 1.0

B_GRID <- if (SMOKE_TEST) c(4L, 12L) else c(4L, 6L, 8L, 12L, 24L)
N_GRID <- if (SMOKE_TEST) c(24L)     else c(24L, 48L, 72L, 96L, 120L, 144L)
# Power-sweep B grid: skip B = 4 because K = 2 is rank-deficient there.
B_GRID_PWR <- B_GRID[B_GRID >= 6]

cat(sprintf("Mode      : %s\n", if (SMOKE_TEST) "SMOKE" else "PRODUCTION"))
cat(sprintf("NGENES    : %d\n", NGENES))
cat(sprintf("NSIMS     : %d\n", NSIMS))
cat(sprintf("B_GRID    : %s\n", paste(B_GRID, collapse = ", ")))
cat(sprintf("N_GRID    : %s\n", paste(N_GRID, collapse = ", ")))
cat(sprintf("ALPHA2    : %g\n", ALPHA2_ANCHOR))

setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
old_wd <- setwd("code"); source("setup.R"); setwd(old_wd)

out_dir_res <- "output/mh_alternative/results"
dir.create(out_dir_res, recursive = TRUE, showWarnings = FALSE)

# ====================================================================
# 1. Baboon KIM pilot (cosinor estimation)
# ====================================================================
cat("\n=== Step 1: load + cosinor-fit Baboon KIM ===\n")
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

bio_kim <- estCircadianParam(mat_sub, tod_kim,
                              min_rhythm_pval = RHYTHM_PVAL,
                              verbose = TRUE)
bio_kim$ngenes <- NGENES
bio_kim$omega_rhythmic <- NULL; bio_kim$alpha_rhythmic <- NULL
bio_kim$omega_dist <- NULL;     bio_kim$alpha_dist <- NULL

# ====================================================================
# 2. Phase-MSE helper
# ====================================================================
circular_diff <- function(phi_hat, phi_truth, period = 24) {
  d <- abs(phi_hat - phi_truth) %% period
  pmin(d, period - d)
}

run_one_cell <- function(B, N) {
  if (N %% B != 0L) return(NULL)
  m <- N %/% B
  cts_active <- rep(seq(0, PERIOD * (1 - 1/B), length.out = B), each = m)

  mse_dcp_all  <- numeric(NSIMS); mse_dcp_hi  <- numeric(NSIMS)
  mse_kfmm_all <- numeric(NSIMS); mse_kfmm_hi <- numeric(NSIMS)
  med_dcp_all  <- numeric(NSIMS); med_dcp_hi  <- numeric(NSIMS)
  med_kfmm_all <- numeric(NSIMS); med_kfmm_hi <- numeric(NSIMS)

  for (s in seq_len(NSIMS)) {
    set.seed(GLOBAL_SEED + 1000L * B + 7L * N + s)
    # MH path: omega = 1, alpha2 = ALPHA2_ANCHOR. simCircadianSingleCohort
    # returns the simulator's phase_g vector (rhythmic genes only).
    sim <- simCircadianSingleCohort(bio_kim, cts_active,
                                     alpha2 = ALPHA2_ANCHOR, alpha3 = 0)
    expr   <- sim$expr
    is_rh  <- as.logical(sim$is_rhythmic)
    # sim$phase_g has length ngenes (non-rhythmic genes are 0); index by is_rh.
    phi_truth_hr <- sim$phase_g[is_rh] %% PERIOD

    r_truth <- sim$r_values[is_rh]
    hi_snr  <- !is.na(r_truth) & r_truth >= R_FLOOR

    # --- cosinor (DCP) phase estimate ---
    cos_fit <- fitCosinorAll_fast(expr, cts_active, period = PERIOD)
    phi_hat_cos <- cos_fit$phi[is_rh]
    d_cos       <- circular_diff(phi_hat_cos, phi_truth_hr, PERIOD)
    mse_dcp_all[s] <- mean(d_cos^2, na.rm = TRUE)
    med_dcp_all[s] <- median(d_cos^2, na.rm = TRUE)
    mse_dcp_hi[s]  <- mean(d_cos[hi_snr]^2, na.rm = TRUE)
    med_dcp_hi[s]  <- median(d_cos[hi_snr]^2, na.rm = TRUE)

    # --- K = 2 OLS phase estimate (correctly specified under MH truth) ---
    omega0 <- 2 * pi / PERIOD
    X_k2 <- cbind(1,
                  cos(omega0 * cts_active),    sin(omega0 * cts_active),
                  cos(2 * omega0 * cts_active), sin(2 * omega0 * cts_active))
    XtX_inv <- tryCatch(solve(crossprod(X_k2)), error = function(e) NULL)
    if (is.null(XtX_inv)) {
      mse_kfmm_all[s] <- NA_real_; med_kfmm_all[s] <- NA_real_
      mse_kfmm_hi[s]  <- NA_real_; med_kfmm_hi[s]  <- NA_real_
    } else {
      coef_k2 <- tcrossprod(expr %*% X_k2, XtX_inv)
      phi_hat_k2 <- (atan2(coef_k2[, 3], coef_k2[, 2]) / omega0) %% PERIOD
      d_k2 <- circular_diff(phi_hat_k2[is_rh], phi_truth_hr, PERIOD)
      mse_kfmm_all[s] <- mean(d_k2^2, na.rm = TRUE)
      med_kfmm_all[s] <- median(d_k2^2, na.rm = TRUE)
      mse_kfmm_hi[s]  <- mean(d_k2[hi_snr]^2, na.rm = TRUE)
      med_kfmm_hi[s]  <- median(d_k2[hi_snr]^2, na.rm = TRUE)
    }
  }

  list(B = B, N = N, R_FLOOR = R_FLOOR,
       mse_dcp_all_raw   = mse_dcp_all,
       mse_dcp_hi_raw    = mse_dcp_hi,
       mse_kfmm_all_raw  = mse_kfmm_all,
       mse_kfmm_hi_raw   = mse_kfmm_hi,
       med_dcp_all_raw   = med_dcp_all,
       med_dcp_hi_raw    = med_dcp_hi,
       med_kfmm_all_raw  = med_kfmm_all,
       med_kfmm_hi_raw   = med_kfmm_hi,
       mse_dcp_all_mean  = mean(mse_dcp_all,  na.rm = TRUE),
       mse_dcp_hi_mean   = mean(mse_dcp_hi,   na.rm = TRUE),
       mse_kfmm_all_mean = mean(mse_kfmm_all, na.rm = TRUE),
       mse_kfmm_hi_mean  = mean(mse_kfmm_hi,  na.rm = TRUE),
       med_dcp_all_mean  = mean(med_dcp_all,  na.rm = TRUE),
       med_dcp_hi_mean   = mean(med_dcp_hi,   na.rm = TRUE),
       med_kfmm_all_mean = mean(med_kfmm_all, na.rm = TRUE),
       med_kfmm_hi_mean  = mean(med_kfmm_hi,  na.rm = TRUE))
}

# ====================================================================
# 3. Sweep
# ====================================================================
cat("\n=== Step 2: phase-MSE sweep (active KIM, MH truth) ===\n")
results <- list()
for (N in N_GRID) {
  for (B in B_GRID) {
    if (N %% B != 0L) next
    cat(sprintf("  Cell B=%2d N=%3d\n", B, N))
    t0 <- Sys.time()
    cell <- run_one_cell(B, N)
    dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (!is.null(cell)) {
      cat(sprintf("    DCP  all: mean=%.2f  med=%.3f  | K=2 med=%.3f  | %.1fs\n",
                  cell$mse_dcp_all_mean,  cell$med_dcp_all_mean,
                  cell$med_kfmm_all_mean, dt))
      results[[sprintf("B%d_N%d", B, N)]] <- cell
    }
  }
}

# ====================================================================
# 4. Additional pass: detection power per (B, N) for DCP and K=2 under MH truth.
#    Powers populate Fig 6 MH Panels A and C; phase MSE populates Panel B.
# ====================================================================
cat("\n=== Step 3: detection power per (B, N), DCP + K=2 (MH truth) ===\n")
analysis_pwr <- CircadianAnalysisOptions(alpha           = 0.05,
                                          p.adjust.method = "BH",
                                          fdr_thresholds  = 0.05)
NSIMS_PWR <- if (SMOKE_TEST) 3L else 20L
B_GRID_PWR <- intersect(B_GRID, c(6L, 8L, 12L, 24L))  # power omits B=4 (K=2 not identifiable)

power_dcp <- matrix(NA_real_, nrow = length(B_GRID_PWR), ncol = length(N_GRID),
                    dimnames = list(paste0("B=", B_GRID_PWR),
                                    paste0("N=", N_GRID)))
power_k2 <- power_dcp

for (i in seq_along(B_GRID_PWR)) {
  for (j in seq_along(N_GRID)) {
    B <- B_GRID_PWR[i]; N <- N_GRID[j]
    if (N %% B != 0L) next
    m <- N %/% B
    cts_design <- rep(seq(0, PERIOD * (1 - 1/B), length.out = B), each = m)
    design <- CircadianDesignOptions(sample_sizes = N, nsims = NSIMS_PWR,
                                      design = "active",
                                      cts = cts_design, B_values = B)
    set.seed(GLOBAL_SEED)
    res_dcp <- runSimsSingleCohort(bio_kim, design, analysis_pwr,
                                    method = "DCP",
                                    harmonics = c(ALPHA2_ANCHOR, 0),
                                    mc.cores = N_CORES, verbose = FALSE)
    set.seed(GLOBAL_SEED)
    res_k2 <- runSimsSingleCohort(bio_kim, design, analysis_pwr,
                                   method = "FMM", K = 2L,
                                   harmonics = c(ALPHA2_ANCHOR, 0),
                                   mc.cores = N_CORES, verbose = FALSE)
    power_dcp[i, j] <- mean(res_dcp$marginal_power, na.rm = TRUE)
    power_k2 [i, j] <- mean(res_k2 $marginal_power, na.rm = TRUE)
  }
  cat(sprintf("  B=%2d: DCP %s | K2 %s\n", B_GRID_PWR[i],
              paste(sprintf("%.2f", power_dcp[i, ]), collapse = " "),
              paste(sprintf("%.2f", power_k2 [i, ]), collapse = " ")))
}

# ====================================================================
# 5. Save
# ====================================================================
out <- list(
  framework      = "MH alternative -- active design phase MSE + power",
  GLOBAL_SEED    = GLOBAL_SEED,
  NGENES         = NGENES, NSIMS = NSIMS, NSIMS_PWR = NSIMS_PWR,
  PERIOD         = PERIOD,
  pilot          = "baboon KIM (active, GSE98965)",
  alpha2_anchor  = ALPHA2_ANCHOR,
  B_GRID         = B_GRID,
  B_GRID_PWR     = B_GRID_PWR,
  N_GRID         = N_GRID,
  results        = results,       # Panel B (phase MSE)
  power_dcp      = power_dcp,     # Panel A
  power_k2       = power_k2       # Panel C
)
ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
rds <- file.path(out_dir_res, sprintf("fig6_mh_phase_mse_%s.rds", ts))
saveRDS(out, rds)
cat(sprintf("\nSaved: %s\n", rds))
