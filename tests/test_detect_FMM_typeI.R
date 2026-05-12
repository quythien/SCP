#' =======================================================================
#' tests/test_detect_FMM_typeI.R
#'
#' Type-I error regression test for detect_FMM (K-harmonic LRT)
#' and detect_DCP (cosinor F-test).
#'
#' Verifies under-null calibration at multiple sample sizes and
#' alpha levels. Empirical rejection rates must match nominal alpha
#' within Monte Carlo tolerance. Also checks p-value uniformity
#' via Kolmogorov-Smirnov. detect_FMM is tested on both the
#' exact-F path (ebayes=FALSE) and the moderated-F path
#' (ebayes=TRUE, the LimoRhyde default).
#'
#' Run with:
#'   Rscript tests/test_detect_FMM_typeI.R
#' =======================================================================

setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
old_wd <- setwd("code"); source("setup.R"); setwd(old_wd)

# Test grid
N_VALUES <- c(12L, 24L, 48L, 72L, 144L)
N_SIMS   <- 2000L         # Monte Carlo SE on T1 at 0.05: ~ 0.005
PERIOD   <- 24
SEED     <- 2026L

# Tolerance bounds: nominal +/- 3 SE-ish at B = 2000
T1_TOL_05 <- c(0.035, 0.065)   # ~ 0.05 +/- 0.015
T1_TOL_01 <- c(0.004, 0.020)   # ~ 0.01 +/- 0.010 (looser for tail)

# ---------------------------------------------------------------------
# Helper: simulate B null gene matrices and run a detector
# ---------------------------------------------------------------------
.run_null <- function(n, n_sim, detector_fn, seed) {
  set.seed(seed)
  times <- seq(0, PERIOD * (1 - 1/n), length.out = n)   # uniform active grid
  null_mat <- matrix(rnorm(n_sim * n), nrow = n_sim, ncol = n)
  pvals <- detector_fn(null_mat, times)
  list(pvals = pvals, times = times)
}

.summary <- function(name, n, pv) {
  pv_clean <- pv[!is.na(pv)]
  if (length(pv_clean) < 50L) {
    cat(sprintf("  %-12s n=%-3d  FAIL (only %d non-NA p-values)\n",
                name, n, length(pv_clean)))
    return(list(pass = FALSE, t1_05 = NA, t1_01 = NA, ks_p = NA))
  }
  t1_05 <- mean(pv_clean < 0.05)
  t1_01 <- mean(pv_clean < 0.01)
  ks_p  <- suppressWarnings(
    ks.test(pv_clean + runif(length(pv_clean)) * 1e-9, "punif")$p.value)
  pass_05 <- t1_05 >= T1_TOL_05[1] && t1_05 <= T1_TOL_05[2]
  pass_01 <- t1_01 >= T1_TOL_01[1] && t1_01 <= T1_TOL_01[2]
  pass    <- pass_05 && pass_01
  cat(sprintf("  %-12s n=%-3d  T1@.05=%.3f  T1@.01=%.3f  KS_p=%.3g  %s\n",
              name, n, t1_05, t1_01, ks_p,
              if (pass) "PASS" else "FAIL"))
  list(pass = pass, t1_05 = t1_05, t1_01 = t1_01, ks_p = ks_p)
}

# ---------------------------------------------------------------------
# 1. detect_DCP type-I error
# ---------------------------------------------------------------------
cat("\n=== detect_DCP type-I error (cosinor F-test) ===\n")
fail_DCP <- 0L
for (n in N_VALUES) {
  r <- .run_null(n, N_SIMS,
                 function(M, t) detect_DCP(M, t, period = PERIOD),
                 seed = SEED + n)
  s <- .summary("DCP", n, r$pvals)
  if (!s$pass) fail_DCP <- fail_DCP + 1L
}

# ---------------------------------------------------------------------
# 2. detect_FMM type-I error at K = 1, both ebayes paths
#    - ebayes=FALSE verifies exact-F calibration (F(2K, n-2K-1) under H0)
#    - ebayes=TRUE  verifies the moderated-F calibration that LimoRhyde uses
# ---------------------------------------------------------------------
cat("\n=== detect_FMM K=1 (equivalent to DCP) ===\n")
fail_K1 <- 0L
for (n in N_VALUES) {
  r_eF <- .run_null(n, N_SIMS,
                 function(M, t) detect_FMM(M, t, period = PERIOD,
                                            K = 1L, ebayes = FALSE,
                                            adjust.method = NULL)$p.value,
                 seed = SEED + 100L + n)
  s_eF <- .summary("FMM K=1 exactF", n, r_eF$pvals)
  if (!s_eF$pass) fail_K1 <- fail_K1 + 1L

  r_eT <- .run_null(n, N_SIMS,
                 function(M, t) detect_FMM(M, t, period = PERIOD,
                                            K = 1L, ebayes = TRUE,
                                            adjust.method = NULL)$p.value,
                 seed = SEED + 150L + n)
  s_eT <- .summary("FMM K=1 eBayes", n, r_eT$pvals)
  if (!s_eT$pass) fail_K1 <- fail_K1 + 1L
}

# ---------------------------------------------------------------------
# 3. detect_FMM type-I error at K = 2 (default; PAPER PRIMARY), both ebayes paths
# ---------------------------------------------------------------------
cat("\n=== detect_FMM K=2 (paper default) ===\n")
fail_K2 <- 0L
for (n in N_VALUES) {
  r_eF <- .run_null(n, N_SIMS,
                 function(M, t) detect_FMM(M, t, period = PERIOD,
                                            K = 2L, ebayes = FALSE,
                                            adjust.method = NULL)$p.value,
                 seed = SEED + 200L + n)
  s_eF <- .summary("FMM K=2 exactF", n, r_eF$pvals)
  if (!s_eF$pass) fail_K2 <- fail_K2 + 1L

  r_eT <- .run_null(n, N_SIMS,
                 function(M, t) detect_FMM(M, t, period = PERIOD,
                                            K = 2L, ebayes = TRUE,
                                            adjust.method = NULL)$p.value,
                 seed = SEED + 250L + n)
  s_eT <- .summary("FMM K=2 eBayes", n, r_eT$pvals)
  if (!s_eT$pass) fail_K2 <- fail_K2 + 1L
}

# ---------------------------------------------------------------------
# 4. detect_FMM type-I error at K = 3, exact-F only
#    (eBayes path uses the same moderated-F across K, already exercised at K=2)
# ---------------------------------------------------------------------
cat("\n=== detect_FMM K=3 exact-F (smaller n excluded for identifiability) ===\n")
fail_K3 <- 0L
for (n in N_VALUES) {
  if (n < 8L) next
  r <- .run_null(n, N_SIMS,
                 function(M, t) detect_FMM(M, t, period = PERIOD,
                                            K = 3L, ebayes = FALSE,
                                            adjust.method = NULL)$p.value,
                 seed = SEED + 300L + n)
  s <- .summary("FMM K=3 exactF", n, r$pvals)
  if (!s$pass) fail_K3 <- fail_K3 + 1L
}

# ---------------------------------------------------------------------
# 5. Rank-deficient design (B < 2K+1): K=2 at B=4 must warn AND be conservative
# ---------------------------------------------------------------------
cat("\n=== detect_FMM K=2 at B=4 (rank-deficient, must be conservative) ===\n")
B_VAL <- 4L
N_RD  <- 24L
m_rd  <- N_RD %/% B_VAL
times_rd <- rep(seq(0, PERIOD * (1 - 1/B_VAL), length.out = B_VAL), each = m_rd)
set.seed(SEED + 400L)
null_mat_rd <- matrix(rnorm(N_SIMS * N_RD), nrow = N_SIMS, ncol = N_RD)
warning_seen <- FALSE
withCallingHandlers(
  {
    res_rd <- detect_FMM(null_mat_rd, times_rd, period = PERIOD,
                         K = 2L, adjust.method = NULL)
  },
  warning = function(w) {
    if (grepl("Nyquist|distinct timepoints", conditionMessage(w))) {
      warning_seen <<- TRUE
    }
    invokeRestart("muffleWarning")
  }
)
t1_rd <- mean(res_rd$p.value < 0.05, na.rm = TRUE)
cat(sprintf("  rank-defc    n=%d (B=%d)  T1@.05=%.3f (must be <= 0.05) ",
            N_RD, B_VAL, t1_rd))
fail_RD <- 0L
if (!warning_seen) {
  cat("FAIL (warning did not fire)\n"); fail_RD <- 1L
} else if (t1_rd > 0.055) {
  cat("FAIL (anti-conservative T1 in rank-deficient regime)\n"); fail_RD <- 1L
} else {
  cat("PASS\n")
}

# ---------------------------------------------------------------------
# 6. Final report
# ---------------------------------------------------------------------
total_fails <- fail_DCP + fail_K1 + fail_K2 + fail_K3 + fail_RD
cat(sprintf("\n=== Summary: %d failures (DCP=%d, K=1=%d, K=2=%d, K=3=%d, RD=%d) ===\n",
            total_fails, fail_DCP, fail_K1, fail_K2, fail_K3, fail_RD))
if (total_fails == 0L) {
  cat("ALL TESTS PASS\n")
  quit(status = 0)
} else {
  cat("SOME TESTS FAILED. See details above.\n")
  quit(status = 1)
}
