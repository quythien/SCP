# cosinor_fast.R
#
# R-side wrappers for the Rcpp implementations in cosinor_fast.cpp.
#
# Usage: source this file AFTER calling Rcpp::sourceCpp() on cosinor_fast.cpp.
# Both functions are drop-in replacements for the R originals; they accept the
# same arguments and return the same types.
#
# Loading:
#   Rcpp::sourceCpp("code/src/cosinor_fast.cpp")
#   source("code/src/cosinor_fast.R")
#
# Then:
#   fitCosinorAll_fast()     replaces fitCosinorAll() in bootstrap_sim.R
#   sim_cosinor_expr_fast()  replaces the for-g loop in simCircadianDiff() and
#                            runSimsSingleCohort() in simulation.R / runner.R


# ---------------------------------------------------------------------------
# fitCosinorAll_fast()
#
# Drop-in replacement for fitCosinorAll() in bootstrap_sim.R.
# Returns a data.frame with the same columns (gene, M, A, phi, sigma, pvalue,
# r, is_rhythmic) but computes the cosinor F-test p-values via C++.
#
# The p-values are numerically identical to one_cosinor_OLS() up to floating-
# point rounding (max difference < 2e-12 in practice; see test results).
#
# Arguments:
#   data           : G x N matrix (genes x samples)
#   times          : numeric vector length N
#   period         : circadian period (default 24)
#   min_rhythm_pval: threshold for is_rhythmic flag (default 0.01)
# ---------------------------------------------------------------------------
fitCosinorAll_fast <- function(data, times, period = 24, min_rhythm_pval = 0.01) {

  if (!exists("cosinor_pvals_cpp", mode = "function"))
    stop("cosinor_pvals_cpp not found. Call Rcpp::sourceCpp('code/src/cosinor_fast.cpp') first.")

  G     <- nrow(data)
  omega <- 2 * pi / period

  # Cosinor p-values via C++ (replaces the lapply/one_cosinor_OLS loop)
  pvals <- cosinor_pvals_cpp(data, times, period)

  # The C++ kernel returns p-values only; we still need M, A, phi, sigma for
  # the data.frame.  These are needed only for pilot fitting (fitCosinorAll),
  # which runs once, not inside the simulation loop.  Use a vectorised OLS
  # pass to recover estimates without re-entering C++ (the fitted values are
  # a cheap matrix operation at this scale).
  #
  # Design matrix (shared across all genes)
  x1     <- cos(omega * times)
  x2     <- sin(omega * times)
  X      <- cbind(1, x1, x2)                 # N x 3
  XtX_inv <- tryCatch(solve(crossprod(X)), error = function(e) NULL)

  if (is.null(XtX_inv)) {
    # Singular design: return NA estimates with p-values from C++
    df <- data.frame(
      gene       = seq_len(G),
      M          = NA_real_,
      A          = NA_real_,
      phi        = NA_real_,
      sigma      = NA_real_,
      pvalue     = pvals,
      r          = NA_real_,
      is_rhythmic = !is.na(pvals) & pvals < min_rhythm_pval,
      stringsAsFactors = FALSE
    )
    return(df)
  }

  # Vectorised OLS: (G x N) %*% (N x 3) -> G x 3 coefficient matrix
  Xty  <- data %*% X                          # G x 3
  coef <- tcrossprod(Xty, XtX_inv)           # G x 3  (= Xty %*% t(XtX_inv))

  M_hat    <- coef[, 1]
  beta1    <- coef[, 2]
  beta2    <- coef[, 3]
  A_hat    <- sqrt(beta1^2 + beta2^2)

  # Phase: atan2(beta2, beta1) / omega, wrapped to [0, period)
  phi_hat  <- (atan2(beta2, beta1) / omega) %% period

  # Residual sigma: sqrt( RSS / (N - 3) )
  yhat     <- tcrossprod(coef, X)             # G x N
  resid    <- data - yhat                     # G x N
  RSS_g    <- rowSums(resid^2)
  df_resid  <- length(times) - 3
  sigma_hat <- if (df_resid > 0) sqrt(RSS_g / df_resid) else rep(NA_real_, length(RSS_g))

  r_hat <- A_hat / pmax(sigma_hat, 1e-6)

  data.frame(
    gene        = seq_len(G),
    M           = M_hat,
    A           = A_hat,
    phi         = phi_hat,
    sigma       = sigma_hat,
    pvalue      = pvals,
    r           = r_hat,
    is_rhythmic = !is.na(pvals) & pvals < min_rhythm_pval & !is.na(A_hat) & A_hat > 0,
    stringsAsFactors = FALSE
  )
}


# ---------------------------------------------------------------------------
# sim_cosinor_expr_fast()
#
# Drop-in replacement for the per-gene for-loop in:
#   - simCircadianDiff()          [simulation.R, lines ~360-380]
#   - runSimsSingleCohort()       [runner.R, lines ~917-924]
#   - simCircadianSingleCohort()  [simulation.R, lines ~506-513]
#
# Arguments match sim_cosinor_expr_cpp() in cosinor_fast.cpp.
# Returns a G x N numeric matrix, identical (up to RNG state) to the R loop.
#
# Note on reproducibility: the C++ function calls R's RNG (R::rnorm) for each
# draw in the same row-major order as the R for-loop, so set.seed() in R
# produces identical results for both implementations.
# ---------------------------------------------------------------------------
sim_cosinor_expr_fast <- function(mesor, amplitude, phase, sigma, times,
                                   period = 24, alpha2 = 0, alpha3 = 0) {
  if (!exists("sim_cosinor_expr_cpp", mode = "function"))
    stop("sim_cosinor_expr_cpp not found. Call Rcpp::sourceCpp('code/src/cosinor_fast.cpp') first.")

  sim_cosinor_expr_cpp(mesor, amplitude, phase, sigma, times,
                        period = period, alpha2 = alpha2, alpha3 = alpha3)
}
