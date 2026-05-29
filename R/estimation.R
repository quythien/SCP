#' Estimate Circadian Parameters from Pilot Data
#'
#' @description
#' Fits a cosinor model to each gene in the pilot dataset, selects the
#' top-K rhythmic genes (ranked by p-value, K = min(300, n_rhythmic)),
#' and returns empirical distributions of mesor, amplitude, phase, and
#' noise as a \code{CircadianBioOptions} object for use in downstream
#' simulation and power analysis.
#'
#' @param data Gene expression matrix (genes x samples).
#' @param times Numeric vector of sample time points (hours), length = \code{ncol(data)}.
#' @param period Circadian period in hours (default 24).
#' @param min_rhythm_pval P-value threshold for classifying a gene as rhythmic (default 0.01).
#' @param verbose Print estimation summary (default TRUE).
#'
#' @return A \code{CircadianBioOptions} S3 object (27 fields) including:
#'   \code{prop_rhythmic}, \code{lBaselineExpr}, \code{lOD}, \code{amplitude},
#'   \code{sigma_rhythmic}, \code{cts}, \code{phase}, \code{prop_DR},
#'   \code{prop_DP}, \code{prop_DM}, \code{phase_diff}, \code{amp_diff},
#'   \code{mesor_diff}, \code{ngenes}, \code{period}, \code{sim.seed}.
#'   Pass directly to \code{runSimsSingleCohort()} or \code{runDifferentialPower()}.
#'
#' @seealso \code{\link{estCircadianParamTwoGroup}} for two-group differential setup,
#'   \code{\link{CircadianBioOptions}} for the options constructor.


#' Wrap Angle to [0, 2pi)
#'
#' @param x Numeric. Angle(s) in radians.
#' @return Numeric. Angle(s) mapped to \code{[0, 2pi)}.
adjust.to.2pi = function(x) x %% (2 * pi)

#' Fit a Cosinor Model to a Single Gene by OLS
#'
#' @description
#' Fits the cosinor model \eqn{y = M + A\cos(\omega t - \phi) + \epsilon}
#' via ordinary least squares by reparameterising as
#' \eqn{y = M + \beta_1\cos(\omega t) + \beta_2\sin(\omega t) + \epsilon}.
#' Returns mesor, amplitude, acrophase, p-value, and R^2.
#'
#' @param tod Numeric vector. Sample time points (hours).
#' @param y Numeric vector. Gene expression values (same length as \code{tod}).
#' @param period Numeric. Rhythm period in hours (default 24).
#' @param compute.phase.CI Logical. Reserved; not yet implemented (default FALSE).
#' @param CI.level Numeric. Confidence level for phase CI (default 0.95).
#'
#' @return Named list with fields: \code{M} (mesor), \code{A} (amplitude),
#'   \code{phi} (acrophase in hours), \code{pvalue} (F-test p-value),
#'   \code{R2} (coefficient of determination). Returns \code{NA} for all
#'   fields if the design matrix is singular or \eqn{n \le 3}.
one_cosinor_OLS = function(tod, y, period = 24, compute.phase.CI = FALSE, CI.level = 0.95) {
  n     = length(tod)
  omega = 2 * pi / period
  x1    = cos(omega * tod)
  x2    = sin(omega * tod)

  mat.S = matrix(c(n, sum(x1), sum(x2),
                   sum(x1), sum(x1^2), sum(x1 * x2),
                   sum(x2), sum(x1 * x2), sum(x2^2)),
                 nrow = 3, byrow = TRUE)
  vec.d = c(sum(y), sum(y * x1), sum(y * x2))

  mat.S.inv = tryCatch(solve(mat.S), error = function(e) NULL)
  if (is.null(mat.S.inv)) return(list(pvalue = NA, M = NA, A = NA, phi = NA))

  est        = mat.S.inv %*% vec.d
  m.hat      = est[1]
  beta1.hat  = est[2]
  beta2.hat  = est[3]
  A.hat      = sqrt(beta1.hat^2 + beta2.hat^2)
  phase.hat  = adjust.to.2pi(atan2(beta2.hat, beta1.hat)) / omega

  TSS  = sum((y - mean(y))^2)
  yhat = m.hat + beta1.hat * x1 + beta2.hat * x2
  RSS  = sum((y - yhat)^2)
  MSS  = TSS - RSS

  if (n <= 3 || RSS < 0) return(list(pvalue = NA, M = m.hat, A = A.hat, phi = phase.hat))

  Fstat = (MSS / 2) / (RSS / (n - 3))
  pval  = stats::pf(Fstat, 2, n - 3, lower.tail = FALSE)

  list(M = m.hat, A = A.hat, phi = phase.hat, pvalue = pval, R2 = MSS / TSS)
}

estimate_circadian_params = function(data, times, period = 24,
                                     min_rhythm_pval = 0.01,
                                     verbose = TRUE) {

  if (verbose) {
    cat("=== Estimating Circadian Parameters from Pilot Data ===\n")
    cat("Genes:", nrow(data), "\n")
    cat("Samples:", ncol(data), "\n")
    cat("Time points:", length(unique(times)), "\n\n")
  }

  G = nrow(data)

  # Fit cosinor to each gene
  fits = lapply(1:G, function(g) {
    tryCatch({
      one_cosinor_OLS(times, data[g, ], period, compute.phase.CI = FALSE)
    }, error = function(e) {
      list(M = NA, A = NA, phi = NA, sigma = NA, pvalue = NA, r = NA)
    })
  })

  # Extract parameters
  M_vals = sapply(fits, function(x) x$M)
  A_vals = sapply(fits, function(x) x$A)
  phi_vals = sapply(fits, function(x) x$phi)
  pvals = sapply(fits, function(x) x$pvalue)

  # Estimate sigma from residuals
  sigma_vals = sapply(1:G, function(g) {
    y = data[g, ]
    fit = fits[[g]]
    if (any(is.na(c(fit$M, fit$A, fit$phi)))) return(NA)

    omega = 2 * pi / period
    yhat = fit$M + fit$A * cos(omega * times - omega * fit$phi)
    sqrt(sum((y - yhat)^2, na.rm = TRUE) / (length(y) - 3))
  })

  # Effect size
  r_vals = A_vals / sigma_vals

  # G_R^cand: all genes passing alpha_pilot (used for prop_rhythmic and is_rhythmic)
  rhythmic_genes = pvals < min_rhythm_pval & !is.na(A_vals) & A_vals > 0
  n_cand <- sum(rhythmic_genes, na.rm = TRUE)

  # G_R: top min(300, |G_R^cand|) by p-value (used for F_{A,sigma} and F_phi only)
  # Ranking by p-value gives the strongest-signal genes, improving distribution stability
  K <- min(300L, n_cand)
  if (K > 0) {
    cand_idx     <- which(rhythmic_genes)
    top_idx      <- cand_idx[order(pvals[cand_idx])][seq_len(K)]
    estim_set    <- rep(FALSE, G)
    estim_set[top_idx] <- TRUE
  } else {
    estim_set <- rhythmic_genes
  }

  if (verbose) {
    cat("Genes passing rhythm filter (p <", min_rhythm_pval, "):",
        n_cand, " |  estimation set K =", K, "\n")
  }

  # Estimate parameter distributions
  params = list(
    # Mesor: all G_0 genes (reflects shared baseline regardless of rhythmic status)
    M_mean = mean(M_vals, na.rm = TRUE),
    M_sd = sd(M_vals, na.rm = TRUE),

    # Amplitude: top-K rhythmic genes only
    A_mean = mean(A_vals[estim_set], na.rm = TRUE),
    A_sd = sd(A_vals[estim_set], na.rm = TRUE),
    A_median = median(A_vals[estim_set], na.rm = TRUE),
    A_q25 = quantile(A_vals[estim_set], 0.25, na.rm = TRUE),
    A_q75 = quantile(A_vals[estim_set], 0.75, na.rm = TRUE),

    # Phase: top-K rhythmic genes only
    phi_mean = circular_mean(phi_vals[estim_set], period),
    phi_concentration = circular_concentration(phi_vals[estim_set], period),
    phi_uniform = FALSE,

    # Noise level (all genes)
    sigma_mean = mean(sigma_vals, na.rm = TRUE),
    sigma_sd = sd(sigma_vals, na.rm = TRUE),
    sigma_median = median(sigma_vals, na.rm = TRUE),

    # Effect size: top-K rhythmic genes only
    r_mean = mean(r_vals[estim_set], na.rm = TRUE),
    r_sd = sd(r_vals[estim_set], na.rm = TRUE),
    r_median = median(r_vals[estim_set], na.rm = TRUE),
    r_q25 = quantile(r_vals[estim_set], 0.25, na.rm = TRUE),
    r_q75 = quantile(r_vals[estim_set], 0.75, na.rm = TRUE),

    # Proportion rhythmic: full candidate set / G_0 (not capped)
    prop_rhythmic = n_cand / G
  )

  # Test for uniform phase distribution using top-K genes
  if (length(na.omit(phi_vals[estim_set])) > 10) {
    rayleigh_test = rayleigh_test_circular(phi_vals[estim_set], period)
    params$phi_uniform = rayleigh_test$pvalue > 0.05
    params$phi_rayleigh_p = rayleigh_test$pvalue
  }

  # Store raw values for flexible simulation. `gene` carries the source gene
  # identifiers (matrix rownames) so the per-gene rhythm table can report which
  # genes are rhythmic; falls back to integer positions if the matrix is unnamed.
  params$raw = list(
    M = M_vals,
    A = A_vals,
    phi = phi_vals,
    sigma = sigma_vals,
    r = r_vals,
    pvalue = pvals,
    gene = if (!is.null(rownames(data))) rownames(data) else as.character(seq_len(nrow(data))),
    is_rhythmic = rhythmic_genes,  # G_R^cand: full set for prop_rhythmic / DR/DP/DM typing
    in_estim_set = estim_set       # G_R: top-K for F_{A,sigma}, F_phi
  )

  if (verbose) {
    cat("\n=== Parameter Estimates ===\n")
    cat(sprintf("Mesor: %.2f ± %.2f\n", params$M_mean, params$M_sd))
    cat(sprintf("Amplitude (rhythmic): %.3f ± %.3f (median: %.3f)\n",
                params$A_mean, params$A_sd, params$A_median))
    cat(sprintf("Phase mean: %.1f h (concentration: %.2f)\n",
                params$phi_mean, params$phi_concentration))
    cat(sprintf("Noise (σ): %.3f ± %.3f\n", params$sigma_mean, params$sigma_sd))
    cat(sprintf("Effect size (r = A/σ): %.3f ± %.3f (median: %.3f)\n",
                params$r_mean, params$r_sd, params$r_median))
    cat(sprintf("Proportion rhythmic: %.1f%%\n", 100 * params$prop_rhythmic))
    cat(sprintf("Phase distribution: %s\n",
                ifelse(params$phi_uniform, "Uniform", "Concentrated")))
  }

  return(params)
}


#' Circular mean calculation
#'
#' @param angles Angles (in hours)
#' @param period Period (default 24)
#'
#' @return Mean angle (in hours)
circular_mean = function(angles, period = 24) {
  angles = angles[!is.na(angles)]
  if (length(angles) == 0) return(NA)

  # Convert to radians
  omega = 2 * pi / period
  radians = omega * angles

  # Circular mean
  mean_angle = atan2(mean(sin(radians)), mean(cos(radians)))

  # Convert back to hours
  (mean_angle / omega) %% period
}


#' Circular concentration (1 - circular variance)
#'
#' @param angles Angles (in hours)
#' @param period Period (default 24)
#'
#' @return Concentration (0 = uniform, 1 = all same)
circular_concentration = function(angles, period = 24) {
  angles = angles[!is.na(angles)]
  if (length(angles) == 0) return(NA)

  # Convert to radians
  omega = 2 * pi / period
  radians = omega * angles

  # Resultant vector length
  R = sqrt(mean(sin(radians))^2 + mean(cos(radians))^2)

  return(R)
}


#' Rayleigh test for uniform circular distribution
#'
#' @param angles Angles (in hours)
#' @param period Period (default 24)
#'
#' @return List with statistic and p-value
rayleigh_test_circular = function(angles, period = 24) {
  angles = angles[!is.na(angles)]
  n = length(angles)
  if (n < 5) return(list(statistic = NA, pvalue = NA))

  # Convert to radians
  omega = 2 * pi / period
  radians = omega * angles

  # Resultant vector length
  R = sqrt(sum(sin(radians))^2 + sum(cos(radians))^2)

  # Rayleigh statistic
  Z = R^2 / n

  # P-value: large-sample exponential approximation. For n < 25 this can be
  # noticeably anti-conservative; warn so callers can apply Mardia's correction
  # or the exact small-n table if needed.
  pvalue = exp(-Z)
  if (n < 25L)
    warning(sprintf(
      "rayleigh_test_circular: n=%d is below 25; p-value approximation unreliable. ",
      n),
      "Consider circular::rayleigh.test() or Mardia's correction for small n.")

  return(list(
    statistic = Z,
    pvalue = pvalue
  ))
}


#' Estimate Circadian Parameters from Pilot Data (Single-Cohort)
#'
#' @description
#' Fits a cosinor model per gene, ranks rhythmic genes by p-value, and returns
#' empirical parameter distributions packaged as a \code{CircadianBioOptions}
#' object suitable for downstream simulation and power analysis.
#'
#' @param data Numeric matrix of pilot expression (genes x samples).
#'   A data.frame is coerced to matrix.
#' @param times Numeric vector of sample collection times in hours,
#'   length \code{ncol(data)}. The function does NOT validate units;
#'   pass hours (0 to 24 scale).
#' @param period Period in the same units as \code{times} (default 24 hours).
#' @param min_rhythm_pval Uncorrected cosinor p-value threshold for the
#'   rhythmic-candidate set (default 0.01). Used to compute prop_rhythmic
#'   and to select the top-K (K = min(300, n_rhythmic)) genes whose
#'   amplitudes, phases, and noise drive simulation truth.
#' @param prop_DR,prop_DP,prop_DM Proportions of differentially-rhythmic,
#'   differentially-phased, and differentially-mesor genes (defaults
#'   0.15, 0.10, 0.00). Capped automatically at prop_rhythmic.
#' @param mesor_diff,phase_diff,amp_diff Numeric length-2 ranges for
#'   the differential effect-size sampling.
#' @param dp_shift_mode One of \code{"fixed"} or \code{"uniform"} (default
#'   "fixed").
#' @param dr_amp_scale,dr_sigma_scale Scaling factors for DR-class amplitude
#'   and noise distributions (default 1).
#' @param paired_sigma Logical, default FALSE. If TRUE, joint (A, sigma)
#'   sampling preserves pilot correlation in downstream simulations.
#' @param sim.seed Integer seed for downstream simulation (default 12345).
#' @param verbose Logical, default TRUE. Print summary at end.
#'
#' @return A \code{CircadianBioOptions} S3 object usable directly with
#'   \code{runSimsSingleCohort()} or \code{runDifferentialPower()}.
#'
#' @seealso \code{\link{estCircadianParamFMM}} for FMM-aware pre-screen,
#'   \code{\link{estCircadianParamTwoGroup}} for two-group setup,
#'   \code{\link{prepCircadianData}} for normalisation upstream.
#'
#' @examples
#' \dontrun{
#'   set.seed(1)
#'   n_genes <- 200; n_samples <- 30
#'   times <- rep(seq(0, 22, by = 2), length.out = n_samples)
#'   expr  <- matrix(rnorm(n_genes * n_samples), nrow = n_genes)
#'   for (g in 1:20) expr[g, ] <- expr[g, ] + 2 * cos(2*pi*times/24)
#'   bio <- estCircadianParam(expr, times, period = 24, verbose = FALSE)
#' }
#'
#' @export
estCircadianParam <- function(data, times, period = 24,
                              min_rhythm_pval = 0.01,
                              prop_DR = 0.15, prop_DP = 0.10,
                              prop_DM = 0.00, mesor_diff = c(0.5, 2.0),
                              phase_diff = c(-6, 6), amp_diff = c(0.5, 2),
                              dp_shift_mode = c("fixed", "uniform"),
                              dr_amp_scale = 1.0,
                              dr_sigma_scale = 1.0,
                              paired_sigma = FALSE,
                              sim.seed = 12345, verbose = TRUE) {

  # Input validation: catch transposed-matrix and wrong-type errors at the
  # earliest opportunity rather than producing silently wrong estimates.
  if (is.data.frame(data)) data <- as.matrix(data)
  if (!is.matrix(data) || !is.numeric(data))
    stop("'data' must be a numeric matrix (genes x samples). Got: ",
         paste(class(data), collapse = "/"))
  if (length(times) != ncol(data))
    stop(sprintf(
      "length(times) (%d) must equal ncol(data) (%d). Did you mean to transpose your matrix?",
      length(times), ncol(data)))

  dp_shift_mode <- match.arg(dp_shift_mode)

  params <- estimate_circadian_params(data, times, period = period,
                                      min_rhythm_pval = min_rhythm_pval,
                                      verbose = verbose)

  ngenes <- nrow(data)
  rhythmic_idx <- params$raw$is_rhythmic

  # Build empirical vectors for CircadianBioOptions
  # lBaselineExpr: use all gene mesors (log scale)
  lBaselineExpr_emp <- params$raw$M[!is.na(params$raw$M)]

  # lOD (sigma): use all gene sigmas, log-transformed
  sigma_valid <- params$raw$sigma[!is.na(params$raw$sigma) & params$raw$sigma > 0]
  lOD_emp <- log(sigma_valid)

  # Amplitude and sigma: top-K estimation set (G_R), paired to preserve A-sigma correlation
  # Phase: top-K estimation set only
  estim_idx   <- params$raw$in_estim_set
  estim_valid <- estim_idx & !is.na(params$raw$A) & params$raw$A > 0 &
                 !is.na(params$raw$sigma) & params$raw$sigma > 0
  amp_emp            <- params$raw$A[estim_valid]
  sigma_rhythmic_emp <- params$raw$sigma[estim_valid]
  phase_emp          <- params$raw$phi[estim_idx & !is.na(params$raw$phi)]

  # Cap differential proportions at the estimated rhythmic budget.
  # CircadianBioOptions requires prop_DR + prop_DP + prop_DM <= prop_rhythmic
  # because every differential gene must be rhythmic in at least one group.
  total_diff <- prop_DR + prop_DP + prop_DM
  if (total_diff > params$prop_rhythmic && total_diff > 0) {
    scale_factor <- params$prop_rhythmic / total_diff
    if (verbose) {
      message(sprintf(
        paste0("estCircadianParam: prop_DR+prop_DP+prop_DM (%.3f) exceeds estimated ",
               "prop_rhythmic (%.3f). Scaling differential props by %.3f to fit budget."),
        total_diff, params$prop_rhythmic, scale_factor))
    }
    prop_DR <- prop_DR * scale_factor
    prop_DP <- prop_DP * scale_factor
    prop_DM <- prop_DM * scale_factor
  }

  obj <- CircadianBioOptions(
    ngenes = ngenes,
    prop_rhythmic = params$prop_rhythmic,
    period = period,
    lBaselineExpr = lBaselineExpr_emp,
    lOD = lOD_emp,
    amplitude = amp_emp,
    sigma_rhythmic = sigma_rhythmic_emp,
    paired_sigma = paired_sigma,
    cts = times,
    phase = phase_emp,
    prop_DR = prop_DR,
    prop_DP = prop_DP,
    prop_DM = prop_DM,
    mesor_diff = mesor_diff,
    phase_diff = phase_diff,
    amp_diff = amp_diff,
    dp_shift_mode = dp_shift_mode,
    dr_amp_scale = dr_amp_scale,
    dr_sigma_scale = dr_sigma_scale,
    sim.seed = sim.seed
  )

  # Attach the per-gene rhythm table so the threshold (alpha_pilot) can be
  # re-selected at load time (see scp_load_pilot()).  Stored capped at
  # `pilot_rhythm_cap` and sorted ascending by p-value, so any alpha <= cap
  # is a prefix slice; ngenes (above) is the prop_rhythmic denominator, so
  # genes with p >= cap can be dropped without affecting the proportion.
  obj$rhythm_fit  <- .build_rhythm_fit(params$raw)
  obj$pilot_cap   <- .pilot_rhythm_cap
  obj$alpha_pilot <- min_rhythm_pval
  obj$pilot_top_k <- 300L     # top-K cap used to form the effect-size distribution
  obj
}

# Cap on the stored per-gene rhythm table.  alpha_pilot is re-selectable in
# (0, .pilot_rhythm_cap]; genes with raw cosinor p >= this are never rhythmic
# at any supported alpha and are dropped to keep the pilot rds small.
.pilot_rhythm_cap <- 0.2

#' Build the capped, p-sorted per-gene rhythm table for a pilot
#'
#' @param raw The \code{$raw} list returned by \code{estimate_circadian_params}
#'   (fields M, A, phi, sigma, pvalue across all genes).
#' @param cap Upper p-value cap for retained genes (default
#'   \code{.pilot_rhythm_cap}).
#' @return A data.frame with columns \code{pvalue}, \code{A}, \code{phi},
#'   \code{sigma}, one row per gene with \code{pvalue < cap} and a finite,
#'   positive amplitude, sorted ascending by \code{pvalue}.
#' @keywords internal
.build_rhythm_fit <- function(raw, cap = .pilot_rhythm_cap) {
  keep <- is.finite(raw$pvalue) & raw$pvalue < cap &
          is.finite(raw$A) & raw$A > 0 &
          is.finite(raw$sigma) & raw$sigma > 0
  gene <- if (!is.null(raw$gene)) raw$gene[keep] else as.character(which(keep))
  # mesor (M) is power-irrelevant (the cosinor F-test is mesor-invariant, and the
  # simulator draws baseline from the lBaselineExpr pool) but is stored per gene so
  # the gene-level cosinor display can show absolute expression level. NA if the
  # fit did not return M.
  mesor <- if (!is.null(raw$M)) raw$M[keep] else rep(NA_real_, sum(keep))
  df <- data.frame(
    gene   = gene,
    pvalue = raw$pvalue[keep],
    A      = raw$A[keep],
    phi    = raw$phi[keep],
    sigma  = raw$sigma[keep],
    mesor  = mesor,
    stringsAsFactors = FALSE
  )
  df[order(df$pvalue), , drop = FALSE]
}

#' Build the capped, p-sorted per-gene rhythm table for a two-harmonic (K=2) pilot
#'
#' @description Two-harmonic analog of \code{\link{.build_rhythm_fit}}. The
#'   threshold column \code{pvalue} holds the K=1 cosinor F-test p-value
#'   (\code{p_K1}) so that \code{alpha_pilot} has identical "which genes are
#'   rhythmic" semantics in K=1 and K=2 modes; the extra columns carry the
#'   gene-paired 2nd-harmonic parameters so the joint
#'   \eqn{(A_1, \phi_1, A_2, \phi_2, \sigma)} draw is preserved when re-sliced.
#' @param p_K1,A1,phi1,A2,phi2,sigma Per-gene vectors (length = ngenes) from the
#'   two-harmonic fit.
#' @param cap Upper p-value cap (default \code{.pilot_rhythm_cap}).
#' @return A data.frame with columns \code{pvalue, A, phi, A2, phi2, sigma},
#'   one row per gene with \code{p_K1 < cap} and a finite positive \code{A1} and
#'   \code{sigma}, sorted ascending by \code{pvalue}. Presence of the \code{A2}
#'   column is the signal that downstream re-slicing should treat the pilot as
#'   two-harmonic.
#' @keywords internal
.build_rhythm_fit2h <- function(p_K1, A1, phi1, A2, phi2, sigma,
                                gene = NULL, mesor = NULL, cap = .pilot_rhythm_cap) {
  keep <- is.finite(p_K1) & p_K1 < cap &
          is.finite(A1) & A1 > 0 &
          is.finite(sigma) & sigma > 0
  g <- if (!is.null(gene)) gene[keep] else as.character(which(keep))
  m <- if (!is.null(mesor)) mesor[keep] else rep(NA_real_, sum(keep))
  df <- data.frame(
    gene   = g,
    pvalue = p_K1[keep],
    A      = A1[keep],
    phi    = phi1[keep],
    A2     = A2[keep],
    phi2   = phi2[keep],
    sigma  = sigma[keep],
    mesor  = m,
    stringsAsFactors = FALSE
  )
  df[order(df$pvalue), , drop = FALSE]
}


#' Estimate CircadianBioOptions Under the Two-Harmonic Cosinor Model
#'
#' @description Two-harmonic analog of \code{\link{estCircadianParam}}. Fits the
#' regression
#' \deqn{y = \mu + a_1 \cos(\omega_0 t) + b_1 \sin(\omega_0 t)
#'             + a_2 \cos(2\omega_0 t) + b_2 \sin(2\omega_0 t) + \epsilon}
#' to each gene by OLS, tests the joint rhythmicity null
#' \eqn{H_0: a_1 = b_1 = a_2 = b_2 = 0} on 4 d.f., and returns a
#' \code{CircadianBioOptions} object whose 5-way joint empirical distribution
#' \eqn{F_{(\mu, A_1, \phi_1, A_2, \phi_2, \sigma)}} is preserved by gene-paired
#' index draws downstream.
#'
#' The returned object carries the standard 1-harmonic fields
#' (\code{amplitude}, \code{phase}, \code{sigma_rhythmic}) holding the 1st-harmonic
#' \eqn{A_{g,1}}, \eqn{\phi_{g,1}}, and \eqn{\sigma_g} respectively, and three
#' extra fields attached to the list after construction:
#' \itemize{
#'   \item \code{amplitude2}: 2nd-harmonic amplitudes \eqn{A_{g,2}}
#'   \item \code{phase2}:     2nd-harmonic acrophases \eqn{\phi_{g,2}} in
#'                            hours \eqn{[0, period/2)}
#'   \item \code{paired_2h}:  \code{TRUE} flag (downstream simulators read this
#'                            to switch on the 2-harmonic generative model).
#' }
#' These extras live alongside \code{$diagnostics}; the constructor is not
#' modified.
#'
#' @section Methodological notes:
#' \enumerate{
#'   \item Pilot summary is a 5-way joint empirical distribution. When the
#'         simulator draws a rhythmic gene, the same index \eqn{j} is used
#'         across \eqn{(A_1, \phi_1, A_2, \phi_2, \sigma)}, preserving the
#'         cross-parameter dependence observed in the pilot.
#'   \item No 2nd-harmonic significance threshold is applied: ALL top-K genes
#'         passing the joint 4-d.f. F-test contribute, even those whose
#'         \eqn{A_2} is essentially noise. The empirical \eqn{F_{A_2}} will
#'         have a noise spike near zero, which is fine.
#'   \item K=2 extension is single-cohort only -- the differential simulator
#'         (\code{simCircadianDiff}) is not extended.
#' }
#'
#' @param data Numeric matrix of pilot expression (genes x samples).
#' @param times Numeric vector of sample collection times in hours,
#'   length \code{ncol(data)}.
#' @param period Period in the same units as \code{times} (default 24 hours).
#' @param min_rhythm_pval P-value threshold for the joint 4-d.f. F-test
#'   (default 0.01). Genes with p < this define the candidate set.
#' @param top_k Cap on number of genes contributing to the joint
#'   \eqn{(A_1, \phi_1, A_2, \phi_2, \sigma)} distribution (default 300),
#'   ranked by F-test p-value.
#' @param prop_DR,prop_DP,prop_DM Differential proportions (defaults 0.15,
#'   0.10, 0.00). Forwarded unchanged to \code{CircadianBioOptions}; note the
#'   K=2 extension itself is single-cohort only.
#' @param mesor_diff,phase_diff,amp_diff Differential effect-size ranges.
#' @param dp_shift_mode "fixed" or "uniform" (default "fixed").
#' @param dr_amp_scale,dr_sigma_scale DR scale factors (default 1).
#' @param paired_sigma Logical (default FALSE). Forwarded to
#'   \code{CircadianBioOptions}.
#' @param sim.seed Integer seed (default 12345).
#' @param verbose Logical (default TRUE). Print diagnostic summary.
#'
#' @return A \code{CircadianBioOptions} object with extra fields
#'   \code{$amplitude2}, \code{$phase2}, \code{$paired_2h}, and a
#'   \code{$diagnostics} list summarising the 2H fit (median A2/A1 ratio,
#'   fraction of genes with significant 2nd-harmonic block, etc.).
#'
#' @seealso \code{\link{estCircadianParam}} for the 1-harmonic case,
#'   \code{\link{simCircadianSingleCohort2H}} for the matching simulator.
#'
#' @export
estCircadianParam2H <- function(data, times, period = 24,
                                 min_rhythm_pval = 0.01,
                                 top_k = 500,
                                 prop_DR = 0.15, prop_DP = 0.10, prop_DM = 0.00,
                                 mesor_diff = c(0.5, 2.0),
                                 phase_diff = c(-6, 6),
                                 amp_diff = c(0.5, 2),
                                 dp_shift_mode = c("fixed", "uniform"),
                                 dr_amp_scale = 1.0,
                                 dr_sigma_scale = 1.0,
                                 paired_sigma = FALSE,
                                 sim.seed = 12345,
                                 verbose = TRUE) {

  # Input validation mirrors estCircadianParam.
  if (is.data.frame(data)) data <- as.matrix(data)
  if (!is.matrix(data) || !is.numeric(data))
    stop("'data' must be a numeric matrix (genes x samples). Got: ",
         paste(class(data), collapse = "/"))
  if (length(times) != ncol(data))
    stop(sprintf(
      "length(times) (%d) must equal ncol(data) (%d). Did you mean to transpose your matrix?",
      length(times), ncol(data)))

  dp_shift_mode <- match.arg(dp_shift_mode)

  ngenes  <- nrow(data)
  N       <- ncol(data)
  omega_0 <- 2 * pi / period
  top_k   <- as.integer(top_k)

  if (verbose) {
    cat("=== Estimating Two-Harmonic Cosinor Parameters from Pilot Data ===\n")
    cat("Genes:", ngenes, "\n")
    cat("Samples:", N, "\n")
    cat("Time points:", length(unique(times)), "\n\n")
  }

  # Build the K=2 cosinor design matrix once: columns
  #   1   cos(w t)   sin(w t)   cos(2w t)   sin(2w t)
  X2 <- cbind(
    1,
    cos(omega_0 * times), sin(omega_0 * times),
    cos(2 * omega_0 * times), sin(2 * omega_0 * times)
  )
  # Reduced (null) design: intercept only.
  X0 <- matrix(1, nrow = N, ncol = 1)

  # Fit per gene by OLS. For numerical robustness use lm.fit (QR).
  # Vectors: M, a1, b1, a2, b2, sigma_hat, K=2 F-test pvalue, K=1 F-test pvalue.
  Mhat  <- rep(NA_real_, ngenes)
  a1hat <- rep(NA_real_, ngenes); b1hat <- rep(NA_real_, ngenes)
  a2hat <- rep(NA_real_, ngenes); b2hat <- rep(NA_real_, ngenes)
  sigma_hat <- rep(NA_real_, ngenes)
  pvals     <- rep(NA_real_, ngenes)   # K=2 joint F-test (diagnostic)
  p_K1      <- rep(NA_real_, ngenes)   # K=1 cosinor F-test (used for top-K ranking)
  # Secondary p-value: test of 2nd-harmonic block alone (a2 = b2 = 0)
  # conditioned on the 1st harmonic, used for diagnostics only.
  p_2h_only <- rep(NA_real_, ngenes)

  df_resid_full <- N - 5L
  df_resid_null <- N - 1L
  if (df_resid_full < 1L)
    stop(sprintf("Two-harmonic cosinor needs N >= 6 samples; got N=%d.", N))

  X1 <- cbind(1,
              cos(omega_0 * times), sin(omega_0 * times))  # 1-harmonic design
  df_resid_1h <- N - 3L

  for (g in seq_len(ngenes)) {
    y <- as.numeric(data[g, ])
    fit2 <- tryCatch(stats::lm.fit(X2, y), error = function(e) NULL)
    if (is.null(fit2) || any(is.na(fit2$coefficients))) next

    coefs <- fit2$coefficients
    res   <- fit2$residuals
    RSS_full <- sum(res * res)
    RSS_null <- sum((y - mean(y))^2)

    # Joint K=2 F-test: H0: a1=b1=a2=b2=0 (4 d.f.); kept for diagnostics
    if (RSS_full <= 0 || df_resid_full < 1L) next
    Fstat <- ((RSS_null - RSS_full) / 4) / (RSS_full / df_resid_full)
    pv    <- stats::pf(Fstat, 4, df_resid_full, lower.tail = FALSE)

    # K=1 cosinor F-test (used for top-K ranking, matching Section 2.1)
    fit1 <- tryCatch(stats::lm.fit(X1, y), error = function(e) NULL)
    pv1  <- NA_real_; pv2 <- NA_real_
    if (!is.null(fit1) && !any(is.na(fit1$coefficients))) {
      RSS_1h <- sum(fit1$residuals^2)
      if (RSS_1h > 0 && df_resid_1h >= 1L) {
        F_1h <- ((RSS_null - RSS_1h) / 2) / (RSS_1h / df_resid_1h)
        pv1  <- stats::pf(F_1h, 2, df_resid_1h, lower.tail = FALSE)
      }
      # Optional diagnostic: 2nd-harmonic-only F-test (2 d.f.)
      if (RSS_full > 0 && df_resid_full >= 1L && RSS_1h >= RSS_full) {
        F_2h <- ((RSS_1h - RSS_full) / 2) / (RSS_full / df_resid_full)
        pv2  <- stats::pf(F_2h, 2, df_resid_full, lower.tail = FALSE)
      }
    }

    Mhat[g]  <- coefs[1L]
    a1hat[g] <- coefs[2L]; b1hat[g] <- coefs[3L]
    a2hat[g] <- coefs[4L]; b2hat[g] <- coefs[5L]
    sigma_hat[g] <- sqrt(RSS_full / df_resid_full)
    pvals[g]     <- pv
    p_K1[g]      <- pv1
    p_2h_only[g] <- pv2
  }

  # (A_k, phi_k) conversion
  A1_g  <- sqrt(a1hat^2 + b1hat^2)
  A2_g  <- sqrt(a2hat^2 + b2hat^2)
  # phi_1 wrapped to [0, period); phi_2 wrapped to [0, period/2)
  phi1_g <- (atan2(b1hat, a1hat) / omega_0) %% period
  phi2_g <- (atan2(b2hat, a2hat) / (2 * omega_0)) %% (period / 2)

  # Candidate rhythmic set + top-K selection (Option B: K=1 F-test ranking).
  # Pre-screen by K=1 cosinor F-test at p < min_rhythm_pval to match the
  # rhythmic-gene definition used in Section 2.1; rank by K=1 F-test p-value
  # (monotone with r1 = A1/sigma at fixed N) to take the strongest-cosinor
  # genes. K=2 fits (A1, A2, phi1, phi2, sigma) are still stored for the
  # joint 5-way pilot summary.
  rhythmic_cand <- !is.na(p_K1) & p_K1 < min_rhythm_pval &
                   !is.na(A1_g) & A1_g > 0 &
                   !is.na(sigma_hat) & sigma_hat > 0
  n_cand <- sum(rhythmic_cand, na.rm = TRUE)
  prop_rhythmic_emp <- n_cand / ngenes

  K_use <- min(top_k, n_cand)
  if (K_use <= 0) {
    stop("No rhythmic genes passed the K=1 cosinor F-test at p < ",
         min_rhythm_pval)
  }
  cand_idx <- which(rhythmic_cand)
  top_idx  <- cand_idx[order(p_K1[cand_idx])][seq_len(K_use)]

  # Empirical distributions
  lBaselineExpr_emp <- Mhat[!is.na(Mhat)]
  sigma_valid       <- sigma_hat[!is.na(sigma_hat) & sigma_hat > 0]
  lOD_emp           <- log(sigma_valid)

  # JOINT 5-way tuples on top-K, paired by gene index
  A1_emp     <- A1_g[top_idx]
  phi1_emp   <- phi1_g[top_idx]
  A2_emp     <- A2_g[top_idx]
  phi2_emp   <- phi2_g[top_idx]
  sigma_emp  <- sigma_hat[top_idx]

  # Cap differential proportions at the rhythmic budget.
  total_diff <- prop_DR + prop_DP + prop_DM
  if (total_diff > prop_rhythmic_emp && total_diff > 0) {
    scale_factor <- prop_rhythmic_emp / total_diff
    if (verbose) {
      message(sprintf(
        paste0("estCircadianParam2H: prop_DR+prop_DP+prop_DM (%.3f) exceeds estimated ",
               "prop_rhythmic (%.3f). Scaling differential props by %.3f to fit budget."),
        total_diff, prop_rhythmic_emp, scale_factor))
    }
    prop_DR <- prop_DR * scale_factor
    prop_DP <- prop_DP * scale_factor
    prop_DM <- prop_DM * scale_factor
  }

  opts <- CircadianBioOptions(
    ngenes         = ngenes,
    prop_rhythmic  = prop_rhythmic_emp,
    period         = period,
    lBaselineExpr  = lBaselineExpr_emp,
    lOD            = lOD_emp,
    amplitude      = A1_emp,
    sigma_rhythmic = sigma_emp,
    paired_sigma   = paired_sigma,
    cts            = times,
    phase          = phi1_emp,
    prop_DR        = prop_DR,
    prop_DP        = prop_DP,
    prop_DM        = prop_DM,
    mesor_diff     = mesor_diff,
    phase_diff     = phase_diff,
    amp_diff       = amp_diff,
    dp_shift_mode  = dp_shift_mode,
    dr_amp_scale   = dr_amp_scale,
    dr_sigma_scale = dr_sigma_scale,
    sim.seed       = sim.seed
  )
  # Restore the joint 5-way pilot pairing.
  # CircadianBioOptions() independently resamples amplitude/phase/sigma down
  # to length n_rhythmic via setAmplitude/setPhase, which destroys the
  # (A1, phi1, A2, phi2, sigma) tuple structure required by the K=2
  # simulator. We overwrite the five paired fields with the raw length-K_use
  # pilot vectors so that simCircadianSingleCohort2H draws a single shared
  # ji index of length n_rhythmic and applies it to all five vectors.
  opts$amplitude      <- A1_emp
  opts$phase          <- phi1_emp
  opts$sigma_rhythmic <- sigma_emp
  opts$amplitude2     <- A2_emp
  opts$phase2         <- phi2_emp
  opts$paired_2h      <- TRUE

  # Diagnostics
  n_2h_signif <- sum(!is.na(p_2h_only[top_idx]) & p_2h_only[top_idx] < 0.05,
                     na.rm = TRUE)
  A2_A1_ratio <- A2_emp / pmax(A1_emp, 1e-12)
  opts$diagnostics <- list(
    n_pre_screen     = n_cand,
    top_k_used       = K_use,
    prop_rhythmic    = prop_rhythmic_emp,
    A1_median        = stats::median(A1_emp, na.rm = TRUE),
    A2_median        = stats::median(A2_emp, na.rm = TRUE),
    A2_over_A1_med   = stats::median(A2_A1_ratio, na.rm = TRUE),
    n_2h_signif_p05  = n_2h_signif,
    prop_2h_signif   = n_2h_signif / max(K_use, 1L),
    sigma_median     = stats::median(sigma_emp, na.rm = TRUE),
    p_2h_only        = p_2h_only[top_idx],
    A1_emp           = A1_emp, A2_emp = A2_emp,
    phi1_emp         = phi1_emp, phi2_emp = phi2_emp,
    sigma_emp        = sigma_emp,
    screen_method    = "joint_2H_F4"
  )

  if (verbose) {
    cat(sprintf("Joint 4-d.f. F-test screen: %d/%d genes pass p < %.2g\n",
                n_cand, ngenes, min_rhythm_pval))
    cat(sprintf("Top-K used for joint (A1, phi1, A2, phi2, sigma): %d\n",
                K_use))
    cat(sprintf("Median A1 = %.3f   Median A2 = %.3f   Median A2/A1 = %.3f\n",
                opts$diagnostics$A1_median,
                opts$diagnostics$A2_median,
                opts$diagnostics$A2_over_A1_med))
    cat(sprintf("Genes with 2nd-harmonic-block p < 0.05: %d/%d (%.1f%%)\n",
                n_2h_signif, K_use,
                100 * opts$diagnostics$prop_2h_signif))
    cat(sprintf("Median residual sigma: %.3f\n", opts$diagnostics$sigma_median))
    cat(sprintf("Proportion rhythmic (joint F-test): %.1f%%\n",
                100 * prop_rhythmic_emp))
  }

  # Attach the wider per-gene rhythm table so alpha_pilot is re-selectable at
  # load time (scp_load_pilot). Thresholded on the K=1 cosinor p (p_K1), the
  # same statistic that defines the candidate set above, so alpha_pilot means
  # the same thing in K=1 and K=2 modes.
  # gene IDs from the source matrix rownames, aligned to the per-gene vectors;
  # length-guarded so a subset/transform never produces a wrong mapping.
  gene_ids_2h <- if (!is.null(rownames(data)) && length(p_K1) == nrow(data)) rownames(data) else NULL
  opts$rhythm_fit  <- .build_rhythm_fit2h(p_K1, A1_g, phi1_g, A2_g, phi2_g, sigma_hat,
                                          gene = gene_ids_2h, mesor = Mhat)
  opts$pilot_cap   <- .pilot_rhythm_cap
  opts$alpha_pilot <- min_rhythm_pval
  opts$pilot_top_k <- top_k

  opts
}


#' Estimate CircadianBioOptions with FMM Per-Gene Parameters
#'
#' @description Like \code{\link{estCircadianParam}}, but additionally fits the
#' FMM (Frequency Modulated Mobius) waveform per top-K rhythmic gene to obtain
#' empirical \eqn{\hat\omega_g} (shape) and \eqn{\hat\alpha_g} (peak location,
#' radians). Returns a \code{CircadianBioOptions} populated with paired
#' \code{omega_rhythmic} and \code{alpha_rhythmic} vectors (paired with
#' amplitude and sigma_rhythmic by gene index). A \code{$diagnostics} list
#' carries empirical fit summaries: Beta(1, beta) MoM estimate, sigma_alpha (standard
#' deviation of alpha in hours), median R^2 of the FMM fits, and the raw
#' per-gene vectors for diagnostic plotting.
#'
#' Used by Fig 5 (B-vs-m grid under empirical FMM truth) and the diagnostic
#' supplementary figure. Fig 4 (sensitivity sweep) does NOT need this -- it
#' uses synthetic \code{omega_dist}/\code{alpha_dist} draws.
#'
#' @param data Gene expression matrix (genes x samples).
#' @param times Sample collection times (length = ncol(data)).
#' @param period Period in hours (default 24).
#' @param min_rhythm_pval detect_FMM pre-screen p-value threshold (default 0.01).
#'   Genes with detect_FMM p < this value are classified as rhythmic candidates
#'   that contribute to the empirical (A, sigma, omega, alpha) distributions
#'   and to \code{prop_rhythmic}. Matches the downstream detect_FMM detector
#'   (so simulation truth and detection use the same model).
#' @param top_k Cap on number of FMM-fitted genes (default 300, matches
#'   estCircadianParam's top-K logic).
#' @param paired_sigma Logical (default TRUE). Forwarded to CircadianBioOptions.
#' @param paired_omega Logical (default TRUE). Forwarded.
#' @param paired_alpha Logical (default TRUE). Forwarded.
#' @param mc.cores Parallel cores for FMM fitting (default 1).
#' @param verbose Print progress (default TRUE).
#' @param ... Additional args passed to CircadianBioOptions().
#'
#' @return CircadianBioOptions object with:
#'   - amplitude, sigma_rhythmic, omega_rhythmic, alpha_rhythmic vectors
#'     of equal length (= number of successfully fitted top-K genes).
#'     \strong{omega_rhythmic} is the FMM shape parameter in [0, 1];
#'     \strong{alpha_rhythmic} is the FMM peak location in radians [0, 2*pi].
#'   - $diagnostics list: beta_hat, sigma_alpha_hat (hours), R2_median,
#'     omega_emp, alpha_emp, n_fitted, screen_method = "detect_FMM".
#'
#' @section FMM scale convention:
#' This function fits \code{FMM::fitFMM()} on \emph{radians} timepoints
#' \eqn{t \in [0, 2\pi]}, because \code{FMM::fitFMM()} does NOT internally
#' rescale \code{timePoints}. Empirically, on a known signal with true
#' \eqn{\omega = 0.3}: fitting in radians recovers \eqn{\hat\omega = 0.30};
#' fitting on \eqn{[0, 1]} gives \eqn{\hat\omega = 0.046} (off by ~7x and
#' meaningless as an FMM shape parameter). Radians is also the convention
#' used by the simulator (\code{simCircadianFMM}: \code{cts_rad = cts * 2*pi/24}
#' and \code{FMM::generateFMM(..., to = 2*pi)}), so the (omega, alpha)
#' returned here are directly usable as simulation truth.
#'
#' Note: the pre-screen below uses the K=2 harmonic LRT (\code{detect_FMM}),
#' which has the same null hypothesis as the downstream detector. This keeps
#' the truth-label generator and the detector consistent. (omega, alpha)
#' from the FMM fit are propagated only as simulation truth, not detection.
#'
#' @export
estCircadianParamFMM <- function(data, times, period = 24,
                                 min_rhythm_pval = 0.01,
                                 top_k = 300L,
                                 screen_K = 2L,
                                 paired_sigma = TRUE,
                                 paired_omega = TRUE,
                                 paired_alpha = TRUE,
                                 prop_DR = 0.0,
                                 prop_DP = 0.0,
                                 prop_DM = 0.0,
                                 mc.cores = 1L,
                                 verbose = TRUE,
                                 ...) {

  if (!requireNamespace("FMM", quietly = TRUE))
    stop("estCircadianParamFMM requires the FMM package. Install with install.packages('FMM').")

  # Step 1a: cosinor fit on all genes — needed for residual sigma estimates
  # (used downstream as F_sigma in simulation). The cosinor p-value is NOT used
  # for screening; the detect_FMM screen below replaces it. This keeps the noise
  # estimate consistent with previous pipelines while letting detect_FMM define
  # which genes count as "rhythmic", internally consistent with the
  # K-harmonic F-test (detect_FMM) detector that downstream Fig 4/5 use.
  if (verbose) cat("Step 1a: cosinor fit (for sigma estimates)...\n")
  pdf_ <- fitCosinorAll(data, times, period = period,
                        min_rhythm_pval = 1.0)   # no cosinor screen here

  # Step 1b: detect_FMM pre-screen (same null as downstream detector).
  # Tests H0: a_1=b_1=...=a_K=b_K=0 via exact F(2K, n-2K-1). screen_K is
  # exposed as a parameter so the pre-screen K can match the analysis K.
  screen_K <- as.integer(screen_K)
  if (verbose) cat(sprintf("Step 1b: K-harmonic pre-screen (K=%d) on %d genes...\n",
                            screen_K, nrow(data)))
  fmm_pvals <- detect_cosinor(data, times, K = screen_K, period = period)

  rhy <- which(!is.na(fmm_pvals) & fmm_pvals < min_rhythm_pval &
               !is.na(pdf_$sigma) & pdf_$sigma > 0)
  if (length(rhy) == 0)
    stop("No rhythmic genes passed detect_FMM pre-screen at p < ", min_rhythm_pval)

  # Top-K selection by detect_FMM p-value
  n_screen_pass <- length(rhy)
  if (length(rhy) > top_k) {
    rhy <- rhy[order(fmm_pvals[rhy])][seq_len(top_k)]
  }
  if (verbose) cat(sprintf("  detect_FMM screen: %d genes pass p<%.2g -> top-%d for FMM-fit summary\n",
                            n_screen_pass, min_rhythm_pval, length(rhy)))

  # Step 2: fit FMM per top-K gene to extract (omega, alpha, A) parameters.
  #
  # SCALE CONVENTION: timePoints in radians [0, 2*pi]. This matters because
  # FMM::fitFMM does NOT internally rescale timePoints — the returned omega
  # and alpha are in the input scale. Verified empirically: on a known
  # signal with true omega=0.3, fitting in radians recovers omega=0.30,
  # while fitting on [0,1] gives omega=0.046 (off by ~7x and meaningless
  # as an FMM shape parameter).
  #
  # The simulator (simCircadianFMM, see simulation.R cts_rad and FMM::generateFMM
  # with to=2*pi) expects omega/alpha in radians, so radians is correct here.
  #
  # Note: the pre-screen uses detect_FMM (linear K-harmonic F-test) which is
  # scale-invariant in time. The FMM model fit below uses radian timepoints
  # to match the simulator's convention. Cross-pipeline consistency holds
  # because alpha/omega from the FMM fit are propagated only as simulation
  # truth, not back into the detector.
  if (verbose) cat(sprintf("Step 2: fit FMM (mc.cores=%d)...\n", mc.cores))
  tod_rad <- (times %% period) * 2 * pi / period

  fit_one <- function(g) {
    y <- as.numeric(data[g, ])
    tryCatch({
      fit <- FMM::fitFMM(y, timePoints = tod_rad,
                         nback           = 1L,             # match detector
                         lengthAlphaGrid = 12L,
                         lengthOmegaGrid = 12L,
                         showProgress    = FALSE)
      list(omega = fit@omega, alpha = fit@alpha,
           A = fit@A, R2 = fit@R2)
    }, error = function(e) NULL)
  }
  fmm_fits <- if (mc.cores > 1L) {
    parallel::mclapply(rhy, fit_one, mc.cores = mc.cores)
  } else {
    lapply(rhy, fit_one)
  }

  ok_idx <- which(!sapply(fmm_fits, is.null))
  if (length(ok_idx) == 0)
    stop("All FMM fits failed; check input data and FMM package.")
  if (verbose) cat(sprintf("  FMM fit: %d/%d succeeded\n",
                            length(ok_idx), length(rhy)))

  # Aligned vectors (only successfully fitted genes)
  rhy_ok <- rhy[ok_idx]
  fits_ok <- fmm_fits[ok_idx]

  amp_emp            <- pdf_$A[rhy_ok]                          # cosinor amplitude
  sigma_rhythmic_emp <- pdf_$sigma[rhy_ok]                      # cosinor residual sd
  phase_emp          <- pdf_$phi[rhy_ok]                        # cosinor acrophase (hours)
  omega_emp          <- vapply(fits_ok, `[[`, numeric(1), "omega")
  alpha_emp          <- vapply(fits_ok, `[[`, numeric(1), "alpha")
  R2_emp             <- vapply(fits_ok, `[[`, numeric(1), "R2")

  # Step 3: empirical diagnostics
  beta_hat        <- (1 - mean(omega_emp)) / mean(omega_emp)    # MoM for Beta(1, β)
  sigma_alpha_hat <- sd(alpha_emp) * 24 / (2 * pi)              # in hours
  R2_median       <- median(R2_emp, na.rm = TRUE)

  if (verbose) {
    cat(sprintf("  Diagnostics: beta_hat=%.3f  sigma_alpha_hat=%.3fh  R2_median=%.3f\n",
                beta_hat, sigma_alpha_hat, R2_median))
  }

  # Step 4: assemble baseline (mesor and lOD distributions across ALL genes,
  # same as estCircadianParam — these aren't FMM-specific)
  lBaselineExpr_emp <- pdf_$M[!is.na(pdf_$M)]
  sigma_valid       <- pdf_$sigma[!is.na(pdf_$sigma) & pdf_$sigma > 0]
  lOD_emp           <- log(sigma_valid)

  # prop_rhythmic from detect_FMM screen (matches detector used downstream)
  # — full screen-pass count, not just top-K
  prop_rhythmic_emp <- mean(!is.na(fmm_pvals) & fmm_pvals < min_rhythm_pval,
                             na.rm = TRUE)

  # Auto-scale differential props if they exceed the rhythmic budget,
  # mirroring the logic in estCircadianParam.
  total_diff <- prop_DR + prop_DP + prop_DM
  if (total_diff > prop_rhythmic_emp && total_diff > 0) {
    scale_factor <- prop_rhythmic_emp / total_diff
    if (verbose) {
      message(sprintf(
        paste0("estCircadianParamFMM: prop_DR+prop_DP+prop_DM (%.3f) exceeds estimated ",
               "prop_rhythmic (%.3f). Scaling differential props by %.3f to fit budget."),
        total_diff, prop_rhythmic_emp, scale_factor))
    }
    prop_DR <- prop_DR * scale_factor
    prop_DP <- prop_DP * scale_factor
    prop_DM <- prop_DM * scale_factor
  }

  opts <- CircadianBioOptions(
    ngenes         = nrow(data),
    prop_rhythmic  = prop_rhythmic_emp,
    period         = period,
    lBaselineExpr  = lBaselineExpr_emp,
    lOD            = lOD_emp,
    amplitude      = amp_emp,
    sigma_rhythmic = sigma_rhythmic_emp,
    omega_rhythmic = omega_emp,
    alpha_rhythmic = alpha_emp,
    paired_sigma   = paired_sigma,
    paired_omega   = paired_omega,
    paired_alpha   = paired_alpha,
    prop_DR        = prop_DR,
    prop_DP        = prop_DP,
    prop_DM        = prop_DM,
    cts            = times,
    phase          = phase_emp,
    ...
  )
  opts$diagnostics <- list(
    beta_hat        = beta_hat,
    sigma_alpha_hat = sigma_alpha_hat,
    R2_median       = R2_median,
    omega_emp       = omega_emp,
    alpha_emp       = alpha_emp,
    R2_emp          = R2_emp,
    n_fitted        = length(ok_idx),
    n_pre_screen    = sum(!is.na(fmm_pvals) & fmm_pvals < min_rhythm_pval, na.rm = TRUE),
    screen_method   = "detect_FMM",
    screen_K        = screen_K,
    top_k_used      = length(rhy)
  )
  opts
}


#' Estimate CircadianBioOptions from Two-Group Pilot Data
#'
#' @description Bridge function for two-group pilot data. Estimates circadian
#' parameters separately for each group, then derives differential simulation
#' hyperparameters (prop_DR, prop_DP, phase_diff, amp_diff) directly from the
#' empirical between-group differences. This provides data-driven starting
#' values for the simulation inputs, replacing the need to specify them by hand.
#'
#' @param data_1 Gene expression matrix (genes x samples) for group 1
#' @param data_2 Gene expression matrix (genes x samples) for group 2
#' @param times_1 Time points for group 1
#' @param times_2 Time points for group 2
#' @param period Circadian period (default 24)
#' @param min_rhythm_pval P-value threshold for rhythmic classification (default 0.01)
#' @param phase_shift_threshold Minimum |delta_phi| in hours to classify a
#'   jointly-rhythmic gene as differentially phased (default 2 h)
#' @param sim.seed Random seed
#' @param verbose Print diagnostic summary
#'
#' @return A CircadianBioOptions object whose prop_DR, prop_DP, phase_diff, and
#'   amp_diff are estimated from pilot data. A $diagnostics list is attached
#'   with all intermediate estimates for user inspection and reporting.
#'
#' @details
#' Estimation logic:
#' \itemize{
#'   \item prop_DR: fraction of all genes rhythmic in exactly one group
#'         (xor of per-gene rhythmicity calls).
#'   \item prop_DP: fraction of all genes that are (i) jointly rhythmic in
#'         both groups AND (ii) have |circular phase difference| > phase_shift_threshold.
#'   \item phase_diff: [Q25, Q75] of the signed empirical phase differences
#'         delta_phi = phi_2 - phi_1 (circular, in hours) among DP genes.
#'         Used as the Unif(phase_diff[1], phase_diff[2]) draw in simulation.
#'   \item amp_diff: [Q25, Q75] of the amplitude ratio A_2/A_1 among jointly
#'         rhythmic genes. Passed as amp_diff bounds.
#' }
#' Baseline distributions (lBaselineExpr, lOD, amplitude, phase) are taken
#' from group 1, consistent with estCircadianParam() for a single group.
circular_difference <- function(phi1, phi2, period = 24) {
  diff <- phi1 - phi2
  ((diff + period/2) %% period) - period/2
}

#' Estimate a two-group circadian pilot for differential power
#'
#' Fits cosinor models to two groups of pilot data and assembles the two-group
#' pilot summary used by the differential power simulator, encoding
#' differential rhythmicity (DR), differential phase (DP), and differential
#' MESOR (DM) structure.
#'
#' @param data_1,data_2 Gene-by-sample expression matrices for the two groups.
#' @param times_1,times_2 Numeric collection times for each group.
#' @param period Period in hours (default 24).
#' @param min_rhythm_pval Marginal p-value screen defining the rhythmic pilot set.
#' @param phase_shift_threshold Minimum phase difference (hours) for the DP class.
#' @param prop_DM Optional target proportion of differential-MESOR genes.
#' @param mesor_diff Optional MESOR difference for the DM component.
#' @param dp_shift_mode Differential-phase shift model, "uniform" or "fixed".
#' @param paired_sigma Preserve the paired (amplitude, noise) relationship across groups.
#' @param sim.seed Random seed.
#' @param verbose Print progress.
#' @return A two-group \code{CircadianBioOptions} pilot summary.
#' @export
estCircadianParamTwoGroup <- function(data_1, data_2, times_1, times_2,
                                      period = 24,
                                      min_rhythm_pval = 0.01,
                                      phase_shift_threshold = 2,
                                      prop_DM = NULL,
                                      mesor_diff = NULL,
                                      dp_shift_mode = c("uniform", "fixed"),
                                      paired_sigma = FALSE,
                                      sim.seed = 12345,
                                      verbose = TRUE) {

  if (nrow(data_1) != nrow(data_2))
    stop("data_1 and data_2 must have the same number of genes (rows).")

  if (verbose) cat("\n--- Fitting group 1 ---\n")
  p1 <- estimate_circadian_params(data_1, times_1, period = period,
                                   min_rhythm_pval = min_rhythm_pval,
                                   verbose = verbose)

  if (verbose) cat("\n--- Fitting group 2 ---\n")
  p2 <- estimate_circadian_params(data_2, times_2, period = period,
                                   min_rhythm_pval = min_rhythm_pval,
                                   verbose = verbose)

  ngenes     <- nrow(data_1)
  rhythmic_1 <- p1$raw$is_rhythmic
  rhythmic_2 <- p2$raw$is_rhythmic

  # prop_DR: rhythmic in exactly one group ------------------------------------
  dr_mask     <- xor(rhythmic_1, rhythmic_2) & !is.na(rhythmic_1) & !is.na(rhythmic_2)
  prop_DR_emp <- mean(dr_mask)

  # Union rhythmic: rhythmic in at least one group ----------------------------
  # This is the correct budget for CircadianBioOptions: every DR/DP gene is
  # rhythmic in at least one group, so the budget is union, not group-1 alone.
  union_mask     <- (rhythmic_1 | rhythmic_2) & !is.na(rhythmic_1) & !is.na(rhythmic_2)
  prop_union_rhy <- mean(union_mask)

  # Jointly rhythmic pool -----------------------------------------------------
  jointly    <- rhythmic_1 & rhythmic_2 & !is.na(rhythmic_1) & !is.na(rhythmic_2)
  prop_joint <- mean(jointly)

  # Phase differences among jointly rhythmic genes ----------------------------
  delta_phi <- circular_difference(
    p2$raw$phi[jointly],
    p1$raw$phi[jointly],
    period
  )

  # prop_DP: fraction of ALL genes that are jointly rhythmic AND |delta_phi| > threshold
  dp_among_joint <- abs(delta_phi) > phase_shift_threshold & !is.na(delta_phi)
  prop_DP_emp    <- if (length(dp_among_joint) > 0)
    prop_joint * mean(dp_among_joint, na.rm = TRUE) else 0

  # prop_DM: jointly rhythmic genes with significant mesor difference -----------
  # SE(M) ≈ sigma/sqrt(n) for OLS intercept; two-sample z-test on mesor difference
  n1 <- ncol(data_1); n2 <- ncol(data_2)
  se_mesor_diff <- sqrt(p1$raw$sigma^2 / n1 + p2$raw$sigma^2 / n2)
  valid_mesor   <- jointly & !is.na(p1$raw$M) & !is.na(p2$raw$M) &
                   !is.na(se_mesor_diff) & se_mesor_diff > 0
  z_DM          <- rep(NA_real_, length(jointly))
  z_DM[valid_mesor] <- (p2$raw$M[valid_mesor] - p1$raw$M[valid_mesor]) /
                        se_mesor_diff[valid_mesor]
  p_DM_g        <- 2 * pnorm(-abs(z_DM))
  dm_mask       <- jointly & !is.na(p_DM_g) & p_DM_g < min_rhythm_pval
  prop_DM_emp   <- if (!is.null(prop_DM)) prop_DM else mean(dm_mask, na.rm = TRUE)
  n_dm_pilot    <- sum(dm_mask, na.rm = TRUE)

  # mesor_diff: IQR of |M2 - M1| among detected DM genes (simulation draw range)
  dm_raw_diffs <- abs(p2$raw$M - p1$raw$M)[dm_mask & !is.na(p2$raw$M - p1$raw$M)]
  if (!is.null(mesor_diff)) {
    mesor_diff_emp <- mesor_diff
  } else if (length(dm_raw_diffs) >= 5) {
    mesor_diff_emp <- c(quantile(dm_raw_diffs, 0.25, names = FALSE),
                        quantile(dm_raw_diffs, 0.75, names = FALSE))
    if (diff(mesor_diff_emp) < 0.1)
      mesor_diff_emp <- c(mesor_diff_emp[1] * 0.9, mesor_diff_emp[2] * 1.1 + 0.1)
  } else {
    mesor_diff_emp <- c(0.3, 1.0)   # conservative fallback when too few DM genes
  }

  # Budget guard: cap prop_DR + prop_DP + prop_DM to fit within union rhythmic budget.
  # Floating-point arithmetic can cause the sum to equal or fractionally exceed
  # prop_union_rhy, which would error in CircadianBioOptions.
  total_diff <- sum(c(prop_DR_emp, prop_DP_emp, prop_DM_emp), na.rm = TRUE)
  if (total_diff >= prop_union_rhy && total_diff > 0) {
    scale_factor <- prop_union_rhy * 0.999 / total_diff
    prop_DR_emp  <- prop_DR_emp  * scale_factor
    prop_DP_emp  <- prop_DP_emp  * scale_factor
    prop_DM_emp  <- prop_DM_emp  * scale_factor
  }

  # phase_diff range: signed IQR among DP genes (fallback to +/- threshold)
  dp_vals <- delta_phi[dp_among_joint & !is.na(delta_phi)]
  if (length(dp_vals) >= 5) {
    phase_diff_emp <- c(quantile(dp_vals, 0.25, names = FALSE),
                        quantile(dp_vals, 0.75, names = FALSE))
    if (diff(phase_diff_emp) < 0.5)
      phase_diff_emp <- c(phase_diff_emp[1] - 0.25, phase_diff_emp[2] + 0.25)
  } else {
    phase_diff_emp <- c(-phase_shift_threshold, phase_shift_threshold)
  }

  # SNR (r = A/sigma) for each group's rhythmic genes -------------------------
  r1_vals <- p1$raw$r[rhythmic_1 & !is.na(p1$raw$r) & is.finite(p1$raw$r)]
  r2_vals <- p2$raw$r[rhythmic_2 & !is.na(p2$raw$r) & is.finite(p2$raw$r)]

  snr_summary <- function(x) {
    if (length(x) < 3) return(c(q25 = NA, median = NA, q75 = NA))
    c(q25    = quantile(x, 0.25, names = FALSE),
      median = median(x),
      q75    = quantile(x, 0.75, names = FALSE))
  }
  r1_snr <- snr_summary(r1_vals)
  r2_snr <- snr_summary(r2_vals)

  # Effective DP SNR: r_dp = 2 * r * |sin(pi * delta_phi / period)| ----------
  # Use group-1 SNR for jointly-rhythmic genes as the baseline r
  r_joint <- p1$raw$r[jointly]
  valid_dp <- dp_among_joint & !is.na(delta_phi) & !is.na(r_joint) & is.finite(r_joint)
  r_dp_eff <- 2 * r_joint[valid_dp] * abs(sin(pi * delta_phi[valid_dp] / period))
  r_dp_snr <- snr_summary(r_dp_eff)

  n_dp_pilot <- sum(valid_dp, na.rm = TRUE)

  # Warn if too few pilot genes for any endpoint to trust power estimates ------
  n_dr_pilot <- sum(dr_mask, na.rm = TRUE)
  min_pilot  <- 50L

  if (n_dr_pilot < min_pilot && verbose) {
    warning(sprintf(
      paste0("Only %d DR genes found in pilot (rhythmic in exactly one group) — ",
             "fewer than the recommended minimum of %d.\n",
             "  DR power will be near zero and estimates unreliable.\n",
             "  --> Remove 'DR' from test_types in CircadianDesignOptions() to skip this endpoint,\n",
             "      or use a larger pilot / relax min_rhythm_pval (currently %.3f)."),
      n_dr_pilot, min_pilot, min_rhythm_pval))
  }
  if (n_dp_pilot < min_pilot && verbose) {
    warning(sprintf(
      paste0("Only %d DP genes found in pilot (jointly rhythmic & |Δφ|>%.1fh) — ",
             "fewer than the recommended minimum of %d.\n",
             "  DP power will be near zero and estimates unreliable.\n",
             "  --> Remove 'DP' from test_types in CircadianDesignOptions() to skip this endpoint,\n",
             "      or reduce phase_shift_threshold (currently %.1fh) / use a larger pilot."),
      n_dp_pilot, phase_shift_threshold, min_pilot, phase_shift_threshold))
  }
  if (n_dm_pilot < min_pilot && verbose) {
    warning(sprintf(
      paste0("Only %d DM genes found in pilot (jointly rhythmic & significant mesor diff) — ",
             "fewer than the recommended minimum of %d.\n",
             "  DM power will be near zero and estimates unreliable.\n",
             "  --> Remove 'DM' from test_types in CircadianDesignOptions() to skip this endpoint,\n",
             "      or use a larger pilot / relax min_rhythm_pval (currently %.3f)."),
      n_dm_pilot, min_pilot, min_rhythm_pval))
  }

  amp_diff_emp <- c(0.5, 2.0)

  # Baseline distributions from group 1 (consistent with estCircadianParam) ---
  lBaselineExpr_emp <- p1$raw$M[!is.na(p1$raw$M)]
  sigma_valid       <- p1$raw$sigma[!is.na(p1$raw$sigma) & p1$raw$sigma > 0]
  lOD_emp           <- log(sigma_valid)
  # F_{A,sigma} and F_phi: top-K estimation set (G_R), paired for joint sampling
  estim_valid_1      <- p1$raw$in_estim_set & !is.na(p1$raw$A) & p1$raw$A > 0 &
                        !is.na(p1$raw$sigma) & p1$raw$sigma > 0
  amp_emp            <- p1$raw$A[estim_valid_1]
  sigma_rhythmic_emp <- p1$raw$sigma[estim_valid_1]
  phase_emp          <- p1$raw$phi[p1$raw$in_estim_set & !is.na(p1$raw$phi)]

  # Group-2 distributions for fully symmetric two-group simulation -----------
  rhythmic_idx_2     <- p2$raw$is_rhythmic
  amp_emp2           <- p2$raw$A[p2$raw$in_estim_set & !is.na(p2$raw$A) & p2$raw$A > 0]
  sigma_valid_2      <- p2$raw$sigma[!is.na(p2$raw$sigma) & p2$raw$sigma > 0]
  lOD_emp2           <- log(sigma_valid_2)
  lBaselineExpr2_emp <- p2$raw$M[!is.na(p2$raw$M)]   # group-2 mesor distribution

  # Diagnostics ---------------------------------------------------------------
  diagnostics <- list(
    # Simulation hyperparameters
    prop_DR_emp           = prop_DR_emp,
    prop_DP_emp           = prop_DP_emp,
    prop_DM_emp           = prop_DM_emp,
    phase_diff_emp        = phase_diff_emp,
    mesor_diff_emp        = mesor_diff_emp,
    phase_shift_threshold = phase_shift_threshold,
    n_DR_genes_pilot      = n_dr_pilot,
    n_DP_genes_pilot      = n_dp_pilot,
    n_DM_genes_pilot      = n_dm_pilot,
    # Rhythmicity summary
    prop_rhythmic_1       = p1$prop_rhythmic,
    prop_rhythmic_2       = p2$prop_rhythmic,
    prop_union_rhythmic   = prop_union_rhy,
    prop_jointly_rhythmic = prop_joint,
    # SNR guidance (key for interpreting power curves)
    r1_snr                = r1_snr,   # A/sigma for group-1 rhythmic genes
    r2_snr                = r2_snr,   # A/sigma for group-2 rhythmic genes
    r_dp_eff_snr          = r_dp_snr  # effective DP SNR = 2r|sin(pi*delta_phi/24)|
  )

  if (verbose) {
    cat("\n=== Two-Group Empirical Differential Parameter Estimates ===\n")
    cat(sprintf("  Group 1 rhythmic: %.1f%%   Group 2 rhythmic: %.1f%%\n",
                100 * p1$prop_rhythmic, 100 * p2$prop_rhythmic))
    cat(sprintf("  Jointly rhythmic: %.1f%%\n", 100 * prop_joint))
    cat(sprintf("  Estimated prop_DR  (rhythmic in exactly one group): %.4f\n",
                prop_DR_emp))
    cat(sprintf("  Estimated prop_DP  (jointly rhythmic & |Δφ|>%.1fh): %.4f\n",
                phase_shift_threshold, prop_DP_emp))
    cat(sprintf("  Estimated phase_diff IQR (DP genes, n=%d): [%.2f, %.2f] h\n",
                n_dp_pilot, phase_diff_emp[1], phase_diff_emp[2]))
    cat(sprintf("  Estimated prop_DM  (jointly rhythmic & sig. mesor diff, p<%.2f): %.4f  (n=%d genes)\n",
                min_rhythm_pval, prop_DM_emp, n_dm_pilot))
    cat(sprintf("  Estimated mesor_diff IQR (DM genes): [%.3f, %.3f] log-CPM units\n",
                mesor_diff_emp[1], mesor_diff_emp[2]))
    cat("\n  --- Signal-to-Noise Ratio guidance (r = A/sigma) ---\n")
    cat(sprintf("  Group 1 rhythmic genes  r: median=%.2f  IQR [%.2f, %.2f]\n",
                r1_snr["median"], r1_snr["q25"], r1_snr["q75"]))
    cat(sprintf("  Group 2 rhythmic genes  r: median=%.2f  IQR [%.2f, %.2f]\n",
                r2_snr["median"], r2_snr["q25"], r2_snr["q75"]))
    cat(sprintf("  Effective DP SNR (2r|sin(πΔφ/24)|): median=%.2f  IQR [%.2f, %.2f]\n",
                r_dp_snr["median"], r_dp_snr["q25"], r_dp_snr["q75"]))
    cat("  --> Compare these SNR values to the stratified power curves to\n")
    cat("      identify which fraction of your pilot genes will be detectable\n")
    cat("      at your planned sample size.\n")
    cat("\n  prop_DR, prop_DP, prop_DM, phase_diff, mesor_diff are used directly in CircadianBioOptions.\n")
    cat("  Inspect $diagnostics to review or override any estimate.\n")
  }

  opts <- CircadianBioOptions(
    ngenes          = ngenes,
    prop_rhythmic   = prop_union_rhy,
    period          = period,
    lBaselineExpr   = lBaselineExpr_emp,
    lBaselineExpr2  = lBaselineExpr2_emp,  # F̂_M2: group-2 mesor distribution
    lOD             = lOD_emp,
    lOD2            = lOD_emp2,            # F̂_σ2: group-2 noise distribution
    amplitude       = amp_emp,
    sigma_rhythmic  = sigma_rhythmic_emp,
    paired_sigma    = paired_sigma,
    cts             = times_1,             # F̂_TOD1: group-1 sampling time distribution
    amplitude2      = amp_emp2,            # F̂_A2: used for g2-only DR genes
    cts2            = times_2,             # F̂_TOD2: group-2 sampling time distribution
    phase           = phase_emp,
    prop_DR         = prop_DR_emp,
    prop_DP         = prop_DP_emp,
    prop_DM         = prop_DM_emp,
    mesor_diff      = mesor_diff_emp,
    phase_diff      = phase_diff_emp,
    amp_diff        = amp_diff_emp,
    dp_shift_mode   = match.arg(dp_shift_mode),
    sim.seed        = sim.seed
  )

  opts$diagnostics <- diagnostics
  opts
}
