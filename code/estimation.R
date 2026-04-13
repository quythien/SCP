#' Estimate Circadian Parameters from Pilot Data
#'
#' @description Estimate distributions of amplitude, phase, noise level, and
#' effect size from real circadian time-series data.
#'
#' @param data Gene expression matrix (genes x samples)
#' @param times Time points for each sample
#' @param period Period (default 24)
#' @param min_rhythm_pval P-value threshold for including rhythmic genes (default 0.1)
#' @param verbose Print progress
#'
#' @return List with parameter distributions
#' 
#' 
#' 
#' #' @description SEMIPARAMETRIC power analysis framework:
#' - Parametric cosinor model for mean function
#' - Empirical parameter distributions from pilot data
#' - Empirical Bayes variance estimation
#' - No assumptions on amplitude, phase, or noise distributions
#'
#' @details The framework is semiparametric because while the circadian
#' model structure (cosinor) is parametric, the parameter distributions
#' can be learned nonparametrically from pilot data via 
#' estimate_circadian_params().


estimate_circadian_params = function(data, times, period = 24,
                                     min_rhythm_pval = 0.1,
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
    sqrt(mean((y - yhat)^2, na.rm = TRUE))
  })

  # Effect size
  r_vals = A_vals / sigma_vals

  # Identify likely rhythmic genes for parameter estimation
  rhythmic_genes = pvals < min_rhythm_pval & !is.na(A_vals) & A_vals > 0

  if (verbose) {
    cat("Genes passing rhythm filter (p <", min_rhythm_pval, "):",
        sum(rhythmic_genes), "\n")
  }

  # Estimate parameter distributions
  params = list(
    # Mesor (baseline expression)
    M_mean = mean(M_vals, na.rm = TRUE),
    M_sd = sd(M_vals, na.rm = TRUE),

    # Amplitude (among rhythmic genes)
    A_mean = mean(A_vals[rhythmic_genes], na.rm = TRUE),
    A_sd = sd(A_vals[rhythmic_genes], na.rm = TRUE),
    A_median = median(A_vals[rhythmic_genes], na.rm = TRUE),
    A_q25 = quantile(A_vals[rhythmic_genes], 0.25, na.rm = TRUE),
    A_q75 = quantile(A_vals[rhythmic_genes], 0.75, na.rm = TRUE),

    # Phase (among rhythmic genes)
    phi_mean = circular_mean(phi_vals[rhythmic_genes], period),
    phi_concentration = circular_concentration(phi_vals[rhythmic_genes], period),
    phi_uniform = FALSE,  # Will be set to TRUE if distribution is uniform

    # Noise level (all genes)
    sigma_mean = mean(sigma_vals, na.rm = TRUE),
    sigma_sd = sd(sigma_vals, na.rm = TRUE),
    sigma_median = median(sigma_vals, na.rm = TRUE),

    # Effect size (among rhythmic genes)
    r_mean = mean(r_vals[rhythmic_genes], na.rm = TRUE),
    r_sd = sd(r_vals[rhythmic_genes], na.rm = TRUE),
    r_median = median(r_vals[rhythmic_genes], na.rm = TRUE),
    r_q25 = quantile(r_vals[rhythmic_genes], 0.25, na.rm = TRUE),
    r_q75 = quantile(r_vals[rhythmic_genes], 0.75, na.rm = TRUE),

    # Proportion rhythmic (empirical estimate)
    prop_rhythmic = mean(rhythmic_genes, na.rm = TRUE)
  )

  # Test for uniform phase distribution
  if (length(na.omit(phi_vals[rhythmic_genes])) > 10) {
    rayleigh_test = rayleigh_test_circular(phi_vals[rhythmic_genes], period)
    params$phi_uniform = rayleigh_test$pvalue > 0.05
    params$phi_rayleigh_p = rayleigh_test$pvalue
  }

  # Store raw values for flexible simulation
  params$raw = list(
    M = M_vals,
    A = A_vals,
    phi = phi_vals,
    sigma = sigma_vals,
    r = r_vals,
    pvalue = pvals,
    is_rhythmic = rhythmic_genes
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

  # P-value (approximation)
  pvalue = exp(-Z)

  return(list(
    statistic = Z,
    pvalue = pvalue
  ))
}


#' Estimate parameters for two-group comparison
#'
#' @description Estimate parameters separately for two conditions
#'
#' @param data_A Data matrix for condition A
#' @param data_B Data matrix for condition B
#' @param times_A Time points for condition A
#' @param times_B Time points for condition B
#' @param ... Other arguments
#'
#' @return List with parameters for both conditions and differences
estimate_two_group_params = function(data_A, data_B, times_A, times_B, ...) {

  params_A = estimate_circadian_params(data_A, times_A, ...)
  params_B = estimate_circadian_params(data_B, times_B, ...)

  # Estimate typical differences
  # Only consider genes rhythmic in at least one condition
  rhythmic_either = params_A$raw$is_rhythmic | params_B$raw$is_rhythmic

  delta_A = params_B$raw$A[rhythmic_either] - params_A$raw$A[rhythmic_either]
  delta_phi = circular_difference(
    params_B$raw$phi[rhythmic_either],
    params_A$raw$phi[rhythmic_either],
    24
  )

  diff_params = list(
    delta_A_mean = mean(abs(delta_A), na.rm = TRUE),
    delta_A_sd = sd(delta_A, na.rm = TRUE),
    delta_phi_mean = circular_mean(abs(delta_phi), 24),
    delta_phi_sd = circular_sd(delta_phi, 24),
    prop_gain = mean(delta_A > 0.1, na.rm = TRUE),  # Amplitude increase
    prop_loss = mean(delta_A < -0.1, na.rm = TRUE), # Amplitude decrease
    prop_phase_shift = mean(abs(delta_phi) > 2, na.rm = TRUE)  # > 2h shift
  )

  return(list(
    params_A = params_A,
    params_B = params_B,
    diff_params = diff_params
  ))
}


#' Circular difference (accounting for wraparound)
#'
#' @param phi1 First angle
#' @param phi2 Second angle
#' @param period Period
#'
#' @return Difference in [-period/2, period/2]
circular_difference = function(phi1, phi2, period = 24) {
  diff = phi1 - phi2
  ((diff + period/2) %% period) - period/2
}


#' Circular standard deviation
#'
#' @param angles Angles (in hours)
#' @param period Period (default 24)
#'
#' @return Circular SD (in hours)
circular_sd = function(angles, period = 24) {
  angles = angles[!is.na(angles)]
  if (length(angles) == 0) return(NA)

  R = circular_concentration(angles, period)

  # Circular variance
  var_circ = 1 - R

  # Circular SD (approximate)
  sd_circ = sqrt(2 * var_circ) * period / (2 * pi)

  return(sd_circ)
}


#' Sample from estimated parameter distributions
#'
#' @description Generate random parameter values for simulation
#'
#' @param params Parameter list from estimate_circadian_params
#' @param n Number of samples
#'
#' @return Data frame with sampled parameters
sample_params = function(params, n = 1000) {

  # Mesor
  M = rnorm(n, params$M_mean, params$M_sd)

  # Amplitude (truncated at 0)
  A = pmax(rnorm(n, params$A_mean, params$A_sd), 0)

  # Phase
  if (params$phi_uniform) {
    phi = runif(n, 0, 24)
  } else {
    # Von Mises-like sampling (approximate with wrapped normal)
    phi = (params$phi_mean + rnorm(n, 0, 1 / params$phi_concentration)) %% 24
  }

  # Noise
  sigma = pmax(rnorm(n, params$sigma_mean, params$sigma_sd), 0.05)

  # Effect size
  r = A / sigma

  return(data.frame(M = M, A = A, phi = phi, sigma = sigma, r = r))
}


#' Estimate CircadianBioOptions from Pilot Data
#'
#' @description Bridge function: estimates circadian parameters from pilot
#' expression data and returns a CircadianBioOptions object with empirical
#' distributions (numeric vectors resampled by set*() helpers).
#'
#' @param data Gene expression matrix (genes x samples) for one group
#' @param times Time points for each sample
#' @param period Circadian period (default 24)
#' @param min_rhythm_pval P-value threshold for rhythmic gene inclusion
#' @param prop_DR Differential rhythmicity proportion (user-specified)
#' @param prop_DP Differential phase proportion (user-specified)
#' @param prop_DA Differential amplitude proportion (user-specified)
#' @param phase_diff Phase shift range for DP genes
#' @param amp_diff Amplitude ratio range for DA genes
#' @param dp_shift_mode "fixed" (use phase_diff[2]) or "uniform" (sample within phase_diff range)
#' @param dr_amp_scale Scale factor for amplitude (A) to adjust DR strength
#' @param dr_sigma_scale Scale factor for sigma to adjust DR strength
#' @param sim.seed Random seed
#' @param verbose Print progress
#'
#' @return CircadianBioOptions with empirical parameter distributions
estCircadianParam <- function(data, times, period = 24,
                              min_rhythm_pval = 0.1,
                              prop_DR = 0.15, prop_DP = 0.10, prop_DA = 0.10,
                              prop_DM = 0.00, mesor_diff = c(0.5, 2.0),
                              phase_diff = c(-6, 6), amp_diff = c(0.5, 2),
                              dp_shift_mode = c("fixed", "uniform"),
                              dr_amp_scale = 1.0,
                              dr_sigma_scale = 1.0,
                              sim.seed = 12345, verbose = TRUE) {

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

  # Amplitude and sigma: use rhythmic genes only, keeping them paired (same gene index)
  # to preserve the empirical A-sigma correlation in downstream joint sampling.
  rhythmic_valid <- rhythmic_idx & !is.na(params$raw$A) & params$raw$A > 0 &
                    !is.na(params$raw$sigma) & params$raw$sigma > 0
  amp_emp         <- params$raw$A[rhythmic_valid]
  sigma_rhythmic_emp <- params$raw$sigma[rhythmic_valid]   # paired with amp_emp

  # Phase: use rhythmic genes
  phase_emp <- params$raw$phi[rhythmic_idx & !is.na(params$raw$phi)]

  # Cap differential proportions at the estimated rhythmic budget.
  # CircadianBioOptions requires prop_DR + prop_DP + prop_DA + prop_DM <= prop_rhythmic
  # because every differential gene must be rhythmic in at least one group.
  total_diff <- prop_DR + prop_DP + prop_DA + prop_DM
  if (total_diff > params$prop_rhythmic && total_diff > 0) {
    scale_factor <- params$prop_rhythmic / total_diff
    if (verbose) {
      message(sprintf(
        paste0("estCircadianParam: prop_DR+prop_DP+prop_DA+prop_DM (%.3f) exceeds estimated ",
               "prop_rhythmic (%.3f). Scaling differential props by %.3f to fit budget."),
        total_diff, params$prop_rhythmic, scale_factor))
    }
    prop_DR <- prop_DR * scale_factor
    prop_DP <- prop_DP * scale_factor
    prop_DA <- prop_DA * scale_factor
    prop_DM <- prop_DM * scale_factor
  }

  CircadianBioOptions(
    ngenes = ngenes,
    prop_rhythmic = params$prop_rhythmic,
    period = period,
    lBaselineExpr = lBaselineExpr_emp,
    lOD = lOD_emp,
    amplitude = amp_emp,
    sigma_rhythmic = sigma_rhythmic_emp,
    phase = phase_emp,
    prop_DR = prop_DR,
    prop_DP = prop_DP,
    prop_DA = prop_DA,
    prop_DM = prop_DM,
    mesor_diff = mesor_diff,
    phase_diff = phase_diff,
    amp_diff = amp_diff,
    dp_shift_mode = dp_shift_mode,
    dr_amp_scale = dr_amp_scale,
    dr_sigma_scale = dr_sigma_scale,
    sim.seed = sim.seed
  )
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
#' @param min_rhythm_pval P-value threshold for rhythmic classification (default 0.1)
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
estCircadianParamTwoGroup <- function(data_1, data_2, times_1, times_2,
                                      period = 24,
                                      min_rhythm_pval = 0.1,
                                      phase_shift_threshold = 2,
                                      prop_DM = 0.00,
                                      mesor_diff = c(0.5, 2.0),
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
  prop_DP_emp    <- prop_joint * mean(dp_among_joint, na.rm = TRUE)

  # Budget guard: cap prop_DR + prop_DP + prop_DM to fit within union rhythmic budget.
  # Floating-point arithmetic can cause the sum to equal or fractionally exceed
  # prop_union_rhy, which would error in CircadianBioOptions.
  total_diff <- prop_DR_emp + prop_DP_emp + prop_DM
  if (total_diff >= prop_union_rhy && total_diff > 0) {
    scale_factor <- prop_union_rhy * 0.999 / total_diff
    prop_DR_emp  <- prop_DR_emp * scale_factor
    prop_DP_emp  <- prop_DP_emp * scale_factor
    prop_DM      <- prop_DM    * scale_factor
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

  # Warn if too few DP genes to trust the phase_diff IQR ---------------------
  n_dp_pilot <- sum(valid_dp, na.rm = TRUE)
  if (n_dp_pilot < 20 && verbose) {
    warning(sprintf(
      paste0("Only %d jointly-rhythmic DP genes found in pilot (|Δφ|>%.1fh).\n",
             "  phase_diff IQR estimate may be unreliable.\n",
             "  Consider using phase_diff = c(-6, 6) as a conservative default."),
      n_dp_pilot, phase_shift_threshold))
  }

  # amp_diff: kept for CircadianBioOptions signature (prop_DA = 0, so unused) -
  amp_diff_emp <- c(0.5, 2.0)

  # Baseline distributions from group 1 (consistent with estCircadianParam) ---
  lBaselineExpr_emp <- p1$raw$M[!is.na(p1$raw$M)]
  sigma_valid       <- p1$raw$sigma[!is.na(p1$raw$sigma) & p1$raw$sigma > 0]
  lOD_emp           <- log(sigma_valid)
  rhythmic_idx      <- p1$raw$is_rhythmic
  # Keep A and sigma paired from the same gene for joint sampling
  rhythmic_valid    <- rhythmic_idx & !is.na(p1$raw$A) & p1$raw$A > 0 &
                       !is.na(p1$raw$sigma) & p1$raw$sigma > 0
  amp_emp            <- p1$raw$A[rhythmic_valid]
  sigma_rhythmic_emp <- p1$raw$sigma[rhythmic_valid]
  phase_emp          <- p1$raw$phi[rhythmic_idx & !is.na(p1$raw$phi)]

  # Group-2 distributions for fully symmetric two-group simulation -----------
  rhythmic_idx_2     <- p2$raw$is_rhythmic
  amp_emp2           <- p2$raw$A[rhythmic_idx_2 & !is.na(p2$raw$A) & p2$raw$A > 0]
  sigma_valid_2      <- p2$raw$sigma[!is.na(p2$raw$sigma) & p2$raw$sigma > 0]
  lOD_emp2           <- log(sigma_valid_2)
  lBaselineExpr2_emp <- p2$raw$M[!is.na(p2$raw$M)]   # group-2 mesor distribution

  # Diagnostics ---------------------------------------------------------------
  diagnostics <- list(
    # Simulation hyperparameters
    prop_DR_emp           = prop_DR_emp,
    prop_DP_emp           = prop_DP_emp,
    phase_diff_emp        = phase_diff_emp,
    phase_shift_threshold = phase_shift_threshold,
    n_DP_genes_pilot      = n_dp_pilot,
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
    cat("\n  prop_DR, prop_DP, phase_diff are used directly in CircadianBioOptions.\n")
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
    amplitude2      = amp_emp2,            # F̂_A2: used for g2-only DR genes
    cts2            = times_2,             # F̂_TOD2: group-2 sampling time distribution
    phase           = phase_emp,
    prop_DR         = prop_DR_emp,
    prop_DP         = prop_DP_emp,
    prop_DA         = 0,
    prop_DM         = prop_DM,             # user-specified; not estimated from pilot
    mesor_diff      = mesor_diff,
    phase_diff      = phase_diff_emp,
    amp_diff        = amp_diff_emp,
    dp_shift_mode   = "uniform",
    sim.seed        = sim.seed
  )

  opts$diagnostics <- diagnostics
  opts
}
