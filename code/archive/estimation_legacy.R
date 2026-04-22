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
    delta_phi_mean = mean(abs(delta_phi), na.rm = TRUE),
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
#' @param phase_diff Phase shift range for DP genes
#' @param amp_diff Unused; retained for interface compatibility
#' @param dp_shift_mode "fixed" (use phase_diff[2]) or "uniform" (sample within phase_diff range)
#' @param dr_amp_scale Scale factor for amplitude (A) to adjust DR strength
#' @param dr_sigma_scale Scale factor for sigma to adjust DR strength
#' @param sim.seed Random seed
#' @param verbose Print progress
#'
#' @return CircadianBioOptions with empirical parameter distributions
