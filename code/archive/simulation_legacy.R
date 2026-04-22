simCircadian <- function(simOptions,
                         times = "even",
                         n = 24,
                         design = c("active", "passive"),
                         cts = NULL,
                         phi = 0,
                         noise_type = "gaussian") {

  design = match.arg(design)

  # Handle time points
  if (design == "active") {
    if (is.character(times) && times == "even") {
      times = seq(0, simOptions$period, length.out = n + 1L)[seq_len(n)]
    }
    times_expanded = times
  } else {
    # Passive design: sample from TOD distribution
    if (is.null(cts)) {
      stop("For passive design, must provide cts (TOD distribution)")
    }
    times_expanded = sampleTimesFromDist(n, cts)
  }

  # Generate data
  ngenes = simOptions$ngenes
  n_rhythmic = round(ngenes * simOptions$prop_rhythmic)

  # Select which genes are rhythmic
  rhythmic_id = sample(1:ngenes, n_rhythmic)
  is_rhythmic = rep(FALSE, ngenes)
  is_rhythmic[rhythmic_id] = TRUE

  # Assign parameters
  mesor = simOptions$lBaselineExpr
  sigma = exp(simOptions$lOD)

  # Amplitude (0 for non-rhythmic genes)
  amplitude = rep(0, ngenes)
  amplitude[rhythmic_id] = simOptions$amplitude

  # Phase
  phase = rep(0, ngenes)
  phase[rhythmic_id] = simOptions$phase

  # Apply phase shift for passive design
  if (design == "passive" && is.numeric(phi)) {
    phase = (phase + phi) %% simOptions$period
  }

  # Generate expression matrix
  expr = matrix(NA, nrow = ngenes, ncol = n)
  omega = 2 * pi / simOptions$period

  for (g in 1:ngenes) {
    mu = mesor[g] + amplitude[g] * cos(omega * times_expanded - omega * phase[g])

    if (noise_type == "gaussian") {
      expr[g, ] = rnorm(n, mu, sigma[g])
    } else if (noise_type == "t") {
      expr[g, ] = mu + sigma[g] * rt(n, df = 5)
    } else if (noise_type == "lognormal") {
      expr[g, ] = exp(rnorm(n, mu, sigma[g]))
    }
  }

  rownames(expr) = paste0("Gene", 1:ngenes)
  colnames(expr) = paste0("Sample", 1:n)

  # Ground truth
  effectsize = amplitude / sigma  # r = A/σ

  ground_truth = data.frame(
    gene = 1:ngenes,
    is_rhythmic = is_rhythmic,
    mesor = mesor,
    amplitude = amplitude,
    phase = phase,
    sigma = sigma,
    effectsize = effectsize,
    stringsAsFactors = FALSE
  )

  return(list(
    expr = expr,
    times = times_expanded,
    ground_truth = ground_truth,
    rhythmic_id = rhythmic_id,
    effectsize = effectsize,
    simOptions = simOptions,
    design = design
  ))
}

#' Sample time points from TOD distribution
#'
#' @param n Number of time points to sample
#' @param cts Observed TOD distribution (kernel density estimated internally)


simulate_circadian_data = function(G = 1000,
                                   G_rhythmic = NULL,
                                   prop_rhythmic = 0.1,
                                   n = 3,
                                   times = seq(0, 48, by = 4),
                                   period = 24,
                                   params = NULL,
                                   noise_type = "gaussian",
                                   waveform = "sinusoid",
                                   seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # Total samples
  N = length(times) * n

  # Determine number of rhythmic genes
  if (is.null(G_rhythmic)) {
    G_rhythmic = floor(G * prop_rhythmic)
  }
  G_non = G - G_rhythmic

  # Default parameters (can be overridden)
  if (is.null(params)) {
    params = list(
      M_mean = 8,           # Mean expression (log scale)
      M_sd = 1,             # Variation in baseline
      A_mean = 0.5,         # Mean amplitude
      A_sd = 0.3,           # Amplitude variation
      phi_mean = NULL,      # NULL = uniform on [0, period)
      sigma_mean = 0.3,     # Mean noise level
      sigma_sd = 0.1        # Noise variation
    )
  }

  # Generate parameters for each gene
  M = rnorm(G, params$M_mean, params$M_sd)
  sigma = pmax(rnorm(G, params$sigma_mean, params$sigma_sd), 0.05)  # Min sigma

  # Non-rhythmic genes: A = 0
  A = rep(0, G)
  phi = rep(0, G)

  # Rhythmic genes: draw A and phi
  if (G_rhythmic > 0) {
    rhythmic_idx = sample(G, G_rhythmic)
    A[rhythmic_idx] = pmax(rnorm(G_rhythmic, params$A_mean, params$A_sd), 0.01)

    # Check if we should use uniform phase distribution
    use_uniform_phase = is.null(params$phi_mean) ||
                        (!is.null(params$phi_uniform) && params$phi_uniform)

    if (use_uniform_phase) {
      # Uniform phase distribution
      phi[rhythmic_idx] = runif(G_rhythmic, 0, period)
    } else {
      # Concentrated around mean
      # Use von Mises-like approach
      phi_sd = if (is.null(params$phi_sd)) 4 else params$phi_sd  # Default SD = 4h
      phi[rhythmic_idx] = params$phi_mean + rnorm(G_rhythmic, 0, phi_sd)
      phi[rhythmic_idx] = phi[rhythmic_idx] %% period
    }
  }

  # Generate data matrix
  data = matrix(NA, nrow = G, ncol = N)
  colnames(data) = paste0("Sample", 1:N)
  rownames(data) = paste0("Gene", 1:G)

  # Expand times for replicates
  times_expanded = rep(times, each = n)

  # Angular frequency
  omega = 2 * pi / period

  # Generate expression for each gene
  for (g in 1:G) {
    # Mean function (rhythmic component)
    if (waveform == "sinusoid") {
      mu = M[g] + A[g] * cos(omega * times_expanded - omega * phi[g])
    } else if (waveform == "damped") {
      # Amplitude decays over time
      damping = exp(-0.01 * times_expanded)
      mu = M[g] + A[g] * damping * cos(omega * times_expanded - omega * phi[g])
    } else if (waveform == "asymmetric") {
      # Add second harmonic for asymmetry
      mu = M[g] + A[g] * (cos(omega * times_expanded - omega * phi[g]) +
                          0.3 * cos(2 * omega * times_expanded - 2 * omega * phi[g]))
    }

    # Add noise
    if (noise_type == "gaussian") {
      y = rnorm(N, mu, sigma[g])
    } else if (noise_type == "t") {
      # t-distribution (df=5 for heavy tails)
      y = mu + sigma[g] * rt(N, df = 5)
    } else if (noise_type == "lognormal") {
      # Log-normal noise
      y = exp(rnorm(N, mu, sigma[g]))
    } else if (noise_type == "negbinom") {
      # Negative binomial (for counts)
      # Convert to count scale
      mu_count = exp(mu)
      size = 10  # Dispersion parameter
      y = rnbinom(N, size = size, mu = mu_count)
    }

    data[g, ] = y
  }

  # Ground truth
  ground_truth = data.frame(
    gene = 1:G,
    is_rhythmic = A > 0,
    M = M,
    A = A,
    phi = phi,
    sigma = sigma,
    r = A / sigma,  # Effect size
    stringsAsFactors = FALSE
  )

  return(list(
    data = data,
    times = times_expanded,
    times_unique = times,
    n_per_time = n,
    ground_truth = ground_truth,
    params = params,
    noise_type = noise_type,
    waveform = waveform,
    period = period
  ))
}


#' Simulate Two-Group Circadian Comparison
#'
#' @description Generate data for comparing circadian rhythms between two conditions
#'
#' @param G Number of genes
#' @param n_A Samples per time point in condition A
#' @param n_B Samples per time point in condition B
#' @param times Time points (same for both conditions)
#' @param diff_type "amplitude", "phase", "rhythmicity", or "mixed"
#' @param delta_A Amplitude difference (for amplitude/mixed type)
#' @param delta_phi Phase difference in hours (for phase/mixed type)
#' @param prop_diff Proportion of genes with difference
#' @param ... Other arguments passed to simulate_circadian_data
#'
#' @return List with data for both conditions and ground truth
simulate_two_group = function(G = 1000,
                              n_A = 3,
                              n_B = 3,
                              times = seq(0, 48, by = 4),
                              diff_type = "amplitude",
                              delta_A = 0.3,
                              delta_phi = 4,
                              prop_diff = 0.1,
                              period = 24,
                              params = NULL,
                              noise_type = "gaussian",
                              seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # Default parameters
  if (is.null(params)) {
    params = list(
      M_mean = 8,
      M_sd = 1,
      A_mean = 0.5,
      A_sd = 0.3,
      phi_mean = NULL,
      sigma_mean = 0.3,
      sigma_sd = 0.1
    )
  }

  N_A = length(times) * n_A
  N_B = length(times) * n_B
  G_diff = floor(G * prop_diff)

  # Generate base parameters
  M = rnorm(G, params$M_mean, params$M_sd)
  sigma_A = pmax(rnorm(G, params$sigma_mean, params$sigma_sd), 0.05)
  sigma_B = sigma_A  # Same noise in both conditions

  # Amplitude and phase
  A_A = pmax(rnorm(G, params$A_mean, params$A_sd), 0)

  # Check if we should use uniform phase distribution
  use_uniform_phase = is.null(params$phi_mean) ||
                      (!is.null(params$phi_uniform) && params$phi_uniform)

  if (use_uniform_phase) {
    phi_A = runif(G, 0, period)
  } else {
    phi_sd = if (is.null(params$phi_sd)) 4 else params$phi_sd  # Default SD = 4h
    phi_A = params$phi_mean + rnorm(G, 0, phi_sd)
    phi_A = phi_A %% period
  }

  # Initialize B parameters (same as A initially)
  A_B = A_A
  phi_B = phi_A

  # Introduce differences
  diff_idx = sample(G, G_diff)
  diff_type_gene = character(G)

  if (diff_type == "amplitude") {
    A_B[diff_idx] = pmax(A_A[diff_idx] + delta_A, 0)
    diff_type_gene[diff_idx] = "amplitude"
  } else if (diff_type == "phase") {
    phi_B[diff_idx] = (phi_A[diff_idx] + delta_phi) %% period
    diff_type_gene[diff_idx] = "phase"
  } else if (diff_type == "rhythmicity") {
    # Rhythmic in A, arrhythmic in B
    A_B[diff_idx] = 0
    diff_type_gene[diff_idx] = "rhythmicity"
  } else if (diff_type == "mixed") {
    # Randomly assign different types
    n_amp = floor(G_diff / 3)
    n_phase = floor(G_diff / 3)
    n_rhyth = G_diff - n_amp - n_phase

    amp_idx = diff_idx[1:n_amp]
    phase_idx = diff_idx[(n_amp + 1):(n_amp + n_phase)]
    rhyth_idx = diff_idx[(n_amp + n_phase + 1):G_diff]

    A_B[amp_idx] = pmax(A_A[amp_idx] + delta_A, 0)
    phi_B[phase_idx] = (phi_A[phase_idx] + delta_phi) %% period
    A_B[rhyth_idx] = 0

    diff_type_gene[amp_idx] = "amplitude"
    diff_type_gene[phase_idx] = "phase"
    diff_type_gene[rhyth_idx] = "rhythmicity"
  }

  diff_type_gene[setdiff(1:G, diff_idx)] = "none"

  # Generate data for each condition
  omega = 2 * pi / period
  times_A = rep(times, each = n_A)
  times_B = rep(times, each = n_B)

  data_A = matrix(NA, nrow = G, ncol = N_A)
  data_B = matrix(NA, nrow = G, ncol = N_B)

  for (g in 1:G) {
    # Condition A
    mu_A = M[g] + A_A[g] * cos(omega * times_A - omega * phi_A[g])
    data_A[g, ] = rnorm(N_A, mu_A, sigma_A[g])

    # Condition B
    mu_B = M[g] + A_B[g] * cos(omega * times_B - omega * phi_B[g])
    data_B[g, ] = rnorm(N_B, mu_B, sigma_B[g])
  }

  colnames(data_A) = paste0("A_Sample", 1:N_A)
  colnames(data_B) = paste0("B_Sample", 1:N_B)
  rownames(data_A) = paste0("Gene", 1:G)
  rownames(data_B) = paste0("Gene", 1:G)

  # Ground truth
  ground_truth = data.frame(
    gene = 1:G,
    is_rhythmic_A = A_A > 0,
    is_rhythmic_B = A_B > 0,
    diff_type = diff_type_gene,
    has_diff = diff_type_gene != "none",
    M = M,
    A_A = A_A,
    A_B = A_B,
    delta_A = A_B - A_A,
    phi_A = phi_A,
    phi_B = phi_B,
    delta_phi = ((phi_B - phi_A + period/2) %% period) - period/2,  # Circular difference
    sigma_A = sigma_A,
    sigma_B = sigma_B,
    r_A = A_A / sigma_A,
    r_B = A_B / sigma_B,
    stringsAsFactors = FALSE
  )

  return(list(
    data_A = data_A,
    data_B = data_B,
    times_A = times_A,
    times_B = times_B,
    times_unique = times,
    n_per_time_A = n_A,
    n_per_time_B = n_B,
    ground_truth = ground_truth,
    diff_type = diff_type,
    params = params,
    period = period
  ))
}


# =====================================================================
# Single-cohort data generator (method-agnostic, harmonic support)
# =====================================================================

#' Simulate single-cohort circadian expression data
#'
#' Internal data generator used by runSingleCohortPower() and
#' runBootstrapDesignGrid(mode="single"). Supports harmonic waveform
#' deviation (alpha2, alpha3) so the same function feeds both standard
#' power analysis and waveform-violation experiments.
#'
#' @param bio.opts   CircadianBioOptions from estCircadianParam()
#' @param cts        Numeric vector of time points, length = N (already expanded)
#' @param alpha2     2nd harmonic relative amplitude (0 = pure cosinor)
#' @param alpha3     3rd harmonic relative amplitude
#' @param seed       Optional random seed
#'
#' @return list(expr = matrix[G x N], is_rhythmic = logical[G], r_values = numeric[G])
