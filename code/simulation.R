#' Simulate Circadian Time-Series Data
#'
#' Generate simulated circadian expression data for power analysis
#'
#' @param simOptions Simulation options from createSimOptions()
#' @param times Time points (vector or "even" for evenly-spaced)
#' @param n Total number of samples
#' @param design "active" or "passive"
#' @param cts For passive design: TOD distribution to sample from
#' @param phi For passive design: phase shift (or "vary" to test all)
#' @param noise_type "gaussian", "t", "lognormal"
#'
#' @return List with simulated data and ground truth

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
      times = seq(0, simOptions$period, length.out = n)
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

sampleTimesFromDist <- function(n, cts) {
  # Map all times to [0, 24) first.  Pilot TOD values can be negative (e.g.,
  # -5.6h = 18.4h in circular time); the original density(from=0,to=24) call
  # silently discarded those, causing the passive simulation to miss evening
  # coverage that the real pilot has.
  cts_circ <- cts %% 24

  # Circular (wrap-around) KDE: replicate the data at -24 and +24 so that
  # boundary contributions near 0h and 24h are handled correctly, then
  # evaluate only on [0, 24].
  cts_ext    <- c(cts_circ - 24, cts_circ, cts_circ + 24)
  tod_density <- density(cts_ext, from = 0, to = 24)

  # Sample n time points from the wrapped density
  sampled <- sample(tod_density$x, size = n, replace = TRUE, prob = tod_density$y)

  # Add small noise to avoid exact ties, then wrap back to [0, 24)
  sampled <- (sampled + rnorm(n, 0, 0.1)) %% 24

  return(sampled)
}
#' Simulate Two-Group Differential Circadian Data
#'
#' Simulate data for comparing circadian rhythms between two groups
#' (e.g., control vs treatment, disease vs healthy).
#'
#' @section Gene expression model:
#' Each gene g follows the cosinor model:
#'   y_g(t) = M_g + A_g * cos(omega * t - omega * phi_g) + epsilon_g
#' where omega = 2*pi/period and epsilon_g ~ N(0, sigma_g^2).
#'
#' @section Default parameter distributions:
#' When lBaselineExpr and lOD are not provided, parameters are drawn from:
#'
#' **Noise (sigma):**
#'   log(sigma_g) ~ N(mu = -1, sd = 0.3)
#'   => sigma_g ~ LogNormal(mu = -1, sd = 0.3)
#'   => median(sigma) = exp(-1) ~ 0.37, IQR ~ [0.27, 0.49]
#'
#' **Amplitude (A):**
#'   A_g ~ max(LogNormal(mu = log(0.4), sd = 0.5), 0.05)
#'   => median(A) = 0.4, IQR ~ [0.24, 0.66]
#'   (Clamped below at 0.05 to prevent near-zero amplitudes)
#'
#' **Signal-to-noise ratio (r = A/sigma):**
#'   Since both A and sigma are log-normal, r = A/sigma is also approximately
#'   log-normal with:
#'     log(r) ~ N(log(0.4) - (-1), sqrt(0.5^2 + 0.3^2))
#'            = N(log(0.4) + 1, 0.58)
#'            ~ N(0.08, 0.58)
#'   => median(r) ~ exp(0.08) ~ 1.08
#'   => Most genes fall in r = 0.3 - 2.0 (right-skewed distribution)
#'
#'   This produces a realistic scenario where the bulk of target genes have
#'   moderate signal strength (r = 0.5 - 1.5), few genes have strong signals
#'   (r > 3), and marginal power is dominated by the crowded low-to-moderate
#'   r strata. See Panel F in the stratified power plots.
#'
#' @section Types of differential expression:
#' \describe{
#'   \item{Type 0}{Non-rhythmic in both groups}
#'   \item{Type 1}{Rhythmic in both, same parameters (control)}
#'   \item{Type 2}{DR: Rhythmic in Group 1 only}
#'   \item{Type 3}{DR: Rhythmic in Group 2 only}
#'   \item{Type 4}{DP: Differential phase (both rhythmic, different peak times)}
#'   \item{Type 5}{DA: Differential amplitude (both rhythmic, different amplitudes)}
#' }
#'
#' @param ngenes Number of genes
#' @param n1 Sample size for group 1
#' @param n2 Sample size for group 2
#' @param lBaselineExpr Log baseline expression. Default: rnorm(ngenes, 5, 2)
#' @param lOD Log over-dispersion (noise). Default: rnorm(ngenes, -1, 0.3),
#'   so sigma = exp(lOD) ~ LogNormal(mu=-1, sd=0.3) with median ~0.37
#' @param prop_rhythmic Overall proportion of rhythmic genes
#' @param prop_DR Proportion with differential rhythmicity
#' @param prop_DP Proportion with differential phase (among rhythmic in both)
#' @param prop_DA Proportion with differential amplitude (among rhythmic in both)
#' @param phase_diff Range of phase shift (hours) for DP genes
#' @param amp_diff Range of amplitude ratio for DA genes (e.g., c(0.5, 2))
#' @param period Period (default 24)
#' @param design "active" or "passive"
#' @param cts TOD distribution for passive design
#' @param sim.seed Random seed
#'
#' @return List with simulated data for both groups and ground truth

simCircadianDiff <- function(ngenes = 5000,
                           n1 = 24,
                           n2 = 24,
                           lBaselineExpr,
                           lOD,
                           lOD2 = NULL,
                           amplitude,
                           amplitude2 = NULL,
                           sigma_rhythmic = NULL,
                           prop_rhythmic = 0.25,
                           prop_DR = 0.1,
                           prop_DP = 0.2,
                           prop_DA = 0.2,
                           prop_DM = 0,
                           phase_diff = c(-6, 6),
                           dp_shift_mode = c("fixed", "uniform"),
                           amp_diff = c(0.5, 2),
                           mesor_diff = c(-1, 1),
                           period = 24,
                           design = c("active", "passive"),
                           cts = NULL,
                           cts2 = NULL,
                           sim.seed = 42,
                           harmonics = NULL) {

  design = match.arg(design)
  dp_shift_mode <- match.arg(dp_shift_mode)

  # Save group-2 distributions BEFORE local variables shadow these parameters.
  amp_dist2 <- amplitude2
  lod_dist2 <- lOD2

  # ---------------------------------------------------------------
  # IMPORTANT: Gene assignment (diff_type) must use a seed that does
  # NOT depend on sample size (n1, n2). Time sampling consumes random
  # numbers proportional to n, so if done first it shifts the RNG
  # state and makes diff_type change across sample sizes.
  # Fix: do gene assignment with sim.seed FIRST, then reseed for
  # n-dependent sampling with a derived seed.
  # ---------------------------------------------------------------
  set.seed(sim.seed)

  # Set baseline expression and noise
  # M_g ~ N(5, 2) on log scale
  if (missing(lBaselineExpr)) {
    lBaselineExpr = rnorm(ngenes, mean = 5, sd = 2)
  } else if (length(lBaselineExpr) != ngenes) {
    # Resample empirical vector to match ngenes
    lBaselineExpr = sample(lBaselineExpr, ngenes, replace = TRUE)
  }
  # log(sigma_g) ~ N(-1, 0.3)  =>  sigma_g ~ LogNormal(mu=-1, sd=0.3)
  # median(sigma) = exp(-1) ~ 0.37,  IQR ~ [0.27, 0.49]
  if (missing(lOD)) {
    lOD = rnorm(ngenes, mean = -1, sd = 0.3)
  } else if (length(lOD) != ngenes) {
    lOD = sample(lOD, ngenes, replace = TRUE)
  }

  # Initialize parameters for both groups
  mesor = lBaselineExpr
  sigma = exp(lOD)   # group-1 noise
  # Group-2 noise: use lOD2 if provided (two-group pilot), else share with group 1
  sigma2_vec <- if (!is.null(lod_dist2)) {
    lod2_res <- if (length(lod_dist2) == ngenes) lod_dist2
                else sample(lod_dist2, ngenes, replace = TRUE)
    exp(lod2_res)
  } else {
    sigma
  }
  omega = 2 * pi / period

  # Initialize classification BEFORE time sampling (n-independent)
  # Categories: 0=non-rhythmic both, 1=rhythmic both (same),
  # 2=DR (only g1 rhythmic), 3=DR (only g2 rhythmic),
  # 4=DP (diff phase), 5=DA (diff amplitude)
  diff_type = rep(0, ngenes)

  # Count differentially rhythmic genes
  n_DR = round(ngenes * prop_DR)
  n_DP = round(ngenes * prop_DP)
  n_DA = round(ngenes * prop_DA)

  # Assign categories
  gene_idx = 1:ngenes

  # Differential Rhythmicity (DR): rhythmic in g1 only or g2 only
  if (n_DR > 0) {
    DR_genes = sample(gene_idx, n_DR)
    DR_g1_only = DR_genes[1:ceiling(n_DR/2)]
    DR_g2_only = DR_genes[(ceiling(n_DR/2)+1):n_DR]
    diff_type[DR_g1_only] = 2
    diff_type[DR_g2_only] = 3
  } else {
    DR_genes = integer(0)
    DR_g1_only = DR_g2_only = integer(0)
  }

  # Rhythmic in both (for DP, DA, or same)
  # prop_rhythmic = total fraction rhythmic in at least one group.
  # DR genes already count toward that budget; rhythmic_both fills the remainder.
  remaining = setdiff(gene_idx, DR_genes)
  n_rhythmic_total <- round(ngenes * prop_rhythmic)
  n_rhythmic_both  <- max(n_DP + n_DA, n_rhythmic_total - n_DR)  # must fit DP + DA
  n_rhythmic_both  <- min(n_rhythmic_both, length(remaining))
  rhythmic_both    <- if (n_rhythmic_both > 0) sample(remaining, n_rhythmic_both) else integer(0)

  # Differential Phase (DP)
  if (n_DP > 0 && length(rhythmic_both) > 0) {
    DP_genes = sample(setdiff(rhythmic_both, DR_genes), min(n_DP, length(rhythmic_both)))
    diff_type[DP_genes] = 4
  } else {
    DP_genes = integer(0)
  }

  # Differential Amplitude (DA)
  if (n_DA > 0 && length(rhythmic_both) > 0) {
    DA_pool = setdiff(rhythmic_both, c(DR_genes, DP_genes))
    n_DA_actual = min(n_DA, length(DA_pool))
    DA_genes = sample(DA_pool, n_DA_actual)
    diff_type[DA_genes] = 5
  } else {
    DA_genes = integer(0)
  }

  # Rhythmic in both, same (control genes)
  rhythmic_same = setdiff(rhythmic_both, c(DP_genes, DA_genes))
  if (length(rhythmic_same) > 0) {
    diff_type[rhythmic_same] = 1
  }

  # Generate base parameters (group 1)
  phase1 = rep(0, ngenes)
  amplitude1 = rep(0, ngenes)

  # All rhythmic genes get base phase and amplitude
  rhythmic_idx = c(rhythmic_same, DP_genes, DA_genes, DR_g1_only)

  if (length(rhythmic_idx) > 0) {
    # phi_g ~ Uniform(0, 24)  (no preferred peak time)
    phase1[rhythmic_idx] = runif(length(rhythmic_idx), 0, period)
    # Amplitude: resample from empirical vector if provided, else parametric default
    if (!missing(amplitude) && !is.null(sigma_rhythmic) &&
        length(sigma_rhythmic) == length(amplitude)) {
      # Joint sampling: draw A and sigma from the same pilot gene to preserve
      # empirical A-sigma correlation (avoids marginal sampling bias in r = A/sigma).
      joint_idx <- sample(length(amplitude), length(rhythmic_idx), replace = TRUE)
      amplitude1[rhythmic_idx] <- pmax(amplitude[joint_idx], 0.05)
      sigma[rhythmic_idx]      <- pmax(sigma_rhythmic[joint_idx], 1e-6)
    } else if (!missing(amplitude)) {
      amplitude1[rhythmic_idx] = pmax(sample(amplitude, length(rhythmic_idx), replace = TRUE), 0.05)
    } else {
      # A_g ~ max(LogNormal(mu=log(0.4), sd=0.5), 0.05)
      # median(A) = 0.4,  IQR ~ [0.24, 0.66]
      # Combined with sigma: r = A/sigma ~ LogNormal(mu~0.08, sd~0.58)
      #   => median(r) ~ 1.08, bulk of genes in r = 0.3 - 2.0
      amplitude1[rhythmic_idx] = pmax(rlnorm(length(rhythmic_idx), log(0.4), 0.5), 0.05)
    }
  }

  # Generate group 2 parameters
  phase2 = phase1
  amplitude2 = amplitude1

  # DR: g1 only - group 2 is arrhythmic
  if (length(DR_g1_only) > 0) {
    amplitude2[DR_g1_only] = 0
  }

  # DR: g2 only - group 1 is arrhythmic
  if (length(DR_g2_only) > 0) {
    amplitude1[DR_g2_only] = 0
    phase2[DR_g2_only] = runif(length(DR_g2_only), 0, period)
    # Use group-2 amplitude distribution (amplitude2) when available;
    # fall back to group-1 distribution (amplitude) otherwise.
    amp_src2 <- if (!is.null(amp_dist2)) amp_dist2 else
                if (!missing(amplitude)) amplitude else NULL
    if (!is.null(amp_src2)) {
      amplitude2[DR_g2_only] = pmax(sample(amp_src2, length(DR_g2_only), replace = TRUE), 0.05)
    } else {
      amplitude2[DR_g2_only] = pmax(rlnorm(length(DR_g2_only), log(0.4), 0.5), 0.05)
    }
  }

  # DP: Shift phase for group 2
  if (length(DP_genes) > 0) {
    if (dp_shift_mode == "uniform") {
      phase_shift = runif(length(DP_genes), phase_diff[1], phase_diff[2])
    } else {
      # Fixed shift = phase_diff[2]
      phase_shift = rep(phase_diff[2], length(DP_genes))
    }
    phase2[DP_genes] = (phase1[DP_genes] + phase_shift) %% period
  }

  # DA: Modify amplitude for group 2
  if (length(DA_genes) > 0) {
    amp_ratio = runif(length(DA_genes), amp_diff[1], amp_diff[2])
    amplitude2[DA_genes] = amplitude1[DA_genes] * amp_ratio
  }

  # ---------------------------------------------------------------
  # Now reseed for n-dependent sampling (time points + expression noise).
  # This ensures gene parameters above are identical across sample sizes.
  # ---------------------------------------------------------------
  set.seed(sim.seed + n1 * 7919L + n2 * 104729L)

  # Generate time points (depends on n1, n2)
  if (design == "active") {
    # If fixed times are provided, use them (enables replicated active designs).
    if (!is.null(cts)) {
      if (length(cts) != n1) stop("Active design: length(cts) must equal n1")
      times1 = cts
      if (!is.null(cts2)) {
        if (length(cts2) != n2) stop("Active design: length(cts2) must equal n2")
        times2 = cts2
      } else {
        times2 = cts
      }
    } else {
      # Exclude endpoint (period ≡ 0 mod period) for true evenly-spaced coverage.
      times1 = seq(0, period, length.out = n1 + 1L)[seq_len(n1)]
      times2 = seq(0, period, length.out = n2 + 1L)[seq_len(n2)]
    }
  } else {
    if (is.null(cts)) {
      stop("For passive design, must provide cts (TOD distribution)")
    }
    times1 = sampleTimesFromDist(n1, cts)
    times2 = sampleTimesFromDist(n2, if (!is.null(cts2)) cts2 else cts)
  }

  # Generate expression data
  expr1 = matrix(NA, nrow = ngenes, ncol = n1)
  expr2 = matrix(NA, nrow = ngenes, ncol = n2)

  for (g in 1:ngenes) {
    # Group 1
    mu1 = mesor[g] + amplitude1[g] * cos(omega * times1 - omega * phase1[g])
    if (!is.null(harmonics) && length(harmonics) >= 1 && harmonics[1] != 0) {
      mu1 = mu1 + amplitude1[g] * harmonics[1] * cos(2 * omega * times1 - omega * phase1[g])
    }
    if (!is.null(harmonics) && length(harmonics) >= 2 && harmonics[2] != 0) {
      mu1 = mu1 + amplitude1[g] * harmonics[2] * cos(3 * omega * times1 - omega * phase1[g])
    }
    expr1[g, ] = rnorm(n1, mu1, sigma[g])

    # Group 2
    mu2 = mesor[g] + amplitude2[g] * cos(omega * times2 - omega * phase2[g])
    if (!is.null(harmonics) && length(harmonics) >= 1 && harmonics[1] != 0) {
      mu2 = mu2 + amplitude2[g] * harmonics[1] * cos(2 * omega * times2 - omega * phase2[g])
    }
    if (!is.null(harmonics) && length(harmonics) >= 2 && harmonics[2] != 0) {
      mu2 = mu2 + amplitude2[g] * harmonics[2] * cos(3 * omega * times2 - omega * phase2[g])
    }
    expr2[g, ] = rnorm(n2, mu2, sigma2_vec[g])
  }

  rownames(expr1) = rownames(expr2) = paste0("Gene", 1:ngenes)
  colnames(expr1) = paste0("G1_", 1:n1)
  colnames(expr2) = paste0("G2_", 1:n2)

  # Create detailed ground truth
  # Calculate phase difference with wrap-around correction
  phase_diff_calc = (phase2 - phase1) %% period
  phase_diff_calc[phase_diff_calc > 12] = phase_diff_calc[phase_diff_calc > 12] - 24

  ground_truth = data.frame(
    gene = 1:ngenes,
    diff_type = diff_type,
    diff_type_label = c("Non-rhythmic both", "Rhythmic both (same)",
                     "DR: G1 only", "DR: G2 only",
                     "DP: Diff phase", "DA: Diff amp")[diff_type + 1],
    mesor = mesor,
    amplitude1 = amplitude1,
    amplitude2 = amplitude2,
    phase1 = phase1,
    phase2 = phase2,
    sigma = sigma,          # group-1 noise (backwards compat)
    sigma2 = sigma2_vec,   # group-2 noise (= sigma when lOD2 not provided)
    phase_diff = phase_diff_calc,
    amp_ratio = ifelse(amplitude1 > 0, amplitude2 / amplitude1, NA),
    is_rhythmic_g1 = amplitude1 > 0,
    is_rhythmic_g2 = amplitude2 > 0,
    stringsAsFactors = FALSE
  )

  # Effect sizes for detection (group-specific sigma when lOD2 is provided)
  effectsize_DR1 = amplitude1 / sigma        # r_{g,1} = A_{g,1} / sigma_{g,1}
  effectsize_DR2 = amplitude2 / sigma2_vec   # r_{g,2} = A_{g,2} / sigma_{g,2}
  effectsize_phase = sqrt(
    amplitude1^2 + amplitude2^2 -
      2*amplitude1*amplitude2*cos(omega*(phase2-phase1))
  ) / sigma   # displacement relative to group-1 noise as reference
  effectsize_amp = abs(amplitude2 - amplitude1) / sigma  # relative to group-1 noise

  return(list(
    expr1 = expr1,
    expr2 = expr2,
    times1 = times1,
    times2 = times2,
    ground_truth = ground_truth,
    effectsize_DR1 = effectsize_DR1,
    effectsize_DR2 = effectsize_DR2,
    effectsize_phase = effectsize_phase,
    effectsize_amp = effectsize_amp,
    simOptions = list(
      ngenes = ngenes,
      n1 = n1,
      n2 = n2,
      prop_rhythmic = prop_rhythmic,
      prop_DR = prop_DR,
      prop_DP = prop_DP,
      prop_DA = prop_DA,
      period = period,
      design = design
    )
  ))
}
#' Simulate Circadian Time-Series Data
#'
#' @description Generate simulated circadian expression data with flexible
#' noise models and waveform shapes.
#'
#' @param G Number of genes
#' @param G_rhythmic Number of rhythmic genes (default: floor(G * prop_rhythmic))
#' @param prop_rhythmic Proportion of rhythmic genes (default 0.1)
#' @param n Number of samples per time point (replicates)
#' @param times Vector of time points (in hours)
#' @param period Period (default 24)
#' @param params List of parameters (M, A, phi, sigma). If NULL, use defaults.
#' @param noise_type "gaussian", "t", "negbinom", "lognormal"
#' @param waveform "sinusoid", "damped", "asymmetric"
#' @param seed Random seed
#'
#' @return List with data matrix, ground truth, and parameters
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
