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
#'   \item{Type 5}{DM: Differential mesor (both rhythmic, same A/phi, different mean level)}
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
#' @param phase_diff Range of phase shift (hours) for DP genes
#' @param amp_diff Unused; retained for interface compatibility
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
                           lBaselineExpr2 = NULL,
                           lOD,
                           lOD2 = NULL,
                           amplitude,
                           amplitude2 = NULL,
                           sigma_rhythmic = NULL,
                           prop_rhythmic = 0.25,
                           prop_DR = 0.1,
                           prop_DP = 0.2,
                           prop_DM = 0,
                           phase_diff = c(-6, 6),
                           dp_shift_mode = c("fixed", "uniform"),
                           amp_diff = c(0.5, 2),
                           mesor_diff = c(0.5, 2.0),
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
  # Group-2 baseline: use lBaselineExpr2 if provided (two-group pilot), else share group-1
  if (!is.null(lBaselineExpr2)) {
    lBaselineExpr2 <- if (length(lBaselineExpr2) == ngenes) lBaselineExpr2
                      else sample(lBaselineExpr2, ngenes, replace = TRUE)
  }
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
  # 4=DP (diff phase), 5=DM (diff mesor, same A/phi)
  diff_type = rep(0, ngenes)

  # Count differentially rhythmic genes
  n_DR = round(ngenes * prop_DR)
  n_DP = round(ngenes * prop_DP)
  n_DM = round(ngenes * prop_DM)

  # Assign categories
  gene_idx = 1:ngenes

  # Differential Rhythmicity (DR): rhythmic in g1 only or g2 only
  if (n_DR > 0) {
    DR_genes  = sample(gene_idx, n_DR)
    split_idx = ceiling(n_DR / 2)
    DR_g1_only = DR_genes[seq_len(split_idx)]
    DR_g2_only = if (split_idx < n_DR) DR_genes[(split_idx + 1L):n_DR] else integer(0)
    diff_type[DR_g1_only] = 2
    diff_type[DR_g2_only] = 3
  } else {
    DR_genes = integer(0)
    DR_g1_only = DR_g2_only = integer(0)
  }

  # Rhythmic in both (for DP, DM, or same)
  # prop_rhythmic = total fraction rhythmic in at least one group.
  # DR genes already count toward that budget; rhythmic_both fills the remainder.
  remaining = setdiff(gene_idx, DR_genes)
  n_rhythmic_total <- round(ngenes * prop_rhythmic)
  n_rhythmic_both  <- max(n_DP + n_DM, n_rhythmic_total - n_DR)  # must fit DP + DM
  n_rhythmic_both  <- min(n_rhythmic_both, length(remaining))
  rhythmic_both    <- if (n_rhythmic_both > 0) sample(remaining, n_rhythmic_both) else integer(0)

  # Differential Phase (DP)
  if (n_DP > 0 && length(rhythmic_both) > 0) {
    DP_genes = sample(setdiff(rhythmic_both, DR_genes), min(n_DP, length(rhythmic_both)))
    diff_type[DP_genes] = 4
  } else {
    DP_genes = integer(0)
  }

  # Differential Mesor (DM) — Type 5: both rhythmic, same A/phi, different mean
  if (n_DM > 0 && length(rhythmic_both) > 0) {
    DM_pool = setdiff(rhythmic_both, c(DR_genes, DP_genes))
    n_DM_actual = min(n_DM, length(DM_pool))
    DM_genes = sample(DM_pool, n_DM_actual)
    diff_type[DM_genes] = 5L
  } else {
    DM_genes = integer(0)
  }

  # Rhythmic in both, same (control genes)
  rhythmic_same = setdiff(rhythmic_both, c(DP_genes, DM_genes))
  if (length(rhythmic_same) > 0) {
    diff_type[rhythmic_same] = 1
  }

  # Generate base parameters (group 1)
  phase1 = rep(0, ngenes)
  amplitude1 = rep(0, ngenes)

  # All rhythmic genes get base phase and amplitude
  rhythmic_idx = c(rhythmic_same, DP_genes, DM_genes, DR_g1_only)

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
      if (!is.null(sigma_rhythmic) && length(sigma_rhythmic) != length(amplitude)) {
        warning(sprintf(
          "sigma_rhythmic length (%d) != amplitude length (%d): falling back to marginal sampling. ",
          length(sigma_rhythmic), length(amplitude)),
          "Empirical A-sigma correlation from pilot will not be preserved.",
          call. = FALSE)
      }
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

  # DM: group-2 mesor is shifted; amplitude and phase unchanged.
  # mesor2 starts as a copy of group-1 mesor (or two-pilot group-2 baseline).
  mesor2 <- if (!is.null(lBaselineExpr2)) lBaselineExpr2 else mesor
  if (length(DM_genes) > 0) {
    # Draw shift magnitude from Uniform[mesor_diff[1], mesor_diff[2]] (log-scale units).
    # Sign is random (up- or down-shift equally likely) to avoid directional bias.
    shift_mag  <- runif(length(DM_genes), mesor_diff[1], mesor_diff[2])
    shift_sign <- sample(c(-1L, 1L), length(DM_genes), replace = TRUE)
    mesor2[DM_genes] <- mesor[DM_genes] + shift_sign * shift_mag
  }

  # ---------------------------------------------------------------
  # Now reseed for n-dependent sampling (time points + expression noise).
  # This ensures gene parameters above are identical across sample sizes.
  # ---------------------------------------------------------------
  set.seed(sim.seed + n1 * 7919L + n2 * 104729L)

  # Generate time points (depends on n1, n2)
  if (design == "active") {
    # If fixed times are provided, use them (enables replicated active designs).
    # If cts is a ZT template shorter than n1 (e.g., a 6-ZT grid passed from
    # runner.R before expansion), expand it by repeating to match n1.
    if (!is.null(cts)) {
      if (length(cts) != n1) {
        cts <- sort(rep_len(cts, n1))
      }
      times1 = cts
      if (!is.null(cts2)) {
        if (length(cts2) != n2) {
          cts2 <- sort(rep_len(cts2, n2))
        }
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
      mu1 = mu1 + amplitude1[g] * harmonics[1] * cos(2 * omega * times1 - 2 * omega * phase1[g])
    }
    if (!is.null(harmonics) && length(harmonics) >= 2 && harmonics[2] != 0) {
      mu1 = mu1 + amplitude1[g] * harmonics[2] * cos(3 * omega * times1 - 3 * omega * phase1[g])
    }
    expr1[g, ] = rnorm(n1, mu1, sigma[g])

    # Group 2 (uses mesor2 which equals mesor unless two-pilot or DM perturbation)
    mu2 = mesor2[g] + amplitude2[g] * cos(omega * times2 - omega * phase2[g])
    if (!is.null(harmonics) && length(harmonics) >= 1 && harmonics[1] != 0) {
      mu2 = mu2 + amplitude2[g] * harmonics[1] * cos(2 * omega * times2 - 2 * omega * phase2[g])
    }
    if (!is.null(harmonics) && length(harmonics) >= 2 && harmonics[2] != 0) {
      mu2 = mu2 + amplitude2[g] * harmonics[2] * cos(3 * omega * times2 - 3 * omega * phase2[g])
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
                     "DP: Diff phase", "DM: Diff mesor")[diff_type + 1],
    mesor1 = mesor,
    mesor2 = mesor2,
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
  effectsize_amp   = abs(amplitude2 - amplitude1) / sigma  # relative to group-1 noise
  effectsize_mesor = abs(mesor2 - mesor) / sigma           # |Δμ| / sigma_g1

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
    effectsize_mesor = effectsize_mesor,
    simOptions = list(
      ngenes = ngenes,
      n1 = n1,
      n2 = n2,
      prop_rhythmic = prop_rhythmic,
      prop_DR = prop_DR,
      prop_DP = prop_DP,
      prop_DM = prop_DM,
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
simCircadianSingleCohort <- function(bio.opts, cts, alpha2 = 0, alpha3 = 0,
                                     seed = NULL) {
  stopifnot(inherits(bio.opts, "CircadianBioOptions"))

  if (!is.null(seed)) set.seed(seed)

  ngenes        <- bio.opts$ngenes
  prop_rhythmic <- bio.opts$prop_rhythmic
  period        <- bio.opts$period %||% 24
  N             <- length(cts)
  omega         <- 2 * pi / period

  has_joint <- !is.null(bio.opts$sigma_rhythmic) &&
               length(bio.opts$sigma_rhythmic) == length(bio.opts$amplitude)

  # --- gene assignment ---
  n_rhythmic  <- round(ngenes * prop_rhythmic)
  rhythmic_id <- sample(ngenes, n_rhythmic)
  is_rhythmic <- logical(ngenes)
  is_rhythmic[rhythmic_id] <- TRUE

  # --- parameter draws ---
  mesor_g <- bio.opts$lBaselineExpr
  sigma_g <- exp(bio.opts$lOD)
  amp_g   <- numeric(ngenes)
  phase_g <- numeric(ngenes)

  if (n_rhythmic > 0) {
    if (has_joint) {
      ji <- sample(length(bio.opts$amplitude), n_rhythmic, replace = TRUE)
      amp_g[rhythmic_id]   <- pmax(bio.opts$amplitude[ji], 0.05)
      sigma_g[rhythmic_id] <- pmax(bio.opts$sigma_rhythmic[ji], 1e-6)
    } else {
      amp_g[rhythmic_id] <- pmax(
        sample(bio.opts$amplitude, n_rhythmic, replace = TRUE), 0.05)
    }
    phase_g[rhythmic_id] <- sample(bio.opts$phase, n_rhythmic, replace = TRUE)
  }

  r_values <- amp_g / sigma_g

  # --- expression matrix ---
  expr <- matrix(NA_real_, nrow = ngenes, ncol = N)
  for (g in seq_len(ngenes)) {
    mu <- mesor_g[g] +
          amp_g[g] * (cos(omega * cts - omega * phase_g[g]) +
                      alpha2 * cos(2 * omega * cts - 2 * omega * phase_g[g]) +
                      alpha3 * cos(3 * omega * cts - 3 * omega * phase_g[g]))
    expr[g, ] <- rnorm(N, mu, sigma_g[g])
  }

  list(expr = expr, is_rhythmic = is_rhythmic, r_values = r_values)
}
