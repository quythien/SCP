#' Sample Time Points from an Empirical TOD Distribution
#'
#' @description
#' Draws \code{n} time-of-day (TOD) values from a circular kernel density
#' estimate of a pilot TOD distribution. Handles wrap-around at 0h/24h
#' using a tripled-data trick so boundary contributions are preserved.
#'
#' @param n Integer. Number of time points to sample.
#' @param cts Numeric vector. Pilot TOD values (hours, any range — wrapped to
#'   \code{[0, 24)} internally).
#'
#' @return Numeric vector of length \code{n} with sampled TOD values in
#'   \code{[0, 24)}.
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
  probs   <- pmax(tod_density$y, 1e-10)   # floor prevents zero-sum crash on tight clusters
  sampled <- sample(tod_density$x, size = n, replace = TRUE, prob = probs)

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
#' @return Named list with the following elements:
#'   \describe{
#'     \item{\code{expr1}}{Gene expression matrix for group 1 (ngenes x n1).}
#'     \item{\code{expr2}}{Gene expression matrix for group 2 (ngenes x n2).}
#'     \item{\code{times1}}{Numeric vector of sample times for group 1 (hours).}
#'     \item{\code{times2}}{Numeric vector of sample times for group 2 (hours).}
#'     \item{\code{ground_truth}}{data.frame (ngenes rows) with columns:
#'       \code{gene}, \code{diff_type} (0–5), \code{diff_type_label},
#'       \code{mesor1}, \code{mesor2}, \code{amplitude1}, \code{amplitude2},
#'       \code{phase1}, \code{phase2}, \code{sigma}, \code{sigma2},
#'       \code{phase_diff}, \code{amp_ratio}, \code{is_rhythmic_g1},
#'       \code{is_rhythmic_g2}.}
#'     \item{\code{effectsize_DR1}}{Per-gene SNR in group 1: \eqn{A_{g,1}/\sigma_{g,1}}.}
#'     \item{\code{effectsize_DR2}}{Per-gene SNR in group 2: \eqn{A_{g,2}/\sigma_{g,2}}.}
#'     \item{\code{effectsize_phase}}{Per-gene phase-displacement effect size.}
#'     \item{\code{effectsize_amp}}{Per-gene amplitude-difference effect size.}
#'     \item{\code{effectsize_mesor}}{Per-gene mesor-difference effect size.}
#'     \item{\code{simOptions}}{List of simulation parameters used (ngenes, n1, n2,
#'       prop_rhythmic, prop_DR, prop_DP, prop_DM, period, design).}
#'   }

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
  a2 <- if (!is.null(harmonics) && length(harmonics) >= 1) harmonics[1] else 0
  a3 <- if (!is.null(harmonics) && length(harmonics) >= 2) harmonics[2] else 0
  if (exists(".CPP_LOADED", inherits = TRUE) && isTRUE(get(".CPP_LOADED", inherits = TRUE)) &&
      exists("sim_cosinor_expr_fast", mode = "function")) {
    expr1 <- sim_cosinor_expr_fast(mesor,  amplitude1, phase1, sigma,      times1, period, a2, a3)
    expr2 <- sim_cosinor_expr_fast(mesor2, amplitude2, phase2, sigma2_vec, times2, period, a2, a3)
  } else {
    expr1 = matrix(NA, nrow = ngenes, ncol = n1)
    expr2 = matrix(NA, nrow = ngenes, ncol = n2)
    for (g in 1:ngenes) {
      mu1 = mesor[g] + amplitude1[g] * cos(omega * times1 - omega * phase1[g])
      if (a2 != 0) mu1 = mu1 + amplitude1[g] * a2 * cos(2 * omega * times1 - 2 * omega * phase1[g])
      if (a3 != 0) mu1 = mu1 + amplitude1[g] * a3 * cos(3 * omega * times1 - 3 * omega * phase1[g])
      expr1[g, ] = rnorm(n1, mu1, sigma[g])
      mu2 = mesor2[g] + amplitude2[g] * cos(omega * times2 - omega * phase2[g])
      if (a2 != 0) mu2 = mu2 + amplitude2[g] * a2 * cos(2 * omega * times2 - 2 * omega * phase2[g])
      if (a3 != 0) mu2 = mu2 + amplitude2[g] * a3 * cos(3 * omega * times2 - 3 * omega * phase2[g])
      expr2[g, ] = rnorm(n2, mu2, sigma2_vec[g])
    }
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
#' Simulate Single-Cohort Circadian Gene Expression Data
#'
#' @description Generates a synthetic gene expression matrix for one cohort
#' by drawing per-gene parameters (mesor, amplitude, phase, noise) from the
#' empirical distributions stored in a \code{CircadianBioOptions} object.
#' Rhythmic genes receive nonzero amplitude and a randomly drawn peak phase;
#' non-rhythmic genes have amplitude fixed at zero.
#'
#' @param bio.opts A \code{CircadianBioOptions} object from \code{estCircadianParam}.
#' @param cts Numeric vector of sample collection times (hours). Length determines
#'   the number of simulated samples (\code{N}).
#' @param alpha2 Second-harmonic coefficient (default 0). Adds a
#'   \eqn{A \cdot \alpha_2 \cos(4\pi t/T - 4\pi\phi/T)} term.
#'   Only used when \code{omega = 1} (cosinor path). Ignored when \code{omega < 1}.
#' @param alpha3 Third-harmonic coefficient (default 0). Same conditions as \code{alpha2}.
#' @param omega FMM waveform shape parameter in \eqn{(0, 1]} (default \code{1.0}).
#'   \code{omega = 1} uses the traditional cosinor generator (with \code{alpha2}/\code{alpha3}).
#'   \code{omega < 1} uses \code{\link{simCircadianFMM}} — a non-sinusoidal Möbius waveform.
#'   The two paths are mutually exclusive: FMM takes precedence when \code{omega < 1}.
#' @param beta FMM orientation parameter (default \eqn{\pi}). Controls the phase of
#'   waveform asymmetry. Only used when \code{omega < 1}. \eqn{\beta}-invariance
#'   ensures DCP detection power is unaffected by \code{beta}.
#' @param seed Optional integer random seed for reproducibility.
#'
#' @return Named list with:
#'   \describe{
#'     \item{\code{expr}}{Gene expression matrix (ngenes x N).}
#'     \item{\code{is_rhythmic}}{Logical vector of length ngenes; TRUE for rhythmic genes.}
#'     \item{\code{r_values}}{Per-gene signal-to-noise ratio \eqn{A_g / \sigma_g}.}
#'   }
#'
#' @seealso \code{\link{estCircadianParam}}, \code{\link{simCircadianDiff}}
simCircadianSingleCohort <- function(bio.opts, cts, alpha2 = 0, alpha3 = 0,
                                     omega = 1.0, beta = pi,
                                     seed = NULL) {
  stopifnot(inherits(bio.opts, "CircadianBioOptions"))
  if (omega <= 0 || omega > 1)
    stop("omega must be in (0, 1]. Use omega = 1 (default) for cosinor simulation.")

  # FMM and Fourier harmonic paths are mutually exclusive.
  # omega < 1 activates the FMM generator; alpha2/alpha3 are then ignored.
  if (omega < 1.0 && (alpha2 != 0 || alpha3 != 0))
    warning("simCircadianSingleCohort: omega < 1 activates FMM simulation; ",
            "alpha2/alpha3 are ignored. Set omega = 1 to use the Fourier harmonic path.")

  if (!is.null(seed)) set.seed(seed)

  # FMM path is also activated when bio.opts has per-gene FMM parameters
  # (omega_dist, alpha_dist, omega_rhythmic, or alpha_rhythmic) — these are
  # only respected by simCircadianFMM, not the cosinor generator.
  has_fmm_per_gene <- !is.null(bio.opts$omega_dist) ||
                       !is.null(bio.opts$alpha_dist) ||
                       (!is.null(bio.opts$omega_rhythmic) &&
                        length(bio.opts$omega_rhythmic) > 0) ||
                       (!is.null(bio.opts$alpha_rhythmic) &&
                        length(bio.opts$alpha_rhythmic) > 0)

  # FMM path: non-sinusoidal waveform when omega < 1 OR per-gene FMM params present
  if (omega < 1.0 || has_fmm_per_gene)
    return(simCircadianFMM(bio.opts, cts, omega = omega, beta = beta))

  ngenes        <- bio.opts$ngenes
  prop_rhythmic <- bio.opts$prop_rhythmic
  period        <- bio.opts$period %||% 24
  N             <- length(cts)
  omega_circ    <- 2 * pi / period   # circular frequency — distinct from FMM omega param

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

  # Resample pilot vectors if ngenes was overridden (e.g. bio$ngenes <- 5000)
  # but lBaselineExpr/lOD still have pilot K length.
  if (length(mesor_g) != ngenes) {
    idx     <- sample(length(mesor_g), ngenes, replace = TRUE)
    mesor_g <- mesor_g[idx]
    sigma_g <- sigma_g[idx]
  }

  amp_g   <- numeric(ngenes)
  phase_g <- numeric(ngenes)

  if (n_rhythmic > 0) {
    if (has_joint) {
      ji <- sample(length(bio.opts$amplitude), n_rhythmic, replace = TRUE)
      amp_g[rhythmic_id]   <- pmax(bio.opts$amplitude[ji], 0.05)
      sigma_g[rhythmic_id] <- pmax(bio.opts$sigma_rhythmic[ji], 1e-6)
    } else {
      amp_src <- bio.opts$amplitude
      amp_g[rhythmic_id] <- pmax(
        if (length(amp_src) == 1L) rep(amp_src, n_rhythmic)
        else sample(amp_src, n_rhythmic, replace = TRUE), 0.05)
    }
    phase_src <- bio.opts$phase
    phase_g[rhythmic_id] <- if (length(phase_src) == 1L)
      rep(phase_src, n_rhythmic)
    else
      sample(phase_src, n_rhythmic, replace = TRUE)
  }

  r_values <- amp_g / sigma_g

  # --- expression matrix ---
  if (exists(".CPP_LOADED", inherits = TRUE) && isTRUE(get(".CPP_LOADED", inherits = TRUE)) &&
      exists("sim_cosinor_expr_fast", mode = "function")) {
    expr <- sim_cosinor_expr_fast(mesor_g, amp_g, phase_g, sigma_g, cts, period, alpha2, alpha3)
  } else {
    expr <- matrix(NA_real_, nrow = ngenes, ncol = N)
    for (g in seq_len(ngenes)) {
      mu <- mesor_g[g] +
            amp_g[g] * (cos(omega_circ * cts - omega_circ * phase_g[g]) +
                        alpha2 * cos(2 * omega_circ * cts - 2 * omega_circ * phase_g[g]) +
                        alpha3 * cos(3 * omega_circ * cts - 3 * omega_circ * phase_g[g]))
      expr[g, ] <- rnorm(N, mu, sigma_g[g])
    }
  }

  list(expr = expr, is_rhythmic = is_rhythmic, r_values = r_values)
}

#' Simulate Single-Cohort Circadian Data Under FMM Waveform
#'
#' @description Generates a synthetic gene expression matrix for one cohort
#' using the Frequency Modulated Möbius (FMM) model.  The waveform shape is
#' controlled by \code{omega}: \code{omega = 1} reproduces a pure sinusoid
#' (identical to DCP cosinor truth), \code{omega = 0} produces a flat signal
#' (no rhythmicity), and intermediate values produce a peaked waveform.
#' Gene-level parameters (mesor, amplitude, phase, noise) are drawn from the
#' empirical distributions in \code{bio.opts} exactly as in
#' \code{\link{simCircadianSingleCohort}}.
#'
#' \strong{Beta-invariance:} The orientation parameter \code{beta} shifts the
#' peak location on the 24-hour cycle but does not alter the amplitude of the
#' fitted cosinor signal.  DCP amplitude-based detection power is therefore
#' invariant to \code{beta}; only phase estimation is affected.
#'
#' @references Rueda C, Rodríguez-Collado A, Peddada SD (2019).
#'   "A single notable oscillation in humans and mammals with strong
#'   implications in chronobiology." Sci Rep 9, 17982.
#'   \doi{10.1038/s41598-019-54569-1}
#'
#' @param bio.opts A \code{CircadianBioOptions} object from
#'   \code{\link{estCircadianParam}}.
#' @param cts Numeric vector of sample collection times in hours \eqn{[0, 24)}.
#'   The length of this vector determines the number of simulated samples
#'   (\eqn{N}).
#' @param omega Numeric scalar in \eqn{[0, 1]}.  FMM waveform shape parameter.
#'   \code{omega = 1} gives a pure sinusoid (cosinor truth);
#'   \code{omega = 0} gives a flat (arrhythmic) signal;
#'   values in \eqn{(0, 1)} give increasingly non-sinusoidal peaked waveforms.
#'   Defaults to \code{1.0}.
#' @param beta Numeric scalar.  FMM orientation (skewness) parameter. Default
#'   \eqn{\pi} gives a symmetric peak.  Does not affect amplitude-based power
#'   (beta-invariance); affects only peak-time (phase) estimation.
#' @param seed Optional integer random seed for reproducibility.
#'
#' @return Named list with:
#'   \describe{
#'     \item{\code{expr}}{Gene expression matrix (\code{ngenes} x \eqn{N}).}
#'     \item{\code{is_rhythmic}}{Logical vector of length \code{ngenes}; TRUE
#'       for rhythmic genes.}
#'     \item{\code{r_values}}{Per-gene signal-to-noise ratio
#'       \eqn{A_g / \sigma_g}.}
#'     \item{\code{omega}}{The \code{omega} value used.}
#'   }
#'
#' @seealso \code{\link{simCircadianSingleCohort}}, \code{\link{estCircadianParam}},
#'   \code{\link{plotFMMViolation}}
simCircadianFMM <- function(bio.opts, cts, omega = 1.0, beta = pi,
                            alpha_fixed = NULL, seed = NULL) {
  stopifnot(inherits(bio.opts, "CircadianBioOptions"))
  if (omega < 0 || omega > 1)
    stop("omega must be in [0, 1]; got ", omega)

  if (!is.null(seed)) set.seed(seed)

  ngenes        <- bio.opts$ngenes
  prop_rhythmic <- bio.opts$prop_rhythmic
  N             <- length(cts)

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

  # Resample pilot vectors if ngenes was overridden
  if (length(mesor_g) != ngenes) {
    idx     <- sample(length(mesor_g), ngenes, replace = TRUE)
    mesor_g <- mesor_g[idx]
    sigma_g <- sigma_g[idx]
  }

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
    if (!is.null(alpha_fixed)) {
      phase_g[rhythmic_id] <- alpha_fixed   # override: all rhythmic genes share this acrophase
    } else {
      phase_g[rhythmic_id] <- sample(bio.opts$phase, n_rhythmic, replace = TRUE)
    }
  }

  r_values <- amp_g / sigma_g

  # --- Resolve per-gene omega_g and alpha_g (radians) for rhythmic genes ---
  # Priority for omega: omega_dist > omega_rhythmic > scalar `omega` arg
  # Priority for alpha: alpha_dist (with optional jitter on alpha_rhythmic) >
  #                     alpha_rhythmic > derived from phase_g (acrophase, radians)
  # All clipped/wrapped to safe ranges.
  if (n_rhythmic > 0) {
    .clip_omega <- function(x) pmin(pmax(x, 1e-4), 1 - 1e-4)
    omega_g <- if (!is.null(bio.opts$omega_dist)) {
      spec <- bio.opts$omega_dist
      raw  <- switch(spec$family,
        beta  = rbeta(n_rhythmic, spec$a, spec$b),
        fixed = rep(spec$value, n_rhythmic),
        stop("Unknown omega_dist$family: ", spec$family))
      .clip_omega(raw)
    } else if (!is.null(bio.opts$omega_rhythmic) &&
               length(bio.opts$omega_rhythmic) == n_rhythmic) {
      .clip_omega(bio.opts$omega_rhythmic)
    } else {
      rep(omega, n_rhythmic)  # backward-compat scalar (no clip — caller validated)
    }

    # Determine the "base" acrophase in radians (from phase_g, in hours, → radians)
    alpha_base_rad <- phase_g[rhythmic_id] * (2 * pi / 24)

    # Phase noise: when sd_hours is small the wrapped-normal approximates
    # Gaussian jitter; for sd_hours >= 3 it becomes effectively circular-uniform.
    # Use von Mises (the principled circular analog of Gaussian noise on a
    # circle) so the alpha sweep is statistically defensible at any σ.
    # rvonmises (Best & Fisher 1979 rejection method, simple implementation).
    .rvonmises <- function(n, mu, kappa) {
      if (kappa <= 0 || !is.finite(kappa)) return(runif(n, 0, 2 * pi))
      a <- 1 + sqrt(1 + 4 * kappa^2)
      b <- (a - sqrt(2 * a)) / (2 * kappa)
      r <- (1 + b^2) / (2 * b)
      out <- numeric(n); k <- 0L
      while (k < n) {
        z <- cos(pi * runif(1))
        f <- (1 + r * z) / (r + z)
        c <- kappa * (r - f)
        u <- runif(1)
        if (u < c * (2 - c) || log(c / u) + 1 - c >= 0) {
          k <- k + 1L
          theta <- if (runif(1) - 0.5 < 0) -acos(f) else acos(f)
          out[k] <- (mu + theta) %% (2 * pi)
        }
      }
      out
    }
    alpha_g <- if (!is.null(bio.opts$alpha_dist)) {
      spec   <- bio.opts$alpha_dist
      sd_hr  <- spec$sd_hours %||% spec$sd %||% 0
      sd_rad <- sd_hr * (2 * pi / 24)   # hours → radians
      mu     <- spec$mean %||% 0
      if (spec$family != "normal")
        stop("Unknown alpha_dist$family: ", spec$family)
      # When sd_rad is small (< pi/6 ≈ 1 hr), Gaussian and von Mises coincide
      # to first order. We use von Mises throughout for principled circular noise.
      kappa <- if (sd_rad <= 0) Inf else 1 / sd_rad^2
      base_alpha <- if (!is.null(bio.opts$alpha_rhythmic) &&
                        length(bio.opts$alpha_rhythmic) == n_rhythmic) {
        bio.opts$alpha_rhythmic
      } else {
        alpha_base_rad
      }
      if (sd_rad == 0) {
        base_alpha
      } else {
        # Add circular noise gene-by-gene around each gene's empirical alpha
        sapply(seq_len(n_rhythmic), function(i)
          .rvonmises(1, mu = base_alpha[i] + mu, kappa = kappa))
      }
    } else if (!is.null(bio.opts$alpha_rhythmic) &&
               length(bio.opts$alpha_rhythmic) == n_rhythmic) {
      bio.opts$alpha_rhythmic
    } else {
      alpha_base_rad
    }

    # Wrap to [0, 2pi) for FMM
    alpha_g <- alpha_g %% (2 * pi)
  } else {
    omega_g <- numeric(0)
    alpha_g <- numeric(0)
  }

  # --- FMM time grid (radians): FMM package uses [0, 2*pi] ---
  cts_rad <- cts * (2 * pi / 24)

  # Per-gene rhythmic index lookup: for gene g, find its position in rhythmic_id
  rhythmic_pos <- integer(ngenes)
  rhythmic_pos[rhythmic_id] <- seq_len(n_rhythmic)

  # --- expression matrix ---
  expr <- matrix(NA_real_, nrow = ngenes, ncol = N)

  for (g in seq_len(ngenes)) {
    if (is_rhythmic[g]) {
      gp        <- rhythmic_pos[g]   # 1..n_rhythmic position
      omega_use <- omega_g[gp]
      alpha_use <- alpha_g[gp]

      if (omega_use > 0) {
        fmm_out <- FMM::generateFMM(
          M          = mesor_g[g],
          A          = amp_g[g],
          alpha      = alpha_use,
          beta       = beta,
          omega      = omega_use,
          from       = 0,
          to         = 2 * pi + 1e-4,
          length.out = 1000L,
          plot       = FALSE,
          outvalues  = TRUE,
          sigmaNoise = 0
        )
        fmm_signal <- approx(fmm_out$t, fmm_out$y, xout = cts_rad,
                             method = "linear", rule = 2)$y
        expr[g, ] <- rnorm(N, fmm_signal, sigma_g[g])
      } else {
        # omega=0 → flat (FMM degenerate)
        expr[g, ] <- rnorm(N, mesor_g[g], sigma_g[g])
      }
    } else {
      expr[g, ] <- rnorm(N, mesor_g[g], sigma_g[g])
    }
  }

  list(expr = expr, is_rhythmic = is_rhythmic, r_values = r_values,
       omega = omega,
       omega_g = if (n_rhythmic > 0) omega_g else NULL,
       alpha_g = if (n_rhythmic > 0) alpha_g else NULL)
}


#' Simulate Two-Group Differential Circadian Data under FMM Waveform
#'
#' @description
#' Extends \code{\link{simCircadianDiff}} by replacing the cosinor generative
#' model with the FMM (Frequency Modulated Möbius) waveform.  Gene-type
#' assignment (DR/DP/DM/null) and parameter draws are identical to
#' \code{simCircadianDiff}; only expression generation uses FMM.
#'
#' Uses a precomputed unit-amplitude FMM template (one per omega value) and
#' per-gene interpolation, avoiding the per-gene \code{FMM::generateFMM()}
#' call and giving a ~1000x speed-up over naive implementation.
#'
#' \strong{Beta-invariance:} The orientation parameter \code{beta} shifts the
#' peak location but does not alter amplitude.  DCP differential-amplitude
#' power (DR, DP) is therefore invariant to \code{beta}.  Only phase
#' estimation (DP peak-time accuracy) is affected by \code{beta}.
#'
#' \strong{omega = 1 identity:} When \code{omega = 1} the FMM waveform reduces
#' exactly to the cosinor model; this function returns the \code{simCircadianDiff}
#' result unmodified for efficiency.
#'
#' @references Rueda C, Rodríguez-Collado A, Peddada SD (2019).
#'   "A single notable oscillation in humans and mammals with strong
#'   implications in chronobiology." Sci Rep 9, 17982.
#'   \doi{10.1038/s41598-019-54569-1}
#'
#' @param omega Numeric in [0, 1]. FMM waveform shape: \code{1} = pure
#'   cosinor (default); \code{0} = flat (arrhythmic); values in \eqn{(0,1)}
#'   give increasingly peaked non-sinusoidal waveforms.
#' @param beta  Numeric. FMM orientation parameter (default \eqn{\pi} gives
#'   a symmetric peak at \code{alpha}).  Does not affect amplitude-based power
#'   (beta-invariance).
#' @param ...   All other arguments passed directly to \code{\link{simCircadianDiff}}.
#'
#' @return Same list structure as \code{\link{simCircadianDiff}}.
#'
#' @seealso \code{\link{simCircadianDiff}}, \code{\link{simCircadianFMM}},
#'   \code{\link{plotFMMDifferential}}
simCircadianDiffFMM <- function(..., omega = 1.0, beta = pi) {

  if (omega < 0 || omega > 1)
    stop("omega must be in [0, 1]; got ", omega)

  # Step 1: run the cosinor simulation to get parameter assignments and times
  sim <- simCircadianDiff(...)
  if (omega == 1.0) return(sim)   # omega=1 IS the cosinor — return as-is

  gt    <- sim$ground_truth
  ngenes <- nrow(gt)
  n1     <- ncol(sim$expr1)
  n2     <- ncol(sim$expr2)

  # Step 2: precompute unit-amplitude FMM template on [0, 2pi)
  # Y_template(x) = cos(beta + 2*atan(omega * tan(x/2))), x in [0, 2pi)
  # Singularity at x=pi is handled by R's atan(Inf) = pi/2 numerically.
  n_pts  <- 2000L
  x_grid <- seq(0, 2 * pi, length.out = n_pts + 1L)[seq_len(n_pts)]
  y_tmpl <- cos(beta + 2 * atan(omega * tan(x_grid / 2)))

  # Step 3: re-generate expression using FMM template via interpolation
  cts_rad1 <- (sim$times1 %% 24) * (2 * pi / 24)
  cts_rad2 <- (sim$times2 %% 24) * (2 * pi / 24)

  fmm_signal <- function(cts_rad, A, alpha_h, M, sigma) {
    if (A == 0 || omega == 0) return(rnorm(length(cts_rad), M, sigma))
    alpha_rad <- (alpha_h %% 24) * (2 * pi / 24)
    x_shifted <- (cts_rad - alpha_rad) %% (2 * pi)
    sig <- approx(x_grid, y_tmpl, xout = x_shifted, method = "linear", rule = 2)$y
    rnorm(length(cts_rad), M + A * sig, sigma)
  }

  expr1 <- matrix(NA_real_, nrow = ngenes, ncol = n1)
  expr2 <- matrix(NA_real_, nrow = ngenes, ncol = n2)

  for (g in seq_len(ngenes)) {
    expr1[g, ] <- fmm_signal(cts_rad1, gt$amplitude1[g], gt$phase1[g],
                              gt$mesor1[g], gt$sigma[g])
    expr2[g, ] <- fmm_signal(cts_rad2, gt$amplitude2[g], gt$phase2[g],
                              gt$mesor2[g], gt$sigma2[g])
  }

  sim$expr1 <- expr1
  sim$expr2 <- expr2
  sim
}
