# =====================================================================
# Built-in Pilot Data Loader
# =====================================================================

# Cache environment to avoid repeated disk reads
.pilot_cache <- new.env(parent = emptyenv())

#' Load built-in pilot dataset by name
#'
#' @param name Dataset name (e.g., "ba11_ba47_younger")
#' @return List with lBaselineExpr, lOD, amplitude, phase, prop_rhythmic
.load_pilot_data <- function(name) {
  # Return cached version if available
  if (exists(name, envir = .pilot_cache)) {
    return(get(name, envir = .pilot_cache))
  }

  # Find data/ directory relative to code/
  # options.R lives in code/, data/ is a sibling directory
  code_dir <- tryCatch(
    normalizePath(sys.frame(1)$ofile %||% ".", mustWork = FALSE),
    error = function(e) "."
  )
  # Try multiple resolution strategies
  candidates <- c(
    file.path(dirname(code_dir), "data", paste0(name, ".rds")),
    file.path(getwd(), "data", paste0(name, ".rds")),
    file.path(getwd(), "..", "data", paste0(name, ".rds"))
  )

  rds_path <- NULL
  for (cand in candidates) {
    if (file.exists(cand)) {
      rds_path <- normalizePath(cand)
      break
    }
  }

  if (is.null(rds_path)) {
    stop(sprintf(
      "Pilot dataset '%s' not found. Expected at data/%s.rds\n  Searched: %s\n  Run data/estimate_pilot_params.R to generate it.",
      name, name, paste(candidates, collapse = "\n           ")
    ))
  }

  pilot <- readRDS(rds_path)
  assign(name, pilot, envir = .pilot_cache)
  pilot
}


# Helper functions for setting parameters

setBaselineExpr <- function(input, ngenes) {
  if (is.numeric(input)) {
    if (length(input) == 1) {
      return(rep(input, ngenes))
    } else if (length(input) != ngenes) {
      return(sample(input, ngenes, replace = TRUE))
    } else {
      return(input)
    }
  } else if (is.function(input)) {
    return(input(ngenes))
  } else if (is.character(input)) {
    pilot <- .load_pilot_data(input)
    return(sample(pilot$lBaselineExpr, ngenes, replace = TRUE))
  } else {
    stop("Unrecognized form of lBaselineExpr")
  }
}

setOD <- function(input, ngenes) {
  if (is.numeric(input)) {
    if (length(input) == 1) {
      return(rep(input, ngenes))
    } else if (length(input) != ngenes) {
      return(sample(input, ngenes, replace = TRUE))
    } else {
      return(input)
    }
  } else if (is.function(input)) {
    return(input(ngenes))
  } else if (is.character(input)) {
    pilot <- .load_pilot_data(input)
    return(sample(pilot$lOD, ngenes, replace = TRUE))
  } else {
    stop("Unrecognized form of lOD")
  }
}

setAmplitude <- function(input, n_rhythmic) {
  if (is.numeric(input)) {
    if (length(input) == 1) {
      return(rep(input, n_rhythmic))
    } else if (length(input) != n_rhythmic) {
      return(sample(input, n_rhythmic, replace = TRUE))
    } else {
      return(input)
    }
  } else if (is.function(input)) {
    return(input(n_rhythmic))
  } else if (is.character(input)) {
    pilot <- .load_pilot_data(input)
    return(sample(pilot$amplitude, n_rhythmic, replace = TRUE))
  } else {
    stop("Unrecognized form of amplitude")
  }
}

setPhase <- function(input, n_rhythmic, period) {
  if (is.character(input) && input == "uniform") {
    return(runif(n_rhythmic, 0, period))
  } else if (is.character(input)) {
    pilot <- .load_pilot_data(input)
    return(sample(pilot$phase, n_rhythmic, replace = TRUE))
  } else if (is.numeric(input)) {
    if (length(input) == 1) {
      return(rep(input, n_rhythmic))
    } else if (length(input) != n_rhythmic) {
      return(sample(input, n_rhythmic, replace = TRUE))
    } else {
      return(input)
    }
  } else if (is.function(input)) {
    return(input(n_rhythmic))
  } else {
    stop("Unrecognized form of phase")
  }
}


# =====================================================================
# Null-coalescing operator
# =====================================================================
`%||%` <- function(a, b) if (!is.null(a)) a else b


# =====================================================================
# Configuration Constructors
# =====================================================================

#' Create Biology + Differential Options
#'
#' @param ngenes Number of genes.
#' @param prop_rhythmic Proportion of rhythmic genes. If \code{NULL} and a pilot
#'   dataset name is supplied to \code{lBaselineExpr}, the pilot's estimated
#'   proportion is used.
#' @param period Circadian period in hours.
#' @param lBaselineExpr Log baseline expression. Can be a numeric scalar (constant
#'   for all genes), a numeric vector (resampled to \code{ngenes}), a function
#'   \code{f(n)} returning \code{n} values, or a character string naming a built-in
#'   pilot dataset stored in \code{data/}. No default, must be supplied.
#' @param lBaselineExpr2 As \code{lBaselineExpr} but for group 2 (differential
#'   analyses only). \code{NULL} uses the same distribution as group 1.
#' @param lOD Log over-dispersion (noise). Same forms accepted as
#'   \code{lBaselineExpr}. No default, must be supplied.
#' @param lOD2 As \code{lOD} but for group 2. \code{NULL} shares group 1 values.
#' @param amplitude Amplitude distribution for rhythmic genes. Same forms
#'   accepted as \code{lBaselineExpr}. No default, must be supplied.
#' @param amplitude2 As \code{amplitude} but for group 2.
#' @param sigma_rhythmic Optional numeric vector of per-gene noise values for
#'   rhythmic genes (same length as \code{amplitude}). When provided, amplitude
#'   and sigma are drawn jointly to preserve pilot A-sigma correlation.
#' @param paired_sigma Logical. If \code{TRUE} and \code{sigma_rhythmic} is
#'   supplied, expand \code{amplitude} and \code{sigma_rhythmic} jointly by the
#'   same resampled gene index so their pilot-estimated pairwise correlation is
#'   preserved; if \code{FALSE} (default), amplitude and sigma are expanded
#'   independently.
#' @param paired_omega Logical. If \code{TRUE} and \code{omega_rhythmic} is
#'   supplied (same length as \code{amplitude}), expand omega using the same
#'   shared gene-index used for amplitude/sigma so all three are paired per gene.
#' @param paired_alpha Logical. Same as \code{paired_omega} for \code{alpha_rhythmic}.
#' @param omega_rhythmic Optional numeric vector of per-gene FMM omega values
#'   from a pilot FMM fit (same length as \code{amplitude}). Used when
#'   \code{paired_omega = TRUE}; ignored when \code{omega_dist} is also set.
#' @param alpha_rhythmic Optional numeric vector of per-gene FMM alpha values
#'   in radians from a pilot FMM fit (same length as \code{amplitude}).
#' @param omega_dist Optional list specifying a distribution from which omega
#'   is drawn at simulation time, e.g. \code{list(family="beta", a=1, b=5)} or
#'   \code{list(family="fixed", value=0.7)}. Overrides \code{omega_rhythmic}
#'   when both are provided.
#' @param alpha_dist Optional list specifying a distribution for alpha,
#'   \code{list(family="normal", mean=0, sd_hours=2)}. Internally
#'   \code{sd_hours} is converted to radians via \code{sd_rad = sd_hours * 2*pi/24}.
#'   When \code{alpha_rhythmic} is also provided and lengths match,
#'   the perturbation is added to each gene's empirical alpha
#'   (alpha_g = alpha_rhythmic[g] + N(0, sd_rad^2)).
#'   When \code{alpha_rhythmic} is missing, alpha_g = N(0, sd_rad^2).
#' @param cts Numeric vector of sample collection times for passive design.
#' @param cts2 As \code{cts} for group 2.
#' @param phase Phase distribution for rhythmic genes. \code{"uniform"} (default)
#'   or the same forms as \code{lBaselineExpr}.
#' @param prop_DR Proportion with differential rhythmicity.
#' @param prop_DP Proportion with differential phase.
#' @param prop_DM Proportion with differential mesor (mean shift, both groups rhythmic).
#' @param phase_diff Range of phase shift for DP genes \code{c(min, max)}.
#' @param amp_diff Range of amplitude ratio (unused; retained for interface compatibility).
#' @param mesor_diff Range of mesor shift for DM genes \code{c(min, max)}
#'   (additive, log-scale units).
#' @param dp_shift_mode \code{"fixed"} (use \code{phase_diff[2]}) or
#'   \code{"uniform"} (sample uniformly within \code{phase_diff} range).
#' @param dr_amp_scale Scale factor for amplitude (A) to adjust DR strength.
#' @param dr_sigma_scale Scale factor for sigma to adjust DR strength.
#' @param alpha_pilot Rhythmicity pre-screen threshold used to define the
#'   pilot's rhythmic gene set (the value at which \code{prop_rhythmic} and the
#'   effect-size distribution were calibrated). Stored on the object for
#'   provenance and reporting; default \code{NULL}.
#' @param sim.seed Random seed.
#'
#' @return Object of class \code{"CircadianBioOptions"}.
#' @export
CircadianBioOptions <- function(ngenes = 5000,
                                prop_rhythmic = NULL,
                                period = 24,
                                lBaselineExpr = NULL,
                                lBaselineExpr2 = NULL,
                                lOD = NULL,
                                lOD2 = NULL,
                                amplitude = NULL,
                                amplitude2 = NULL,
                                sigma_rhythmic = NULL,
                                paired_sigma = FALSE,
                                paired_omega = FALSE,
                                paired_alpha = FALSE,
                                omega_rhythmic = NULL,
                                alpha_rhythmic = NULL,
                                omega_dist = NULL,
                                alpha_dist = NULL,
                                cts = NULL,
                                cts2 = NULL,
                                phase = "uniform",
                                prop_DR = 0.15,
                                prop_DP = 0.10,
                                prop_DM = 0.00,
                                phase_diff = c(-6, 6),
                                amp_diff = c(0.5, 2),
                                mesor_diff = c(0.5, 2.0),
                                dp_shift_mode = c("fixed", "uniform"),
                                dr_amp_scale = 1.0,
                                dr_sigma_scale = 1.0,
                                alpha_pilot = NULL,
                                sim.seed = 12345) {

  dp_shift_mode <- match.arg(dp_shift_mode)

  # Require the three core distribution parameters, no hard-coded defaults.
  # Users must supply their own pilot data name (character string) or numeric
  # vector/scalar/function.  estCircadianParam() is the recommended entry point.
  if (is.null(lBaselineExpr))
    stop("lBaselineExpr is required. Supply a numeric vector, scalar, function, ",
         "or a character string naming a built-in pilot dataset.")
  if (is.null(lOD))
    stop("lOD is required. Supply a numeric vector, scalar, function, ",
         "or a character string naming a built-in pilot dataset.")
  if (is.null(amplitude))
    stop("amplitude is required. Supply a numeric vector, scalar, function, ",
         "or a character string naming a built-in pilot dataset.")

  # If prop_rhythmic is NULL, use pilot estimate
  if (is.null(prop_rhythmic)) {
    # Try to get from pilot data if any distribution arg is a pilot string
    pilot_name <- NULL
    for (arg in list(lBaselineExpr, lOD, amplitude)) {
      if (is.character(arg) && arg != "uniform") {
        pilot_name <- arg
        break
      }
    }
    if (!is.null(pilot_name)) {
      pilot <- .load_pilot_data(pilot_name)
      prop_rhythmic <- pilot$prop_rhythmic
    } else {
      prop_rhythmic <- 0.30  # fallback default
    }
  }

  # Validate proportions
  stopifnot(prop_rhythmic >= 0, prop_rhythmic <= 1)
  stopifnot(prop_DR >= 0, prop_DR <= 1)
  stopifnot(prop_DP >= 0, prop_DP <= 1)
  stopifnot(prop_DM >= 0, prop_DM <= 1)
  total_diff <- prop_DR + prop_DP + prop_DM
  if (total_diff > 1) stop("prop_DR + prop_DP + prop_DM must be <= 1")

  # All differential genes are rhythmic in at least one group, so prop_rhythmic
  # cannot be smaller than their combined fraction without contradicting the truth.
  rhythmic_required <- prop_DR + prop_DP + prop_DM
  if (prop_rhythmic < rhythmic_required - 1e-9) {
    stop(sprintf(
      "prop_rhythmic (%.3f) must be >= prop_DR + prop_DP + prop_DM (%.3f).\n%s",
      prop_rhythmic, rhythmic_required,
      "All DR/DP/DM genes are rhythmic in at least one group and count toward the rhythmic budget."
    ))
  }

  stopifnot(ngenes > 0, period > 0)
  stopifnot(length(phase_diff) == 2, length(amp_diff) == 2)
  stopifnot(dr_amp_scale > 0, dr_sigma_scale > 0)

  rng_state <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    .GlobalEnv$.Random.seed else NULL
  set.seed(sim.seed)
  on.exit(
    if (!is.null(rng_state))
      assign(".Random.seed", rng_state, envir = .GlobalEnv),
    add = TRUE
  )
  n_rhythmic <- round(ngenes * prop_rhythmic)

  # Resolve distributions using existing set*() helpers, store specs for re-seeding
  lBaselineExpr_resolved  <- setBaselineExpr(lBaselineExpr, ngenes)
  lBaselineExpr2_resolved <- if (!is.null(lBaselineExpr2)) setBaselineExpr(lBaselineExpr2, ngenes) else NULL
  lOD_resolved  <- setOD(lOD,  ngenes)
  lOD2_resolved <- if (!is.null(lOD2)) setOD(lOD2, ngenes) else NULL
  # paired_sigma=TRUE: expand amplitude and sigma_rhythmic jointly using a shared
  # index so each simulated rhythmic gene draws (A, sigma) from the same pilot gene,
  # preserving the empirical r_tilde = A/sigma distribution.
  # paired_sigma=FALSE (default): original behaviour, amplitude expanded independently
  # via setAmplitude, sigma drawn from lOD in simulation (has_joint=FALSE).
  #
  # If paired_omega = TRUE / paired_alpha = TRUE and the corresponding *_rhythmic
  # vector matches the amplitude length, we reuse the SAME ji index so that
  # (A, sigma, omega, alpha) for each simulated gene all come from the same pilot gene.
  amp_len_orig   <- if (is.numeric(amplitude)) length(amplitude) else NA_integer_
  paired_pairing <- !is.na(amp_len_orig) && amp_len_orig >= 2L

  # Compute the shared ji once if any pairing is requested with sufficient
  # amplitude length. Used for sigma, omega, and alpha jointly.
  any_pair_request <- (paired_sigma && !is.null(sigma_rhythmic) &&
                       is.numeric(amplitude) && length(amplitude) > 1L &&
                       length(amplitude) == length(sigma_rhythmic) &&
                       length(amplitude) != n_rhythmic) ||
                      (paired_omega && !is.null(omega_rhythmic) &&
                       length(omega_rhythmic) == amp_len_orig &&
                       length(omega_rhythmic) != n_rhythmic) ||
                      (paired_alpha && !is.null(alpha_rhythmic) &&
                       length(alpha_rhythmic) == amp_len_orig &&
                       length(alpha_rhythmic) != n_rhythmic)
  ji <- if (any_pair_request && paired_pairing)
          sample(amp_len_orig, n_rhythmic, replace = TRUE) else NULL

  if (paired_sigma &&
      !is.null(sigma_rhythmic) &&
      is.numeric(amplitude) && length(amplitude) > 1L &&
      length(amplitude) == length(sigma_rhythmic) &&
      length(amplitude) != n_rhythmic) {
    amplitude_resolved <- amplitude[ji]
    sigma_rhythmic     <- sigma_rhythmic[ji]
  } else {
    amplitude_resolved <- setAmplitude(amplitude, n_rhythmic)
    # Independent fallback resample when sigma_rhythmic is supplied without
    # pairing (or with a length mismatch the paired branch couldn't handle).
    # Without this, a length(sigma_rhythmic) != n_rhythmic value would silently
    # propagate to the simulator and either crash or recycle incorrectly.
    if (!is.null(sigma_rhythmic) && length(sigma_rhythmic) != n_rhythmic) {
      sigma_rhythmic <- if (length(sigma_rhythmic) == 1L)
        rep(sigma_rhythmic, n_rhythmic)
      else
        sample(sigma_rhythmic, n_rhythmic, replace = TRUE)
    }
  }

  # Expand omega_rhythmic with the shared ji when paired_omega requested.
  # If the paired expansion cannot proceed (length mismatch or amplitude is
  # not a numeric vector), warn so the caller knows pairing is silently off.
  if (paired_omega) {
    if (!is.null(omega_rhythmic) &&
        length(omega_rhythmic) == amp_len_orig &&
        length(omega_rhythmic) != n_rhythmic &&
        !is.null(ji)) {
      omega_rhythmic <- omega_rhythmic[ji]
    } else if (!is.null(omega_rhythmic) &&
               length(omega_rhythmic) != n_rhythmic) {
      warning(sprintf(
        "paired_omega=TRUE but omega_rhythmic length (%d) does not match amplitude length (%s) or n_rhythmic (%d); pairing skipped.",
        length(omega_rhythmic),
        if (is.na(amp_len_orig)) "NA (amplitude is non-numeric)" else as.character(amp_len_orig),
        n_rhythmic))
    }
  }

  # Same for alpha.
  if (paired_alpha) {
    if (!is.null(alpha_rhythmic) &&
        length(alpha_rhythmic) == amp_len_orig &&
        length(alpha_rhythmic) != n_rhythmic &&
        !is.null(ji)) {
      alpha_rhythmic <- alpha_rhythmic[ji]
    } else if (!is.null(alpha_rhythmic) &&
               length(alpha_rhythmic) != n_rhythmic) {
      warning(sprintf(
        "paired_alpha=TRUE but alpha_rhythmic length (%d) does not match amplitude length (%s) or n_rhythmic (%d); pairing skipped.",
        length(alpha_rhythmic),
        if (is.na(amp_len_orig)) "NA (amplitude is non-numeric)" else as.character(amp_len_orig),
        n_rhythmic))
    }
  }

  amplitude2_resolved <- if (!is.null(amplitude2)) setAmplitude(amplitude2, n_rhythmic) else NULL
  phase_resolved <- setPhase(phase, n_rhythmic, period)

  opts <- list(
    ngenes = ngenes,
    prop_rhythmic = prop_rhythmic,
    period = period,
    lBaselineExpr = lBaselineExpr_resolved,
    lBaselineExpr_spec = lBaselineExpr,
    lBaselineExpr2 = lBaselineExpr2_resolved,
    lOD = lOD_resolved,
    lOD_spec = lOD,
    lOD2 = lOD2_resolved,
    amplitude = amplitude_resolved,
    amplitude_spec = amplitude,
    amplitude2 = amplitude2_resolved,
    sigma_rhythmic = sigma_rhythmic,
    paired_omega   = paired_omega,
    paired_alpha   = paired_alpha,
    omega_rhythmic = omega_rhythmic,
    alpha_rhythmic = alpha_rhythmic,
    omega_dist     = omega_dist,
    alpha_dist     = alpha_dist,
    cts = cts,
    cts2 = cts2,
    phase = phase_resolved,
    phase_spec = phase,
    prop_DR = prop_DR,
    prop_DP = prop_DP,
    prop_DM = prop_DM,
    phase_diff = phase_diff,
    amp_diff = amp_diff,
    mesor_diff = mesor_diff,
    dp_shift_mode = dp_shift_mode,
    dr_amp_scale = dr_amp_scale,
    dr_sigma_scale = dr_sigma_scale,
    alpha_pilot = alpha_pilot,
    sim.seed = sim.seed
  )
  class(opts) <- "CircadianBioOptions"
  opts
}


#' Create Study Design Options
#'
#' @param sample_sizes Vector of sample sizes per group
#' @param nsims Number of simulation replicates
#' @param design "active" or "passive"
#' @param cts Time-of-day distribution for passive design (numeric vector)
#' @param B_values Optional vector of time-point counts (B) to sweep, in
#'   addition to \code{sample_sizes}; NULL runs a single B implied by \code{cts}.
#' @param test_types Character vector of tests to run ("DR", "DP", "DM")
#' @param omega FMM waveform shape parameter in (0, 1]. \code{omega = 1}
#'   (default) is the pure cosinor path; \code{omega < 1} simulates a
#'   non-sinusoidal waveform via \code{simCircadianFMM()}. \code{omega = 0} is
#'   not allowed.
#' @param beta FMM orientation parameter (peak-location offset, default
#'   \code{pi}); only used when \code{omega < 1}, ignored for the cosinor path.
#'
#' @return Object of class "CircadianDesignOptions"
#' @export
CircadianDesignOptions <- function(sample_sizes = c(10, 20, 40, 60, 80, 100),
                                   nsims      = 100,
                                   design     = c("active", "passive"),
                                   cts        = NULL,
                                   B_values   = NULL,
                                   test_types = c("DR", "DP", "DM"),
                                   omega      = 1.0,
                                   beta       = pi) {
  # omega: FMM waveform shape parameter.
  #   omega = 1  (default) -> pure cosinor simulation (traditional path).
  #   omega < 1  -> FMM non-sinusoidal simulation via simCircadianFMM().
  #   omega = 0  is not allowed (degenerate flat waveform).
  # beta: FMM orientation parameter (peak location offset, default pi).
  #   Only used when omega < 1; ignored for cosinor simulation.

  design <- match.arg(design)
  stopifnot(all(sample_sizes > 0), nsims > 0)
  if (design == "passive" && is.null(cts))
    stop("cts (time-of-day vector) is required for passive design")
  if (omega <= 0 || omega > 1)
    stop("omega must be in (0, 1]. Use omega = 1 (default) for the cosinor path.")

  opts <- list(
    sample_sizes = sample_sizes,
    nsims        = nsims,
    design       = design,
    cts          = cts,
    B_values     = B_values,
    test_types   = test_types,
    omega        = omega,
    beta         = beta
  )
  class(opts) <- "CircadianDesignOptions"
  opts
}


# Helper to auto-generate strata labels from breakpoints
.make_strata_labels <- function(r_strata) {
  n <- length(r_strata) - 1
  labels <- character(n)
  for (i in seq_len(n)) {
    lo <- r_strata[i]
    hi <- r_strata[i + 1]
    if (is.infinite(hi)) {
      labels[i] <- sprintf(">%g", lo)
    } else {
      labels[i] <- sprintf("(%g,%g]", lo, hi)
    }
  }
  labels
}


#' Create Analysis & Reporting Options
#'
#' @param alpha Significance level for DCP pipeline
#' @param p.adjust.method Multiple testing correction method
#' @param parallel.ncores Number of cores for parallel DCP
#' @param amp.cutoff Amplitude cutoff for DCP_Rhythmicity
#' @param target_effect Minimum effect size to count as "interesting"
#' @param fdr_thresholds FDR thresholds for power curves
#' @param reference_n Reference sample size for Panel B / phase shift plots
#' @param r_strata Breakpoints for A/sigma stratification
#' @param strata_labels Labels for strata (auto-generated if NULL)
#' @param phase_shifts Phase shift magnitudes to sweep (for sensitivity analysis)
#' @param DCmethod Differential-rhythmicity detector to report against,
#'   one of "DCP" (default) or "CircaCompare".
#'
#' @return Object of class "CircadianAnalysisOptions"
#' @export
CircadianAnalysisOptions <- function(alpha = 0.05,
                                     p.adjust.method = "BH",
                                     parallel.ncores = 1,
                                     amp.cutoff = 0,
                                     target_effect = 0.1,
                                     fdr_thresholds = c(0.01, 0.05, 0.10, 0.20),
                                     reference_n = 60,
                                     r_strata = c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5,
                                                  1.75, 2, 2.5, 3, 3.5, 4, 4.5, 5, Inf),
                                     strata_labels = NULL,
                                     phase_shifts = c(0, 0.5, 1, 2, 4, 6, 8, 10, 12),
                                     DCmethod = "DCP") {

  DCmethod <- match.arg(DCmethod, c("DCP", "CircaCompare"))

  stopifnot(alpha > 0, alpha < 1)
  stopifnot(all(fdr_thresholds > 0), all(fdr_thresholds < 1))
  stopifnot(reference_n > 0)

  if (is.null(strata_labels)) {
    strata_labels <- .make_strata_labels(r_strata)
  } else {
    stopifnot(length(strata_labels) == length(r_strata) - 1)
  }

  opts <- list(
    alpha = alpha,
    p.adjust.method = p.adjust.method,
    parallel.ncores = parallel.ncores,
    amp.cutoff = amp.cutoff,
    target_effect = target_effect,
    fdr_thresholds = fdr_thresholds,
    reference_n = reference_n,
    r_strata = r_strata,
    strata_labels = strata_labels,
    phase_shifts = phase_shifts,
    DCmethod = DCmethod
  )
  class(opts) <- "CircadianAnalysisOptions"
  opts
}


#' Create Bootstrap Design Grid Options
#'
#' @param design_vector Ordered time points to sweep over (length >= max(B_values))
#' @param B_values Number of distinct time points to test
#' @param m_values Replicates per time point (if NULL, derived as round(N/B))
#' @param N_values Total per-group sample sizes to test (if NULL, derived as B*m)
#' @param nboot Number of bootstrap parameter draws (outer uncertainty loop)
#' @param nsims_inner Simulations per bootstrap draw
#' @param design "active" or "passive"
#' @param seed Random seed
#'
#' @return Object of class "CircadianBootstrapOptions"
#' @export
CircadianBootstrapOptions <- function(design_vector,
                                      B_values     = c(4, 6, 8, 12, 24),
                                      m_values     = NULL,
                                      N_values     = NULL,
                                      nboot        = 200,
                                      nsims_inner  = 20,
                                      design       = c("active", "passive"),
                                      seed         = 42) {
  design <- match.arg(design)

  if (missing(design_vector) || is.null(design_vector)) {
    stop("design_vector is required for CircadianBootstrapOptions")
  }

  if (design == "active" && length(design_vector) < max(B_values)) {
    stop(sprintf(
      "For active design, length(design_vector) >= max(B_values). Got %d < %d",
      length(design_vector), max(B_values)
    ))
  }

  stopifnot(nboot > 0, nsims_inner > 0)
  stopifnot(all(B_values > 0))

  # Cosinor model requires >= 3 distinct time points for identifiable fitting
  # (2 parameters: amplitude + phase, plus intercept = minimum 3 time points).
  if (any(B_values < 3)) {
    stop(sprintf(
      "All B_values must be >= 3 (cosinor needs at least 3 time points for identifiable fitting). Got: %s",
      paste(B_values[B_values < 3], collapse = ", ")
    ))
  }

  # Passive design does not vary by B: the TOD distribution (design_vector) is
  # fixed regardless of B, so all B values would produce identical simulations.
  # A B-vs-m grid is only identified under active design.
  if (design == "passive" && length(B_values) > 1) {
    stop(paste(
      "Passive design cannot identify a B effect: subject times are sampled from",
      "design_vector regardless of B, so all B values produce identical results.",
      "\n  For passive mode: set B_values to a single placeholder value (e.g., B_values = 4)",
      "and interpret power_mean as a function of N only."
    ))
  }

  # Derive N_values from B_values and m_values, or vice versa
  if (is.null(N_values) && !is.null(m_values)) {
    N_values <- sort(unique(as.integer(outer(B_values, m_values, FUN = "*"))))
  }
  if (is.null(N_values)) {
    stop("Provide N_values (or m_values to derive N_values = B * m)")
  }

  opts <- list(
    design_vector = design_vector,
    B_values      = sort(B_values),
    m_values      = m_values,
    N_values      = sort(N_values),
    nboot         = nboot,
    nsims_inner   = nsims_inner,
    design        = design,
    seed          = seed
  )
  class(opts) <- "CircadianBootstrapOptions"
  opts
}


# =====================================================================
# S3 Print Methods
# =====================================================================

#' Print a CircadianBioOptions object
#'
#' @param x A \code{CircadianBioOptions} object.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @exportS3Method
print.CircadianBioOptions <- function(x, ...) {
  cat("CircadianBioOptions\n")
  cat(sprintf("  ngenes:         %d\n", x$ngenes))
  if (!is.null(x$alpha_pilot)) {
    adj <- if (identical(x$adjust_pilot, "BH")) "BH-FDR q" else "raw p"
    cat(sprintf("  prop_rhythmic:  %.0f%% (alpha_pilot = %g, %s)\n",
                100 * x$prop_rhythmic, x$alpha_pilot, adj))
  } else {
    cat(sprintf("  prop_rhythmic:  %.0f%%\n", 100 * x$prop_rhythmic))
  }
  cat(sprintf("  period:         %g h\n", x$period))
  cat(sprintf("  prop_DR:        %.0f%%\n", 100 * x$prop_DR))
  cat(sprintf("  prop_DP:        %.0f%%\n", 100 * x$prop_DP))
  cat(sprintf("  prop_DM:        %.0f%%\n", 100 * (x$prop_DM %||% 0)))
  cat(sprintf("  phase_diff:     [%g, %g] h\n", x$phase_diff[1], x$phase_diff[2]))
  cat(sprintf("  amp_diff:       [%g, %g]\n", x$amp_diff[1], x$amp_diff[2]))
  cat(sprintf("  dp_shift_mode:  %s\n", x$dp_shift_mode))
  cat(sprintf("  dr_amp_scale:   %g\n", x$dr_amp_scale))
  cat(sprintf("  dr_sigma_scale: %g\n", x$dr_sigma_scale))
  cat(sprintf("  sim.seed:       %d\n", x$sim.seed))
  # Show distribution type for key params
  .print_spec <- function(name, spec) {
    # Map pilot dataset names to friendly labels
    .pilot_labels <- c(
      ba11_ba47_younger = "empirical (BA11+BA47 younger)"
    )
    if (is.function(spec)) {
      cat(sprintf("  %-15s function\n", paste0(name, ":")))
    } else if (is.character(spec) && spec %in% names(.pilot_labels)) {
      cat(sprintf("  %-15s %s\n", paste0(name, ":"), .pilot_labels[spec]))
    } else if (is.character(spec)) {
      cat(sprintf("  %-15s \"%s\"\n", paste0(name, ":"), spec))
    } else if (is.numeric(spec) && length(spec) == 1) {
      cat(sprintf("  %-15s %g (constant)\n", paste0(name, ":"), spec))
    } else if (is.numeric(spec)) {
      cat(sprintf("  %-15s vector[%d]\n", paste0(name, ":"), length(spec)))
    }
  }
  .print_spec("lBaselineExpr", x$lBaselineExpr_spec)
  .print_spec("lOD", x$lOD_spec)
  .print_spec("amplitude", x$amplitude_spec)
  .print_spec("phase", x$phase_spec)
  invisible(x)
}

#' Print a CircadianDesignOptions object
#'
#' @param x A \code{CircadianDesignOptions} object.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @exportS3Method
print.CircadianDesignOptions <- function(x, ...) {
  cat("CircadianDesignOptions\n")
  cat(sprintf("  sample_sizes:   %s\n", paste(x$sample_sizes, collapse = ", ")))
  cat(sprintf("  nsims:          %d\n", x$nsims))
  cat(sprintf("  design:         %s\n", x$design))
  cat(sprintf("  test_types:     %s\n", paste(x$test_types, collapse = ", ")))
  if (!is.null(x$B_values)) {
    cat(sprintf("  B_values:       %s\n", paste(x$B_values, collapse = ", ")))
  }
  if (!is.null(x$cts)) {
    cat(sprintf("  cts:            %d time points\n", length(x$cts)))
  }
  omega <- x$omega %||% 1.0
  if (omega < 1.0) {
    cat(sprintf("  sim model:      FMM  (omega = %.2f, beta = %.4f)\n",
                omega, x$beta %||% pi))
  } else {
    cat("  sim model:      cosinor  (omega = 1, traditional)\n")
  }
  invisible(x)
}

#' Print a CircadianAnalysisOptions object
#'
#' @param x A \code{CircadianAnalysisOptions} object.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @exportS3Method
print.CircadianAnalysisOptions <- function(x, ...) {
  cat("CircadianAnalysisOptions\n")
  cat(sprintf("  alpha:          %g\n", x$alpha))
  cat(sprintf("  p.adjust:       %s\n", x$p.adjust.method))
  cat(sprintf("  parallel.ncores:%d\n", x$parallel.ncores))
  cat(sprintf("  amp.cutoff:     %g\n", x$amp.cutoff))
  cat(sprintf("  target_effect:  %g\n", x$target_effect))
  cat(sprintf("  fdr_thresholds: %s\n", paste(x$fdr_thresholds, collapse = ", ")))
  cat(sprintf("  reference_n:    %d\n", x$reference_n))
  cat(sprintf("  r_strata:       %d strata\n", length(x$r_strata) - 1))
  cat(sprintf("  phase_shifts:   %s h\n", paste(x$phase_shifts, collapse = ", ")))
  cat(sprintf("  DCmethod:       %s\n", x$DCmethod))
  invisible(x)
}

#' Print a CircadianBootstrapOptions object
#'
#' @param x A \code{CircadianBootstrapOptions} object.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @exportS3Method
print.CircadianBootstrapOptions <- function(x, ...) {
  cat("CircadianBootstrapOptions\n")
  cat(sprintf("  design:         %s\n", x$design))
  cat(sprintf("  design_vector:  %d time points [%g to %g]\n",
              length(x$design_vector), min(x$design_vector), max(x$design_vector)))
  cat(sprintf("  B_values:       %s\n", paste(x$B_values, collapse = ", ")))
  cat(sprintf("  N_values:       %s\n", paste(x$N_values, collapse = ", ")))
  cat(sprintf("  nboot:          %d\n", x$nboot))
  cat(sprintf("  nsims_inner:    %d\n", x$nsims_inner))
  cat(sprintf("  seed:           %d\n", x$seed))
  invisible(x)
}
