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

#' Create Biology and Differential Options
#'
#' @description
#' Builds the \code{CircadianBioOptions} object that describes the simulated
#' transcriptome: the number of genes, the proportion that are rhythmic, and
#' the gene-level distributions of baseline expression, noise, amplitude, and
#' phase. For two-group studies it also records the proportions of genes that
#' differ in rhythmicity, phase, or mesor between groups. Most users obtain
#' this object indirectly from \code{\link{estCircadianParam}} or
#' \code{\link{scp_load_pilot}}; the constructor is exposed for building a
#' configuration by hand.
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
#' @param lOD Log residual noise: \code{log(sigma)}, the per-gene noise around
#'   the cosinor fit (\code{sigma = exp(lOD)}). Same forms accepted as
#'   \code{lBaselineExpr}. No default, must be supplied.
#' @param lOD2 As \code{lOD} (the residual noise) but for group 2. \code{NULL}
#'   shares group 1 values.
#' @param amplitude Amplitude distribution for rhythmic genes, given in any
#'   of the forms accepted for \code{lBaselineExpr}. This argument is
#'   required and has no default.
#' @param amplitude2 As \code{amplitude} but for group 2.
#' @param sigma_rhythmic Optional numeric vector of per-gene noise values for
#'   rhythmic genes (same length as \code{amplitude}). When provided, amplitude
#'   and sigma are drawn jointly to preserve pilot A-sigma correlation.
#' @param paired_sigma Logical. If \code{TRUE} and \code{sigma_rhythmic} is
#'   supplied, expand \code{amplitude} and \code{sigma_rhythmic} jointly by the
#'   same resampled gene index so their pilot-estimated pairwise correlation is
#'   preserved; if \code{FALSE} (default), amplitude and sigma are expanded
#'   independently.
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
#' @examples
#' # Build a synthetic biology configuration without pilot data.
#' bio <- CircadianBioOptions(ngenes = 300L, prop_rhythmic = 0.3,
#'                            lBaselineExpr = 5, lOD = -1, amplitude = 1)
#' bio
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
  # preserving the empirical r = A/sigma distribution.
  # paired_sigma=FALSE (default): original behaviour, amplitude expanded independently
  # via setAmplitude, sigma drawn from lOD in simulation (has_joint=FALSE).
  if (paired_sigma &&
      !is.null(sigma_rhythmic) &&
      is.numeric(amplitude) && length(amplitude) > 1L &&
      length(amplitude) == length(sigma_rhythmic) &&
      length(amplitude) != n_rhythmic) {
    ji                 <- sample(length(amplitude), n_rhythmic, replace = TRUE)
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
#' @description
#' Collects the sampling-design choices for a power analysis: the sample
#' sizes to evaluate, the number of simulation replicates per size, and
#' whether collection times are placed on an even grid (active design) or
#' drawn from the pilot's time-of-day distribution (passive design). The
#' resulting \code{CircadianDesignOptions} object is passed to the
#' \code{runSims*} and \code{recommendDesign} functions.
#'
#' @param sample_sizes Numeric vector of sample sizes to evaluate.
#' @param nsims Number of simulation replicates run at each sample size.
#' @param design Sampling design: \code{"active"} places collection times on
#'   an even grid, \code{"passive"} draws them from \code{cts}.
#' @param cts Numeric vector giving the time-of-day distribution used for a
#'   passive design.
#' @param B_values Optional vector of time-point counts (B) to sweep, in
#'   addition to \code{sample_sizes}. \code{NULL} runs the single B implied
#'   by \code{cts}.
#' @param test_types Character vector of differential tests to run, any of
#'   \code{"DR"}, \code{"DP"}, and \code{"DM"}.
#'
#' @return Object of class "CircadianDesignOptions"
#' @examples
#' # Active design over two small sample sizes with few simulations.
#' dopt <- CircadianDesignOptions(sample_sizes = c(24, 48), nsims = 5)
#' dopt
#' @export
CircadianDesignOptions <- function(sample_sizes = c(10, 20, 40, 60, 80, 100),
                                   nsims      = 100,
                                   design     = c("active", "passive"),
                                   cts        = NULL,
                                   B_values   = NULL,
                                   test_types = c("DR", "DP", "DM")) {

  design <- match.arg(design)
  stopifnot(all(sample_sizes > 0), nsims > 0)
  if (design == "passive" && is.null(cts))
    stop("cts (time-of-day vector) is required for passive design")

  opts <- list(
    sample_sizes = sample_sizes,
    nsims        = nsims,
    design       = design,
    cts          = cts,
    B_values     = B_values,
    test_types   = test_types
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


#' Create Analysis and Reporting Options
#'
#' @description
#' Collects the settings that govern how each simulated dataset is analysed
#' and how power is reported: the per-gene significance level, the
#' multiple-testing correction, the FDR thresholds at which power curves are
#' drawn, and the effect-size strata used to break power down by signal
#' strength. The resulting \code{CircadianAnalysisOptions} object is passed
#' to the \code{runSims*} and plotting functions.
#'
#' @usage
#' CircadianAnalysisOptions(
#'   alpha = 0.05, p.adjust.method = "BH", parallel.ncores = 1,
#'   amp.cutoff = 0, target_effect = 0.1,
#'   fdr_thresholds = c(0.01, 0.05, 0.10, 0.20), reference_n = 60,
#'   r_strata = c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75,
#'                2, 2.5, 3, 3.5, 4, 4.5, 5, Inf),
#'   strata_labels = NULL,
#'   phase_shifts = c(0, 0.5, 1, 2, 4, 6, 8, 10, 12),
#'   DCmethod = "DCP")
#'
#' @param alpha Per-gene significance level for the cosinor rhythmicity test.
#' @param p.adjust.method Multiple-testing correction applied across genes,
#'   passed to \code{stats::p.adjust} (default Benjamini-Hochberg).
#' @param parallel.ncores Number of cores used to fit the per-gene tests.
#' @param amp.cutoff Minimum amplitude a gene must exceed to enter the
#'   cosinor rhythmicity screen.
#' @param target_effect Minimum effect size a gene must reach to count as a
#'   meaningful discovery.
#' @param fdr_thresholds FDR thresholds at which power curves are reported.
#' @param reference_n Reference sample size marked on the panel-B and
#'   phase-shift plots.
#' @param r_strata Breakpoints, in \eqn{r = A / \sigma} units, that define
#'   the effect-size strata used to report power by signal strength.
#' @param strata_labels Labels for the strata; generated automatically from
#'   \code{r_strata} when \code{NULL}.
#' @param phase_shifts Phase-shift magnitudes (hours) swept in the
#'   phase-sensitivity analysis.
#' @param DCmethod Differential-rhythmicity detector to report against,
#'   either \code{"DCP"} (the cosinor F-test, default) or
#'   \code{"CircaCompare"}. \code{"DCP"} is a retained legacy alias for the
#'   cosinor detector.
#'
#' @return Object of class "CircadianAnalysisOptions"
#' @examples
#' # Reporting options: 5 percent alpha, power at two FDR thresholds.
#' aopt <- CircadianAnalysisOptions(alpha = 0.05,
#'                                  fdr_thresholds = c(0.05, 0.20))
#' aopt
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


#' Create Bootstrap Grid Options
#'
#' @description
#' Collects the settings for a bootstrapped design-grid analysis, which
#' quantifies how much the finite pilot contributes to uncertainty in
#' power across combinations of time points (B) and replicates (m). It
#' defines the grid of designs to evaluate and the number of bootstrap
#' draws and inner simulations used to form the confidence bands. The
#' resulting object is passed to \code{\link{runBootstrapDesignGrid}}.
#'
#' @param design_vector Ordered vector of candidate time points to sample
#'   from; its length must be at least \code{max(B_values)}.
#' @param B_values Numbers of distinct time points to evaluate.
#' @param m_values Replicates per time point. When \code{NULL}, derived as
#'   \code{round(N / B)}.
#' @param N_values Total per-group sample sizes to evaluate. When
#'   \code{NULL}, derived as \code{B * m}.
#' @param nboot Number of bootstrap parameter draws in the outer
#'   uncertainty loop.
#' @param nsims_inner Number of simulations run per bootstrap draw.
#' @param design Sampling design, either \code{"active"} or \code{"passive"}.
#' @param seed Random seed.
#'
#' @return Object of class "CircadianBootstrapOptions"
#' @examples
#' # A small active-design B-vs-m grid over two sample sizes.
#' bopt <- CircadianBootstrapOptions(design_vector = seq(0, 22, by = 2),
#'                                   B_values = 4, N_values = c(24, 48),
#'                                   nboot = 3, nsims_inner = 4)
#' bopt
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
