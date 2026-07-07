#' Bootstrap Design Grid for Circadian Power Analysis
#'
#' PowerSim performs pilot-informed power analysis for cross-sectional circadian
#' studies. Gene-level circadian parameter distributions are estimated from pilot
#' data, and uncertainty in those distributions is propagated through an outer
#' bootstrap. Conditional on a chosen parameter distribution, independent subjects
#' are simulated under a population-average cosinor model for candidate designs.
#'
#' The bootstrap quantifies uncertainty in the pilot-implied parameter distribution
#' (mesor, amplitude, phase, noise). It does NOT create subject-level variability
#' within a single simulation; that is captured by sigma_g (residual SD from the
#' cosinor fit), which represents pooled between-subject spread conditional on TOD.
#' sigma_g does NOT model subject-specific random effects in M, A, or phi, the
#' model is a population-average cosinor, which is a defensible first-order
#' baseline for cross-sectional post-mortem studies.
#'
#' Key notation:
#'   B = number of distinct time bins (subset of design_vector; must be >= 3)
#'   m = subjects per time bin  -> N = B * m per-group sample size (exact, no rounding)
#'   nboot = outer bootstrap draws for pilot-estimation uncertainty


# =====================================================================
# Pilot Parameter Fitting
# =====================================================================

#' Fit cosinor model to all genes in pilot data
#'
#' @param data Genes x samples matrix
#' @param times Sample time points (length = ncol(data))
#' @param period Circadian period in hours (default 24)
#' @param min_rhythm_pval P-value threshold for classifying rhythmic genes
#'
#' @return Data frame: gene, M, A, phi, sigma, pvalue, is_rhythmic
#' @export
fitCosinorAll <- function(data, times, period = 24, min_rhythm_pval = 0.01) {
  if (exists(".CPP_LOADED", inherits = TRUE) && isTRUE(get(".CPP_LOADED", inherits = TRUE)) &&
      exists("fitCosinorAll_fast", mode = "function"))
    return(fitCosinorAll_fast(data, times, period = period, min_rhythm_pval = min_rhythm_pval))

  G <- nrow(data)
  omega <- 2 * pi / period

  results <- lapply(seq_len(G), function(g) {
    y <- data[g, ]
    tryCatch({
      fit   <- one_cosinor_OLS(times, y, period, compute.phase.CI = FALSE)
      n_obs <- sum(!is.na(y))
      if (n_obs <= 3L || is.na(fit$pvalue))
        return(list(gene=g, M=fit$M, A=fit$A, phi=fit$phi,
                    sigma=NA_real_, pvalue=fit$pvalue, r=NA_real_))
      yhat      <- fit$M + fit$A * cos(omega * times - omega * fit$phi)
      sigma_hat <- sqrt(sum((y - yhat)^2, na.rm = TRUE) / (n_obs - 3L))
      list(gene = g, M = fit$M, A = fit$A, phi = fit$phi, sigma = sigma_hat,
           pvalue = fit$pvalue, r = fit$A / max(sigma_hat, 1e-6))
    }, error = function(e) {
      list(gene = g, M = NA, A = NA, phi = NA, sigma = NA, pvalue = NA, r = NA)
    })
  })

  df <- do.call(rbind, lapply(results, as.data.frame))
  df$is_rhythmic <- !is.na(df$pvalue) & df$pvalue < min_rhythm_pval & !is.na(df$A) & df$A > 0
  df
}


# =====================================================================
# Bootstrap Parameter Resampling
# =====================================================================

#' Bootstrap parameter sets from pilot fit
#'
#' Resamples G rows (genes) with replacement from param_df to generate
#' nboot independent parameter draws.
#'
#' @param param_df Data frame from fitCosinorAll() (G rows)
#' @param nboot Number of bootstrap draws
#' @param seed Random seed
#'
#' @return List of length nboot; each element is a data frame like param_df
bootstrapParams <- function(param_df, nboot, seed = 42) {
  rng_state <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    .GlobalEnv$.Random.seed else NULL
  set.seed(seed)
  on.exit(
    if (!is.null(rng_state))
      assign(".Random.seed", rng_state, envir = .GlobalEnv),
    add = TRUE
  )
  G <- nrow(param_df)
  lapply(seq_len(nboot), function(b) {
    idx <- sample(G, G, replace = TRUE)
    param_df[idx, ]
  })
}


#' Reproducible subject-resample indices for the subject bootstrap
#'
#' Draws \code{nboot} sets of \code{n_pilot} subject indices with replacement
#' (Efron nonparametric bootstrap over subjects). The RNG state is saved and
#' restored so the draws are reproducible from \code{seed} without disturbing
#' the caller's random stream.
#'
#' @param n_pilot Number of pilot subjects (columns of the pilot matrix).
#' @param nboot Number of bootstrap draws.
#' @param seed Random seed.
#'
#' @return List of length \code{nboot}; each element is an integer index vector
#'   of length \code{n_pilot}.
.resampleSubjectIndices <- function(n_pilot, nboot, seed = 42) {
  rng_state <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    .GlobalEnv$.Random.seed else NULL
  set.seed(seed)
  on.exit(
    if (!is.null(rng_state))
      assign(".Random.seed", rng_state, envir = .GlobalEnv),
    add = TRUE
  )
  lapply(seq_len(nboot), function(b) sample(n_pilot, n_pilot, replace = TRUE))
}


# =====================================================================
# Internal helper: build CircadianBioOptions from bootstrap draw
# =====================================================================

#' Build CircadianBioOptions from a Bootstrap Parameter Draw
#'
#' @description
#' Constructs a \code{CircadianBioOptions} object from a single bootstrap
#' resample of pilot-gene cosinor fits. Differential parameters
#' (prop_DR, prop_DP, phase_diff, etc.) are inherited from \code{bio_diff.opts}.
#' Amplitude and sigma are kept jointly sampled to preserve the empirical
#' A-sigma correlation.
#'
#' @param boot_df Data frame. One bootstrap resample of \code{fitCosinorAll()}
#'   output (columns: M, A, phi, sigma, is_rhythmic).
#' @param bio_diff.opts List or \code{CircadianBioOptions}. Source for
#'   differential-design parameters (prop_DR, prop_DP, phase_diff, etc.).
#'
#' @return A \code{CircadianBioOptions} object ready for \code{simCircadianSingleCohort()}.
.buildBioFromBoot <- function(boot_df, bio_diff.opts) {
  # Use the bootstrap draw for base distributions; differential params from bio_diff.opts
  lBaselineExpr_vec <- boot_df$M
  lOD_vec           <- log(pmax(boot_df$sigma, 1e-6))
  rhythmic_mask     <- !is.na(boot_df$is_rhythmic) & boot_df$is_rhythmic
  prop_rhythmic     <- mean(rhythmic_mask, na.rm = TRUE)
  # Guard against degenerate bootstrap draws
  if (prop_rhythmic < 0.01) prop_rhythmic <- 0.01
  if (prop_rhythmic > 0.99) prop_rhythmic <- 0.99

  # Keep A and sigma paired from the same gene (joint sampling preserves A-sigma correlation).
  # Apply the same top-K=min(300, n_rhythmic) selection as estCircadianParam so that the
  # bootstrap bio mirrors the plug-in's pipeline exactly: only the top-300 rhythmic genes
  # ranked by p-value are used for the amplitude/sigma distribution.
  rhythmic_valid <- rhythmic_mask & !is.na(boot_df$A) & boot_df$A > 0 &
                    !is.na(boot_df$sigma) & boot_df$sigma > 0 & !is.na(boot_df$pvalue)
  if (sum(rhythmic_valid) > 300L) {
    # Rank by p-value ascending; keep top-300 (strongest signal, matching estCircadianParam)
    rv_idx   <- which(rhythmic_valid)
    pv_order <- order(boot_df$pvalue[rv_idx])
    keep     <- rv_idx[pv_order[seq_len(300L)]]
    rhythmic_valid <- logical(nrow(boot_df))
    rhythmic_valid[keep] <- TRUE
  }
  amplitude_vec      <- boot_df$A[rhythmic_valid]
  sigma_rhythmic_vec <- boot_df$sigma[rhythmic_valid]

  if (length(amplitude_vec) == 0 || all(is.na(amplitude_vec))) {
    amplitude_vec      <- pmax(rlnorm(max(10, nrow(boot_df)), log(0.4), 0.5), 0.05)
    sigma_rhythmic_vec <- NULL
  }
  amplitude_vec <- pmax(amplitude_vec[!is.na(amplitude_vec)], 0.05)

  phase_vec <- boot_df$phi[rhythmic_mask]
  if (length(phase_vec) == 0 || all(is.na(phase_vec))) {
    phase_vec <- runif(max(10, nrow(boot_df)), 0, 24)
  }
  phase_vec <- phase_vec[!is.na(phase_vec)]

  # Build options, pass resolved vectors directly
  opts <- list(
    ngenes        = nrow(boot_df),
    prop_rhythmic = prop_rhythmic,
    period        = bio_diff.opts$period %||% 24,
    lBaselineExpr = lBaselineExpr_vec,
    lBaselineExpr_spec = lBaselineExpr_vec,
    lOD           = lOD_vec,
    lOD_spec      = lOD_vec,
    lOD2          = NULL,
    amplitude      = amplitude_vec,
    amplitude_spec = amplitude_vec,
    sigma_rhythmic = sigma_rhythmic_vec,
    amplitude2     = NULL,
    cts2          = NULL,
    phase         = phase_vec,
    phase_spec    = phase_vec,
    prop_DR        = bio_diff.opts$prop_DR,
    prop_DP        = bio_diff.opts$prop_DP,
    prop_DM        = bio_diff.opts$prop_DM   %||% 0,
    mesor_diff     = bio_diff.opts$mesor_diff %||% c(0.5, 2.0),
    lBaselineExpr2 = bio_diff.opts$lBaselineExpr2,
    phase_diff     = bio_diff.opts$phase_diff,
    amp_diff       = bio_diff.opts$amp_diff,
    dp_shift_mode  = bio_diff.opts$dp_shift_mode,
    dr_amp_scale   = bio_diff.opts$dr_amp_scale %||% 1.0,
    dr_sigma_scale = bio_diff.opts$dr_sigma_scale %||% 1.0,
    sim.seed       = sample.int(.Machine$integer.max, 1)
  )
  class(opts) <- "CircadianBioOptions"
  opts
}


# =====================================================================
# Internal helper: select B evenly-spaced time points from design_vector
# =====================================================================

#' Select B Evenly-Spaced Time Points from a Design Vector
#'
#' @param design_vector Numeric vector. Full set of candidate ZT time points.
#' @param B Integer. Number of distinct time bins to select.
#'
#' @return Numeric vector of length \code{min(B, length(design_vector))} with
#'   evenly-spaced time points drawn from \code{design_vector}.
.selectTimePoints <- function(design_vector, B) {
  n_full <- length(design_vector)
  if (B >= n_full) return(design_vector)
  idx <- round(seq(1, n_full, length.out = B))
  idx <- pmin(pmax(idx, 1), n_full)
  design_vector[unique(idx)]
}


# =====================================================================
# Main Bootstrap Design Grid Runner
# =====================================================================

#' Run bootstrap design grid simulation
#'
#' Estimates uncertainty in power across (N, B) design configurations by
#' bootstrapping parameter sets from pilot data.
#'
#' @param pilot_data Genes x samples matrix (pilot experiment).
#' @param pilot_times Numeric vector of sample time points, length = ncol(pilot_data).
#' @param boot.opts \code{CircadianBootstrapOptions} from \code{CircadianBootstrapOptions()}.
#'   Must set \code{B_values} to a single value for passive designs.
#' @param analysis.opts \code{CircadianAnalysisOptions}.
#' @param bio_diff.opts \code{CircadianBioOptions} from \code{estCircadianParam()}. Required
#'   even in \code{mode = "single"}, differential proportion fields are read from here.
#' @param pilot_data_2 Optional group-2 pilot matrix for \code{mode = "differential"}.
#' @param pilot_times_2 Optional group-2 pilot times.
#' @param mode \code{"single"} or \code{"differential"}.
#' @param methods Detection method(s): \code{"DCP"} (default), \code{"JTK"}, \code{"RAIN"}.
#' @param test_types Endpoints for differential mode: any of \code{"DR"}, \code{"DP"}, \code{"DM"}.
#' @param alpha2 Second-harmonic coefficient (default 0).
#' @param alpha3 Third-harmonic coefficient (default 0).
#' @param min_rhythm_pval P-value threshold for pilot rhythmicity (default 0.01).
#' @param resample Resampling scheme for the outer bootstrap. \code{"subject"}
#'   (default) is the Efron nonparametric bootstrap over pilot subjects: each
#'   draw resamples subject columns with replacement (carrying their collection
#'   times) and refits the pilot summary, matching the Section 2.3 description.
#'   \code{"gene"} is the legacy scheme that fits once and resamples gene rows of
#'   the fit; retained for backward compatibility and diagnostics only.
#' @param verbose Print progress (default TRUE).
#' @param mc.cores Parallel cores for inner simulation loop (default 1).
#'
#' @return Named list with 15 elements:
#'   \describe{
#'     \item{\code{N_values}}{Candidate sample sizes tested.}
#'     \item{\code{B_values}}{Candidate time-point counts tested.}
#'     \item{\code{m_matrix}}{Matrix of total samples (N x B).}
#'     \item{\code{nboot}}{Number of bootstrap replicates.}
#'     \item{\code{design}}{Design type (\code{"active"} or \code{"passive"}).}
#'     \item{\code{test_types}}{Endpoints evaluated.}
#'     \item{\code{fdr_threshold}}{FDR threshold used.}
#'     \item{\code{boot_power}}{4-D array \code{[nboot x n_N x n_B x n_tests]}.}
#'     \item{\code{power_mean}}{Mean power array \code{[n_N x n_B x n_tests]}.}
#'     \item{\code{power_se}}{Bootstrap SE array, same shape as \code{power_mean}.}
#'     \item{\code{power_ci_lo}}{2.5th percentile CI array.}
#'     \item{\code{power_ci_hi}}{97.5th percentile CI array.}
#'     \item{\code{optimal_B}}{Integer vector (length n_N): B with highest mean power at each N.}
#'     \item{\code{optimal_B_ci_lo}}{Bootstrap 2.5th percentile of optimal B.}
#'     \item{\code{optimal_B_ci_hi}}{Bootstrap 97.5th percentile of optimal B.}
#'   }
#' @export
runBootstrapDesignGrid <- function(pilot_data,
                                   pilot_times,
                                   boot.opts,
                                   analysis.opts,
                                   bio_diff.opts,
                                   pilot_data_2    = NULL,
                                   pilot_times_2   = NULL,
                                   mode            = NULL,
                                   methods         = "DCP",
                                   test_types      = c("DR", "DP", "DM"),
                                   alpha2          = 0,
                                   alpha3          = 0,
                                   min_rhythm_pval = 0.01,
                                   resample        = c("subject", "gene"),
                                   verbose         = TRUE,
                                   mc.cores        = 1L) {

  stopifnot(inherits(boot.opts, "CircadianBootstrapOptions"))
  stopifnot(inherits(analysis.opts, "CircadianAnalysisOptions"))
  stopifnot(inherits(bio_diff.opts, "CircadianBioOptions"))
  stopifnot(ncol(pilot_data) == length(pilot_times))

  # Auto-detect mode from pilot data and bio.opts
  if (is.null(mode)) {
    mode <- if (!is.null(pilot_data_2)) "differential" else "single"
  }
  mode <- match.arg(mode, c("single", "differential"))

  N_values    <- boot.opts$N_values
  B_values    <- boot.opts$B_values
  nboot       <- boot.opts$nboot
  nsims_inner <- boot.opts$nsims_inner
  design      <- boot.opts$design
  seed        <- boot.opts$seed
  design_vector <- boot.opts$design_vector

  n_N    <- length(N_values)
  n_B    <- length(B_values)
  period <- bio_diff.opts$period %||% 24

  # Single-cohort mode uses one test type ("rhythmic"); differential uses DR/DP/DM
  if (mode == "single") {
    all_tests <- "rhythmic"
  } else {
    all_tests <- c()
    if (bio_diff.opts$prop_DR > 0) all_tests <- c(all_tests, "DR")
    if (bio_diff.opts$prop_DP > 0) all_tests <- c(all_tests, "DP")
    if (!is.null(bio_diff.opts$prop_DM) && bio_diff.opts$prop_DM > 0) all_tests <- c(all_tests, "DM")
    if (length(all_tests) == 0) all_tests <- "DR"
  }
  n_tests <- length(all_tests)

  # m matrix: m[i, j] = N / B (exact integer division only).
  # NA marks invalid (N, B) pairs where B does not divide N evenly or B > N.
  # Simulating with B > N or round(N/B) would produce a different actual N than
  # requested, mislabelling the power surface; those cells are skipped instead.
  m_matrix <- outer(N_values, B_values, FUN = function(N, B) {
    ifelse(B > N | (N %% B) != 0L, NA_integer_, as.integer(N / B))
  })

  # Warn about (N, B) pairs where B does not divide N exactly
  invalid_pairs <- which(is.na(m_matrix), arr.ind = TRUE)
  if (nrow(invalid_pairs) > 0 && verbose) {
    msg <- apply(invalid_pairs, 1, function(idx)
      sprintf("N=%d/B=%d", N_values[idx[1]], B_values[idx[2]]))
    warning(sprintf(
      "Skipping %d (N, B) pair(s) where B does not divide N exactly: %s",
      nrow(invalid_pairs), paste(msg, collapse = ", ")))
  }

  if (verbose) {
    message(paste(
      "",
      "=== Bootstrap Design Grid ===",
      sprintf("  Pilot data: %d genes x %d samples", nrow(pilot_data), ncol(pilot_data)),
      sprintf("  N_values:   %s", paste(N_values, collapse = ", ")),
      sprintf("  B_values:   %s", paste(B_values, collapse = ", ")),
      sprintf("  nboot:      %d", nboot),
      sprintf("  nsims_inner:%d", nsims_inner),
      sep = "\n"))
  }

  # Step 1-2: Prepare the resampled pilot draws.
  #  resample = "subject" (default): Efron nonparametric bootstrap over pilot
  #    SUBJECTS. Each draw resamples subject columns with replacement (carrying
  #    their collection times) and REFITS the cosinor inside the worker loop,
  #    matching the Section 2.3 description ("resample subjects with replacement,
  #    recompute the pilot summary"). This is the pilot-uncertainty bootstrap.
  #  resample = "gene": legacy scheme that fits once and resamples gene rows;
  #    retained for backward compatibility and diagnostics only.
  resample <- match.arg(resample)
  if (identical(resample, "subject")) {
    if (verbose) cat("\nResampling pilot subjects (subject bootstrap)...\n")
    subj_idx  <- .resampleSubjectIndices(ncol(pilot_data), nboot, seed = seed)
    boot_list <- NULL
  } else {
    if (verbose) cat("\nFitting cosinor to pilot data (gene bootstrap)...\n")
    param_df  <- fitCosinorAll(pilot_data, pilot_times, period = period,
                               min_rhythm_pval = min_rhythm_pval)
    boot_list <- bootstrapParams(param_df, nboot, seed = seed)
    subj_idx  <- NULL
  }

  # Step 3: Sweep over (b, N, B), bootstrap draws parallelized via mclapply.
  # Each draw is independent (different seed); results combined into boot_power_arr.
  fdr_threshold <- min(analysis.opts$fdr_thresholds)

  if (verbose) cat(sprintf("Running %d bootstrap draws (mc.cores=%d)...\n", nboot, mc.cores))

  boot_results <- parallel::mclapply(seq_len(nboot), function(b) {
    if (identical(resample, "subject")) {
      j          <- subj_idx[[b]]
      refit_df   <- fitCosinorAll(pilot_data[, j, drop = FALSE], pilot_times[j],
                                  period = period, min_rhythm_pval = min_rhythm_pval)
      bio_b      <- .buildBioFromBoot(refit_df, bio_diff.opts)
      boot_times <- pilot_times[j]     # passive TOD is resampled jointly with the subjects
    } else {
      bio_b      <- .buildBioFromBoot(boot_list[[b]], bio_diff.opts)
      boot_times <- design_vector
    }
    bio_b$sim.seed <- seed + 1000L + b  # reproducible inner-simulation seed per draw
    result_b <- array(NA_real_, dim = c(n_N, n_B, n_tests))

    for (n_idx in seq_len(n_N)) {
      N <- N_values[n_idx]

      for (B_idx in seq_len(n_B)) {
        B <- B_values[B_idx]
        m <- m_matrix[n_idx, B_idx]
        if (is.na(m)) next
        actual_N <- B * m

        if (design == "active") {
          time_pts <- .selectTimePoints(design_vector, B)
          cts_exp  <- rep(time_pts, each = m)
          if (length(cts_exp) != actual_N) {
            cts_exp <- rep(time_pts, times = ceiling(actual_N / length(time_pts)))[seq_len(actual_N)]
          }
        } else {
          cts_exp <- boot_times
        }

        bio_b_n        <- bio_b
        # Use the user-specified ngenes (from bio_diff.opts) rather than the raw gene count
        # from the bootstrap draw. This ensures the bootstrap simulation uses the same
        # number of genes as the plug-in for consistent FDR calibration.
        bio_b_n$ngenes <- bio_diff.opts$ngenes %||% nrow(boot_list[[b]])

        iter_design <- CircadianDesignOptions(
          sample_sizes = actual_N,
          nsims        = nsims_inner,
          design       = design,
          cts          = cts_exp,
          B_values     = B,
          test_types   = all_tests
        )

        tryCatch({
          if (mode == "single") {
            sim_out <- runSingleCohortGrid(bio_b_n, iter_design, analysis.opts,
                         methods = "DCP", alpha2 = 0,
                         mc.cores = 1L, verbose = FALSE)
            result_b[n_idx, B_idx, 1L] <- sim_out$power_df$power[1L]
          } else {
            sim_out <- runSimsDiff(bio_b_n, iter_design, analysis.opts)

            for (t_idx in seq_len(n_tests)) {
              tt      <- all_tests[t_idx]
              fdr_key <- paste0("fdr_", tt)
              fdr_arr <- sim_out[[fdr_key]]
              if (is.null(fdr_arr)) next

              powers <- sapply(seq_len(nsims_inner), function(sim_i) {
                fdr_vec       <- fdr_arr[, 1, sim_i]
                diff_type_vec <- sim_out$diff_type[[sim_i]]
                target_idx <- if (tt == "DR") {
                  diff_type_vec %in% c(2, 3)
                } else if (tt == "DP") {
                  diff_type_vec == 4
                } else if (tt == "DM") {
                  diff_type_vec == 5
                } else {
                  rep(FALSE, length(fdr_vec))
                }
                n_target <- sum(target_idx)
                if (n_target == 0) return(NA_real_)
                sum(fdr_vec[target_idx] <= fdr_threshold, na.rm = TRUE) / n_target
              })
              result_b[n_idx, B_idx, t_idx] <- mean(powers, na.rm = TRUE)
            }
          }
        }, error = function(e) {          # NA stays in result_b on error
          warning(sprintf("runBootstrapDesignGrid: cell (n_idx=%d, B_idx=%d) failed (left NA): %s",
                          n_idx, B_idx, conditionMessage(e)))
          NULL
        })
      }
    }
    result_b
  }, mc.cores = mc.cores)

  # Combine per-draw results into boot_power_arr[b, n_idx, B_idx, test_idx]
  boot_power_arr <- array(NA_real_,
    dim = c(nboot, n_N, n_B, n_tests),
    dimnames = list(
      paste0("boot", seq_len(nboot)),
      paste0("N", N_values),
      paste0("B", B_values),
      all_tests
    )
  )
  for (b in seq_len(nboot)) {
    res <- boot_results[[b]]
    if (!inherits(res, "error") && !is.null(res)) {
      boot_power_arr[b, , , ] <- res
    }
  }

  # Step 4: Aggregate summaries over bootstrap dimension
  # Aggregate over dim 1 (bootstrap) manually to guarantee [n_N x n_B x n_tests]
  # shape even when n_B=1 or n_tests=1 (apply() drops singleton dimensions).
  power_mean  <- array(apply(boot_power_arr, c(2,3,4), mean,     na.rm=TRUE), dim=c(n_N,n_B,n_tests))
  power_se    <- array(apply(boot_power_arr, c(2,3,4), sd,       na.rm=TRUE), dim=c(n_N,n_B,n_tests))
  power_ci_lo <- array(apply(boot_power_arr, c(2,3,4), function(x) as.numeric(quantile(x, 0.025, na.rm=TRUE))), dim=c(n_N,n_B,n_tests))
  power_ci_hi <- array(apply(boot_power_arr, c(2,3,4), function(x) as.numeric(quantile(x, 0.975, na.rm=TRUE))), dim=c(n_N,n_B,n_tests))

  # Optimal B for each N (highest mean power), using first test type
  optimal_B     <- integer(n_N)
  optimal_B_ci_lo <- integer(n_N)
  optimal_B_ci_hi <- integer(n_N)

  for (n_idx in seq_len(n_N)) {
    mean_by_B  <- power_mean[n_idx, , 1]
    best_idx   <- which.max(mean_by_B)
    optimal_B[n_idx] <- if (length(best_idx) > 0L) B_values[best_idx] else NA_integer_

    # Bootstrap distribution of optimal B
    # Explicitly reshape to [nboot, n_B] to avoid dimension-drop when n_B == 1
    pm_by_B    <- matrix(boot_power_arr[, n_idx, , 1], nrow = nboot, ncol = n_B)
    opt_B_boot <- apply(pm_by_B, 1, function(row) {
      if (all(is.na(row))) return(NA_integer_)
      B_values[which.max(row)]
    })
    q <- quantile(opt_B_boot, c(0.025, 0.975), na.rm = TRUE)
    optimal_B_ci_lo[n_idx] <- q[1]
    optimal_B_ci_hi[n_idx] <- q[2]
  }

  if (verbose) cat("\nBootstrap design grid complete.\n")

  list(
    N_values        = N_values,
    B_values        = B_values,
    m_matrix        = m_matrix,
    nboot           = nboot,
    design          = design,
    test_types      = all_tests,
    fdr_threshold   = fdr_threshold,
    boot_power      = boot_power_arr,
    power_mean      = power_mean,
    power_se        = power_se,
    power_ci_lo     = power_ci_lo,
    power_ci_hi     = power_ci_hi,
    optimal_B       = optimal_B,
    optimal_B_ci_lo = optimal_B_ci_lo,
    optimal_B_ci_hi = optimal_B_ci_hi
  )
}


# =====================================================================
# Summary
# =====================================================================

#' Summarize bootstrap design grid results
#'
#' @param result Output from runBootstrapDesignGrid()
#' @param test_type Test type to summarize ("DR", "DP", "DM")
#' @param fdr_threshold FDR threshold (used for display label only)
#' @param verbose Print table
#'
#' @return Invisible data frame with summary
summaryBootstrapDesignGrid <- function(result,
                                       test_type     = "DR",
                                       fdr_threshold = NULL,
                                       verbose       = TRUE) {
  t_idx <- match(test_type, result$test_types)
  if (is.na(t_idx)) {
    warning(sprintf("test_type '%s' not found; using first: '%s'",
                    test_type, result$test_types[1]))
    t_idx <- 1
    test_type <- result$test_types[1]
  }

  fdr_thr <- fdr_threshold %||% result$fdr_threshold

  rows <- list()
  for (n_idx in seq_along(result$N_values)) {
    N <- result$N_values[n_idx]
    for (B_idx in seq_along(result$B_values)) {
      B <- result$B_values[B_idx]
      m <- result$m_matrix[n_idx, B_idx]
      pmean  <- result$power_mean[n_idx, B_idx, t_idx]
      plo    <- result$power_ci_lo[n_idx, B_idx, t_idx]
      phi    <- result$power_ci_hi[n_idx, B_idx, t_idx]
      opt    <- result$optimal_B[n_idx] == B

      rows[[length(rows) + 1]] <- data.frame(
        N            = N,
        B            = B,
        m            = m,
        Power_mean   = round(pmean * 100, 1),
        Power_CI_lo  = round(plo  * 100, 1),
        Power_CI_hi  = round(phi  * 100, 1),
        Optimal      = ifelse(opt, "*", ""),
        stringsAsFactors = FALSE
      )
    }
  }

  df <- do.call(rbind, rows)

  if (verbose) {
    cat(sprintf("\nBootstrap Design Grid Summary (%s test at FDR %.0f%%)\n",
                test_type, 100 * fdr_thr))
    cat(sprintf("  nboot = %d, design = %s\n\n", result$nboot, result$design))
    print(df, row.names = FALSE)
    cat("\n* = Optimal B for this N (highest mean power)\n")
  }

  invisible(df)
}


# =====================================================================
# Plotting
# =====================================================================

#' Plot bootstrap design grid results
#'
#' Panel A: Power vs N, one line per B, with bootstrap CI ribbon.
#' Panel B: Heatmap of mean power over (N, B) grid.
#'
#' @param result Output from runBootstrapDesignGrid()
#' @param test_type Test type to plot
#' @param fdr_threshold FDR threshold (label only)
#' @param panels Which panel(s) to draw: "A" (power vs N) and/or "B" (heatmap)
#'   (default "A").
#' @param output_file Path for PDF output (NULL = screen)
plotBootstrapDesignGrid <- function(result,
                                    test_type     = "DR",
                                    fdr_threshold = NULL,
                                    panels        = "A",
                                    output_file   = NULL) {
  panels <- match.arg(panels, choices = c("A", "B"), several.ok = TRUE)

  t_idx <- match(test_type, result$test_types)
  if (is.na(t_idx)) {
    t_idx <- 1
    test_type <- result$test_types[1]
  }

  fdr_thr <- fdr_threshold %||% result$fdr_threshold

  N_values <- result$N_values
  B_values <- result$B_values
  n_N <- length(N_values)
  n_B <- length(B_values)

  power_mean  <- matrix(result$power_mean[, , t_idx],  nrow = n_N, ncol = n_B)
  power_ci_lo <- matrix(result$power_ci_lo[, , t_idx], nrow = n_N, ncol = n_B)
  power_ci_hi <- matrix(result$power_ci_hi[, , t_idx], nrow = n_N, ncol = n_B)

  cols <- rainbow(n_B, s = 0.7, v = 0.85)

  n_panels <- length(panels)
  if (!is.null(output_file)) {
    pdf(output_file, width = if (n_panels == 1) 7 else 12, height = 5)
    on.exit(dev.off())
  }

  old_par <- par(mfrow = c(1, n_panels), mar = c(4.5, 4.5, 3, 1.5))
  on.exit(par(old_par), add = TRUE)

  # -----------------------------------------------------------
  # Panel A: Power vs N, one line per B, with CI ribbon
  # -----------------------------------------------------------
  if (!("A" %in% panels)) {
    # skip to Panel B
  } else {
  y_max <- min(1, max(power_ci_hi, na.rm = TRUE) * 1.05)
  plot(N_values, power_mean[, 1],
       type  = "n",
       xlim  = range(N_values),
       ylim  = c(0, y_max),
       xlab  = "N per group",
       ylab  = "Power",
       main  = sprintf("Bootstrap Power vs N (%s, FDR <= %.0f%%)\nline: bootstrap mean, band: 95%% bootstrap CI", test_type, 100 * fdr_thr),
       las   = 1, xaxt = "n")
  axis(1, at = N_values)
  abline(h = 0.80, lty = 2, col = "gray40")
  text(min(N_values), 0.82, "80%", col = "gray40", cex = 0.75, adj = 0)
  abline(h  = seq(0.2, 0.8, by = 0.2), col = "gray85", lty = 3)

  for (B_idx in seq_len(n_B)) {
    pm  <- power_mean[, B_idx]
    plo <- power_ci_lo[, B_idx]
    phi <- power_ci_hi[, B_idx]

    # CI ribbon
    polygon(c(N_values, rev(N_values)),
            c(plo, rev(phi)),
            col    = adjustcolor(cols[B_idx], alpha.f = 0.15),
            border = NA)
    lines(N_values, pm,  col = cols[B_idx], lwd = 2)
    points(N_values, pm, col = cols[B_idx], pch = 19, cex = 0.7)
  }

  legend("bottomright",
         legend = paste0("B=", B_values),
         col    = cols,
         lwd    = 2, lty = 1,
         title  = "# Time Points",
         cex    = 0.8, bty = "n")
  } # end Panel A

  if (!("B" %in% panels)) return(invisible(NULL))

  # -----------------------------------------------------------
  # Panel B: Heatmap of mean power over (N, B)
  # -----------------------------------------------------------
  col_ramp <- colorRampPalette(c("white", "lightyellow", "orange", "red", "darkred"))(100)
  pm_pct <- power_mean * 100
  image(seq_len(n_N), seq_len(n_B), pm_pct,
        col    = col_ramp,
        zlim   = c(0, 100),
        xaxt   = "n", yaxt = "n",
        xlab   = "N per group",
        ylab   = "B (# time points)",
        main   = sprintf("Mean Power (%%) Grid\n(%s test, FDR <= %.0f%%)",
                         test_type, 100 * fdr_thr))
  axis(1, at = seq_len(n_N), labels = N_values)
  axis(2, at = seq_len(n_B), labels = B_values, las = 1)

  # Cell annotations, explicitly mark infeasible (N, B) pairs
  for (n_idx in seq_len(n_N)) {
    for (B_idx in seq_len(n_B)) {
      val <- pm_pct[n_idx, B_idx]
      if (!is.na(val)) {
        txt_col <- if (val > 60) "white" else "black"
        text(n_idx, B_idx, sprintf("%.0f", val),
             col = txt_col, cex = 0.75)
      } else {
        # Infeasible cell: B does not divide N exactly
        rect(n_idx - 0.48, B_idx - 0.48, n_idx + 0.48, B_idx + 0.48,
             col = "gray80", border = NA)
        text(n_idx, B_idx, "\u2014", col = "gray40", cex = 1.0)  # em-dash
      }
    }
  }

  # Mark optimal B per N
  for (n_idx in seq_len(n_N)) {
    opt_B_idx <- which(B_values == result$optimal_B[n_idx])
    if (length(opt_B_idx) > 0) {
      rect(n_idx - 0.45, opt_B_idx - 0.45, n_idx + 0.45, opt_B_idx + 0.45,
           border = "blue", lwd = 2)
    }
  }
  mtext("Blue border = optimal B (highest mean power)", side = 1,
        line = 3.5, cex = 0.7, col = "blue")

  invisible(NULL)
}
