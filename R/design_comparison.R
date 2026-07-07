#' Design Comparison: Two-Stage vs Bootstrap-Simulation Power
#'
#' Compares two approaches to pilot-informed power analysis for a FIXED design:
#'
#' 1. Two-stage (point estimate): Estimate circadian parameters from pilot once,
#'    simulate with those fixed estimates -> one power curve (no uncertainty).
#'
#' 2. Bootstrap simulation: Resample pilot gene rows (nboot times) to get a
#'    distribution of parameter sets, simulate for each -> power +/- CI.
#'
#' The comparison isolates pilot-estimation uncertainty: both approaches use the
#' same design (same B and m); the only difference is whether the parameter
#' distribution is treated as fixed or as an uncertain quantity.
#'
#' IMPORTANT: Pass bootstrap results from a single fixed B value. Averaging over
#' multiple candidate B values produces a quantity with no actionable scientific
#' meaning, because the user must choose one specific design, not the mean of
#' several designs with different B/m splits.


# =====================================================================
# Two-Stage Power Runner
# =====================================================================

#' Run two-stage power analysis (point estimate approach)
#'
#' Estimates circadian parameters from pilot data once, then runs power
#' simulations using those point-estimate parameters.
#'
#' @param pilot_data      Genes x samples matrix
#' @param pilot_times     Sample time points (length = ncol(pilot_data))
#' @param design.opts     CircadianDesignOptions
#' @param analysis.opts   CircadianAnalysisOptions
#' @param bio_diff.opts   CircadianBioOptions (for differential params only)
#' @param min_rhythm_pval Threshold for calling rhythmic genes in pilot
#' @param test_type       "DR", "DP", or "DM"
#' @param verbose         Print progress
#'
#' @return List with marginal_power matrix [n_sizes x nsims] and sample_sizes
runTwoStagePower <- function(pilot_data,
                              pilot_times,
                              design.opts,
                              analysis.opts,
                              bio_diff.opts,
                              min_rhythm_pval = 0.01,
                              test_type       = "DR",
                              verbose         = TRUE,
                              mc.cores        = 1L) {

  stopifnot(inherits(design.opts,   "CircadianDesignOptions"))
  stopifnot(inherits(analysis.opts, "CircadianAnalysisOptions"))
  stopifnot(inherits(bio_diff.opts, "CircadianBioOptions"))
  stopifnot(ncol(pilot_data) == length(pilot_times))

  if (verbose) cat("\nTwo-stage power analysis\n")

  # Step 1: Estimate parameters from pilot data
  bio_pilot <- estCircadianParam(
    data             = pilot_data,
    times            = pilot_times,
    period           = bio_diff.opts$period %||% 24,
    prop_DR          = bio_diff.opts$prop_DR,
    prop_DP          = bio_diff.opts$prop_DP,
    prop_DM          = bio_diff.opts$prop_DM   %||% 0,
    mesor_diff       = bio_diff.opts$mesor_diff %||% c(0.5, 2.0),
    min_rhythm_pval  = min_rhythm_pval,
    verbose          = verbose
  )

  # Step 2: Override remaining differential params from bio_diff.opts
  if (!is.null(bio_diff.opts$phase_diff)) {
    bio_pilot <- updateBioOptions(bio_pilot, phase_diff = bio_diff.opts$phase_diff)
  }
  if (!is.null(bio_diff.opts$amp_diff)) {
    bio_pilot <- updateBioOptions(bio_pilot, amp_diff = bio_diff.opts$amp_diff)
  }
  bio_pilot <- updateBioOptions(bio_pilot,
                                 dp_shift_mode  = bio_diff.opts$dp_shift_mode,
                                 dr_amp_scale   = bio_diff.opts$dr_amp_scale %||% 1.0,
                                 dr_sigma_scale = bio_diff.opts$dr_sigma_scale %||% 1.0)

  # Step 3: Run power analysis using point estimates
  sample_sizes  <- design.opts$sample_sizes
  nsims         <- design.opts$nsims
  fdr_threshold <- min(analysis.opts$fdr_thresholds)
  n_sizes       <- length(sample_sizes)

  marginal_power <- matrix(NA_real_, nrow = n_sizes, ncol = nsims)

  run_one_n <- function(j) {
    n <- sample_sizes[j]
    if (verbose) cat(sprintf("  n = %d\n", n))

    # Expand cts to per-sample schedule: active designs require length(cts)==n1
    cts_n <- if (!is.null(design.opts$cts) && design.opts$design == "active") {
      B_pts <- design.opts$cts
      m     <- n / length(B_pts)
      rep(B_pts, each = round(m))
    } else {
      design.opts$cts
    }

    iter_design <- CircadianDesignOptions(
      sample_sizes = n,
      nsims        = nsims,
      design       = design.opts$design,
      cts          = cts_n,
      test_types   = c(test_type)
    )

    pow_row <- rep(NA_real_, nsims)
    tryCatch({
      sim_out <- runSimsDiff(bio_pilot, iter_design, analysis.opts)
      pow_row <- .computeMarginalPower(sim_out, test_type, fdr_threshold, nsims)
    }, error = function(e) {
      if (verbose) warning(sprintf("  Two-stage n=%d failed: %s", n, e$message))
    })
    pow_row
  }

  results <- parallel::mclapply(seq_len(n_sizes), run_one_n, mc.cores = mc.cores)
  for (j in seq_len(n_sizes)) {
    if (!is.null(results[[j]])) marginal_power[j, ] <- results[[j]]
  }

  if (verbose) cat("Two-stage power analysis complete.\n")

  list(
    sample_sizes   = sample_sizes,
    marginal_power = marginal_power,
    power_mean     = rowMeans(marginal_power, na.rm = TRUE),
    power_se       = apply(marginal_power, 1, function(x) {
                      x <- x[!is.na(x)]
                      if (length(x) < 2) NA_real_ else sd(x) / sqrt(length(x))
                    }),
    test_type      = test_type,
    fdr_threshold  = fdr_threshold
  )
}


# =====================================================================
# Compare Approaches
# =====================================================================

#' Compare two-stage vs bootstrap power estimates
#'
#' Merges per-N results from both approaches into a comparison data frame.
#' For bootstrap results, averages over B values to get marginal power per N.
#'
#' @param two_stage_result  Output from runTwoStagePower()
#' @param bootstrap_result  Output from runBootstrapDesignGrid()
#' @param test_type         Test type to extract from bootstrap result
#' @param target_power      Target power level for n80 calculation
#'
#' @return List with comparison data frame and n_80 summaries
compareDesignApproaches <- function(two_stage_result,
                                    bootstrap_result,
                                    test_type    = "DR",
                                    target_power = 0.80) {

  # --- Two-stage ---
  ts_sizes <- two_stage_result$sample_sizes
  ts_power <- two_stage_result$power_mean
  ts_se    <- two_stage_result$power_se

  # --- Bootstrap: extract per-N power ---
  # IMPORTANT: bootstrap_result should come from a single fixed B value.
  # If multiple B values are present, marginal power is averaged over B,
  # which conflates design variation with parameter uncertainty.
  t_idx <- match(test_type, bootstrap_result$test_types)
  if (is.na(t_idx)) t_idx <- 1

  if (length(bootstrap_result$B_values) > 1) {
    stop(paste(
      "bootstrap_result contains multiple B values.",
      "compareDesignApproaches() requires a single fixed-B bootstrap result so that",
      "CI bands reflect parameter-estimation uncertainty only, not design variation.",
      "\n  Run runBootstrapDesignGrid() with B_values of length 1 and pass that result here.",
      "\n  B values found:", paste(bootstrap_result$B_values, collapse = ", ")
    ))
  }

  boot_sizes   <- bootstrap_result$N_values
  boot_mean    <- rowMeans(bootstrap_result$power_mean[, , t_idx,  drop = FALSE], na.rm = TRUE)
  boot_ci_lo   <- rowMeans(bootstrap_result$power_ci_lo[, , t_idx, drop = FALSE], na.rm = TRUE)
  boot_ci_hi   <- rowMeans(bootstrap_result$power_ci_hi[, , t_idx, drop = FALSE], na.rm = TRUE)

  # Align on common N values
  common_N <- intersect(ts_sizes, boot_sizes)

  ts_idx   <- match(common_N, ts_sizes)
  boot_idx <- match(common_N, boot_sizes)

  comparison <- data.frame(
    n               = common_N,
    two_stage_power = ts_power[ts_idx],
    two_stage_se    = ts_se[ts_idx],
    boot_power_mean = boot_mean[boot_idx],
    boot_ci_lo      = boot_ci_lo[boot_idx],
    boot_ci_hi      = boot_ci_hi[boot_idx],
    stringsAsFactors = FALSE
  )

  # --- n₈₀ estimates ---
  # Two-stage: smallest N where power_mean >= target_power
  n80_ts_idx <- which(ts_power >= target_power)
  n80_ts <- if (length(n80_ts_idx) > 0) ts_sizes[min(n80_ts_idx)] else NA_integer_

  # Bootstrap: across nboot draws, compute distribution of n₈₀
  # For each bootstrap draw b, average over B → marginal power vs N curve
  nboot <- bootstrap_result$nboot

  n80_boot_vec <- sapply(seq_len(nboot), function(b) {
    # Average over B for each N: apply over dim 2 (N) of the [1 x N x B x 1] slice
    tmp <- bootstrap_result$boot_power[b, , , t_idx, drop = FALSE]
    boot_pm_b <- apply(tmp, 2, mean, na.rm = TRUE)
    idx <- which(boot_pm_b >= target_power)
    if (length(idx) == 0) return(NA_real_)
    bootstrap_result$N_values[min(idx)]
  })
  n80_boot_vec <- n80_boot_vec[!is.na(n80_boot_vec)]

  n80_boot_median  <- if (length(n80_boot_vec) > 0) median(n80_boot_vec) else NA
  n80_boot_lo      <- if (length(n80_boot_vec) > 0) quantile(n80_boot_vec, 0.025) else NA
  n80_boot_hi      <- if (length(n80_boot_vec) > 0) quantile(n80_boot_vec, 0.975) else NA

  list(
    comparison       = comparison,
    test_type        = test_type,
    target_power     = target_power,
    n80_two_stage    = n80_ts,
    n80_boot_median  = n80_boot_median,
    n80_boot_lo      = n80_boot_lo,
    n80_boot_hi      = n80_boot_hi
  )
}


# =====================================================================
# Plotting
# =====================================================================

#' Plot design approach comparison
#'
#' Panel A: Power curves, two-stage (solid blue) vs bootstrap (dashed orange + CI band).
#' Panel B: Recommended N at target_power, bar chart comparing methods.
#'
#' @param comparison    Output from compareDesignApproaches()
#' @param test_type     Test type label
#' @param target_power  Target power (horizontal reference line)
#' @param output_file   Path for PDF output (NULL = screen)
plotDesignComparison <- function(comparison,
                                  test_type    = NULL,
                                  target_power = NULL,
                                  panels       = c("A", "B"),
                                  output_file  = NULL) {

  test_type    <- test_type    %||% comparison$test_type    %||% "DR"
  target_power <- target_power %||% comparison$target_power %||% 0.80

  df       <- comparison$comparison
  n_values <- df$n
  fdr_thr  <- NULL  # not stored; use label from test_type

  panels <- match.arg(panels, choices = c("A", "B"), several.ok = TRUE)
  n_panels <- length(panels)

  if (!is.null(output_file)) {
    pdf(output_file, width = if (n_panels == 1) 7 else 12, height = 5)
    on.exit(dev.off())
  }

  old_par <- par(mfrow = c(1, n_panels), mar = c(4.5, 4.5, 3, 1.5))
  on.exit(par(old_par), add = TRUE)

  # -----------------------------------------------------------
  # Panel A: Power curves
  # -----------------------------------------------------------
  if ("A" %in% panels) {
  y_max <- min(1.05, max(c(df$two_stage_power, df$boot_ci_hi), na.rm = TRUE) * 1.05)

  plot(n_values, df$two_stage_power,
       type = "n",
       xlim = range(n_values),
       ylim = c(0, y_max),
       xlab = "N per group",
       ylab = "Power",
       main = sprintf("Power Curves: Two-Stage vs Bootstrap\n(%s test)", test_type),
       las  = 1, xaxt = "n")
  axis(1, at = n_values)

  # Bootstrap CI ribbon
  polygon(c(n_values, rev(n_values)),
          c(df$boot_ci_lo, rev(df$boot_ci_hi)),
          col    = adjustcolor("darkorange", alpha.f = 0.15),
          border = NA)

  # Bootstrap mean line
  lines(n_values, df$boot_power_mean, col = "darkorange", lwd = 2, lty = 2)

  # Two-stage line + SE bars
  lines(n_values, df$two_stage_power, col = "steelblue", lwd = 2.5, lty = 1)
  arrows(n_values, df$two_stage_power - df$two_stage_se,
         n_values, df$two_stage_power + df$two_stage_se,
         length = 0.05, angle = 90, code = 3, col = "steelblue", lwd = 1.2)

  # Target power line
  abline(h = target_power, lty = 2, col = "gray40")
  text(min(n_values), target_power + 0.02,
       sprintf("%.0f%%", 100 * target_power),
       col = "gray40", cex = 0.75, adj = 0)

  # Vertical ticks at n₈₀
  n80_ts   <- comparison$n80_two_stage
  n80_boot <- comparison$n80_boot_median
  if (!is.na(n80_ts)) {
    abline(v = n80_ts, col = "steelblue", lty = 3, lwd = 1.5)
  }
  if (!is.null(n80_boot) && !is.na(n80_boot)) {
    abline(v = n80_boot, col = "darkorange", lty = 3, lwd = 1.5)
  }

  legend("bottomright",
         legend = c("Two-stage (point est.)",
                    "Bootstrap (mean)",
                    "Bootstrap 95% CI"),
         col    = c("steelblue", "darkorange", adjustcolor("darkorange", 0.3)),
         lty    = c(1, 2, NA),
         lwd    = c(2.5, 2, NA),
         pch    = c(NA, NA, 15),
         pt.cex = c(NA, NA, 2),
         cex    = 0.8, bty = "n")
  } # end Panel A

  # -----------------------------------------------------------
  # Panel B: n₈₀ point + 95% interval
  #
  # Two rows, not four bars:
  #   Row 1, Two-stage:  a single point (no CI, point estimate only)
  #   Row 2, Bootstrap:  median point + 95% CI bar
  #
  # This correctly represents bootstrap n₈₀ as ONE uncertain recommendation,
  # not three separate competing recommendations.
  # -----------------------------------------------------------
  if (!("B" %in% panels)) return(invisible(NULL))

  n80_ts   <- comparison$n80_two_stage
  n80_med  <- comparison$n80_boot_median
  n80_lo   <- comparison$n80_boot_lo
  n80_hi   <- comparison$n80_boot_hi

  any_valid <- any(!is.na(c(n80_ts, n80_med, n80_lo, n80_hi)))

  if (!any_valid) {
    plot.new()
    title(main = sprintf("N at %.0f%% Power\n(%s test)", 100 * target_power, test_type))
    text(0.5, 0.5,
         sprintf("Power did not reach %.0f%%\nfor any N tested.\nIncrease N range or nsims.",
                 100 * target_power),
         cex = 1.1, col = "gray40", adj = c(0.5, 0.5))
    return(invisible(NULL))
  }

  all_vals <- c(n80_ts, n80_med, n80_lo, n80_hi)
  x_min <- min(all_vals, na.rm = TRUE) * 0.85
  x_max <- max(all_vals, na.rm = TRUE) * 1.15

  plot(NA, NA,
       xlim = c(x_min, x_max),
       ylim = c(0.5, 2.5),
       xlab = "Recommended N per group",
       ylab = "",
       main = sprintf("N at %.0f%% Power\n(%s test)", 100 * target_power, test_type),
       yaxt = "n", las = 1)
  axis(2, at = c(1, 2), labels = c("Bootstrap", "Two-stage"), las = 1, tick = FALSE)
  abline(v = seq(x_min, x_max, length.out = 6), col = "gray90", lty = 3)

  # Two-stage: single point
  if (!is.na(n80_ts)) {
    points(n80_ts, 2, pch = 18, col = "steelblue", cex = 2.2)
    text(n80_ts, 2.25, sprintf("N = %g", n80_ts),
         col = "steelblue", cex = 0.85, adj = 0.5)
  } else {
    text(mean(c(x_min, x_max)), 2, "N/A", col = "gray50", cex = 0.9, adj = 0.5)
  }

  # Bootstrap: median + CI bar
  if (!is.na(n80_med)) {
    if (!is.na(n80_lo) && !is.na(n80_hi)) {
      lines(c(n80_lo, n80_hi), c(1, 1), col = "darkorange", lwd = 2.5)
      # End caps
      for (xv in c(n80_lo, n80_hi)) {
        lines(c(xv, xv), c(0.88, 1.12), col = "darkorange", lwd = 2.5)
      }
      text(n80_lo - (x_max - x_min) * 0.02, 1.28,
           sprintf("%.0f%%\nN=%g", 2.5, n80_lo),
           col = "darkorange", cex = 0.72, adj = 1)
      text(n80_hi + (x_max - x_min) * 0.02, 1.28,
           sprintf("%.0f%%\nN=%g", 97.5, n80_hi),
           col = "darkorange", cex = 0.72, adj = 0)
    }
    points(n80_med, 1, pch = 19, col = "darkorange", cex = 2.2)
    text(n80_med, 0.72, sprintf("med.\nN=%g", n80_med),
         col = "darkorange", cex = 0.72, adj = 0.5)
  } else {
    text(mean(c(x_min, x_max)), 1, "N/A", col = "gray50", cex = 0.9, adj = 0.5)
  }

  mtext("CI width = pilot-estimation uncertainty (bootstrap outer loop)",
        side = 1, line = 3.5, cex = 0.7, col = "gray30")

  invisible(NULL)
}


# =====================================================================
# Ground Truth Calibration
# =====================================================================

#' Internal helper: compute marginal power from runSimsDiff() output
#'
#' @param sim_out Output from runSimsDiff() with sample_sizes of length 1
#' @param test_type "DR", "DP", or "DM"
#' @param fdr_threshold Numeric threshold
#' @param nsims Number of simulations (3rd dim of fdr array)
#' @return Numeric vector of length nsims (power per simulation)
.computeMarginalPower <- function(sim_out, test_type, fdr_threshold, nsims) {
  fdr_key <- paste0("fdr_", test_type)
  fdr_arr <- sim_out[[fdr_key]]   # [ngenes, 1, nsims]
  if (is.null(fdr_arr)) return(rep(NA_real_, nsims))

  sapply(seq_len(nsims), function(sim_i) {
    fdr_vec   <- fdr_arr[, 1, sim_i]
    diff_type <- sim_out$diff_type[[sim_i]]

    target_idx <- switch(test_type,
      DR = diff_type %in% c(2, 3),
      DP = diff_type == 4,
      DM = diff_type == 5,
      rep(FALSE, length(fdr_vec))
    )

    n_target <- sum(target_idx)
    if (n_target == 0) return(NA_real_)
    sum(fdr_vec[target_idx] <= fdr_threshold, na.rm = TRUE) / n_target
  })
}


#' Generate synthetic pilot data from known true parameters
#'
#' Simulates a single-group cross-sectional dataset (genes x n_pilot) with
#' subject times sampled from pilot_times. Used as a stand-in for real pilot
#' data when the ground truth is known.
#'
#' @param true_bio.opts CircadianBioOptions with known true parameters
#' @param n_pilot       Number of pilot subjects
#' @param pilot_times   TOD distribution to sample subject times from
#' @param seed          Random seed
#'
#' @return List: data (genes x n_pilot matrix), times (length n_pilot)
#' @export
generatePilotData <- function(true_bio.opts, n_pilot, pilot_times, seed = 42) {
  stopifnot(inherits(true_bio.opts, "CircadianBioOptions"))
  stopifnot(n_pilot > 0, length(pilot_times) > 0)

  sim_out <- simCircadianDiff(
    ngenes        = true_bio.opts$ngenes,
    n1            = n_pilot,
    n2            = 1,                      # dummy second group
    lBaselineExpr = true_bio.opts$lBaselineExpr,
    lOD           = true_bio.opts$lOD,
    amplitude     = true_bio.opts$amplitude,
    prop_rhythmic = true_bio.opts$prop_rhythmic,
    prop_DR       = 0,
    prop_DP       = 0,
    period        = true_bio.opts$period %||% 24,
    design        = "passive",
    cts           = pilot_times,
    sim.seed      = seed
  )

  list(data = sim_out$expr1, times = sim_out$times1)
}


#' Ground truth calibration: compare two-stage vs bootstrap when truth is known
#'
#' Generates synthetic pilot data from known parameters, then evaluates how
#' well two-stage (point estimate) and bootstrap (distribution) approaches
#' recover the true power curve.
#'
#' @param true_bio.opts  CircadianBioOptions, the known ground truth biology
#'                       (includes both base params AND differential params)
#' @param design.opts    CircadianDesignOptions, N values and design to evaluate
#' @param analysis.opts  CircadianAnalysisOptions
#' @param n_pilot        Number of synthetic pilot subjects
#' @param pilot_times    TOD distribution for pilot time sampling
#' @param nboot          Bootstrap draws for uncertainty quantification
#' @param nsims_inner    Simulations per bootstrap draw
#' @param test_type      "DR", "DP", or "DM"
#' @param seed           Random seed
#' @param verbose        Print progress
#'
#' @return List with true_power, two-stage power, and bootstrap power +/- CI
runGroundTruthComparison <- function(true_bio.opts,
                                      design.opts,
                                      analysis.opts,
                                      n_pilot,
                                      pilot_times,
                                      nboot        = 50,
                                      nsims_inner  = 10,
                                      test_type    = "DR",
                                      seed         = 42,
                                      verbose      = TRUE) {

  stopifnot(inherits(true_bio.opts, "CircadianBioOptions"))
  stopifnot(inherits(design.opts,   "CircadianDesignOptions"))
  stopifnot(inherits(analysis.opts, "CircadianAnalysisOptions"))

  fdr_threshold <- min(analysis.opts$fdr_thresholds)
  sample_sizes  <- design.opts$sample_sizes
  nsims         <- design.opts$nsims
  n_sizes       <- length(sample_sizes)
  period        <- true_bio.opts$period %||% 24

  # ------------------------------------------------------------------
  # Step 1: Generate synthetic pilot from known true params
  # ------------------------------------------------------------------
  if (verbose) cat("\nGround truth calibration\n")
  if (verbose) cat(sprintf("Generating synthetic pilot (n=%d)...\n", n_pilot))

  pilot <- generatePilotData(true_bio.opts, n_pilot, pilot_times, seed = seed)

  # ------------------------------------------------------------------
  # Step 2: Oracle (true) power, use exact known params
  # ------------------------------------------------------------------
  if (verbose) cat("Computing oracle power (known truth)...\n")

  true_power_mean <- rep(NA_real_, n_sizes)
  true_power_se   <- rep(NA_real_, n_sizes)

  for (j in seq_len(n_sizes)) {
    n <- sample_sizes[j]
    if (verbose) cat(sprintf("  Oracle n = %d\n", n))

    iter_design <- CircadianDesignOptions(
      sample_sizes = n,
      nsims        = nsims,
      design       = design.opts$design,
      cts          = design.opts$cts,
      test_types   = c(test_type)
    )

    tryCatch({
      sim_out <- runSimsDiff(true_bio.opts, iter_design, analysis.opts)
      powers  <- .computeMarginalPower(sim_out, test_type, fdr_threshold, nsims)
      true_power_mean[j] <- mean(powers, na.rm = TRUE)
      true_power_se[j]   <- sd(powers,   na.rm = TRUE)
    }, error = function(e) {
      if (verbose) warning(sprintf("  Oracle n=%d failed: %s", n, e$message))
    })
  }

  # ------------------------------------------------------------------
  # Step 3: Two-stage, estimate from pilot, then simulate
  # ------------------------------------------------------------------
  if (verbose) cat("Running two-stage on synthetic pilot...\n")

  ts_result <- runTwoStagePower(
    pilot_data    = pilot$data,
    pilot_times   = pilot$times,
    design.opts   = design.opts,
    analysis.opts = analysis.opts,
    bio_diff.opts = true_bio.opts,
    test_type     = test_type,
    verbose       = verbose
  )

  # ------------------------------------------------------------------
  # Step 4: Bootstrap, resample pilot params, simulate for each draw
  # ------------------------------------------------------------------
  if (verbose) cat("Running bootstrap on synthetic pilot...\n")

  param_df  <- fitCosinorAll(pilot$data, pilot$times, period = period)
  boot_list <- bootstrapParams(param_df, nboot, seed = seed + 1)

  boot_power_mat <- matrix(NA_real_, nrow = nboot, ncol = n_sizes)

  for (b in seq_len(nboot)) {
    if (verbose && b %% max(1L, floor(nboot / 5L)) == 0L) {
      message(sprintf("  Bootstrap draw %d / %d", b, nboot))
    }

    bio_b <- .buildBioFromBoot(boot_list[[b]], true_bio.opts)

    for (j in seq_len(n_sizes)) {
      n <- sample_sizes[j]

      iter_design <- CircadianDesignOptions(
        sample_sizes = n,
        nsims        = nsims_inner,
        design       = design.opts$design,
        cts          = design.opts$cts,
        test_types   = c(test_type)
      )

      tryCatch({
        sim_out <- runSimsDiff(bio_b, iter_design, analysis.opts)
        powers  <- .computeMarginalPower(sim_out, test_type, fdr_threshold, nsims_inner)
        boot_power_mat[b, j] <- mean(powers, na.rm = TRUE)
      }, error = function(e) NULL)
    }
  }

  boot_mean  <- colMeans(boot_power_mat, na.rm = TRUE)
  boot_ci_lo <- apply(boot_power_mat, 2, quantile, probs = 0.025, na.rm = TRUE)
  boot_ci_hi <- apply(boot_power_mat, 2, quantile, probs = 0.975, na.rm = TRUE)

  # Coverage: proportion of N values where CI contains true power
  coverage <- mean(
    true_power_mean >= boot_ci_lo & true_power_mean <= boot_ci_hi,
    na.rm = TRUE
  )

  if (verbose) {
    message(sprintf("\nCoverage (CI contains truth): %.0f%% of N values",
                    100 * coverage))
  }

  list(
    sample_sizes     = sample_sizes,
    true_power_mean  = true_power_mean,
    true_power_se    = true_power_se,
    ts_power_mean    = ts_result$power_mean,
    ts_power_se      = ts_result$power_se,
    boot_power_mean  = boot_mean,
    boot_ci_lo       = boot_ci_lo,
    boot_ci_hi       = boot_ci_hi,
    boot_power_mat   = boot_power_mat,
    coverage         = coverage,
    pilot_n          = n_pilot,
    nboot            = nboot,
    test_type        = test_type,
    fdr_threshold    = fdr_threshold
  )
}


#' Plot ground truth calibration results
#'
#' Panel A: Power vs N, oracle (black), two-stage (blue), bootstrap mean +/- CI (orange).
#' Panel B: Per-N calibration, point estimates vs CI bars vs true value.
#'
#' @param gt_result    Output from runGroundTruthComparison()
#' @param test_type    Test type label (uses gt_result$test_type if NULL)
#' @param target_power Reference power line (default 0.80)
#' @param output_file  Path for PDF output (NULL = screen)
plotGroundTruthComparison <- function(gt_result,
                                       test_type    = NULL,
                                       target_power = 0.80,
                                       output_file  = NULL) {

  test_type    <- test_type %||% gt_result$test_type %||% "DR"
  N            <- gt_result$sample_sizes
  n_sizes      <- length(N)
  fdr_thr      <- gt_result$fdr_threshold

  if (!is.null(output_file)) {
    pdf(output_file, width = 12, height = 5)
    on.exit(dev.off())
  }

  old_par <- par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.5, 1.5))
  on.exit(par(old_par), add = TRUE)

  # -----------------------------------------------------------
  # Panel A: Power curves, oracle vs two-stage vs bootstrap
  # -----------------------------------------------------------
  y_max <- min(1.05,
    max(c(gt_result$true_power_mean + gt_result$true_power_se,
          gt_result$boot_ci_hi), na.rm = TRUE) * 1.08)

  plot(N, gt_result$true_power_mean,
       type = "n",
       xlim = range(N),
       ylim = c(0, y_max),
       xlab = "N per group",
       ylab = "Power",
       main = sprintf("Ground Truth Calibration\n(%s test, pilot n = %d, FDR <= %.0f%%)",
                      test_type, gt_result$pilot_n, 100 * fdr_thr),
       las = 1, xaxt = "n")
  axis(1, at = N)
  abline(h = target_power, lty = 2, col = "gray40")
  text(min(N), target_power + 0.02,
       sprintf("%.0f%%", 100 * target_power),
       col = "gray40", cex = 0.75, adj = 0)
  abline(h = seq(0.2, 0.8, 0.2), col = "gray88", lty = 3)

  # Bootstrap CI ribbon
  polygon(c(N, rev(N)),
          c(gt_result$boot_ci_lo, rev(gt_result$boot_ci_hi)),
          col    = adjustcolor("darkorange", alpha.f = 0.15),
          border = NA)

  # Lines
  lines(N, gt_result$boot_power_mean,
        col = "darkorange", lwd = 2, lty = 2)
  lines(N, gt_result$ts_power_mean,
        col = "steelblue",  lwd = 2, lty = 2)
  lines(N, gt_result$true_power_mean,
        col = "black", lwd = 3, lty = 1)

  # Oracle SE band, represents Monte Carlo variability from finite nsims,
  # NOT inferential uncertainty about the true power. Rendered very faintly
  # to avoid confusion with the bootstrap CI (which is inferential).
  polygon(c(N, rev(N)),
          c(gt_result$true_power_mean - gt_result$true_power_se,
            rev(gt_result$true_power_mean + gt_result$true_power_se)),
          col = adjustcolor("black", 0.04), border = NA)
  text(max(N), gt_result$true_power_mean[length(N)] + gt_result$true_power_se[length(N)],
       "(MC ±SE)", col = "gray50", cex = 0.65, adj = c(1, -0.3))

  legend("bottomright",
         legend = c("Oracle (true params)",
                    "Two-stage estimate",
                    "Bootstrap mean",
                    "Bootstrap 95% CI",
                    "Oracle ±SE (MC error)"),
         col    = c("black", "steelblue", "darkorange",
                    adjustcolor("darkorange", 0.3),
                    adjustcolor("black", 0.25)),
         lty    = c(1, 2, 2, NA, NA),
         lwd    = c(3, 2, 2, NA, NA),
         pch    = c(NA, NA, NA, 15, 15),
         pt.cex = c(NA, NA, NA, 2, 2),
         cex    = 0.75, bty = "n")

  # -----------------------------------------------------------
  # Panel B: Per-N calibration with CI bars
  # -----------------------------------------------------------
  x   <- seq_len(n_sizes)
  off <- 0.18   # horizontal offset for jittering points

  plot(x, gt_result$true_power_mean * 100,
       type = "n",
       xlim = c(0.5, n_sizes + 0.5),
       ylim = c(0, 100),
       xlab = "N per group",
       ylab = "Power (%)",
       main = sprintf("Calibration by N\n(bootstrap CI coverage = %.0f%%)",
                      100 * gt_result$coverage),
       las = 1, xaxt = "n")
  axis(1, at = x, labels = N)
  abline(h = seq(0, 100, 20), col = "gray88", lty = 3)
  abline(h = target_power * 100, lty = 2, col = "gray40")

  # Bootstrap CI bars
  for (j in seq_len(n_sizes)) {
    lines(c(x[j] + off, x[j] + off),
          c(gt_result$boot_ci_lo[j] * 100, gt_result$boot_ci_hi[j] * 100),
          col = "darkorange", lwd = 2.5)
  }

  # Points
  points(x - off, gt_result$ts_power_mean   * 100,
         col = "steelblue",  pch = 17, cex = 1.3)
  points(x + off, gt_result$boot_power_mean * 100,
         col = "darkorange", pch = 19, cex = 1.3)
  points(x,       gt_result$true_power_mean * 100,
         col = "black",      pch = 15, cex = 1.5)

  legend("bottomright",
         legend = c("Oracle (truth)", "Two-stage", "Bootstrap mean ± CI"),
         col    = c("black", "steelblue", "darkorange"),
         pch    = c(15, 17, 19),
         cex    = 0.85, bty = "n")

  invisible(NULL)
}
