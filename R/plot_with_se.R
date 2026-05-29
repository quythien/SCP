#' =======================================================================
#' Plotting Functions with Standard Error Bars
#' =======================================================================
#'
#' SE plotting functions for DR/DP stratified power and phase shift analyses.
#' SE = sd(power) / sqrt(nsims) computed across simulation replicates.
#'
#' Functions:
#'   add_se_bars()          - Add +/-1 SE error bars to a plot
#'   plotWithSE()           - 6-panel DR/DP power figure with SE bars
#'   plotPhaseShiftWithSE() - 6-panel phase shift figure with SE bars


# =====================================================================
# HELPER: add error bars to a plot
# =====================================================================
#' Add +/-1 SE Error Bars to a Base-R Plot
#'
#' @param x Numeric vector. x-coordinates of the bar centres.
#' @param y Numeric vector. y-coordinates (point estimates).
#' @param se Numeric vector. Standard errors (bar half-width).
#' @param col Colour passed to \code{arrows()}.
#' @param bar_width Numeric. Relative width of the bar caps (default 0.3).
#'
#' @return Invisibly returns \code{NULL}. Called for its side-effect of
#'   drawing arrows on the active plot device.
add_se_bars <- function(x, y, se, col, bar_width = 0.3) {
  valid <- !is.na(y) & !is.na(se) & se > 0
  if (!any(valid)) return(invisible())
  arrows(x[valid], (y - se)[valid], x[valid], (y + se)[valid],
         angle = 90, code = 3, length = bar_width * 0.15,
         col = col, lwd = 1.2)
}


# =====================================================================
# Plot DR/DP with SE bars (replaces plotDRPowerStratified)
# =====================================================================
#' @title Plot differential power with standard-error bands
#' @noRd
#' @param results_file Path to saved .rds results (from run_pipeline.R)
#' @param output_file  Path for output PDF
#' @param test_name    "DR" or "DP"
#' @return Invisible list with marginal_mean and marginal_se
plotWithSE <- function(results_file, output_file, test_name = "DR",
                       analysis.opts = NULL) {

  load(results_file)
  if (exists("dp_power_raw")) dr_power_raw <- dp_power_raw

  sample_sizes  <- dr_power_raw$sample_sizes
  nsims         <- dr_power_raw$nsims
  strata_labels <- dr_power_raw$strata_labels
  r_strata      <- dr_power_raw$r_strata

  # Use analysis.opts if provided, otherwise defaults
  if (!is.null(analysis.opts) && inherits(analysis.opts, "CircadianAnalysisOptions")) {
    fdr_thresholds <- analysis.opts$fdr_thresholds
    reference_n    <- analysis.opts$reference_n
  } else {
    fdr_thresholds <- c(0.01, 0.05, 0.10, 0.20)
    reference_n    <- 60
  }
  n_sizes      <- length(sample_sizes)
  n_thresholds <- length(fdr_thresholds)
  n_strata     <- length(strata_labels)

  # Detect nested ground truth
  nested_gt <- is.list(dr_power_raw$is_target_list[[1]]) && length(dr_power_raw$is_target_list) == n_sizes

  # ---------------------------------------------------------------
  # Recompute power at each FDR threshold: [size, stratum, threshold, sim]
  # ---------------------------------------------------------------
  power_arr <- array(NA, dim = c(n_sizes, n_strata, n_thresholds, nsims))
  TD_arr    <- array(NA, dim = c(n_sizes, n_strata, n_thresholds, nsims))
  FD_arr    <- array(NA, dim = c(n_sizes, n_strata, n_thresholds, nsims))

  for (j in 1:n_sizes) {
    for (s in 1:nsims) {
      pvals <- dr_power_raw$pvalues[j, , s]

      if (nested_gt) {
        r_vec     <- dr_power_raw$r_values_list[[j]][[s]]
        is_target <- dr_power_raw$is_target_list[[j]][[s]]
        is_null   <- dr_power_raw$is_null_list[[j]][[s]]
      } else {
        r_vec     <- dr_power_raw$r_values_list[[s]]
        is_target <- dr_power_raw$is_target_list[[s]]
        is_null   <- dr_power_raw$is_null_list[[s]]
      }

      xgr <- cut(r_vec, breaks = r_strata, include.lowest = TRUE, labels = FALSE)

      qvals <- rep(1, length(pvals))
      tested <- pvals < 1
      if (sum(tested) > 0) qvals[tested] <- p.adjust(pvals[tested], method = "BH")

      for (t in 1:n_thresholds) {
        discoveries <- qvals <= fdr_thresholds[t]
        for (k in 1:n_strata) {
          in_stratum <- xgr == k
          td <- sum(discoveries & is_target & in_stratum, na.rm = TRUE)
          fd <- sum(discoveries & is_null & in_stratum, na.rm = TRUE)
          nt <- sum(is_target & in_stratum, na.rm = TRUE)
          TD_arr[j, k, t, s] <- td
          FD_arr[j, k, t, s] <- fd
          power_arr[j, k, t, s] <- if (nt > 0) td / nt else NA
        }
      }
    }
  }

  # ---------------------------------------------------------------
  # Per-simulation marginal power: [size, threshold, sim]
  # ---------------------------------------------------------------
  marginal_sim <- array(NA, dim = c(n_sizes, n_thresholds, nsims))
  TD_sim       <- array(NA, dim = c(n_sizes, n_thresholds, nsims))

  for (j in 1:n_sizes) {
    for (s in 1:nsims) {
      if (nested_gt) {
        total_targets_s <- sum(dr_power_raw$is_target_list[[j]][[s]], na.rm = TRUE)
      } else {
        total_targets_s <- sum(dr_power_raw$is_target_list[[1]], na.rm = TRUE)
      }
      for (t in 1:n_thresholds) {
        td_total <- sum(TD_arr[j, , t, s], na.rm = TRUE)
        TD_sim[j, t, s] <- td_total
        marginal_sim[j, t, s] <- if (total_targets_s > 0) td_total / total_targets_s else NA
      }
    }
  }

  # Means and SEs
  marginal_mean <- apply(marginal_sim, c(1, 2), mean, na.rm = TRUE)
  marginal_se   <- apply(marginal_sim, c(1, 2), function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
  TD_mean       <- apply(TD_sim, c(1, 2), mean, na.rm = TRUE)
  TD_se         <- apply(TD_sim, c(1, 2), function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

  # ---------------------------------------------------------------
  # GENERATE 6-PANEL FIGURE
  # ---------------------------------------------------------------
  pdf(output_file, width = 16, height = 10)
  par(mfrow = c(2, 3), mai = c(0.9, 1.0, 0.5, 0.2), mgp = c(3, 0.5, 0), oma = c(0, 0, 2, 0))

  threshold_colors <- c("darkgreen", "steelblue", "orange", "red")
  threshold_labels <- c("FDR 1%", "FDR 5%", "FDR 10%", "FDR 20%")
  size_colors <- rainbow(n_sizes, s = 0.6, v = 0.8)

  # --- Panel A: Power by r (multiple n lines) at FDR 5% ---
  mean_power_A <- apply(dr_power_raw$strat_power, c(1, 2), mean, na.rm = TRUE)
  se_power_A   <- apply(dr_power_raw$strat_power, c(1, 2), function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

  matplot(1:n_strata, t(100 * mean_power_A),
          type = "l", lwd = 2, col = size_colors, lty = 1,
          xlim = c(0.5, n_strata + 0.5), ylim = c(0, 100),
          xlab = "r = A/sigma", ylab = "Stratified Power (%)",
          main = sprintf("%s Power by r (FDR 5%%)", test_name), xaxt = "n")
  axis(1, at = 1:n_strata, labels = strata_labels, las = 2, cex.axis = 0.6)
  for (j in 1:n_sizes) {
    points(1:n_strata, 100 * mean_power_A[j, ], pch = 19, col = size_colors[j], cex = 0.6)
    add_se_bars(1:n_strata, 100 * mean_power_A[j, ], 100 * se_power_A[j, ], col = size_colors[j])
  }
  abline(h = 80, lty = 2, col = "gray")
  grid()
  legend("bottomright", paste0("n=", sample_sizes), col = size_colors, lty = 1, lwd = 2, cex = 0.6)
  mtext("A", side = 3, at = -0.02, font = 2)

  # --- Panel B: Power by r at different FDR thresholds (reference n) ---
  idx_n60 <- which.min(abs(sample_sizes - reference_n))
  mean_power_B <- matrix(NA, nrow = n_strata, ncol = n_thresholds)
  se_power_B   <- matrix(NA, nrow = n_strata, ncol = n_thresholds)
  for (t in 1:n_thresholds) {
    for (k in 1:n_strata) {
      vals <- power_arr[idx_n60, k, t, ]
      mean_power_B[k, t] <- mean(vals, na.rm = TRUE)
      se_power_B[k, t]   <- sd(vals, na.rm = TRUE) / sqrt(sum(!is.na(vals)))
    }
  }

  matplot(1:n_strata, 100 * mean_power_B,
          type = "b", pch = 19, lwd = 2, col = threshold_colors, lty = 1,
          xlim = c(0.5, n_strata + 0.5), ylim = c(0, 100),
          xlab = "r = A/sigma", ylab = "Stratified Power (%)",
          main = sprintf("%s Power by r at Different FDR (n=%d)", test_name, sample_sizes[idx_n60]),
          xaxt = "n")
  axis(1, at = 1:n_strata, labels = strata_labels, las = 2, cex.axis = 0.6)
  for (t in 1:n_thresholds) {
    add_se_bars(1:n_strata, 100 * mean_power_B[, t], 100 * se_power_B[, t], col = threshold_colors[t])
  }
  abline(h = 80, lty = 2, col = "gray")
  grid()
  legend("bottomright", threshold_labels, col = threshold_colors, lty = 1, pch = 19, lwd = 2, cex = 0.7)
  mtext("B", side = 3, at = -0.02, font = 2)

  # --- Panel C: Marginal Power vs Sample Size (FDR 5%, clean) ---
  idx_fdr5 <- 2
  plot(sample_sizes, 100 * marginal_mean[, idx_fdr5],
       type = "b", pch = 19, lwd = 2, col = "steelblue",
       xlim = c(0, max(sample_sizes) * 1.1), ylim = c(0, 100),
       xlab = "Sample Size (per group)", ylab = "Marginal Power (%)",
       main = sprintf("%s Marginal Power vs Sample Size (FDR 5%%)", test_name))
  add_se_bars(sample_sizes, 100 * marginal_mean[, idx_fdr5], 100 * marginal_se[, idx_fdr5], col = "steelblue")
  abline(h = 80, lty = 2, col = "gray")
  abline(v = 60, lty = 3, col = "darkgreen", lwd = 1.5)
  grid()
  text(sample_sizes, 100 * marginal_mean[, idx_fdr5] + 4,
       sprintf("%.1f%%", 100 * marginal_mean[, idx_fdr5]), cex = 0.65, col = "steelblue")
  text(63, 5, "n=60", cex = 0.6, col = "darkgreen", pos = 4)
  mtext("C", side = 3, at = -0.02, font = 2)

  # --- Panel D: Marginal Power vs Sample Size (all FDR thresholds) ---
  matplot(sample_sizes, 100 * marginal_mean,
          type = "b", pch = 19, lwd = 2, col = threshold_colors, lty = 1,
          xlim = c(0, max(sample_sizes) * 1.1), ylim = c(0, 100),
          xlab = "Sample Size (per group)", ylab = "Marginal Power (%)",
          main = sprintf("%s Power vs Sample Size by FDR", test_name))
  for (t in 1:n_thresholds) {
    add_se_bars(sample_sizes, 100 * marginal_mean[, t], 100 * marginal_se[, t], col = threshold_colors[t])
  }
  abline(h = 80, lty = 2, col = "gray")
  grid()
  legend("bottomright", threshold_labels, col = threshold_colors, lty = 1, pch = 19, lwd = 2, cex = 0.7)
  mtext("D", side = 3, at = -0.02, font = 2)

  # --- Panel E: True Discoveries by r stratum (FDR 5%, lines per n) ---
  mean_TD_E <- apply(dr_power_raw$strat_TD, c(1, 2), mean, na.rm = TRUE)
  se_TD_E   <- apply(dr_power_raw$strat_TD, c(1, 2), function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

  y_max_E <- max(mean_TD_E + se_TD_E, na.rm = TRUE) * 1.1
  matplot(1:n_strata, t(mean_TD_E),
          type = "l", lwd = 2, col = size_colors, lty = 1,
          xlim = c(0.5, n_strata + 0.5), ylim = c(0, y_max_E),
          xlab = "r = A/sigma", ylab = "Mean True Discoveries (per sim)",
          main = sprintf("%s True Discoveries by r (FDR 5%%)", test_name), xaxt = "n")
  axis(1, at = 1:n_strata, labels = strata_labels, las = 2, cex.axis = 0.6)
  for (j in 1:n_sizes) {
    points(1:n_strata, mean_TD_E[j, ], pch = 19, col = size_colors[j], cex = 0.6)
    add_se_bars(1:n_strata, mean_TD_E[j, ], se_TD_E[j, ], col = size_colors[j])
  }
  grid()
  legend("topright", paste0("n=", sample_sizes), col = size_colors, lty = 1, lwd = 2, cex = 0.6)
  mtext("E", side = 3, at = -0.02, font = 2)

  # --- Panel F: Distribution of target genes across r strata ---
  if (nested_gt) {
    r_vec     <- dr_power_raw$r_values_list[[1]][[1]]
    is_target <- dr_power_raw$is_target_list[[1]][[1]]
  } else {
    r_vec     <- dr_power_raw$r_values_list[[1]]
    is_target <- dr_power_raw$is_target_list[[1]]
  }
  xgr <- cut(r_vec[is_target], breaks = r_strata, include.lowest = TRUE, labels = FALSE)
  target_counts <- tabulate(xgr, nbins = n_strata)

  bp <- barplot(target_counts, names.arg = strata_labels,
                col = "lightblue", border = "gray40",
                xlab = "r = A/sigma", ylab = "Number of Target Genes",
                main = sprintf("Distribution of %s Target Genes by r", test_name),
                las = 2, cex.names = 0.6)
  mtext("F", side = 3, at = -0.02, font = 2)

  mtext(sprintf("%s Power Analysis - Stratified by Signal-to-Noise Ratio", test_name),
        outer = TRUE, font = 2, cex = 1.2)
  dev.off()

  message(sprintf("Figure saved: %s", output_file))

  # Print summary
  cat(sprintf("\n%s MARGINAL POWER (mean +/- SE across %d simulations)\n", test_name, nsims))
  cat(sprintf("%-10s |", "n"))
  for (t in 1:n_thresholds) cat(sprintf(" %-14s |", threshold_labels[t]))
  cat("\n")
  cat(paste0(rep("-", 15 + n_thresholds * 18), collapse = ""), "\n")
  for (j in 1:n_sizes) {
    cat(sprintf("n = %-6d |", sample_sizes[j]))
    for (t in 1:n_thresholds) {
      cat(sprintf(" %5.1f%% +/- %4.1f%% |", 100 * marginal_mean[j, t], 100 * marginal_se[j, t]))
    }
    cat("\n")
  }
  cat("\n")

  invisible(list(marginal_mean = marginal_mean, marginal_se = marginal_se))
}


# =====================================================================
# Plot Phase Shift with SE bars
# =====================================================================
#' @title Plot phase-shift power with standard-error bands
#' @noRd
#' @param results_file Path to saved phase shift .rds results
#' @param output_file  Path for output PDF
#' @return Invisible list with marginal_mean and marginal_se
plotPhaseShiftWithSE <- function(results_file, output_file,
                                 analysis.opts = NULL) {

  load(results_file)

  sample_sizes  <- dp_phase_results$sample_sizes
  phase_shifts  <- dp_phase_results$phase_shifts
  nsims         <- dp_phase_results$nsims
  r_strata      <- dp_phase_results$r_strata
  strata_labels <- dp_phase_results$strata_labels
  strat_power   <- dp_phase_results$strat_power
  strat_TD      <- dp_phase_results$strat_TD
  strat_n_targets <- dp_phase_results$strat_n_targets

  # Drop 0-hour phase shift (null; undefined power)
  keep_phase <- phase_shifts != 0
  if (any(keep_phase)) {
    phase_shifts <- phase_shifts[keep_phase]
    strat_power <- strat_power[keep_phase, , , , drop = FALSE]
    strat_TD <- strat_TD[keep_phase, , , , drop = FALSE]
    strat_n_targets <- strat_n_targets[keep_phase, , , , drop = FALSE]
  }

  # Use analysis.opts if provided, otherwise defaults
  if (!is.null(analysis.opts) && inherits(analysis.opts, "CircadianAnalysisOptions")) {
    reference_n <- analysis.opts$reference_n
  } else {
    reference_n <- 60
  }

  n_phase    <- length(phase_shifts)
  n_size     <- length(sample_sizes)
  n_r_strata <- length(r_strata) - 1

  idx_n60 <- which.min(abs(sample_sizes - reference_n))
  show_phase_idx    <- which(phase_shifts %in% c(2, 4, 6, 8, 10, 12))
  show_phase_labels <- phase_shifts[show_phase_idx]

  show_r_idx    <- c(4, 6, 8, 10, 12, 14)
  show_r_labels <- strata_labels[show_r_idx]
  show_r_colors <- rainbow(6, s = 0.6, v = 0.8)

  # --- Pre-compute means and SEs ---

  # Panel A: Power vs r by Phase Shift (at n=60)
  mean_power_n60 <- array(NA, dim = c(n_phase, n_r_strata))
  se_power_n60   <- array(NA, dim = c(n_phase, n_r_strata))
  for (p in 1:n_phase) {
    for (k in 1:n_r_strata) {
      vals <- strat_power[p, idx_n60, k, ]
      mean_power_n60[p, k] <- mean(vals, na.rm = TRUE)
      se_power_n60[p, k]   <- sd(vals, na.rm = TRUE) / sqrt(sum(!is.na(vals)))
    }
  }

  # Panel C: Power vs Phase Shift by r stratum (at n=60)
  power_by_phase    <- matrix(NA, nrow = n_phase, ncol = length(show_r_idx))
  se_by_phase       <- matrix(NA, nrow = n_phase, ncol = length(show_r_idx))
  for (p in 1:n_phase) {
    for (ri in seq_along(show_r_idx)) {
      k <- show_r_idx[ri]
      vals <- strat_power[p, idx_n60, k, ]
      power_by_phase[p, ri] <- mean(vals, na.rm = TRUE)
      se_by_phase[p, ri]    <- sd(vals, na.rm = TRUE) / sqrt(sum(!is.na(vals)))
    }
  }

  # Panel D: Per-simulation marginal power: [size, phase, sim]
  marginal_sim <- array(NA, dim = c(n_size, n_phase, nsims))
  for (j in 1:n_size) {
    for (p in 1:n_phase) {
      for (s in 1:nsims) {
        td  <- sum(strat_TD[p, j, , s], na.rm = TRUE)
        tgt <- sum(strat_n_targets[p, j, , s], na.rm = TRUE)
        marginal_sim[j, p, s] <- if (tgt > 0) td / tgt else NA
      }
    }
  }
  marginal_mean <- apply(marginal_sim, c(1, 2), mean, na.rm = TRUE)
  marginal_se   <- apply(marginal_sim, c(1, 2), function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

  # --- Create 6-panel PDF ---

  pdf(output_file, width = 16, height = 10)
  par(mfrow = c(2, 3), mai = c(0.9, 1.0, 0.5, 0.2), mgp = c(3, 0.5, 0), oma = c(0, 0, 2, 0))

  phase_colors <- rainbow(length(show_phase_idx), s = 0.6, v = 0.8)

  # --- Panel A: Power vs r by Phase Shift (n=60) ---
  matplot(1:n_r_strata, 100 * t(mean_power_n60[show_phase_idx, ]),
          type = "b", pch = 19, lwd = 2,
          col = phase_colors, lty = 1,
          xlim = c(0.5, n_r_strata + 0.5), ylim = c(0, 100),
          xlab = "r = A/sigma (Signal-to-Noise Ratio)",
          ylab = "Power (%)",
          main = "Power vs r by Phase Shift (n=60)",
          xaxt = "n")
  axis(1, at = 1:n_r_strata, labels = strata_labels, las = 2, cex.axis = 0.6)
  for (ii in seq_along(show_phase_idx)) {
    p <- show_phase_idx[ii]
    add_se_bars(1:n_r_strata, 100 * mean_power_n60[p, ], 100 * se_power_n60[p, ], col = phase_colors[ii])
  }
  abline(h = 80, lty = 2, col = "gray")
  grid()
  legend("bottomright", paste0(show_phase_labels, "h"),
         col = phase_colors, lty = 1, pch = 19, lwd = 2, cex = 0.7)
  mtext("A", side = 3, at = -0.02, font = 2, line = 0.5)

  # --- Panel B: Heatmap of Power (phase shift x r) at n=60 ---
  par(mai = c(0.9, 1.0, 0.5, 0.5))
  z_matrix <- 100 * t(mean_power_n60)
  image(z_matrix,
        xlab = "Phase Shift (hours)", ylab = "r = A/sigma",
        main = "Power Heatmap (n=60)",
        col = colorRampPalette(c("white", "yellow", "orange", "red", "darkred"))(100),
        xaxt = "n", yaxt = "n")
  axis(1, at = seq(0, 1, length.out = n_phase), labels = phase_shifts)
  axis(2, at = seq(0, 1, length.out = n_r_strata), labels = strata_labels, las = 2, cex.axis = 0.6)
  box()
  mtext("B", side = 3, at = -0.02, font = 2, line = 0.5)

  # --- Panel C: Power vs Phase Shift by r stratum (n=60) ---
  par(mai = c(0.9, 1.0, 0.5, 0.2))
  matplot(phase_shifts, 100 * power_by_phase,
          type = "b", pch = 19, lwd = 2,
          col = show_r_colors, lty = 1,
          xlim = c(0, 13), ylim = c(0, 100),
          xlab = "Phase Shift (hours)",
          ylab = "Power (%)",
          main = "Power vs Phase Shift by r (n=60)")
  for (ri in seq_along(show_r_idx)) {
    add_se_bars(phase_shifts, 100 * power_by_phase[, ri], 100 * se_by_phase[, ri], col = show_r_colors[ri])
  }
  abline(h = 80, lty = 2, col = "gray")
  abline(v = 6, lty = 3, col = "darkgreen", lwd = 1.5)
  grid()
  legend("bottomright", show_r_labels,
         col = show_r_colors, lty = 1, pch = 19, lwd = 2, cex = 0.6)
  mtext("C", side = 3, at = -0.02, font = 2, line = 0.5)

  # --- Panel D: Marginal Power vs Sample Size by Phase Shift ---
  par(mai = c(0.9, 1.0, 0.5, 0.2))
  matplot(sample_sizes, 100 * marginal_mean[, show_phase_idx],
          type = "b", pch = 19, lwd = 2,
          col = phase_colors, lty = 1,
          xlim = c(0, 110), ylim = c(0, 100),
          xlab = "Sample Size (per group)",
          ylab = "Marginal Power (%)",
          main = "Marginal Power vs Sample Size")
  for (ii in seq_along(show_phase_idx)) {
    p <- show_phase_idx[ii]
    add_se_bars(sample_sizes, 100 * marginal_mean[, p], 100 * marginal_se[, p], col = phase_colors[ii])
  }
  abline(h = 80, lty = 2, col = "gray")
  abline(v = 60, lty = 3, col = "darkgreen", lwd = 1.5)
  grid()
  legend("bottomright", paste0(show_phase_labels, "h"),
         col = phase_colors, lty = 1, pch = 19, lwd = 2, cex = 0.6)
  mtext("D", side = 3, at = -0.02, font = 2, line = 0.5)

  # --- Panel E: Heatmap of Power (sample size x phase shift) ---
  par(mai = c(0.9, 1.0, 0.5, 0.5))
  z_matrix2 <- 100 * t(marginal_mean)
  image(z_matrix2,
        xlab = "Sample Size (per group)", ylab = "Phase Shift (hours)",
        main = "Power Heatmap (Marginal)",
        col = colorRampPalette(c("white", "yellow", "orange", "red", "darkred"))(100),
        xaxt = "n", yaxt = "n")
  axis(1, at = seq(0, 1, length.out = n_size), labels = sample_sizes)
  axis(2, at = seq(0, 1, length.out = n_phase), labels = phase_shifts)
  box()
  mtext("E", side = 3, at = -0.02, font = 2, line = 0.5)

  # --- Panel F: Minimum detectable phase shift by r-stratum (at reference n) ---
  # For each r-stratum, find the smallest phase shift achieving >=80% power at n=reference_n
  par(mai = c(0.9, 1.0, 0.5, 0.2))

  min_shift <- rep(NA, n_r_strata)
  for (k in 1:n_r_strata) {
    for (p in 1:n_phase) {
      if (phase_shifts[p] == 0) next
      pwr <- mean(strat_power[p, idx_n60, k, ], na.rm = TRUE)
      if (!is.na(pwr) && pwr >= 0.80) {
        # Interpolate between this and previous phase shift
        if (p > 1) {
          pwr_prev <- mean(strat_power[p - 1, idx_n60, k, ], na.rm = TRUE)
          if (!is.na(pwr_prev) && pwr > pwr_prev) {
            frac <- (0.80 - pwr_prev) / (pwr - pwr_prev)
            min_shift[k] <- phase_shifts[p - 1] + frac * (phase_shifts[p] - phase_shifts[p - 1])
          } else {
            min_shift[k] <- phase_shifts[p]
          }
        } else {
          min_shift[k] <- phase_shifts[p]
        }
        break
      }
    }
  }

  # Only show strata that have at least some target genes (non-NA power somewhere)
  has_data <- sapply(1:n_r_strata, function(k) {
    any(!is.na(strat_power[, idx_n60, k, ]), na.rm = TRUE)
  })

  bar_vals <- min_shift
  bar_vals[is.na(bar_vals)] <- max(phase_shifts)  # cap at 12 (max possible)
  bar_colors_f <- ifelse(is.na(min_shift), "gray80",
                         ifelse(min_shift <= 2, "steelblue",
                                ifelse(min_shift <= 6, "orange", "darkred")))
  bar_colors_f[!has_data] <- "white"

  bp <- barplot(bar_vals, names.arg = strata_labels,
                col = bar_colors_f, border = "gray40",
                ylim = c(0, max(phase_shifts) + 1),
                xlab = "r = A/sigma", ylab = "Min Phase Shift (hours)",
                main = sprintf("Min Detectable Shift for 80%% Power (n=%d)", sample_sizes[idx_n60]),
                las = 2, cex.names = 0.6)

  # Add labels
  for (i in 1:n_r_strata) {
    if (!has_data[i]) next
    if (is.na(min_shift[i])) {
      text(bp[i], bar_vals[i] + 0.3, "N/R", cex = 0.55, font = 3, col = "gray40")
    } else {
      text(bp[i], bar_vals[i] + 0.3, sprintf("%.1fh", min_shift[i]), cex = 0.55, font = 2)
    }
  }

  abline(h = 6, lty = 2, col = "darkgreen", lwd = 1.5)
  text(bp[n_r_strata], 6.3, "6h shift", cex = 0.6, col = "darkgreen", pos = 2)
  abline(h = 2, lty = 3, col = "steelblue", lwd = 1.5)
  text(bp[n_r_strata], 2.3, "2h shift", cex = 0.6, col = "steelblue", pos = 2)
  grid()
  legend("topright",
         c("<= 2h", "2-6h", "> 6h", "Not reached"),
         fill = c("steelblue", "orange", "darkred", "gray80"),
         cex = 0.6, border = "gray40")
  mtext("F", side = 3, at = -0.02, font = 2, line = 0.5)

  mtext("Differential Phase Power Analysis: Effect of Phase Shift Magnitude",
        outer = TRUE, font = 2, cex = 1.2)

  dev.off()
  message(sprintf("Phase shift figure saved: %s", output_file))

  invisible(list(marginal_mean = marginal_mean, marginal_se = marginal_se))
}
