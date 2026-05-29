#' DR Power Analysis Plotting Functions
#'
#' Functions for generating stratified power plots
#' from saved DR simulation results

# (First copy removed - see plotDRPowerStratified below)


#' Quick Replot DR Power (convenience wrapper)
#'
#' @param test_name "DR" or "DP"
#' @param results_dir Directory containing results files
#' @param output_file Path to save PDF (optional)
#'
#' @return Invisible list with calculated power statistics
replotDRPower <- function(test_name = "DR",
                          results_dir = "output/dr_power_stratified",
                          output_file = NULL) {

  results_file <- file.path(results_dir, sprintf("%s_power_raw_pvalues.rds", tolower(test_name)))

  if (!file.exists(results_file)) {
    stop(sprintf("Results file not found: %s\n", results_file))
  }

  if (is.null(output_file)) {
    output_file <- file.path("output/figures", sprintf("%s_power_stratified.pdf", tolower(test_name)))
  }

  plotDRPowerStratified(results_file, output_file, test_name)
}


#' Generate Sampling Density Plot
#'
#' @param results_file Path to saved density results
#' @param output_file Path to save PDF figure
#'
#' @return Invisible density results
plotDensityResults <- function(results_file, output_file = NULL) {

  load(results_file)

  if (is.null(output_file)) {
    output_file <- sprintf("output/figures/%s_sampling_density.pdf", tolower(density_results$test_type))
  }

  n_subjects <- density_results$n_subjects
  n_time_points <- density_results$n_time_points
  power_avg <- density_results$power_avg
  TD_avg <- density_results$TD_avg
  FDC_avg <- density_results$FDC_avg
  test_type <- density_results$test_type

  pdf(output_file, width = 14, height = 10)
  par(mfrow = c(2, 2), mai = c(0.9, 1.0, 0.6, 0.3), mgp = c(3, 0.5, 0))

  colors <- rainbow(length(n_time_points), s = 0.6, v = 0.8)

  # Panel A: Power vs Sample Size for different time point densities
  matplot(n_subjects, 100 * power_avg,
          type = "b", pch = 19, lwd = 2,
          col = colors, lty = 1,
          xlim = c(0, max(n_subjects) * 1.1), ylim = c(0, 100),
          xlab = "Number of Subjects (per group)",
          ylab = "Power (%)",
          main = sprintf("%s Power: Subjects vs Time Points", test_type))
  abline(h = 80, lty = 2, col = "gray")
  grid()
  legend("bottomright", paste0(n_time_points, " time points"),
         col = colors, lty = 1, pch = 19, lwd = 2, cex = 0.7)
  mtext("A", side = 3, at = -0.02, font = 2)

  # Panel B: True Discoveries vs Sample Size
  matplot(n_subjects, TD_avg,
          type = "b", pch = 19, lwd = 2,
          col = colors, lty = 1,
          xlim = c(0, max(n_subjects) * 1.1),
          ylim = c(0, max(TD_avg) * 1.1),
          xlab = "Number of Subjects (per group)",
          ylab = "Mean # True Discoveries",
          main = sprintf("%s Discoveries: Subjects vs Time Points", test_type))
  grid()
  mtext("B", side = 3, at = -0.02, font = 2)

  # Panel C: Power Heatmap
  par(mai = c(0.9, 1.0, 0.6, 0.5))
  z_matrix <- 100 * t(power_avg)
  image(z_matrix,
        xlab = "Number of Subjects", ylab = "Time Points",
        main = sprintf("%s Power Heatmap (%%)", test_type),
        col = colorRampPalette(c("white", "yellow", "orange", "red", "darkred"))(100),
        xaxt = "n", yaxt = "n")
  axis(1, at = seq(0, 1, length.out = length(n_subjects)), labels = n_subjects)
  axis(2, at = seq(0, 1, length.out = length(n_time_points)), labels = n_time_points)
  box()
  mtext("C", side = 3, at = -0.02, font = 2)

  # Panel D: FDC Heatmap
  par(mai = c(0.9, 1.0, 0.6, 0.5))
  z_matrix2 <- t(FDC_avg)
  z_matrix2[is.na(z_matrix2) | is.infinite(z_matrix2)] <- 5
  image(z_matrix2,
        xlab = "Number of Subjects", ylab = "Time Points",
        main = sprintf("%s FDC Heatmap", test_type),
        col = colorRampPalette(c("white", "yellow", "orange", "red", "darkred"))(100),
        xaxt = "n", yaxt = "n")
  axis(1, at = seq(0, 1, length.out = length(n_subjects)), labels = n_subjects)
  axis(2, at = seq(0, 1, length.out = length(n_time_points)), labels = n_time_points)
  box()
  mtext("D", side = 3, at = -0.02, font = 2)

  dev.off()

  message(sprintf("Figure saved: %s", output_file))

  invisible(density_results)
}


#' Print Density Analysis Summary Table
#'
#' @description
#' Prints a formatted text table of power, true discoveries, false discoveries,
#' and false discovery cost from a density analysis result, stratified by sample
#' size and number of time points.
#'
#' @param density_results List. Output of \code{runSimsDensity()} containing
#'   \code{n_subjects}, \code{n_time_points}, \code{power_avg}, \code{TD_avg},
#'   \code{FD_avg}, \code{FDC_avg}, and \code{test_type}.
#'
#' @return Invisibly returns \code{NULL}. Called for its side-effect of printing
#'   to the console.
printDensitySummary <- function(density_results) {

  n_subjects <- density_results$n_subjects
  n_time_points <- density_results$n_time_points
  power_avg <- density_results$power_avg
  TD_avg <- density_results$TD_avg
  FD_avg <- density_results$FD_avg
  FDC_avg <- density_results$FDC_avg

  cat("\n====================================================================\n")
  cat(sprintf("%s POWER: SUBJECTS VS TIME POINTS\n", density_results$test_type))
  cat("====================================================================\n\n")

  for (j in seq_along(n_time_points)) {
    nt <- n_time_points[j]
    cat(sprintf("Time points: %d\n", nt))
    cat(sprintf("%-10s | %-8s | %-8s | %-8s | %-8s\n",
                "Subjects", "Power", "TD", "FD", "FDC"))
    cat(paste0(rep("-", 50), collapse = ""), "\n")

    for (i in seq_along(n_subjects)) {
      cat(sprintf("%-10d | %-8.1f%% | %-8.1f | %-8.1f | %-8.3f\n",
                  n_subjects[i],
                  100 * power_avg[i, j],
                  TD_avg[i, j],
                  FD_avg[i, j],
                  FDC_avg[i, j]))
    }
    cat("\n")
  }
}
#' Plot DR Power Stratified by r-Stratum
#'
#' @description
#' Loads a saved power analysis result and produces a multi-panel plot showing
#' detection power stratified by signal-to-noise ratio (r = A/sigma) for each
#' sample size and time-point design.
#'
#' @param results_file Character. Path to a saved \code{.RData} results file
#'   from a power analysis run (must contain \code{dr_power_raw} or
#'   \code{dp_power_raw}).
#' @param output_file Character or NULL. Path for the output PDF. If NULL,
#'   the plot is drawn to the current device.
#' @param test_name Character. Label for the endpoint being plotted (default
#'   \code{"DR"}; can also be \code{"DP"}).
#'
#' @return Invisibly returns \code{NULL}. Called for its side-effect of
#'   producing a plot.
plotDRPowerStratified <- function(results_file, output_file = NULL, test_name = "DR") {

  # Load saved results
  load(results_file)
  
  # Handle both dr_power_raw and dp_power_raw naming
  if (exists("dp_power_raw")) {
    dr_power_raw <- dp_power_raw
  }
  
  sample_sizes <- dr_power_raw$sample_sizes
  nsims <- dr_power_raw$nsims
  strata_labels <- dr_power_raw$strata_labels
  r_strata <- dr_power_raw$r_strata

  # Set default output filename if not provided
  if (is.null(output_file)) {
    output_file <- sprintf("output/figures/%s_power_stratified.pdf", tolower(test_name))
  }

  # Calculate power at FDR thresholds from p-values
  fdr_thresholds <- c(0.01, 0.05, 0.10, 0.20)
  n_sizes <- length(sample_sizes)
  n_thresholds <- length(fdr_thresholds)
  n_strata <- length(strata_labels)

  # Detect ground truth structure: nested [[j]][[s]] or flat [[s]]
  # Nested: each sample size has its own ground truth (independent simulations)
  # Flat: all sample sizes share the same ground truth (legacy format)
  nested_gt <- is.list(dr_power_raw$is_target_list[[1]]) && length(dr_power_raw$is_target_list) == n_sizes

  # Apply FDR correction on-the-fly and calculate power at each threshold
  # Storage: [sample_size, stratum, threshold, sim]
  power_by_threshold <- list(
    power = array(NA, dim = c(n_sizes, n_strata, n_thresholds, nsims)),
    TD = array(NA, dim = c(n_sizes, n_strata, n_thresholds, nsims)),
    FD = array(NA, dim = c(n_sizes, n_strata, n_thresholds, nsims))
  )

  for (j in 1:n_sizes) {
    for (s in 1:nsims) {
      pvals <- dr_power_raw$pvalues[j, , s]

      # Get ground truth for this sample size and simulation
      if (nested_gt) {
        r_vec <- dr_power_raw$r_values_list[[j]][[s]]
        is_target <- dr_power_raw$is_target_list[[j]][[s]]
        is_null <- dr_power_raw$is_null_list[[j]][[s]]
      } else {
        r_vec <- dr_power_raw$r_values_list[[s]]
        is_target <- dr_power_raw$is_target_list[[s]]
        is_null <- dr_power_raw$is_null_list[[s]]
      }

      xgr <- cut(r_vec, breaks = r_strata, include.lowest = TRUE, labels = FALSE)

      # FDR correction: only on tested genes (p < 1), matching runner.R logic
      qvals <- rep(1, length(pvals))
      tested <- pvals < 1
      if (sum(tested) > 0) {
        qvals[tested] <- p.adjust(pvals[tested], method = "BH")
      }

      for (t in 1:n_thresholds) {
        discoveries <- qvals <= fdr_thresholds[t]

        for (k in 1:n_strata) {
          in_stratum <- xgr == k
          power_by_threshold$TD[j, k, t, s] <- sum(discoveries & is_target & in_stratum, na.rm = TRUE)
          power_by_threshold$FD[j, k, t, s] <- sum(discoveries & is_null & in_stratum, na.rm = TRUE)
          n_targets <- sum(is_target & in_stratum, na.rm = TRUE)
          power_by_threshold$power[j, k, t, s] <- if (n_targets > 0) {
            power_by_threshold$TD[j, k, t, s] / n_targets
          } else NA
        }
      }
    }
  }

  # Calculate marginal power per sample size
  # Total targets may differ per j if ground truth is nested (independent sims per n)
  marginal_power <- matrix(NA, nrow = n_sizes, ncol = n_thresholds)
  TD_total <- matrix(NA, nrow = n_sizes, ncol = n_thresholds)

  for (j in 1:n_sizes) {
    # Count total targets for this sample size
    if (nested_gt) {
      total_targets_j <- sum(sapply(1:nsims, function(s) sum(dr_power_raw$is_target_list[[j]][[s]], na.rm = TRUE)))
    } else {
      total_targets_j <- sum(dr_power_raw$is_target_list[[1]], na.rm = TRUE) * nsims
    }
    for (t in 1:n_thresholds) {
      TD_total[j, t] <- sum(power_by_threshold$TD[j, , t, ], na.rm = TRUE)
      marginal_power[j, t] <- if (total_targets_j > 0) TD_total[j, t] / total_targets_j else NA
    }
  }

  # Generate 6-panel figure (2x3)
  pdf(output_file, width = 16, height = 10)
  par(mfrow = c(2, 3), mai = c(0.9, 1.0, 0.5, 0.2), mgp = c(3, 0.5, 0), oma = c(0, 0, 2, 0))

  threshold_colors <- c("darkgreen", "steelblue", "orange", "red")
  threshold_labels <- c("FDR 1%", "FDR 5%", "FDR 10%", "FDR 20%")
  size_colors <- rainbow(n_sizes, s = 0.6, v = 0.8)

  #---------------------------------------------------------------------------
  # Panel A: Power by r (multiple n lines) at FDR 5%
  #---------------------------------------------------------------------------
  mean_power <- apply(dr_power_raw$strat_power, c(1, 2), mean, na.rm = TRUE)
  matplot(1:length(strata_labels), t(100 * mean_power),
          type = "l", lwd = 2, col = size_colors, lty = 1,
          xlim = c(0.5, length(strata_labels) + 0.5), ylim = c(0, 100),
          xlab = "r = A/sigma", ylab = "Stratified Power (%)",
          main = sprintf("%s Power by r (FDR 5%%)", test_name), xaxt = "n")
  axis(1, at = 1:length(strata_labels), labels = strata_labels, las = 2, cex.axis = 0.6)
  abline(h = 80, lty = 2, col = "gray")
  grid()
  legend("bottomright", paste0("n=", sample_sizes), col = size_colors, lty = 1, lwd = 2, cex = 0.6)
  mtext("A", side = 3, at = -0.02, font = 2)

  #---------------------------------------------------------------------------
  # Panel B: Power by r at different FDR thresholds (n=60)
  #---------------------------------------------------------------------------
  idx_n60 <- which.min(abs(sample_sizes - 60))
  power_by_r_threshold <- matrix(NA, nrow = length(strata_labels), ncol = n_thresholds)
  for (t in 1:n_thresholds) {
    for (k in 1:length(strata_labels)) {
      power_by_r_threshold[k, t] <- mean(power_by_threshold$power[idx_n60, k, t, ], na.rm = TRUE)
    }
  }

  matplot(1:length(strata_labels), 100 * power_by_r_threshold,
          type = "b", pch = 19, lwd = 2, col = threshold_colors, lty = 1,
          xlim = c(0.5, length(strata_labels) + 0.5), ylim = c(0, 100),
          xlab = "r = A/sigma", ylab = "Stratified Power (%)",
          main = sprintf("%s Power by r at Different FDR (n=%d)", test_name, sample_sizes[idx_n60]),
          xaxt = "n")
  axis(1, at = 1:length(strata_labels), labels = strata_labels, las = 2, cex.axis = 0.6)
  abline(h = 80, lty = 2, col = "gray")
  grid()
  legend("bottomright", threshold_labels, col = threshold_colors, lty = 1, pch = 19, lwd = 2, cex = 0.7)
  mtext("B", side = 3, at = -0.02, font = 2)

  #---------------------------------------------------------------------------
  # Panel C: Marginal Power vs Sample Size (FDR 5% only, clean)
  #---------------------------------------------------------------------------
  idx_fdr5 <- 2  # FDR 5% is the 2nd threshold
  plot(sample_sizes, 100 * marginal_power[, idx_fdr5],
       type = "b", pch = 19, lwd = 2, col = "steelblue",
       xlim = c(0, max(sample_sizes) * 1.1), ylim = c(0, 100),
       xlab = "Sample Size (per group)", ylab = "Marginal Power (%)",
       main = sprintf("%s Marginal Power vs Sample Size (FDR 5%%)", test_name))
  abline(h = 80, lty = 2, col = "gray")
  abline(v = 60, lty = 3, col = "darkgreen", lwd = 1.5)
  grid()
  # Annotate power values
  text(sample_sizes, 100 * marginal_power[, idx_fdr5] + 4,
       sprintf("%.1f%%", 100 * marginal_power[, idx_fdr5]), cex = 0.65, col = "steelblue")
  text(63, 5, "n=60", cex = 0.6, col = "darkgreen", pos = 4)
  mtext("C", side = 3, at = -0.02, font = 2)

  #---------------------------------------------------------------------------
  # Panel D: Marginal Power vs Sample Size (all FDR thresholds)
  #---------------------------------------------------------------------------
  matplot(sample_sizes, 100 * marginal_power,
          type = "b", pch = 19, lwd = 2, col = threshold_colors, lty = 1,
          xlim = c(0, max(sample_sizes) * 1.1), ylim = c(0, 100),
          xlab = "Sample Size (per group)", ylab = "Marginal Power (%)",
          main = sprintf("%s Power vs Sample Size by FDR", test_name))
  abline(h = 80, lty = 2, col = "gray")
  grid()
  legend("bottomright", threshold_labels, col = threshold_colors, lty = 1, pch = 19, lwd = 2, cex = 0.7)
  mtext("D", side = 3, at = -0.02, font = 2)

  #---------------------------------------------------------------------------
  # Panel E: True Discoveries vs Sample Size
  #---------------------------------------------------------------------------
  matplot(sample_sizes, TD_total,
          type = "b", pch = 19, lwd = 2, col = threshold_colors, lty = 1,
          xlim = c(0, max(sample_sizes) * 1.1), ylim = c(0, max(TD_total) * 1.1),
          xlab = "Sample Size (per group)", ylab = "Total True Discoveries",
          main = sprintf("%s True Discoveries vs Sample Size", test_name))
  grid()
  legend("bottomright", threshold_labels, col = threshold_colors, lty = 1, pch = 19, lwd = 2, cex = 0.7)
  mtext("E", side = 3, at = -0.02, font = 2)

  #---------------------------------------------------------------------------
  # Panel F: Distribution of target genes across r strata
  #---------------------------------------------------------------------------
  if (nested_gt) {
    r_vec <- dr_power_raw$r_values_list[[1]][[1]]
    is_target <- dr_power_raw$is_target_list[[1]][[1]]
  } else {
    r_vec <- dr_power_raw$r_values_list[[1]]
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

  # Overall title
  mtext(sprintf("%s Power Analysis - Stratified by Signal-to-Noise Ratio", test_name),
        outer = TRUE, font = 2, cex = 1.2)

  dev.off()

  # Print summary table
  cat("\n====================================================================\n")
  cat(sprintf("%s MARGINAL POWER AT DIFFERENT FDR THRESHOLDS\n", test_name))
  cat("(FDR applied only to tested genes, matching runner.R logic)\n")
  cat("====================================================================\n\n")

  cat(sprintf("%-10s |", "n"))
  for (t in 1:n_thresholds) cat(sprintf(" %-8s |", threshold_labels[t]))
  cat("\n")
  cat(paste0(rep("-", 15 + n_thresholds * 12), collapse = ""), "\n")

  for (j in 1:n_sizes) {
    cat(sprintf("n = %-7d |", sample_sizes[j]))
    for (t in 1:n_thresholds) cat(sprintf(" %7.1f%% |", 100 * marginal_power[j, t]))
    cat("\n")
  }

  cat(sprintf("\nFigure saved: %s\n", output_file))
  cat("Note: FDR correction applied only to tested genes (p < 1), consistent with runner.R.\n\n")

  invisible(list(
    marginal_power = marginal_power,
    TD = TD_total,
    power_by_threshold = power_by_threshold,
    thresholds = fdr_thresholds,
    sample_sizes = sample_sizes
  ))
}
