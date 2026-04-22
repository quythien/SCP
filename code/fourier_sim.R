#' Fourier Deviation Simulation
#'
#' Tests how power degrades when the true waveform has higher harmonics
#' (2nd, 3rd) that DCP's cosinor model doesn't fit.
#'
#' The cosinor model assumes:
#'   y(t) = M + A*cos(ωt - ωφ) + ε
#'
#' But the true waveform may be:
#'   y(t) = M + A*[cos(ωt - ωφ) + α₂*cos(2ωt - ωφ) + α₃*cos(3ωt - ωφ)] + ε
#'
#' where α₂, α₃ are relative harmonic amplitudes. This module sweeps over
#' (α₂, α₃) combinations and measures resulting power degradation.


# =====================================================================
# Main Fourier Deviation Power Runner
# =====================================================================

#' Run Fourier deviation power analysis
#'
#' Measures how power degrades as higher harmonics are added to the true waveform
#' while the detection model remains a pure cosinor (fundamental only).
#'
#' @param bio.opts      CircadianBioOptions
#' @param design.opts   CircadianDesignOptions
#' @param analysis.opts CircadianAnalysisOptions
#' @param harmonic_grid data.frame with columns alpha2, alpha3 (default: 4x3 grid)
#' @param test_type     "DR", "DP", or "DM"
#' @param verbose       Print progress
#'
#' @return List with power arrays and summaries across harmonic levels
runFourierDeviationPower <- function(bio.opts,
                                     design.opts,
                                     analysis.opts,
                                     harmonic_grid = NULL,
                                     test_type     = "DR",
                                     verbose       = TRUE,
                                     mc.cores      = 1L) {

  stopifnot(inherits(bio.opts,      "CircadianBioOptions"))
  stopifnot(inherits(design.opts,   "CircadianDesignOptions"))
  stopifnot(inherits(analysis.opts, "CircadianAnalysisOptions"))

  if (is.null(harmonic_grid)) {
    harmonic_grid <- expand.grid(
      alpha2 = c(0, 0.25, 0.5, 0.75),
      alpha3 = c(0, 0.25, 0.5)
    )
  }

  harmonic_grid <- as.data.frame(harmonic_grid)
  stopifnot(all(c("alpha2", "alpha3") %in% names(harmonic_grid)))

  # Add readable label
  harmonic_grid$label <- sprintf("(%.2g, %.2g)", harmonic_grid$alpha2, harmonic_grid$alpha3)

  n_harm        <- nrow(harmonic_grid)
  sample_sizes  <- design.opts$sample_sizes
  nsims         <- design.opts$nsims
  n_sizes       <- length(sample_sizes)
  fdr_threshold <- min(analysis.opts$fdr_thresholds)

  if (verbose) {
    cat(sprintf("\n=== Fourier Deviation Power Analysis ===\n"))
    cat(sprintf("  test_type:    %s\n", test_type))
    cat(sprintf("  sample_sizes: %s\n", paste(sample_sizes, collapse = ", ")))
    cat(sprintf("  harmonic combinations: %d\n", n_harm))
    cat(sprintf("  nsims:        %d\n", nsims))
    cat(sprintf("  FDR threshold: %.0f%%\n", 100 * fdr_threshold))
    cat(sprintf("  mc.cores:     %d\n", mc.cores))
  }

  # Flatten (h_idx, s_idx) into a single combo grid for parallel dispatch
  combo_grid <- expand.grid(h_idx = seq_len(n_harm), s_idx = seq_len(n_sizes))

  run_one_combo <- function(ci) {
    h_idx <- combo_grid$h_idx[ci]
    s_idx <- combo_grid$s_idx[ci]
    alpha2        <- harmonic_grid$alpha2[h_idx]
    alpha3        <- harmonic_grid$alpha3[h_idx]
    harmonics_vec <- c(alpha2, alpha3)
    n             <- sample_sizes[s_idx]

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

    tryCatch({
      sim_out <- runSimsDiff(bio.opts, iter_design, analysis.opts,
                             harmonics = harmonics_vec)
      fdr_key <- paste0("fdr_", test_type)
      fdr_arr <- sim_out[[fdr_key]]
      if (is.null(fdr_arr)) return(rep(NA_real_, nsims))

      vapply(seq_len(nsims), function(sim_i) {
        fdr_vec   <- fdr_arr[, 1, sim_i]
        diff_type <- sim_out$diff_type[[sim_i]]
        target_idx <- switch(test_type,
          DR = diff_type %in% c(2, 3),
          DP = diff_type == 4,
          DM = diff_type == 5,
          rep(FALSE, length(fdr_vec))
        )
        n_target <- sum(target_idx)
        if (n_target == 0) NA_real_
        else sum(fdr_vec[target_idx] <= fdr_threshold, na.rm = TRUE) / n_target
      }, numeric(1))

    }, error = function(e) {
      if (verbose) warning(sprintf("  Harmonic %d, N=%d failed: %s", h_idx, n, e$message))
      rep(NA_real_, nsims)
    })
  }

  combo_results <- parallel::mclapply(
    seq_len(nrow(combo_grid)), run_one_combo,
    mc.cores = mc.cores, mc.set.seed = TRUE
  )

  # Reassemble into power_arr[harm_idx, size_idx, sim_idx]
  power_arr <- array(NA_real_,
    dim      = c(n_harm, n_sizes, nsims),
    dimnames = list(harmonic_grid$label, paste0("N", sample_sizes), NULL)
  )
  for (ci in seq_len(nrow(combo_grid))) {
    power_arr[combo_grid$h_idx[ci], combo_grid$s_idx[ci], ] <- combo_results[[ci]]
  }

  power_mean <- apply(power_arr, c(1, 2), mean, na.rm = TRUE)
  power_se   <- apply(power_arr, c(1, 2), sd,   na.rm = TRUE)

  if (verbose) cat("Fourier deviation analysis complete.\n")

  list(
    harmonic_grid = harmonic_grid,
    sample_sizes  = sample_sizes,
    test_type     = test_type,
    fdr_threshold = fdr_threshold,
    power_arr     = power_arr,
    power_mean    = power_mean,
    power_se      = power_se,
    reference_n   = analysis.opts$reference_n
  )
}


# =====================================================================
# Plotting
# =====================================================================

#' Plot Fourier deviation power results
#'
#' Panel A: Heatmap of mean power at reference_n over (α₂, α₃) grid.
#' Panel B: Power vs N lines, one per harmonic combination.
#'
#' @param result      Output from runFourierDeviationPower()
#' @param test_type   Test type label (uses result$test_type if NULL)
#' @param reference_n Reference N for Panel A (uses result$reference_n if NULL)
#' @param output_file Path for PDF output (NULL = screen)
plotFourierDeviation <- function(result,
                                 test_type   = NULL,
                                 reference_n = NULL,
                                 panels      = "B",
                                 output_file = NULL) {
  panels <- match.arg(panels, choices = c("A", "B"), several.ok = TRUE)

  test_type   <- test_type   %||% result$test_type
  reference_n <- reference_n %||% result$reference_n

  hg          <- result$harmonic_grid
  sample_sizes <- result$sample_sizes
  power_mean  <- result$power_mean   # [n_harm x n_sizes]
  fdr_thr     <- result$fdr_threshold

  # Find reference_n index (closest)
  ref_idx <- which.min(abs(sample_sizes - reference_n))

  alpha2_vals <- sort(unique(hg$alpha2))
  alpha3_vals <- sort(unique(hg$alpha3))
  n_harm      <- nrow(hg)

  n_panels <- length(panels)
  if (!is.null(output_file)) {
    pdf(output_file, width = if (n_panels == 1) 7 else 12, height = 5)
    on.exit(dev.off())
  }

  old_par <- par(mfrow = c(1, n_panels), mar = c(4.5, 4.5, 3, 1.5))
  on.exit(par(old_par), add = TRUE)

  col_ramp <- colorRampPalette(c("white", "lightyellow", "orange", "red", "darkred"))(100)

  # -----------------------------------------------------------
  # Panel A: Heatmap at reference_n
  # -----------------------------------------------------------
  if (!("A" %in% panels)) {
    # skip heatmap
  } else if (length(alpha2_vals) > 1 && length(alpha3_vals) > 1) {
    # Build 2D power matrix at reference_n
    pm_grid <- matrix(NA_real_, nrow = length(alpha2_vals), ncol = length(alpha3_vals))
    for (h_idx in seq_len(n_harm)) {
      i2 <- which(alpha2_vals == hg$alpha2[h_idx])
      i3 <- which(alpha3_vals == hg$alpha3[h_idx])
      if (length(i2) == 1 && length(i3) == 1) {
        pm_grid[i2, i3] <- power_mean[h_idx, ref_idx]
      }
    }

    image(seq_along(alpha2_vals), seq_along(alpha3_vals), pm_grid * 100,
          col   = col_ramp, zlim = c(0, 100),
          xaxt  = "n", yaxt = "n",
          xlab  = expression(alpha[2] ~ "(2nd harmonic relative amplitude)"),
          ylab  = expression(alpha[3] ~ "(3rd harmonic relative amplitude)"),
          main  = sprintf("Power (%%) at N=%d\n(%s test, FDR <= %.0f%%)",
                          sample_sizes[ref_idx], test_type, 100 * fdr_thr))
    axis(1, at = seq_along(alpha2_vals), labels = alpha2_vals)
    axis(2, at = seq_along(alpha3_vals), labels = alpha3_vals, las = 1)

    # Cell text annotations
    for (i2 in seq_along(alpha2_vals)) {
      for (i3 in seq_along(alpha3_vals)) {
        val <- pm_grid[i2, i3] * 100
        if (!is.na(val)) {
          txt_col <- if (val > 60) "white" else "black"
          text(i2, i3, sprintf("%.0f", val), col = txt_col, cex = 0.85)
        }
      }
    }

    # Contour lines at 20%, 40%, 60%, 80%
    contour(seq_along(alpha2_vals), seq_along(alpha3_vals), pm_grid * 100,
            levels = c(20, 40, 60, 80),
            labels = c("20%", "40%", "60%", "80%"),
            add = TRUE, col = "gray30", lwd = 1.2)
  } else {
    # Fallback: bar chart when grid is 1D
    pm_at_ref <- power_mean[, ref_idx]
    barplot(pm_at_ref * 100,
            names.arg = hg$label,
            col       = col_ramp[round(pm_at_ref * 99) + 1],
            ylim      = c(0, 100),
            xlab      = "Harmonic (α₂, α₃)",
            ylab      = "Power (%)",
            main      = sprintf("Power at N=%d\n(%s test, FDR <= %.0f%%)",
                                sample_sizes[ref_idx], test_type, 100 * fdr_thr),
            las       = 2, cex.names = 0.75)
    abline(h = 80, lty = 2, col = "gray40")
  } # end Panel A

  if (!("B" %in% panels)) return(invisible(NULL))

  # -----------------------------------------------------------
  # Panel B: Power vs N lines, one per harmonic combination
  # -----------------------------------------------------------
  # Identify pure cosinor (α₂=α₃=0) index
  pure_idx <- which(hg$alpha2 == 0 & hg$alpha3 == 0)

  cols_h <- rainbow(n_harm, s = 0.7, v = 0.85)

  plot(sample_sizes, power_mean[1, ],
       type = "n",
       xlim = range(sample_sizes),
       ylim = c(0, min(1, max(power_mean, na.rm = TRUE) * 1.08)),
       xlab = "N per group",
       ylab = "Power",
       main = sprintf("Power vs N by Harmonic\n(%s test, FDR <= %.0f%%)",
                      test_type, 100 * fdr_thr),
       las  = 1, xaxt = "n")
  axis(1, at = sample_sizes)
  abline(h = 0.80, lty = 2, col = "gray40")
  text(min(sample_sizes), 0.82, "80%", col = "gray40", cex = 0.75, adj = 0)
  abline(h = seq(0.2, 0.8, by = 0.2), col = "gray88", lty = 3)

  for (h_idx in seq_len(n_harm)) {
    is_pure <- length(pure_idx) > 0 && h_idx == pure_idx[1]
    lwd_h   <- if (is_pure) 3 else 1.5
    col_h   <- if (is_pure) "black" else cols_h[h_idx]
    lty_h   <- if (is_pure) 1 else 2
    lines(sample_sizes, power_mean[h_idx, ],
          col = col_h, lwd = lwd_h, lty = lty_h)
    points(sample_sizes, power_mean[h_idx, ],
           col = col_h, pch = if (is_pure) 16 else 1, cex = 0.8)
  }

  legend_cols <- cols_h
  legend_lty  <- rep(2, n_harm)
  legend_lwd  <- rep(1.5, n_harm)
  if (length(pure_idx) > 0) {
    legend_cols[pure_idx[1]] <- "black"
    legend_lty[pure_idx[1]]  <- 1
    legend_lwd[pure_idx[1]]  <- 3
  }

  legend("bottomright",
         legend = sprintf("(α₂,α₃)=%s", hg$label),
         col    = legend_cols,
         lty    = legend_lty,
         lwd    = legend_lwd,
         title  = "Harmonics",
         cex    = 0.7, bty = "n")

  invisible(NULL)
}
