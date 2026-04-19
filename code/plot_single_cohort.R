# Single-cohort Figure 1: 3-panel power summary
# (add_se_bars is defined in plot_with_se.R; local fallback below for standalone use)
if (!exists("add_se_bars")) {
  add_se_bars <- function(x, y, se, col, bar_width = 0.3) {
    valid <- !is.na(y) & !is.na(se) & se > 0
    if (!any(valid)) return(invisible())
    arrows(x[valid], (y - se)[valid], x[valid], (y + se)[valid],
           angle = 90, code = 3, length = bar_width * 0.15,
           col = col, lwd = 1.2)
  }
}
# Panel A: marginal power vs n, lines = FDR thresholds
# Panel B: true discoveries by r-stratum (lines per n) + overlay of gene distribution
# Panel C: stratified power by r-stratum, lines per n

#' Plot single-cohort Figure 1 (3 panels)
#'
#' @param res        Output from \code{runSimsSingleCohort()}.
#' @param out_pdf    Path for PDF output. If NULL, plots to current device.
#' @param title      Overall figure title (default empty).
#' @param fdr_thresholds FDR levels for panel A (default c(0.01,0.05,0.10,0.20)).
#' @param p.adjust.method Correction method (default "BH").
#' @param reference_n    Reference sample size for annotations (default NULL = none).
#' @param width      PDF width in inches (default 15).
#' @param height     PDF height in inches (default 5.5).
#' @return Invisibly returns list of data used per panel.
#' @export
plotSingleCohortFig1 <- function(res, out_pdf = NULL, title = "",
                                 fdr_thresholds = c(0.01, 0.05, 0.10, 0.20),
                                 p.adjust.method = "BH",
                                 reference_n = NULL,
                                 width = 15, height = 5.5) {

  sample_sizes  <- res$sample_sizes
  nsims         <- res$nsims
  r_strata      <- res$r_strata
  strata_labels <- res$strata_labels
  n_sizes       <- length(sample_sizes)
  n_strata      <- length(strata_labels)
  n_thresh      <- length(fdr_thresholds)

  thresh_cols   <- c("darkgreen", "steelblue", "orange", "red")[seq_len(n_thresh)]
  thresh_labels <- paste0("FDR ", round(100 * fdr_thresholds), "%")
  size_colors   <- rainbow(n_sizes, s = 0.6, v = 0.8)

  # ------------------------------------------------------------------
  # Recompute power and TD at each FDR threshold from raw pvalues
  # power_arr[size, stratum, threshold, sim]
  # TD_arr   [size, stratum, threshold, sim]
  # ------------------------------------------------------------------
  power_arr <- array(NA_real_, dim = c(n_sizes, n_strata, n_thresh, nsims))
  TD_arr    <- array(NA_real_, dim = c(n_sizes, n_strata, n_thresh, nsims))

  for (j in seq_len(n_sizes)) {
    for (s in seq_len(nsims)) {
      pvals    <- res$pvalues[j, , s]
      r_values <- res$r_values_list[[j]][[s]]

      is_rhythmic <- r_values > 0   # non-rhythmic genes have r = 0 by construction

      # r-stratum per gene (NA for non-rhythmic)
      r_for_strat <- ifelse(is_rhythmic, r_values, 0)
      xgr <- cut(r_for_strat, breaks = r_strata, include.lowest = TRUE, labels = FALSE)
      xgr[!is_rhythmic] <- NA

      # BH correction
      qvals <- rep(1, length(pvals))
      pvals[is.na(pvals)] <- 1
      qvals <- p.adjust(pvals, method = p.adjust.method)

      for (t in seq_len(n_thresh)) {
        disc <- qvals <= fdr_thresholds[t]
        for (k in seq_len(n_strata)) {
          in_k  <- !is.na(xgr) & xgr == k
          td    <- sum(disc & is_rhythmic & in_k, na.rm = TRUE)
          n_tgt <- sum(is_rhythmic & in_k, na.rm = TRUE)
          TD_arr[j, k, t, s]    <- td
          power_arr[j, k, t, s] <- if (n_tgt > 0) td / n_tgt else NA_real_
        }
      }
    }
  }

  # ------------------------------------------------------------------
  # Panel A: marginal power vs n, lines = FDR thresholds
  # marginal_sim[size, threshold, sim]
  # ------------------------------------------------------------------
  marginal_sim <- array(NA_real_, dim = c(n_sizes, n_thresh, nsims))
  for (j in seq_len(n_sizes)) {
    for (s in seq_len(nsims)) {
      r_values    <- res$r_values_list[[j]][[s]]
      is_rhythmic <- r_values > 0
      n_tgt       <- sum(is_rhythmic)
      for (t in seq_len(n_thresh)) {
        td_total <- sum(TD_arr[j, , t, s], na.rm = TRUE)
        marginal_sim[j, t, s] <- if (n_tgt > 0) td_total / n_tgt else NA_real_
      }
    }
  }
  marginal_mean <- apply(marginal_sim, c(1, 2), mean, na.rm = TRUE)
  marginal_se   <- apply(marginal_sim, c(1, 2),
                         function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

  # ------------------------------------------------------------------
  # Panel B: TD by r-stratum at FDR 5% (use idx_fdr5 = closest to 0.05)
  # ------------------------------------------------------------------
  idx_fdr5 <- which.min(abs(fdr_thresholds - 0.05))
  mean_TD <- apply(TD_arr[, , idx_fdr5, , drop = FALSE], c(1, 2),
                   mean, na.rm = TRUE)  # [size, stratum]
  se_TD   <- apply(TD_arr[, , idx_fdr5, , drop = FALSE], c(1, 2),
                   function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

  # Gene distribution by r (use first sim of first n for representative counts)
  r_sample  <- res$r_values_list[[1]][[1]]
  is_rhy    <- r_sample > 0
  xgr_dist  <- cut(r_sample[is_rhy], breaks = r_strata,
                   include.lowest = TRUE, labels = FALSE)
  gene_counts <- tabulate(xgr_dist, nbins = n_strata)

  # ------------------------------------------------------------------
  # Panel C: stratified power at FDR 5% by r-stratum, lines per n
  # ------------------------------------------------------------------
  mean_pow <- apply(power_arr[, , idx_fdr5, , drop = FALSE], c(1, 2),
                    mean, na.rm = TRUE)  # [size, stratum]
  se_pow   <- apply(power_arr[, , idx_fdr5, , drop = FALSE], c(1, 2),
                    function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

  # ------------------------------------------------------------------
  # Open device
  # ------------------------------------------------------------------
  if (!is.null(out_pdf)) pdf(out_pdf, width = width, height = height)
  par(mfrow = c(1, 3), mai = c(0.9, 1.0, 0.5, 0.2),
      mgp = c(3, 0.5, 0), oma = c(0, 0, 2, 0))
  on.exit({ if (!is.null(out_pdf)) dev.off() }, add = TRUE)

  # ---- Panel A: marginal power vs n, multiple FDR lines ----
  matplot(sample_sizes, 100 * marginal_mean,
          type = "b", pch = 19, lwd = 2,
          col = thresh_cols, lty = 1,
          xlim = c(0, max(sample_sizes) * 1.05), ylim = c(0, 100),
          xlab = "Sample size (n)", ylab = "Marginal Power (%)",
          main = "Marginal Power vs Sample Size")
  for (t in seq_len(n_thresh)) {
    add_se_bars(sample_sizes, 100 * marginal_mean[, t],
                100 * marginal_se[, t], col = thresh_cols[t])
  }
  abline(h = 80, lty = 2, col = "gray")
  if (!is.null(reference_n)) {
    abline(v = reference_n, lty = 3, col = "darkgreen", lwd = 1.5)
  }
  grid()
  legend("bottomright", thresh_labels,
         col = thresh_cols, lty = 1, pch = 19, lwd = 2, cex = 0.7)
  mtext("A", side = 3, at = par("usr")[1], font = 2, line = 0.5)

  # ---- Panel B: TD by r-stratum, n lines + overlaid gene distribution ----
  y_max_TD <- max(mean_TD + se_TD, na.rm = TRUE) * 1.15
  # Scale gene_counts to share axis
  count_scale <- y_max_TD / max(gene_counts + 1)
  scaled_counts <- gene_counts * count_scale

  plot(seq_len(n_strata), rep(0, n_strata),
       type = "n",
       xlim = c(0.5, n_strata + 0.5), ylim = c(0, y_max_TD),
       xlab = "r = A/sigma", ylab = "Mean True Discoveries (per sim)",
       main = "True Discoveries by r (FDR 5%)", xaxt = "n")
  axis(1, at = seq_len(n_strata), labels = strata_labels, las = 2, cex.axis = 0.6)

  # Smoothed gene-count overlay (loess on stratum index, scaled to TD axis)
  lo_counts <- tryCatch(
    predict(loess(scaled_counts ~ seq_len(n_strata), span = 0.5),
            newdata = data.frame(seq_len(n_strata))),
    error = function(e) scaled_counts
  )
  lo_counts <- pmax(lo_counts, 0)
  polygon(c(seq_len(n_strata), rev(seq_len(n_strata))),
          c(lo_counts, rep(0, n_strata)),
          col = "#cccccc55", border = NA)
  lines(seq_len(n_strata), lo_counts, col = "grey60", lwd = 1.5, lty = 2)

  # Secondary axis for gene counts
  axis(4, at = pretty(c(0, y_max_TD), 4),
       labels = round(pretty(c(0, y_max_TD), 4) / count_scale),
       las = 1, cex.axis = 0.65, col.axis = "grey50")

  for (j in seq_len(n_sizes)) {
    lines(seq_len(n_strata), mean_TD[j, ], col = size_colors[j], lwd = 2)
    points(seq_len(n_strata), mean_TD[j, ], pch = 19, col = size_colors[j], cex = 0.6)
    add_se_bars(seq_len(n_strata), mean_TD[j, ], se_TD[j, ], col = size_colors[j])
  }
  grid()
  legend("topright",
         c(paste0("n=", sample_sizes), "Gene density"),
         col = c(size_colors, "grey60"), lty = c(rep(1, n_sizes), 2),
         lwd = c(rep(2, n_sizes), 1.5), cex = 0.6)
  mtext("B", side = 3, at = par("usr")[1], font = 2, line = 0.5)

  # ---- Panel C: stratified power by r, lines per n ----
  matplot(seq_len(n_strata), 100 * t(mean_pow),
          type = "l", lwd = 2, col = size_colors, lty = 1,
          xlim = c(0.5, n_strata + 0.5), ylim = c(0, 100),
          xlab = "r = A/sigma", ylab = "Power (%)",
          main = "Power by r (FDR 5%)", xaxt = "n")
  axis(1, at = seq_len(n_strata), labels = strata_labels, las = 2, cex.axis = 0.6)
  for (j in seq_len(n_sizes)) {
    points(seq_len(n_strata), 100 * mean_pow[j, ],
           pch = 19, col = size_colors[j], cex = 0.6)
    add_se_bars(seq_len(n_strata), 100 * mean_pow[j, ],
                100 * se_pow[j, ], col = size_colors[j])
  }
  abline(h = 80, lty = 2, col = "gray")
  grid()
  legend("bottomright", paste0("n=", sample_sizes),
         col = size_colors, lty = 1, lwd = 2, cex = 0.6)
  mtext("C", side = 3, at = par("usr")[1], font = 2, line = 0.5)

  if (nchar(title) > 0)
    mtext(title, outer = TRUE, cex = 1.0, font = 2)

  invisible(list(marginal_mean = marginal_mean, marginal_se = marginal_se,
                 mean_TD = mean_TD, mean_pow = mean_pow))
}
