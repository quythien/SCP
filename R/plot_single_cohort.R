# Single-cohort Figure 1: 3-panel power summary
# (add_se_bars is defined in plot_with_se.R; local fallback below for standalone use)
if (!exists("add_se_bars")) {
  add_se_bars <- function(x, y, se, col, bar_width = 0.3) {
    valid <- !is.na(y) & !is.na(se) & se > 1e-6
    if (!any(valid)) return(invisible())
    arrows(x[valid], (y - se)[valid], x[valid], (y + se)[valid],
           angle = 90, code = 3, length = bar_width * 0.15,
           col = col, lwd = 1.2)
  }
}
# Panel A: marginal power vs n, lines = FDR thresholds
# Panel B: true discoveries by r-stratum (lines per n) + overlay of gene distribution
# Panel C: stratified power by r-stratum, lines per n

#' Plot single-cohort circadian power (3 panels)
#'
#' @param res        Output from \code{runSimsSingleCohort()}.
#' @param panel_a_res Optional alternate result object used for the Panel A
#'   marginal-power computation only. If NULL (default), Panel A is computed
#'   from \code{res}. Useful when the marginal summary should reflect a
#'   different empirical assumption than the effect-size-stratified panels.
#' @param out_pdf    Path for PDF output. If NULL, plots to current device.
#' @param title      Overall figure title (default empty).
#' @param fdr_thresholds FDR levels for panel A (default c(0.01,0.05,0.10,0.20)).
#' @param p.adjust.method Correction method (default "BH").
#' @param reference_n    Reference sample size for annotations (default NULL = none).
#' @param width      PDF width in inches (default 15).
#' @param height     PDF height in inches (default 5.5).
#' @return Invisibly returns list of data used per panel.
#' @export
plotSingleCohortPower <- function(res, out_pdf = NULL, title = "",
                                 panel_a_res = NULL,
                                 panels = c("A", "B", "C"),
                                 fdr_thresholds = c(0.01, 0.05, 0.10, 0.20),
                                 panel_fdr = 0.05,
                                 vline_power = 0.80,
                                 vline_fdr   = 0.20,
                                 fdr = NULL,
                                 p.adjust.method = "BH",
                                 reference_n = NULL,
                                 display_sizes = NULL,
                                 r_max = 5,
                                 cex_main = 1.50, cex_lab = 1.60, cex_axis = 1.45,
                                 width = 15, height = 5.5) {
  # `fdr` is a single-value alias for panel_fdr + vline_fdr (back-compat)
  if (!is.null(fdr) && is.numeric(fdr) && length(fdr) == 1L) {
    panel_fdr <- fdr
    vline_fdr <- fdr
  }
  # Adaptive r_max: pick the 95th percentile of observed r-tilde values across
  # all simulations so panels B and C only show populated strata (don't waste
  # space on r in (3, 5] when no gene has r > 1.5).
  if (is.null(r_max)) {
    r_pool <- unlist(lapply(res$r_values_list, function(rv) {
      unlist(lapply(rv, function(v) v[v > 0]), use.names = FALSE)
    }), use.names = FALSE)
    if (length(r_pool) > 0L) {
      r_max <- max(stats::quantile(r_pool, 0.95, na.rm = TRUE), 1.0)
      r_max <- ceiling(r_max * 2) / 2   # round up to nearest 0.5
    } else {
      r_max <- 5
    }
  }

  sample_sizes  <- res$sample_sizes
  nsims         <- res$nsims
  r_strata      <- res$r_strata
  strata_labels <- res$strata_labels
  n_sizes       <- length(sample_sizes)
  n_strata      <- length(strata_labels)
  n_thresh      <- length(fdr_thresholds)

  base_cols     <- c("darkgreen", "steelblue", "orange", "red", "purple", "brown", "magenta")
  thresh_cols   <- if (n_thresh <= length(base_cols)) base_cols[seq_len(n_thresh)]
                   else grDevices::rainbow(n_thresh)
  thresh_labels <- paste0("FDR ", round(100 * fdr_thresholds, 1), "%")
  size_colors   <- rainbow(n_sizes, s = 0.6, v = 0.8)

  # Which sample sizes to show in panels B and C (all by default)
  disp_idx <- if (!is.null(display_sizes))
                which(sample_sizes %in% display_sizes)
              else seq_len(n_sizes)

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

  vline_n <- npower(res, target_power = vline_power, fdr = vline_fdr)$n

  # If a separate result object is supplied for Panel A, override the marginal
  # summary (mean, SE, N target) computed above. The effect-size strata in
  # Panels B and C continue to use the primary `res`.
  if (!is.null(panel_a_res)) {
    pa <- panel_a_res
    pa_nsizes  <- length(pa$sample_sizes)
    pa_marg_sim <- array(NA_real_, dim = c(pa_nsizes, n_thresh, pa$nsims))
    for (jj in seq_len(pa_nsizes)) {
      for (ss in seq_len(pa$nsims)) {
        rv2   <- pa$r_values_list[[jj]][[ss]]
        rhy2  <- rv2 > 0
        n_t2  <- sum(rhy2)
        pv2   <- pa$pvalues[jj, , ss]
        pv2[is.na(pv2)] <- 1
        qv2   <- p.adjust(pv2, method = p.adjust.method)
        for (tt in seq_len(n_thresh)) {
          disc2 <- qv2 <= fdr_thresholds[tt]
          pa_marg_sim[jj, tt, ss] <- if (n_t2 > 0) sum(disc2 & rhy2) / n_t2 else NA_real_
        }
      }
    }
    marginal_mean <- apply(pa_marg_sim, c(1, 2), mean, na.rm = TRUE)
    marginal_se   <- apply(pa_marg_sim, c(1, 2),
                           function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
    sample_sizes_a <- pa$sample_sizes
    disp_idx_a     <- if (!is.null(display_sizes))
                        which(sample_sizes_a %in% display_sizes)
                      else seq_along(sample_sizes_a)
    vline_n <- npower(pa, target_power = vline_power, fdr = vline_fdr)$n
  } else {
    sample_sizes_a <- sample_sizes
    disp_idx_a     <- disp_idx
  }

  # ------------------------------------------------------------------
  # Panels B & C: stratum-level results at panel_fdr (closest in fdr_thresholds)
  # ------------------------------------------------------------------
  idx_fdr5  <- which.min(abs(fdr_thresholds - panel_fdr))
  fdr_label <- sprintf("FDR %g%%", panel_fdr * 100)
  mean_TD <- apply(TD_arr[, , idx_fdr5, , drop = FALSE], c(1, 2),
                   mean, na.rm = TRUE)  # [size, stratum]
  se_TD   <- apply(TD_arr[, , idx_fdr5, , drop = FALSE], c(1, 2),
                   function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

  # Gene distribution by r: average counts across all sims and sample sizes
  count_mat <- matrix(0L, nrow = n_sizes * nsims, ncol = n_strata)
  row_i <- 1L
  for (j in seq_len(n_sizes)) {
    for (s in seq_len(nsims)) {
      rv  <- res$r_values_list[[j]][[s]]
      rhy <- rv > 0
      xg  <- cut(rv[rhy], breaks = r_strata, include.lowest = TRUE, labels = FALSE)
      count_mat[row_i, ] <- tabulate(xg, nbins = n_strata)
      row_i <- row_i + 1L
    }
  }
  gene_counts <- colMeans(count_mat)

  # ------------------------------------------------------------------
  # Panel C: stratified power at FDR 5% by r-stratum, lines per n
  # ------------------------------------------------------------------
  mean_pow <- apply(power_arr[, , idx_fdr5, , drop = FALSE], c(1, 2),
                    mean, na.rm = TRUE)  # [size, stratum]
  se_pow   <- apply(power_arr[, , idx_fdr5, , drop = FALSE], c(1, 2),
                    function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

  # ------------------------------------------------------------------
  # Collapse strata with r >= r_max into a single ">r_max" bin
  # ------------------------------------------------------------------
  left_bounds <- r_strata[-length(r_strata)]   # left boundary of each bin
  collapse_from <- which(left_bounds >= r_max)

  if (length(collapse_from) >= 2) {
    keep <- seq_len(min(collapse_from) - 1)

    gene_counts_plt <- c(gene_counts[keep], sum(gene_counts[collapse_from]))

    mean_TD_plt <- cbind(mean_TD[, keep, drop = FALSE],
                         rowSums(mean_TD[, collapse_from, drop = FALSE]))
    se_TD_plt   <- cbind(se_TD[, keep, drop = FALSE],
                         sqrt(rowSums(se_TD[, collapse_from, drop = FALSE]^2)))

    gc_merged <- sum(gene_counts[collapse_from])
    if (gc_merged > 0) {
      pow_merged <- rowSums(mean_TD[, collapse_from, drop = FALSE]) / gc_merged
    } else {
      pow_merged <- rep(0, n_sizes)
    }
    mean_pow_plt <- cbind(mean_pow[, keep, drop = FALSE], pow_merged)
    se_pow_plt   <- cbind(se_pow[, keep, drop = FALSE],
                          matrix(NA_real_, n_sizes, 1))

    strata_labels_plt <- c(strata_labels[keep], sprintf(">%g", r_max))
    n_strata_plt      <- length(strata_labels_plt)
  } else {
    gene_counts_plt   <- gene_counts
    mean_TD_plt       <- mean_TD;  se_TD_plt  <- se_TD
    mean_pow_plt      <- mean_pow; se_pow_plt <- se_pow
    strata_labels_plt <- strata_labels
    n_strata_plt      <- n_strata
  }

  # ------------------------------------------------------------------
  # Open device
  # ------------------------------------------------------------------
  panels <- intersect(c("A", "B", "C"), panels)
  if (length(panels) == 0L) panels <- c("A", "B", "C")
  if (!is.null(out_pdf)) pdf(out_pdf, width = width, height = height)
  # Scale left/bottom margins with the label font so the "Power (%)" y-title
  # is not clipped at large cex (e.g. the enlarged Shiny fonts). Anchored at
  # ~1.55 so the manuscript figures (default cex) are essentially unchanged.
  extra_lab <- max(0, cex_lab - 1.55)
  mai_left  <- 0.85 + 1.25 * extra_lab
  mai_bot   <- 0.85 + 0.60 * extra_lab
  # Outer top oma needs ~3 lines for the cex=1.5 title to clear the top edge.
  par(mfrow = c(1, length(panels)), mai = c(mai_bot, mai_left, 0.55, 0.15),
      mgp = c(3.2, 0.65, 0), oma = c(0, 0, 2.4, 0),
      cex.axis = cex_axis, cex.lab = cex_lab, cex.main = cex_main, font.main = 2)
  on.exit({ if (!is.null(out_pdf)) dev.off() }, add = TRUE)

  sample_legend_n <- length(disp_idx)
  sample_legend_cols <- if (sample_legend_n > 24L) 3L else if (sample_legend_n > 12L) 2L else 1L
  sample_legend_cex <- if (sample_legend_n > 24L) 0.72 else if (sample_legend_n > 12L) 0.85 else 1.05

  if ("A" %in% panels) {
  # ---- Panel A: marginal power vs n, multiple FDR lines ----
  ss_disp_a <- sample_sizes_a[disp_idx_a]
  matplot(ss_disp_a, 100 * marginal_mean[disp_idx_a, , drop = FALSE],
          type = "b", pch = 19, lwd = 2,
          col = thresh_cols, lty = 1,
          xlim = c(0, max(ss_disp_a) * 1.05), ylim = c(0, 100),
          xlab = "Sample size (n)", ylab = "Power (%)",
          main = "")
  title(main = "A   Power vs Sample Size",
        adj = 0.5, font.main = 2, cex.main = cex_main, line = 0.5)
  for (t in seq_len(n_thresh)) {
    add_se_bars(ss_disp_a, 100 * marginal_mean[disp_idx_a, t],
                100 * marginal_se[disp_idx_a, t], col = thresh_cols[t])
  }
  if (!is.null(reference_n)) {
    abline(v = reference_n, lty = 3, col = "darkgreen", lwd = 1.5)
  }
  abline(h = 80, lty = 2, col = "grey50", lwd = 1.2)
  grid()
  legend("topleft", title = "FDR", legend = thresh_labels,
         col = thresh_cols, lty = 1, pch = 19, lwd = 2,
         cex = 1.05, bty = "o", bg = "white", box.col = "grey50", box.lwd = 0.8,
         text.col = "black", title.col = "black",
         inset = c(0.02, 0.02), y.intersp = 1.0, title.adj = 0,
         seg.len = 1.4, x.intersp = 0.7)
  if (!is.na(vline_n)) {
    abline(v = vline_n, lty = 2, col = adjustcolor("steelblue", 0.7), lwd = 1.5)
    # Place "n=NNN" to the RIGHT of the dashed line, low in the plot (y=35).
    # Above the bottom-right legend (~y=0-25) but below the curves' high-N
    # asymptote, so it is clearly visible in white space.
    text(vline_n, 35, sprintf("n=%d", vline_n),
         col = "steelblue", cex = 1.05, adj = c(-0.10, 0.5), font = 2)
  }
  }  # end Panel A

  if ("B" %in% panels) {
  # ---- Panel B: stratified power by r, lines per n ----
  # Increase mgp[1] so r̃ label clears the vertical (las=2) tick labels
  par(mgp = c(4.6, 0.6, 0))
  matplot(seq_len(n_strata_plt), 100 * t(mean_pow_plt[disp_idx, , drop = FALSE]),
          type = "l", lwd = 2, col = size_colors[disp_idx], lty = 1,
          xlim = c(0.5, n_strata_plt + 0.5), ylim = c(0, 100), bty = "l",
          xlab = expression(tilde(r) == A/sigma), ylab = "Power (%)",
          main = "",
          xaxt = "n")
  title(main = bquote(bold("B   ") * bold("Stratified Power by") ~
                       bold(tilde(r)) ~ bold(.(sprintf("(%s)", fdr_label)))),
        adj = 0.5, font.main = 2, cex.main = cex_main, line = 0.5)
  axis(1, at = seq_len(n_strata_plt), labels = strata_labels_plt, las = 2, cex.axis = max(1.15, cex_axis * 0.8))
  for (j in disp_idx) {
    points(seq_len(n_strata_plt), 100 * mean_pow_plt[j, ],
           pch = 19, col = size_colors[j], cex = 0.95)
    add_se_bars(seq_len(n_strata_plt), 100 * mean_pow_plt[j, ],
                100 * se_pow_plt[j, ], col = size_colors[j])
  }
  grid()
  legend("topleft", title = "Sample size", legend = paste0("n = ", sample_sizes[disp_idx]),
         col = size_colors[disp_idx], lty = 1, lwd = 2,
         cex = sample_legend_cex, bty = "o", bg = "white",
         box.col = "grey50", box.lwd = 0.8,
         text.col = "black", title.col = "black",
         inset = c(0.02, 0.02), y.intersp = 0.95, title.adj = 0,
         seg.len = 1.2, x.intersp = 0.6, ncol = sample_legend_cols)
  }  # end Panel B

  if ("C" %in% panels) {
  # ---- Panel C: TD by r-stratum, n lines + overlaid gene distribution ----
  y_max_TD      <- max(gene_counts_plt) * 1.15
  scaled_counts <- gene_counts_plt

  plot(seq_len(n_strata_plt), rep(0, n_strata_plt),
       type = "n", bty = "l",
       xlim = c(0.5, n_strata_plt + 0.5), ylim = c(0, y_max_TD),
       xlab = expression(tilde(r) == A/sigma), ylab = "# True Discoveries",
       main = "",
       xaxt = "n")
  title(main = bquote(bold("C   ") * bold("True Discoveries by") ~
                       bold(tilde(r)) ~ bold(.(sprintf("(%s)", fdr_label)))),
        adj = 0.5, font.main = 2, cex.main = cex_main, line = 0.5)
  axis(1, at = seq_len(n_strata_plt), labels = strata_labels_plt, las = 2, cex.axis = max(1.15, cex_axis * 0.8))

  step_x <- rep(seq(0.5, n_strata_plt + 0.5, by = 1), each = 2)
  step_y <- c(0, rep(scaled_counts, each = 2), 0)
  polygon(step_x, step_y, col = "#cccccc55", border = NA)
  lines(step_x, step_y, col = "grey60", lwd = 1.5, lty = 2)

  for (j in disp_idx) {
    lines(seq_len(n_strata_plt), mean_TD_plt[j, ], col = size_colors[j], lwd = 2)
    points(seq_len(n_strata_plt), mean_TD_plt[j, ], pch = 19, col = size_colors[j], cex = 0.95)
    add_se_bars(seq_len(n_strata_plt), mean_TD_plt[j, ], se_TD_plt[j, ], col = size_colors[j])
  }
  grid()
  legend("topleft", title = "Sample size",
         legend = c(paste0("n = ", sample_sizes[disp_idx]), "Target gene count"),
         col = c(size_colors[disp_idx], "grey60"), lty = c(rep(1, length(disp_idx)), 2),
         lwd = c(rep(2, length(disp_idx)), 1.5),
         cex = sample_legend_cex, bty = "o", bg = "white",
         box.col = "grey50", box.lwd = 0.8,
         text.col = "black", title.col = "black",
         inset = c(0.02, 0.02), y.intersp = 0.95, title.adj = 0,
         seg.len = 1.2, x.intersp = 0.6, ncol = sample_legend_cols)
  par(mgp = c(3.0, 0.6, 0))   # restore
  }  # end Panel C

  if (nchar(title) > 0)
    mtext(title, outer = TRUE, side = 3, line = 0.2, cex = 1.65, font = 2)

  invisible(list(marginal_mean = marginal_mean, marginal_se = marginal_se,
                 mean_TD = mean_TD, mean_pow = mean_pow))
}
