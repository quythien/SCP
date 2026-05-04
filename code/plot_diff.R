# Differential power plot: n_ep * n_comps rows x 3 columns
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
# Differential power Figure 2: 6-row x 3-column (18 panels)
# Row order: DR comp1, DR comp2, DP comp1, DP comp2, DM comp1, DM comp2
# Columns: (A) marginal power vs n | (B) TD by r-stratum | (C) power by r-stratum
#
# Input: res_list — named list of two runSimsDiff() outputs (raw format)
#        comp_labels — character(2) labels for the two comparisons

# ----------------------------------------------------------------------
# Internal: convert raw runSimsDiff() output to stratified arrays
# ----------------------------------------------------------------------
#' Prepare Differential Power Data for Stratified Plotting
#'
#' @description
#' Converts a raw \code{runSimsDiff()} result into stratified power and
#' true-discovery arrays indexed by sample size, r-stratum, FDR threshold,
#' and simulation replicate.
#'
#' @param res List. Output of \code{runSimsDiff()} (contains \code{fdr_DR},
#'   \code{fdr_DP}, \code{fdr_DM}, \code{diff_type}, \code{effectsize}).
#' @param ep Character. Endpoint: one of \code{"DR"}, \code{"DP"}, \code{"DM"}.
#' @param fdr_thresholds Numeric vector. FDR thresholds to evaluate power at.
#' @param r_breaks Numeric vector. Breakpoints defining r-strata (e.g.
#'   \code{c(0, 0.5, 1, 2, Inf)}).
#'
#' @return Named list with arrays \code{power_arr}, \code{TD_arr},
#'   \code{marginal_sim}, and matrix \code{count_mat}.
.prepDiffStratified <- function(res, ep, fdr_thresholds, r_breaks) {
  ngenes   <- res$ngenes
  n_sizes  <- length(res$sample_sizes)
  nsims    <- res$nsims
  n_thresh <- length(fdr_thresholds)
  n_strata <- length(r_breaks) - 1

  target_types <- switch(ep, DR = c(2L, 3L), DP = 4L, DM = 5L, integer(0))
  fdr_mat <- res[[paste0("fdr_", ep)]]   # [ngenes, n_sizes, nsims]

  power_arr    <- array(NA_real_, dim = c(n_sizes, n_strata, n_thresh, nsims))
  TD_arr       <- array(NA_real_, dim = c(n_sizes, n_strata, n_thresh, nsims))
  marginal_sim <- array(NA_real_, dim = c(n_sizes, n_thresh, nsims))
  count_mat    <- matrix(0L, nrow = nsims, ncol = n_strata)

  for (s in seq_len(nsims)) {
    dt <- res$diff_type[[s]]
    es <- res$effectsize[[s]]

    is_target <- dt %in% target_types
    n_tgt_total <- sum(is_target)

    # r value per gene for stratification
    r_full <- rep(0, ngenes)
    if (n_tgt_total > 0) {
      r_full[is_target] <- switch(ep,
        DR = pmax(es$DR1[is_target], es$DR2[is_target]),
        DP = pmin(es$DR1[is_target], es$DR2[is_target]),
        DM = pmin(es$DR1[is_target], es$DR2[is_target]),
        rep(0, n_tgt_total)
      )
    }

    xgr_full <- rep(NA_integer_, ngenes)
    if (n_tgt_total > 0)
      xgr_full[is_target] <- cut(r_full[is_target], breaks = r_breaks,
                                 include.lowest = TRUE, labels = FALSE)

    count_mat[s, ] <- tabulate(xgr_full[is_target], nbins = n_strata)

    for (j in seq_len(n_sizes)) {
      fdr_g <- fdr_mat[, j, s]

      for (t in seq_len(n_thresh)) {
        disc <- !is.na(fdr_g) & fdr_g <= fdr_thresholds[t]
        marginal_sim[j, t, s] <- if (n_tgt_total > 0)
                                    sum(disc & is_target, na.rm = TRUE) / n_tgt_total
                                  else NA_real_

        for (k in seq_len(n_strata)) {
          in_k  <- !is.na(xgr_full) & xgr_full == k
          td    <- sum(disc & is_target & in_k, na.rm = TRUE)
          n_k   <- sum(is_target & in_k, na.rm = TRUE)
          TD_arr[j, k, t, s]    <- td
          power_arr[j, k, t, s] <- if (n_k > 0) td / n_k else NA_real_
        }
      }
    }
  }

  gene_counts <- colMeans(count_mat)

  list(
    power_arr    = power_arr,
    TD_arr       = TD_arr,
    marginal_sim = marginal_sim,
    gene_counts  = gene_counts
  )
}

# ----------------------------------------------------------------------
# Internal: collapse strata with left boundary >= threshold into one bin
# ----------------------------------------------------------------------
#' Collapse High-r Tail Strata into a Single Bin
#'
#' @description
#' Merges all r-strata whose left boundary is at or above \code{threshold}
#' into a single "r ≥ threshold" bin. This prevents sparse high-r strata
#' from cluttering stratified differential power plots.
#'
#' @param r_breaks Numeric vector. Stratum breakpoints (length = n_strata + 1).
#' @param strata_labels Character vector. Labels for each stratum.
#' @param gene_counts Integer vector. Gene counts per stratum per simulation.
#' @param mean_TD Numeric vector. Mean true discoveries per stratum.
#' @param se_TD Numeric vector. SE of true discoveries per stratum.
#' @param mean_pow Numeric vector. Mean power per stratum.
#' @param se_pow Numeric vector. SE of power per stratum.
#' @param threshold Numeric. r value above which strata are collapsed (default 5).
#'
#' @return Named list with collapsed vectors: \code{strata_labels},
#'   \code{gene_counts}, \code{mean_TD}, \code{se_TD}, \code{mean_pow},
#'   \code{se_pow}, and scalar \code{n_strata}.
.collapseTail <- function(r_breaks, strata_labels, gene_counts,
                          mean_TD, se_TD, mean_pow, se_pow,
                          threshold = 5) {
  left_bounds   <- r_breaks[-length(r_breaks)]
  collapse_from <- which(left_bounds >= threshold)

  if (length(collapse_from) < 2) {
    return(list(strata_labels = strata_labels, gene_counts = gene_counts,
                mean_TD = mean_TD, se_TD = se_TD,
                mean_pow = mean_pow, se_pow = se_pow,
                n_strata = length(strata_labels)))
  }

  keep <- seq_len(min(collapse_from) - 1)
  gc   <- c(gene_counts[keep], sum(gene_counts[collapse_from]))

  mTD <- cbind(mean_TD[, keep, drop = FALSE],
               rowSums(mean_TD[, collapse_from, drop = FALSE]))
  sTD <- cbind(se_TD[, keep, drop = FALSE],
               sqrt(rowSums(se_TD[, collapse_from, drop = FALSE]^2)))

  gc_merged <- sum(gene_counts[collapse_from])
  pow_m     <- if (gc_merged > 0)
                 rowSums(mean_TD[, collapse_from, drop = FALSE]) / gc_merged
               else rep(0, nrow(mean_TD))
  mPow <- cbind(mean_pow[, keep, drop = FALSE], pow_m)
  sPow <- cbind(se_pow[, keep, drop = FALSE],
                matrix(NA_real_, nrow(mean_pow), 1))

  list(strata_labels = c(strata_labels[keep], sprintf(">%g", threshold)),
       gene_counts   = gc,
       mean_TD       = mTD, se_TD  = sTD,
       mean_pow      = mPow, se_pow = sPow,
       n_strata      = length(keep) + 1L)
}

#' Plot differential circadian power (18 panels: 6 rows x 3 columns)
#'
#' @param res_list    Named list of length 2, each element output of runSimsDiff().
#' @param comp_labels Character(2) — display labels for the two comparisons.
#' @param endpoints   Which endpoints to include (default c("DR","DP","DM")).
#' @param fdr_thresholds FDR levels for panel A (default c(0.01,0.05,0.10,0.20)).
#' @param display_sizes  Sample sizes to show in all panels (NULL = all).
#' @param r_break_width  Width of r-strata bins (default 0.25).
#' @param out_pdf     Path for PDF. If NULL, plots to current device.
#' @param width       PDF width in inches (default 15).
#' @param height      PDF height in inches (default 30).
#' @return Invisibly, list of panel data per comparison/endpoint.
#' @export
plotDiffPower <- function(res_list,
                         comp_labels    = names(res_list),
                         endpoints      = c("DR", "DP", "DM"),
                         fdr_thresholds = c(0.01, 0.05, 0.10, 0.20),
                         panel_fdr      = 0.05,
                         vline_power    = 0.80,
                         vline_fdr      = 0.20,
                         display_sizes  = NULL,
                         r_break_width  = 0.25,
                         r_max          = 5,
                         r_display_max  = NULL,
                         out_pdf        = NULL,
                         width          = 15,
                         height         = 30) {

  n_comps <- length(res_list)
  if (n_comps < 1) stop("res_list must have at least 1 element.")
  if (is.null(comp_labels) || length(comp_labels) != n_comps)
    comp_labels <- paste("Comparison", seq_len(n_comps))

  n_ep     <- length(endpoints)
  n_thresh <- length(fdr_thresholds)
  thresh_cols   <- c("darkgreen", "steelblue", "orange", "red")[seq_len(n_thresh)]
  thresh_labels <- paste0("FDR ", round(100 * fdr_thresholds), "%")

  r_breaks      <- c(seq(0, 5, by = r_break_width), Inf)
  n_strata_raw  <- length(r_breaks) - 1
  strata_labels_raw <- c(
    sprintf("(%g,%g]", head(r_breaks, -2), r_breaks[seq(2, n_strata_raw)]),
    sprintf(">%g", r_breaks[n_strata_raw])
  )

  # Resolve panel_fdr to a named vector over endpoints (scalar = same for all)
  if (length(panel_fdr) == 1) {
    panel_fdr_vec <- setNames(rep(panel_fdr, n_ep), endpoints)
  } else {
    panel_fdr_vec <- panel_fdr[endpoints]
    if (any(is.na(panel_fdr_vec)))
      stop("panel_fdr must be a scalar or named vector with names matching endpoints.")
  }

  # ------------------------------------------------------------------
  # Precompute n80_max per comparison: max n across endpoints to reach 80% at n80_fdr
  # ------------------------------------------------------------------
  n80_by_comp <- lapply(seq_along(res_list), function(ci) {
    res_ci <- res_list[[ci]]
    n80_ep <- vapply(endpoints, function(ep) {
      np <- npower(res_ci, target_power = vline_power, fdr = vline_fdr, endpoint = ep)
      if (is.na(np$n)) NA_real_ else np$n
    }, numeric(1))
    suppressWarnings(max(n80_ep, na.rm = TRUE))
  })

  # ------------------------------------------------------------------
  # Open device: n_ep * n_comps rows x 3 columns
  # ------------------------------------------------------------------
  if (!is.null(out_pdf)) pdf(out_pdf, width = width, height = height)
  par(mfrow = c(n_ep * n_comps, 3),
      mai   = c(0.9, 0.9, 0.45, 0.2),
      mgp   = c(3.0, 0.5, 0),
      oma   = c(0, 0, 2, 0))
  on.exit({ if (!is.null(out_pdf)) dev.off() }, add = TRUE)

  panel_data       <- list()
  panel_letter_idx <- 1L

  for (ep in endpoints) {
    for (ci in seq_along(res_list)) {

      res          <- res_list[[ci]]
      comp_label   <- comp_labels[ci]
      sample_sizes <- res$sample_sizes
      nsims        <- res$nsims
      n_sizes      <- length(sample_sizes)
      size_colors  <- rainbow(n_sizes, s = 0.6, v = 0.8)

      disp_idx <- if (!is.null(display_sizes))
                    which(sample_sizes %in% display_sizes)
                  else seq_len(n_sizes)
      ss_disp  <- sample_sizes[disp_idx]

      # ------------------------------------------------------------------
      # Preprocess: stratified arrays from raw runSimsDiff() output
      # ------------------------------------------------------------------
      ep_fdr    <- panel_fdr_vec[ep]
      idx_fdr5  <- which.min(abs(fdr_thresholds - ep_fdr))
      fdr_label <- sprintf("FDR %g%%", ep_fdr * 100)

      prep <- .prepDiffStratified(res, ep, fdr_thresholds, r_breaks)

      marginal_mean <- apply(prep$marginal_sim, c(1, 2), mean, na.rm = TRUE)
      marginal_se   <- apply(prep$marginal_sim, c(1, 2),
                             function(x) sd(x, na.rm=TRUE) / sqrt(sum(!is.na(x))))

      mean_TD <- apply(prep$TD_arr[, , idx_fdr5, , drop=FALSE], c(1,2),
                       mean, na.rm = TRUE)
      se_TD   <- apply(prep$TD_arr[, , idx_fdr5, , drop=FALSE], c(1,2),
                       function(x) sd(x, na.rm=TRUE) / sqrt(sum(!is.na(x))))
      mean_pow <- apply(prep$power_arr[, , idx_fdr5, , drop=FALSE], c(1,2),
                        mean, na.rm = TRUE)
      se_pow   <- apply(prep$power_arr[, , idx_fdr5, , drop=FALSE], c(1,2),
                        function(x) sd(x, na.rm=TRUE) / sqrt(sum(!is.na(x))))

      coll <- .collapseTail(r_breaks, strata_labels_raw, prep$gene_counts,
                            mean_TD, se_TD, mean_pow, se_pow, threshold = r_max)

      n_strata_plt  <- coll$n_strata
      strata_lbl    <- coll$strata_labels
      gene_counts   <- coll$gene_counts
      mean_TD_plt   <- coll$mean_TD;  se_TD_plt  <- coll$se_TD
      mean_pow_plt  <- coll$mean_pow; se_pow_plt <- coll$se_pow

      # Cap display at r_display_max: drop strata with upper bound > threshold
      if (!is.null(r_display_max)) {
        upper_b <- sapply(strata_lbl, function(lbl) {
          if (grepl("^>", lbl)) Inf
          else as.numeric(sub(".*,([^]]+)].*", "\\1", lbl))
        })
        keep <- which(upper_b <= r_display_max)
        if (length(keep) > 0) {
          strata_lbl   <- strata_lbl[keep]
          gene_counts  <- gene_counts[keep]
          mean_TD_plt  <- mean_TD_plt[, keep, drop=FALSE]
          se_TD_plt    <- se_TD_plt[,  keep, drop=FALSE]
          mean_pow_plt <- mean_pow_plt[, keep, drop=FALSE]
          se_pow_plt   <- se_pow_plt[,  keep, drop=FALSE]
          n_strata_plt <- length(keep)
        }
      }

      row_label <- sprintf("%s — %s", ep, comp_label)

      # ---- Panel A: marginal power vs n, multiple FDR thresholds ----
      matplot(ss_disp, 100 * marginal_mean[disp_idx, , drop = FALSE],
              type = "b", pch = 19, lwd = 2,
              col = thresh_cols, lty = 1,
              xlim = c(0, max(ss_disp) * 1.05), ylim = c(0, 100),
              xlab = "Sample size (n)", ylab = "Power (%)",
              main = sprintf("%s\nGenome-wide Power vs Sample Size", row_label),
              cex.main = 0.75)
      for (t in seq_len(n_thresh)) {
        add_se_bars(ss_disp, 100 * marginal_mean[disp_idx, t],
                    100 * marginal_se[disp_idx, t], col = thresh_cols[t])
      }
      vline_n <- n80_by_comp[[ci]]
      abline(h = 80, lty = 2, col = "grey50", lwd = 1.2)
      if (!is.na(vline_n) && is.finite(vline_n)) {
        abline(v = vline_n, lty = 2, col = adjustcolor("steelblue", 0.7), lwd = 1.5)
        text(vline_n, 15, sprintf("n=%d", vline_n), col = "steelblue", cex = 0.65, adj = -0.1)
      }
      grid()
      legend("bottomright", thresh_labels,
             col = thresh_cols, lty = 1, pch = 19, lwd = 2, cex = 0.6)

      # ---- Panel B: power by r-stratum at panel_fdr, lines per n ----
      matplot(seq_len(n_strata_plt),
              100 * t(mean_pow_plt[disp_idx, , drop = FALSE]),
              type = "l", lwd = 2, col = size_colors[disp_idx], lty = 1,
              xlim = c(0.5, n_strata_plt + 0.5), ylim = c(0, 100), bty = "l",
              xaxt = "n",
              xlab = expression(tilde(r) == A/sigma), ylab = "Power (%)",
              main = bquote("Stratified Power by" ~ tilde(r) ~ .(sprintf("(%s)", fdr_label))),
              cex.main = 0.75)
      axis(1, at = seq_len(n_strata_plt), labels = strata_lbl,
           las = 2, cex.axis = 0.6)
      for (j in disp_idx) {
        points(seq_len(n_strata_plt), 100 * mean_pow_plt[j, ],
               pch = 19, col = size_colors[j], cex = 0.6)
        add_se_bars(seq_len(n_strata_plt), 100 * mean_pow_plt[j, ],
                    100 * se_pow_plt[j, ], col = size_colors[j])
      }
      grid()
      legend("bottomright", paste0("n=", ss_disp),
             col = size_colors[disp_idx], lty = 1, lwd = 2, cex = 0.55)

      # ---- Panel C: TD by r-stratum at panel_fdr, lines per n ----
      # gene_counts and mean_TD are both in units of genes — use the same natural axis.
      y_max_TD      <- max(gene_counts) * 1.15
      scaled_counts <- gene_counts

      plot(seq_len(n_strata_plt), rep(0, n_strata_plt),
           type = "n", bty = "l", xaxt = "n",
           xlim = c(0.5, n_strata_plt + 0.5), ylim = c(0, y_max_TD),
           xlab = expression(tilde(r) == A/sigma), ylab = sprintf("# True %s Discoveries", ep),
           main = bquote("True Discoveries by" ~ tilde(r) ~ .(sprintf("(%s)", fdr_label))),
           cex.main = 0.75)
      axis(1, at = seq_len(n_strata_plt), labels = strata_lbl,
           las = 2, cex.axis = 0.6)

      step_x <- rep(seq(0.5, n_strata_plt + 0.5, by = 1), each = 2)
      step_y <- c(0, rep(scaled_counts, each = 2), 0)
      polygon(step_x, step_y, col = "#cccccc55", border = NA)
      lines(step_x, step_y, col = "grey60", lwd = 1.5, lty = 2)

      for (j in disp_idx) {
        lines(seq_len(n_strata_plt), mean_TD_plt[j, ], col = size_colors[j], lwd = 2)
        points(seq_len(n_strata_plt), mean_TD_plt[j, ], pch = 19,
               col = size_colors[j], cex = 0.6)
        add_se_bars(seq_len(n_strata_plt), mean_TD_plt[j, ],
                    se_TD_plt[j, ], col = size_colors[j])
      }
      grid()
      legend("topright",
             c(paste0("n=", ss_disp), "# Target Discoveries"),
             col = c(size_colors[disp_idx], "grey60"),
             lty = c(rep(1, length(disp_idx)), 2),
             lwd = c(rep(2, length(disp_idx)), 1.5), cex = 0.55)

      panel_data[[paste(ep, ci, sep = "_")]] <- list(
        marginal_mean = marginal_mean, marginal_se = marginal_se,
        mean_TD  = mean_TD_plt,  se_TD  = se_TD_plt,
        mean_pow = mean_pow_plt, se_pow = se_pow_plt,
        strata_labels = strata_lbl, gene_counts = gene_counts
      )
      panel_letter_idx <- panel_letter_idx + 1L
    }
  }

  mtext("Differential Circadian Power Analysis",
        outer = TRUE, cex = 1.1, font = 2)

  invisible(panel_data)
}
