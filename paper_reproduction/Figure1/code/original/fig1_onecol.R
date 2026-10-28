#!/usr/bin/env Rscript
# =============================================================================
# Figure 1, one-column (86 mm) version for the Bioinformatics port.
# 3 rows x 2 columns: left = GTEx Adrenal Gland (A-C), right = GTEx Liver (D-F).
# Plot-only: reads the same cached runs as fig1_fig2_replot.R. The power and
# true-discovery computation below is copied verbatim from
# plotSingleCohortPower() (SCP 0.4.57) so the curves match the full-width
# figure; keep the two in sync if the package function changes.
# Run from the package root: Rscript examples/publication/onecol/fig1_onecol.R
# =============================================================================
setwd(Sys.getenv("SCP_ROOT", unset = "."))
suppressMessages(library(SCP))
source("examples/publication/onecol/onecol_common.R")

DISPLAY_N   <- c(20, 40, 60, 80, 100, 120, 150, 200)
FDR_GRID    <- c(0.01, 0.05, 0.10, 0.20)
PANEL_FDR   <- 0.05
VLINE_POWER <- 0.80
VLINE_FDR   <- 0.05
R_MAX       <- 3

prep_sc <- function(res, panel_a_res, fdr_thresholds = FDR_GRID,
                    panel_fdr = PANEL_FDR, vline_power = VLINE_POWER,
                    vline_fdr = VLINE_FDR, display_sizes = DISPLAY_N,
                    r_max = R_MAX, p.adjust.method = "BH") {
  sample_sizes  <- res$sample_sizes
  nsims         <- res$nsims
  r_strata      <- res$r_strata
  strata_labels <- res$strata_labels
  n_sizes       <- length(sample_sizes)
  n_strata      <- length(strata_labels)
  n_thresh      <- length(fdr_thresholds)
  base_cols     <- c("darkgreen", "steelblue", "orange", "red", "purple", "brown", "magenta")
  thresh_cols   <- base_cols[seq_len(n_thresh)]
  thresh_labels <- paste0("FDR ", round(100 * fdr_thresholds, 1), "%")
  size_colors   <- rainbow(n_sizes, s = 0.6, v = 0.8)
  disp_idx <- which(sample_sizes %in% display_sizes)

  power_arr <- array(NA_real_, dim = c(n_sizes, n_strata, n_thresh, nsims))
  TD_arr    <- array(NA_real_, dim = c(n_sizes, n_strata, n_thresh, nsims))
  for (j in seq_len(n_sizes)) {
    for (s in seq_len(nsims)) {
      pvals    <- res$pvalues[j, , s]
      r_values <- res$r_values_list[[j]][[s]]
      is_rhythmic <- r_values > 0
      r_for_strat <- ifelse(is_rhythmic, r_values, 0)
      xgr <- cut(r_for_strat, breaks = r_strata, include.lowest = TRUE, labels = FALSE)
      xgr[!is_rhythmic] <- NA
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

  # Panel A: marginal power from the companion (panel_a_res) run
  pa <- panel_a_res
  pa_nsizes  <- length(pa$sample_sizes)
  pa_marg_sim <- array(NA_real_, dim = c(pa_nsizes, n_thresh, pa$nsims))
  for (jj in seq_len(pa_nsizes)) {
    for (ss in seq_len(pa$nsims)) {
      rv2  <- pa$r_values_list[[jj]][[ss]]
      rhy2 <- rv2 > 0
      n_t2 <- sum(rhy2)
      pv2  <- pa$pvalues[jj, , ss]
      pv2[is.na(pv2)] <- 1
      qv2  <- p.adjust(pv2, method = p.adjust.method)
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
  disp_idx_a     <- which(sample_sizes_a %in% display_sizes)
  vline_n <- npower(pa, target_power = vline_power, fdr = vline_fdr)$n

  # Panels B & C at panel_fdr
  idx_fdr5  <- which.min(abs(fdr_thresholds - panel_fdr))
  fdr_label <- sprintf("FDR %g%%", panel_fdr * 100)
  mean_TD <- apply(TD_arr[, , idx_fdr5, , drop = FALSE], c(1, 2), mean, na.rm = TRUE)
  se_TD   <- apply(TD_arr[, , idx_fdr5, , drop = FALSE], c(1, 2),
                   function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
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
  mean_pow <- apply(power_arr[, , idx_fdr5, , drop = FALSE], c(1, 2), mean, na.rm = TRUE)
  se_pow   <- apply(power_arr[, , idx_fdr5, , drop = FALSE], c(1, 2),
                    function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

  # Collapse strata with r >= r_max into one ">r_max" bin
  left_bounds   <- r_strata[-length(r_strata)]
  collapse_from <- which(left_bounds >= r_max)
  if (length(collapse_from) >= 2) {
    keep <- seq_len(min(collapse_from) - 1)
    gene_counts_plt <- c(gene_counts[keep], sum(gene_counts[collapse_from]))
    mean_TD_plt <- cbind(mean_TD[, keep, drop = FALSE],
                         rowSums(mean_TD[, collapse_from, drop = FALSE]))
    se_TD_plt   <- cbind(se_TD[, keep, drop = FALSE],
                         sqrt(rowSums(se_TD[, collapse_from, drop = FALSE]^2)))
    gc_merged <- sum(gene_counts[collapse_from])
    pow_merged <- if (gc_merged > 0)
      rowSums(mean_TD[, collapse_from, drop = FALSE]) / gc_merged else rep(0, n_sizes)
    mean_pow_plt <- cbind(mean_pow[, keep, drop = FALSE], pow_merged)
    se_pow_plt   <- cbind(se_pow[, keep, drop = FALSE], matrix(NA_real_, n_sizes, 1))
    strata_labels_plt <- c(strata_labels[keep], sprintf(">%g", r_max))
  } else {
    gene_counts_plt <- gene_counts
    mean_TD_plt <- mean_TD;  se_TD_plt  <- se_TD
    mean_pow_plt <- mean_pow; se_pow_plt <- se_pow
    strata_labels_plt <- strata_labels
  }

  list(ss_a = sample_sizes_a[disp_idx_a],
       marg_mean = marginal_mean[disp_idx_a, , drop = FALSE],
       marg_se   = marginal_se[disp_idx_a, , drop = FALSE],
       vline_n = vline_n, thresh_cols = thresh_cols, thresh_labels = thresh_labels,
       fdr_label = fdr_label, strata = strata_labels_plt,
       ss_bc = sample_sizes[disp_idx], size_cols = size_colors[disp_idx],
       pow = mean_pow_plt[disp_idx, , drop = FALSE],
       pow_se = se_pow_plt[disp_idx, , drop = FALSE],
       td = mean_TD_plt[disp_idx, , drop = FALSE],
       td_se = se_TD_plt[disp_idx, , drop = FALSE],
       gene_counts = gene_counts_plt)
}

draw_A <- function(d, letter) {
  matplot(d$ss_a, 100 * d$marg_mean, type = "b", pch = 19, lwd = LWD, cex = PT_CEX,
          col = d$thresh_cols, lty = 1, xlim = c(0, max(d$ss_a) * 1.05),
          ylim = c(0, 100), xlab = "Sample size (N)", ylab = "", main = "",
          yaxt = "n")
  axis(2, at = seq(0, 100, 20))
  title(ylab = "Power (%)", line = YLAB_LINE)
  panel_title(letter, "Power vs sample size")
  for (t in seq_along(d$thresh_cols))
    se_bars(d$ss_a, 100 * d$marg_mean[, t], 100 * d$marg_se[, t], d$thresh_cols[t])
  abline(h = 80, lty = 2, col = "grey50", lwd = 0.9)
  grid(lwd = 0.5)
  if (!is.na(d$vline_n)) {
    abline(v = d$vline_n, lty = 2, col = adjustcolor("steelblue", 0.7), lwd = 1.0)
    right <- d$vline_n <= 0.6 * max(d$ss_a)      # label on the roomier side of the line
    text(d$vline_n, if (right) 35 else 10, sprintf("N=%d", d$vline_n), col = "steelblue",
         cex = CEX_LEG, adj = c(if (right) -0.1 else 1.1, 0.5), font = 2)
  }
}

strata_axis <- function(d) {
  # ticks at every bin, labels at every other bin (the last bin always
  # labelled). Labelled positions are passed alone: R drops labels it judges
  # overlapping, and with blank alternates it can drop the visible ones.
  k <- length(d$strata)
  axis(1, at = seq_len(k), labels = d$strata, las = 2, cex.axis = CEX_AXIS)
  title(xlab = expression(r == A/sigma), line = 4.0)
}

draw_B <- function(d, letter, ylab = TRUE) {
  k <- length(d$strata)
  matplot(seq_len(k), 100 * t(d$pow), type = "l", lwd = LWD, col = d$size_cols,
          lty = 1, xlim = c(0.5, k + 0.5), ylim = c(0, 100),
          xaxt = "n", yaxt = if (ylab) "s" else "n", xlab = "", ylab = "", main = "")
  if (!ylab) axis(2, labels = FALSE) else title(ylab = "Power (%)", line = YLAB_LINE)
  panel_title(letter, "Stratified power")
  strata_axis(d)
  grid(lwd = 0.5)
  for (j in seq_len(nrow(d$pow))) {
    points(seq_len(k), 100 * d$pow[j, ], pch = 19, col = d$size_cols[j], cex = PT_CEX)
    se_bars(seq_len(k), 100 * d$pow[j, ], 100 * d$pow_se[j, ], d$size_cols[j])
  }
}

draw_C <- function(d, letter) {
  k <- length(d$strata)
  plot(seq_len(k), rep(0, k), type = "n", xaxt = "n",
       xlim = c(0.5, k + 0.5), ylim = c(0, max(d$gene_counts) * 1.15),
       xlab = "", ylab = "", main = "")
  title(ylab = "# True discoveries", line = YLAB_LINE)
  panel_title(letter, "True discoveries")
  strata_axis(d)
  grid(lwd = 0.5)
  step_x <- rep(seq(0.5, k + 0.5, by = 1), each = 2)
  step_y <- c(0, rep(d$gene_counts, each = 2), 0)
  polygon(step_x, step_y, col = "#cccccc55", border = NA)
  lines(step_x, step_y, col = "grey60", lwd = 1.0, lty = 2)
  for (j in seq_len(nrow(d$td))) {
    lines(seq_len(k), d$td[j, ], col = d$size_cols[j], lwd = LWD)
    points(seq_len(k), d$td[j, ], pch = 19, col = d$size_cols[j], cex = PT_CEX)
    se_bars(seq_len(k), d$td[j, ], d$td_se[j, ], d$size_cols[j])
  }
}

load_pair <- function(primary_rx, panelA_rx) {
  rds_p <- pick_latest("output/single_cohort/results", primary_rx)
  rds_a <- pick_latest("output/single_cohort/results", panelA_rx)
  cat(sprintf("  B/C: %s\n  A  : %s\n", basename(rds_p), basename(rds_a)))
  prep_sc(readRDS(rds_p), readRDS(rds_a))
}
cat("GTEx Adrenal Gland\n")
adr <- load_pair("^single_cohort_power_GTEx_AdrenalGland_2026.*\\.rds$",
                 "^single_cohort_power_GTEx_AdrenalGland_paired_.*\\.rds$")
cat("GTEx Liver\n")
liv <- load_pair("^single_cohort_power_GTEx_Liver_2026.*\\.rds$",
                 "^single_cohort_power_GTEx_Liver_paired_.*\\.rds$")

OUT <- file.path(OUT_DIR, "Fig1_single_cohort_1col.pdf")
stopifnot(identical(adr$ss_bc, liv$ss_bc), identical(adr$size_cols, liv$size_cols))
# 3 x 2 grid: Adrenal Gland A-C in the left column, Liver D-F in the right
# column; FDR legend under row A/D, target-N legend under the last row.
pdf(OUT, width = COL_W, height = 5.0, pointsize = PT)
layout(matrix(c(1, 2, 3, 3, 4, 5, 6, 7, 8, 8), nrow = 5, byrow = TRUE),   # cells numbered in drawing order
       heights = c(1.24, 0.22, 1.6, 1.6, 0.34))
par(oma = c(0, 0, 1.3, 0))
YLAB_LINE <- 2.5                       # y titles aligned across rows (row C has 4-digit ticks)
LEFT <- 0.54                           # one left margin for every row: equal plot widths
base_par(c(0.34, LEFT, 0.20, 0.05))
draw_A(adr, "A"); cx_adr <- plot_centre()   # plot-box centre of A
draw_A(liv, "D"); cx_liv <- plot_centre()   # plot-box centre of D
cx_mid <- mean(c(cx_adr, cx_liv))           # legends centred on the pair of columns
legend_strip(adr$thresh_labels, adr$thresh_cols, pch = 19, ncol = 4, at = cx_mid)      # A, D
base_par(c(0.66, LEFT, 0.20, 0.05))
draw_B(adr, "B"); draw_B(liv, "E")
draw_C(adr, "C"); draw_C(liv, "F")
legend_strip(paste0("N=", adr$ss_bc), adr$size_cols, pch = NA, ncol = 4, at = cx_mid)  # B, C, E, F
# column headers centred on the plot areas of A and D (device NDC)
mtext("GTEx Adrenal Gland", side = 3, outer = TRUE, at = cx_adr, line = 0.1, font = 2, cex = 1)
mtext("GTEx Liver",         side = 3, outer = TRUE, at = cx_liv, line = 0.1, font = 2, cex = 1)
dev.off()
cat(sprintf("Saved %s (N80 at FDR %g%%: Adrenal %s, Liver %s)\n", OUT, 100 * VLINE_FDR,
            adr$vline_n, liv$vline_n))
