#' ============================================================
#' fig1_hybrid_replot.R - Fig 1A/1B dual-RDS replot (standalone variant)
#'
#' Renders single-cohort Fig 1A (Adrenal) and Fig 1B (Liver) by combining
#' the marginal-power Panel A from an alternate result file with the
#' effect-size-stratified Panels B-C from the primary cached result file.
#' This script duplicates the dual-RDS rendering that `fig1_fig2_replot.R`
#' now performs via `plotSingleCohortPower(..., panel_a_res = ...)`.
#'
#' Output:
#'   submission/figures/Fig1A_single_cohort_AdrenalGland.pdf
#'   submission/figures/Fig1B_single_cohort_Liver.pdf
#'   (and mirrored under output/main_figures/)
#' ============================================================

setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
old_wd <- setwd("code"); source("setup.R"); setwd(old_wd)

pick_latest <- function(dir, pattern) {
  cands <- list.files(dir, pattern, full.names = TRUE)
  if (length(cands) == 0) stop("No match for ", pattern, " in ", dir)
  cands[order(file.mtime(cands), decreasing = TRUE)[1]]
}

mirror_pdf <- function(src, dest_paths) {
  for (d in dest_paths) {
    dir.create(dirname(d), recursive = TRUE, showWarnings = FALSE)
    file.copy(src, d, overwrite = TRUE)
    cat(sprintf("Saved: %s\n", d))
  }
}

DISPLAY_N   <- c(20, 40, 60, 80, 100, 120, 150, 200)
R_MAX       <- 3
FDR_VEC     <- c(0.01, 0.05, 0.10, 0.20)
PANEL_FDR   <- 0.05
VLINE_POWER <- 0.80
VLINE_FDR   <- 0.05
THRESH_COLS <- c("darkgreen", "steelblue", "orange", "red")
THRESH_LAB  <- paste0("FDR ", round(100 * FDR_VEC), "%")

# Helper: compute marginal mean and se from a single-cohort RDS
compute_marginal <- function(res) {
  n_sizes  <- length(res$sample_sizes)
  nsims    <- res$nsims
  n_thresh <- length(FDR_VEC)
  marg_sim <- array(NA_real_, dim = c(n_sizes, n_thresh, nsims))
  for (j in seq_len(n_sizes)) {
    for (s in seq_len(nsims)) {
      rv  <- res$r_values_list[[j]][[s]]
      rhy <- rv > 0
      n_t <- sum(rhy)
      fdr_g <- res$fdr_list[[j]][[s]]
      for (t in seq_len(n_thresh)) {
        disc <- !is.na(fdr_g) & fdr_g <= FDR_VEC[t]
        marg_sim[j, t, s] <- if (n_t > 0) sum(disc & rhy) / n_t else NA_real_
      }
    }
  }
  list(
    mean = apply(marg_sim, c(1, 2), mean, na.rm = TRUE),
    se   = apply(marg_sim, c(1, 2),
                 function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))),
    n80  = npower(res, target_power = VLINE_POWER, fdr = VLINE_FDR)$n,
    sample_sizes = res$sample_sizes
  )
}

# Helper: compute strata mean_TD, mean_pow, gene_counts from a single-cohort RDS
compute_strata <- function(res) {
  n_sizes  <- length(res$sample_sizes)
  nsims    <- res$nsims
  r_strata <- res$r_strata
  n_strata <- length(r_strata) - 1L
  strata_labels <- c(
    sprintf("(%g,%g]", head(r_strata, -2), r_strata[seq(2, n_strata)]),
    sprintf(">%g", r_strata[n_strata])
  )

  TD_arr    <- array(NA_real_, dim = c(n_sizes, n_strata, length(FDR_VEC), nsims))
  power_arr <- array(NA_real_, dim = c(n_sizes, n_strata, length(FDR_VEC), nsims))
  count_mat <- matrix(0L, nrow = nsims * n_sizes, ncol = n_strata)

  row_i <- 1L
  for (j in seq_len(n_sizes)) {
    for (s in seq_len(nsims)) {
      rv  <- res$r_values_list[[j]][[s]]
      rhy <- rv > 0
      xg  <- cut(rv, breaks = r_strata, include.lowest = TRUE, labels = FALSE)
      fdr_g <- res$fdr_list[[j]][[s]]
      for (t in seq_along(FDR_VEC)) {
        disc <- !is.na(fdr_g) & fdr_g <= FDR_VEC[t]
        for (k in seq_len(n_strata)) {
          in_k <- !is.na(xg) & xg == k
          td   <- sum(disc & rhy & in_k, na.rm = TRUE)
          n_t  <- sum(rhy & in_k, na.rm = TRUE)
          TD_arr[j, k, t, s] <- td
          power_arr[j, k, t, s] <- if (n_t > 0) td / n_t else NA_real_
        }
      }
      count_mat[row_i, ] <- tabulate(xg[rhy], nbins = n_strata)
      row_i <- row_i + 1L
    }
  }
  count_mat <- count_mat[seq_len(row_i - 1L), , drop = FALSE]

  idx_fdr <- which.min(abs(FDR_VEC - PANEL_FDR))
  mean_TD <- apply(TD_arr[, , idx_fdr, , drop = FALSE], c(1, 2), mean, na.rm = TRUE)
  se_TD   <- apply(TD_arr[, , idx_fdr, , drop = FALSE], c(1, 2),
                   function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
  mean_pow <- apply(power_arr[, , idx_fdr, , drop = FALSE], c(1, 2),
                    mean, na.rm = TRUE)
  se_pow   <- apply(power_arr[, , idx_fdr, , drop = FALSE], c(1, 2),
                    function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
  gene_counts <- colMeans(count_mat)

  # Collapse strata with left bound >= R_MAX
  left_bounds <- r_strata[-length(r_strata)]
  cf <- which(left_bounds >= R_MAX)
  if (length(cf) >= 2) {
    keep <- seq_len(min(cf) - 1)
    gc_m <- sum(gene_counts[cf])
    gene_counts <- c(gene_counts[keep], gc_m)
    mean_TD <- cbind(mean_TD[, keep, drop = FALSE],
                     rowSums(mean_TD[, cf, drop = FALSE]))
    se_TD   <- cbind(se_TD[, keep, drop = FALSE],
                     sqrt(rowSums(se_TD[, cf, drop = FALSE]^2)))
    pow_m <- if (gc_m > 0) rowSums(mean_TD[, ncol(mean_TD), drop = FALSE]) / gc_m
             else rep(0, n_sizes)
    mean_pow <- cbind(mean_pow[, keep, drop = FALSE], pow_m)
    se_pow   <- cbind(se_pow[, keep, drop = FALSE], matrix(NA_real_, n_sizes, 1))
    strata_labels <- c(strata_labels[keep], sprintf(">%g", R_MAX))
  }

  list(
    mean_TD = mean_TD, se_TD = se_TD,
    mean_pow = mean_pow, se_pow = se_pow,
    gene_counts = gene_counts,
    strata_labels = strata_labels,
    sample_sizes = res$sample_sizes
  )
}

add_se_bars <- function(x, y, se, col = "black") {
  arrows(x, y - 1.96 * se, x, y + 1.96 * se,
         code = 3, angle = 90, length = 0.03, col = col, lwd = 1.2)
}

render_hybrid <- function(panel_a, strata, title, out_pdf,
                          width = 11.5, height = 4.2) {
  ss      <- panel_a$sample_sizes
  ss_s    <- strata$sample_sizes
  disp_a  <- which(ss   %in% DISPLAY_N)
  disp_s  <- which(ss_s %in% DISPLAY_N)
  ss_a    <- ss[disp_a]
  ss_b    <- ss_s[disp_s]
  size_cols <- rainbow(length(ss_s), s = 0.6, v = 0.8)

  n_strata <- length(strata$strata_labels)
  vline_n  <- panel_a$n80
  fdr_lab  <- sprintf("FDR %g%%", PANEL_FDR * 100)

  pdf(out_pdf, width = width, height = height)
  par(mfrow = c(1, 3), mai = c(1.15, 1.05, 0.70, 0.12),
      mgp = c(3.4, 0.7, 0), oma = c(0, 0, 2.6, 0),
      cex.axis = 1.35, cex.lab = 1.55, cex.main = 1.45, font.main = 2)

  # ---- Panel A: marginal power vs sample size ----
  matplot(ss_a, 100 * panel_a$mean[disp_a, , drop = FALSE],
          type = "b", pch = 19, lwd = 2,
          col = THRESH_COLS, lty = 1,
          xlim = c(0, max(ss_a) * 1.05), ylim = c(0, 100),
          xlab = "Sample size (n)", ylab = "Power (%)", main = "")
  title(main = "A   Power vs Sample Size",
        adj = 0, font.main = 2, cex.main = 1.45, line = 0.6)
  for (t in seq_along(FDR_VEC)) {
    add_se_bars(ss_a, 100 * panel_a$mean[disp_a, t],
                100 * panel_a$se[disp_a, t], col = THRESH_COLS[t])
  }
  abline(h = 80, lty = 2, col = "grey50", lwd = 1.2)
  if (!is.na(vline_n)) {
    abline(v = vline_n, lty = 2, col = adjustcolor("steelblue", 0.7), lwd = 1.5)
    text(vline_n, 15, sprintf("n=%d", vline_n),
         col = "steelblue", cex = 0.78, adj = -0.1)
  }
  grid()
  legend("bottomright", THRESH_LAB,
         col = THRESH_COLS, lty = 1, pch = 19, lwd = 2,
         cex = 0.62, bty = "o", box.col = "grey70", box.lwd = 0.5,
         inset = 0.01, y.intersp = 0.85)

  # ---- Panel B: effect-size-stratified power ----
  par(mgp = c(4.6, 0.6, 0))
  matplot(seq_len(n_strata),
          100 * t(strata$mean_pow[disp_s, , drop = FALSE]),
          type = "l", lwd = 2, col = size_cols[disp_s], lty = 1,
          xlim = c(0.5, n_strata + 0.5), ylim = c(0, 100), bty = "l",
          xlab = expression(tilde(r) == A/sigma), ylab = "Power (%)",
          main = "", xaxt = "n")
  title(main = bquote("B   " * bold("Stratified Power by") ~
                       bold(tilde(r)) ~ bold(.(sprintf("(%s)", fdr_lab)))),
        adj = 0, font.main = 2, cex.main = 1.45, line = 0.6)
  axis(1, at = seq_len(n_strata), labels = strata$strata_labels,
       las = 2, cex.axis = 0.95)
  for (j in disp_s) {
    points(seq_len(n_strata), 100 * strata$mean_pow[j, ],
           pch = 19, col = size_cols[j], cex = 0.65)
    add_se_bars(seq_len(n_strata), 100 * strata$mean_pow[j, ],
                100 * strata$se_pow[j, ], col = size_cols[j])
  }
  grid()
  legend("bottomright", paste0("n=", ss_b),
         col = size_cols[disp_s], lty = 1, lwd = 2,
         cex = 0.58, bty = "o", box.col = "grey70", box.lwd = 0.5,
         inset = 0.01, y.intersp = 0.85)

  # ---- Panel C: effect-size-stratified true discoveries ----
  y_max <- max(strata$gene_counts) * 1.15
  plot(seq_len(n_strata), rep(0, n_strata), type = "n", bty = "l",
       xlim = c(0.5, n_strata + 0.5), ylim = c(0, y_max),
       xlab = expression(tilde(r) == A/sigma),
       ylab = "# True Discoveries",
       main = "", xaxt = "n")
  title(main = bquote("C   " * bold("True Discoveries by") ~
                       bold(tilde(r)) ~ bold(.(sprintf("(%s)", fdr_lab)))),
        adj = 0, font.main = 2, cex.main = 1.45, line = 0.6)
  axis(1, at = seq_len(n_strata), labels = strata$strata_labels,
       las = 2, cex.axis = 0.95)
  step_x <- rep(seq(0.5, n_strata + 0.5, by = 1), each = 2)
  step_y <- c(0, rep(strata$gene_counts, each = 2), 0)
  polygon(step_x, step_y, col = "#cccccc55", border = NA)
  lines(step_x, step_y, col = "grey60", lwd = 1.5, lty = 2)
  for (j in disp_s) {
    lines(seq_len(n_strata), strata$mean_TD[j, ], col = size_cols[j], lwd = 2)
    points(seq_len(n_strata), strata$mean_TD[j, ], pch = 19,
           col = size_cols[j], cex = 0.65)
    add_se_bars(seq_len(n_strata), strata$mean_TD[j, ],
                strata$se_TD[j, ], col = size_cols[j])
  }
  grid()
  legend("topright",
         c(paste0("n=", ss_b), "# Target Discoveries"),
         col = c(size_cols[disp_s], "grey60"),
         lty = c(rep(1, length(disp_s)), 2),
         lwd = c(rep(2, length(disp_s)), 1.5),
         cex = 0.58, bty = "o", box.col = "grey70", box.lwd = 0.5,
         inset = 0.01, y.intersp = 0.85)
  par(mgp = c(3.0, 0.6, 0))

  if (nchar(title) > 0)
    mtext(title, outer = TRUE, cex = 1.5, font = 2)
  dev.off()
}

# ----- Fig 1A: Adrenal Gland -----
cat("=== Fig 1A (Adrenal Gland) ===\n")
unp_adr <- pick_latest("output/single_cohort/results",
                       "^single_cohort_power_GTEx_AdrenalGland_2026.*\\.rds$")
pai_adr <- pick_latest("output/single_cohort/results",
                       "^single_cohort_power_GTEx_AdrenalGland_paired_.*\\.rds$")
cat(sprintf("  primary (Panel B/C): %s\n", basename(unp_adr)))
cat(sprintf("  alt     (Panel A)  : %s\n", basename(pai_adr)))
adr_a  <- compute_marginal(readRDS(pai_adr))
adr_bc <- compute_strata(readRDS(unp_adr))
cat(sprintf("  N80 (FDR %g%%): %s\n", VLINE_FDR * 100,
            ifelse(is.na(adr_a$n80), "NA", adr_a$n80)))

tmp_adr <- tempfile(fileext = ".pdf")
render_hybrid(adr_a, adr_bc,
              title   = "Single-Cohort Power Analysis (GTEx Adrenal Gland)",
              out_pdf = tmp_adr)
mirror_pdf(tmp_adr,
           c("output/main_figures/Fig1A_single_cohort_AdrenalGland.pdf",
             "submission/figures/Fig1A_single_cohort_AdrenalGland.pdf"))

# ----- Fig 1B: Liver -----
cat("\n=== Fig 1B (Liver) ===\n")
unp_liv <- pick_latest("output/single_cohort/results",
                       "^single_cohort_power_GTEx_Liver_2026.*\\.rds$")
pai_liv <- pick_latest("output/single_cohort/results",
                       "^single_cohort_power_GTEx_Liver_paired_.*\\.rds$")
cat(sprintf("  primary (Panel B/C): %s\n", basename(unp_liv)))
cat(sprintf("  alt     (Panel A)  : %s\n", basename(pai_liv)))
liv_a  <- compute_marginal(readRDS(pai_liv))
liv_bc <- compute_strata(readRDS(unp_liv))
cat(sprintf("  N80 (FDR %g%%): %s\n", VLINE_FDR * 100,
            ifelse(is.na(liv_a$n80), "NA", liv_a$n80)))

tmp_liv <- tempfile(fileext = ".pdf")
render_hybrid(liv_a, liv_bc,
              title   = "Single-Cohort Power Analysis (GTEx Liver)",
              out_pdf = tmp_liv)
mirror_pdf(tmp_liv,
           c("output/main_figures/Fig1B_single_cohort_Liver.pdf",
             "submission/figures/Fig1B_single_cohort_Liver.pdf"))

cat("\n=== Done ===\n")
