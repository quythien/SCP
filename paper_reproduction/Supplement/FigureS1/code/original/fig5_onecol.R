#!/usr/bin/env Rscript
# =============================================================================
# Figure 5, one-column (86 mm) version for the Bioinformatics port: panels A
# (K=2 power vs N) and B (power by first-harmonic effect size) side by side.
# Plot-only: the summaries are computed from the cached simulation results
# exactly as in two_harmonic/fig5_twoharm_framework.R (RUN_SIMS = FALSE path);
# keep the two in sync if that script changes.
# Run from the package root: Rscript examples/publication/onecol/fig5_onecol.R
# =============================================================================
setwd(Sys.getenv("SCP_ROOT", unset = "."))
source("examples/publication/onecol/onecol_common.R")

RES_PATH    <- "output/two_harmonic/results/fig5_panelAB_GTExLiver_topK300_unpaired.rds"
RES_PATH_A  <- "output/two_harmonic/results/fig5_panelA_GTExLiver_topK300_paired.rds"
ALPHA_GRID  <- c(0.01, 0.05, 0.10, 0.20)
R_BREAKS    <- c(seq(0, 2, by = 0.25), Inf)
R_LABELS    <- c(sprintf("(%g,%g]", head(R_BREAKS, -2), R_BREAKS[seq(2, length(R_BREAKS) - 1)]), ">2")
VLINE_FDR   <- 0.05
VLINE_POWER <- 0.80

res_by_N   <- readRDS(RES_PATH)
res_panelA <- if (file.exists(RES_PATH_A)) readRDS(RES_PATH_A) else res_by_N
DISPLAY_N_A <- if (!identical(res_panelA, res_by_N)) c(20, 40, 60, 80, 100, 120, 150, 200, 250, 300) else c(20, 40, 60, 80, 100, 120, 150, 200)
DISPLAY_N_B <- c(20, 40, 60, 80, 100, 120, 150, 200)
n_N_A <- length(DISPLAY_N_A); n_N_B <- length(DISPLAY_N_B); n_a <- length(ALPHA_GRID)

pwr_mean <- matrix(NA_real_, n_N_A, n_a); pwr_se <- matrix(NA_real_, n_N_A, n_a)
for (j in seq_along(DISPLAY_N_A)) {
  sims <- res_panelA[[j]]; sims <- sims[vapply(sims, is.list, logical(1))]
  if (length(sims) == 0) next
  per_rep <- vapply(sims, function(s) {
    n_tgt <- sum(s$is_rhythmic)
    sapply(ALPHA_GRID, function(a) { q <- p.adjust(s$pvalue, "BH"); sum(q <= a & s$is_rhythmic) / max(1, n_tgt) })
  }, FUN.VALUE = numeric(n_a))
  pwr_mean[j, ] <- rowMeans(per_rep, na.rm = TRUE)
  pwr_se[j, ]   <- apply(per_rep, 1, sd, na.rm = TRUE) / sqrt(length(sims))
}
pwr_strat_mean <- matrix(NA_real_, n_N_B, length(R_LABELS)); pwr_strat_se <- pwr_strat_mean
for (j in seq_along(DISPLAY_N_B)) {
  sims <- res_by_N[[j]]; sims <- sims[vapply(sims, is.list, logical(1))]
  if (length(sims) == 0) next
  per_rep <- vapply(sims, function(s) {
    q <- p.adjust(s$pvalue, "BH"); disc <- q <= 0.05 & s$is_rhythmic
    bins <- cut(s$r_tilde, breaks = R_BREAKS, include.lowest = TRUE, labels = FALSE)
    sapply(seq_along(R_LABELS), function(k) {
      in_k <- s$is_rhythmic & !is.na(bins) & bins == k; n_tgt <- sum(in_k)
      if (n_tgt == 0) NA_real_ else sum(disc & in_k) / n_tgt })
  }, FUN.VALUE = numeric(length(R_LABELS)))
  pwr_strat_mean[j, ] <- rowMeans(per_rep, na.rm = TRUE)
  pwr_strat_se[j, ]   <- apply(per_rep, 1, sd, na.rm = TRUE) / sqrt(length(sims))
}
interp_n80 <- function(N_vals, p_vals, target = VLINE_POWER) {
  ok <- !is.na(p_vals) & is.finite(N_vals); N_vals <- N_vals[ok]; p_vals <- p_vals[ok]
  if (length(N_vals) < 2) return(NA_real_)
  if (max(p_vals) < target) return(NA_real_)
  if (min(p_vals) >= target) return(N_vals[1])
  i <- min(which(p_vals >= target)); if (i == 1) return(N_vals[1])
  N_vals[i - 1] + (target - p_vals[i - 1]) * (N_vals[i] - N_vals[i - 1]) / (p_vals[i] - p_vals[i - 1])
}
vline_n <- interp_n80(DISPLAY_N_A, pwr_mean[, which.min(abs(ALPHA_GRID - VLINE_FDR))])

OUT <- file.path(OUT_DIR, "Fig5_twoharm_framework_1col.pdf")
pdf(OUT, width = COL_W, height = 2.6, pointsize = PT)
layout(matrix(1:4, nrow = 2, byrow = TRUE), heights = c(2.1, 0.48))
# ---- Panel A ----
base_par(c(0.66, 0.44, 0.24, 0.05))
alpha_cols <- c("darkgreen", "steelblue", "orange", "red")
matplot(DISPLAY_N_A, 100 * pwr_mean, type = "o", pch = 19, lwd = LWD, col = alpha_cols,
        lty = 1, cex = PT_CEX, xlim = c(0, max(DISPLAY_N_A) * 1.05), ylim = c(0, 100),
        xlab = "Sample size (N)", ylab = "Power (%)", main = "")
panel_title("A", "Power vs sample size"); cx_A <- plot_centre()
for (a in seq_along(ALPHA_GRID)) se_bars(DISPLAY_N_A, 100 * pwr_mean[, a], 100 * pwr_se[, a], alpha_cols[a])
abline(h = 80, lty = 2, col = "grey50", lwd = 0.9)
grid(lwd = 0.5)
if (is.finite(vline_n)) {
  abline(v = vline_n, lty = 2, col = adjustcolor("steelblue", 0.7), lwd = 1.0)
  text(vline_n, 35, sprintf("N=%d", round(vline_n)), col = "steelblue", cex = CEX_LEG, adj = c(1.1, 0.5), font = 2)
}
# ---- Panel B ----
size_colors <- rainbow(n_N_B, s = 0.6, v = 0.8)
matplot(seq_along(R_LABELS), 100 * t(pwr_strat_mean), type = "l", lwd = LWD, col = size_colors,
        lty = 1, xlim = c(0.5, length(R_LABELS) + 0.5), ylim = c(0, 100),
        xlab = "", ylab = "Power (%)", main = "", xaxt = "n")
panel_title("B", "Stratified power"); cx_B <- plot_centre()
axis(1, at = seq_along(R_LABELS), labels = R_LABELS, las = 2, cex.axis = CEX_AXIS)
title(xlab = expression(r[1] == A[1]/sigma), line = 4.0)
grid(lwd = 0.5)
for (j in seq_len(n_N_B)) {
  points(seq_along(R_LABELS), 100 * pwr_strat_mean[j, ], pch = 19, col = size_colors[j], cex = PT_CEX)
  se_bars(seq_along(R_LABELS), 100 * pwr_strat_mean[j, ], 100 * pwr_strat_se[j, ], size_colors[j])
}
legend_strip(paste0(round(100 * ALPHA_GRID), "%"), alpha_cols, pch = 19, ncol = 4,
             title = "FDR threshold", at = cx_A)                               # under A
legend_strip(as.character(DISPLAY_N_B), size_colors, pch = NA, ncol = 4,
             title = "Target sample size N", at = cx_B)                        # under B
dev.off()
cat(sprintf("Saved %s (N80 at FDR %g%%: %s)\n", OUT, 100 * VLINE_FDR,
            ifelse(is.na(vline_n), "out of range", sprintf("%d", round(vline_n)))))
