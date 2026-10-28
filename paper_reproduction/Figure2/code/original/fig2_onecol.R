#!/usr/bin/env Rscript
# =============================================================================
# Figure 2, one-column (86 mm) version for the Bioinformatics port: the three
# differential endpoints (DR, DP, DM) in one row of square panels that share
# the power axis, with the FDR legend below. Plot-only: the marginal power
# curves come from plotDiffPower() itself (its invisible return value, drawn
# into a scratch PDF), so the numbers match the full-width figure exactly.
# Run from the package root: Rscript examples/publication/onecol/fig2_onecol.R
# =============================================================================
setwd(Sys.getenv("SCP_ROOT", unset = "."))
suppressMessages(library(SCP))
source("examples/publication/onecol/onecol_common.R")

DISPLAY_N <- c(20, 40, 60, 80, 100, 120, 150, 200)
FDR_GRID  <- c(0.01, 0.05, 0.10, 0.20)
EP        <- c("DR", "DP", "DM")
EP_NAME   <- c(DR = "Rhythmicity (DR)", DP = "Phase (DP)", DM = "Mesor (DM)")
VLINE_FDR <- 0.20

rds_diff <- tryCatch(
  pick_latest("output/differential/results", "^diff_power_ADR_vs_LIV_paired_.*\\.rds$"),
  error = function(e) pick_latest("output/differential/results",
                                  "^diff_power_ADR_vs_LIV_2026.*\\.rds$"))
cat(sprintf("  RDS: %s\n", basename(rds_diff)))
res <- readRDS(rds_diff)

tmp <- tempfile(fileext = ".pdf")
pd <- plotDiffPower(res_list = list(res), comp_labels = "GTEx Adrenal Gland vs Liver",
                    endpoints = EP, fdr_thresholds = FDR_GRID, panel_fdr = 0.05,
                    vline_power = 0.80, vline_fdr = VLINE_FDR,
                    display_sizes = DISPLAY_N, out_pdf = tmp, width = 12, height = 4.7)
unlink(tmp)
n80 <- vapply(EP, function(ep)
  npower(res, target_power = 0.80, fdr = VLINE_FDR, endpoint = ep)$n, numeric(1))
disp_idx <- which(res$sample_sizes %in% DISPLAY_N)
ss_disp  <- res$sample_sizes[disp_idx]
thresh_cols   <- c("darkgreen", "steelblue", "orange", "red")
thresh_labels <- paste0("FDR ", round(100 * FDR_GRID), "%")

OUT <- file.path(OUT_DIR, "Fig2_differential_ADR_vs_LIV_1col.pdf")
# row height chosen so each plot box is square (width is fixed by the column)
pdf(OUT, width = COL_W, height = 1.76, pointsize = PT)
layout(matrix(c(1, 2, 3, 4, 4, 4), nrow = 2, byrow = TRUE),
       widths = c(1.38, 1, 1), heights = c(1.48, 0.27))
cx <- numeric(0)
for (i in seq_along(EP)) {
  ep <- EP[i]
  first <- i == 1
  base_par(if (first) c(0.36, 0.46, 0.24, 0.04) else c(0.36, 0.08, 0.24, 0.04))
  mm <- pd[[paste0(ep, "_1")]]$marginal_mean[disp_idx, , drop = FALSE]
  se <- pd[[paste0(ep, "_1")]]$marginal_se[disp_idx, , drop = FALSE]
  matplot(ss_disp, 100 * mm, type = "b", pch = 19, lwd = LWD, cex = PT_CEX,
          col = thresh_cols, lty = 1, xlim = c(0, max(ss_disp) * 1.05),
          ylim = c(0, 100), xlab = "N per group", ylab = if (first) "Power (%)" else "",
          yaxt = if (first) "s" else "n", main = "")
  if (!first) axis(2, labels = FALSE)
  panel_title(LETTERS[i], ep)
  for (t in seq_along(thresh_cols))
    se_bars(ss_disp, 100 * mm[, t], 100 * se[, t], thresh_cols[t])
  abline(h = 80, lty = 2, col = "grey50", lwd = 0.9)
  grid(lwd = 0.5)
  if (is.finite(n80[ep])) {
    abline(v = n80[ep], lty = 2, col = adjustcolor("steelblue", 0.7), lwd = 1.0)
    text(n80[ep], 92, sprintf("N=%d", n80[ep]), col = "steelblue",
         cex = CEX_LEG, adj = c(1.08, 0.5), font = 2)
  }
  cx <- c(cx, plot_centre())
}
legend_strip(thresh_labels, thresh_cols, pch = 19, ncol = 4, at = mean(cx))
dev.off()
cat(sprintf("Saved %s (N80 at FDR %g%%: %s)\n", OUT, 100 * VLINE_FDR,
            paste(sprintf("%s=%d", EP, n80), collapse = ", ")))
