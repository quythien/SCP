#!/usr/bin/env Rscript
# =============================================================================
# Figure 6, one-column (86 mm) version for the Bioinformatics port: the five
# panels of fig6_active_design_5panel_gtex.R in two rows (A B over C D E, the
# bottom row sharing one power axis) with the B legend below; set TWO_ROWS to
# FALSE for the three-row variant (A B / C D / E centred). Plot-only: same
# caches, same B -> colour map, same convention that the non-identifiable K=2
# designs B=3,4 are drawn flat at 0.
# Run from the package root: Rscript examples/publication/onecol/fig6_onecol.R
# =============================================================================
setwd(Sys.getenv("SCP_ROOT", unset = "."))
source("examples/publication/onecol/onecol_common.R")

AB_CACHE   <- "output/two_harmonic/results/fig6AB_liver_active.rds"
DIFF_CACHE <- "output/diagnostics/Fig6_differential_Bsweep_fdr20_fine.rds"
ab   <- readRDS(AB_CACHE)
diff <- readRDS(DIFF_CACHE)

B_ALL   <- c(3L, 4L, 6L, 8L, 12L, 24L)
COL_ALL <- c("#440154", "#3b528b", "#21918c", "#5ec962", "#fde725", "#7d028c")
names(COL_ALL) <- B_ALL

nN    <- length(ab$N)
pwr_B <- matrix(0, nN, length(B_ALL), dimnames = list(NULL, paste0("B=", B_ALL)))
se_B  <- matrix(0, nN, length(B_ALL), dimnames = list(NULL, paste0("B=", B_ALL)))
for (j in seq_along(ab$B_K2)) {
  bcol <- match(ab$B_K2[j], B_ALL)
  pwr_B[, bcol] <- ab$pwr_K2[, j]
  se_B[,  bcol] <- ab$se_K2[, j]
}

CX <- numeric(0)   # plot-box centres, for centring the legend bar
draw_panel <- function(N, M, SE, B_grid, letter, main_text, jitter_step, ylab = TRUE) {
  cols <- COL_ALL[as.character(B_grid)]
  plot(NA, xlim = c(0, max(N) * 1.04), ylim = c(0, 100),
       xlab = "Sample size (N)", ylab = if (ylab) "Power (%)" else "",
       yaxt = if (ylab) "s" else "n", main = "")
  if (!ylab) axis(2, labels = FALSE)
  panel_title(letter, main_text)
  abline(h = 80, lty = 2, col = "grey50", lwd = 0.9); grid(lwd = 0.5)
  for (k in seq_along(B_grid)) {
    xj <- N + (k - (length(B_grid) + 1) / 2) * jitter_step
    y  <- 100 * M[, k]; s <- 100 * SE[, k]
    lo <- pmax(0, y - s); hi <- pmin(100, y + s)
    ok <- is.finite(lo) & is.finite(hi) & s > 0
    if (any(ok))
      arrows(xj[ok], lo[ok], xj[ok], hi[ok], code = 3, angle = 90,
             length = 0.015, col = cols[k], lwd = 0.7)
    lines(xj, y, type = "o", pch = 19, lwd = 1.1, col = cols[k], cex = 0.35)
  }
  CX <<- c(CX, plot_centre())
}

OUT <- file.path(OUT_DIR, "Fig6_active_design_5panel_1col.pdf")
TWO_ROWS <- TRUE   # FALSE: A B / C D / E centred (three rows); TRUE: A B / C D E (Fig 2 style)
if (TWO_ROWS) {
  # Row 2 shares one y axis, so C gets the axis margin and D, E do not; the six
  # column widths below make the three plot areas equal and rows 1 and 2 align.
  x <- (COL_W - 0.36) / 3
  widths <- c((x + 0.36) / 2, (x + 0.36) / 2, (x - 0.36) / 2, (x + 0.36) / 2, x / 2, x / 2)
  pdf(OUT, width = COL_W, height = 3.35, pointsize = PT)
  layout(matrix(c(1, 1, 1, 2, 2, 2,
                  3, 3, 4, 4, 5, 5,
                  6, 6, 6, 6, 6, 6), nrow = 3, byrow = TRUE),
         widths = widths, heights = c(1.4, 1.5, 0.42))
  base_par(c(0.34, 0.44, 0.24, 0.05))
  draw_panel(ab$N, ab$pwr_K1, ab$se_K1, ab$B_K1, "A", "Single-harmonic (K=1)", 1.2)
  draw_panel(ab$N, pwr_B, se_B, B_ALL, "B", "Two-harmonic (K=2)", 1.2)
  base_par(c(0.36, 0.44, 0.24, 0.03))
  draw_panel(diff$N_GRID, diff$power$DR, diff$se$DR, diff$B_GRID, "C", "DR", 2.2)
  base_par(c(0.36, 0.08, 0.24, 0.03))
  draw_panel(diff$N_GRID, diff$power$DP, diff$se$DP, diff$B_GRID, "D", "DP", 2.2, ylab = FALSE)
  draw_panel(diff$N_GRID, diff$power$DM, diff$se$DM, diff$B_GRID, "E", "DM", 2.2, ylab = FALSE)
} else {
  pdf(OUT, width = COL_W, height = 4.4, pointsize = PT)
  layout(matrix(c(1, 1, 2, 2,
                  3, 3, 4, 4,
                  0, 5, 5, 0,
                  6, 6, 6, 6), nrow = 4, byrow = TRUE),
         heights = c(1.33, 1.33, 1.33, 0.41))
  base_par(c(0.34, 0.44, 0.24, 0.05))
  draw_panel(ab$N, ab$pwr_K1, ab$se_K1, ab$B_K1, "A", "Single-harmonic (K=1)", 1.2)
  draw_panel(ab$N, pwr_B, se_B, B_ALL, "B", "Two-harmonic (K=2)", 1.2)
  draw_panel(diff$N_GRID, diff$power$DR, diff$se$DR, diff$B_GRID, "C", "DR (rhythmicity)", 2.2)
  draw_panel(diff$N_GRID, diff$power$DP, diff$se$DP, diff$B_GRID, "D", "DP (phase)", 2.2)
  draw_panel(diff$N_GRID, diff$power$DM, diff$se$DM, diff$B_GRID, "E", "DM (mesor)", 2.2)
}
# one legend bar under all panels
legend_strip(paste0("B=", B_ALL), COL_ALL, pch = 19, ncol = length(B_ALL),
             title = "Time points per cycle (B)", lwd = 1.2, at = mean(CX))
dev.off()
cat(sprintf("Saved %s (K=1 N80~%.0f, K=2 N80~%.0f)\n", OUT, mean(ab$n80_K1), mean(ab$n80_K2)))
