#' Fig 6 generator (GTEx recalibration): B-invariance of power under a balanced
#' active design. Five panels:
#'   (A,B) single-cohort biomarker detection, K=1 and K=2, on the GTEx Liver pilot;
#'   (C,D,E) two-group differential power for DR, DP and DM on GTEx Adrenal vs Liver.
#'
#' This is a PLOT-ONLY render (Stage 3). Both calibrations come from pre-computed
#' caches; no simulations are run here.
#'   * A/B biomarker panels: output/two_harmonic/results/fig6AB_liver_active.rds
#'       (K=1 over B={3,4,6,8,12,24}; K=2 over B={6,8,12,24}). Power is a
#'       proportion in [0,1]; multiply by 100 for %.
#'   * C/D/E differential panels: output/diagnostics/Fig6_differential_Bsweep.rds
#'       (DR/DP/DM over N_GRID x B_GRID={4,6,8,12,24}). Power is a proportion.
#'
#' Convention preserved from the earlier version: for the two-harmonic panel (B),
#' the non-identifiable designs B=3 and B=4 (K=2 needs B>=5) are drawn as flat
#' lines at 0 power, and the B->colour map is shared across every panel so the
#' legend reads left to right.

setwd(Sys.getenv("SCP_ROOT", unset = "."))  # run from the package root

# ---------------------------------------------------------------------------
# Caches
# ---------------------------------------------------------------------------
AB_CACHE   <- "output/two_harmonic/results/fig6AB_liver_active.rds"   # A/B (GTEx Liver)
DIFF_CACHE <- "output/diagnostics/Fig6_differential_Bsweep_fdr20_fine.rds"  # C/D/E (Adrenal vs Liver, BH-FDR 20%, fine N grid to match Fig 2)
OUT_PDF    <- "submission/figures/Fig6_active_design_5panel.pdf"
OUT_PDF2   <- "figures/Fig6_active_design_5panel.pdf"

ab   <- readRDS(AB_CACHE)
diff <- readRDS(DIFF_CACHE)

# B -> colour map, consistent across every panel (viridis; last = magenta as shipped)
B_ALL   <- c(3L, 4L, 6L, 8L, 12L, 24L)
COL_ALL <- c("#440154", "#3b528b", "#21918c", "#5ec962", "#fde725", "#7d028c")
names(COL_ALL) <- B_ALL

# ---------------------------------------------------------------------------
# Panel B: assemble a full [N x 6] matrix over B_ALL. Real K=2 curves for
# B={6,8,12,24}; B={3,4} are non-identifiable and forced flat at 0 (no bars).
# ---------------------------------------------------------------------------
nN     <- length(ab$N)
pwr_B  <- matrix(0, nN, length(B_ALL), dimnames = list(NULL, paste0("B=", B_ALL)))
se_B   <- matrix(0, nN, length(B_ALL), dimnames = list(NULL, paste0("B=", B_ALL)))
for (j in seq_along(ab$B_K2)) {
  bcol <- match(ab$B_K2[j], B_ALL)
  pwr_B[, bcol] <- ab$pwr_K2[, j]
  se_B[,  bcol] <- ab$se_K2[, j]
}

# ---------------------------------------------------------------------------
# Drawing helpers
# ---------------------------------------------------------------------------
CEX_MAIN <- 1.30
LEG_CEX  <- 0.80

draw_panel <- function(N, M, SE, B_grid, letter, main_text,
                       jitter_step = 1.2, show_legend = FALSE) {
  cols <- COL_ALL[as.character(B_grid)]
  plot(NA, xlim = c(0, max(N) * 1.04), ylim = c(0, 100),
       xlab = "Sample size (N)", ylab = "Power (%)", main = "")
  title(main = sprintf("%s   %s", letter, main_text),
        adj = 0.5, font.main = 2, cex.main = CEX_MAIN, line = 0.35)
  abline(h = 80, lty = 2, col = "grey50", lwd = 1.3); grid()
  for (k in seq_along(B_grid)) {
    xj <- N + (k - (length(B_grid) + 1) / 2) * jitter_step
    y  <- 100 * M[, k]; s <- 100 * SE[, k]
    lo <- pmax(0, y - s); hi <- pmin(100, y + s)
    ok <- is.finite(lo) & is.finite(hi) & s > 0
    if (any(ok))
      arrows(xj[ok], lo[ok], xj[ok], hi[ok], code = 3, angle = 90,
             length = 0.025, col = cols[k], lwd = 1.0)
    lines(xj, y, type = "o", pch = 19, lwd = 1.7, col = cols[k], cex = 0.6)
  }
  if (show_legend)
    legend("bottomright", paste0("B=", B_grid), col = cols, lty = 1, pch = 19,
           lwd = 1.5, cex = LEG_CEX, bty = "o", box.col = "grey70", box.lwd = 0.5,
           inset = 0.01, y.intersp = 0.84, bg = "white", seg.len = 1.4)
}

# ---------------------------------------------------------------------------
# Render
# ---------------------------------------------------------------------------
dir.create(dirname(OUT_PDF),  recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(OUT_PDF2), recursive = TRUE, showWarnings = FALSE)

pdf(OUT_PDF, width = 10.5, height = 7.0)
layout(matrix(c(1, 1, 1, 2, 2, 2,
                3, 3, 4, 4, 5, 5), nrow = 2, byrow = TRUE))
par(mai = c(0.62, 0.62, 0.34, 0.10), mgp = c(2.1, 0.55, 0),
    oma = c(0, 0, 3.1, 0), cex.axis = 1.15, cex.lab = 1.35, font.main = 2)

# A / B: GTEx Liver biomarker detection (legend on A)
draw_panel(ab$N, ab$pwr_K1, ab$se_K1, ab$B_K1, "A", "Single-harmonic (K=1)",
           jitter_step = 1.2, show_legend = TRUE)
draw_panel(ab$N, pwr_B, se_B, B_ALL, "B", "Two-harmonic (K=2)",
           jitter_step = 1.2, show_legend = FALSE)

# C / D / E: GTEx Adrenal vs Liver differential (legend on C)
draw_panel(diff$N_GRID, diff$power$DR, diff$se$DR, diff$B_GRID, "C",
           "Differential rhythmicity (DR)", jitter_step = 2.2, show_legend = TRUE)
draw_panel(diff$N_GRID, diff$power$DP, diff$se$DP, diff$B_GRID, "D",
           "Differential phase (DP)", jitter_step = 2.2, show_legend = FALSE)
draw_panel(diff$N_GRID, diff$power$DM, diff$se$DM, diff$B_GRID, "E",
           "Differential mesor (DM)", jitter_step = 2.2, show_legend = FALSE)

mtext("Active-design B vs m trade-off in circadian study design",
      outer = TRUE, side = 3, line = 1.7, font = 2, cex = 1.25)
mtext("Biomarker detection (A, B, GTEx Liver); differential endpoints (C-E, Adrenal Gland vs Liver)",
      outer = TRUE, side = 3, line = 0.25, font = 1, cex = 0.95)
dev.off()

file.copy(OUT_PDF, OUT_PDF2, overwrite = TRUE)
cat(sprintf("Saved %s (K=1 N80~%.0f, K=2 N80~%.0f)\n",
            OUT_PDF, mean(ab$n80_K1), mean(ab$n80_K2)))
cat(sprintf("Copied to %s\n", OUT_PDF2))
