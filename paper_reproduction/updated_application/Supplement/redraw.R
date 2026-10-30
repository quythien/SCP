# Supplementary Figure S1: two-harmonic detection power.

# 1. Set the plotting style ------------------------------------------------
FIGURE_DIR <- "Supplement"
source("code/onecol_common.R")

# 2. Read the supplied panel values ---------------------------------------
plot_data <- readRDS(file.path(FIGURE_DIR, "data", "panel_data.rds"))
DISPLAY_N_A <- plot_data$DISPLAY_N_A
DISPLAY_N_B <- plot_data$DISPLAY_N_B
ALPHA_GRID <- plot_data$ALPHA_GRID
R_BREAKS <- plot_data$R_BREAKS
R_LABELS <- plot_data$R_LABELS
VLINE_FDR <- plot_data$VLINE_FDR
VLINE_POWER <- plot_data$VLINE_POWER
pwr_mean <- plot_data$pwr_mean
pwr_se <- plot_data$pwr_se
pwr_strat_mean <- plot_data$pwr_strat_mean
pwr_strat_se <- plot_data$pwr_strat_se
vline_n <- ceiling(plot_data$vline_n)
n_N_A <- length(DISPLAY_N_A)
n_N_B <- length(DISPLAY_N_B)
n_a <- length(ALPHA_GRID)

# 3. Draw the figure ------------------------------------------------------
OUT <- file.path(OUT_DIR, "Supplementary_harmonic_framework.pdf")
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
