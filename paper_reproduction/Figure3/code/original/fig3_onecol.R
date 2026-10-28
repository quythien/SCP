#!/usr/bin/env Rscript
# =============================================================================
# Figure 3, one-column (86 mm) version for the Bioinformatics port, drawn in
# the shared style of the other redrawn figures (7-8 pt text, boxed panels,
# legend strip below). The plotted values come from figures/onecol/
# fig3_values.csv, recovered from the shipped vector PDF by
# fig3_extract_values.py (the generator fig3_bootstrap_subject.R needs the
# GSE160521 and GTEx CPM matrices on the server and saves no numbers).
# Recovered values reproduce the ones quoted in Section 3.4 (caudate N=40:
# plug-in 21%, bootstrap mean 72%, CI [34%, 97%]). Rerun the generator on the
# server with these settings to replace the CSV with computed values.
# Run from the package root: Rscript examples/publication/onecol/fig3_onecol.R
# =============================================================================
setwd(Sys.getenv("SCP_ROOT", unset = "."))
source("examples/publication/onecol/onecol_common.R")

v <- read.csv("figures/onecol/fig3_values.csv")
CP <- "#0072B2"; CB <- "#D55E00"        # plug-in blue, bootstrap orange (as shipped)

draw_boot <- function(d, letter, title, ylab = TRUE) {
  d <- d[order(d$N), ]
  plot(NA, xlim = c(0, max(d$N) * 1.04), ylim = c(0, 100), xlab = "Sample size (N)",
       ylab = if (ylab) "Power (%)" else "", yaxt = if (ylab) "s" else "n", main = "")
  if (!ylab) axis(2, labels = FALSE)
  panel_title(letter, title)
  abline(h = 80, lty = 2, col = "grey50", lwd = 0.9); grid(lwd = 0.5)
  lines(d$N, 100 * d$plugin, lwd = 1.6, col = CP)
  lines(d$N, 100 * d$mean, lwd = 1.0, lty = 2, col = CB)
  ok <- is.finite(d$lo) & is.finite(d$hi) & (d$hi - d$lo) > 1e-6   # zero-length bars are not drawn (as shipped)
  arrows(d$N[ok], 100 * d$lo[ok], d$N[ok], 100 * d$hi[ok], code = 3, angle = 90,
         length = 0.025, lwd = 1.0, col = CB)
  points(d$N, 100 * d$mean, pch = 19, col = CB, cex = PT_CEX + 0.05)
}

OUT <- file.path(OUT_DIR, "Fig3_bootstrap_redraw_1col.pdf")
# cell widths chosen so both plot boxes are the same width: A carries the
# y-axis margin (0.46 in), B only a hairline (0.08 in)
mai_A <- c(0.36, 0.46, 0.24, 0.04); mai_B <- c(0.36, 0.08, 0.24, 0.04)
cell_B <- (COL_W - ((mai_A[2] + mai_A[4]) - (mai_B[2] + mai_B[4]))) / 2
pdf(OUT, width = COL_W, height = 2.05, pointsize = PT)
layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE), widths = c(COL_W - cell_B, cell_B), heights = c(1.8, 0.25))
base_par(mai_A); draw_boot(v[v$panel == "A", ], "A", "Caudate (n = 59)"); cx_A <- plot_centre()
base_par(mai_B); draw_boot(v[v$panel == "B", ], "B", "GTEx muscle (n = 748)", ylab = FALSE); cx_B <- plot_centre()
par(mai = c(0.02, 0.03, 0.02, 0.03), xpd = NA); plot.new()
legend(grconvertX(mean(c(cx_A, cx_B)), "ndc", "user"), 0.5, xjust = 0.5, yjust = 0.5,
       c("Projected power (plug-in)", "Bootstrap mean and 95% interval"),
       col = c(CP, CB), lty = c(1, 2), lwd = c(1.6, 1.0), pch = c(NA, 19), pt.cex = PT_CEX + 0.05,
       ncol = 2, bty = "n", cex = 0.8, seg.len = 1.6, x.intersp = 0.5)
par(xpd = FALSE)
dev.off()
cat(sprintf("Saved %s\n", OUT))
