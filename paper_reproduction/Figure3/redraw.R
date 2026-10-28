# Redraw Figure 3 from the supplied plotting values.

# 1. Load the shared figure style and the input table ----------------------
FIGURE_DIR <- "Figure3"
source("code/onecol_common.R")
v <- read.csv(file.path(FIGURE_DIR, "data", "fig3_values.csv"))
stopifnot(all(c("panel", "N", "plugin", "mean", "lo", "hi") %in% names(v)),
          setequal(v$panel, c("A", "B")), nrow(v) == 16L,
          all(v$plugin >= 0 & v$plugin <= 1), all(v$mean >= 0 & v$mean <= 1))

# 2. Use the original panel drawing function ------------------------------
CP <- "#0072B2"; CB <- "#D55E00"        # plug-in blue, bootstrap orange

draw_boot <- function(d, letter, title, ylab = TRUE) {
  d <- d[order(d$N), ]
  plot(NA, xlim = c(0, max(d$N) * 1.04), ylim = c(0, 100), xlab = "Sample size (N)",
       ylab = if (ylab) "Power (%)" else "", yaxt = if (ylab) "s" else "n", main = "")
  if (!ylab) axis(2, labels = FALSE)
  panel_title(letter, title)
  abline(h = 80, lty = 2, col = "grey50", lwd = 0.9); grid(lwd = 0.5)
  lines(d$N, 100 * d$plugin, lwd = 1.6, col = CP)
  lines(d$N, 100 * d$mean, lwd = 1.0, lty = 2, col = CB)
  ok <- is.finite(d$lo) & is.finite(d$hi) & (d$hi - d$lo) > 1e-6   # zero-length bars are not drawn
  arrows(d$N[ok], 100 * d$lo[ok], d$N[ok], 100 * d$hi[ok], code = 3, angle = 90,
         length = 0.025, lwd = 1.0, col = CB)
  points(d$N, 100 * d$mean, pch = 19, col = CB, cex = PT_CEX + 0.05)
}

OUT <- file.path(OUT_DIR, "Figure3.pdf")
# 3. Draw both panels and their shared legend -----------------------------
# Equal plot widths account for the y-axis margin in panel A.
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
