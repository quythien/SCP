# =====================================================================
# Replot Fig 3 (bootstrap single-cohort) from cached panels
# (output/bootstrap_sc/results/all_panels.rds). The original generator
# fig5_bootstrap_sc.R needs controlled-access pilot data, so the published
# 2-panel figure is regenerated here from the cache.
#
# Advisor preferences (May 2026 review) encoded here:
#  - 2 panels: A = Putamen SCZ, B = GTEx Thyroid (Thyroid capped at N=320).
#  - x-axis labels EVERY simulated N (gap.axis = -1 disables R's auto-thinning)
#    so e.g. 120 is never skipped; smaller tick font keeps them from colliding.
#  - Panel letter prefixed into the centred title ("A   <label>") so letter and
#    title align on one line, like Fig 1/2/4.
#  - Suptitle centred (adj = 0.5), trimmed; compact bottom-right legend.
#  - Wong palette: point estimate #0072B2 (blue), bootstrap #D55E00 (orange).
# Run: Rscript examples/publication/fig3_bootstrap_replot.R
# =====================================================================
setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")

panels_all <- readRDS("output/bootstrap_sc/results/all_panels.rds")
# The published figure shows Putamen SCZ (A) and GTEx Thyroid (B).
sel <- c(1L, 4L)
panels <- panels_all[sel]

col_plugin <- "#0072B2"   # Wong blue   - point estimate
col_boot   <- "#D55E00"   # Wong orange - bootstrap mean + 95% CI

dests <- c("figures/Fig3_bootstrap_singlecohort.pdf",
           "submission/figures/Fig3_bootstrap_singlecohort.pdf",
           "output/main_figures/Fig3_bootstrap_singlecohort.pdf")

tmp <- tempfile(fileext = ".pdf")
cairo_pdf(tmp, width = 7.6, height = 3.6)
par(mfrow = c(1, 2), mar = c(4.0, 4.2, 2.4, 1.6),
    mgp = c(2.5, 0.7, 0), oma = c(0, 0, 2.0, 0),
    cex.axis = 1.0, cex.lab = 1.15, font.main = 2)

letters_seq <- c("A", "B")
# Per-panel display cap on N (Thyroid is already saturated by 320, so the
# 400 point only stretches the axis).
n_cap <- c(Inf, 320)
for (pi in seq_along(panels)) {
  p      <- panels[[pi]]
  keep   <- p$N_grid <= n_cap[pi]
  N_grid <- p$N_grid[keep]
  tp     <- (100 * apply(p$plugin$marginal_power, 1, mean, na.rm = TRUE))[keep]
  bt_mn  <- (100 * p$boot$power_mean[, 1, 1])[keep]
  bt_lo  <- (100 * p$boot$power_ci_lo[, 1, 1])[keep]
  bt_hi  <- (100 * p$boot$power_ci_hi[, 1, 1])[keep]

  plot(N_grid, tp, type = "l", lwd = 2.0, col = col_plugin,
       ylim = c(0, 100), xlim = range(N_grid),
       xlab = "Sample size (n)", ylab = "Power (%)",
       main = "", xaxt = "n")
  # Explicit ticks at EVERY simulated N (interval is whatever the panel used),
  # so 120 and all other sampled sizes are labelled. gap.axis = -1 disables R's
  # automatic label thinning; a smaller tick font keeps them from colliding.
  axis(1, at = N_grid, labels = N_grid, cex.axis = 0.82, gap.axis = -1)
  # Panel letter prefixed into the centred title so the letter and the title
  # sit on the same line at the same size (matches Fig 1/2/4).
  title(main = sprintf("%s   %s", letters_seq[pi], p$label),
        line = 0.5, cex.main = 1.02, font.main = 2)
  abline(h = 80, lty = 2, col = "grey55", lwd = 1.2)
  lines(N_grid, bt_mn, lwd = 1.2, col = col_boot, lty = 2)
  arrows(N_grid, bt_lo, N_grid, bt_hi, code = 3, angle = 90,
         length = 0.04, lwd = 1.6, col = col_boot)
  points(N_grid, bt_mn, pch = 19, col = col_boot, cex = 0.7)
  if (pi == 2) {
    legend("bottomright",
           legend = c("Point estimate", "Bootstrap mean ± 95% CI"),
           col = c(col_plugin, col_boot), lwd = c(2.0, 1.6),
           lty = c(1, NA), pch = c(NA, 19),
           bty = "o", box.col = "grey55", box.lwd = 0.8, bg = "white",
           cex = 0.58, inset = 0.015, seg.len = 0.8,
           y.intersp = 0.85, x.intersp = 0.35)
  }
}
mtext("Bootstrap Uncertainty in Single-Cohort Power Estimates",
      outer = TRUE, side = 3, line = 0.25, font = 2, cex = 1.10, adj = 0.5)
dev.off()

for (d in dests) {
  dir.create(dirname(d), recursive = TRUE, showWarnings = FALSE)
  file.copy(tmp, d, overwrite = TRUE)
  cat("Saved:", d, "\n")
}
cat("Panel A N:", paste(panels[[1]]$N_grid, collapse = ","), "\n")
cat("Panel B N:", paste(panels[[2]]$N_grid, collapse = ","), "\n")
