#!/usr/bin/env Rscript
# Render Fig 6 both Options A and B for visual comparison
setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")

render_fig6 <- function(rds_path, out_pdf, title_suffix) {
  res <- readRDS(rds_path)
  pdf(out_pdf, width = 10.5, height = 3.7)
  par(mfrow = c(1, 3), mai = c(0.88, 0.78, 0.50, 0.10),
      mgp = c(2.6, 0.55, 0), oma = c(0, 0, 1.7, 0),
      cex.axis = 0.85, cex.lab = 1.00, cex.main = 0.95, font.main = 2,
      tcl = -0.35)

  N <- res$N
  cols_K1 <- viridis::viridis(length(res$B_K1))
  if (!requireNamespace("viridis", quietly = TRUE))
    cols_K1 <- c("#440154","#3b528b","#21918c","#5ec962","#fde725")[seq_along(res$B_K1)]
  cols_K2 <- c("#440154","#3b528b","#21918c","#5ec962","#fde725")[seq_along(res$B_K2)]

  # ---- Panel A: K=1 power ----
  matplot(N, 100 * res$pwr_K1, type = "b", pch = 19, lwd = 1.8,
          col = cols_K1, lty = 1, ylim = c(0, 100),
          xlim = c(0, max(N) * 1.05),
          xlab = "Sample size (N)", ylab = "K=1 Power (%)",
          main = "")
  title(main = "A   K=1 Power", adj = 0, font.main = 2,
        cex.main = 0.95, line = 0.4)
  abline(h = 80, lty = 2, col = "grey50", lwd = 1.2)
  grid()
  legend("bottomright", paste0("B=", res$B_K1),
         col = cols_K1, lty = 1, pch = 19, lwd = 1.5, cex = 0.55,
         bty = "o", box.col = "grey70", inset = 0.01, y.intersp = 0.78,
         bg = "white")

  # ---- Panel B: K=2 power ----
  matplot(N, 100 * res$pwr_K2, type = "b", pch = 19, lwd = 1.8,
          col = cols_K2, lty = 1, ylim = c(0, 100),
          xlim = c(0, max(N) * 1.05),
          xlab = "Sample size (N)", ylab = "K=2 Power (%)",
          main = "")
  title(main = "B   K=2 Power", adj = 0, font.main = 2,
        cex.main = 0.95, line = 0.4)
  abline(h = 80, lty = 2, col = "grey50", lwd = 1.2)
  grid()
  legend("bottomright", paste0("B=", res$B_K2),
         col = cols_K2, lty = 1, pch = 19, lwd = 1.5, cex = 0.55,
         bty = "o", box.col = "grey70", inset = 0.01, y.intersp = 0.78,
         bg = "white")

  # ---- Panel C: K=2 vs K=1 advantage (paired comparison) ----
  # Take the B values shared between K=1 and K=2 (B in {6,8,12,24})
  shared_B <- intersect(res$B_K1, res$B_K2)
  idx_K1 <- match(shared_B, res$B_K1)
  idx_K2 <- match(shared_B, res$B_K2)
  diff_mat <- 100 * (res$pwr_K2[, idx_K2] - res$pwr_K1[, idx_K1])
  matplot(N, diff_mat, type = "b", pch = 19, lwd = 1.8,
          col = cols_K2[idx_K2], lty = 1,
          ylim = range(diff_mat, na.rm = TRUE, finite = TRUE) + c(-5, 5),
          xlim = c(0, max(N) * 1.05),
          xlab = "Sample size (N)", ylab = "K=2 minus K=1 (pp)",
          main = "")
  title(main = "C   K=2 advantage", adj = 0, font.main = 2,
        cex.main = 0.95, line = 0.4)
  abline(h = 0, lty = 2, col = "grey50", lwd = 1.2)
  grid()
  legend("topright", paste0("B=", shared_B),
         col = cols_K2[idx_K2], lty = 1, pch = 19, lwd = 1.5, cex = 0.55,
         bty = "o", box.col = "grey70", inset = 0.01, y.intersp = 0.78,
         bg = "white")

  mtext(sprintf("Active design B-invariance: %s", title_suffix),
        outer = TRUE, side = 3, line = 0.2, font = 2, cex = 0.95)
  dev.off()
  cat(sprintf("Saved: %s\n", out_pdf))
}

if (file.exists("output/two_harmonic/results/fig6_optionA_Hughes_smallN.rds")) {
  render_fig6("output/two_harmonic/results/fig6_optionA_Hughes_smallN.rds",
              "submission/figures/Fig6_optionA_Hughes.pdf",
              "Hughes 2009 mouse liver (small-N)")
}
if (file.exists("output/two_harmonic/results/fig6_optionB_BaboonKIM.rds")) {
  render_fig6("output/two_harmonic/results/fig6_optionB_BaboonKIM.rds",
              "submission/figures/Fig6_optionB_BaboonKIM.pdf",
              "Baboon KIM")
}
