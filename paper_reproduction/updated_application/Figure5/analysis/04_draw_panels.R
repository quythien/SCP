# Original-style A/B review after production: Rscript 04_draw_panels.R REPORT_DIRECTORY NEW_FIGURE_DIRECTORY

# 1. Read saved numbers and require a new figure directory -----------------
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
report <- normalizePath(args[1], mustWork = TRUE)
destination <- normalizePath(args[2], mustWork = FALSE)
stopifnot(!dir.exists(destination))
cache <- readRDS(file.path(report, "figure_2_5_AB_cache.rds"))
stopifnot(isTRUE(cache$balanced_B_curves_coincide), isTRUE(cache$nonidentifiable_designs_not_included))
dir.create(destination, recursive = TRUE)

# 2. Retain the original B colors, axes and reference line without jitter --
colors <- c("3" = "#440154", "4" = "#3b528b", "6" = "#21918c",
            "8" = "#5ec962", "12" = "#fde725", "24" = "#7d028c")
draw_panels <- function() {
  par(mfrow = c(1, 2), mar = c(3.8, 3.5, 2.3, .6), mgp = c(2, .6, 0),
      tcl = -.25, cex.axis = .82, cex.lab = .85)
  for (panel in c("A", "B")) {
    B <- if (panel == "A") cache$B_K1 else cache$B_K2
    M <- if (panel == "A") cache$pwr_K1 else cache$pwr_K2
    SE <- if (panel == "A") cache$se_K1 else cache$se_K2
    plot(NA, xlim = c(0, max(cache$N) * 1.04), ylim = c(0, 100),
      xlab = "Sample size (N)", ylab = if (panel == "A") "Power (%)" else "",
      main = paste(panel, if (panel == "A") "Single-harmonic (K=1)" else "Two-harmonic (K=2)"),
      cex.main = .9)
    grid(lwd = .5, col = "grey85")
    abline(h = 80, lty = 2, col = "grey50", lwd = .9)
    for (j in seq_along(B)) {
      y <- 100 * M[, j]
      se <- 100 * SE[, j]
      ok <- se > 0
      arrows(cache$N[ok], pmax(0, y[ok] - se[ok]), cache$N[ok], pmin(100, y[ok] + se[ok]),
             code = 3, angle = 90, length = .015, col = colors[as.character(B[j])], lwd = .7)
      lines(cache$N, y, type = "o", pch = 19, cex = .35, lwd = 1.1,
            col = colors[as.character(B[j])])
    }
    legend("bottomright", paste0("B=", B), col = colors[as.character(B)],
           lty = 1, pch = 19, pt.cex = .5, cex = .64, bty = "n")
    mtext("Balanced B curves overlap", side = 3, line = .25, cex = .66)
  }
}
pdf(file.path(destination, "Figure_2_5_AB_review.pdf"), width = 86 / 25.4,
    height = 1.8, pointsize = 8, family = "Helvetica")
draw_panels()
dev.off()
png(file.path(destination, "Figure_2_5_AB_review.png"), width = 1016, height = 540,
    res = 300, type = if (capabilities("cairo")) "cairo" else getOption("bitmapType"))
draw_panels()
dev.off()

# 3. Explain exact overlap and the two separate truth calibrations ---------
writeLines(c("The B curves coincide under these exactly balanced Gaussian designs: the harmonic information matrix depends on total N and is identical across all displayed B values.",
  "Parameters and canonical Gaussian draws are shared across B within each panel. Each B cell has 600 replicates. Error bars show one Monte Carlo standard error for a fixed pilot.",
  "Panel A uses Figure 2.1's pure first-harmonic Liver calibration. Panel B uses the verified two-harmonic Liver calibration. Comparing panels does not isolate a detector advantage.",
  "No coordinates were jittered to separate the overlapping curves. Non-identifiable K2 designs B=3 and B=4 are unavailable and are not drawn as zero power."),
  file.path(destination, "Figure_notes.txt"))
cat("Saved original-style A/B review figures from production numbers only.\n")
