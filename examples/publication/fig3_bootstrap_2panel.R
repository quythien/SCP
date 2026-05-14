#' =======================================================================
#' fig3_bootstrap_2panel.R  -  Trim Fig 3 to the two extreme-pilot panels
#'
#' Reads the cached panel list (Putamen SCZ n=28, Seney ACC n=60,
#' GTEx Pancreas n=249, GTEx Thyroid n=416) and emits a 2-panel figure
#' showing only the smallest (Putamen, n=28) and largest (Thyroid,
#' n=416) pilots, removing the Seney panel (in-house data with no
#' public accession) and the Pancreas panel (kept the largest-contrast
#' framing). Replaces the previous 4-panel Fig 3.
#' =======================================================================

setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
source("examples/publication/_pub_style.R")

panels_all <- readRDS("output/bootstrap_sc/results/all_panels.rds")
stopifnot(length(panels_all) == 4)

panels <- panels_all[c(1, 4)]   # Putamen + Thyroid

out_dir <- "output/bootstrap_sc"
fig_path_local <- file.path(out_dir, "figures", "fig_bootstrap_sc_2panel.pdf")
fig_path_main  <- "output/main_figures/Fig3_bootstrap_singlecohort.pdf"
fig_path_sub   <- "submission/figures/Fig3_bootstrap_singlecohort.pdf"
dir.create(dirname(fig_path_local), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(fig_path_main),  recursive = TRUE, showWarnings = FALSE)

pal_DT     <- pub_palette_detector()
col_plugin <- unname(pal_DT["DCP"])
col_boot   <- unname(pal_DT["FMM"])

for (out_pdf in c(fig_path_local, fig_path_main, fig_path_sub)) {
  cairo_pdf(out_pdf, width = 7.2, height = 3.4)
  pub_par(mfrow = c(1, 2), mar = c(4.0, 4.2, 2.4, 1.0))
  letters_seq <- c("a", "b")

  for (pi in seq_along(panels)) {
    p      <- panels[[pi]]
    N_grid <- p$N_grid
    tp     <- apply(p$plugin$marginal_power, 1, mean, na.rm = TRUE)
    bt_mn  <- p$boot$power_mean[, 1, 1]
    bt_lo  <- p$boot$power_ci_lo[, 1, 1]
    bt_hi  <- p$boot$power_ci_hi[, 1, 1]

    plot(N_grid, tp, type = "l", lwd = 1.8, col = col_plugin,
         ylim = c(0, 1), xlab = "N (total samples)", ylab = "Power",
         main = p$label)
    panel_label(letters_seq[pi])
    abline_80pct()
    lines(N_grid, bt_mn, lwd = 1.0, col = col_boot, lty = 2)
    arrows(N_grid, bt_lo, N_grid, bt_hi,
           code = 3, angle = 90, length = 0.04, lwd = 1.6, col = col_boot)
    points(N_grid, bt_mn, pch = 19, col = col_boot, cex = 0.7)
    if (pi == 1) {
      pub_legend("bottomright",
                 legend = c("Plug-in", "Bootstrap mean + 95% CI"),
                 col = c(col_plugin, col_boot),
                 lwd = c(1.8, 1.6), lty = c(1, NA),
                 pch = c(NA, 19))
    }
  }

  dev.off()
  cat(sprintf("Saved: %s\n", out_pdf))
}
