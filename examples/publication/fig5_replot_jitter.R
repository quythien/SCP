#' Replot Fig 5 with per-B linetype + small x-jitter so the four B curves
#' separate visually instead of overlapping into one band.
setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
source("examples/publication/_pub_style.R")

cache_path <- rev(list.files(
  "output/active_vs_passive/results",
  "^fig5_v4_active_passive_.*\\.rds$", full.names = TRUE))[1]
r <- readRDS(cache_path)
B_GRID   <- r$B_GRID_ACTIVE
N_GRID   <- r$N_GRID_ACTIVE
K_HARM   <- r$K_HARM

# Per-B visual encoding: color (already sequential), linetype (1..4), small x-jitter
n_B       <- length(B_GRID)
pal_B     <- pub_palette_sequential(n_B)
lty_B     <- c(1, 2, 3, 4)[seq_len(n_B)]
x_jitter  <- seq(-1.5, 1.5, length.out = n_B)
pch_B     <- c(19, 17, 15, 18)[seq_len(n_B)]

fig_paths <- c(
  "output/active_vs_passive/figures/fig5_v4_active_vs_passive.pdf",
  "output/main_figures/Fig5_active_vs_passive.pdf",
  "submission/figures/Fig5_active_vs_passive.pdf"
)

draw_panel <- function(results_block, main_label, panel_letter,
                       legend_in_this_panel = FALSE) {
  plot(NA, xlim = range(N_GRID), ylim = c(0, 1),
       xlab = "N (total samples)", ylab = "Power",
       main = main_label)
  panel_label(panel_letter)
  abline_80pct()
  for (k in seq_along(B_GRID)) {
    df <- results_block[[as.character(B_GRID[k])]]
    xv <- df$N + x_jitter[k]
    lines(xv, df$power, col = pal_B[k], lwd = 1.8, lty = lty_B[k])
    points(xv, df$power, pch = pch_B[k], col = pal_B[k], cex = 0.7)
  }
  if (legend_in_this_panel) {
    pub_legend("bottomright",
               legend = sprintf("%d", B_GRID),
               col   = pal_B,
               lty   = lty_B,
               pch   = pch_B,
               lwd   = 1.6,
               title = "B (timepoints)")
  }
}

for (out_pdf in fig_paths) {
  cairo_pdf(out_pdf, width = 7.2, height = 7.0)
  pub_par(mfrow = c(2, 1), mar = c(4.0, 4.2, 2.4, 1.0))
  draw_panel(r$results_active$DCP, "Active, DCP",                 "a", TRUE)
  draw_panel(r$results_active$FMM, sprintf("Active, K = %d", K_HARM), "b", FALSE)
  dev.off()
  cat(sprintf("Saved: %s\n", out_pdf))
}
