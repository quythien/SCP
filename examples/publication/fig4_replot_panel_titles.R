#' Redraw Fig 4 from cached RDS with the panel-title beta=0 dropped.
setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
source("examples/publication/_pub_style.R")

rds <- readRDS("output/sensitivity/results/fig4_sensitivity_20260514_034908.rds")

beta_hat        <- rds$beta_hat
sigma_alpha_hat <- rds$sigma_alpha_hat
N_GRID          <- rds$N_grid
BETA_GRID       <- rds$beta_grid
SDHR_GRID       <- rds$sdhr_grid
power_A         <- rds$power_A
power_B         <- rds$power_B

fmm_curve <- function(t_h, omega, alpha_h, beta = 0) {
  t_rad <- t_h * (2 * pi / 24)
  alpha_rad <- alpha_h * (2 * pi / 24)
  cos(beta + 2 * atan(omega * tan((t_rad - alpha_rad) / 2)))
}

OMEGA_DEMO <- c(1.00, 0.60, 0.30, 0.15)
ALPHA_DEMO <- c(0, 6, 12, 18)
T_GRID_FINE <- seq(0, 24, length.out = 401)

fig_paths <- c(
  "output/sensitivity/figures/fig4_sensitivity.pdf",
  "output/main_figures/Fig4_sensitivity.pdf",
  "submission/figures/Fig4_sensitivity.pdf"
)

for (out_pdf in fig_paths) {
  cairo_pdf(out_pdf, width = 7.2, height = 6.6)
  pub_par(mfrow = c(2, 2), mar = c(4.0, 4.2, 2.4, 1.0))

  pal_OMEGA <- pub_palette_sequential(length(OMEGA_DEMO))
  pal_ALPHA <- pub_palette_sequential(length(ALPHA_DEMO))
  pal_A     <- pub_palette_sequential(length(BETA_GRID))
  pal_B     <- pub_palette_sequential(length(SDHR_GRID))

  # Panel a: omega illustration
  plot(NA, xlim = c(0, 24), ylim = c(-1.2, 1.2),
       xlab = "Time (h)", ylab = "FMM signal y(t)",
       main = expression(omega ~ "shape (alpha = 12 h)"),
       xaxs = "i")
  panel_label("a")
  abline(h = 0, col = "grey85", lwd = 0.6)
  for (k in seq_along(OMEGA_DEMO)) {
    yv <- fmm_curve(T_GRID_FINE, omega = OMEGA_DEMO[k], alpha_h = 12)
    lines(T_GRID_FINE, yv, col = pal_OMEGA[k], lwd = 1.8)
  }
  pub_legend("bottomright",
             legend = sprintf("%.2f", OMEGA_DEMO),
             col = pal_OMEGA, lwd = 1.6,
             title = expression(omega),
             cex = 0.7)

  # Panel a (right half): alpha illustration
  plot(NA, xlim = c(0, 24), ylim = c(-1.2, 1.2),
       xlab = "Time (h)", ylab = "FMM signal y(t)",
       main = expression(alpha ~ "peak shift (omega = 0.4)"),
       xaxs = "i")
  abline(h = 0, col = "grey85", lwd = 0.6)
  for (k in seq_along(ALPHA_DEMO)) {
    yv <- fmm_curve(T_GRID_FINE, omega = 0.4, alpha_h = ALPHA_DEMO[k])
    lines(T_GRID_FINE, yv, col = pal_ALPHA[k], lwd = 1.8)
  }
  pub_legend("bottomright",
             legend = sprintf("%d h", ALPHA_DEMO),
             col = pal_ALPHA, lwd = 1.6,
             title = expression(alpha),
             cex = 0.7)

  # Panel b: omega sweep
  plot(NA, xlim = range(N_GRID), ylim = c(0, 1),
       xlab = "N (total samples)", ylab = "Power",
       main = expression(omega ~ "sweep"))
  panel_label("b")
  abline_80pct()
  for (k in seq_along(BETA_GRID)) {
    lty_k <- if (abs(BETA_GRID[k] - beta_hat) < 0.05) 1 else 2
    lwd_k <- if (abs(BETA_GRID[k] - beta_hat) < 0.05) 2.2 else 1.5
    lines(N_GRID, power_A[, k], col = pal_A[k], lwd = lwd_k, lty = lty_k)
    points(N_GRID, power_A[, k], pch = 19, col = pal_A[k], cex = 0.7)
  }
  pub_legend("bottomright",
             legend = sprintf("%g", BETA_GRID),
             col = pal_A, lwd = 1.6,
             title = expression(eta),
             cex = 0.7)

  # Panel c: sigma_alpha sweep
  plot(NA, xlim = range(N_GRID), ylim = c(0, 1),
       xlab = "N (total samples)", ylab = "Power",
       main = expression(sigma[alpha] ~ "(peak-phase dispersion) sweep"))
  panel_label("c")
  abline_80pct()
  for (k in seq_along(SDHR_GRID)) {
    lines(N_GRID, power_B[, k], col = pal_B[k], lwd = 1.6)
    points(N_GRID, power_B[, k], pch = 19, col = pal_B[k], cex = 0.7)
  }
  pub_legend("bottomright",
             legend = sprintf("%g", SDHR_GRID),
             col = pal_B, lwd = 1.6,
             title = expression(sigma[alpha] ~ "(h)"),
             cex = 0.7)

  dev.off()
  cat(sprintf("Saved: %s\n", out_pdf))
}
