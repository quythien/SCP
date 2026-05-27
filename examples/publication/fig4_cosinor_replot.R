#' =======================================================================
#' fig4_cosinor_replot.R — Render Fig 4 from cached cosinor-sensitivity RDS
#' =======================================================================
#'
#' 2 x 2 layout:
#'   Top row, Panel A
#'     Left  — FMM waveform shape vs omega (alpha fixed at 12 h)
#'     Right — FMM peak shift vs alpha (omega fixed at 0.4)
#'   Bottom row
#'     Panel B — eta sweep: DCP power vs N at varying omega ~ Beta(1, eta)
#'     Panel C — sigma_alpha sweep: DCP power vs N at varying peak dispersion
#'
#' INPUT:  output/sensitivity/results/fig4_cosinor_sensitivity_*.rds (latest)
#' OUTPUT: submission/figures/Fig4_sensitivity.pdf (+ mirrors)
#'
#' USAGE:  Rscript examples/publication/fig4_cosinor_replot.R

setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
source("examples/publication/_pub_style.R")

cache_dir  <- "output/sensitivity/results"
cache_path <- rev(sort(list.files(
  cache_dir, "^fig4_cosinor_sensitivity_.*\\.rds$", full.names = TRUE)))[1]
if (is.na(cache_path)) stop("No fig4_cosinor_sensitivity_*.rds in ", cache_dir)
cat(sprintf("Loading: %s\n", cache_path))
r <- readRDS(cache_path)

N_GRID    <- r$N_grid
BETA_GRID <- r$beta_grid           # eta values for Panel B
SDHR_GRID <- r$sdhr_grid           # sigma_alpha values (hours) for Panel C
eta_hat       <- r$beta_hat
sd_alpha_hat  <- r$sigma_alpha_hat
power_eta <- r$power_A             # [N x eta]
power_sda <- r$power_B             # [N x sigma_alpha]

# FMM waveform helper: y(t) = cos(beta + 2*atan(omega * tan((t_rad - alpha_rad)/2))).
fmm_curve <- function(t_h, omega, alpha_h, beta = 0) {
  t_rad     <- t_h * 2 * pi / 24
  alpha_rad <- alpha_h * 2 * pi / 24
  cos(beta + 2 * atan(omega * tan((t_rad - alpha_rad) / 2)))
}

OMEGA_DEMO  <- c(1.00, 0.60, 0.30, 0.15)   # omega varies, alpha = 12
ALPHA_DEMO  <- c(0, 6, 12, 18)             # alpha varies, omega = 0.4
T_GRID_FINE <- seq(0, 24, length.out = 401)

out_paths <- c(
  "output/sensitivity/figures/Fig4_sensitivity.pdf",
  "output/main_figures/Fig4_sensitivity.pdf",
  "submission/figures/Fig4_sensitivity.pdf"
)
for (p in out_paths) dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)

pal_OMEGA <- pub_palette_sequential(length(OMEGA_DEMO))
pal_ALPHA <- pub_palette_sequential(length(ALPHA_DEMO))
pal_eta   <- pub_palette_sequential(length(BETA_GRID))
pal_sda   <- pub_palette_sequential(length(SDHR_GRID))

idx_eta_anchor <- which.min(abs(BETA_GRID  - eta_hat))
idx_sda_anchor <- which.min(abs(SDHR_GRID  - sd_alpha_hat))

draw_fig4 <- function(out_pdf) {
  cairo_pdf(out_pdf, width = 9.0, height = 7.2)
  pub_par(mfrow = c(2, 2), mar = c(4.0, 4.5, 2.6, 1.0))

  # --- Panel A-left: omega illustration ---
  plot(NA, xlim = c(0, 24), ylim = c(-1.25, 1.25),
       xlab = "Time (h)", ylab = "FMM signal y(t)",
       main = expression(bold("Vary " * omega * "  (peak at " * alpha == 12 ~ "h)")),
       xaxs = "i")
  panel_label("A")
  abline(h = 0, col = "grey90", lwd = 0.6)
  for (k in seq_along(OMEGA_DEMO)) {
    yv <- fmm_curve(T_GRID_FINE, omega = OMEGA_DEMO[k], alpha_h = 12)
    lines(T_GRID_FINE, yv, col = pal_OMEGA[k], lwd = 2.0)
  }
  axis(1, at = c(0, 6, 12, 18, 24))
  pub_legend("bottomright",
             legend = sapply(OMEGA_DEMO, function(o)
               as.expression(bquote(omega == .(format(o, nsmall = 2))))),
             col   = pal_OMEGA, lwd = 1.8, title = NULL)

  # --- Panel A-right: alpha illustration ---
  plot(NA, xlim = c(0, 24), ylim = c(-1.25, 1.25),
       xlab = "Time (h)", ylab = "FMM signal y(t)",
       main = expression(bold("Vary " * alpha * "  (shape " * omega == 0.4 * ")")),
       xaxs = "i")
  abline(h = 0, col = "grey90", lwd = 0.6)
  for (k in seq_along(ALPHA_DEMO)) {
    yv <- fmm_curve(T_GRID_FINE, omega = 0.4, alpha_h = ALPHA_DEMO[k])
    lines(T_GRID_FINE, yv, col = pal_ALPHA[k], lwd = 2.0)
  }
  axis(1, at = c(0, 6, 12, 18, 24))
  pub_legend("bottomright",
             legend = sapply(ALPHA_DEMO, function(a)
               as.expression(bquote(alpha == .(a) ~ "h"))),
             col   = pal_ALPHA, lwd = 1.8, title = NULL)

  # --- Panel B: eta sweep ---
  plot(NA, xlim = range(N_GRID), ylim = c(0, 1),
       xlab = expression(N ~ "(total samples)"),
       ylab = "Power (DCP, FDR = 0.05)",
       main = expression(bold(omega[g] ~ "~" ~ Beta(1, eta) * "  (waveform shape sweep)")))
  panel_label("B")
  abline_80pct()
  for (k in seq_along(BETA_GRID)) {
    lwd_k <- if (k == idx_eta_anchor) 2.6 else 1.7
    lty_k <- if (BETA_GRID[k] == 0) 2 else 1
    lines(N_GRID, power_eta[, k], col = pal_eta[k], lwd = lwd_k, lty = lty_k)
    points(N_GRID, power_eta[, k], col = pal_eta[k], pch = 19, cex = 0.6)
  }
  pub_legend("bottomright",
             legend = sapply(BETA_GRID, function(b)
               if (b == 0) expression(eta == 0 ~ "(cosinor)")
               else as.expression(bquote(eta == .(b)))),
             col   = pal_eta, lwd = 1.7,
             lty   = ifelse(BETA_GRID == 0, 2, 1),
             title = NULL)

  # --- Panel C: sigma_alpha sweep ---
  plot(NA, xlim = range(N_GRID), ylim = c(0, 1),
       xlab = expression(N ~ "(total samples)"),
       ylab = "Power (DCP, FDR = 0.05)",
       main = expression(bold(alpha[g] ~ "~" ~ vonMises(0, kappa[alpha]) * "  (peak-dispersion sweep)")))
  panel_label("C")
  abline_80pct()
  for (k in seq_along(SDHR_GRID)) {
    lwd_k <- if (k == idx_sda_anchor) 2.6 else 1.7
    lines(N_GRID, power_sda[, k], col = pal_sda[k], lwd = lwd_k)
    points(N_GRID, power_sda[, k], col = pal_sda[k], pch = 19, cex = 0.6)
  }
  pub_legend("bottomright",
             legend = sapply(SDHR_GRID, function(s)
               as.expression(bquote(sigma[alpha] == .(s) ~ "h"))),
             col   = pal_sda, lwd = 1.7, title = NULL)

  dev.off()
  cat(sprintf("Saved: %s\n", out_pdf))
}

for (p in out_paths) draw_fig4(p)
