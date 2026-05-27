#' Render Fig 4 MH from cache. 2x3 layout:
#'   Top row (illustrations of the MH signal):
#'     A1 -- vary alpha2 (alpha3 = 0, peak at noon)
#'     A2 -- vary alpha3 (alpha2 = anchor, peak at noon)
#'     A3 -- vary phi  (alpha2 = anchor, alpha3 = 0)
#'   Bottom row (DCP power sweeps under MH truth):
#'     B  -- alpha2 sweep (alpha3 = 0)
#'     C  -- alpha3 sweep (alpha2 = anchor)
#'     D  -- phi-dispersion sweep (alpha2 = anchor, alpha3 = 0)
#'
#' Output: submission/figures/mh/Fig4_mh_sensitivity.pdf (+ mirror)

setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
source("examples/publication/_pub_style.R")

cache_dir  <- "output/mh_alternative/results"
cache_path <- rev(sort(list.files(cache_dir,
                                  "^fig4_mh_sensitivity_.*\\.rds$",
                                  full.names = TRUE)))[1]
cat(sprintf("Loading: %s\n", cache_path))
r <- readRDS(cache_path)

N_GRID       <- r$N_grid
ETA_A2_GRID  <- r$eta_a2_grid
ETA_A3_GRID  <- r$eta_a3_grid
SDHR_GRID    <- r$sdhr_grid
eta_a2_anc   <- r$eta_a2_anchor
power_a2     <- r$power_A
power_a3     <- r$power_B
power_phi    <- r$power_C
mean_a2_anc  <- 1 / (1 + eta_a2_anc)   # for illustration anchor

# Illustration helper: y(t) = cos(omega0 (t - phi)) + a2 cos(2 omega0 (t - phi)) + a3 cos(3 omega0 (t - phi))
mh_curve <- function(t_h, alpha2 = 0, alpha3 = 0, phi_h = 12) {
  omega0 <- 2 * pi / 24
  cos(omega0 * (t_h - phi_h)) +
    alpha2 * cos(2 * omega0 * (t_h - phi_h)) +
    alpha3 * cos(3 * omega0 * (t_h - phi_h))
}

ALPHA2_DEMO <- c(0.00, 0.10, 0.30, 0.50)
ALPHA3_DEMO <- c(0.00, 0.10, 0.30, 0.50)
PHI_DEMO    <- c(0, 6, 12, 18)
T_GRID_FINE <- seq(0, 24, length.out = 401)

out_paths <- c(
  "output/mh_alternative/figures/Fig4_mh_sensitivity.pdf",
  "submission/figures/mh/Fig4_mh_sensitivity.pdf"
)
for (p in out_paths) dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)

pal_A2_demo <- pub_palette_sequential(length(ALPHA2_DEMO))
pal_A3_demo <- pub_palette_sequential(length(ALPHA3_DEMO))
pal_PHI     <- pub_palette_sequential(length(PHI_DEMO))
pal_a2      <- pub_palette_sequential(length(ETA_A2_GRID))
pal_a3      <- pub_palette_sequential(length(ETA_A3_GRID))
pal_sda     <- pub_palette_sequential(length(SDHR_GRID))
idx_a2_anchor  <- which.min(abs(ETA_A2_GRID - eta_a2_anc))
idx_a3_anchor  <- which.max(ETA_A3_GRID)         # large eta_a3 = near-cosinor baseline
idx_sda_anchor <- 1L

draw <- function(out_pdf) {
  cairo_pdf(out_pdf, width = 13.5, height = 8.5)
  par(mfrow = c(2, 3),
      mai = c(1.20, 1.10, 0.85, 0.25),
      mgp = c(3.6, 0.75, 0),
      oma = c(0, 0, 2.6, 0),
      family = "Helvetica", las = 1, tcl = -0.3,
      cex.axis = 1.20, cex.lab = 1.45, cex.main = 1.45, font.main = 2)

  # ------- TOP ROW: signal illustrations -------
  plot(NA, xlim = c(0, 24), ylim = c(-1.8, 1.8),
       xlab = "Time (h)", ylab = "MH signal y(t)",
       main = bquote(bold("Vary " ~ alpha[2] ~ "  (" * alpha[3] ~ "= 0, " * phi ~ "= 12 h)")),
       xaxs = "i")
  panel_label("A")
  abline(h = 0, col = "grey90", lwd = 0.6)
  for (k in seq_along(ALPHA2_DEMO))
    lines(T_GRID_FINE,
          mh_curve(T_GRID_FINE, alpha2 = ALPHA2_DEMO[k], alpha3 = 0, phi_h = 12),
          col = pal_A2_demo[k], lwd = 2.0)
  axis(1, at = c(0, 6, 12, 18, 24))
  pub_legend("bottomright",
             legend = sapply(ALPHA2_DEMO, function(a)
               as.expression(bquote(alpha[2] == .(format(a, nsmall = 2))))),
             col = pal_A2_demo, lwd = 1.8, title = NULL)

  plot(NA, xlim = c(0, 24), ylim = c(-1.8, 1.8),
       xlab = "Time (h)", ylab = "MH signal y(t)",
       main = bquote(bold("Vary " ~ alpha[3] ~ "  (mean " * alpha[2] * "=" * .(format(mean_a2_anc, digits = 2)) * ", " * phi ~ "= 12 h)")),
       xaxs = "i")
  panel_label("B")
  abline(h = 0, col = "grey90", lwd = 0.6)
  for (k in seq_along(ALPHA3_DEMO))
    lines(T_GRID_FINE,
          mh_curve(T_GRID_FINE, alpha2 = mean_a2_anc, alpha3 = ALPHA3_DEMO[k], phi_h = 12),
          col = pal_A3_demo[k], lwd = 2.0)
  axis(1, at = c(0, 6, 12, 18, 24))
  pub_legend("bottomright",
             legend = sapply(ALPHA3_DEMO, function(a)
               as.expression(bquote(alpha[3] == .(format(a, nsmall = 2))))),
             col = pal_A3_demo, lwd = 1.8, title = NULL)

  plot(NA, xlim = c(0, 24), ylim = c(-1.8, 1.8),
       xlab = "Time (h)", ylab = "MH signal y(t)",
       main = bquote(bold("Vary " ~ phi ~ "  (mean " * alpha[2] * "=" * .(format(mean_a2_anc, digits = 2)) * ", " * alpha[3] ~ "= 0)")),
       xaxs = "i")
  panel_label("C")
  abline(h = 0, col = "grey90", lwd = 0.6)
  for (k in seq_along(PHI_DEMO))
    lines(T_GRID_FINE,
          mh_curve(T_GRID_FINE, alpha2 = mean_a2_anc, alpha3 = 0, phi_h = PHI_DEMO[k]),
          col = pal_PHI[k], lwd = 2.0)
  axis(1, at = c(0, 6, 12, 18, 24))
  pub_legend("bottomright",
             legend = sapply(PHI_DEMO, function(p)
               as.expression(bquote(phi == .(p) ~ "h"))),
             col = pal_PHI, lwd = 1.8, title = NULL)

  # ------- BOTTOM ROW: DCP power sweeps -------
  plot(NA, xlim = range(N_GRID), ylim = c(0, 1),
       xlab = expression(N ~ "(total samples)"),
       ylab = "Power (DCP, FDR = 0.05)",
       main = bquote(bold(alpha[2*"g"] ~ "~ Beta(1, " ~ eta[alpha[2]] ~ ")")))
  panel_label("D")
  abline_80pct()
  for (k in seq_along(ETA_A2_GRID)) {
    lwd_k <- if (k == idx_a2_anchor) 2.6 else 1.7
    lines(N_GRID, power_a2[, k], col = pal_a2[k], lwd = lwd_k)
    points(N_GRID, power_a2[, k], col = pal_a2[k], pch = 19, cex = 0.6)
  }
  legend("bottomright",
         legend = sapply(ETA_A2_GRID, function(b)
           as.expression(bquote(eta[alpha[2]] == .(b)))),
         col = pal_a2, lwd = 1.7,
         bty = "o", bg = "white", cex = 0.85)

  plot(NA, xlim = range(N_GRID), ylim = c(0, 1),
       xlab = expression(N ~ "(total samples)"),
       ylab = "Power (DCP, FDR = 0.05)",
       main = bquote(bold(alpha[3*"g"] ~ "~ Beta(1, " ~ eta[alpha[3]] ~ ")")))
  panel_label("E")
  abline_80pct()
  for (k in seq_along(ETA_A3_GRID)) {
    lwd_k <- if (k == idx_a3_anchor) 2.6 else 1.7
    lines(N_GRID, power_a3[, k], col = pal_a3[k], lwd = lwd_k)
    points(N_GRID, power_a3[, k], col = pal_a3[k], pch = 19, cex = 0.6)
  }
  legend("bottomright",
         legend = sapply(ETA_A3_GRID, function(b)
           as.expression(bquote(eta[alpha[3]] == .(b)))),
         col = pal_a3, lwd = 1.7,
         bty = "o", bg = "white", cex = 0.85)

  plot(NA, xlim = range(N_GRID), ylim = c(0, 1),
       xlab = expression(N ~ "(total samples)"),
       ylab = "Power (DCP, FDR = 0.05)",
       main = bquote(bold(phi[g] ~ "~" ~ vonMises(0, kappa[phi]) ~ "sweep")))
  panel_label("F")
  abline_80pct()
  for (k in seq_along(SDHR_GRID)) {
    lwd_k <- if (k == idx_sda_anchor) 2.6 else 1.7
    lines(N_GRID, power_phi[, k], col = pal_sda[k], lwd = lwd_k)
    points(N_GRID, power_phi[, k], col = pal_sda[k], pch = 19, cex = 0.6)
  }
  legend("bottomright",
         legend = sapply(SDHR_GRID, function(s)
           as.expression(bquote(sigma[phi] == .(s) ~ "h"))),
         col = pal_sda, lwd = 1.7,
         bty = "o", bg = "white", cex = 0.95)

  mtext("MH Sensitivity Analysis: Cosinor-Based DCP under Multi-Harmonic Truth (Baboon LUN, Active)",
        outer = TRUE, cex = 1.4, font = 2)
  dev.off()
  cat(sprintf("Saved: %s\n", out_pdf))
}

for (p in out_paths) draw(p)
