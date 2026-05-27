#' Render Fig 6 MH from cache. 1x3 layout mirroring fig6_active_replot.R:
#'   Panel A -- DCP detection power on active KIM (data from fig5_v4 cache, FMM truth)
#'   Panel B -- median phase MSE under MH truth, vs N with curves per B + 95% CI
#'   Panel C -- K = 2 detection power on active KIM (data from fig5_v4 cache, FMM truth)
#'
#' Note: Panels A/C still use FMM-truth power data because no MH-truth
#' active-power sweep was run; only the phase MSE panel (B) uses MH truth.
#' This mirrors the existing fig6_active_replot.R layout and lets the
#' reader compare phase MSE under MH vs FMM truth.
#'
#' Output: submission/figures/mh/Fig6_mh_active_design.pdf (+ mirrors)

setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
source("examples/publication/_pub_style.R")

# All three panels: load the latest MH cache (phase MSE + power_dcp + power_k2)
mse_cache <- rev(sort(list.files(
  "output/mh_alternative/results",
  "^fig6_mh_phase_mse_.*\\.rds$", full.names = TRUE)))[1]
cat(sprintf("Loading MH cache: %s\n", mse_cache))
m <- readRDS(mse_cache)
B_GRID_MSE  <- m$B_GRID
N_GRID_MSE  <- m$N_GRID
B_GRID_AC   <- m$B_GRID_PWR              # power skips B = 4 (K=2 not identifiable)
N_GRID_AC   <- m$N_GRID
power_dcp_M <- m$power_dcp               # [B x N]
power_k2_M  <- m$power_k2                # [B x N]

# Mean + 95% CI of per-replicate median squared phase error
med_mse <- matrix(NA_real_, nrow = length(B_GRID_MSE), ncol = length(N_GRID_MSE),
                  dimnames = list(paste0("B=", B_GRID_MSE),
                                  paste0("N=", N_GRID_MSE)))
med_lo  <- med_mse; med_hi <- med_mse
for (i in seq_along(B_GRID_MSE)) {
  for (j in seq_along(N_GRID_MSE)) {
    cell <- m$results[[sprintf("B%d_N%d", B_GRID_MSE[i], N_GRID_MSE[j])]]
    if (is.null(cell)) next
    raw  <- cell$med_dcp_all_raw
    if (!is.null(raw) && sum(is.finite(raw)) >= 2) {
      mu  <- mean(raw, na.rm = TRUE)
      se  <- sd(raw, na.rm = TRUE) / sqrt(sum(is.finite(raw)))
      med_mse[i, j] <- mu
      med_lo[i, j]  <- mu - 1.96 * se
      med_hi[i, j]  <- mu + 1.96 * se
    } else {
      med_mse[i, j] <- cell$med_dcp_all_mean
    }
  }
}

pal_B_AC  <- pub_palette_sequential(length(B_GRID_AC))
pal_B_MSE <- pub_palette_sequential(length(B_GRID_MSE))

out_paths <- c(
  "output/mh_alternative/figures/Fig6_mh_active_design.pdf",
  "submission/figures/mh/Fig6_mh_active_design.pdf"
)
for (q in out_paths) dir.create(dirname(q), recursive = TRUE, showWarnings = FALSE)

draw <- function(out_pdf) {
  cairo_pdf(out_pdf, width = 11.5, height = 4.9)
  par(mfrow = c(1, 3),
      mai = c(1.30, 1.15, 0.85, 0.18),
      mgp = c(3.6, 0.75, 0),
      oma = c(0, 0, 2.4, 0),
      family = "Helvetica", las = 1, tcl = -0.3,
      cex.axis = 1.35, cex.lab = 1.55, cex.main = 1.55, font.main = 2)

  # Panel A: DCP power on active KIM (MH truth)
  plot(NA, xlim = range(N_GRID_AC), ylim = c(0, 1),
       xlab = expression(N ~ "(total samples)"),
       ylab = "Power (DCP, FDR = 0.05)",
       main = expression(bold("Active KIM,  cosinor F-test")))
  panel_label("A")
  abline_80pct()
  for (k in seq_along(B_GRID_AC)) {
    yk <- power_dcp_M[k, ]
    lines(N_GRID_AC, yk, col = pal_B_AC[k], lwd = 1.8)
    points(N_GRID_AC, yk, col = pal_B_AC[k], pch = 19, cex = 0.55)
  }
  legend("bottomright",
         legend = sapply(B_GRID_AC, function(b) as.expression(bquote(B == .(b)))),
         col = pal_B_AC, lwd = 1.8, pch = 19,
         bty = "o", bg = "white", cex = 0.95)

  # Panel B: phase MSE under MH truth
  y_max <- max(med_hi, med_mse, na.rm = TRUE) * 1.05
  plot(NA, xlim = range(N_GRID_MSE), ylim = c(0, y_max),
       xlab = expression(N ~ "(total samples)"),
       ylab = expression("Median phase MSE (h"^2 * ")"),
       main = expression(bold("Phase estimation under MH truth")))
  panel_label("B")
  for (i in seq_along(B_GRID_MSE)) {
    yi  <- med_mse[i, ]; ylo <- med_lo[i, ]; yhi <- med_hi[i, ]
    if (all(is.na(yi))) next
    lines(N_GRID_MSE, yi, col = pal_B_MSE[i], lwd = 1.8)
    points(N_GRID_MSE, yi, col = pal_B_MSE[i], pch = 19, cex = 0.7)
    keep <- is.finite(yi) & is.finite(ylo) & is.finite(yhi)
    if (any(keep)) {
      arrows(N_GRID_MSE[keep], ylo[keep], N_GRID_MSE[keep], yhi[keep],
             angle = 90, code = 3, length = 0.03,
             col = pal_B_MSE[i], lwd = 1.0)
    }
  }
  legend("topright",
         legend = sapply(B_GRID_MSE, function(b) as.expression(bquote(B == .(b)))),
         col = pal_B_MSE, lwd = 1.8, pch = 19,
         bty = "o", bg = "white", cex = 0.95)

  # Panel C: K = 2 power on active KIM (MH truth, correctly specified)
  plot(NA, xlim = range(N_GRID_AC), ylim = c(0, 1),
       xlab = expression(N ~ "(total samples)"),
       ylab = "Power (K = 2, FDR = 0.05)",
       main = expression(bold("Active KIM,  K = 2 harmonic F-test")))
  panel_label("C")
  abline_80pct()
  for (k in seq_along(B_GRID_AC)) {
    yk <- power_k2_M[k, ]
    lines(N_GRID_AC, yk, col = pal_B_AC[k], lwd = 1.8)
    points(N_GRID_AC, yk, col = pal_B_AC[k], pch = 19, cex = 0.55)
  }
  legend("bottomright",
         legend = sapply(B_GRID_AC, function(b) as.expression(bquote(B == .(b)))),
         col = pal_B_AC, lwd = 1.8, pch = 19,
         bty = "o", bg = "white", cex = 0.95)

  mtext("Active Design Analysis: MH Phase MSE (Baboon KIM, Active)",
        outer = TRUE, cex = 1.5, font = 2)
  dev.off()
  cat(sprintf("Saved: %s\n", out_pdf))
}

for (q in out_paths) draw(q)
