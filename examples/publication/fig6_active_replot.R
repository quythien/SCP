#' =======================================================================
#' fig6_active_replot.R — Render Fig 6 (active design: DCP / phase MSE / K=2)
#' =======================================================================
#'
#' Reads cached results from two sources and renders the 3-panel PDF:
#'   - Panel A (DCP power on active KIM) and
#'   - Panel C (K = 2 power on active KIM)
#'     come from output/active_vs_passive/results/fig5_v4_active_passive_*.rds
#'   - Panel B (median phase MSE under cosinor, B-sweep)
#'     comes from output/phase_mse/results/fig6_phase_mse_*.rds
#'
#' OUTPUT:
#'   output/main_figures/Fig6_active_design.pdf
#'   submission/figures/Fig6_active_design.pdf
#'
#' USAGE:
#'   Rscript examples/publication/fig6_active_replot.R

setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
source("examples/publication/_pub_style.R")

# ----- Load Panels A & C (active KIM, fig5_v4) -----
pwr_dir   <- "output/active_vs_passive/results"
pwr_cache <- rev(sort(list.files(
  pwr_dir, "^fig5_v4_active_passive_.*\\.rds$", full.names = TRUE)))[1]
if (is.na(pwr_cache)) stop("No fig5_v4_active_passive_*.rds in ", pwr_dir)
cat(sprintf("Loading power cache: %s\n", pwr_cache))
p <- readRDS(pwr_cache)

B_GRID_AC <- p$B_GRID_ACTIVE                # e.g. c(6, 8, 12, 24)
N_GRID_AC <- p$N_GRID_ACTIVE
res_dcp   <- p$results_active$DCP           # list keyed by B (character)
res_fmm   <- p$results_active$FMM

# ----- Load Panel B (phase MSE) -----
mse_dir   <- "output/phase_mse/results"
mse_cache <- rev(sort(list.files(
  mse_dir, "^fig6_phase_mse_.*\\.rds$", full.names = TRUE)))[1]
if (is.na(mse_cache)) stop("No fig6_phase_mse_*.rds in ", mse_dir)
cat(sprintf("Loading phase-MSE cache: %s\n", mse_cache))
m <- readRDS(mse_cache)

B_GRID_MSE <- m$B_GRID                       # e.g. c(4, 6, 8, 12, 24)
N_GRID_MSE <- m$N_GRID                       # e.g. c(48, 96)

# Extract median squared error (h^2) as a (B x N) matrix, with 95% CI on
# the mean of per-replicate medians: mean +/- 1.96 * SE,
# SE = sd(per-replicate medians) / sqrt(n_replicates).
med_mse <- matrix(NA_real_, nrow = length(B_GRID_MSE), ncol = length(N_GRID_MSE),
                  dimnames = list(paste0("B=", B_GRID_MSE),
                                  paste0("N=", N_GRID_MSE)))
med_lo  <- med_mse                  # mean - 1.96 * SE
med_hi  <- med_mse                  # mean + 1.96 * SE
for (i in seq_along(B_GRID_MSE)) {
  for (j in seq_along(N_GRID_MSE)) {
    key <- sprintf("B%d_N%d", B_GRID_MSE[i], N_GRID_MSE[j])
    cell <- m$results[[key]]
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

# ----- Aesthetics -----
pal_B_AC  <- pub_palette_sequential(length(B_GRID_AC))
pal_N_MSE <- pub_palette_sequential(length(N_GRID_MSE))

# ----- Output paths -----
out_paths <- c(
  "output/main_figures/Fig6_active_design.pdf",
  "submission/figures/Fig6_active_design.pdf"
)
for (q in out_paths) dir.create(dirname(q), recursive = TRUE, showWarnings = FALSE)

# ----- Draw -----
draw_fig6 <- function(out_pdf) {
  cairo_pdf(out_pdf, width = 11.5, height = 4.5)
  par(mfrow    = c(1, 3),
      mai      = c(1.30, 1.15, 0.85, 0.18),
      mgp      = c(3.6, 0.75, 0),
      oma      = c(0, 0, 0.4, 0),
      family   = "Helvetica", las = 1, tcl = -0.3,
      cex.axis = 1.35, cex.lab = 1.55, cex.main = 1.55, font.main = 2)

  # --- Panel A: DCP power vs N, curves per B (active KIM) ---
  plot(NA, xlim = range(N_GRID_AC), ylim = c(0, 1),
       xlab = expression(N ~ "(total samples)"),
       ylab = "Power (DCP, FDR = 0.05)",
       main = expression(bold("Active KIM,  cosinor F-test")))
  panel_label("A")
  abline_80pct()
  for (k in seq_along(B_GRID_AC)) {
    df <- res_dcp[[as.character(B_GRID_AC[k])]]
    lines(df$N, df$power, col = pal_B_AC[k], lwd = 1.8)
    points(df$N, df$power, col = pal_B_AC[k], pch = 19, cex = 0.55)
  }
  pub_legend("bottomright",
             legend = sapply(B_GRID_AC, function(b)
               as.expression(bquote(B == .(b)))),
             col   = pal_B_AC, lwd = 1.8,
             title = NULL)

  # --- Panel B: median phase MSE vs N, curves per B + 95% MC error bars ---
  pal_B_MSE <- pub_palette_sequential(length(B_GRID_MSE))
  y_max <- max(med_hi, med_mse, na.rm = TRUE) * 1.05
  plot(NA, xlim = range(N_GRID_MSE), ylim = c(0, y_max),
       xlab = expression(N ~ "(total samples)"),
       ylab = expression("Median phase MSE (h" ^ 2 * ")"),
       main = expression(bold("Phase estimation  (" * hat(phi)[g]^{cos} * ")")))
  panel_label("B")
  # Mean lines + vertical 95% MC error bars (whiskers from 2.5% to 97.5% across replicates)
  for (i in seq_along(B_GRID_MSE)) {
    yi  <- med_mse[i, ]
    ylo <- med_lo[i, ]
    yhi <- med_hi[i, ]
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
         legend = sapply(B_GRID_MSE, function(b)
           as.expression(bquote(B == .(b)))),
         col   = pal_B_MSE, lwd = 1.8, pch = 19,
         bty   = "o", bg = "white", cex = 0.85)

  # --- Panel C: K=2 power vs N, curves per B (active KIM) ---
  plot(NA, xlim = range(N_GRID_AC), ylim = c(0, 1),
       xlab = expression(N ~ "(total samples)"),
       ylab = "Power (K = 2, FDR = 0.05)",
       main = expression(bold("Active KIM,  K = 2 harmonic F-test")))
  panel_label("C")
  abline_80pct()
  for (k in seq_along(B_GRID_AC)) {
    df <- res_fmm[[as.character(B_GRID_AC[k])]]
    lines(df$N, df$power, col = pal_B_AC[k], lwd = 1.8)
    points(df$N, df$power, col = pal_B_AC[k], pch = 19, cex = 0.55)
  }
  pub_legend("bottomright",
             legend = sapply(B_GRID_AC, function(b)
               as.expression(bquote(B == .(b)))),
             col   = pal_B_AC, lwd = 1.8,
             title = NULL)

  dev.off()
  cat(sprintf("Saved: %s\n", out_pdf))
}

for (q in out_paths) draw_fig6(q)
