#' =======================================================================
#' fig5_fmm_replot.R — Render Fig 5 from cached full-FMM-framework RDS
#' =======================================================================
#'
#' 3-panel layout (mirrors Fig 1 styling for Panels A and B):
#'   Panel A — marginal power vs N, one line per BH-FDR alpha
#'   Panel B — power stratified by r-tilde at FDR = 0.05, one line per N
#'   Panel C — eta sweep: power vs N at varying omega_g ~ Beta(1, eta)
#'
#' INPUT:  output/fmm_framework/results/fig5_fmm_framework_*.rds (latest)
#' OUTPUT: submission/figures/Fig5_fmm_framework.pdf (+ mirrors)
#'
#' USAGE:  Rscript examples/publication/fig5_fmm_replot.R

setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
source("examples/publication/_pub_style.R")

cache_dir  <- "output/fmm_framework/results"
cache_path <- rev(sort(list.files(
  cache_dir, "^fig5_fmm_framework_.*\\.rds$", full.names = TRUE)))[1]
if (is.na(cache_path)) stop("No fig5_fmm_framework_*.rds in ", cache_dir)
cat(sprintf("Loading: %s\n", cache_path))
r <- readRDS(cache_path)

N_GRID    <- r$N_GRID
ETA_GRID  <- r$ETA_GRID
FDR_GRID  <- r$FDR_GRID
eta_hat   <- r$eta_hat
strata    <- r$strata_labels
power_marg  <- r$power_marg            # [N x FDR x sim]
power_strat <- r$power_strat           # [N x stratum x FDR x sim]
power_eta   <- r$power_eta             # [N x eta]

# Mean + SE over sims for plotting (SE = sd / sqrt(n) per cell)
marg_mean  <- apply(power_marg,  c(1, 2),    mean, na.rm = TRUE)
marg_se    <- apply(power_marg,  c(1, 2),
                    function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
strat_mean <- apply(power_strat, c(1, 2, 3), mean, na.rm = TRUE)
strat_se   <- apply(power_strat, c(1, 2, 3),
                    function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

# Panel B uses FDR = 0.05 by default
idx_fdr_pB <- which(FDR_GRID == 0.05)

# Drop strata with no rhythmic genes at any N (NaN columns)
strat_at_05    <- strat_mean[, , idx_fdr_pB]      # [N x stratum]
strat_at_05_se <- strat_se  [, , idx_fdr_pB]
strata_ok      <- apply(strat_at_05, 2, function(x) any(is.finite(x)))
strat_pB       <- strat_at_05   [, strata_ok, drop = FALSE]
strat_pB_se    <- strat_at_05_se[, strata_ok, drop = FALSE]
strata_lbl     <- strata[strata_ok]
# Rewrite any "(X, Inf]" stratum label as ">X" to match Fig 1 convention.
strata_lbl     <- sub("\\(([^,]+),\\s*Inf\\]$", ">\\1", strata_lbl)

# Panel B: down-select N to a representative subset (mirrors Fig 1's
# display_sizes pattern). 5 N values keep the legend compact and visually
# separate without flattening the trend.
disp_N    <- c(24L, 48L, 72L, 96L, 144L)
disp_idx  <- which(N_GRID %in% disp_N)

# Helper for vertical SE error bars (matches Fig 1/2 style)
add_se_bars <- function(x, y, se, col, bar_width = 0.3) {
  valid <- !is.na(y) & !is.na(se) & se > 0
  if (!any(valid)) return(invisible())
  arrows(x[valid], (y - se)[valid], x[valid], (y + se)[valid],
         angle = 90, code = 3, length = bar_width * 0.05,
         col = col, lwd = 1.2)
}

# Output paths
out_paths <- c(
  "output/fmm_framework/figures/Fig5_fmm_framework.pdf",
  "output/main_figures/Fig5_fmm_framework.pdf",
  "submission/figures/Fig5_fmm_framework.pdf"
)
for (p in out_paths) dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)

# Aesthetics — match Fig 1 conventions
thresh_cols <- c("darkgreen", "steelblue", "orange", "red")[seq_along(FDR_GRID)]
thresh_lbls <- paste0("FDR ", round(100 * FDR_GRID), "%")
N_cols      <- rainbow(length(N_GRID), s = 0.6, v = 0.8)
pal_eta     <- pub_palette_sequential(length(ETA_GRID))
idx_eta_anchor <- which.min(abs(ETA_GRID - eta_hat))

# ----- Draw -----
draw_fig5 <- function(out_pdf) {
  # Match Fig 1's font sizes + canvas dimensions so panel titles and axis
  # labels render at the same effective size after manuscript scaling.
  cairo_pdf(out_pdf, width = 11.5, height = 4.9)
  par(mfrow    = c(1, 3),
      mai      = c(1.60, 1.15, 0.85, 0.18),
      mgp      = c(3.6, 0.75, 0),
      oma      = c(0, 0, 2.4, 0),
      family   = "Helvetica", las = 1, tcl = -0.3,
      cex.axis = 1.35, cex.lab = 1.55, cex.main = 1.55, font.main = 2)

  # --- Panel A: marginal power vs N, lines per FDR + SE bars + boxed legend ---
  plot(NA, xlim = range(N_GRID), ylim = c(0, 1),
       xlab = expression(N ~ "(total samples)"),
       ylab = "Power",
       main = "Power vs Sample Size")
  panel_label("A")
  abline_80pct()
  for (a_i in seq_along(FDR_GRID)) {
    lines(N_GRID, marg_mean[, a_i], col = thresh_cols[a_i], lwd = 2.2)
    points(N_GRID, marg_mean[, a_i], col = thresh_cols[a_i], pch = 19, cex = 0.7)
    add_se_bars(N_GRID, marg_mean[, a_i], marg_se[, a_i], col = thresh_cols[a_i])
  }
  legend("bottomright",
         legend = thresh_lbls,
         col    = thresh_cols, lwd = 2.2, pch = 19,
         bty    = "o", bg = "white", cex = 0.95)

  # --- Panel B: power stratified by r-tilde at FDR = 0.05 ---
  # mgp[1] = 5.5 pushes xlab below the rotated strata tick labels (matches Fig 1).
  n_strata <- length(strata_lbl)
  par(mgp = c(5.5, 0.7, 0))
  plot(NA, xlim = c(0.5, n_strata + 0.5), ylim = c(0, 1),
       xlab = expression(tilde(r) == A / sigma ~ " stratum"),
       ylab = "Power (FDR = 0.05)",
       main = bquote(bold("Stratified Power by") ~ bold(tilde(r))),
       xaxt = "n")
  panel_label("B")
  abline_80pct()
  axis(1, at = seq_len(n_strata), labels = strata_lbl,
       las = 2, cex.axis = 0.95)
  for (j in disp_idx) {
    yj <- strat_pB[j, ]
    sj <- strat_pB_se[j, ]
    lines(seq_len(n_strata), yj, col = N_cols[j], lwd = 2.0)
    points(seq_len(n_strata), yj, col = N_cols[j], pch = 19, cex = 0.7)
    add_se_bars(seq_len(n_strata), yj, sj, col = N_cols[j])
  }
  legend("bottomright",
         legend = paste0("N = ", N_GRID[disp_idx]),
         col    = N_cols[disp_idx], lwd = 2.0, pch = 19,
         bty    = "o", bg = "white", cex = 0.95)
  par(mgp = c(3.4, 0.7, 0))

  # --- Panel C: eta sweep ---
  plot(NA, xlim = range(N_GRID), ylim = c(0, 1),
       xlab = expression(N ~ "(total samples)"),
       ylab = "Power (FDR = 0.05)",
       main = bquote(bold(omega[g]) ~ bold("~") ~ bold(Beta(1, eta))))
  panel_label("C")
  abline_80pct()
  for (k in seq_along(ETA_GRID)) {
    lwd_k <- if (k == idx_eta_anchor) 2.6 else 1.7
    lty_k <- if (ETA_GRID[k] == 0) 2 else 1
    lines(N_GRID, power_eta[, k], col = pal_eta[k], lwd = lwd_k, lty = lty_k)
    points(N_GRID, power_eta[, k], col = pal_eta[k], pch = 19, cex = 0.6)
  }
  legend("bottomright",
         legend = sapply(ETA_GRID, function(b)
           if (b == 0) expression(eta == 0 ~ "(cosinor)")
           else as.expression(bquote(eta == .(b)))),
         col   = pal_eta, lwd = 1.8, pch = NA,
         lty   = ifelse(ETA_GRID == 0, 2, 1),
         bty   = "o", bg = "white", cex = 0.95)

  mtext("Full FMM Framework Power Analysis (Baboon LUN, Active)",
        outer = TRUE, cex = 1.5, font = 2)

  dev.off()
  cat(sprintf("Saved: %s\n", out_pdf))
}

for (p in out_paths) draw_fig5(p)
