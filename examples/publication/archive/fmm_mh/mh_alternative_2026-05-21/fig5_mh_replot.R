#' Render Fig 5 MH from cache. Mirrors fig5_fmm_replot.R structure:
#'   Panel A -- marginal power vs N, one line per FDR
#'   Panel B -- power stratified by r-tilde at FDR = 0.05, one line per N
#'   Panel C -- alpha2 sweep (replaces eta sweep in FMM version)
#' Output: submission/figures/mh/Fig5_mh_framework.pdf (+ mirrors)

setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
source("examples/publication/_pub_style.R")

cache_dir  <- "output/mh_alternative/results"
cache_path <- rev(sort(list.files(cache_dir,
                                  "^fig5_mh_framework_.*\\.rds$",
                                  full.names = TRUE)))[1]
cat(sprintf("Loading: %s\n", cache_path))
r <- readRDS(cache_path)

N_GRID       <- r$N_GRID
ETA_A2_GRID  <- r$ETA_A2_GRID
FDR_GRID     <- r$FDR_GRID
eta_a2_anc   <- r$eta_a2_anchor
strata       <- r$strata_labels
power_marg   <- r$power_marg                # [N x FDR x sim]
power_strat  <- r$power_strat               # [N x stratum x FDR x sim]
power_a2     <- r$power_eta                 # [N x eta_a2]

marg_mean  <- apply(power_marg,  c(1, 2),    mean, na.rm = TRUE)
marg_se    <- apply(power_marg,  c(1, 2),
                    function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
strat_mean <- apply(power_strat, c(1, 2, 3), mean, na.rm = TRUE)
strat_se   <- apply(power_strat, c(1, 2, 3),
                    function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
idx_fdr_pB <- which(FDR_GRID == 0.05)
strat_at_05    <- strat_mean[, , idx_fdr_pB]
strat_at_05_se <- strat_se  [, , idx_fdr_pB]
strata_ok      <- apply(strat_at_05, 2, function(x) any(is.finite(x)))
strat_pB       <- strat_at_05   [, strata_ok, drop = FALSE]
strat_pB_se    <- strat_at_05_se[, strata_ok, drop = FALSE]
strata_lbl     <- strata[strata_ok]
strata_lbl     <- sub("\\(([^,]+),\\s*Inf\\]$", ">\\1", strata_lbl)

disp_N   <- c(24L, 48L, 72L, 96L, 144L)
disp_idx <- which(N_GRID %in% disp_N)

add_se_bars <- function(x, y, se, col, bar_width = 0.3) {
  valid <- !is.na(y) & !is.na(se) & se > 0
  if (!any(valid)) return(invisible())
  arrows(x[valid], (y - se)[valid], x[valid], (y + se)[valid],
         angle = 90, code = 3, length = bar_width * 0.05,
         col = col, lwd = 1.2)
}

out_paths <- c(
  "output/mh_alternative/figures/Fig5_mh_framework.pdf",
  "submission/figures/mh/Fig5_mh_framework.pdf"
)
for (p in out_paths) dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)

thresh_cols <- c("darkgreen", "steelblue", "orange", "red")[seq_along(FDR_GRID)]
thresh_lbls <- paste0("FDR ", round(100 * FDR_GRID), "%")
N_cols      <- rainbow(length(N_GRID), s = 0.6, v = 0.8)
pal_a2        <- pub_palette_sequential(length(ETA_A2_GRID))
idx_a2_anchor <- which.min(abs(ETA_A2_GRID - eta_a2_anc))

draw <- function(out_pdf) {
  cairo_pdf(out_pdf, width = 11.5, height = 4.9)
  par(mfrow = c(1, 3),
      mai = c(1.60, 1.15, 0.85, 0.18),
      mgp = c(3.6, 0.75, 0),
      oma = c(0, 0, 2.4, 0),
      family = "Helvetica", las = 1, tcl = -0.3,
      cex.axis = 1.35, cex.lab = 1.55, cex.main = 1.55, font.main = 2)

  # Panel A: marginal power vs N, lines per FDR
  plot(NA, xlim = range(N_GRID), ylim = c(0, 1),
       xlab = expression(N ~ "(total samples)"), ylab = "Power",
       main = "Power vs Sample Size")
  panel_label("A")
  abline_80pct()
  for (a_i in seq_along(FDR_GRID)) {
    lines(N_GRID, marg_mean[, a_i], col = thresh_cols[a_i], lwd = 2.2)
    points(N_GRID, marg_mean[, a_i], col = thresh_cols[a_i], pch = 19, cex = 0.7)
    add_se_bars(N_GRID, marg_mean[, a_i], marg_se[, a_i], col = thresh_cols[a_i])
  }
  legend("bottomright", legend = thresh_lbls,
         col = thresh_cols, lwd = 2.2, pch = 19,
         bty = "o", bg = "white", cex = 0.95)

  # Panel B: stratified power by r-tilde at FDR = 0.05
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
    yj <- strat_pB[j, ]; sj <- strat_pB_se[j, ]
    lines(seq_len(n_strata), yj, col = N_cols[j], lwd = 2.0)
    points(seq_len(n_strata), yj, col = N_cols[j], pch = 19, cex = 0.7)
    add_se_bars(seq_len(n_strata), yj, sj, col = N_cols[j])
  }
  legend("bottomright", legend = paste0("N = ", N_GRID[disp_idx]),
         col = N_cols[disp_idx], lwd = 2.0, pch = 19,
         bty = "o", bg = "white", cex = 0.95)
  par(mgp = c(3.6, 0.75, 0))

  # Panel C: alpha2 sweep
  plot(NA, xlim = range(N_GRID), ylim = c(0, 1),
       xlab = expression(N ~ "(total samples)"),
       ylab = "Power (FDR = 0.05)",
       main = bquote(bold(alpha[2*"g"] ~ "~ Beta(1, " ~ eta[alpha[2]] ~ ")")))
  panel_label("C")
  abline_80pct()
  for (k in seq_along(ETA_A2_GRID)) {
    lwd_k <- if (k == idx_a2_anchor) 2.6 else 1.7
    lines(N_GRID, power_a2[, k], col = pal_a2[k], lwd = lwd_k)
    points(N_GRID, power_a2[, k], col = pal_a2[k], pch = 19, cex = 0.6)
  }
  legend("bottomright",
         legend = sapply(ETA_A2_GRID, function(b)
           as.expression(bquote(eta[alpha[2]] == .(b)))),
         col = pal_a2, lwd = 1.8, pch = NA,
         bty = "o", bg = "white", cex = 0.85)

  mtext("Full MH Framework Power Analysis (Baboon LUN, Active)",
        outer = TRUE, cex = 1.5, font = 2)
  dev.off()
  cat(sprintf("Saved: %s\n", out_pdf))
}

for (p in out_paths) draw(p)
