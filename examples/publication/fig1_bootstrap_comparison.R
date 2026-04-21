#' =======================================================================
#' fig1_bootstrap_comparison.R — Figure 1: Two-Stage vs Bootstrap Power
#' =======================================================================
#'
#' Generates a 3-panel figure comparing two-stage (plug-in) vs bootstrap
#' power curves for three pilot datasets. Each panel shows:
#'   - Solid line:  two-stage power (SE as T-bars)
#'   - Dashed line: bootstrap mean power
#'   - Shaded band: bootstrap 95% CI
#'   - 80% reference line
#'
#' Datasets:
#'   A. Baboon LUN vs CER  (n=12, active, r~1.72, DR~41%)
#'   B. Mouse D1 vs D2     (n=45, active, r~0.65, DR~20%)
#'   C. Human PFC aging    (n=60, passive, r~0.56, DR~14%)
#'
#' Output: paper/PowerSim/figures/bootstrap_summary.pdf
#'
#' Usage:
#'   Rscript examples/publication/fig1_bootstrap_comparison.R

POWERSIM_ROOT <- {
  env <- Sys.getenv("POWERSIM_ROOT", unset = "")
  if (nzchar(env)) env else
    "/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim"
}
setwd(POWERSIM_ROOT)

# -----------------------------------------------------------------------
# Dataset registry
# -----------------------------------------------------------------------
RUN_TAG_08 <- Sys.getenv("RUN_TAG", "20260330_1126")

datasets <- list(
  list(
    rds     = file.path("output", "08_two_stage_vs_bootstrap_realdata_20260330_1126",
                        "s1_baboon_comparison.rds"),
    label   = "Baboon LUN vs CER\n(n=12 pilot, active, r~1.72, DR~41%)",
    col     = "darkorange",
    panel   = "A"
  ),
  list(
    rds     = file.path("output", paste0("08_two_stage_vs_bootstrap_", RUN_TAG_08),
                        "s2_d1d2_comparison.rds"),
    label   = "Mouse D1 vs D2\n(n=45 pilot, active, r~0.65, DR~20%)",
    col     = "forestgreen",
    panel   = "B"
  ),
  list(
    rds     = file.path("output", paste0("08_two_stage_vs_bootstrap_", RUN_TAG_08),
                        "s3_seney_comparison.rds"),
    label   = "Human PFC aging\n(n=60 pilot, passive, r~0.56, DR~14%)",
    col     = "firebrick",
    panel   = "C"
  )
)

# -----------------------------------------------------------------------
# Load RDS files
# -----------------------------------------------------------------------
panels <- lapply(datasets, function(d) {
  if (!file.exists(d$rds)) {
    warning(sprintf("Missing: %s", d$rds))
    return(NULL)
  }
  dat <- readRDS(d$rds)
  cmp <- dat$comparison$comparison
  list(
    n       = cmp$n,
    ts      = cmp$two_stage_power,
    ts_se   = cmp$two_stage_se,
    bm      = cmp$boot_power_mean,
    bm_lo   = cmp$boot_ci_lo,
    bm_hi   = cmp$boot_ci_hi,
    n80_ts  = dat$comparison$n80_two_stage,
    n80_bm  = dat$comparison$n80_boot_median
  )
})

# -----------------------------------------------------------------------
# Build 3-panel figure
# -----------------------------------------------------------------------
out_file <- "paper/PowerSim/figures/bootstrap_summary.pdf"
pdf(out_file, width = 15, height = 5)
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3.5, 1.5), oma = c(0, 0, 0, 0))

for (i in seq_along(datasets)) {
  d <- datasets[[i]]
  p <- panels[[i]]

  if (is.null(p)) {
    plot.new()
    title(main = d$label, sub = "(data unavailable)")
    mtext(d$panel, side = 3, adj = 0, font = 2, line = 0.5)
    next
  }

  N    <- p$n
  ts   <- p$ts  * 100
  ts_se <- (p$ts_se %||% rep(0, length(ts))) * 100
  bm   <- p$bm  * 100
  lo   <- p$bm_lo * 100
  hi   <- p$bm_hi * 100

  xlim <- c(0, max(N) * 1.08)
  ylim <- c(0, 100)

  plot(NA, xlim = xlim, ylim = ylim, las = 1, xaxt = "n",
       xlab = "N per group", ylab = "DR power (FDR 5%)",
       main = d$label)
  axis(1, at = N)
  abline(h = c(20, 40, 60, 80), lty = 3, col = "gray85")

  # Bootstrap CI band
  polygon(c(N, rev(N)), c(lo, rev(hi)),
          col = adjustcolor(d$col, 0.15), border = NA)

  # Bootstrap mean (dashed)
  lines(N, bm, col = d$col, lwd = 2, lty = 2)
  points(N, bm, col = d$col, pch = 1, cex = 0.85)

  # Two-stage (solid) with SE T-bars
  lines(N, ts, col = d$col, lwd = 2.5, lty = 1)
  points(N, ts, col = d$col, pch = 19, cex = 0.85)
  if (any(ts_se > 0, na.rm = TRUE)) {
    arrows(N, ts - ts_se, N, ts + ts_se,
           angle = 90, code = 3, length = 0.04, col = d$col, lwd = 1.2)
  }

  # n80 markers
  if (!is.na(p$n80_ts)) {
    abline(v = p$n80_ts, lty = 2, col = adjustcolor(d$col, 0.5))
    text(p$n80_ts + max(N) * 0.02, 10,
         sprintf("n80~%d\n(TS)", p$n80_ts),
         cex = 0.68, col = d$col, adj = 0, font = 2)
  }
  if (!is.na(p$n80_bm) && (is.na(p$n80_ts) || p$n80_bm != p$n80_ts)) {
    abline(v = p$n80_bm, lty = 3, col = adjustcolor(d$col, 0.5))
    text(p$n80_bm + max(N) * 0.02, 22,
         sprintf("n80~%d\n(BS)", p$n80_bm),
         cex = 0.68, col = d$col, adj = 0, font = 2)
  }

  legend("bottomright",
         legend = c("Two-stage (plug-in)", "Bootstrap", "Bootstrap 95% CI"),
         col    = c(d$col, d$col, adjustcolor(d$col, 0.4)),
         lwd    = c(2.5, 2, 6),
         lty    = c(1, 2, 1),
         pch    = c(19, 1, NA),
         pt.cex = 0.85,
         bty = "n", cex = 0.75)

  mtext(d$panel, side = 3, adj = 0, font = 2, line = 0.5)
}

dev.off()
cat(sprintf("Figure saved: %s\n", out_file))

# Print n80 summary
cat("\n--- n80 summary ---\n")
for (i in seq_along(datasets)) {
  d <- datasets[[i]]
  p <- panels[[i]]
  if (is.null(p)) next
  cat(sprintf("  %s: TS n80=%s  BS n80=%s\n",
              sub("\n.*", "", d$label),
              ifelse(is.na(p$n80_ts), ">max(N)", p$n80_ts),
              ifelse(is.na(p$n80_bm), ">max(N)", p$n80_bm)))
}
