#' =======================================================================
#' fig3_multi_dataset_dr_power.R — Figure 3: DR Power Across Datasets
#' =======================================================================
#'
#' Generates a 4-panel marginal DR power vs N figure combining all four
#' pilot datasets. Each panel plots mean power ± SE at FDR 5%, with the
#' 80% reference line.  Panels are ordered by effect size (r) descending
#' so the figure communicates the r → n80 relationship at a glance.
#'
#' Datasets (in panel order):
#'   1. Mouse LIV vs CER  (r~2.9,  active, GSE54651)     — 03d run
#'   2. Baboon LUN vs CER (r~1.72, active)                — 03b run
#'   3. Mouse D1 vs D2    (r~0.65, active)                — 03c run
#'   4. Human PFC aging   (r~0.56, passive)               — key run
#'
#' Output: paper/PowerSim/figures/multi_dataset_dr_power.pdf
#'
#' Usage:
#'   Rscript examples/publication/fig3_multi_dataset_dr_power.R

POWERSIM_ROOT <- {
  env <- Sys.getenv("POWERSIM_ROOT", unset = "")
  if (nzchar(env)) env else
    "/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim"
}
setwd(POWERSIM_ROOT)
source_dir <- file.path(getwd(), "code")
old_wd <- setwd(source_dir); source("setup.R"); setwd(old_wd)

# -----------------------------------------------------------------------
# Dataset registry — each entry: rds path, label, color, r, prop_DR
# -----------------------------------------------------------------------
datasets <- list(
  list(
    file    = "output/03d_power_core_mouse_gse_20260402_0023/dr_power_raw_pvalues.rds",
    label   = "Mouse LIV vs CER\n(GSE54651, r~2.9, DR~30%)",
    short   = "Mouse LIV vs CER",
    col     = "steelblue",
    r_med   = 2.9,
    prop_dr = 0.30
  ),
  list(
    file    = "output/03b_power_core_active_20260331_1159/dr_power_raw_pvalues.rds",
    label   = "Baboon LUN vs CER\n(r~1.72, DR~41%)",
    short   = "Baboon LUN vs CER",
    col     = "darkorange",
    r_med   = 1.72,
    prop_dr = 0.41
  ),
  list(
    file    = "output/03c_power_core_mouse_20260331_1215/dr_power_raw_pvalues.rds",
    label   = "Mouse D1 vs D2\n(r~0.65, DR~20%)",
    short   = "Mouse D1 vs D2",
    col     = "forestgreen",
    r_med   = 0.65,
    prop_dr = 0.20
  ),
  list(
    file    = "output/run_20260330_1126/dr_power_raw_pvalues.rds",
    label   = "Human PFC aging\n(r~0.56, DR~14%, passive)",
    short   = "Human PFC aging",
    col     = "firebrick",
    r_med   = 0.56,
    prop_dr = 0.14
  )
)

# -----------------------------------------------------------------------
# Extract marginal power mean + SE at FDR 5% directly from raw pvalues
# -----------------------------------------------------------------------
.extract_marginal <- function(rds_file) {
  nm  <- load(rds_file)
  obj <- get(nm)
  if (exists("dp_power_raw")) obj <- dp_power_raw   # safety alias

  n_sizes <- length(obj$sample_sizes)
  nsims   <- obj$nsims
  FDR_THR <- 0.05

  # Detect nested ground truth format
  nested_gt <- is.list(obj$is_target_list[[1]]) && length(obj$is_target_list) == n_sizes

  power_sim <- matrix(NA_real_, nrow = n_sizes, ncol = nsims)
  for (j in seq_len(n_sizes)) {
    for (s in seq_len(nsims)) {
      pvals <- obj$pvalues[j, , s]
      if (nested_gt) {
        is_target <- obj$is_target_list[[j]][[s]]
      } else {
        is_target <- obj$is_target_list[[s]]
      }
      qvals <- rep(1, length(pvals))
      tested <- !is.na(pvals) & pvals < 1
      if (sum(tested) > 0) qvals[tested] <- p.adjust(pvals[tested], method = "BH")
      n_tgt <- sum(is_target, na.rm = TRUE)
      power_sim[j, s] <- if (n_tgt > 0) sum(qvals <= FDR_THR & is_target, na.rm = TRUE) / n_tgt else NA
    }
  }

  list(
    N    = obj$sample_sizes,
    mean = rowMeans(power_sim, na.rm = TRUE),
    se   = apply(power_sim, 1, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
  )
}

cat("Extracting marginal power for each dataset...\n")
panels <- lapply(datasets, function(d) {
  cat(sprintf("  %s\n", d$short))
  tryCatch(
    .extract_marginal(d$file),
    error = function(e) { warning(sprintf("Failed %s: %s", d$short, e$message)); NULL }
  )
})

# -----------------------------------------------------------------------
# Build 4-panel figure
# -----------------------------------------------------------------------
out_file <- "paper/PowerSim/figures/multi_dataset_dr_power.pdf"
pdf(out_file, width = 14, height = 4)
par(mfrow = c(1, 4), mar = c(4.5, 4.5, 3.5, 1.5), oma = c(0, 0, 2, 0))
panel_labels <- c("A", "B", "C", "D")

for (i in seq_along(datasets)) {
  d <- datasets[[i]]
  p <- panels[[i]]
  if (is.null(p)) {
    plot.new(); title(main = d$short, sub = "(data unavailable)")
    next
  }

  N    <- p$N
  pm   <- p$mean * 100
  se   <- p$se   * 100

  xlim <- c(0, max(N) * 1.08)
  ylim <- c(0, 100)

  plot(N, pm, type = "b", pch = 19, lwd = 2, col = d$col,
       xlim = xlim, ylim = ylim, las = 1,
       xlab = "N per group", ylab = "DR power (FDR 5%)",
       main = d$label, xaxt = "n")
  axis(1, at = N)
  abline(h = 80, lty = 2, col = "gray40")
  abline(h = c(20, 40, 60), lty = 3, col = "gray85")

  # SE bars
  arrows(N, pm - se, N, pm + se,
         angle = 90, code = 3, length = 0.05, col = d$col, lwd = 1.2)

  # Annotate power values
  text(N, pm + 5, sprintf("%.0f%%", pm), cex = 0.65, col = d$col)

  # n80 marker: interpolate
  n80 <- NA
  if (any(pm >= 80, na.rm = TRUE)) {
    idx <- which(pm >= 80)[1]
    if (idx > 1) {
      frac <- (80 - pm[idx - 1]) / (pm[idx] - pm[idx - 1])
      n80  <- N[idx - 1] + frac * (N[idx] - N[idx - 1])
    } else {
      n80 <- N[idx]
    }
    abline(v = n80, lty = 2, col = adjustcolor(d$col, 0.6))
    mtext(sprintf("n80~%d", round(n80)), side = 1, at = n80, line = 2.2,
          cex = 0.65, col = d$col)
  }

  mtext(panel_labels[i], side = 3, adj = 0, font = 2, line = 0.5)
}

mtext("Differential Rhythmicity Power by Dataset (FDR 5%, mean ± SE across 50 simulations)",
      outer = TRUE, font = 2, cex = 0.95)
dev.off()
cat(sprintf("\nFigure saved: %s\n", out_file))
