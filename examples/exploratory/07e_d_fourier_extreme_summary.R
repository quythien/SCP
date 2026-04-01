#' =======================================================================
#' 07e_d_fourier_extreme_summary.R — Extreme Harmonics: Summary Figure
#' =======================================================================
#'
#' Loads RDS outputs from 07e_a/07e_b/07e_c and generates:
#'   - s4_extreme_comparison.pdf  (3-panel cross-dataset figure)
#'   - extreme_summary.txt
#'
#' Run AFTER 07e_a, 07e_b, 07e_c have all completed:
#'   POWERSIM_ROOT=$ROOT RUN_TAG=$TAG Rscript examples/exploratory/07e_d_fourier_extreme_summary.R

POWERSIM_ROOT <- {
  env <- Sys.getenv("POWERSIM_ROOT", unset = "")
  if (nzchar(env)) env else
    "/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim"
}
setwd(POWERSIM_ROOT)

RUN_TAG <- Sys.getenv("RUN_TAG", format(Sys.Date(), "%Y%m%d"))
out_dir <- file.path("output", paste0("07_fourier_extreme_", RUN_TAG))

if (!dir.exists(out_dir)) {
  stop(sprintf("Output dir not found: %s\nDid you set RUN_TAG correctly?", out_dir))
}

source_dir <- file.path(getwd(), "code")
old_wd <- setwd(source_dir); source("setup.R"); setwd(old_wd)

cat(sprintf("Loading results from: %s/\n\n", out_dir))

rds_map <- list(
  mouse  = file.path(out_dir, "s1_mouse_extreme.rds"),
  baboon = file.path(out_dir, "s2_baboon_extreme.rds"),
  d1d2   = file.path(out_dir, "s3_d1d2_extreme.rds")
)

datasets <- list()
for (nm in names(rds_map)) {
  f <- rds_map[[nm]]
  if (!file.exists(f)) {
    warning(sprintf("Missing RDS: %s — skipping %s", f, nm)); next
  }
  datasets[[nm]] <- readRDS(f)
  cat(sprintf("  Loaded: %s\n", nm))
}

if (length(datasets) == 0) stop("No datasets loaded.")

# Reconstruct HARM_GRID alpha2 values from first available result
first_res <- datasets[[1]]$lo
a2_vals   <- sort(unique(first_res$harmonic_grid$alpha2))

# =======================================================================
# SECTION 4: CROSS-DATASET COMPARISON FIGURE
# =======================================================================
cat("\n====================================================================\n")
cat("SECTION 4: Cross-dataset extreme harmonic comparison\n")
cat("====================================================================\n\n")

fig_s4 <- file.path(out_dir, "s4_extreme_comparison.pdf")
pdf(fig_s4, width = 15, height = 5)
par(mfrow = c(1, length(datasets)), mar = c(4.5, 4.5, 3.5, 1.5))

summary_lines <- c(
  "Fourier Extreme Waveform Robustness",
  sprintf("RUN_TAG: %s", RUN_TAG),
  "alpha3 = 0 throughout (isolates 2nd harmonic effect)",
  ""
)

for (nm in names(datasets)) {
  ds <- datasets[[nm]]
  if (is.null(ds$lo) || is.null(ds$hi)) { plot.new(); next }

  ref_idx_lo <- which.min(abs(ds$lo$sample_sizes - ds$ref_n))
  ref_idx_hi <- which.min(abs(ds$hi$sample_sizes - ds$ref_n))
  pow_lo <- ds$lo$power_mean[, ref_idx_lo]
  pow_hi <- ds$hi$power_mean[, ref_idx_hi]
  se_lo  <- ds$lo$power_se[,  ref_idx_lo]
  se_hi  <- ds$hi$power_se[,  ref_idx_hi]
  ylim   <- c(0, min(1, max(pow_lo + se_lo, pow_hi + se_hi, na.rm = TRUE) * 1.1))

  plot(a2_vals, pow_lo, type = "b", pch = 16, lwd = 2.5, col = ds$col,
       ylim = ylim, xlab = expression(alpha[2]),
       ylab = "Power (FDR 5%)", las = 1,
       main = sprintf("%s\nN=%d", ds$label, ds$ref_n))
  polygon(c(a2_vals, rev(a2_vals)), c(pow_lo - se_lo, rev(pow_lo + se_lo)),
          col = adjustcolor(ds$col, 0.12), border = NA)
  lines(a2_vals, pow_hi, type = "b", pch = 17, lwd = 2.5, col = ds$col, lty = 2)
  polygon(c(a2_vals, rev(a2_vals)), c(pow_hi - se_hi, rev(pow_hi + se_hi)),
          col = adjustcolor(ds$col, 0.08), border = NA)
  abline(h = 0.80, lty = 3, col = "gray50")
  abline(v = 0.5,  lty = 2, col = "gray70")
  legend("topright", legend = sprintf("B=%d", c(ds$B_low, ds$B_high)),
         col = ds$col, lwd = 2.5, pch = c(16, 17), lty = c(1, 2), bty = "n", cex = 0.85)

  for (res_nm in c("lo", "hi")) {
    res <- ds[[res_nm]]
    B   <- if (res_nm == "lo") ds$B_low else ds$B_high
    if (!is.null(res)) {
      ref_idx <- which.min(abs(res$sample_sizes - ds$ref_n))
      p0 <- res$power_mean[a2_vals == 0,        ref_idx]
      p1 <- res$power_mean[which.max(a2_vals),  ref_idx]
      summary_lines <- c(summary_lines,
        sprintf("%-10s B=%2d (ref_n=%d): power %.0f%% (a2=0) -> %.0f%% (a2=1.0)  delta=%.0f pp",
                toupper(nm), B, ds$ref_n, 100 * p0, 100 * p1, 100 * (p0 - p1)))
    }
  }
  summary_lines <- c(summary_lines, "")
}
dev.off()
cat(sprintf("Figure: %s\n", fig_s4))

summary_lines <- c(summary_lines,
  "Key interpretation:",
  "  - alpha2=0.5 is the typical upper bound for bulk RNA-seq circadian data",
  "  - alpha2=0.75-1.0 approximates single-cell or spike-like waveforms",
  "  - Small delta (<10 pp) across a2=0 to 1.0 => cosinor-DCP is robust",
  "  - Growing B advantage at extreme alpha2 => more time points help for non-sinusoidal data"
)

writeLines(summary_lines, file.path(out_dir, "extreme_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
cat(sprintf("\nOutput: %s/\nDone.\n", out_dir))
