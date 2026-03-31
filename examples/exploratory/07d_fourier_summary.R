#' =======================================================================
#' 07d_fourier_summary.R — Fourier Robustness: Cross-Dataset Summary
#' =======================================================================
#' Loads RDS outputs from 07a/07b/07c and generates the Section 4
#' cross-dataset B-protection summary figure + fourier_summary.txt.
#'
#' Run AFTER 07a, 07b, 07c have all completed:
#'   RUN_TAG=20260401 Rscript examples/exploratory/07d_fourier_summary.R
#'
#' Uses same RUN_TAG as the per-dataset scripts (default: today's date).

POWERSIM_ROOT <- {
  env <- Sys.getenv("POWERSIM_ROOT", unset = "")
  if (nzchar(env)) env else
    "/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim"
}
setwd(POWERSIM_ROOT)

RUN_TAG <- Sys.getenv("RUN_TAG", format(Sys.Date(), "%Y%m%d"))
out_dir <- file.path("output", paste0("07_fourier_robustness_", RUN_TAG))

if (!dir.exists(out_dir)) {
  stop(sprintf("Output dir not found: %s\nDid you set RUN_TAG correctly?", out_dir))
}

cat(sprintf("Loading results from: %s/\n\n", out_dir))

# -----------------------------------------------------------------------
# Load per-dataset RDS files
# -----------------------------------------------------------------------
rds_files <- list(
  mouse  = file.path(out_dir, "s1_mouse_fourier.rds"),
  baboon = file.path(out_dir, "s2_baboon_fourier.rds"),
  d1d2   = file.path(out_dir, "s3_d1d2_fourier.rds")
)

sec4_results <- list()
for (nm in names(rds_files)) {
  f <- rds_files[[nm]]
  if (!file.exists(f)) {
    warning(sprintf("Missing RDS: %s — skipping %s", f, nm))
    next
  }
  dat <- readRDS(f)
  # Reconstruct the sec4 entry expected by the plotting code
  # Each RDS contains: B_low/B_high RDS data + label/col/B_low/B_high/ref_n
  b_names <- names(dat)[!names(dat) %in% c("label", "col", "B_low", "B_high", "ref_n")]
  sec4_results[[nm]] <- list(
    label  = dat$label,
    col    = dat$col,
    B_low  = dat$B_low,
    B_high = dat$B_high,
    ref_n  = dat$ref_n,
    res_lo = dat[[b_names[1]]],
    res_hi = dat[[b_names[2]]]
  )
  cat(sprintf("  Loaded: %s (%s)\n", nm, dat$label))
}

if (length(sec4_results) == 0) stop("No datasets loaded — cannot produce summary.")

# -----------------------------------------------------------------------
# Section 4: Cross-dataset B protection figure
# -----------------------------------------------------------------------
cat("\n====================================================================\n")
cat("SECTION 4: Cross-dataset B protection summary\n")
cat("====================================================================\n\n")

fig_s4 <- file.path(out_dir, "s4_b_protection.pdf")
pdf(fig_s4, width = 5 * length(sec4_results), height = 5)
par(mfrow = c(1, length(sec4_results)), mar = c(4.5, 4.5, 3.5, 1.5))

SMOKE <- nzchar(Sys.getenv("POWERSIM_SMOKE"))
summary_lines <- c(
  "Fourier Robustness — B Protection Summary",
  sprintf("RUN_TAG: %s", RUN_TAG),
  ""
)

for (nm in names(sec4_results)) {
  sc     <- sec4_results[[nm]]
  res_lo <- sc$res_lo
  res_hi <- sc$res_hi
  ref_n  <- sc$ref_n

  if (is.null(res_lo) || is.null(res_hi)) {
    plot.new(); title(sprintf("%s\nresults unavailable", sc$label))
    next
  }

  hg      <- res_lo$harmonic_grid
  a3_0    <- hg$alpha3 == 0
  a2_vals <- sort(unique(hg$alpha2[a3_0]))

  ref_idx_lo <- which.min(abs(res_lo$sample_sizes - ref_n))
  ref_idx_hi <- which.min(abs(res_hi$sample_sizes - ref_n))

  pow_lo <- res_lo$power_mean[a3_0, ref_idx_lo]
  pow_hi <- res_hi$power_mean[a3_0, ref_idx_hi]
  se_lo  <- res_lo$power_se[a3_0,  ref_idx_lo]
  se_hi  <- res_hi$power_se[a3_0,  ref_idx_hi]

  pow_lo_0   <- pow_lo[a2_vals == 0]
  pow_hi_0   <- pow_hi[a2_vals == 0]
  pow_lo_max <- pow_lo[which.max(a2_vals)]
  pow_hi_max <- pow_hi[which.max(a2_vals)]
  delta_lo <- if (length(pow_lo_0) && length(pow_lo_max)) pow_lo_0 - pow_lo_max else NA
  delta_hi <- if (length(pow_hi_0) && length(pow_hi_max)) pow_hi_0 - pow_hi_max else NA

  summary_lines <- c(summary_lines,
    sprintf("%s (ref_n=%d):", sc$label, ref_n),
    sprintf("  B=%d: power %.0f%% (alpha2=0) -> %.0f%% (alpha2=%.2g)   delta=%.0f pp",
            sc$B_low,  100*pow_lo_0, 100*pow_lo_max, max(a2_vals), 100*delta_lo),
    sprintf("  B=%d: power %.0f%% (alpha2=0) -> %.0f%% (alpha2=%.2g)   delta=%.0f pp",
            sc$B_high, 100*pow_hi_0, 100*pow_hi_max, max(a2_vals), 100*delta_hi),
    sprintf("  B advantage: B=%d loses %.0f pp less than B=%d",
            sc$B_high, 100*(delta_lo - delta_hi), sc$B_low),
    ""
  )

  plot(a2_vals, pow_lo * 100, type = "b", pch = 16, lwd = 2.5, col = sc$col,
       ylim = c(0, 100), las = 1,
       xlab = expression(alpha[2] ~ "(2nd harmonic, a3=0)"),
       ylab = "Power (%, FDR 5%)",
       main = sprintf("%s\nN=%d: does B protect?", sc$label, ref_n))
  arrows(a2_vals, (pow_lo - se_lo) * 100, a2_vals, (pow_lo + se_lo) * 100,
         length = 0.05, angle = 90, code = 3, col = sc$col, lwd = 1.2)
  lines(a2_vals, pow_hi * 100, type = "b", pch = 17, lwd = 2.5, col = sc$col, lty = 2)
  arrows(a2_vals, (pow_hi - se_hi) * 100, a2_vals, (pow_hi + se_hi) * 100,
         length = 0.05, angle = 90, code = 3, col = sc$col, lwd = 1.2)
  abline(h = 80, lty = 3, col = "gray50")
  text(min(a2_vals), 82, "80%", col = "gray50", cex = 0.8, adj = 0)
  legend("bottomleft", legend = sprintf("B=%d", c(sc$B_low, sc$B_high)),
         col = sc$col, lwd = 2.5, pch = c(16, 17), lty = c(1, 2), bty = "n", cex = 0.85)
}

dev.off()
cat(sprintf("Figure: %s\n", fig_s4))

writeLines(summary_lines, file.path(out_dir, "fourier_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
cat(sprintf("\nOutput: %s/\nDone.\n", out_dir))
