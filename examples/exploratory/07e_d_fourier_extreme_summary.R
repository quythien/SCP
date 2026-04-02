#' =======================================================================
#' 07e_d_fourier_extreme_summary.R — Extreme Harmonics: Summary Figures
#' =======================================================================
#'
#' Loads RDS outputs from 07e_a/07e_b/07e_c and generates:
#'   - s4_heatmap.pdf          (3-dataset x 2-B heatmaps across full N x alpha2 grid)
#'   - s4_extreme_comparison.pdf  (original 3-panel slice at actual pilot N)
#'   - extreme_summary.txt
#'
#' Run AFTER 07e_a, 07e_b, 07e_c have all completed.

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

# Actual pilot N per dataset (n subjects per group in the real study)
pilot_n <- list(mouse = 24L, baboon = 12L, d1d2 = 45L)

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

a2_vals <- sort(unique(datasets[[1]]$lo$harmonic_grid$alpha2))

power_palette <- colorRampPalette(
  c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b")
)(100)

# =======================================================================
# FIGURE 1: Full N x alpha2 heatmaps (3 datasets x 2 B values)
# =======================================================================
cat("\n====================================================================\n")
cat("FIGURE 1: Full grid heatmaps\n")
cat("====================================================================\n\n")

fig_heat <- file.path(out_dir, "s4_heatmap.pdf")
n_ds <- length(datasets)
pdf(fig_heat, width = 5 * 2, height = 4.5 * n_ds)
par(mfrow = c(n_ds, 2), mar = c(4, 4.5, 3, 3.5))

for (nm in names(datasets)) {
  ds     <- datasets[[nm]]
  pn     <- pilot_n[[nm]]

  for (b_slot in c("lo", "hi")) {
    res <- ds[[b_slot]]
    B   <- if (b_slot == "lo") ds$B_low else ds$B_high
    if (is.null(res)) { plot.new(); next }

    # power_mean is [n_a2 x n_N]; transpose to [n_N x n_a2] for image (rows=N, cols=a2)
    pm  <- t(res$power_mean)   # [n_N x n_a2]
    Ns  <- res$sample_sizes
    a2s <- a2_vals

    image(x = a2s, y = Ns,
          z = t(pm),           # image: z[i,j] -> x[i], y[j]
          col  = power_palette,
          zlim = c(0, 1),
          xlab = expression(alpha[2] ~ "(2nd harmonic)"),
          ylab = "N per group",
          main = sprintf("%s  |  B = %d", ds$label, B),
          las  = 1)

    # Contour at 80%
    contour(x = a2s, y = Ns, z = t(pm),
            levels = 0.80, labels = "80%",
            lwd = 2, col = "white", add = TRUE, labcex = 0.8)

    # Mark actual pilot N
    abline(h = pn, lty = 2, col = "red", lwd = 1.5)
    mtext(sprintf("pilot n=%d", pn), side = 4, at = pn,
          col = "red", cex = 0.7, las = 1, line = 0.3)

    # Mark alpha2 = 0.5 practical limit
    abline(v = 0.5, lty = 3, col = "white", lwd = 1.2)

    # Color bar (manual)
    cat(sprintf("  %s B=%d: N grid %s | a2 grid %s\n",
                nm, B,
                paste(Ns,  collapse = ","),
                paste(a2s, collapse = ",")))
  }
}
dev.off()
cat(sprintf("Figure: %s\n", fig_heat))

# =======================================================================
# FIGURE 2: Original-style power vs alpha2 at actual pilot N
# (replaces old ref_n — now uses true pilot size per dataset)
# =======================================================================
cat("\n====================================================================\n")
cat("FIGURE 2: Power vs alpha2 at actual pilot N\n")
cat("====================================================================\n\n")

fig_s4 <- file.path(out_dir, "s4_extreme_comparison.pdf")
pdf(fig_s4, width = 5 * n_ds, height = 5)
par(mfrow = c(1, n_ds), mar = c(4.5, 4.5, 3.5, 1.5))

summary_lines <- c(
  "Fourier Extreme Waveform Robustness",
  sprintf("RUN_TAG: %s", RUN_TAG),
  "alpha3 = 0 throughout (isolates 2nd harmonic effect)",
  "ref_n = actual pilot size per dataset",
  ""
)

for (nm in names(datasets)) {
  ds <- datasets[[nm]]
  pn <- pilot_n[[nm]]
  if (is.null(ds$lo) || is.null(ds$hi)) { plot.new(); next }

  ref_idx_lo <- which.min(abs(ds$lo$sample_sizes - pn))
  ref_idx_hi <- which.min(abs(ds$hi$sample_sizes - pn))
  pow_lo <- ds$lo$power_mean[, ref_idx_lo]
  pow_hi <- ds$hi$power_mean[, ref_idx_hi]
  se_lo  <- ds$lo$power_se[,  ref_idx_lo]
  se_hi  <- ds$hi$power_se[,  ref_idx_hi]
  ylim   <- c(0, min(1, max(pow_lo + se_lo, pow_hi + se_hi, na.rm = TRUE) * 1.1))

  plot(a2_vals, pow_lo, type = "b", pch = 16, lwd = 2.5, col = ds$col,
       ylim = ylim,
       xlab = expression(alpha[2] ~ "(2nd harmonic relative amplitude)"),
       ylab = "Power (FDR 5%)", las = 1,
       main = sprintf("%s\npilot N = %d per group", ds$label, pn))
  polygon(c(a2_vals, rev(a2_vals)), c(pow_lo - se_lo, rev(pow_lo + se_lo)),
          col = adjustcolor(ds$col, 0.12), border = NA)
  lines(a2_vals, pow_hi, type = "b", pch = 17, lwd = 2.5, col = ds$col, lty = 2)
  polygon(c(a2_vals, rev(a2_vals)), c(pow_hi - se_hi, rev(pow_hi + se_hi)),
          col = adjustcolor(ds$col, 0.08), border = NA)
  abline(h = 0.80, lty = 3, col = "gray50")
  abline(v = 0.5,  lty = 2, col = "gray70")
  text(0.51, ylim[2] * 0.95, "bulk\nmax", col = "gray50", cex = 0.7, adj = 0)
  legend("topright",
         legend = sprintf("B=%d", c(ds$B_low, ds$B_high)),
         col = ds$col, lwd = 2.5, pch = c(16, 17), lty = c(1, 2),
         bty = "n", cex = 0.85)

  for (b_slot in c("lo", "hi")) {
    res <- ds[[b_slot]]
    B   <- if (b_slot == "lo") ds$B_low else ds$B_high
    if (!is.null(res)) {
      ref_idx <- which.min(abs(res$sample_sizes - pn))
      p0 <- res$power_mean[a2_vals == 0,       ref_idx]
      p1 <- res$power_mean[which.max(a2_vals), ref_idx]
      summary_lines <- c(summary_lines,
        sprintf("%-10s B=%2d (pilot_n=%d): power %.0f%% (a2=0) -> %.0f%% (a2=1.0)  delta=%.0f pp",
                toupper(nm), B, pn, 100 * p0, 100 * p1, 100 * (p0 - p1)))
    }
  }
  summary_lines <- c(summary_lines, "")
}
dev.off()
cat(sprintf("Figure: %s\n", fig_s4))

# =======================================================================
# FIGURE 3: Power vs alpha2 — all N values on one plot per dataset
# =======================================================================
cat("\n====================================================================\n")
cat("FIGURE 3: Power vs alpha2 across all N values\n")
cat("====================================================================\n\n")

fig_allN <- file.path(out_dir, "s4_allN_comparison.pdf")
pdf(fig_allN, width = 5 * n_ds, height = 5)
par(mfrow = c(1, n_ds), mar = c(4.5, 4.5, 3.5, 1.5))

for (nm in names(datasets)) {
  ds  <- datasets[[nm]]
  pn  <- pilot_n[[nm]]
  res <- ds$lo   # use sparse B only to avoid overplotting; add hi as dashed
  if (is.null(res)) { plot.new(); next }

  Ns       <- res$sample_sizes
  n_cols   <- colorRampPalette(c("gray80", ds$col))(length(Ns))
  ylim     <- c(0, 1)

  plot(NA, xlim = range(a2_vals), ylim = ylim,
       xlab = expression(alpha[2]),
       ylab = "Power (FDR 5%)", las = 1,
       main = sprintf("%s\nsolid=B=%d, dashed=B=%d",
                      ds$label, ds$B_low, ds$B_high))
  abline(h = 0.80, lty = 3, col = "gray70")
  abline(v = 0.50, lty = 3, col = "gray70")

  for (ni in seq_along(Ns)) {
    lines(a2_vals, res$power_mean[, ni],
          col = n_cols[ni], lwd = 2, lty = 1)
    if (!is.null(ds$hi)) {
      lines(a2_vals, ds$hi$power_mean[, ni],
            col = n_cols[ni], lwd = 2, lty = 2)
    }
    # mark pilot N
    if (Ns[ni] == pn) {
      points(a2_vals, res$power_mean[, ni], pch = 16, col = n_cols[ni], cex = 1.2)
    }
  }
  legend("topright",
         legend = paste0("N=", Ns),
         col = n_cols, lwd = 2, bty = "n", cex = 0.7)
}
dev.off()
cat(sprintf("Figure: %s\n", fig_allN))

# =======================================================================
# Summary text
# =======================================================================
summary_lines <- c(summary_lines,
  "Key interpretation:",
  "  - ref_n = actual pilot size (Mouse=24, Baboon=12, D1D2=45)",
  "  - alpha2=0.5 is the practical upper bound for bulk RNA-seq",
  "  - alpha2=0.75-1.0 approximates single-cell or spike-like waveforms",
  "  - Red dashed line in heatmap = actual pilot N",
  "  - White contour in heatmap = 80% power boundary"
)

writeLines(summary_lines, file.path(out_dir, "extreme_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
cat(sprintf("\nOutput: %s/\nDone.\n", out_dir))
