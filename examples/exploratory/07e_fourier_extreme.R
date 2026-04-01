#' =======================================================================
#' 07e_fourier_extreme.R — Extreme Waveform Deviation Robustness
#' =======================================================================
#'
#' PURPOSE
#'   Extends 07_fourier_robustness.R to more extreme non-sinusoidal shapes.
#'   The standard analysis (07a/b/c) tests α₂ up to 0.5 — moderate deviation.
#'   This script pushes to α₂ = 1.0 (2nd harmonic as strong as fundamental),
#'   approaching square-wave and spike-like waveforms.
#'
#'   Motivating question: at what α₂ threshold does cosinor-based DCP start
#'   losing substantial power, and does higher B help at extreme deviations?
#'
#' HARMONIC SHAPES TESTED
#'   α₂ = 0.00 — pure sinusoid (cosinor assumption satisfied)
#'   α₂ = 0.25 — mild deviation (asymmetric peak)
#'   α₂ = 0.50 — moderate deviation (bimodal tendency)
#'   α₂ = 0.75 — strong deviation (two-peak waveform)
#'   α₂ = 1.00 — extreme (2nd harmonic = fundamental; approximates square wave)
#'   α₃ = 0 throughout (isolates 2nd harmonic effect cleanly)
#'
#'   For context: bulk RNA-seq circadian data typically shows α₂ ≤ 0.3.
#'   Single-cell or spike-like expression patterns may reach α₂ ~ 0.5–0.75.
#'   α₂ = 1.0 is near the theoretical worst case for a 24h oscillation.
#'
#' DATASETS
#'   All 3 active datasets — so we can see if signal strength (r) modifies
#'   the sensitivity to waveform misspecification.
#'
#' OUTPUTS
#'   output/07_fourier_extreme_<RUN_TAG>/
#'     s1_mouse_extreme.rds / .pdf
#'     s2_baboon_extreme.rds / .pdf
#'     s3_d1d2_extreme.rds / .pdf
#'     extreme_summary.txt
#'     s4_extreme_comparison.pdf  — all 3 datasets on one plot
#'
#' USAGE
#'   POWERSIM_ROOT=$ROOT Rscript examples/exploratory/07e_fourier_extreme.R
#'   POWERSIM_SMOKE=1   Rscript examples/exploratory/07e_fourier_extreme.R

# -----------------------------------------------------------------------
# 0. Path + settings
# -----------------------------------------------------------------------
POWERSIM_ROOT <- {
  env <- Sys.getenv("POWERSIM_ROOT", unset = "")
  if (nzchar(env)) env else
    "/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim"
}
setwd(POWERSIM_ROOT)

SMOKE       <- nzchar(Sys.getenv("POWERSIM_SMOKE"))
NGENES      <- if (SMOKE) 300L else 3000L
NSIMS       <- if (SMOKE) 5L  else 20L
RHYTHM_PVAL <- 0.05

# Extended harmonic grid: α₂ up to 1.0, α₃ = 0 only (isolate 2nd harmonic)
HARM_GRID <- if (SMOKE) {
  expand.grid(alpha2 = c(0, 0.5, 1.0), alpha3 = 0)
} else {
  expand.grid(alpha2 = c(0, 0.25, 0.5, 0.75, 1.0), alpha3 = 0)
}

# N grids — same as 07a/b/c
N_GRID_MOUSE  <- if (SMOKE) c(16L, 24L, 32L) else c(8L, 16L, 24L, 48L, 72L)
N_GRID_BABOON <- if (SMOKE) c(24L, 36L, 48L) else c(12L, 24L, 36L, 48L, 72L, 96L)
N_GRID_D1D2   <- if (SMOKE) c(24L, 48L, 72L) else c(24L, 48L, 72L, 96L, 120L)

# B values — compare low vs high B at extreme deviations
B_VALS_MOUSE  <- c(4L, 8L)
B_VALS_BABOON <- c(4L, 12L)
B_VALS_D1D2   <- c(3L, 6L)

REF_N_MOUSE  <- 24L
REF_N_BABOON <- 48L
REF_N_D1D2   <- 72L

DATA_MOUSE      <- "data/mice_GSE54651_CPM.RData"
DATA_BABOON     <- "data/CAMO_PRC_hmb.RData"
DATA_D1D2_PHENO <- "data/mouse_clinicalinfo_03082021_rmOutliers.csv"
DATA_D1D2_EXPR  <- "data/mouse_D1D2_logCPMfiltered_counts.csv"

RUN_TAG <- Sys.getenv("RUN_TAG", format(Sys.Date(), "%Y%m%d"))
out_dir <- file.path("output", paste0("07_fourier_extreme_", RUN_TAG))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Mode: %s | NGENES=%d | NSIMS=%d | alpha2_max=1.0 | harm_combos=%d\n",
            if (SMOKE) "SMOKE" else "PRODUCTION", NGENES, NSIMS, nrow(HARM_GRID)))
cat(sprintf("Output -> %s/\n\n", out_dir))

source_dir <- file.path(getwd(), "code")
old_wd <- setwd(source_dir); source("setup.R"); setwd(old_wd)

opts_analysis <- CircadianAnalysisOptions(
  alpha = 0.05, p.adjust.method = "BH",
  fdr_thresholds = c(0.05), reference_n = 48L
)

# -----------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------
.run_fourier_one_B <- function(bio_opts, N_grid, B, design_type = "active", label = "") {
  cts_B <- seq(0, 24, length.out = B + 1)[seq_len(B)]
  design_opts <- CircadianDesignOptions(
    sample_sizes = N_grid, nsims = NSIMS, design = design_type,
    cts = cts_B, test_types = "DR"
  )
  cat(sprintf("    B=%d: %s\n", B, label))
  tryCatch(
    runFourierDeviationPower(bio_opts, design_opts, opts_analysis,
                             harmonic_grid = HARM_GRID, test_type = "DR",
                             verbose = FALSE),
    error = function(e) { warning(sprintf("B=%d failed: %s", B, e$message)); NULL }
  )
}

.plot_extreme <- function(res_low, res_high, B_low, B_high,
                           dataset_label, ref_n, col, output_file) {
  pdf(output_file, width = 7, height = 5)
  par(mar = c(4.5, 4.5, 3.5, 1.5))

  a2_vals <- sort(unique(HARM_GRID$alpha2))

  if (!is.null(res_low) && !is.null(res_high)) {
    ref_idx_lo <- which.min(abs(res_low$sample_sizes  - ref_n))
    ref_idx_hi <- which.min(abs(res_high$sample_sizes - ref_n))

    # all rows are alpha3=0 since HARM_GRID only has alpha3=0
    pow_lo <- res_low$power_mean[,  ref_idx_lo]
    pow_hi <- res_high$power_mean[, ref_idx_hi]
    se_lo  <- res_low$power_se[,   ref_idx_lo]
    se_hi  <- res_high$power_se[,  ref_idx_hi]

    ylim <- c(0, min(1, max(c(pow_lo, pow_hi) + c(se_lo, se_hi), na.rm = TRUE) * 1.1))

    plot(a2_vals, pow_lo, type = "b", pch = 16, lwd = 2.5, col = col,
         ylim = ylim, xlab = expression(alpha[2] ~ "(2nd harmonic relative amplitude)"),
         ylab = "Power (FDR 5%)", las = 1,
         main = sprintf("%s — extreme harmonics\nN=%d; solid=B=%d, dashed=B=%d",
                        dataset_label, ref_n, B_low, B_high))
    polygon(c(a2_vals, rev(a2_vals)),
            c(pow_lo - se_lo, rev(pow_lo + se_lo)),
            col = adjustcolor(col, 0.12), border = NA)
    lines(a2_vals, pow_hi, type = "b", pch = 17, lwd = 2.5, col = col, lty = 2)
    polygon(c(a2_vals, rev(a2_vals)),
            c(pow_hi - se_hi, rev(pow_hi + se_hi)),
            col = adjustcolor(col, 0.08), border = NA)
    abline(h = 0.80, lty = 3, col = "gray50")
    # Vertical reference lines at α₂ breakpoints
    abline(v = 0.5, lty = 2, col = "gray70")
    text(0.51, ylim[2] * 0.95, "standard\nmax", col = "gray50", cex = 0.7, adj = 0)
    legend("topright",
           legend = sprintf("B=%d", c(B_low, B_high)),
           col = col, lwd = 2.5, pch = c(16, 17), lty = c(1, 2), bty = "n")
  } else {
    plot.new(); text(0.5, 0.5, "Results unavailable", cex = 1.2)
  }
  dev.off()
  cat(sprintf("  Figure: %s\n", output_file))
}

# Storage for summary
results_all  <- list()
summary_lines <- c(
  "Fourier Extreme Waveform Robustness",
  sprintf("Mode: %s | NGENES=%d | NSIMS=%d | alpha2: 0 to 1.0",
          if (SMOKE) "SMOKE" else "PRODUCTION", NGENES, NSIMS),
  "alpha3 = 0 throughout (isolates 2nd harmonic effect)",
  ""
)

# =======================================================================
# SECTION 1: MOUSE GSE54651
# =======================================================================
cat("====================================================================\n")
cat("SECTION 1: Mouse GSE54651 — extreme harmonics\n")
cat("====================================================================\n\n")

dat_mouse <- readRDS(DATA_MOUSE)
ct_liv    <- dat_mouse$tod[["LIV"]]
prep_liv  <- prepCircadianData(dat_mouse$count_clean[["LIV"]], times = ct_liv, input_type = "log2")
mat_liv   <- prep_liv$data[rowSums(prep_liv$data > 0) >= 4, , drop = FALSE]
set.seed(1)
mat_liv_s <- mat_liv[sample(nrow(mat_liv), min(NGENES, nrow(mat_liv))), , drop = FALSE]

prep_cer <- prepCircadianData(dat_mouse$count_clean[["CER"]], times = dat_mouse$tod[["CER"]], input_type = "log2")
mat_cer  <- prep_cer$data[rowSums(prep_cer$data > 0) >= 4, , drop = FALSE]
common_m <- intersect(rownames(mat_liv_s), rownames(mat_cer))
fit_liv  <- fitCosinorAll(mat_liv_s[common_m, ], times = ct_liv,       period = 24)
fit_cer  <- fitCosinorAll(mat_cer[common_m, ],   times = prep_cer$times, period = 24)
rhy_liv  <- !is.na(fit_liv$pvalue) & fit_liv$pvalue < RHYTHM_PVAL
rhy_cer  <- !is.na(fit_cer$pvalue) & fit_cer$pvalue < RHYTHM_PVAL

bio_mouse <- estCircadianParam(
  data = mat_liv_s, times = ct_liv, period = 24,
  prop_DR = mean(xor(rhy_liv, rhy_cer)), prop_DP = 0, prop_DA = 0,
  min_rhythm_pval = RHYTHM_PVAL, verbose = FALSE
)
cat(sprintf("Mouse: r_med=%.2f, prop_DR=%.1f%%\n",
            median(as.numeric(fit_liv$A[rhy_liv]) / as.numeric(fit_liv$sigma[rhy_liv]), na.rm=TRUE),
            100 * mean(xor(rhy_liv, rhy_cer))))

s1_lo <- .run_fourier_one_B(bio_mouse, N_GRID_MOUSE, B_VALS_MOUSE[1], label = "sparse B=4")
s1_hi <- .run_fourier_one_B(bio_mouse, N_GRID_MOUSE, B_VALS_MOUSE[2], label = "dense B=8")
.plot_extreme(s1_lo, s1_hi, B_VALS_MOUSE[1], B_VALS_MOUSE[2],
              "GSE54651 (r~2.9)", REF_N_MOUSE, "steelblue",
              file.path(out_dir, "s1_mouse_extreme.pdf"))
saveRDS(list(lo = s1_lo, hi = s1_hi, B_low = B_VALS_MOUSE[1], B_high = B_VALS_MOUSE[2],
             ref_n = REF_N_MOUSE, label = "GSE54651 (r~2.9)", col = "steelblue"),
        file.path(out_dir, "s1_mouse_extreme.rds"))

# Summary stats
a2_vals <- sort(unique(HARM_GRID$alpha2))
for (nm in c("s1_lo", "s1_hi")) {
  res <- get(nm); B <- if (nm == "s1_lo") B_VALS_MOUSE[1] else B_VALS_MOUSE[2]
  if (!is.null(res)) {
    ref_idx <- which.min(abs(res$sample_sizes - REF_N_MOUSE))
    p0 <- res$power_mean[a2_vals == 0,    ref_idx]
    p1 <- res$power_mean[which.max(a2_vals), ref_idx]
    summary_lines <- c(summary_lines,
      sprintf("GSE54651 B=%d (ref_n=%d): power %.0f%% (a2=0) -> %.0f%% (a2=1.0)  delta=%.0f pp",
              B, REF_N_MOUSE, 100*p0, 100*p1, 100*(p0-p1)))
  }
}
summary_lines <- c(summary_lines, "")


# =======================================================================
# SECTION 2: BABOON
# =======================================================================
cat("\n====================================================================\n")
cat("SECTION 2: Baboon LUN vs CER — extreme harmonics\n")
cat("====================================================================\n\n")

load(DATA_BABOON)
bab_expr <- baboon_withTOD$baboon; bab_tod <- baboon_withTOD$tod
prep_lun  <- prepCircadianData(bab_expr[["LUN"]], times = bab_tod[["LUN"]], input_type = "cpm")
prep_cerb <- prepCircadianData(bab_expr[["CER"]], times = bab_tod[["CER"]], input_type = "cpm")
mat_lun   <- prep_lun$data[rowSums(prep_lun$data > 0) >= 6, , drop = FALSE]
mat_cerb  <- prep_cerb$data[rowSums(prep_cerb$data > 0) >= 6, , drop = FALSE]
tod_lun   <- prep_lun$times; tod_cerb <- prep_cerb$times
common_b  <- intersect(rownames(mat_lun), rownames(mat_cerb))
set.seed(2); g_idx_b <- sample(common_b, min(NGENES, length(common_b)))
mat_lun_s <- mat_lun[g_idx_b, , drop = FALSE]; mat_cerb_s <- mat_cerb[g_idx_b, , drop = FALSE]
fit_lun   <- fitCosinorAll(mat_lun_s,  times = tod_lun,  period = 24)
fit_cerb  <- fitCosinorAll(mat_cerb_s, times = tod_cerb, period = 24)
rhy_lun   <- !is.na(fit_lun$pvalue)  & fit_lun$pvalue  < RHYTHM_PVAL
rhy_cerb  <- !is.na(fit_cerb$pvalue) & fit_cerb$pvalue < RHYTHM_PVAL

bio_bab <- estCircadianParam(
  data = mat_lun_s, times = tod_lun, period = 24,
  prop_DR = mean(xor(rhy_lun, rhy_cerb)), prop_DP = 0, prop_DA = 0,
  min_rhythm_pval = RHYTHM_PVAL, verbose = FALSE
)
cat(sprintf("Baboon: r_med=%.2f, prop_DR=%.1f%%\n",
            median(as.numeric(fit_lun$A[rhy_lun]) / as.numeric(fit_lun$sigma[rhy_lun]), na.rm=TRUE),
            100 * mean(xor(rhy_lun, rhy_cerb))))

s2_lo <- .run_fourier_one_B(bio_bab, N_GRID_BABOON, B_VALS_BABOON[1], label = "sparse B=4")
s2_hi <- .run_fourier_one_B(bio_bab, N_GRID_BABOON, B_VALS_BABOON[2], label = "dense B=12")
.plot_extreme(s2_lo, s2_hi, B_VALS_BABOON[1], B_VALS_BABOON[2],
              "Baboon (r~1.7)", REF_N_BABOON, "darkorange",
              file.path(out_dir, "s2_baboon_extreme.pdf"))
saveRDS(list(lo = s2_lo, hi = s2_hi, B_low = B_VALS_BABOON[1], B_high = B_VALS_BABOON[2],
             ref_n = REF_N_BABOON, label = "Baboon (r~1.7)", col = "darkorange"),
        file.path(out_dir, "s2_baboon_extreme.rds"))

for (nm in c("s2_lo", "s2_hi")) {
  res <- get(nm); B <- if (nm == "s2_lo") B_VALS_BABOON[1] else B_VALS_BABOON[2]
  if (!is.null(res)) {
    ref_idx <- which.min(abs(res$sample_sizes - REF_N_BABOON))
    p0 <- res$power_mean[a2_vals == 0,    ref_idx]
    p1 <- res$power_mean[which.max(a2_vals), ref_idx]
    summary_lines <- c(summary_lines,
      sprintf("Baboon   B=%d (ref_n=%d): power %.0f%% (a2=0) -> %.0f%% (a2=1.0)  delta=%.0f pp",
              B, REF_N_BABOON, 100*p0, 100*p1, 100*(p0-p1)))
  }
}
summary_lines <- c(summary_lines, "")


# =======================================================================
# SECTION 3: MOUSE D1D2
# =======================================================================
cat("\n====================================================================\n")
cat("SECTION 3: Mouse D1 vs D2 — extreme harmonics\n")
cat("====================================================================\n\n")

pheno_d1d2 <- read.csv(DATA_D1D2_PHENO, row.names = 1)
prep_d1d2  <- prepCircadianData(DATA_D1D2_EXPR, times = "time", input_type = "counts",
                                pheno = pheno_d1d2, sample_col = "sample")
log_d1d2   <- prep_d1d2$data
d1_samp <- pheno_d1d2$sample[pheno_d1d2$cell == "D1"]
mat_d1  <- log_d1d2[, colnames(log_d1d2) %in% d1_samp, drop = FALSE]
mat_d2  <- log_d1d2[, colnames(log_d1d2) %in% pheno_d1d2$sample[pheno_d1d2$cell == "D2"], drop = FALSE]
tod_d1  <- pheno_d1d2$time[match(colnames(mat_d1), pheno_d1d2$sample)]
keep_d1 <- rowSums(mat_d1 > 1) >= floor(ncol(mat_d1) / 2)
mat_d1  <- mat_d1[keep_d1, , drop = FALSE]; mat_d2 <- mat_d2[keep_d1, , drop = FALSE]
set.seed(4); g_idx_d <- sample(nrow(mat_d1), min(NGENES, nrow(mat_d1)))
mat_d1_s <- mat_d1[g_idx_d, , drop = FALSE]; mat_d2_s <- mat_d2[g_idx_d, , drop = FALSE]
fit_d1   <- fitCosinorAll(mat_d1_s, times = tod_d1, period = 24)
rhy_d1   <- !is.na(fit_d1$pvalue) & fit_d1$pvalue < RHYTHM_PVAL
fit_d2   <- fitCosinorAll(mat_d2_s, times = pheno_d1d2$time[match(colnames(mat_d2_s), pheno_d1d2$sample)], period = 24)
rhy_d2   <- !is.na(fit_d2$pvalue) & fit_d2$pvalue < RHYTHM_PVAL

bio_d1d2 <- estCircadianParam(
  data = mat_d1_s, times = tod_d1, period = 24,
  prop_DR = mean(xor(rhy_d1, rhy_d2)), prop_DP = 0, prop_DA = 0,
  min_rhythm_pval = RHYTHM_PVAL, verbose = FALSE
)
cat(sprintf("D1D2: r_med=%.2f, prop_DR=%.1f%%\n",
            median(as.numeric(fit_d1$A[rhy_d1]) / as.numeric(fit_d1$sigma[rhy_d1]), na.rm=TRUE),
            100 * mean(xor(rhy_d1, rhy_d2))))

s3_lo <- .run_fourier_one_B(bio_d1d2, N_GRID_D1D2, B_VALS_D1D2[1], label = "sparse B=3")
s3_hi <- .run_fourier_one_B(bio_d1d2, N_GRID_D1D2, B_VALS_D1D2[2], label = "dense B=6")
.plot_extreme(s3_lo, s3_hi, B_VALS_D1D2[1], B_VALS_D1D2[2],
              "D1D2 (r~0.66)", REF_N_D1D2, "forestgreen",
              file.path(out_dir, "s3_d1d2_extreme.pdf"))
saveRDS(list(lo = s3_lo, hi = s3_hi, B_low = B_VALS_D1D2[1], B_high = B_VALS_D1D2[2],
             ref_n = REF_N_D1D2, label = "D1D2 (r~0.66)", col = "forestgreen"),
        file.path(out_dir, "s3_d1d2_extreme.rds"))

for (nm in c("s3_lo", "s3_hi")) {
  res <- get(nm); B <- if (nm == "s3_lo") B_VALS_D1D2[1] else B_VALS_D1D2[2]
  if (!is.null(res)) {
    ref_idx <- which.min(abs(res$sample_sizes - REF_N_D1D2))
    p0 <- res$power_mean[a2_vals == 0,    ref_idx]
    p1 <- res$power_mean[which.max(a2_vals), ref_idx]
    summary_lines <- c(summary_lines,
      sprintf("D1D2     B=%d (ref_n=%d): power %.0f%% (a2=0) -> %.0f%% (a2=1.0)  delta=%.0f pp",
              B, REF_N_D1D2, 100*p0, 100*p1, 100*(p0-p1)))
  }
}
summary_lines <- c(summary_lines, "")


# =======================================================================
# SECTION 4: CROSS-DATASET COMPARISON FIGURE
# =======================================================================
cat("\n====================================================================\n")
cat("SECTION 4: Cross-dataset extreme harmonic comparison\n")
cat("====================================================================\n\n")

datasets <- list(
  list(lo = s1_lo, hi = s1_hi, B_low = B_VALS_MOUSE[1],  B_high = B_VALS_MOUSE[2],
       ref_n = REF_N_MOUSE,  label = "GSE54651 (r~2.9)", col = "steelblue"),
  list(lo = s2_lo, hi = s2_hi, B_low = B_VALS_BABOON[1], B_high = B_VALS_BABOON[2],
       ref_n = REF_N_BABOON, label = "Baboon (r~1.7)",   col = "darkorange"),
  list(lo = s3_lo, hi = s3_hi, B_low = B_VALS_D1D2[1],   B_high = B_VALS_D1D2[2],
       ref_n = REF_N_D1D2,   label = "D1D2 (r~0.66)",   col = "forestgreen")
)

fig_s4 <- file.path(out_dir, "s4_extreme_comparison.pdf")
pdf(fig_s4, width = 15, height = 5)
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3.5, 1.5))

for (ds in datasets) {
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
}
dev.off()
cat(sprintf("Figure: %s\n", fig_s4))

# =======================================================================
# WRAP-UP
# =======================================================================
summary_lines <- c(summary_lines,
  "Key interpretation:",
  "  - alpha2=0.5 is the typical upper bound for bulk RNA-seq circadian data",
  "  - alpha2=0.75-1.0 approximates single-cell or spike-like waveforms",
  "  - If power loss is small even at alpha2=1.0, cosinor-DCP is robust",
  "  - If B advantage grows at extreme alpha2, more time points help for non-sinusoidal data"
)
writeLines(summary_lines, file.path(out_dir, "extreme_summary.txt"))
cat("\n====================================================================\n")
cat("07e_fourier_extreme COMPLETE\n")
cat("====================================================================\n\n")
cat(paste(summary_lines, collapse = "\n"), "\n")
cat(sprintf("\nOutput: %s/\nDone.\n", out_dir))
