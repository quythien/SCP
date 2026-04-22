#' =======================================================================
#' 07c_fourier_d1d2.R — Fourier Robustness: Mouse D1 vs D2
#' =======================================================================
#' Split from 07_fourier_robustness.R. Runs Section 3 (D1D2) only.
#' See 07a_fourier_mouse.R header for parallel launch instructions.

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

HARM_GRID <- if (SMOKE) {
  expand.grid(alpha2 = c(0, 0.5), alpha3 = c(0, 0.5))
} else {
  expand.grid(alpha2 = c(0, 0.25, 0.5), alpha3 = c(0, 0.25, 0.5))
}

N_GRID_D1D2     <- if (SMOKE) c(24L, 48L, 72L) else c(24L, 48L, 72L, 96L, 120L)
B_VALS_D1D2     <- c(3L, 6L)
REF_N_D1D2      <- 72L
DATA_D1D2_PHENO <- "data/mouse_clinicalinfo_03082021_rmOutliers.csv"
DATA_D1D2_EXPR  <- "data/mouse_D1D2_logCPMfiltered_counts.csv"

RUN_TAG <- Sys.getenv("RUN_TAG", format(Sys.Date(), "%Y%m%d"))
out_dir <- file.path("output", paste0("07_fourier_robustness_", RUN_TAG))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Mode: %s | NGENES=%d | NSIMS=%d | harm_combos=%d\n",
            if (SMOKE) "SMOKE" else "PRODUCTION", NGENES, NSIMS, nrow(HARM_GRID)))
cat(sprintf("Output -> %s/\n\n", out_dir))

source_dir <- file.path(getwd(), "code")
old_wd <- setwd(source_dir); source("setup.R"); setwd(old_wd)

opts_analysis <- CircadianAnalysisOptions(
  alpha = 0.05, p.adjust.method = "BH",
  fdr_thresholds = c(0.05), reference_n = 48L
)

.run_fourier_one_B <- function(bio_opts, N_grid, B, pilot_tod_full,
                                design_type = "active", label = "") {
  cts_B <- seq(0, 24, length.out = B + 1)[seq_len(B)]
  design_opts <- CircadianDesignOptions(
    sample_sizes = N_grid, nsims = NSIMS, design = design_type,
    cts = if (design_type == "active") cts_B else pilot_tod_full,
    test_types = "DR"
  )
  cat(sprintf("    B=%d: %s\n", B, label))
  tryCatch(
    runFourierDeviationPower(bio_opts, design_opts, opts_analysis,
                             harmonic_grid = HARM_GRID, test_type = "DR",
                             verbose = FALSE),
    error = function(e) { warning(sprintf("B=%d failed: %s", B, e$message)); NULL }
  )
}

.plot_dataset_fourier <- function(res_low, res_high, B_low, B_high,
                                   dataset_label, ref_n, output_file) {
  pdf(output_file, width = 6, height = 5)
  par(mfrow = c(1, 1), mar = c(4.5, 4.5, 3.5, 1.5))
  if (!is.null(res_low) && !is.null(res_high)) {
    hg <- res_low$harmonic_grid
    a3_0 <- hg$alpha3 == 0
    a2_vals <- sort(unique(hg$alpha2[a3_0]))
    ref_idx_lo <- which.min(abs(res_low$sample_sizes  - ref_n))
    ref_idx_hi <- which.min(abs(res_high$sample_sizes - ref_n))
    pow_lo <- res_low$power_mean[a3_0,  ref_idx_lo]
    pow_hi <- res_high$power_mean[a3_0, ref_idx_hi]
    se_lo  <- res_low$power_se[a3_0,   ref_idx_lo]
    se_hi  <- res_high$power_se[a3_0,  ref_idx_hi]
    plot(a2_vals, pow_lo, type = "b", pch = 16, lwd = 2, col = "steelblue",
         ylim = c(0, 1), xlab = expression(alpha[2] ~ "(2nd harmonic)"),
         ylab = "Power (FDR 5%)", las = 1,
         main = sprintf("%s\nB vs harmonics at N=%d (a3=0)", dataset_label, ref_n))
    arrows(a2_vals, pow_lo - se_lo, a2_vals, pow_lo + se_lo,
           length = 0.05, angle = 90, code = 3, col = "steelblue", lwd = 1.2)
    lines(a2_vals, pow_hi, type = "b", pch = 17, lwd = 2, col = "darkorange", lty = 2)
    arrows(a2_vals, pow_hi - se_hi, a2_vals, pow_hi + se_hi,
           length = 0.05, angle = 90, code = 3, col = "darkorange", lwd = 1.2)
    abline(h = 0.80, lty = 3, col = "gray50")
    legend("bottomleft", legend = sprintf("B=%d", c(B_low, B_high)),
           col = c("steelblue", "darkorange"), lwd = 2,
           pch = c(16, 17), lty = c(1, 2), bty = "n")
  } else {
    plot.new(); text(0.5, 0.5, "Comparison unavailable", cex = 1.2)
  }
  dev.off()
  cat(sprintf("  Figure: %s\n", output_file))
}

# =======================================================================
# SECTION 3: MOUSE D1D2 — D1 pilot
# =======================================================================
cat("====================================================================\n")
cat("SECTION 3: Mouse D1 vs D2 Fourier robustness\n")
cat("====================================================================\n\n")

pheno_d1d2 <- read.csv(DATA_D1D2_PHENO, row.names = 1)
prep_d1d2  <- prepCircadianData(DATA_D1D2_EXPR, times = "time", input_type = "counts",
                                pheno = pheno_d1d2, sample_col = "sample")
log_d1d2   <- prep_d1d2$data

d1_samp <- pheno_d1d2$sample[pheno_d1d2$cell == "D1"]
d2_samp <- pheno_d1d2$sample[pheno_d1d2$cell == "D2"]
mat_d1  <- log_d1d2[, colnames(log_d1d2) %in% d1_samp, drop = FALSE]
mat_d2  <- log_d1d2[, colnames(log_d1d2) %in% d2_samp, drop = FALSE]
tod_d1  <- pheno_d1d2$time[match(colnames(mat_d1), pheno_d1d2$sample)]
tod_d2  <- pheno_d1d2$time[match(colnames(mat_d2), pheno_d1d2$sample)]

keep_d1  <- rowSums(mat_d1 > 1) >= floor(ncol(mat_d1) / 2)
mat_d1   <- mat_d1[keep_d1, , drop = FALSE]
mat_d2   <- mat_d2[keep_d1, , drop = FALSE]
set.seed(4)
g_idx_d  <- sample(nrow(mat_d1), min(NGENES, nrow(mat_d1)))
mat_d1_s <- mat_d1[g_idx_d, , drop = FALSE]
mat_d2_s <- mat_d2[g_idx_d, , drop = FALSE]

fit_d1   <- fitCosinorAll(mat_d1_s, times = tod_d1, period = 24)
fit_d2   <- fitCosinorAll(mat_d2_s, times = tod_d2, period = 24)
rhy_d1   <- !is.na(fit_d1$pvalue) & fit_d1$pvalue < RHYTHM_PVAL
rhy_d2   <- !is.na(fit_d2$pvalue) & fit_d2$pvalue < RHYTHM_PVAL
prop_DR_d <- mean(xor(rhy_d1, rhy_d2))

bio_d1d2 <- estCircadianParam(
  data = mat_d1_s, times = tod_d1, period = 24,
  prop_DR = prop_DR_d, prop_DP = 0,
  min_rhythm_pval = RHYTHM_PVAL, verbose = FALSE
)
cat(sprintf("D1D2 pilot: n=%d, prop_DR=%.1f%%, r_med=%.2f\n",
            ncol(mat_d1_s), 100 * prop_DR_d,
            median(as.numeric(fit_d1$A[rhy_d1]) /
                   as.numeric(fit_d1$sigma[rhy_d1]), na.rm = TRUE)))

cat("Running Fourier analysis (B=3 and B=6)...\n")
s3_B3 <- .run_fourier_one_B(bio_d1d2, N_GRID_D1D2, B_VALS_D1D2[1], tod_d1, label = "3 ZT (every 8h)")
s3_B6 <- .run_fourier_one_B(bio_d1d2, N_GRID_D1D2, B_VALS_D1D2[2], tod_d1, label = "6 ZT (every 4h)")

.plot_dataset_fourier(s3_B3, s3_B6, B_VALS_D1D2[1], B_VALS_D1D2[2],
                       "Mouse D1 vs D2", REF_N_D1D2,
                       file.path(out_dir, "s3_d1d2_fourier.pdf"))

saveRDS(list(B3 = s3_B3, B6 = s3_B6,
             label = "D1D2 (r~0.66)", col = "forestgreen",
             B_low = B_VALS_D1D2[1], B_high = B_VALS_D1D2[2], ref_n = REF_N_D1D2),
        file.path(out_dir, "s3_d1d2_fourier.rds"))

cat(sprintf("\nDone. Output: %s/\n", out_dir))
