#' =======================================================================
#' 07b_fourier_baboon.R — Fourier Robustness: Baboon LUN vs CER
#' =======================================================================
#' Split from 07_fourier_robustness.R. Runs Section 2 (Baboon) only.
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

N_GRID_BABOON <- if (SMOKE) c(24L, 36L, 48L) else c(12L, 24L, 36L, 48L, 72L, 96L)
B_VALS_BABOON <- c(4L, 12L)
REF_N_BABOON  <- 48L
DATA_BABOON   <- "data/CAMO_PRC_hmb.RData"

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
# SECTION 2: BABOON — LUN pilot
# =======================================================================
cat("====================================================================\n")
cat("SECTION 2: Baboon LUN vs CER Fourier robustness\n")
cat("====================================================================\n\n")

load(DATA_BABOON)
bab_expr  <- baboon_withTOD$baboon
bab_tod   <- baboon_withTOD$tod

prep_lun  <- prepCircadianData(bab_expr[["LUN"]], times = bab_tod[["LUN"]], input_type = "cpm")
prep_cerb <- prepCircadianData(bab_expr[["CER"]], times = bab_tod[["CER"]], input_type = "cpm")
mat_lun   <- prep_lun$data[rowSums(prep_lun$data > 0) >= 6, , drop = FALSE]
mat_cerb  <- prep_cerb$data[rowSums(prep_cerb$data > 0) >= 6, , drop = FALSE]
tod_lun   <- prep_lun$times
tod_cerb  <- prep_cerb$times

common_b  <- intersect(rownames(mat_lun), rownames(mat_cerb))
set.seed(2)
g_idx_b   <- sample(common_b, min(NGENES, length(common_b)))
mat_lun_s  <- mat_lun[g_idx_b, , drop = FALSE]
mat_cerb_s <- mat_cerb[g_idx_b, , drop = FALSE]

fit_lun  <- fitCosinorAll(mat_lun_s,  times = tod_lun,  period = 24)
fit_cerb <- fitCosinorAll(mat_cerb_s, times = tod_cerb, period = 24)
rhy_lun  <- !is.na(fit_lun$pvalue)  & fit_lun$pvalue  < RHYTHM_PVAL
rhy_cerb <- !is.na(fit_cerb$pvalue) & fit_cerb$pvalue < RHYTHM_PVAL
prop_DR_b <- mean(xor(rhy_lun, rhy_cerb))

bio_bab <- estCircadianParam(
  data = mat_lun_s, times = tod_lun, period = 24,
  prop_DR = prop_DR_b, prop_DP = 0,
  min_rhythm_pval = RHYTHM_PVAL, verbose = FALSE
)
cat(sprintf("Baboon pilot: n=%d, prop_DR=%.1f%%, r_med=%.2f\n",
            ncol(mat_lun_s), 100 * prop_DR_b,
            median(as.numeric(fit_lun$A[rhy_lun]) /
                   as.numeric(fit_lun$sigma[rhy_lun]), na.rm = TRUE)))

cat("Running Fourier analysis (B=4 and B=12)...\n")
s2_B4  <- .run_fourier_one_B(bio_bab, N_GRID_BABOON, B_VALS_BABOON[1], tod_lun, label = "4 ZT (every 6h)")
s2_B12 <- .run_fourier_one_B(bio_bab, N_GRID_BABOON, B_VALS_BABOON[2], tod_lun, label = "12 ZT (full)")

.plot_dataset_fourier(s2_B4, s2_B12, B_VALS_BABOON[1], B_VALS_BABOON[2],
                       "Baboon LUN vs CER", REF_N_BABOON,
                       file.path(out_dir, "s2_baboon_fourier.pdf"))

saveRDS(list(B4 = s2_B4, B12 = s2_B12,
             label = "Baboon (r~1.7)", col = "darkorange",
             B_low = B_VALS_BABOON[1], B_high = B_VALS_BABOON[2], ref_n = REF_N_BABOON),
        file.path(out_dir, "s2_baboon_fourier.rds"))

cat(sprintf("\nDone. Output: %s/\n", out_dir))
