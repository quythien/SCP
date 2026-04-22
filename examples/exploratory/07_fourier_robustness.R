#' =======================================================================
#' 07_fourier_robustness.R — Fourier Waveform Robustness Analysis
#' =======================================================================
#'
#' PURPOSE
#'   Tests how DR power degrades when the true waveform has higher harmonics
#'   that the cosinor detection model cannot fit. The DCP pipeline assumes a
#'   pure sinusoid; real circadian gene expression often has:
#'     - Asymmetric peaks (morning-biased vs evening-biased)
#'     - Sharp transitions (more square-wave-like)
#'     - Two-component oscillations
#'   These are captured by 2nd and 3rd harmonics (α₂, α₃).
#'
#'   The simulation asks two questions for each dataset:
#'
#'   Q1: Power degradation — how much does DR power fall as harmonics increase
#'       (α₂ and α₃ grow from 0 toward 0.75)?
#'       Expectation: moderate degradation at α₂ ≤ 0.25; substantial at α₂ ≥ 0.5.
#'
#'   Q2: B vs m under harmonics — for a fixed N, does using more time points
#'       (higher B) protect against power loss from harmonics?
#'       Expectation: yes, because denser time coverage captures harmonic peaks
#'       that sparse designs miss; sparse designs (low B) lose more power.
#'
#' DATASETS (all active designs; passive Seney omitted — B not identifiable)
#'   1. Mouse GSE54651 — LIV pilot (r~2.9, B∈{4,8})
#'   2. Baboon LUN — pilot (r~1.72, B∈{4,12})
#'   3. Mouse D1D2 — D1 pilot (r~0.66, B∈{3,6})
#'
#' OUTPUTS
#'   output/07_fourier_robustness_<timestamp>/
#'     s1_mouse_fourier.pdf       — Panel A: heatmap, Panel B: power vs N per harmonic
#'     s2_baboon_fourier.pdf
#'     s3_d1d2_fourier.pdf
#'     s4_b_protection.pdf        — power vs α₂ at fixed N, one line per (dataset, B)
#'     fourier_summary.txt        — key numbers
#'
#' USAGE
#'   Rscript examples/exploratory/07_fourier_robustness.R
#'   POWERSIM_SMOKE=1 Rscript examples/exploratory/07_fourier_robustness.R
#'
#' @author Thien Pham

# -----------------------------------------------------------------------
# 0. Path configuration
# -----------------------------------------------------------------------
POWERSIM_ROOT <- {
  env <- Sys.getenv("POWERSIM_ROOT", unset = "")
  if (nzchar(env)) env else
    "/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim"
}
setwd(POWERSIM_ROOT)

# -----------------------------------------------------------------------
# 1. Settings
# -----------------------------------------------------------------------
SMOKE       <- nzchar(Sys.getenv("POWERSIM_SMOKE"))

NGENES      <- if (SMOKE) 300L  else 3000L
NSIMS       <- if (SMOKE) 5L    else 20L
RHYTHM_PVAL <- 0.05

# Harmonic grid: α₂ (2nd) and α₃ (3rd harmonic relative amplitudes)
# Smoke: 2x2 = 4 combinations; Production: 4x3 = 12 combinations
HARM_GRID <- if (SMOKE) {
  expand.grid(alpha2 = c(0, 0.5), alpha3 = c(0, 0.5))
} else {
  expand.grid(alpha2 = c(0, 0.25, 0.5), alpha3 = c(0, 0.25, 0.5))
}

# N grids — must be divisible by max(B) for each dataset
# Mouse:  LCM(4, 8)  = 8  → c(8,16,24,48,72)
# Baboon: LCM(4, 12) = 12 → c(12,24,36,48,72,96)
# D1D2:   LCM(3, 6)  = 6  → c(24,48,72,96,120)
N_GRID_MOUSE  <- if (SMOKE) c(16L, 24L, 32L)     else c(8L, 16L, 24L, 48L, 72L)
N_GRID_BABOON <- if (SMOKE) c(24L, 36L, 48L)     else c(12L, 24L, 36L, 48L, 72L, 96L)
N_GRID_D1D2   <- if (SMOKE) c(24L, 48L, 72L)     else c(24L, 48L, 72L, 96L, 120L)

# B values to compare per dataset (low B = sparse; high B = dense coverage)
B_VALS_MOUSE  <- c(4L, 8L)    # 4 = {ZT 0,6,12,18}; 8 = full 8 CT
B_VALS_BABOON <- c(4L, 12L)   # 4 = every 6h subset; 12 = all ZT
B_VALS_D1D2   <- c(3L, 6L)    # 3 = every 8h; 6 = every 4h

# Reference N for the "B protection" summary plot (Section 4)
REF_N_MOUSE  <- 24L
REF_N_BABOON <- 48L
REF_N_D1D2   <- 72L

DATA_MOUSE  <- "data/mice_GSE54651_CPM.RData"
DATA_BABOON <- "data/CAMO_PRC_hmb.RData"
DATA_D1D2_PHENO <- "data/mouse_clinicalinfo_03082021_rmOutliers.csv"
DATA_D1D2_EXPR  <- "data/mouse_D1D2_logCPMfiltered_counts.csv"

cat(sprintf("Mode: %s | NGENES=%d | NSIMS=%d | harm_combos=%d\n",
            if (SMOKE) "SMOKE" else "PRODUCTION",
            NGENES, NSIMS, nrow(HARM_GRID)))

# -----------------------------------------------------------------------
# 2. Setup
# -----------------------------------------------------------------------
source_dir <- file.path(getwd(), "code")
old_wd <- setwd(source_dir)
source("setup.R")
setwd(old_wd)

ts      <- format(Sys.time(), "%Y%m%d_%H%M")
out_dir <- file.path("output", paste0("07_fourier_robustness_", ts))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("Output -> %s/\n\n", out_dir))

opts_analysis <- CircadianAnalysisOptions(
  alpha           = 0.05,
  p.adjust.method = "BH",
  fdr_thresholds  = c(0.05),
  reference_n     = 48L
)

# Storage for Section 4 cross-dataset comparison
sec4_results <- list()

# -----------------------------------------------------------------------
# Helper: run Fourier robustness for one B value, return power_mean matrix
# -----------------------------------------------------------------------
.run_fourier_one_B <- function(bio_opts, N_grid, B, pilot_tod_full,
                                design_type = "active", label = "") {
  # Build evenly-spaced time points for this B
  cts_B <- seq(0, 24, length.out = B + 1)[seq_len(B)]

  design_opts <- CircadianDesignOptions(
    sample_sizes = N_grid,
    nsims        = NSIMS,
    design       = design_type,
    cts          = if (design_type == "active") cts_B else pilot_tod_full,
    test_types   = "DR"
  )

  cat(sprintf("    B=%d: %s\n", B, label))
  tryCatch(
    runFourierDeviationPower(bio_opts, design_opts, opts_analysis,
                             harmonic_grid = HARM_GRID,
                             test_type     = "DR",
                             verbose       = FALSE),
    error = function(e) {
      warning(sprintf("  B=%d failed: %s", B, e$message))
      NULL
    }
  )
}

# -----------------------------------------------------------------------
# Helper: plot per-dataset Fourier figure (Panel A heatmap + Panel B lines)
# -----------------------------------------------------------------------
.plot_dataset_fourier <- function(res_low, res_high, B_low, B_high,
                                   dataset_label, ref_n, output_file) {
  pdf(output_file, width = 6, height = 5)
  par(mfrow = c(1, 1), mar = c(4.5, 4.5, 3.5, 1.5))

  # Power vs alpha2 at ref_n, one line per B
  # (α₃ = 0 slice for clarity)
  if (!is.null(res_low) && !is.null(res_high)) {
    hg <- res_low$harmonic_grid
    a3_0 <- hg$alpha3 == 0
    a2_vals <- sort(unique(hg$alpha2[a3_0]))

    ref_idx_lo <- which.min(abs(res_low$sample_sizes  - ref_n))
    ref_idx_hi <- which.min(abs(res_high$sample_sizes - ref_n))

    pow_lo    <- res_low$power_mean[a3_0,  ref_idx_lo]
    pow_hi    <- res_high$power_mean[a3_0, ref_idx_hi]
    se_lo     <- res_low$power_se[a3_0,   ref_idx_lo]
    se_hi     <- res_high$power_se[a3_0,  ref_idx_hi]

    plot(a2_vals, pow_lo, type = "b", pch = 16, lwd = 2,
         col = "steelblue",
         ylim = c(0, 1), xlab = expression(alpha[2] ~ "(2nd harmonic)"),
         ylab = "Power (FDR 5%)", las = 1,
         main = sprintf("%s\nB vs harmonics at N=%d (a3=0)", dataset_label, ref_n))
    arrows(a2_vals, pow_lo - se_lo, a2_vals, pow_lo + se_lo,
           length = 0.05, angle = 90, code = 3, col = "steelblue", lwd = 1.2)
    lines(a2_vals, pow_hi, type = "b", pch = 17, lwd = 2, col = "darkorange", lty = 2)
    arrows(a2_vals, pow_hi - se_hi, a2_vals, pow_hi + se_hi,
           length = 0.05, angle = 90, code = 3, col = "darkorange", lwd = 1.2)
    abline(h = 0.80, lty = 3, col = "gray50")
    legend("bottomleft",
           legend = sprintf("B=%d", c(B_low, B_high)),
           col = c("steelblue", "darkorange"),
           lwd = 2, pch = c(16, 17), lty = c(1, 2), bty = "n")
  } else {
    plot.new()
    text(0.5, 0.5, "Comparison unavailable", cex = 1.2)
  }

  dev.off()
  cat(sprintf("  Figure: %s\n", output_file))
}


# =======================================================================
# SECTION 1: MOUSE GSE54651 — LIV pilot
# =======================================================================
cat("\n====================================================================\n")
cat("SECTION 1: Mouse GSE54651 — LIV vs CER Fourier robustness\n")
cat("====================================================================\n\n")

dat_mouse    <- readRDS(DATA_MOUSE)
ct_liv       <- dat_mouse$tod[["LIV"]]
prep_liv     <- prepCircadianData(dat_mouse$count_clean[["LIV"]], times = ct_liv,
                                  input_type = "log2")
mat_liv      <- prep_liv$data
mat_liv      <- mat_liv[rowSums(mat_liv > 0) >= 4, , drop = FALSE]
set.seed(1)
g_idx_m      <- sample(nrow(mat_liv), min(NGENES, nrow(mat_liv)))
mat_liv_s    <- mat_liv[g_idx_m, , drop = FALSE]

# Estimate prop_DR (LIV vs CER)
prep_cer     <- prepCircadianData(dat_mouse$count_clean[["CER"]], times = dat_mouse$tod[["CER"]],
                                  input_type = "log2")
mat_cer      <- prep_cer$data
mat_cer      <- mat_cer[rowSums(mat_cer > 0) >= 4, , drop = FALSE]
common_m     <- intersect(rownames(mat_liv_s), rownames(mat_cer))
fit_liv      <- fitCosinorAll(mat_liv_s[common_m, ], times = ct_liv, period = 24)
fit_cer      <- fitCosinorAll(mat_cer[common_m, prep_cer$times == prep_cer$times[1] |
                                TRUE, drop = FALSE],
                              times = prep_cer$times, period = 24)
rhy_liv      <- !is.na(fit_liv$pvalue) & fit_liv$pvalue < RHYTHM_PVAL
rhy_cer      <- !is.na(fit_cer$pvalue) & fit_cer$pvalue < RHYTHM_PVAL
prop_DR_m    <- mean(xor(rhy_liv, rhy_cer))

bio_mouse <- estCircadianParam(
  data = mat_liv_s, times = ct_liv, period = 24,
  prop_DR = prop_DR_m, prop_DP = 0,
  min_rhythm_pval = RHYTHM_PVAL, verbose = FALSE
)
cat(sprintf("Mouse pilot: n=%d, prop_DR=%.1f%%, r_med=%.2f\n",
            ncol(mat_liv_s), 100 * prop_DR_m,
            median(as.numeric(fit_liv$A[rhy_liv]) /
                   as.numeric(fit_liv$sigma[rhy_liv]), na.rm = TRUE)))

cat("Running Fourier analysis (B=4 and B=8)...\n")
s1_B4  <- .run_fourier_one_B(bio_mouse, N_GRID_MOUSE, B_VALS_MOUSE[1], ct_liv,
                               label = "sparse (every 6h)")
s1_B8  <- .run_fourier_one_B(bio_mouse, N_GRID_MOUSE, B_VALS_MOUSE[2], ct_liv,
                               label = "full 8 CT")

.plot_dataset_fourier(s1_B4, s1_B8, B_VALS_MOUSE[1], B_VALS_MOUSE[2],
                       "GSE54651 LIV vs CER", REF_N_MOUSE,
                       file.path(out_dir, "s1_mouse_fourier.pdf"))

saveRDS(list(B4 = s1_B4, B8 = s1_B8), file.path(out_dir, "s1_mouse_fourier.rds"))
sec4_results[["mouse"]] <- list(
  label  = "GSE54651 (r~2.9)", col = "steelblue",
  B_low  = B_VALS_MOUSE[1], B_high = B_VALS_MOUSE[2],
  ref_n  = REF_N_MOUSE,
  res_lo = s1_B4, res_hi = s1_B8
)


# =======================================================================
# SECTION 2: BABOON — LUN pilot
# =======================================================================
cat("\n====================================================================\n")
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
mat_lun_s <- mat_lun[g_idx_b, , drop = FALSE]
mat_cerb_s <- mat_cerb[g_idx_b, , drop = FALSE]

fit_lun   <- fitCosinorAll(mat_lun_s,  times = tod_lun,  period = 24)
fit_cerb  <- fitCosinorAll(mat_cerb_s, times = tod_cerb, period = 24)
rhy_lun   <- !is.na(fit_lun$pvalue)  & fit_lun$pvalue  < RHYTHM_PVAL
rhy_cerb  <- !is.na(fit_cerb$pvalue) & fit_cerb$pvalue < RHYTHM_PVAL
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
s2_B4  <- .run_fourier_one_B(bio_bab, N_GRID_BABOON, B_VALS_BABOON[1], tod_lun,
                               label = "4 ZT (every 6h)")
s2_B12 <- .run_fourier_one_B(bio_bab, N_GRID_BABOON, B_VALS_BABOON[2], tod_lun,
                               label = "12 ZT (full coverage)")

.plot_dataset_fourier(s2_B4, s2_B12, B_VALS_BABOON[1], B_VALS_BABOON[2],
                       "Baboon LUN vs CER", REF_N_BABOON,
                       file.path(out_dir, "s2_baboon_fourier.pdf"))

saveRDS(list(B4 = s2_B4, B12 = s2_B12), file.path(out_dir, "s2_baboon_fourier.rds"))
sec4_results[["baboon"]] <- list(
  label  = "Baboon (r~1.7)", col = "darkorange",
  B_low  = B_VALS_BABOON[1], B_high = B_VALS_BABOON[2],
  ref_n  = REF_N_BABOON,
  res_lo = s2_B4, res_hi = s2_B12
)


# =======================================================================
# SECTION 3: MOUSE D1D2 — D1 pilot
# =======================================================================
cat("\n====================================================================\n")
cat("SECTION 3: Mouse D1 vs D2 Fourier robustness\n")
cat("====================================================================\n\n")

pheno_d1d2 <- read.csv(DATA_D1D2_PHENO, row.names = 1)
prep_d1d2  <- prepCircadianData(DATA_D1D2_EXPR, times = "time", input_type = "counts",
                                pheno = pheno_d1d2, sample_col = "sample")
log_d1d2   <- prep_d1d2$data

d1_samp  <- pheno_d1d2$sample[pheno_d1d2$cell == "D1"]
d2_samp  <- pheno_d1d2$sample[pheno_d1d2$cell == "D2"]
mat_d1   <- log_d1d2[, colnames(log_d1d2) %in% d1_samp, drop = FALSE]
mat_d2   <- log_d1d2[, colnames(log_d1d2) %in% d2_samp, drop = FALSE]
tod_d1   <- pheno_d1d2$time[match(colnames(mat_d1), pheno_d1d2$sample)]
tod_d2   <- pheno_d1d2$time[match(colnames(mat_d2), pheno_d1d2$sample)]

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
s3_B3  <- .run_fourier_one_B(bio_d1d2, N_GRID_D1D2, B_VALS_D1D2[1], tod_d1,
                               label = "3 ZT (every 8h)")
s3_B6  <- .run_fourier_one_B(bio_d1d2, N_GRID_D1D2, B_VALS_D1D2[2], tod_d1,
                               label = "6 ZT (every 4h)")

.plot_dataset_fourier(s3_B3, s3_B6, B_VALS_D1D2[1], B_VALS_D1D2[2],
                       "Mouse D1 vs D2", REF_N_D1D2,
                       file.path(out_dir, "s3_d1d2_fourier.pdf"))

saveRDS(list(B3 = s3_B3, B6 = s3_B6), file.path(out_dir, "s3_d1d2_fourier.rds"))
sec4_results[["d1d2"]] <- list(
  label  = "D1D2 (r~0.66)", col = "forestgreen",
  B_low  = B_VALS_D1D2[1], B_high = B_VALS_D1D2[2],
  ref_n  = REF_N_D1D2,
  res_lo = s3_B3, res_hi = s3_B6
)


# =======================================================================
# SECTION 4: CROSS-DATASET SUMMARY — B protection against harmonics
# =======================================================================
# For each dataset: power vs α₂ at reference_n, B_low vs B_high (α₃=0 slice).
# Shows whether B improvement is universal or dataset-specific.

cat("\n====================================================================\n")
cat("SECTION 4: Cross-dataset B protection summary\n")
cat("====================================================================\n\n")

fig_s4 <- file.path(out_dir, "s4_b_protection.pdf")
pdf(fig_s4, width = 5 * length(sec4_results), height = 5)
par(mfrow = c(1, length(sec4_results)), mar = c(4.5, 4.5, 3.5, 1.5))

summary_lines <- c(
  "Fourier Robustness — B Protection Summary",
  sprintf("Mode: %s | NGENES=%d | NSIMS=%d", if (SMOKE) "SMOKE" else "PRODUCTION", NGENES, NSIMS),
  sprintf("Harmonic combinations: %d", nrow(HARM_GRID)),
  ""
)

for (nm in names(sec4_results)) {
  sc <- sec4_results[[nm]]
  res_lo <- sc$res_lo; res_hi <- sc$res_hi
  ref_n  <- sc$ref_n

  if (is.null(res_lo) || is.null(res_hi)) {
    plot.new(); title(sprintf("%s\nresults unavailable", sc$label))
    next
  }

  hg <- res_lo$harmonic_grid
  a3_0 <- hg$alpha3 == 0
  a2_vals <- sort(unique(hg$alpha2[a3_0]))

  ref_idx_lo <- which.min(abs(res_lo$sample_sizes - ref_n))
  ref_idx_hi <- which.min(abs(res_hi$sample_sizes - ref_n))

  pow_lo <- res_lo$power_mean[a3_0, ref_idx_lo]
  pow_hi <- res_hi$power_mean[a3_0, ref_idx_hi]
  se_lo  <- res_lo$power_se[a3_0,  ref_idx_lo]
  se_hi  <- res_hi$power_se[a3_0,  ref_idx_hi]

  # Power at pure cosinor (α₂=0) and at max harmonic
  pow_lo_0 <- pow_lo[a2_vals == 0]
  pow_hi_0 <- pow_hi[a2_vals == 0]
  pow_lo_max <- pow_lo[which.max(a2_vals)]
  pow_hi_max <- pow_hi[which.max(a2_vals)]

  delta_lo <- if (length(pow_lo_0) && length(pow_lo_max)) pow_lo_0 - pow_lo_max else NA
  delta_hi <- if (length(pow_hi_0) && length(pow_hi_max)) pow_hi_0 - pow_hi_max else NA

  summary_lines <- c(summary_lines,
    sprintf("%s (ref_n=%d):", sc$label, ref_n),
    sprintf("  B=%d: power %.0f%% (α₂=0) → %.0f%% (α₂=%.2g)   Δ=%.0f pp",
            sc$B_low,  100*pow_lo_0, 100*pow_lo_max, max(a2_vals), 100*delta_lo),
    sprintf("  B=%d: power %.0f%% (α₂=0) → %.0f%% (α₂=%.2g)   Δ=%.0f pp",
            sc$B_high, 100*pow_hi_0, 100*pow_hi_max, max(a2_vals), 100*delta_hi),
    sprintf("  B advantage: B=%d loses %.0f pp less than B=%d",
            sc$B_high, 100*(delta_lo - delta_hi), sc$B_low),
    ""
  )

  plot(a2_vals, pow_lo * 100, type = "b", pch = 16, lwd = 2.5,
       col = sc$col,
       ylim = c(0, 100), las = 1,
       xlab = expression(alpha[2] ~ "(2nd harmonic, a3=0)"),
       ylab = "Power (%, FDR 5%)",
       main = sprintf("%s\nN=%d: does B protect?", sc$label, ref_n))
  arrows(a2_vals, (pow_lo - se_lo) * 100, a2_vals, (pow_lo + se_lo) * 100,
         length = 0.05, angle = 90, code = 3, col = sc$col, lwd = 1.2)
  lines(a2_vals, pow_hi * 100, type = "b", pch = 17, lwd = 2.5,
        col = sc$col, lty = 2)
  arrows(a2_vals, (pow_hi - se_hi) * 100, a2_vals, (pow_hi + se_hi) * 100,
         length = 0.05, angle = 90, code = 3, col = sc$col, lwd = 1.2)
  abline(h = 80, lty = 3, col = "gray50")
  text(min(a2_vals), 82, "80%", col = "gray50", cex = 0.8, adj = 0)
  legend("bottomleft",
         legend = sprintf("B=%d", c(sc$B_low, sc$B_high)),
         col = sc$col, lwd = 2.5, pch = c(16, 17), lty = c(1, 2),
         bty = "n", cex = 0.85)
}

dev.off()
cat(sprintf("Figure: %s\n", fig_s4))


# =======================================================================
# WRAP-UP
# =======================================================================
writeLines(summary_lines, file.path(out_dir, "fourier_summary.txt"))
cat("\n====================================================================\n")
cat("07_fourier_robustness COMPLETE\n")
cat("====================================================================\n\n")
cat(paste(summary_lines, collapse = "\n"), "\n")
cat(sprintf("Output: %s/\n", out_dir))
cat("Done.\n")
