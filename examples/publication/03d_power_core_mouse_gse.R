#' =======================================================================
#' 03d_power_core_mouse_gse.R
#'    Core Power Analysis — Active Design (Mouse LIV vs CER, GSE54651)
#' =======================================================================
#'
#' PURPOSE
#'   Companion to 03b (Baboon) and 03c (D1D2). Runs DR + DP power pipeline
#'   and B vs m tradeoff bootstrap for Mouse LIV vs CER (Zhang et al. 2014,
#'   GSE54651). Strong signal (r~2.9) provides the high-effect-size anchor
#'   for cross-dataset comparison.
#'
#' DATASET
#'   Mouse Liver (LIV) vs Cerebellum (CER) — Zhang et al. 2014 (GSE54651)
#'     data/mice_GSE54651_CPM.RData
#'     - Active: B=8 ZT time points, m=3 per ZT
#'     - r_median(LIV) ~ 2.9,  prop_DR(LIV vs CER) ~ 27%
#'     - Expected n80 for DR: ~12-24 (very strong signal)
#'
#' OUTPUTS
#'   output/03d_power_core_mouse_gse_<ts>/
#'     dr_power_raw_pvalues.rds
#'     dp_power_raw_pvalues.rds
#'     bm_boot_grid.rds
#'     signal_summary.txt
#'     figures/
#'       dr_power.pdf
#'       dp_power.pdf
#'       bm_tradeoff.pdf
#'
#' USAGE
#'   Rscript examples/publication/03d_power_core_mouse_gse.R
#'   POWERSIM_SMOKE=1 Rscript examples/publication/03d_power_core_mouse_gse.R

# =====================================================================
# SECTION 1: SETUP & CONFIGURATION
# =====================================================================

SMOKE        <- nzchar(Sys.getenv("POWERSIM_SMOKE"))
NSIMS_CORE   <- if (SMOKE) 5L    else 50L
NGENES_CORE  <- if (SMOKE) 500L  else 5000L
NBOOT        <- if (SMOKE) 5L    else 50L
NSIMS_INNER  <- if (SMOKE) 5L    else 20L
NCORES       <- if (SMOKE) 2L    else 20L
RHYTHM_PVAL  <- 0.05

# DR/DP core power: strong signal, n80 expected ~12-24
N_GRID_CORE <- if (SMOKE) c(12L, 24L, 36L) else
                 c(12L, 24L, 36L, 48L, 60L, 72L)

# B vs m tradeoff: B sweep over {4, 6, 8}; pilot has B=8, so B<=8.
# N must be divisible by LCM(4,6,8)=24.
B_VALS    <- c(4L, 6L, 8L)
N_GRID_BM <- if (SMOKE) c(24L, 48L) else
               c(24L, 48L, 72L, 96L)

POWERSIM_ROOT <- {
  env <- Sys.getenv("POWERSIM_ROOT", unset = "")
  if (nzchar(env)) env else
    "/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim"
}
setwd(POWERSIM_ROOT)

source_dir <- file.path(getwd(), "code")
old_wd <- setwd(source_dir)
source("setup.R")
setwd(old_wd)
source("code/plot_with_se.R")

cat(sprintf("Mode: %s | NGENES=%d | NSIMS=%d | NBOOT=%d | NCORES=%d\n",
            if (SMOKE) "SMOKE" else "PRODUCTION",
            NGENES_CORE, NSIMS_CORE, NBOOT, NCORES))

run_tag  <- format(Sys.time(), "%Y%m%d_%H%M")
base_out <- file.path("output", paste0("03d_power_core_mouse_gse_", run_tag))
fig_dir  <- file.path(base_out, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("Output: %s/\n\n", base_out))

opts_analysis <- CircadianAnalysisOptions(
  alpha           = 0.05,
  p.adjust.method = "BH",
  parallel.ncores = 1,
  amp.cutoff      = 0,
  target_effect   = 0.1,
  fdr_thresholds  = c(0.01, 0.05, 0.10, 0.20),
  reference_n     = 24L,
  phase_shifts    = c(0, 0.5, 1, 2, 4, 6, 8, 10, 12)
)

# =====================================================================
# SECTION 2: LOAD MOUSE GSE54651 PILOT DATA
# =====================================================================

cat("Loading Mouse LIV vs CER (GSE54651)...\n")
dat_mouse <- readRDS("data/mice_GSE54651_CPM.RData")

ct_liv <- dat_mouse$tod[["LIV"]]
ct_cer <- dat_mouse$tod[["CER"]]

prep_liv <- prepCircadianData(dat_mouse$count_clean[["LIV"]],
                               times = ct_liv, input_type = "log2")
prep_cer <- prepCircadianData(dat_mouse$count_clean[["CER"]],
                               times = ct_cer, input_type = "log2")
rm(dat_mouse)

keep_liv <- rowSums(prep_liv$data > 0) >= 4
keep_cer <- rowSums(prep_cer$data > 0) >= 4
common_g <- intersect(rownames(prep_liv$data)[keep_liv],
                       rownames(prep_cer$data)[keep_cer])

set.seed(7)
g_idx   <- sample(common_g, min(NGENES_CORE, length(common_g)))
mat_liv <- prep_liv$data[g_idx, , drop = FALSE]
mat_cer <- prep_cer$data[g_idx, , drop = FALSE]
tod_liv <- prep_liv$times
tod_cer <- prep_cer$times
rm(prep_liv, prep_cer)

cat(sprintf("  LIV: %d genes x %d samples\n", nrow(mat_liv), ncol(mat_liv)))
cat(sprintf("  CER: %d genes x %d samples\n", nrow(mat_cer), ncol(mat_cer)))
cat(sprintf("  ZT points (LIV): %s\n\n",
            paste(sort(unique(tod_liv)), collapse = ", ")))

# =====================================================================
# SECTION 3: ESTIMATE PARAMETERS
# =====================================================================

cat("Estimating circadian parameters from Mouse LIV vs CER (two-group)...\n")

opts_bio <- estCircadianParamTwoGroup(
  data_1          = mat_liv,
  data_2          = mat_cer,
  times_1         = tod_liv,
  times_2         = tod_cer,
  period          = 24,
  min_rhythm_pval = RHYTHM_PVAL,
  verbose         = FALSE
)
opts_bio <- updateBioOptions(opts_bio, ngenes = NGENES_CORE)

prop_DR_liv <- opts_bio$prop_DR
fit_liv     <- fitCosinorAll(mat_liv, times = tod_liv, period = 24)
rhy_liv     <- !is.na(fit_liv$pvalue) & fit_liv$pvalue < RHYTHM_PVAL
r_liv_med   <- median(as.numeric(fit_liv$A[rhy_liv]) /
                      as.numeric(fit_liv$sigma[rhy_liv]), na.rm = TRUE)
rm(fit_liv)

cat(sprintf("  LIV: prop_rhythmic=%.1f%%  r_median=%.2f\n",
            mean(rhy_liv) * 100, r_liv_med))
cat(sprintf("  prop_DR (LIV vs CER, two-group): %.1f%%\n\n", prop_DR_liv * 100))

writeLines(c(
  "Mouse LIV vs CER (GSE54651) — circadian signal summary",
  sprintf("Design: active, B=%d ZT points, m=%d per ZT",
          length(unique(tod_liv)),
          round(ncol(mat_liv) / length(unique(tod_liv)))),
  sprintf("LIV: prop_rhythmic=%.1f%%  r_median=%.2f  (p<%.2f)",
          mean(rhy_liv)*100, r_liv_med, RHYTHM_PVAL),
  sprintf("prop_DR (LIV vs CER, two-group): %.1f%%", prop_DR_liv*100),
  "",
  "Context (active trio):",
  sprintf("  Mouse LIV vs CER (active, this):    r~%.1f,  prop_DR~%.0f%%",
          r_liv_med, prop_DR_liv*100),
  "  Baboon LUN vs CER (active m=1):     r~1.7,  prop_DR~41%",
  "  Mouse D1 vs D2 (active m~8):        r~0.65, prop_DR~20%"
), file.path(base_out, "signal_summary.txt"))

design_vec_full <- sort(unique(tod_liv))

# =====================================================================
# SECTION 4: DR POWER ANALYSIS
# =====================================================================

cat("\n====================================================================\n")
cat("SECTION 4: Differential Rhythmicity (DR) — Mouse LIV vs CER\n")
cat("====================================================================\n\n")

opts_bio_DR <- updateBioOptions(opts_bio,
  prop_DR    = prop_DR_liv,
  prop_DP    = 0.00,
  prop_DA    = 0.00,
  phase_diff = c(0, 0),
  amp_diff   = c(1, 1)
)

opts_design_DR <- CircadianDesignOptions(
  sample_sizes = N_GRID_CORE,
  nsims        = NSIMS_CORE,
  design       = "active",
  cts          = design_vec_full,
  test_types   = "DR"
)

dr_power_raw <- runPowerAnalysis(opts_bio_DR, opts_design_DR, opts_analysis,
                                  test_type = "DR")

dr_file <- file.path(base_out, "dr_power_raw_pvalues.rds")
save(dr_power_raw, file = dr_file)
cat(sprintf("DR results saved: %s\n", dr_file))

plotWithSE(dr_file, file.path(fig_dir, "dr_power.pdf"),
           test_name = "DR", analysis.opts = opts_analysis)
cat(sprintf("Figure: %s/figures/dr_power.pdf\n", base_out))

# =====================================================================
# SECTION 5: DP POWER ANALYSIS
# =====================================================================

cat("\n====================================================================\n")
cat("SECTION 5: Differential Phase (DP) — Mouse LIV vs CER\n")
cat("====================================================================\n\n")

opts_bio_DP <- updateBioOptions(opts_bio,
  prop_DR    = 0.00,
  prop_DP    = 0.10,
  prop_DA    = 0.00,
  phase_diff = c(-6, 6),
  amp_diff   = c(1, 1)
)

opts_design_DP <- CircadianDesignOptions(
  sample_sizes = N_GRID_CORE,
  nsims        = NSIMS_CORE,
  design       = "active",
  cts          = design_vec_full,
  test_types   = "DP"
)

dp_power_raw <- runPowerAnalysis(opts_bio_DP, opts_design_DP, opts_analysis,
                                  test_type = "DP")

dp_file <- file.path(base_out, "dp_power_raw_pvalues.rds")
save(dp_power_raw, file = dp_file)
cat(sprintf("DP results saved: %s\n", dp_file))

plotWithSE(dp_file, file.path(fig_dir, "dp_power.pdf"),
           test_name = "DP", analysis.opts = opts_analysis)
cat(sprintf("Figure: %s/figures/dp_power.pdf\n", base_out))

# =====================================================================
# SECTION 6: B vs m TRADEOFF (bootstrap design grid)
# =====================================================================

cat("\n====================================================================\n")
cat("SECTION 6: B vs m tradeoff — Mouse LIV (sweep B=4, 6, 12)\n")
cat("====================================================================\n\n")

cat(sprintf("B values: %s\n", paste(B_VALS, collapse = ", ")))
cat(sprintf("N grid:   %s\n\n", paste(N_GRID_BM, collapse = ", ")))

boot_opts <- CircadianBootstrapOptions(
  design        = "active",
  design_vector = design_vec_full,
  B_values      = B_VALS,
  N_values      = N_GRID_BM,
  nboot         = NBOOT,
  nsims_inner   = NSIMS_INNER,
  seed          = 42L
)

bio_bm <- updateBioOptions(opts_bio,
  prop_DR    = prop_DR_liv,
  prop_DP    = 0.00,
  prop_DA    = 0.00,
  phase_diff = c(0, 0),
  amp_diff   = c(1, 1)
)

cat("Running bootstrap design grid...\n")
bm_boot <- runBootstrapDesignGrid(
  pilot_data    = mat_liv,
  pilot_times   = tod_liv,
  boot.opts     = boot_opts,
  analysis.opts = opts_analysis,
  bio_diff.opts = bio_bm,
  verbose       = TRUE,
  mc.cores      = NCORES
)

saveRDS(bm_boot, file.path(base_out, "bm_boot_grid.rds"))
cat(sprintf("Bootstrap grid saved: %s/bm_boot_grid.rds\n", base_out))

summaryBootstrapDesignGrid(bm_boot, test_type = "DR", fdr_threshold = 0.05)

# B vs m tradeoff plot
n_N <- length(bm_boot$N_values)
n_B <- length(bm_boot$B_values)
pm  <- matrix(bm_boot$power_mean[,,1], nrow = n_N, ncol = n_B)
lo  <- matrix(bm_boot$power_ci_lo[,,1], nrow = n_N, ncol = n_B)
hi  <- matrix(bm_boot$power_ci_hi[,,1], nrow = n_N, ncol = n_B)

cols_b <- c("steelblue", "darkorange", "forestgreen")
fig_bm <- file.path(fig_dir, "bm_tradeoff.pdf")
pdf(fig_bm, width = 7, height = 5)
par(mar = c(4.5, 4.5, 3.5, 1.5))
plot(NA, xlim = range(bm_boot$N_values), ylim = c(0, 1),
     xlab = "N per group", ylab = "Power (FDR 5%)",
     main = sprintf("Mouse LIV vs CER: B vs m tradeoff (DR)\nprop_DR=%.0f%%, r~%.2f",
                    prop_DR_liv * 100, r_liv_med),
     las = 1)
abline(h = 0.80, lty = 2, col = "gray40")
abline(h = c(0.2, 0.4, 0.6), lty = 3, col = "gray85")
for (b in seq_len(n_B)) {
  polygon(c(bm_boot$N_values, rev(bm_boot$N_values)),
          c(lo[, b], rev(hi[, b])),
          col = adjustcolor(cols_b[b], 0.15), border = NA)
  lines(bm_boot$N_values,  pm[, b], col = cols_b[b], lwd = 2, lty = b)
  points(bm_boot$N_values, pm[, b], col = cols_b[b], pch = 15 + b, cex = 0.9)
}
m_approx <- pmax(1L, round(min(N_GRID_BM) / bm_boot$B_values))
legend("bottomright",
       legend = sprintf("B=%d (m~%d repl/ZT)", bm_boot$B_values, m_approx),
       col = cols_b, lwd = 2, lty = seq_len(n_B), pch = 15 + seq_len(n_B),
       bty = "n", cex = 0.85)
dev.off()
cat(sprintf("Figure: %s\n", fig_bm))

# n80 per B
cat("\n--- n80 per B value (bootstrap median, FDR 5%) ---\n")
for (b in seq_len(n_B)) {
  n80_b  <- N_GRID_BM[which(pm[, b]  >= 0.80)[1]]
  n80_lo <- N_GRID_BM[which(hi[, b]  >= 0.80)[1]]
  n80_hi <- N_GRID_BM[which(lo[, b]  >= 0.80)[1]]
  cat(sprintf("  B=%2d: n80=%s  [95%% CI: %s, %s]\n",
              bm_boot$B_values[b],
              ifelse(is.na(n80_b),  ">max", n80_b),
              ifelse(is.na(n80_lo), ">max", n80_lo),
              ifelse(is.na(n80_hi), ">max", n80_hi)))
}

cat("\n====================================================================\n")
cat("03d_power_core_mouse_gse COMPLETE\n")
cat("====================================================================\n")
cat(sprintf("Output: %s/\n", base_out))
cat("Done.\n")
