#' =======================================================================
#' 03b_power_core_active.R
#'    Core Power Analysis — Active Design (Baboon LUN vs CER)
#' =======================================================================
#'
#' PURPOSE
#'   Companion to 03_power_core.R (passive / human aging). Runs the same
#'   DR + DP power pipeline on an active-design pilot (Baboon LUN vs CER,
#'   B=12 ZT points, r~1.7) and adds a B vs m tradeoff section.
#'
#'   Together with 03_power_core.R this answers:
#'     Q1: How does active vs passive design affect DR/DP power curves?
#'         (active with strong signal needs far fewer subjects)
#'     Q2: For an active design, is it better to spread subjects across
#'         more time points (B↑, m↓) or concentrate replicates (B↓, m↑)?
#'
#' DATASET
#'   Baboon LUN vs CER (Mure et al. 2018)
#'     data/CAMO_PRC_hmb.RData  →  baboon_withTOD$baboon / $tod
#'     - Active: B=12 ZT time points (every ~2h), m=1 per ZT per tissue
#'     - r_median(LUN) ~ 1.72,  prop_DR(LUN vs CER) ~ 41%
#'     - Expected n80 for DR: ~24-36 (strong signal)
#'
#' OUTPUTS
#'   output/03b_power_core_active_<ts>/
#'     dr_power_raw_pvalues.rds
#'     dp_power_raw_pvalues.rds
#'     signal_summary.txt
#'     figures/
#'       dr_power.pdf
#'       dp_power.pdf
#'       bm_tradeoff.pdf
#'
#' USAGE
#'   Rscript examples/publication/03b_power_core_active.R
#'   POWERSIM_SMOKE=1 Rscript examples/publication/03b_power_core_active.R
#'
#' @author Thien Pham

# =====================================================================
# SECTION 1: SETUP & CONFIGURATION
# =====================================================================

SMOKE        <- nzchar(Sys.getenv("POWERSIM_SMOKE"))
NSIMS_CORE   <- if (SMOKE) 5L    else 50L
NGENES_CORE  <- if (SMOKE) 500L  else 5000L
NBOOT        <- if (SMOKE) 5L    else 50L
NSIMS_INNER  <- if (SMOKE) 5L    else 20L
RHYTHM_PVAL  <- 0.05

# DR/DP core power: N grid (active, strong signal — n80 expected ~24-36)
N_GRID_CORE <- if (SMOKE) c(12L, 24L, 36L) else
                 c(12L, 24L, 36L, 48L, 60L, 72L)

# B vs m tradeoff: B sweep over {4, 6, 12}; N must be divisible by 12
B_VALS      <- c(4L, 6L, 12L)
N_GRID_BM   <- if (SMOKE) c(24L, 36L, 48L) else
                 c(24L, 36L, 48L, 60L, 72L, 96L)

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

cat(sprintf("Mode: %s | NGENES=%d | NSIMS=%d | NBOOT=%d\n",
            if (SMOKE) "SMOKE" else "PRODUCTION",
            NGENES_CORE, NSIMS_CORE, NBOOT))

run_tag  <- format(Sys.time(), "%Y%m%d_%H%M")
base_out <- file.path("output", paste0("03b_power_core_active_", run_tag))
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
  reference_n     = 36L,
  phase_shifts    = c(0, 0.5, 1, 2, 4, 6, 8, 10, 12)
)

# =====================================================================
# SECTION 2: LOAD BABOON PILOT DATA
# =====================================================================

cat("Loading Baboon LUN vs CER pilot data...\n")
load("data/CAMO_PRC_hmb.RData")
bab_expr <- baboon_withTOD$baboon
bab_tod  <- baboon_withTOD$tod
rm(baboon_withTOD)

prep_lun  <- prepCircadianData(bab_expr[["LUN"]], times = bab_tod[["LUN"]], input_type = "cpm")
prep_cer  <- prepCircadianData(bab_expr[["CER"]], times = bab_tod[["CER"]], input_type = "cpm")
rm(bab_expr, bab_tod)

# Filter: keep genes expressed in ≥ 6 samples in each tissue
keep_lun  <- rowSums(prep_lun$data > 0) >= 6
keep_cer  <- rowSums(prep_cer$data > 0) >= 6
common_g  <- intersect(rownames(prep_lun$data)[keep_lun],
                        rownames(prep_cer$data)[keep_cer])

set.seed(7)
g_idx     <- sample(common_g, min(NGENES_CORE, length(common_g)))
mat_lun   <- prep_lun$data[g_idx, , drop = FALSE]
mat_cer   <- prep_cer$data[g_idx, , drop = FALSE]
tod_lun   <- prep_lun$times
tod_cer   <- prep_cer$times
rm(prep_lun, prep_cer)

cat(sprintf("  LUN: %d genes x %d samples\n", nrow(mat_lun), ncol(mat_lun)))
cat(sprintf("  CER: %d genes x %d samples\n", nrow(mat_cer), ncol(mat_cer)))
cat(sprintf("  ZT points (LUN): %s\n\n",
            paste(sort(unique(tod_lun)), collapse = ", ")))


# =====================================================================
# SECTION 3: ESTIMATE PARAMETERS
# =====================================================================

cat("Estimating circadian parameters from Baboon LUN vs CER (two-group)...\n")

# Use two-group estimation so prop_rhythmic is computed from the union of both
# tissues, correctly accommodating genes that are rhythmic in CER but not LUN.
opts_bio <- estCircadianParamTwoGroup(
  data_1          = mat_lun,
  data_2          = mat_cer,
  times_1         = tod_lun,
  times_2         = tod_cer,
  period          = 24,
  min_rhythm_pval = RHYTHM_PVAL,
  verbose         = FALSE
)
opts_bio <- updateBioOptions(opts_bio, ngenes = NGENES_CORE)

prop_DR_b <- opts_bio$prop_DR
fit_lun   <- fitCosinorAll(mat_lun, times = tod_lun, period = 24)
rhy_lun   <- !is.na(fit_lun$pvalue) & fit_lun$pvalue < RHYTHM_PVAL
r_lun_med <- median(as.numeric(fit_lun$A[rhy_lun]) /
                    as.numeric(fit_lun$sigma[rhy_lun]), na.rm = TRUE)
rm(fit_lun)

cat(sprintf("  LUN: prop_rhythmic=%.1f%%  r_median=%.2f\n",
            mean(rhy_lun) * 100, r_lun_med))
cat(sprintf("  prop_DR (LUN vs CER, two-group): %.1f%%\n\n", prop_DR_b * 100))

# Signal summary
writeLines(c(
  "Baboon LUN vs CER — circadian signal summary",
  sprintf("Design: active, B=%d ZT points, m=1 per ZT", length(unique(tod_lun))),
  sprintf("LUN: prop_rhythmic=%.1f%%  r_median=%.2f  (p<%.2f)",
          mean(rhy_lun)*100, r_lun_med, RHYTHM_PVAL),
  sprintf("prop_DR (LUN vs CER, two-group): %.1f%%", prop_DR_b*100),
  "",
  "Context (active vs passive):",
  "  Baboon LUN vs CER (active, this):    r~1.7,  prop_DR~41%, n80 expected ~24-36",
  "  Human aging young vs old (passive):  r~0.5,  prop_DR~14%, n80 expected >300"
), file.path(base_out, "signal_summary.txt"))

# Design vector: full B=12 ZT coverage from pilot
design_vec_full <- sort(unique(tod_lun))


# =====================================================================
# SECTION 4: DR POWER ANALYSIS
# =====================================================================

cat("\n====================================================================\n")
cat("SECTION 4: Differential Rhythmicity (DR) — active design, B=12\n")
cat("====================================================================\n\n")

opts_bio_DR <- updateBioOptions(opts_bio,
  prop_DR    = prop_DR_b,
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
cat("SECTION 5: Differential Phase (DP) — active design, B=12\n")
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
# Question: At fixed N, is it better to spread subjects across more ZT
# time points (larger B, fewer replicates per ZT) or concentrate them
# (smaller B, more replicates per ZT)?
# B=4 → ~every 6h (sparse),  B=6 → every 4h,  B=12 → every 2h (dense)

cat("\n====================================================================\n")
cat("SECTION 6: B vs m tradeoff — Baboon LUN (sweep B=4, 6, 12)\n")
cat("====================================================================\n\n")

cat(sprintf("B values: %s\n", paste(B_VALS, collapse=", ")))
cat(sprintf("N grid:   %s  (all divisible by 12)\n\n", paste(N_GRID_BM, collapse=", ")))

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
  prop_DR    = prop_DR_b,
  prop_DP    = 0.00,
  prop_DA    = 0.00,
  phase_diff = c(0, 0),
  amp_diff   = c(1, 1)
)

cat("Running bootstrap design grid...\n")
bm_boot <- runBootstrapDesignGrid(
  pilot_data    = mat_lun,
  pilot_times   = tod_lun,
  boot.opts     = boot_opts,
  analysis.opts = opts_analysis,
  bio_diff.opts = bio_bm,
  verbose       = TRUE,
  mc.cores      = 32L
)

saveRDS(bm_boot, file.path(base_out, "bm_boot_grid.rds"))
cat(sprintf("Bootstrap grid saved: %s/bm_boot_grid.rds\n", base_out))

summaryBootstrapDesignGrid(bm_boot, test_type = "DR", fdr_threshold = 0.05)

# B vs m tradeoff plot
n_N   <- length(bm_boot$N_values)
n_B   <- length(bm_boot$B_values)
fdr_idx <- which.min(abs(bm_boot$fdr_thresholds - 0.05))

pm <- matrix(bm_boot$power_mean[, , fdr_idx], nrow = n_N, ncol = n_B)
lo <- matrix(bm_boot$power_ci_lo[, , fdr_idx], nrow = n_N, ncol = n_B)
hi <- matrix(bm_boot$power_ci_hi[, , fdr_idx], nrow = n_N, ncol = n_B)

cols_b <- c("steelblue", "darkorange", "forestgreen")

fig_bm <- file.path(fig_dir, "bm_tradeoff.pdf")
pdf(fig_bm, width = 7, height = 5)
par(mar = c(4.5, 4.5, 3.5, 1.5))

plot(NA, xlim = range(bm_boot$N_values), ylim = c(0, 1),
     xlab = "N per group", ylab = "Power (FDR 5%)",
     main = sprintf("Baboon LUN vs CER: B vs m tradeoff (DR)\nprop_DR=%.0f%%, r~%.2f",
                    prop_DR_b * 100, r_lun_med),
     las = 1)
abline(h = 0.80, lty = 2, col = "gray40")
abline(h = c(0.2, 0.4, 0.6), lty = 3, col = "gray85")

for (b in seq_len(n_B)) {
  col_b <- cols_b[b]
  polygon(c(bm_boot$N_values, rev(bm_boot$N_values)),
          c(lo[, b], rev(hi[, b])),
          col = adjustcolor(col_b, 0.15), border = NA)
  lines(bm_boot$N_values, pm[, b], col = col_b, lwd = 2, lty = b)
  points(bm_boot$N_values, pm[, b], col = col_b, pch = 15 + b, cex = 0.9)
}

# Label: B value and approximate m (replicates per ZT) at minimum N
m_approx <- pmax(1L, round(min(N_GRID_BM) / bm_boot$B_values))
legend("bottomright",
       legend = sprintf("B=%d (m~%d repl/ZT)", bm_boot$B_values, m_approx),
       col = cols_b, lwd = 2, lty = seq_len(n_B), pch = 15 + seq_len(n_B),
       bty = "n", cex = 0.85)

dev.off()
cat(sprintf("Figure: %s\n", fig_bm))


# =====================================================================
# WRAP-UP
# =====================================================================

# n80 per B value
cat("\n--- n80 per B value (bootstrap median, FDR 5%) ---\n")
for (b in seq_len(n_B)) {
  n80_b <- N_GRID_BM[which(pm[, b] >= 0.80)[1]]
  n80_lo <- N_GRID_BM[which(hi[, b] >= 0.80)[1]]
  n80_hi <- N_GRID_BM[which(lo[, b] >= 0.80)[1]]
  cat(sprintf("  B=%2d: n80=%s  [95%% CI: %s, %s]\n",
              bm_boot$B_values[b],
              ifelse(is.na(n80_b),  ">max", n80_b),
              ifelse(is.na(n80_lo), ">max", n80_lo),
              ifelse(is.na(n80_hi), ">max", n80_hi)))
}

cat("\n====================================================================\n")
cat("03b_power_core_active COMPLETE\n")
cat("====================================================================\n\n")
cat(sprintf("Output: %s/\n", base_out))
cat("Done.\n")
