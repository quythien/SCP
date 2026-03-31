#!/usr/bin/env Rscript
#' =======================================================================
#' 03c_power_core_mouse.R
#'    Core Power Analysis — Active Multi-Replicate Design (Mouse D1 vs D2)
#' =======================================================================
#'
#' PURPOSE
#'   Companion to 03_power_core.R (passive / human aging) and
#'   03b_power_core_active.R (active single-replicate / Baboon).
#'   Runs the DR + DP power pipeline on an active multi-replicate pilot
#'   (Mouse D1 vs D2 cell types, B=6 ZT points, m~7-8 per ZT) and adds
#'   a B vs m tradeoff section.
#'
#'   Together the 03 trio answers:
#'     03:  Passive design (human aging) — baseline cost of uncontrolled TOD
#'     03b: Active, m=1 (Baboon) — benefit of temporal control
#'     03c: Active, m=8 (Mouse D1D2) — benefit of biological replicates per ZT
#'          and the B vs m tradeoff question
#'
#' DATASET
#'   Mouse D1 vs D2 striatal cell types
#'     data/mouse_clinicalinfo_03082021_rmOutliers.csv  — sample metadata
#'     data/mouse_D1D2_logCPMfiltered_counts.csv        — raw counts
#'     - Active: B=6 ZT (2,6,10,14,18,22h every 4h), m~7-8 per ZT per group
#'     - r_median(D1) ~ 0.65,  prop_DR(D1 vs D2) ~ 18-22%
#'     - Expected n80 for DR: ~80-120
#'
#' OUTPUTS
#'   output/03c_power_core_mouse_<ts>/
#'     signal_summary.txt
#'     dr_power_raw_pvalues.rds
#'     dp_power_raw_pvalues.rds
#'     bm_boot_grid.rds
#'     figures/
#'       dr_power.pdf
#'       dp_power.pdf
#'       bm_tradeoff.pdf
#'
#' USAGE
#'   Rscript examples/publication/03c_power_core_mouse.R
#'   POWERSIM_SMOKE=1 Rscript examples/publication/03c_power_core_mouse.R
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

# DR/DP core power: N grid (active, moderate signal — n80 expected ~80-120)
N_GRID_CORE <- if (SMOKE) c(24L, 48L, 72L) else
                 c(24L, 48L, 72L, 96L, 120L, 150L)

# B vs m tradeoff: B=3 (every 8h) vs B=6 (every 4h, full pilot coverage)
# N must be divisible by 6
B_VALS    <- c(3L, 6L)
N_GRID_BM <- if (SMOKE) c(24L, 48L, 72L) else
               c(24L, 48L, 72L, 96L, 120L, 150L)

DATA_PHENO <- "data/mouse_clinicalinfo_03082021_rmOutliers.csv"
DATA_EXPR  <- "data/mouse_D1D2_logCPMfiltered_counts.csv"

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
base_out <- file.path("output", paste0("03c_power_core_mouse_", run_tag))
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
  reference_n     = 96L,
  phase_shifts    = c(0, 0.5, 1, 2, 4, 6, 8, 10, 12)
)


# =====================================================================
# SECTION 2: LOAD MOUSE D1 vs D2 PILOT DATA
# =====================================================================

cat("Loading Mouse D1 vs D2 pilot data...\n")

pheno <- read.csv(DATA_PHENO, row.names = 1)

# Raw counts (despite filename); prepCircadianData normalizes to log2(CPM+1)
prep    <- prepCircadianData(DATA_EXPR, times = "time", input_type = "counts",
                             pheno = pheno, sample_col = "sample")
log_mat <- prep$data

# Split by cell type
d1_samp <- pheno$sample[pheno$cell == "D1"]
d2_samp <- pheno$sample[pheno$cell == "D2"]
mat_d1  <- log_mat[, colnames(log_mat) %in% d1_samp, drop = FALSE]
mat_d2  <- log_mat[, colnames(log_mat) %in% d2_samp, drop = FALSE]
tod_d1  <- pheno$time[match(colnames(mat_d1), pheno$sample)]
tod_d2  <- pheno$time[match(colnames(mat_d2), pheno$sample)]

# Filter: keep genes expressed in ≥ half of D1 samples
keep   <- rowSums(mat_d1 > 1) >= floor(ncol(mat_d1) / 2)
mat_d1 <- mat_d1[keep, , drop = FALSE]
mat_d2 <- mat_d2[keep, , drop = FALSE]

# Subsample to NGENES_CORE
set.seed(42)
g_idx  <- sample(nrow(mat_d1), min(NGENES_CORE, nrow(mat_d1)))
mat_d1 <- mat_d1[g_idx, , drop = FALSE]
mat_d2 <- mat_d2[g_idx, , drop = FALSE]

cat(sprintf("  D1: %d genes x %d samples\n", nrow(mat_d1), ncol(mat_d1)))
cat(sprintf("  D2: %d genes x %d samples\n", nrow(mat_d2), ncol(mat_d2)))
cat(sprintf("  ZT points (D1): %s\n\n",
            paste(sort(unique(tod_d1)), collapse = ", ")))


# =====================================================================
# SECTION 3: ESTIMATE PARAMETERS
# =====================================================================

cat("Estimating circadian parameters from Mouse D1 vs D2 (two-group)...\n")

opts_bio <- estCircadianParamTwoGroup(
  data_1          = mat_d1,
  data_2          = mat_d2,
  times_1         = tod_d1,
  times_2         = tod_d2,
  period          = 24,
  min_rhythm_pval = RHYTHM_PVAL,
  verbose         = FALSE
)
opts_bio <- updateBioOptions(opts_bio, ngenes = NGENES_CORE)

prop_DR_d1d2 <- opts_bio$prop_DR
fit_d1       <- fitCosinorAll(mat_d1, times = tod_d1, period = 24)
rhy_d1       <- !is.na(fit_d1$pvalue) & fit_d1$pvalue < RHYTHM_PVAL
r_d1_med     <- median(as.numeric(fit_d1$A[rhy_d1]) /
                       as.numeric(fit_d1$sigma[rhy_d1]), na.rm = TRUE)
rm(fit_d1)

cat(sprintf("  D1: prop_rhythmic=%.1f%%  r_median=%.2f\n",
            mean(rhy_d1) * 100, r_d1_med))
cat(sprintf("  prop_DR (D1 vs D2, two-group): %.1f%%\n\n", prop_DR_d1d2 * 100))

writeLines(c(
  "Mouse D1 vs D2 — circadian signal summary",
  sprintf("Design: active, B=%d ZT points, m~%d per ZT",
          length(unique(tod_d1)), round(ncol(mat_d1) / length(unique(tod_d1)))),
  sprintf("D1: prop_rhythmic=%.1f%%  r_median=%.2f  (p<%.2f)",
          mean(rhy_d1)*100, r_d1_med, RHYTHM_PVAL),
  sprintf("prop_DR (D1 vs D2, two-group): %.1f%%", prop_DR_d1d2*100),
  "",
  "Context (active trio):",
  "  Baboon LUN vs CER (active m=1):     r~1.7,  prop_DR~41%, n80 expected ~24-36",
  "  Mouse D1 vs D2 (active m~8, this):  r~0.65, prop_DR~20%, n80 expected ~80-120",
  "  Human aging young vs old (passive): r~0.5,  prop_DR~14%, n80 expected >300"
), file.path(base_out, "signal_summary.txt"))

# Design vector: full B=6 ZT coverage from pilot
design_vec_full <- sort(unique(tod_d1))


# =====================================================================
# SECTION 4: DR POWER ANALYSIS
# =====================================================================

cat("\n====================================================================\n")
cat("SECTION 4: Differential Rhythmicity (DR) — active design, B=6\n")
cat("====================================================================\n\n")

opts_bio_DR <- updateBioOptions(opts_bio,
  prop_DR    = prop_DR_d1d2,
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
cat("SECTION 5: Differential Phase (DP) — active design, B=6\n")
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
# time points (B=6, fewer replicates per ZT) or concentrate replicates
# (B=3, more replicates per ZT)?
# B=3 → every 8h (sparse),  B=6 → every 4h (full pilot coverage)

cat("\n====================================================================\n")
cat("SECTION 6: B vs m tradeoff — Mouse D1 vs D2 (sweep B=3, 6)\n")
cat("====================================================================\n\n")

cat(sprintf("B values: %s\n", paste(B_VALS, collapse=", ")))
cat(sprintf("N grid:   %s  (all divisible by 6)\n\n", paste(N_GRID_BM, collapse=", ")))

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
  prop_DR    = prop_DR_d1d2,
  prop_DP    = 0.00,
  prop_DA    = 0.00,
  phase_diff = c(0, 0),
  amp_diff   = c(1, 1)
)

cat("Running bootstrap design grid (mc.cores=32)...\n")
bm_boot <- runBootstrapDesignGrid(
  pilot_data    = mat_d1,
  pilot_times   = tod_d1,
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
n_N     <- length(bm_boot$N_values)
n_B     <- length(bm_boot$B_values)
fdr_idx <- which.min(abs(bm_boot$fdr_thresholds - 0.05))

pm <- matrix(bm_boot$power_mean[, , fdr_idx], nrow = n_N, ncol = n_B)
lo <- matrix(bm_boot$power_ci_lo[, , fdr_idx], nrow = n_N, ncol = n_B)
hi <- matrix(bm_boot$power_ci_hi[, , fdr_idx], nrow = n_N, ncol = n_B)

cols_b <- c("steelblue", "darkorange")

fig_bm <- file.path(fig_dir, "bm_tradeoff.pdf")
pdf(fig_bm, width = 7, height = 5)
par(mar = c(4.5, 4.5, 3.5, 1.5))

plot(NA, xlim = range(bm_boot$N_values), ylim = c(0, 1),
     xlab = "N per group", ylab = "Power (FDR 5%)",
     main = sprintf("Mouse D1 vs D2: B vs m tradeoff (DR)\nprop_DR=%.0f%%, r~%.2f",
                    prop_DR_d1d2 * 100, r_d1_med),
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

cat("\n--- n80 per B value (bootstrap median, FDR 5%) ---\n")
for (b in seq_len(n_B)) {
  n80_b  <- N_GRID_BM[which(pm[, b] >= 0.80)[1]]
  n80_lo <- N_GRID_BM[which(hi[, b] >= 0.80)[1]]
  n80_hi <- N_GRID_BM[which(lo[, b] >= 0.80)[1]]
  cat(sprintf("  B=%2d: n80=%s  [95%% CI: %s, %s]\n",
              bm_boot$B_values[b],
              ifelse(is.na(n80_b),  ">max", n80_b),
              ifelse(is.na(n80_lo), ">max", n80_lo),
              ifelse(is.na(n80_hi), ">max", n80_hi)))
}

cat("\n====================================================================\n")
cat("03c_power_core_mouse COMPLETE\n")
cat("====================================================================\n\n")
cat(sprintf("Output: %s/\n", base_out))
cat("Done.\n")
