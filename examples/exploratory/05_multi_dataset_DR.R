#' =======================================================================
#' 05_multi_dataset_DR.R — Multi-Dataset DR Power Comparison
#' =======================================================================
#'
#' PURPOSE
#'   Compare differential rhythmicity (DR) power across four biological
#'   scenarios spanning a wide range of signal strength and design types:
#'
#'   1. Mouse GSE54651 — LIV vs CER (strongest: r~2.9, active CT, m=1)
#'   2. Baboon cross-tissue — LUN vs CER (moderate: r~1.6, active ZT, m=1)
#'   3. Seney MDD vs Control ACC (weak: r~0.55, passive post-mortem) [PRIMARY]
#'      Alt: Human aging young vs old BA11/BA47 (r~0.5, passive) — see
#'           "SECTION 3 ALTERNATIVE" block below; swap by toggling if(FALSE)
#'   4. Mouse D1 vs D2 cell types (moderate: r~0.66, active ZT, m~8)
#'
#'   For each scenario the bootstrap design grid answers:
#'     Q1: What N achieves 80% power? (bootstrap median n80)
#'     Q2: How uncertain is that n80? (bootstrap 95% CI)
#'     Q3: For fixed N, is more time-point coverage (B up) or more
#'         replicates (m up) better? (power vs N, stratified by B)
#'
#'   Q2 subsumes Q1: the bootstrap median IS the point estimate, with
#'   added uncertainty quantification from pilot data resampling.
#'
#' OUTPUTS
#'   output/05_multi_dataset_DR_<timestamp>/
#'     s0_ct_justification.txt          (mouse: 4 vs 8 CT quality table)
#'     s1_mouse_boot_grid.rds
#'     s2_baboon_boot_grid.rds
#'     s3_human_boot_grid.rds
#'     figures/s1_mouse_bootstrap.pdf
#'     figures/s2_baboon_bootstrap.pdf
#'     figures/s3_human_bootstrap.pdf
#'     figures/s4_three_way_comparison.pdf
#'     figures/s5_bm_tradeoff.pdf
#'
#' DATA
#'   Mouse:  mice_GSE54651_CPM.RData  (readRDS; count_clean; ENSMUSG IDs)
#'   Baboon: CAMO_PRC_hmb.RData       (load; baboon_withTOD; ENSG IDs)
#'   Human:  combined_data.rds        (readRDS; expr + pheno; ENSG IDs)
#'
#' USAGE
#'   Rscript examples/exploratory/05_multi_dataset_DR.R
#'   POWERSIM_SMOKE=1 Rscript examples/exploratory/05_multi_dataset_DR.R
#'
#' @author Thien Pham

# =====================================================================
# SETTINGS
# =====================================================================

SMOKE       <- nzchar(Sys.getenv("POWERSIM_SMOKE"))

NGENES      <- if (SMOKE) 500L  else 5000L
NBOOT       <- if (SMOKE) 3L    else 100L
NSIMS_INNER <- if (SMOKE) 3L    else 20L

N_GRID_MOUSE  <- if (SMOKE) c(8L, 12L, 24L)          else c(8L, 12L, 24L, 36L, 48L, 60L)
N_GRID_BABOON <- if (SMOKE) c(12L, 24L, 36L)          else c(12L, 24L, 36L, 48L, 60L, 72L, 96L)
N_GRID_HUMAN  <- if (SMOKE) c(40L, 60L, 100L)         else c(40L, 60L, 100L, 160L, 200L, 300L)
# Note: all N must be divisible by B=4 (passive design uses B=4 as fixed placeholder;
# framework requires N %% B == 0). Values 40,60,100,160,200,300 all satisfy this.

# D1D2 mouse: active, B=6 ZT (2,6,10,14,18,22), m~8 replicates/ZT
# All N must be divisible by max(B_VALS_D1D2)=6
N_GRID_D1D2   <- if (SMOKE) c(24L, 36L, 48L)          else c(24L, 36L, 48L, 60L, 72L, 96L, 120L)
B_VALS_D1D2   <- c(3L, 6L)   # B=3 (every 8h) vs B=6 (every 4h)

DATA_D1D2_PHENO <- "data/mouse_clinicalinfo_03082021_rmOutliers.csv"
DATA_D1D2_EXPR  <- "data/mouse_D1D2_logCPMfiltered_counts.csv"

# B sweep: number of distinct time points to test (active designs only)
# Mouse:  4 = {ZT 4,10,16,22}; 8 = treat all 8 CT as distinct (full dataset)
# Baboon: 4,6,8,12 = subsets of the 12 ZT time points
B_VALS_MOUSE  <- c(4L, 8L)
B_VALS_BABOON <- c(4L, 6L, 8L, 12L)

RHYTHM_PVAL <- 0.05   # p-value threshold for rhythmic classification

DATA_MOUSE  <- "data/mice_GSE54651_CPM.RData"
DATA_BABOON <- "data/CAMO_PRC_hmb.RData"
DATA_HUMAN  <- "data/combined_data.rds"   # kept for reference; Seney used as passive example
DATA_SENEY_META <- "data/MD5_MetaData_1-15-25.xlsx"
DATA_SENEY_TOD  <- "data/TOD.xlsx"
DATA_SENEY_EXPR <- "data/ACC_RNA_filtered_normalized.csv"

cat(sprintf("Mode: %s | NGENES=%d | NBOOT=%d | NSIMS_INNER=%d\n",
            if (SMOKE) "SMOKE" else "PRODUCTION",
            NGENES, NBOOT, NSIMS_INNER))

# =====================================================================
# PATH CONFIGURATION
# =====================================================================
# Set POWERSIM_ROOT as an environment variable for portability.
# This lets you run the script on any machine (local or server) without
# editing the code — just set the variable before launching R:
#
#   Linux/Mac shell:  export POWERSIM_ROOT=/path/to/PowerSim
#   R console:        Sys.setenv(POWERSIM_ROOT="/path/to/PowerSim")
#   Slurm/PBS header: #SBATCH --export=POWERSIM_ROOT=/path/to/PowerSim
#
# If POWERSIM_ROOT is not set, the local development path below is used.

POWERSIM_ROOT <- {
  env <- Sys.getenv("POWERSIM_ROOT", unset = "")
  if (nzchar(env)) env else
    "/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim"
}

# =====================================================================
# SETUP
# =====================================================================

POWERSIM_DIR <- POWERSIM_ROOT
setwd(POWERSIM_DIR)
source_dir <- file.path(getwd(), "code")
old_wd <- setwd(source_dir)
source("setup.R")
setwd(old_wd)

out_dir <- file.path("output", sprintf("05_multi_dataset_DR_%s", format(Sys.time(), "%Y%m%d_%H%M")))
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("Output: %s/\n\n", out_dir))

# Shared analysis options
opts_analysis <- CircadianAnalysisOptions(
  alpha           = 0.05,
  p.adjust.method = "BH",
  fdr_thresholds  = c(0.01, 0.05, 0.10),
  reference_n     = 60
)

# =====================================================================
# SECTION 0: MOUSE CT JUSTIFICATION
# =====================================================================
# Compare parameter estimation quality: all 8 CT vs day1 (4 CT) vs day2 (4 CT).
# With only 4 time points and 3 cosinor parameters, residual df = 1 → very unstable.
# Conclusion: always use all 8 CT together (5 residual df).

cat("====================================================================\n")
cat("SECTION 0: Mouse CT justification (all 8 vs day-by-day)\n")
cat("====================================================================\n\n")

dat_mouse <- readRDS(DATA_MOUSE)
ct_all <- dat_mouse$tod[["LIV"]]   # 22,28,34,40,46,52,58,64
prep_liv_full <- prepCircadianData(dat_mouse$count_clean[["LIV"]], times = ct_all,
                                   input_type = "log2", verbose = FALSE)
mat_liv_full  <- prep_liv_full$data

day1_idx <- ct_all <= 40   # CT 22,28,34,40
day2_idx <- ct_all >= 46   # CT 46,52,58,64

cat(sprintf("CT all   (n=%d): %s\n", sum(!logical(length(ct_all))), paste(ct_all, collapse=",")))
cat(sprintf("Day1     (n=%d): %s\n", sum(day1_idx), paste(ct_all[day1_idx], collapse=",")))
cat(sprintf("Day2     (n=%d): %s\n\n", sum(day2_idx), paste(ct_all[day2_idx], collapse=",")))

# Filter and subsample
mat_liv_full <- mat_liv_full[rowSums(mat_liv_full > 0) >= 4, , drop=FALSE]
set.seed(42)
scan_idx <- sample(nrow(mat_liv_full), min(2000L, nrow(mat_liv_full)))
mat_scan <- mat_liv_full[scan_idx, , drop=FALSE]

.ct_summary <- function(mat, tod, label, n_params=3) {
  fit <- fitCosinorAll(mat, times=tod, period=24)
  rhy <- !is.na(fit$pvalue) & fit$pvalue < RHYTHM_PVAL
  r   <- as.numeric(fit$A[rhy]) / as.numeric(fit$sigma[rhy])
  df  <- ncol(mat) - n_params
  cat(sprintf("  %-12s: n=%d  resid_df=%d  prop_rhy=%.1f%%  r_median=%.2f  [df=1 → unstable]\n",
              label, ncol(mat), df, mean(rhy)*100, median(r, na.rm=TRUE)))
  invisible(list(prop_rhy=mean(rhy), r_median=median(r,na.rm=TRUE), df=df))
}

cat("LIV cosinor fit quality:\n")
q_all  <- .ct_summary(mat_scan,               ct_all,           "All 8 CT")
q_day1 <- .ct_summary(mat_scan[, day1_idx],   ct_all[day1_idx], "Day1 (n=4)")
q_day2 <- .ct_summary(mat_scan[, day2_idx],   ct_all[day2_idx], "Day2 (n=4)")

cat("\nConclusion: all 8 CT used throughout this script (residual df=5 vs 1).\n")

# Save justification text
s0_lines <- c(
  "Mouse CT justification: all 8 CT vs day-by-day",
  sprintf("All 8 CT: n=8, resid_df=5, prop_rhy=%.1f%%, r_median=%.2f", q_all$prop_rhy*100, q_all$r_median),
  sprintf("Day1 (n=4): resid_df=1, prop_rhy=%.1f%%, r_median=%.2f", q_day1$prop_rhy*100, q_day1$r_median),
  sprintf("Day2 (n=4): resid_df=1, prop_rhy=%.1f%%, r_median=%.2f", q_day2$prop_rhy*100, q_day2$r_median),
  "Decision: use all 8 CT (df=5 >> df=1)"
)
writeLines(s0_lines, file.path(out_dir, "s0_ct_justification.txt"))

# =====================================================================
# SECTION 1: MOUSE — LIV vs CER
# =====================================================================

cat("\n====================================================================\n")
cat("SECTION 1: Mouse GSE54651 — LIV vs CER (bootstrap design grid)\n")
cat("====================================================================\n\n")

# Load tissues (already log2-scaled; prepCircadianData validates + coerces to matrix)
prep_liv <- prepCircadianData(dat_mouse$count_clean[["LIV"]], times = dat_mouse$tod[["LIV"]],
                              input_type = "log2")
prep_cer <- prepCircadianData(dat_mouse$count_clean[["CER"]], times = dat_mouse$tod[["CER"]],
                              input_type = "log2")
mat_liv <- prep_liv$data; ct_liv <- prep_liv$times
mat_cer <- prep_cer$data; ct_cer <- prep_cer$times

# Filter expressed genes
mat_liv <- mat_liv[rowSums(mat_liv > 0) >= 4, , drop=FALSE]
mat_cer <- mat_cer[rowSums(mat_cer > 0) >= 4, , drop=FALSE]

# Common genes, subsample
common_m <- intersect(rownames(mat_liv), rownames(mat_cer))
cat(sprintf("Common genes (LIV ∩ CER): %d\n", length(common_m)))
set.seed(1)
g_idx_m <- sample(common_m, min(NGENES, length(common_m)))

mat_liv_s <- mat_liv[g_idx_m, , drop=FALSE]
mat_cer_s <- mat_cer[g_idx_m, , drop=FALSE]

# Estimate prop_DR from pilot cosinor fits
fit_liv <- fitCosinorAll(mat_liv_s, times=ct_liv, period=24)
fit_cer <- fitCosinorAll(mat_cer_s, times=ct_cer, period=24)
rhy_liv <- !is.na(fit_liv$pvalue) & fit_liv$pvalue < RHYTHM_PVAL
rhy_cer <- !is.na(fit_cer$pvalue) & fit_cer$pvalue < RHYTHM_PVAL
dr_mask_m <- xor(rhy_liv, rhy_cer)
prop_DR_mouse <- mean(dr_mask_m)
r_liv_med <- median(as.numeric(fit_liv$A[rhy_liv]) / as.numeric(fit_liv$sigma[rhy_liv]), na.rm=TRUE)

cat(sprintf("LIV: prop_rhy=%.1f%%  r_median=%.2f\n", mean(rhy_liv)*100, r_liv_med))
cat(sprintf("CER: prop_rhy=%.1f%%\n", mean(rhy_cer)*100))
cat(sprintf("prop_DR (LIV vs CER): %.1f%%\n\n", prop_DR_mouse*100))

# Bootstrap design grid
boot_opts_mouse <- CircadianBootstrapOptions(
  design        = "active",
  design_vector = seq(0, 24, length.out = max(B_VALS_MOUSE) + 1)[seq_len(max(B_VALS_MOUSE))],
  B_values      = B_VALS_MOUSE,
  N_values      = N_GRID_MOUSE,
  nboot         = NBOOT,
  nsims_inner   = NSIMS_INNER,
  seed          = 42
)

bio_mouse <- estCircadianParam(
  data = mat_liv_s, times = ct_liv, period = 24,
  prop_DR = prop_DR_mouse, prop_DP = 0,
  min_rhythm_pval = RHYTHM_PVAL, verbose = FALSE
)

cat("Running bootstrap design grid (mouse)...\n")
s1_boot <- runBootstrapDesignGrid(
  pilot_data    = mat_liv_s,
  pilot_times   = ct_liv,
  boot.opts     = boot_opts_mouse,
  analysis.opts = opts_analysis,
  bio_diff.opts = bio_mouse,
  verbose       = TRUE
)
cat("\n--- Bootstrap Summary (Mouse DR, FDR 5%) ---\n")
summaryBootstrapDesignGrid(s1_boot, test_type = "DR", fdr_threshold = 0.05)

plotBootstrapDesignGrid(s1_boot, test_type = "DR", fdr_threshold = 0.05,
                        output_file = file.path(fig_dir, "s1_mouse_bootstrap.pdf"))
saveRDS(s1_boot, file.path(out_dir, "s1_mouse_boot_grid.rds"))
cat(sprintf("Figure: %s/figures/s1_mouse_bootstrap.pdf\n", out_dir))

# =====================================================================
# SECTION 2: BABOON — LUN vs CER
# =====================================================================

cat("\n====================================================================\n")
cat("SECTION 2: Baboon — LUN vs CER (bootstrap design grid)\n")
cat("====================================================================\n\n")

load(DATA_BABOON)   # loads: baboon_withTOD, gtex, mice (the older mice object)
bab_expr <- baboon_withTOD$baboon
bab_tod  <- baboon_withTOD$tod

# Load and transform tissues (raw CPM data.frame → log2(CPM+1) matrix)
prep_lun  <- prepCircadianData(bab_expr[["LUN"]], times = bab_tod[["LUN"]], input_type = "cpm")
prep_cerb <- prepCircadianData(bab_expr[["CER"]], times = bab_tod[["CER"]], input_type = "cpm")
mat_lun   <- prep_lun$data;  tod_lun <- prep_lun$times
mat_cer_b <- prep_cerb$data; tod_cer <- prep_cerb$times

# Filter expressed genes
mat_lun   <- mat_lun[rowSums(mat_lun > 0) >= 6, , drop=FALSE]
mat_cer_b <- mat_cer_b[rowSums(mat_cer_b > 0) >= 6, , drop=FALSE]

# Common genes, subsample
common_b <- intersect(rownames(mat_lun), rownames(mat_cer_b))
cat(sprintf("Common genes (LUN ∩ CER): %d\n", length(common_b)))
set.seed(2)
g_idx_b <- sample(common_b, min(NGENES, length(common_b)))

mat_lun_s   <- mat_lun[g_idx_b, , drop=FALSE]
mat_cer_bs  <- mat_cer_b[g_idx_b, , drop=FALSE]

# Estimate prop_DR
fit_lun <- fitCosinorAll(mat_lun_s,  times=tod_lun, period=24)
fit_cerb <- fitCosinorAll(mat_cer_bs, times=tod_cer, period=24)
rhy_lun  <- !is.na(fit_lun$pvalue)  & fit_lun$pvalue  < RHYTHM_PVAL
rhy_cerb <- !is.na(fit_cerb$pvalue) & fit_cerb$pvalue < RHYTHM_PVAL
dr_mask_b   <- xor(rhy_lun, rhy_cerb)
prop_DR_bab <- mean(dr_mask_b)
r_lun_med   <- median(as.numeric(fit_lun$A[rhy_lun]) / as.numeric(fit_lun$sigma[rhy_lun]), na.rm=TRUE)

cat(sprintf("LUN: prop_rhy=%.1f%%  r_median=%.2f\n", mean(rhy_lun)*100, r_lun_med))
cat(sprintf("CER: prop_rhy=%.1f%%\n", mean(rhy_cerb)*100))
cat(sprintf("prop_DR (LUN vs CER): %.1f%%\n\n", prop_DR_bab*100))

# Bootstrap design grid
bab_zt_full <- tod_lun   # ZT 0,2,...,22 (12 points)
boot_opts_bab <- CircadianBootstrapOptions(
  design        = "active",
  design_vector = bab_zt_full,
  B_values      = B_VALS_BABOON,
  N_values      = N_GRID_BABOON,
  nboot         = NBOOT,
  nsims_inner   = NSIMS_INNER,
  seed          = 42
)

bio_bab <- estCircadianParam(
  data = mat_lun_s, times = tod_lun, period = 24,
  prop_DR = prop_DR_bab, prop_DP = 0,
  min_rhythm_pval = RHYTHM_PVAL, verbose = FALSE
)

cat("Running bootstrap design grid (baboon)...\n")
s2_boot <- runBootstrapDesignGrid(
  pilot_data    = mat_lun_s,
  pilot_times   = tod_lun,
  boot.opts     = boot_opts_bab,
  analysis.opts = opts_analysis,
  bio_diff.opts = bio_bab,
  verbose       = TRUE
)
cat("\n--- Bootstrap Summary (Baboon DR, FDR 5%) ---\n")
summaryBootstrapDesignGrid(s2_boot, test_type = "DR", fdr_threshold = 0.05)

plotBootstrapDesignGrid(s2_boot, test_type = "DR", fdr_threshold = 0.05,
                        output_file = file.path(fig_dir, "s2_baboon_bootstrap.pdf"))
saveRDS(s2_boot, file.path(out_dir, "s2_baboon_boot_grid.rds"))
cat(sprintf("Figure: %s/figures/s2_baboon_bootstrap.pdf\n", out_dir))

# =====================================================================
# SECTION 3: SENEY MDD vs CONTROL (ACC, passive post-mortem)
# =====================================================================
# Passive design example: MDD vs unaffected control, anterior cingulate cortex (ACC).
# TOD = official post-mortem pronounced time — irregular, uncontrolled.
#
# SCIENTIFIC STORY — sex stratification:
#   Sex is a well-documented source of heterogeneity in circadian gene expression
#   in the ACC and in MDD biology (Seney et al. 2018). When male and female subjects
#   are pooled, two rhythm profiles with different phase alignments and amplitudes
#   are averaged into a single cosinor estimate. This sex-driven heterogeneity:
#     (1) Reduces apparent amplitude A — phase incoherence between sexes causes
#         partial cancellation of the rhythmic signal across subjects
#     (2) Inflates residual variance σ — cosinor residuals absorb between-sex variation
#   Together these reduce r = A/σ from ~0.82 (male-only pilot) to ~0.55 (combined),
#   a ~50% drop. Since power scales as r² × N, this requires ~3× more subjects.
#   Stratifying by sex restores within-sex phase coherence: n80 ~200 (male-only)
#   vs ~300 (combined). See 07_seney_sex_DR.R for the full sex-stratified comparison.
#
# SIGNAL (p<0.05, combined):  r~0.55, prop_DR~19%, n80 ~280-300
# SIGNAL (p<0.05, male-only): r~0.82, prop_DR~19%, n80 ~200-220
# HUMAN AGING ALTERNATIVE: see "SECTION 3 ALTERNATIVE" block below

cat("\n====================================================================\n")
cat("SECTION 3: Seney MDD vs Control ACC (bootstrap, passive design)\n")
cat("====================================================================\n\n")

suppressPackageStartupMessages(library(readxl))
meta_seney <- read_excel(DATA_SENEY_META)
tod_seney  <- read_excel(DATA_SENEY_TOD)
expr_seney <- as.matrix(read.csv(DATA_SENEY_EXPR, row.names = 1, check.names = FALSE))
# (prepCircadianData applied below after metadata join + ok_s filter)

# Match samples: strip trailing letter from colnames (e.g. "13219A" -> "13219")
col_ids_s  <- gsub("[A-Za-z]+$", "", colnames(expr_seney))
meta_idx_s <- match(col_ids_s, as.character(meta_seney$HU_NUM))
tod_idx_s  <- match(col_ids_s, as.character(tod_seney$HU_NUM))

# TOD = hour of day from official pronounced time
tod_hour_s <- as.numeric(format(tod_seney$OFFC_TIME[tod_idx_s], "%H")) +
              as.numeric(format(tod_seney$OFFC_TIME[tod_idx_s], "%M")) / 60

disease_s <- meta_seney$Disease[meta_idx_s]   # 1=Control, 2=MDD

# Filter to complete cases; use prepCircadianData for standard validation + coercion
ok_s <- !is.na(disease_s) & !is.na(tod_hour_s)
prep_seney <- prepCircadianData(expr_seney[, ok_s], times = tod_hour_s[ok_s], input_type = "log2")
expr_seney <- prep_seney$data
tod_hour_s <- prep_seney$times
disease_s  <- disease_s[ok_s]

ctrl_idx <- disease_s == 1
mdd_idx  <- disease_s == 2
tod_ctrl <- tod_hour_s[ctrl_idx]
tod_mdd  <- tod_hour_s[mdd_idx]

cat(sprintf("Control: n=%d  TOD range=[%.1f, %.1f]h\n", sum(ctrl_idx), min(tod_ctrl,na.rm=TRUE), max(tod_ctrl,na.rm=TRUE)))
cat(sprintf("MDD:     n=%d  TOD range=[%.1f, %.1f]h\n", sum(mdd_idx),  min(tod_mdd,na.rm=TRUE),  max(tod_mdd,na.rm=TRUE)))

# Subsample genes
set.seed(3)
g_idx_s  <- sample(nrow(expr_seney), min(NGENES, nrow(expr_seney)))
mat_ctrl <- expr_seney[g_idx_s, ctrl_idx, drop=FALSE]
mat_mdd  <- expr_seney[g_idx_s, mdd_idx,  drop=FALSE]

# Estimate parameters from both groups (p<0.05; estCircadianParamTwoGroup uses
# union rhythmic budget so prop_DR is never clipped by group-1 rhythmicity alone)
cat("Estimating parameters (two-group)...\n")
bio_seney <- estCircadianParamTwoGroup(
  data_1          = mat_ctrl,
  data_2          = mat_mdd,
  times_1         = tod_ctrl,
  times_2         = tod_mdd,
  period          = 24,
  min_rhythm_pval = RHYTHM_PVAL,
  verbose         = FALSE
)
prop_DR_seney <- bio_seney$prop_DR
r_ctrl_med    <- bio_seney$diagnostics$r1_snr["median"]

cat(sprintf("Control: prop_rhy=%.1f%%  r_median=%.2f\n",
            100*bio_seney$diagnostics$prop_rhythmic_1, r_ctrl_med))
cat(sprintf("MDD:     prop_rhy=%.1f%%\n", 100*bio_seney$diagnostics$prop_rhythmic_2))
cat(sprintf("prop_DR (union budget): %.1f%%\n\n", 100*prop_DR_seney))

# Bootstrap options (passive: B fixed at 4L — sweeping B has no meaning in passive)
boot_opts_seney <- CircadianBootstrapOptions(
  design        = "passive",
  design_vector = tod_ctrl,     # control TOD as reference distribution
  B_values      = 4L,
  N_values      = N_GRID_HUMAN,
  nboot         = NBOOT,
  nsims_inner   = NSIMS_INNER,
  seed          = 42
)

cat("Running bootstrap design grid (Seney MDD vs Control)...\n")
s3_boot <- runBootstrapDesignGrid(
  pilot_data    = mat_ctrl,
  pilot_times   = tod_ctrl,
  boot.opts     = boot_opts_seney,
  analysis.opts = opts_analysis,
  bio_diff.opts = bio_seney,
  verbose       = TRUE
)
cat("\n--- Bootstrap Summary (Seney DR, FDR 5%) ---\n")
summaryBootstrapDesignGrid(s3_boot, test_type = "DR", fdr_threshold = 0.05)

plotBootstrapDesignGrid(s3_boot, test_type = "DR", fdr_threshold = 0.05,
                        output_file = file.path(fig_dir, "s3_seney_bootstrap.pdf"))
saveRDS(s3_boot, file.path(out_dir, "s3_seney_boot_grid.rds"))
cat(sprintf("Figure: %s/figures/s3_seney_bootstrap.pdf\n", out_dir))

# =====================================================================
# SECTION 3 ALTERNATIVE: HUMAN AGING — young vs old BA11/BA47
# =====================================================================
# This block is commented out. To use human aging instead of Seney:
#   1. Comment out or remove Section 3 (Seney) above
#   2. Uncomment this block
#   3. Replace s3_boot / n80_s / pow_s references in Sections 5-6 with
#      the human aging equivalents below
#
# WHY SENEY IS PRIMARY:
#   Seney adds a sex-stratification dimension (07_seney_sex_DR.R) showing
#   that controlling a biological confounder reduces n80 from ~300 → ~200.
#   Human aging tells only the "passive is costly" story with no design lever.
#
# WHY KEEP AGING:
#   Well-known dataset (published), directly mirrors the 03_power_core.R
#   pipeline, and aging circadian disruption is a major research area.
#
# SIGNAL COMPARISON:
#   Human aging: r~0.50, prop_DR~14%, passive → n80 >300
#   Seney combined: r~0.55, prop_DR~19%, passive → n80 ~280-300
#   Both tell the same "passive is expensive" story; Seney is slightly richer.

# if (FALSE) {   # <-- change to if (TRUE) { to activate
#
# cat("\n====================================================================\n")
# cat("SECTION 3 (ALT): Human aging young vs old BA11/BA47\n")
# cat("====================================================================\n\n")
#
# dat_human <- readRDS(DATA_HUMAN)
# expr_human <- dat_human$expr
# pheno_human <- dat_human$pheno
#
# # Key columns: AgeGroup ("younger"/"older"), TOD.y (time of death, h)
# tod_col  <- if ("TOD.x" %in% colnames(pheno_human)) "TOD.x" else "TOD.y"
# age_col  <- if ("AgeGroup" %in% colnames(pheno_human)) "AgeGroup" else "age_group"
# valid_h  <- !is.na(pheno_human[[age_col]]) & !is.na(pheno_human[[tod_col]]) &
#             pheno_human[[age_col]] %in% c("younger", "older")
# pheno_h  <- pheno_human[valid_h, ]
# expr_h   <- expr_human[, valid_h]
# tod_h    <- pheno_h[[tod_col]]
# young_idx <- pheno_h[[age_col]] == "younger"   # n=62 pilot
# old_idx   <- pheno_h[[age_col]] == "older"     # n=74
#
# g_idx_h <- sample(nrow(expr_h), min(NGENES, nrow(expr_h)))
# mat_young <- expr_h[g_idx_h, young_idx, drop=FALSE]
# mat_old   <- expr_h[g_idx_h, old_idx,   drop=FALSE]
# tod_young <- tod_h[young_idx]
# tod_old   <- tod_h[old_idx]
#
# # Estimate parameters using both groups (p<0.10 works better for weak passive signal)
# bio_aging <- estCircadianParamTwoGroup(
#   mat_young, mat_old, tod_young, tod_old,
#   period=24, min_rhythm_pval=0.10, verbose=FALSE)
#
# cat(sprintf("Young rhythmic: %.1f%%  Old rhythmic: %.1f%%  prop_DR: %.1f%%  r: %.3f\n",
#   100*bio_aging$diagnostics$prop_rhythmic_1,
#   100*bio_aging$diagnostics$prop_rhythmic_2,
#   100*bio_aging$prop_DR,
#   bio_aging$diagnostics$r1_snr["median"]))
#
# boot_opts_aging <- CircadianBootstrapOptions(
#   design       = "passive",
#   design_vector = tod_young,
#   B_values     = 4L,          # passive: B fixed, sweeping has no meaning
#   N_values     = N_GRID_HUMAN,
#   nboot        = NBOOT,
#   nsims_inner  = NSIMS_INNER
# )
# s3_boot <- runBootstrapDesignGrid(
#   pilot_data    = mat_young,
#   pilot_times   = tod_young,
#   boot.opts     = boot_opts_aging,
#   bio_diff.opts = bio_aging,
#   test_type     = "DR"
# )
# summaryBootstrapDesignGrid(s3_boot, test_type="DR", fdr_threshold=0.05)
# plotBootstrapDesignGrid(s3_boot, test_type="DR", fdr_threshold=0.05,
#   output_file=file.path(fig_dir, "s3_aging_bootstrap.pdf"))
# saveRDS(s3_boot, file.path(out_dir, "s3_aging_boot_grid.rds"))
#
# }  # end if (FALSE)

# =====================================================================
# SECTION 4: MOUSE D1 vs D2 CELL TYPES
# =====================================================================
# Active design, m~8 replicates/ZT — demonstrates B vs m tradeoff
# with a well-replicated pilot; effect size r~0.66 (moderate, between baboon and human)

cat("\n====================================================================\n")
cat("SECTION 4: Mouse D1 vs D2 cell types (bootstrap design grid)\n")
cat("====================================================================\n\n")

# Raw counts (despite filename); prepCircadianData normalizes to log2(CPM+1)
# and aligns sample order via pheno$sample column
pheno_d1d2 <- read.csv(DATA_D1D2_PHENO, row.names = 1)
prep_d1d2  <- prepCircadianData(DATA_D1D2_EXPR, times = "time", input_type = "counts",
                                pheno = pheno_d1d2, sample_col = "sample")
log_d1d2   <- prep_d1d2$data

d1_samp <- pheno_d1d2$sample[pheno_d1d2$cell == "D1"]
d2_samp <- pheno_d1d2$sample[pheno_d1d2$cell == "D2"]
mat_d1  <- log_d1d2[, colnames(log_d1d2) %in% d1_samp, drop=FALSE]
mat_d2  <- log_d1d2[, colnames(log_d1d2) %in% d2_samp, drop=FALSE]
tod_d1  <- pheno_d1d2$time[match(colnames(mat_d1), pheno_d1d2$sample)]
tod_d2  <- pheno_d1d2$time[match(colnames(mat_d2), pheno_d1d2$sample)]

cat(sprintf("D1: n=%d  D2: n=%d  ZT=%s\n", ncol(mat_d1), ncol(mat_d2),
            paste(sort(unique(tod_d1)), collapse=",")))

# Expression filter + subsample
keep_d1 <- rowSums(mat_d1 > 1) >= floor(ncol(mat_d1)/2)
mat_d1  <- mat_d1[keep_d1, , drop=FALSE]
mat_d2  <- mat_d2[keep_d1, , drop=FALSE]
set.seed(4)
g_idx_d <- sample(nrow(mat_d1), min(NGENES, nrow(mat_d1)))
mat_d1_s <- mat_d1[g_idx_d, , drop=FALSE]
mat_d2_s <- mat_d2[g_idx_d, , drop=FALSE]

# Estimate prop_DR
fit_d1  <- fitCosinorAll(mat_d1_s, times=tod_d1, period=24)
fit_d2  <- fitCosinorAll(mat_d2_s, times=tod_d2, period=24)
rhy_d1  <- !is.na(fit_d1$pvalue) & fit_d1$pvalue < RHYTHM_PVAL
rhy_d2  <- !is.na(fit_d2$pvalue) & fit_d2$pvalue < RHYTHM_PVAL
dr_mask_d   <- xor(rhy_d1, rhy_d2)
prop_DR_d1d2 <- mean(dr_mask_d)
r_d1_med    <- median(as.numeric(fit_d1$A[rhy_d1]) / as.numeric(fit_d1$sigma[rhy_d1]), na.rm=TRUE)

cat(sprintf("D1: prop_rhy=%.1f%%  r_median=%.2f\n", mean(rhy_d1)*100, r_d1_med))
cat(sprintf("D2: prop_rhy=%.1f%%\n", mean(rhy_d2)*100))
cat(sprintf("prop_DR (D1 vs D2): %.1f%%\n\n", prop_DR_d1d2*100))

# Bootstrap design grid
d1_zt_full <- sort(unique(tod_d1))   # 2,6,10,14,18,22
boot_opts_d1d2 <- CircadianBootstrapOptions(
  design        = "active",
  design_vector = d1_zt_full,
  B_values      = B_VALS_D1D2,
  N_values      = N_GRID_D1D2,
  nboot         = NBOOT,
  nsims_inner   = NSIMS_INNER,
  seed          = 42
)

bio_d1d2 <- estCircadianParam(
  data            = mat_d1_s,
  times           = tod_d1,
  period          = 24,
  prop_DR         = prop_DR_d1d2,
  prop_DP         = 0,
  min_rhythm_pval = RHYTHM_PVAL,
  verbose         = FALSE
)

cat("Running bootstrap design grid (D1D2)...\n")
s4_boot <- runBootstrapDesignGrid(
  pilot_data    = mat_d1_s,
  pilot_times   = tod_d1,
  boot.opts     = boot_opts_d1d2,
  analysis.opts = opts_analysis,
  bio_diff.opts = bio_d1d2,
  verbose       = TRUE
)
cat("\n--- Bootstrap Summary (D1D2 DR, FDR 5%) ---\n")
summaryBootstrapDesignGrid(s4_boot, test_type = "DR", fdr_threshold = 0.05)

plotBootstrapDesignGrid(s4_boot, test_type = "DR", fdr_threshold = 0.05,
                        output_file = file.path(fig_dir, "s4_d1d2_bootstrap.pdf"))
saveRDS(s4_boot, file.path(out_dir, "s4_d1d2_boot_grid.rds"))
cat(sprintf("Figure: %s/figures/s4_d1d2_bootstrap.pdf\n", out_dir))

# =====================================================================
# SECTION 5: FOUR-WAY COMPARISON FIGURE
# =====================================================================
# Panel A: Power curves (bootstrap median + 95% CI) vs N per group
# Panel B: n80 ± CI horizontal bar chart

cat("\n====================================================================\n")
cat("SECTION 5: Four-way comparison figure\n")
cat("====================================================================\n\n")

# Helper: extract marginal power (averaged over B) from bootstrap result
.extract_power_margin <- function(boot_res, fdr_thr = 0.05) {
  fdr_idx <- which.min(abs(boot_res$fdr_thresholds - fdr_thr))
  if (length(fdr_idx) == 0) fdr_idx <- 1L
  n_N <- length(boot_res$N_values)
  n_B <- length(boot_res$B_values)
  pm <- boot_res$power_mean[, , fdr_idx, drop=FALSE]
  lo <- boot_res$power_ci_lo[, , fdr_idx, drop=FALSE]
  hi <- boot_res$power_ci_hi[, , fdr_idx, drop=FALSE]
  data.frame(
    N          = boot_res$N_values,
    power_mean = rowMeans(matrix(pm, nrow=n_N, ncol=n_B)),
    ci_lo      = rowMeans(matrix(lo, nrow=n_N, ncol=n_B)),
    ci_hi      = rowMeans(matrix(hi, nrow=n_N, ncol=n_B))
  )
}

.n80 <- function(N_vals, power_vec, ci_lo_vec, ci_hi_vec) {
  list(med = N_vals[which(power_vec  >= 0.80)[1]],
       lo  = N_vals[which(ci_hi_vec >= 0.80)[1]],
       hi  = N_vals[which(ci_lo_vec >= 0.80)[1]])
}

pow_m <- .extract_power_margin(s1_boot)
pow_b <- .extract_power_margin(s2_boot)
pow_s <- .extract_power_margin(s3_boot)
pow_d <- .extract_power_margin(s4_boot)

n80_m <- .n80(pow_m$N, pow_m$power_mean, pow_m$ci_lo, pow_m$ci_hi)
n80_b <- .n80(pow_b$N, pow_b$power_mean, pow_b$ci_lo, pow_b$ci_hi)
n80_s <- .n80(pow_s$N, pow_s$power_mean, pow_s$ci_lo, pow_s$ci_hi)
n80_d <- .n80(pow_d$N, pow_d$power_mean, pow_d$ci_lo, pow_d$ci_hi)

.fmt <- function(x) ifelse(is.na(x), ">max", as.character(x))
cat("n80 summary:\n")
cat(sprintf("  GSE54651 LIV vs CER:    median=%s  CI=[%s, %s]\n", .fmt(n80_m$med), .fmt(n80_m$lo), .fmt(n80_m$hi)))
cat(sprintf("  Baboon   LUN vs CER:    median=%s  CI=[%s, %s]\n", .fmt(n80_b$med), .fmt(n80_b$lo), .fmt(n80_b$hi)))
cat(sprintf("  Seney    MDD vs Ctrl:   median=%s  CI=[%s, %s]\n", .fmt(n80_s$med), .fmt(n80_s$lo), .fmt(n80_s$hi)))
cat(sprintf("  Mouse    D1 vs D2:      median=%s  CI=[%s, %s]\n", .fmt(n80_d$med), .fmt(n80_d$lo), .fmt(n80_d$hi)))

cols <- c(gse="steelblue", baboon="darkorange", seney="firebrick", d1d2="forestgreen")

.draw_power_ribbon <- function(pdat, col) {
  polygon(c(pdat$N, rev(pdat$N)), c(pdat$ci_lo, rev(pdat$ci_hi)),
          col=adjustcolor(col, 0.15), border=NA)
  lines(pdat$N, pdat$power_mean, col=col, lwd=2)
  points(pdat$N, pdat$power_mean, col=col, pch=16, cex=0.7)
}

fig_compare <- file.path(fig_dir, "s5_four_way_comparison.pdf")
pdf(fig_compare, width=14, height=5)
par(mfrow=c(1,2), mar=c(4,4,3,1))

# Panel A: power curves (log x-axis to handle wide N range)
all_N <- sort(unique(c(pow_m$N, pow_b$N, pow_s$N, pow_d$N)))
plot(NA, xlim=range(all_N), ylim=c(0,1),
     xlab="N per group", ylab="Power",
     main="DR Power: Bootstrap median + 95% CI\n(4 datasets)")
abline(h=0.80, lty=2, col="gray40")
abline(h=c(0.2,0.4,0.6), lty=3, col="gray85")
.draw_power_ribbon(pow_m, cols["gse"])
.draw_power_ribbon(pow_b, cols["baboon"])
.draw_power_ribbon(pow_d, cols["d1d2"])
.draw_power_ribbon(pow_s, cols["seney"])
legend("bottomright",
       legend=c("GSE54651 LIV vs CER (active, m=1, r~2.9)",
                "Baboon LUN vs CER   (active, m=1, r~1.7)",
                "Mouse D1 vs D2      (active, m~8, r~0.66)",
                "Seney MDD vs Ctrl   (passive, r~0.55)"),
       col=cols, lwd=2, bty="n", cex=0.72)

# Panel B: n80 bar chart
datasets <- c("GSE54651\nLIV vs CER", "Baboon\nLUN vs CER", "Mouse\nD1 vs D2", "Seney\nMDD vs Ctrl")
N_maxes  <- c(max(N_GRID_MOUSE), max(N_GRID_BABOON), max(N_GRID_D1D2), max(N_GRID_HUMAN))
n80_meds <- mapply(function(n80, nm) ifelse(is.na(n80), nm+10, n80),
                   list(n80_m$med, n80_b$med, n80_d$med, n80_s$med), as.list(N_maxes))
n80_los  <- mapply(function(n80, nm) ifelse(is.na(n80), nm+10, n80),
                   list(n80_m$lo,  n80_b$lo,  n80_d$lo,  n80_s$lo),  as.list(N_maxes))
n80_his  <- mapply(function(n80, nm) ifelse(is.na(n80), nm+10, n80),
                   list(n80_m$hi,  n80_b$hi,  n80_d$hi,  n80_s$hi),  as.list(N_maxes))

bp <- barplot(unlist(n80_meds), names.arg=datasets, col=unname(cols), horiz=FALSE,
              las=1, ylim=c(0, max(unlist(n80_his), na.rm=TRUE)*1.2),
              ylab="n80 (N per group for 80% power)",
              main="n80 Recommendation (bootstrap median +- CI)",
              cex.names=0.75)
arrows(bp, unlist(n80_los), bp, unlist(n80_his),
       angle=90, code=3, length=0.07, col="gray30")
text(bp, unlist(n80_meds) + max(unlist(n80_his),na.rm=TRUE)*0.04,
     ifelse(unlist(n80_meds) > N_maxes, ">max(N)", as.character(round(unlist(n80_meds)))),
     cex=0.75, col="black")
dev.off()
cat(sprintf("Figure: %s\n", fig_compare))

# =====================================================================
# SECTION 5: B vs m TRADEOFF (active designs: mouse and baboon)
# =====================================================================
# For each active dataset, show power vs N stratified by B (time-point count).
# This directly answers Q3: for fixed N, is higher B (more coverage) better?

cat("\n====================================================================\n")
cat("SECTION 6: B vs m tradeoff — active designs (GSE54651, baboon, D1D2)\n")
cat("====================================================================\n\n")

.plot_bm_tradeoff <- function(boot_res, title, fdr_thr=0.05, col_palette=NULL) {
  fdr_idx <- which.min(abs(boot_res$fdr_thresholds - fdr_thr))
  if (length(fdr_idx) == 0) fdr_idx <- 1L

  N_vals <- boot_res$N_values
  B_vals <- boot_res$B_values
  n_B    <- length(B_vals)
  n_N    <- length(N_vals)

  if (is.null(col_palette)) col_palette <- rainbow(n_B, s=0.7, v=0.85)

  # Extract power per (N, B) — matrix [n_N x n_B]
  pm <- matrix(boot_res$power_mean[, , fdr_idx], nrow=n_N, ncol=n_B)
  lo <- matrix(boot_res$power_ci_lo[, , fdr_idx], nrow=n_N, ncol=n_B)
  hi <- matrix(boot_res$power_ci_hi[, , fdr_idx], nrow=n_N, ncol=n_B)

  plot(NA, xlim=range(N_vals), ylim=c(0,1),
       xlab="N per group", ylab="Power",
       main=title)
  abline(h=0.80, lty=2, col="gray40")
  abline(h=c(0.2,0.4,0.6), lty=3, col="gray85")

  for (b in seq_len(n_B)) {
    col_b <- col_palette[b]
    polygon(c(N_vals, rev(N_vals)),
            c(lo[, b], rev(hi[, b])),
            col=adjustcolor(col_b, 0.15), border=NA)
    lines(N_vals, pm[, b], col=col_b, lwd=2, lty=b)
    points(N_vals, pm[, b], col=col_b, pch=15+b, cex=0.8)
  }
  legend("bottomright",
         legend=sprintf("B=%d (%d repl/ZT)", B_vals, pmax(1L, round(min(N_vals)/B_vals))),
         col=col_palette, lwd=2, lty=seq_len(n_B), pch=15+seq_len(n_B),
         bty="n", cex=0.75)
}

fig_bm <- file.path(fig_dir, "s6_bm_tradeoff.pdf")
pdf(fig_bm, width=18, height=5)
par(mfrow=c(1,3), mar=c(4,4,3,1))

.plot_bm_tradeoff(s1_boot,
  title=sprintf("GSE54651 LIV vs CER\nB vs m tradeoff (prop_DR=%.0f%%)", prop_DR_mouse*100))
.plot_bm_tradeoff(s2_boot,
  title=sprintf("Baboon LUN vs CER\nB vs m tradeoff (prop_DR=%.0f%%)", prop_DR_bab*100))
.plot_bm_tradeoff(s4_boot,
  title=sprintf("Mouse D1 vs D2\nB vs m tradeoff (prop_DR=%.0f%%)", prop_DR_d1d2*100))

dev.off()
cat(sprintf("Figure: %s\n", fig_bm))

# =====================================================================
# WRAP-UP
# =====================================================================

cat("\n====================================================================\n")
cat("05_multi_dataset_DR COMPLETE\n")
cat("====================================================================\n\n")

cat("--- Summary ---\n")
cat(sprintf("GSE54651 LIV vs CER:  prop_DR=%.1f%%  r_median(LIV)=%.2f    n80=%s\n",
            prop_DR_mouse*100, r_liv_med,   .fmt(n80_m$med)))
cat(sprintf("Baboon   LUN vs CER:  prop_DR=%.1f%%  r_median(LUN)=%.2f    n80=%s\n",
            prop_DR_bab*100,   r_lun_med,   .fmt(n80_b$med)))
cat(sprintf("Seney    MDD vs Ctrl: prop_DR=%.1f%%  r_median(Ctrl)=%.2f   n80=%s\n",
            prop_DR_seney*100, r_ctrl_med,  .fmt(n80_s$med)))
cat(sprintf("Mouse    D1 vs D2:    prop_DR=%.1f%%  r_median(D1)=%.2f     n80=%s\n",
            prop_DR_d1d2*100,  r_d1_med,    .fmt(n80_d$med)))

cat(sprintf("\nOutput directory: %s/\n", out_dir))
cat(sprintf("Figures: %s/\n", fig_dir))
