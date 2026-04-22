#' =======================================================================
#' 08g_bootstrap_diff.R — Differential Bootstrap: Plug-in vs Bootstrap CI
#' =======================================================================
#'
#' Differential analogue of 08f_bootstrap_sc.R. Shows that the plug-in
#' differential power curve overpredicts when the pilot is small, and that
#' bootstrap CIs widen as pilot size shrinks — for DR endpoint.
#'
#' Three panels ordered by pilot size:
#'   A. Baboon LUN vs CER      n_pilot = 12   active  B=4   (small pilot → wide CI)
#'   B. Mouse D1 vs D2         n_pilot = 45   active  B=4   (moderate pilot)
#'   C. Seney CTL vs MDD ACC   n_pilot = 60   passive B=4   (larger pilot → narrow CI)
#'
#' USAGE:
#'   Rscript examples/publication/08g_bootstrap_diff.R
#'   SMOKE_TEST=true Rscript examples/publication/08g_bootstrap_diff.R

SMOKE_TEST <- identical(Sys.getenv("SMOKE_TEST"), "true")

NGENES      <- if (SMOKE_TEST) 300L  else 5000L
NBOOT       <- if (SMOKE_TEST) 5L    else 50L
NSIMS       <- if (SMOKE_TEST) 5L    else 30L
NSIMS_INNER <- if (SMOKE_TEST) 5L    else 25L
N_CORES     <- as.integer(Sys.getenv("MC_CORES", unset = "6"))
RHYTHM_PVAL <- 0.01  # alpha_pilot per paper (SCP.tex §2.1); must match bootstrap's fitCosinorAll threshold
B_VAL       <- 4L
GLOBAL_SEED <- 2025L

cat(sprintf("Mode        : %s\n", if (SMOKE_TEST) "SMOKE" else "PRODUCTION"))
cat(sprintf("NGENES      : %d\n", NGENES))
cat(sprintf("NBOOT       : %d\n", NBOOT))
cat(sprintf("NSIMS       : %d\n", NSIMS))
cat(sprintf("NSIMS_INNER : %d\n", NSIMS_INNER))
cat(sprintf("MC_CORES    : %d\n\n", N_CORES))

old_wd <- setwd("code"); source("setup.R"); setwd(old_wd)
suppressPackageStartupMessages(library(readxl))

out_dir <- "output/bootstrap_diff"
dir.create(file.path(out_dir, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

analysis <- CircadianAnalysisOptions(
  alpha           = 0.05,
  p.adjust.method = "BH",
  fdr_thresholds  = 0.05
)

# -----------------------------------------------------------------------
# Helper: run plug-in + bootstrap for one two-group dataset and save
# -----------------------------------------------------------------------
.run_diff_comparison <- function(mat_1, tod_1, mat_2, tod_2,
                                  bio_diff, N_grid, B_val,
                                  design_type, label, out_prefix) {
  cat(sprintf("\n--- %s ---\n", label))
  cat(sprintf("    n1=%d  n2=%d  ngenes=%d  B=%d  design=%s\n",
              ncol(mat_1), ncol(mat_2), nrow(mat_1), B_val, design_type))

  design_vec <- if (design_type == "active") {
    seq(0, 24 * (1 - 1/B_val), length.out = B_val)
  } else {
    tod_1
  }

  # Plug-in
  design_opts <- CircadianDesignOptions(
    sample_sizes = N_grid,
    nsims        = NSIMS,
    design       = design_type,
    cts          = design_vec,
    B_values     = B_val
  )
  cat("  Running plug-in (runDifferentialPower)...\n")
  set.seed(GLOBAL_SEED)
  pi_result <- tryCatch(
    runDifferentialPower(bio_diff, design_opts, analysis,
                         methods    = "DCP",
                         test_types = "DR",
                         alpha2     = 0,
                         mc.cores   = N_CORES,
                         plot       = FALSE,
                         verbose    = FALSE),
    error = function(e) { warning(sprintf("Plug-in failed: %s", e$message)); NULL }
  )

  # Bootstrap
  boot_opts <- CircadianBootstrapOptions(
    design_vector = design_vec,
    B_values      = B_val,
    N_values      = N_grid,
    nboot         = NBOOT,
    nsims_inner   = NSIMS_INNER,
    design        = design_type,
    seed          = GLOBAL_SEED
  )
  cat(sprintf("  Running bootstrap (%d draws, mc.cores=%d)...\n", NBOOT, N_CORES))
  set.seed(GLOBAL_SEED)
  bt_result <- tryCatch(
    runBootstrapDesignGrid(
      pilot_data      = mat_1,
      pilot_times     = tod_1,
      pilot_data_2    = mat_2,
      pilot_times_2   = tod_2,
      bio_diff.opts   = bio_diff,
      boot.opts       = boot_opts,
      analysis.opts   = analysis,
      mode            = "differential",
      min_rhythm_pval = RHYTHM_PVAL,
      mc.cores        = N_CORES,
      verbose         = FALSE
    ),
    error = function(e) { warning(sprintf("Bootstrap failed: %s", e$message)); NULL }
  )

  if (is.null(pi_result) || is.null(bt_result)) {
    cat("  FAILED — skipping\n"); return(NULL)
  }

  # DR is test index 1 in the result arrays
  dr_idx <- 1L

  # Summary table
  pi_pwr <- pi_result$power_df$power[pi_result$power_df$test_type == "DR"]
  cat(sprintf("  %-5s  plug-in  bootstrap  gap(pp)  CI-width\n", "N"))
  for (i in seq_along(N_grid)) {
    tp <- if (i <= length(pi_pwr)) pi_pwr[i] else NA_real_
    bp <- bt_result$power_mean[i, 1, dr_idx]
    lo <- bt_result$power_ci_lo[i, 1, dr_idx]
    hi <- bt_result$power_ci_hi[i, 1, dr_idx]
    cat(sprintf("  N=%-4d  %.1f%%     %.1f%%      %+.1fpp   %.1fpp\n",
        N_grid[i], 100*tp, 100*bp, 100*(tp - bp), 100*(hi - lo)))
  }

  out <- list(label = label, n_pilot = ncol(mat_1), N_grid = N_grid,
              pi_result = pi_result, bt_result = bt_result,
              pi_pwr_dr = pi_pwr, dr_idx = dr_idx)
  saveRDS(out, paste0(out_prefix, ".rds"))
  cat(sprintf("  Saved: %s.rds\n", basename(out_prefix)))
  out
}


# =======================================================================
# PANEL A: Baboon LUN vs CER (n=12, active, small pilot)
# =======================================================================
cat("================================================================\n")
cat("PANEL A: Baboon LUN vs CER (n=12, active)\n")
cat("================================================================\n")

N_GRID_BAB <- if (SMOKE_TEST) c(12L, 24L, 36L) else c(12L, 24L, 36L, 48L, 60L, 72L, 96L)

load("data/CAMO_PRC_hmb.RData")
prep_lun  <- prepCircadianData(baboon_withTOD$baboon[["LUN"]],
                                times = baboon_withTOD$tod[["LUN"]], input_type = "cpm")
prep_cer  <- prepCircadianData(baboon_withTOD$baboon[["CER"]],
                                times = baboon_withTOD$tod[["CER"]], input_type = "cpm")
mat_lun   <- prep_lun$data[rowSums(prep_lun$data > 0) >= 6, , drop = FALSE]
mat_cer   <- prep_cer$data[rowSums(prep_cer$data > 0) >= 6, , drop = FALSE]
tod_lun   <- prep_lun$times
tod_cer   <- prep_cer$times
common_b  <- intersect(rownames(mat_lun), rownames(mat_cer))
set.seed(GLOBAL_SEED)
g_bab     <- sample(common_b, min(NGENES, length(common_b)))
mat_lun   <- mat_lun[g_bab, , drop = FALSE]
mat_cer   <- mat_cer[g_bab, , drop = FALSE]
rm(baboon_withTOD, gtex, mice, prep_lun, prep_cer)

bio_bab <- estCircadianParamTwoGroup(
  data_1 = mat_lun, data_2 = mat_cer,
  times_1 = tod_lun, times_2 = tod_cer,
  period = 24, min_rhythm_pval = RHYTHM_PVAL, verbose = FALSE
)
bio_bab$ngenes <- NGENES
cat(sprintf("Baboon: prop_DR=%.1f%%\n", 100 * bio_bab$prop_DR))

pA <- .run_diff_comparison(mat_lun, tod_lun, mat_cer, tod_cer,
                            bio_bab, N_GRID_BAB, B_VAL,
                            design_type = "active",
                            label = sprintf("Baboon LUN/CER (n=%d)", ncol(mat_lun)),
                            out_prefix = file.path(out_dir, "results", "panelA_baboon"))


# =======================================================================
# PANEL B: Mouse D1 vs D2 (n=45, active, moderate pilot)
# =======================================================================
cat("\n================================================================\n")
cat("PANEL B: Mouse D1 vs D2 (n=45, active)\n")
cat("================================================================\n")

N_GRID_D1 <- if (SMOKE_TEST) c(40L, 80L, 120L) else c(40L, 60L, 80L, 120L, 160L, 200L, 250L, 300L)

pheno   <- read.csv("data/mouse_clinicalinfo_03082021_rmOutliers.csv", row.names = 1)
prep_d1d2 <- prepCircadianData("data/mouse_D1D2_logCPMfiltered_counts.csv",
                                times = "time", input_type = "counts",
                                pheno = pheno, sample_col = "sample")
log_d1d2  <- prep_d1d2$data
d1_samp   <- pheno$sample[pheno$cell == "D1"]
d2_samp   <- pheno$sample[pheno$cell == "D2"]
mat_d1    <- log_d1d2[, colnames(log_d1d2) %in% d1_samp, drop = FALSE]
mat_d2    <- log_d1d2[, colnames(log_d1d2) %in% d2_samp, drop = FALSE]
tod_d1    <- pheno$time[match(colnames(mat_d1), pheno$sample)]
tod_d2    <- pheno$time[match(colnames(mat_d2), pheno$sample)]
keep      <- rowSums(mat_d1 > 1) >= floor(ncol(mat_d1) / 2)
mat_d1    <- mat_d1[keep, , drop = FALSE]
mat_d2    <- mat_d2[keep, , drop = FALSE]
set.seed(GLOBAL_SEED)
g_d1      <- sample(nrow(mat_d1), min(NGENES, nrow(mat_d1)))
mat_d1    <- mat_d1[g_d1, , drop = FALSE]
mat_d2    <- mat_d2[g_d1, , drop = FALSE]
rm(pheno, prep_d1d2, log_d1d2)

bio_d1d2 <- estCircadianParamTwoGroup(
  data_1 = mat_d1, data_2 = mat_d2,
  times_1 = tod_d1, times_2 = tod_d2,
  period = 24, min_rhythm_pval = RHYTHM_PVAL, verbose = FALSE
)
bio_d1d2$ngenes <- NGENES
cat(sprintf("Mouse D1/D2: prop_DR=%.1f%%\n", 100 * bio_d1d2$prop_DR))

pB <- .run_diff_comparison(mat_d1, tod_d1, mat_d2, tod_d2,
                            bio_d1d2, N_GRID_D1, B_VAL,
                            design_type = "active",
                            label = sprintf("Mouse D1/D2 (n=%d)", ncol(mat_d1)),
                            out_prefix = file.path(out_dir, "results", "panelB_d1d2"))


# =======================================================================
# PANEL C: Seney CTL vs MDD ACC (n=60, passive, larger pilot)
# =======================================================================
cat("\n================================================================\n")
cat("PANEL C: Seney CTL vs MDD ACC (n=60, passive)\n")
cat("================================================================\n")

N_GRID_SEN <- if (SMOKE_TEST) c(40L, 80L, 120L) else c(40L, 80L, 120L, 160L, 200L, 250L, 300L)

meta_s   <- read_excel("data/MD5_MetaData_1-15-25.xlsx")
tod_s    <- read_excel("data/TOD.xlsx")
expr_raw <- as.matrix(read.csv("data/ACC_RNA_filtered_normalized.csv",
                                row.names = 1, check.names = FALSE))
col_ids  <- gsub("[A-Za-z]+$", "", colnames(expr_raw))
meta_idx <- match(col_ids, as.character(meta_s$HU_NUM))
tod_idx  <- match(col_ids, as.character(tod_s$HU_NUM))
tod_hour <- as.numeric(format(tod_s$OFFC_TIME[tod_idx], "%H")) +
            as.numeric(format(tod_s$OFFC_TIME[tod_idx], "%M")) / 60
disease  <- meta_s$Disease[meta_idx]
ok       <- !is.na(disease) & !is.na(tod_hour)
prep_sen <- prepCircadianData(expr_raw[, ok], times = tod_hour[ok], input_type = "log2")
ctrl_idx <- disease[ok] == 1
mdd_idx  <- disease[ok] == 2
mat_ctrl <- prep_sen$data[, ctrl_idx, drop = FALSE]
mat_mdd  <- prep_sen$data[, mdd_idx,  drop = FALSE]
tod_ctrl <- prep_sen$times[ctrl_idx]
tod_mdd  <- prep_sen$times[mdd_idx]
set.seed(GLOBAL_SEED)
g_sen    <- sample(nrow(mat_ctrl), min(NGENES, nrow(mat_ctrl)))
mat_ctrl <- mat_ctrl[g_sen, , drop = FALSE]
mat_mdd  <- mat_mdd[g_sen,  , drop = FALSE]
rm(meta_s, tod_s, expr_raw, prep_sen)

bio_sen <- estCircadianParamTwoGroup(
  data_1 = mat_ctrl, data_2 = mat_mdd,
  times_1 = tod_ctrl, times_2 = tod_mdd,
  period = 24, min_rhythm_pval = RHYTHM_PVAL, verbose = FALSE
)
bio_sen$ngenes <- NGENES
cat(sprintf("Seney CTL/MDD: prop_DR=%.1f%%\n", 100 * bio_sen$prop_DR))

pC <- .run_diff_comparison(mat_ctrl, tod_ctrl, mat_mdd, tod_mdd,
                            bio_sen, N_GRID_SEN, B_VAL,
                            design_type = "passive",
                            label = sprintf("Seney CTL/MDD (n=%d)", ncol(mat_ctrl)),
                            out_prefix = file.path(out_dir, "results", "panelC_seney"))


# =======================================================================
# Combined figure
# =======================================================================
cat("\n=== Generating combined figure ===\n")
panels <- Filter(Negate(is.null), list(pA, pB, pC))

if (length(panels) > 0) {
  pdf(file.path(out_dir, "figures", "fig_bootstrap_diff.pdf"),
      width = 14, height = 5)
  par(mfrow = c(1, length(panels)), mar = c(4, 4, 3, 1), las = 1)

  panel_letters <- LETTERS[seq_along(panels)]
  for (pi in seq_along(panels)) {
    p      <- panels[[pi]]
    N_grid <- p$N_grid
    ts_pwr <- p$pi_pwr_dr
    dr_idx <- p$dr_idx
    bt_mn  <- p$bt_result$power_mean[, 1, dr_idx]
    bt_lo  <- p$bt_result$power_ci_lo[, 1, dr_idx]
    bt_hi  <- p$bt_result$power_ci_hi[, 1, dr_idx]

    plot(N_grid, ts_pwr, type = "l", lwd = 2, col = "steelblue",
         ylim = c(0, 1), xlab = "N per group", ylab = "Power (DR)",
         main = sprintf("(%s) %s", panel_letters[pi], p$label))
    polygon(c(N_grid, rev(N_grid)),
            c(bt_lo, rev(bt_hi)),
            col = adjustcolor("tomato", 0.25), border = NA)
    lines(N_grid, bt_mn, lwd = 2, col = "tomato")
    abline(h = 0.80, lty = 2, col = "grey50")
    if (pi == 1) {
      legend("bottomright", bty = "n", lwd = 2,
             col = c("steelblue", "tomato"),
             legend = c("Plug-in", "Bootstrap mean + 95% CI"))
    }
  }
  dev.off()
  cat(sprintf("Saved: %s/figures/fig_bootstrap_diff.pdf\n", out_dir))
}

saveRDS(panels, file.path(out_dir, "results", "all_panels.rds"))
cat("\n=== Done ===\n")
