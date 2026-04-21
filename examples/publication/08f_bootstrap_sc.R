#' =======================================================================
#' 08f_bootstrap_sc.R — Single-Cohort Bootstrap: Plug-in vs Bootstrap CI
#' =======================================================================
#'
#' Shows that the plug-in (two-stage) single-cohort power curve overpredicts
#' when the pilot is small, and that bootstrap CIs widen as pilot size shrinks.
#'
#' Three panels ordered by pilot size:
#'   A. Baboon LUN          n_pilot = 12   active  B=4   r~1.72  (small pilot → wide CI)
#'   B. Mouse D1 striatum   n_pilot = 45   passive B=4   r~0.65  (moderate pilot)
#'   C. Seney CTL ACC       n_pilot = 60   passive B=4   r~0.79  (larger pilot → narrow CI)
#'
#' Story: overprediction gap (plug-in minus bootstrap mean) and CI width both
#' shrink as n_pilot grows, validating that the bootstrap correction matters
#' most when pilot data are scarce.
#'
#' USAGE:
#'   Rscript examples/publication/08f_bootstrap_sc.R
#'   SMOKE_TEST=true Rscript examples/publication/08f_bootstrap_sc.R

SMOKE_TEST <- identical(Sys.getenv("SMOKE_TEST"), "true")

NGENES      <- if (SMOKE_TEST) 300L  else 3000L
NBOOT       <- if (SMOKE_TEST) 5L    else 50L
NSIMS       <- if (SMOKE_TEST) 5L    else 30L
NSIMS_INNER <- if (SMOKE_TEST) 5L    else 20L
N_CORES     <- as.integer(Sys.getenv("MC_CORES", unset = "6"))
RHYTHM_PVAL <- 0.05
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

out_dir <- "output/bootstrap_sc"
dir.create(file.path(out_dir, "results"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"),  recursive = TRUE, showWarnings = FALSE)

analysis <- CircadianAnalysisOptions(
  alpha           = 0.05,
  p.adjust.method = "BH",
  fdr_thresholds  = 0.05
)

# -----------------------------------------------------------------------
# Helper: run plug-in + bootstrap for one dataset and save
# -----------------------------------------------------------------------
.run_sc_comparison <- function(mat, tod, bio_sc, N_grid, B_val,
                                label, out_prefix) {
  cat(sprintf("\n--- %s ---\n", label))
  cat(sprintf("    n_pilot=%d  ngenes=%d  B=%d\n",
              ncol(mat), nrow(mat), B_val))

  design_vec <- seq(0, 24 * (1 - 1/B_val), length.out = B_val)

  # Plug-in
  design_ts <- CircadianDesignOptions(
    sample_sizes = N_grid, nsims = NSIMS, design = "active",
    cts = design_vec, B_values = B_val
  )
  cat("  Running plug-in (runSingleCohortPower)...\n")
  set.seed(GLOBAL_SEED)
  ts_result <- tryCatch(
    runSingleCohortPower(bio_sc, design_ts, analysis, methods = "DCP",
                         alpha2 = 0, mc.cores = N_CORES,
                         plot = FALSE, verbose = FALSE),
    error = function(e) { warning(sprintf("Plug-in failed: %s", e$message)); NULL }
  )

  # Bootstrap
  boot_opts <- CircadianBootstrapOptions(
    design_vector = design_vec, B_values = B_val, N_values = N_grid,
    nboot = NBOOT, nsims_inner = NSIMS_INNER, design = "active", seed = GLOBAL_SEED
  )
  cat(sprintf("  Running bootstrap (%d draws, mc.cores=%d)...\n", NBOOT, N_CORES))
  set.seed(GLOBAL_SEED)
  bt_result <- tryCatch(
    runBootstrapDesignGrid(
      pilot_data    = mat,
      pilot_times   = tod,
      bio_diff.opts = bio_sc,
      boot.opts     = boot_opts,
      analysis.opts = analysis,
      mode          = "single",
      mc.cores      = N_CORES,
      verbose       = FALSE
    ),
    error = function(e) { warning(sprintf("Bootstrap failed: %s", e$message)); NULL }
  )

  if (is.null(ts_result) || is.null(bt_result)) {
    cat("  FAILED — skipping\n"); return(NULL)
  }

  # Summary table
  cat(sprintf("  %-5s  plug-in  bootstrap  gap(pp)  CI-width\n", "N"))
  for (i in seq_along(N_grid)) {
    tp <- ts_result$power_df$power[i]
    bp <- bt_result$power_mean[i, 1, 1]
    lo <- bt_result$power_ci_lo[i, 1, 1]
    hi <- bt_result$power_ci_hi[i, 1, 1]
    cat(sprintf("  N=%-4d  %.1f%%     %.1f%%      %+.1fpp   %.1fpp\n",
        N_grid[i], 100*tp, 100*bp, 100*(tp-bp), 100*(hi-lo)))
  }

  out <- list(label = label, n_pilot = ncol(mat), N_grid = N_grid,
              ts_result = ts_result, bt_result = bt_result)
  saveRDS(out, paste0(out_prefix, ".rds"))
  cat(sprintf("  Saved: %s.rds\n", basename(out_prefix)))
  out
}


# =======================================================================
# PANEL A: Baboon LUN (n=12, active pilot, r~1.72)
# =======================================================================
cat("================================================================\n")
cat("PANEL A: Baboon LUN (n=12, active, r~1.72)\n")
cat("================================================================\n")

N_GRID_BAB <- if (SMOKE_TEST) c(12L, 24L, 36L) else c(12L, 24L, 36L, 48L, 60L, 72L, 96L)
# All divisible by B_VAL=4: 12%4=0, 24%4=0, ... ✓

load("data/CAMO_PRC_hmb.RData")
prep_lun <- prepCircadianData(baboon_withTOD$baboon[["LUN"]],
                               times = baboon_withTOD$tod[["LUN"]],
                               input_type = "cpm")
mat_bab  <- prep_lun$data[rowSums(prep_lun$data > 0) >= 6, , drop = FALSE]
tod_bab  <- prep_lun$times
set.seed(GLOBAL_SEED)
g_bab    <- sample(nrow(mat_bab), min(NGENES, nrow(mat_bab)))
mat_bab  <- mat_bab[g_bab, , drop = FALSE]
rm(baboon_withTOD, gtex, mice, prep_lun)

bio_bab <- estCircadianParam(mat_bab, times = tod_bab, period = 24,
                              min_rhythm_pval = RHYTHM_PVAL, verbose = FALSE)
bio_bab$ngenes <- NGENES

pA <- .run_sc_comparison(mat_bab, tod_bab, bio_bab, N_GRID_BAB, B_VAL,
                          sprintf("Baboon LUN (n=%d, r~1.72)", ncol(mat_bab)),
                          file.path(out_dir, "results", "panelA_baboon"))


# =======================================================================
# PANEL B: Mouse D1 striatum (n=45, passive, r~0.65)
# =======================================================================
cat("\n================================================================\n")
cat("PANEL B: Mouse D1 striatum (n=45, passive, r~0.65)\n")
cat("================================================================\n")

N_GRID_D1 <- if (SMOKE_TEST) c(40L, 80L, 120L) else c(40L, 60L, 80L, 120L, 160L, 200L)

pheno   <- read.csv("data/mouse_clinicalinfo_03082021_rmOutliers.csv", row.names = 1)
prep_d1 <- prepCircadianData("data/mouse_D1D2_logCPMfiltered_counts.csv",
                              times      = "time",
                              input_type = "counts",
                              pheno      = pheno,
                              sample_col = "sample")
d1_samp <- pheno$sample[pheno$cell == "D1"]
mat_d1  <- prep_d1$data[, colnames(prep_d1$data) %in% d1_samp, drop = FALSE]
tod_d1  <- pheno$time[match(colnames(mat_d1), pheno$sample)]
mat_d1  <- mat_d1[rowSums(mat_d1 > 1) >= floor(ncol(mat_d1) / 2), , drop = FALSE]
set.seed(GLOBAL_SEED)
g_d1    <- sample(nrow(mat_d1), min(NGENES, nrow(mat_d1)))
mat_d1  <- mat_d1[g_d1, , drop = FALSE]
rm(pheno, prep_d1)

bio_d1 <- estCircadianParam(mat_d1, times = tod_d1, period = 24,
                              min_rhythm_pval = RHYTHM_PVAL, verbose = FALSE)
bio_d1$ngenes <- NGENES

pB <- .run_sc_comparison(mat_d1, tod_d1, bio_d1, N_GRID_D1, B_VAL,
                          sprintf("Mouse D1 (n=%d, r~0.65)", ncol(mat_d1)),
                          file.path(out_dir, "results", "panelB_d1"))


# =======================================================================
# PANEL C: Seney CTL ACC (n=60, passive, r~0.79)
# =======================================================================
cat("\n================================================================\n")
cat("PANEL C: Seney CTL ACC (n=60, passive, r~0.79)\n")
cat("================================================================\n")

N_GRID_SEN <- if (SMOKE_TEST) c(40L, 80L, 120L) else c(40L, 80L, 120L, 160L, 200L)

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
mat_sen  <- prep_sen$data[, ctrl_idx, drop = FALSE]
tod_sen  <- prep_sen$times[ctrl_idx]
set.seed(GLOBAL_SEED)
g_sen    <- sample(nrow(mat_sen), min(NGENES, nrow(mat_sen)))
mat_sen  <- mat_sen[g_sen, , drop = FALSE]
rm(meta_s, tod_s, expr_raw, prep_sen)

bio_sen <- estCircadianParam(mat_sen, times = tod_sen, period = 24,
                              min_rhythm_pval = RHYTHM_PVAL, verbose = FALSE)
bio_sen$ngenes <- NGENES

pC <- .run_sc_comparison(mat_sen, tod_sen, bio_sen, N_GRID_SEN, B_VAL,
                          sprintf("Seney CTL (n=%d, r~0.79)", ncol(mat_sen)),
                          file.path(out_dir, "results", "panelC_seney"))


# =======================================================================
# Combined figure
# =======================================================================
cat("\n=== Generating combined figure ===\n")
panels <- Filter(Negate(is.null), list(pA, pB, pC))

if (length(panels) > 0) {
  pdf(file.path(out_dir, "figures", "fig_bootstrap_sc.pdf"),
      width = 14, height = 5)
  par(mfrow = c(1, length(panels)), mar = c(4, 4, 3, 1), las = 1)

  panel_letters <- LETTERS[seq_along(panels)]
  for (pi in seq_along(panels)) {
    p      <- panels[[pi]]
    N_grid <- p$N_grid
    ts_pwr <- p$ts_result$power_df$power
    bt_mn  <- p$bt_result$power_mean[, 1, 1]
    bt_lo  <- p$bt_result$power_ci_lo[, 1, 1]
    bt_hi  <- p$bt_result$power_ci_hi[, 1, 1]

    plot(N_grid, ts_pwr, type = "l", lwd = 2, col = "steelblue",
         ylim = c(0, 1), xlab = "N (total samples)", ylab = "Power",
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
  cat(sprintf("Saved: %s/figures/fig_bootstrap_sc.pdf\n", out_dir))
}

saveRDS(panels, file.path(out_dir, "results", "all_panels.rds"))
cat("\n=== Done ===\n")
