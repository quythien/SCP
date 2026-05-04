#' =======================================================================
#' fig4b_fmm_differential.R — FMM Waveform Robustness: Differential Power
#' =======================================================================
#'
#' Companion to fig4_fmm_violation.R (single-cohort). Shows how DCP
#' differential power (DR, DP, DM) degrades as the FMM waveform departs
#' from a pure sinusoid.
#'
#' Dataset: GTEx Adrenal Gland vs Liver (PASSIVE design — post-mortem GTEx).
#' Same dataset as Figure 2.
#'
#' Design: passive — collection times KDE-sampled from pilot TOD distributions
#' (bio$cts for ADR group 1, bio$cts2 for LIV group 2).
#'
#' USAGE:
#'   Rscript examples/publication/fig4b_fmm_differential.R
#'   SMOKE_TEST=true Rscript examples/publication/fig4b_fmm_differential.R

SMOKE_TEST  <- identical(Sys.getenv("SMOKE_TEST"), "true")
GLOBAL_SEED <- 2025L
set.seed(GLOBAL_SEED)

old_wd <- setwd("code"); source("setup.R"); setwd(old_wd)
library(FMM)
library(ggplot2)

OMEGA_VALS  <- if (SMOKE_TEST) c(0.0, 0.5, 1.0) else c(0.0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0)
BETA_VALS   <- if (SMOKE_TEST) c(0, pi/2, pi)  else c(0, pi/4, pi/2, 3*pi/4, pi, 3*pi/2)
OMEGA_FIXED <- 0.5   # fixed omega for beta sweep
N_GRID      <- if (SMOKE_TEST) c(30L, 60L, 100L) else
               c(20L, 30L, 40L, 50L, 60L, 80L, 100L, 120L, 150L, 200L, 250L, 300L)
NSIMS       <- if (SMOKE_TEST) 5L   else 20L
NGENES     <- if (SMOKE_TEST) 200L else 5000L
FDR_THRESH <- 0.05
N_CORES    <- as.integer(Sys.getenv("MC_CORES", unset = "40"))

cat(sprintf("Mode     : %s\n", if (SMOKE_TEST) "SMOKE" else "PRODUCTION"))
cat(sprintf("omega    : %s\n", paste(OMEGA_VALS, collapse = ", ")))
cat(sprintf("N grid   : %s\n", paste(N_GRID, collapse = ", ")))
cat(sprintf("nsims    : %d\n", NSIMS))
cat(sprintf("ngenes   : %d\n", NGENES))
cat(sprintf("mc.cores : %d\n\n", N_CORES))

out_dir <- "output/fmm_differential"
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "results"),  recursive = TRUE, showWarnings = FALSE)

# =====================================================================
# 1. Load GTEx ADR vs LIV pilot (passive design)
# =====================================================================
cat("--- Loading GTEx ADR vs LIV pilot ---\n")
pilot_rds <- "data/gtex_adr_vs_liv_pilot.rds"
bio <- readRDS(pilot_rds)
bio$ngenes <- NGENES

cat(sprintf("ADR n_pilot=%d  LIV n_pilot=%d\n",
            length(bio$cts), length(bio$cts2)))
cat(sprintf("prop_DR=%.3f  prop_DP=%.3f  prop_DM=%.3f\n",
            bio$prop_DR, bio$prop_DP, bio$prop_DM %||% 0))

analysis <- CircadianAnalysisOptions(
  alpha           = 0.05,
  p.adjust.method = "BH",
  fdr_thresholds  = FDR_THRESH
)

# =====================================================================
# 2. Power simulation helper
# =====================================================================
get_power <- function(fdr_arr, diff_type_vec, type_codes, nsims, fdr_thresh) {
  mean(sapply(seq_len(nsims), function(i) {
    fdr <- fdr_arr[, 1, i]
    gt  <- diff_type_vec[[i]]
    idx <- if (length(type_codes) == 2) gt %in% type_codes else gt == type_codes
    if (!any(idx, na.rm = TRUE)) return(NA_real_)
    sum(fdr[idx] <= fdr_thresh, na.rm = TRUE) / sum(idx)
  }), na.rm = TRUE)
}

# =====================================================================
# 3. Sweep omega x N
# =====================================================================
all_rows <- list()

for (omega_val in OMEGA_VALS) {
  cat(sprintf("\nomega = %.1f\n", omega_val))
  omega_tag <- sprintf("%.1f", omega_val)

  for (N_val in N_GRID) {
    cat(sprintf("  N=%d ... ", N_val))

    # Passive design: sample collection times from pilot TOD KDE per sim
    set.seed(GLOBAL_SEED + round(omega_val * 100) * 1000L + N_val)

    sim_results <- parallel::mclapply(seq_len(NSIMS), function(s) {
      # KDE-sample passive times for this simulation replicate
      cts1 <- sampleTimesFromDist(N_val, bio$cts)
      cts2 <- sampleTimesFromDist(N_val, bio$cts2)

      sim_args <- list(
        ngenes        = NGENES,
        n1            = N_val, n2 = N_val,
        lBaselineExpr = bio$lBaselineExpr,
        lOD           = bio$lOD,
        lOD2          = bio$lOD2,
        amplitude     = bio$amplitude,
        amplitude2    = bio$amplitude2,
        sigma_rhythmic = bio$sigma_rhythmic,
        lBaselineExpr2 = bio$lBaselineExpr2,
        prop_rhythmic = bio$prop_rhythmic,
        prop_DR       = bio$prop_DR,
        prop_DP       = bio$prop_DP,
        prop_DM       = bio$prop_DM %||% 0,
        phase_diff    = bio$phase_diff,
        dp_shift_mode = bio$dp_shift_mode %||% "uniform",
        mesor_diff    = bio$mesor_diff %||% c(0.5, 2.0),
        design        = "passive",
        cts           = cts1,
        cts2          = cts2,
        sim.seed      = GLOBAL_SEED + s * 1000L,
        omega         = omega_val
      )

      sim <- do.call(simCircadianDiffFMM, sim_args)

      # Apply DCP differential detection — same pipeline as runSimsDiff / Fig 2
      gt         <- sim$ground_truth$diff_type
      gene_names <- paste0("Gene", seq_len(NGENES))
      ngenes     <- NGENES
      pval_DR    <- rep(1, ngenes)
      pval_DP    <- rep(1, ngenes)
      pval_DM    <- rep(1, ngenes)

      tryCatch({
        x1 <- format_for_DCP(sim$expr1, sim$times1, gene_names)
        x2 <- format_for_DCP(sim$expr2, sim$times2, gene_names)

        rhy <- DCP_Rhythmicity(x1 = x1, x2 = x2, method = "Sidak_FS",
                               period = 24, amp.cutoff = 0,
                               alpha = FDR_THRESH, CI = FALSE,
                               p.adjust.method = "BH",
                               parallel.ncores = 1L)

        # DR test
        n_dr <- sum(rhy$rhythm.joint$TOJR != "arrhy")
        if (n_dr > 0) {
          dr <- DCP_DiffR2(rhy, method = "LR", TOJR = NULL,
                           alpha = FDR_THRESH, p.adjust.method = "BH",
                           parallel.ncores = 1L)
          idx <- match(gene_names, dr$gname)
          pval_DR[!is.na(idx)] <- dr$p.R2[idx[!is.na(idx)]]
        }

        # DP + DM test
        n_both <- sum(rhy$rhythm.joint$TOJR == "both")
        if (n_both > 0) {
          dp <- DCP_DiffPar(rhy, Par = "A&phase&M", TOJR = NULL,
                            alpha = FDR_THRESH, p.adjust.method = "BH",
                            parallel.ncores = 1L)
          idx <- match(gene_names, dp$gname)
          if ("p.delta.peak" %in% colnames(dp))
            pval_DP[!is.na(idx)] <- dp$p.delta.peak[idx[!is.na(idx)]]
          if ("p.delta.M" %in% colnames(dp))
            pval_DM[!is.na(idx)] <- dp$p.delta.M[idx[!is.na(idx)]]
        }
      }, error = function(e) NULL)

      # BH correction across all genes (same as runSimsDiff)
      fdr_DR <- p.adjust(pval_DR, "BH")
      fdr_DP <- p.adjust(pval_DP, "BH")
      fdr_DM <- p.adjust(pval_DM, "BH")

      calc_pwr <- function(fdr, type_codes) {
        idx <- if (length(type_codes) == 2) gt %in% type_codes else gt == type_codes
        if (!any(idx, na.rm = TRUE)) return(NA_real_)
        sum(fdr[idx] <= FDR_THRESH, na.rm = TRUE) / sum(idx)
      }

      list(
        DR = calc_pwr(fdr_DR, c(2L, 3L)),
        DP = calc_pwr(fdr_DP, 4L),
        DM = calc_pwr(fdr_DM, 5L)
      )
    }, mc.cores = N_CORES)

    # Aggregate
    valid  <- Filter(Negate(is.null), sim_results)
    pwr_DR <- mean(sapply(valid, function(x) x$DR), na.rm = TRUE)
    pwr_DP <- mean(sapply(valid, function(x) x$DP), na.rm = TRUE)
    pwr_DM <- mean(sapply(valid, function(x) x$DM), na.rm = TRUE)
    se_DR  <- sd(sapply(valid, function(x) x$DR), na.rm = TRUE) / sqrt(length(valid))
    se_DP  <- sd(sapply(valid, function(x) x$DP), na.rm = TRUE) / sqrt(length(valid))
    se_DM  <- sd(sapply(valid, function(x) x$DM), na.rm = TRUE) / sqrt(length(valid))

    cat(sprintf("DR=%.2f DP=%.2f DM=%.2f\n", pwr_DR, pwr_DP, pwr_DM))

    all_rows <- c(all_rows, list(data.frame(
      omega = omega_val, N = N_val,
      DR = pwr_DR, DP = pwr_DP, DM = pwr_DM,
      SE_DR = se_DR, SE_DP = se_DP, SE_DM = se_DM,
      stringsAsFactors = FALSE
    )))
  }
  saveRDS(do.call(rbind, all_rows),
          file.path(out_dir, "results",
                    sprintf("results_FMM_diff_omega_%s.rds", omega_tag)))
}

# =====================================================================
# 4. Beta sweep: vary beta at fixed omega (same pipeline as omega loop)
# =====================================================================
cat("\n=== BETA SWEEP (omega fixed at", OMEGA_FIXED, ") ===\n")
beta_rows <- list()

for (beta_val in BETA_VALS) {
  beta_tag <- sprintf("%.4f", round(beta_val, 4))
  cat(sprintf("\nbeta = %.4f\n", beta_val))

  for (N_val in N_GRID) {
    cat(sprintf("  N=%d ... ", N_val))
    set.seed(GLOBAL_SEED + round(beta_val * 100) * 777L + N_val)

    sim_results <- parallel::mclapply(seq_len(NSIMS), function(s) {
      cts1 <- sampleTimesFromDist(N_val, bio$cts)
      cts2 <- sampleTimesFromDist(N_val, bio$cts2)
      sim_args <- list(
        ngenes = NGENES, n1 = N_val, n2 = N_val,
        lBaselineExpr = bio$lBaselineExpr, lOD = bio$lOD, lOD2 = bio$lOD2,
        amplitude = bio$amplitude, amplitude2 = bio$amplitude2,
        sigma_rhythmic = bio$sigma_rhythmic, lBaselineExpr2 = bio$lBaselineExpr2,
        prop_rhythmic = bio$prop_rhythmic, prop_DR = bio$prop_DR,
        prop_DP = bio$prop_DP, prop_DM = bio$prop_DM %||% 0,
        phase_diff = bio$phase_diff, dp_shift_mode = bio$dp_shift_mode %||% "uniform",
        mesor_diff = bio$mesor_diff %||% c(0.5, 2.0),
        design = "passive", cts = cts1, cts2 = cts2,
        sim.seed = GLOBAL_SEED + s * 1000L,
        omega = OMEGA_FIXED, beta = beta_val
      )
      sim <- do.call(simCircadianDiffFMM, sim_args)
      gt  <- sim$ground_truth$diff_type
      gene_names <- paste0("Gene", seq_len(NGENES))
      pval_DR <- rep(1, NGENES); pval_DP <- rep(1, NGENES); pval_DM <- rep(1, NGENES)
      tryCatch({
        x1 <- format_for_DCP(sim$expr1, sim$times1, gene_names)
        x2 <- format_for_DCP(sim$expr2, sim$times2, gene_names)
        rhy <- DCP_Rhythmicity(x1=x1, x2=x2, method="Sidak_FS", period=24,
                               amp.cutoff=0, alpha=FDR_THRESH, CI=FALSE,
                               p.adjust.method="BH", parallel.ncores=1L)
        n_dr <- sum(rhy$rhythm.joint$TOJR != "arrhy")
        if (n_dr > 0) {
          dr <- DCP_DiffR2(rhy, method="LR", TOJR=NULL, alpha=FDR_THRESH,
                           p.adjust.method="BH", parallel.ncores=1L)
          idx <- match(gene_names, dr$gname)
          pval_DR[!is.na(idx)] <- dr$p.R2[idx[!is.na(idx)]]
        }
        n_both <- sum(rhy$rhythm.joint$TOJR == "both")
        if (n_both > 0) {
          dp <- DCP_DiffPar(rhy, Par="A&phase&M", TOJR=NULL, alpha=FDR_THRESH,
                            p.adjust.method="BH", parallel.ncores=1L)
          idx <- match(gene_names, dp$gname)
          if ("p.delta.peak" %in% colnames(dp))
            pval_DP[!is.na(idx)] <- dp$p.delta.peak[idx[!is.na(idx)]]
          if ("p.delta.M" %in% colnames(dp))
            pval_DM[!is.na(idx)] <- dp$p.delta.M[idx[!is.na(idx)]]
        }
      }, error = function(e) NULL)
      fdr_DR <- p.adjust(pval_DR,"BH"); fdr_DP <- p.adjust(pval_DP,"BH"); fdr_DM <- p.adjust(pval_DM,"BH")
      calc_pwr <- function(fdr, tc) {
        idx <- if (length(tc)==2) gt %in% tc else gt==tc
        if (!any(idx,na.rm=TRUE)) return(NA_real_)
        sum(fdr[idx]<=FDR_THRESH,na.rm=TRUE)/sum(idx)
      }
      list(DR=calc_pwr(fdr_DR,c(2L,3L)), DP=calc_pwr(fdr_DP,4L), DM=calc_pwr(fdr_DM,5L))
    }, mc.cores = N_CORES)

    valid <- Filter(Negate(is.null), sim_results)
    pwr_DR <- mean(sapply(valid,function(x)x$DR),na.rm=TRUE)
    pwr_DP <- mean(sapply(valid,function(x)x$DP),na.rm=TRUE)
    pwr_DM <- mean(sapply(valid,function(x)x$DM),na.rm=TRUE)
    se_DR  <- sd(sapply(valid,function(x)x$DR),na.rm=TRUE)/sqrt(length(valid))
    se_DP  <- sd(sapply(valid,function(x)x$DP),na.rm=TRUE)/sqrt(length(valid))
    se_DM  <- sd(sapply(valid,function(x)x$DM),na.rm=TRUE)/sqrt(length(valid))
    cat(sprintf("DR=%.2f DP=%.2f DM=%.2f\n", pwr_DR, pwr_DP, pwr_DM))
    beta_rows <- c(beta_rows, list(data.frame(
      beta=beta_val, N=N_val, DR=pwr_DR, DP=pwr_DP, DM=pwr_DM,
      SE_DR=se_DR, SE_DP=se_DP, SE_DM=se_DM, stringsAsFactors=FALSE
    )))
  }
  saveRDS(do.call(rbind, beta_rows),
          file.path(out_dir, "results", sprintf("results_FMM_diff_beta_%s.rds", beta_tag)))
}

omega_df <- do.call(rbind, all_rows)
beta_df  <- do.call(rbind, beta_rows)
saveRDS(omega_df, file.path(out_dir, "results", "results_FMM_diff_all.rds"))
saveRDS(beta_df,  file.path(out_dir, "results", "results_FMM_diff_beta_all.rds"))
cat("\nSaved all results.\n")
full_df <- omega_df  # for plotFMMDifferential (omega sweep)

# =====================================================================
# 5. Plot using API function
# =====================================================================
fig4b_path <- file.path(out_dir, "figures", "fig4b_fmm_differential.pdf")
plotFMMDifferential(full_df, nsims = NSIMS,
                    dataset_label = "GTEx Adrenal Gland vs Liver (passive, r~1.0 vs r~0.4)",
                    output_file = fig4b_path, width = 14, height = 5)
cat(sprintf("Saved: %s\n", fig4b_path))
cat("\n=== Done ===\n")
