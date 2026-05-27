#' =======================================================================
#' 08_two_stage_vs_bootstrap_realdata.R
#'    Two-Stage vs Bootstrap Design on Real Pilot Data
#' =======================================================================
#'
#' PURPOSE
#'   Extends the synthetic calibration in 02_calibration.R to three real
#'   pilot datasets with very different pilot sizes:
#'
#'     Dataset         Pilot n   Design     Expected CI width
#'     Baboon LUN      n=12      Active     WIDE   (small pilot → high uncertainty)
#'     Mouse D1D2 D1   n=45      Active     Medium
#'     Seney CTL ACC   n=60      Passive    Narrower
#'
#'   For each dataset both methods are run at a fixed single B value:
#'     Two-stage:   estimate once from pilot → one power curve, no CI
#'     Bootstrap:   resample pilot genes NBOOT times → power ± 95% CI
#'
#'   The comparison directly answers:
#'     Q1: Do two-stage and bootstrap agree on the central n80?
#'         (they should, since bootstrap median ≈ two-stage point estimate)
#'     Q2: How much does CI width differ across pilot sizes?
#'         (Baboon CI >> D1D2 CI >> Seney CI, demonstrating that
#'          small pilots give false precision in the two-stage approach)
#'
#' WHY SINGLE B?
#'   compareDesignApproaches() requires a single fixed-B bootstrap result.
#'   CI bands then reflect parameter-estimation uncertainty only, not
#'   design variation. We use the full B of each pilot (best possible design).
#'
#' OUTPUTS
#'   output/08_two_stage_vs_bootstrap_realdata_<ts>/
#'     s1_baboon_comparison.pdf
#'     s2_d1d2_comparison.pdf
#'     s3_seney_comparison.pdf
#'     s4_ci_width_summary.pdf    — CI width vs pilot size across all 3
#'     comparison_summary.txt
#'
#' USAGE
#'   Rscript examples/publication/08_two_stage_vs_bootstrap_realdata.R
#'   POWERSIM_SMOKE=1 Rscript examples/publication/08_two_stage_vs_bootstrap_realdata.R
#'
#' @author Thien Quy Pham

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
NBOOT       <- if (SMOKE) 5L    else 50L
NSIMS       <- if (SMOKE) 5L    else 30L    # two-stage inner sims
NSIMS_INNER <- if (SMOKE) 5L    else 20L    # bootstrap inner sims per draw
RHYTHM_PVAL <- 0.05

# Fixed single-B per dataset (full pilot coverage → best point estimate)
# Note: N grids must be divisible by B
B_BABOON <- 12L   # all 12 ZT points; N divisible by 12
B_D1D2   <- 6L    # all 6 ZT points; N divisible by 6
B_SENEY  <- 4L    # passive placeholder; N divisible by 4

# N grids: smoke = 3 values; production = 6-7 values covering expected n80
N_GRID_BABOON <- if (SMOKE) c(12L, 24L, 36L) else
                   c(12L, 24L, 36L, 48L, 60L, 72L, 96L)
N_GRID_D1D2   <- if (SMOKE) c(24L, 48L, 72L) else
                   c(24L, 36L, 48L, 60L, 72L, 96L, 120L)
N_GRID_SENEY  <- if (SMOKE) c(40L, 100L, 200L) else
                   c(40L, 80L, 120L, 160L, 200L, 300L)   # all div by 4

DATA_BABOON     <- "data/CAMO_PRC_hmb.RData"
DATA_D1D2_PHENO <- "data/mouse_clinicalinfo_03082021_rmOutliers.csv"
DATA_D1D2_EXPR  <- "data/mouse_D1D2_logCPMfiltered_counts.csv"
DATA_SENEY_META <- "data/MD5_MetaData_1-15-25.xlsx"
DATA_SENEY_TOD  <- "data/TOD.xlsx"
DATA_SENEY_EXPR <- "data/ACC_RNA_filtered_normalized.csv"

cat(sprintf("Mode: %s | NGENES=%d | NBOOT=%d | NSIMS=%d | NSIMS_INNER=%d\n",
            if (SMOKE) "SMOKE" else "PRODUCTION",
            NGENES, NBOOT, NSIMS, NSIMS_INNER))

# -----------------------------------------------------------------------
# 2. Setup
# -----------------------------------------------------------------------
suppressPackageStartupMessages(library(readxl))

source_dir <- file.path(getwd(), "code")
old_wd <- setwd(source_dir)
source("setup.R")
setwd(old_wd)

ts      <- format(Sys.time(), "%Y%m%d_%H%M")
out_dir <- file.path("output", paste0("08_two_stage_vs_bootstrap_realdata_", ts))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("Output -> %s/\n\n", out_dir))

opts_analysis <- CircadianAnalysisOptions(
  alpha           = 0.05,
  p.adjust.method = "BH",
  fdr_thresholds  = c(0.05),
  reference_n     = 60L
)

# Storage for cross-dataset summary
results_all <- list()
summary_lines <- c(
  "Two-Stage vs Bootstrap Real-Data Comparison",
  sprintf("Mode: %s | NGENES=%d | NBOOT=%d | NSIMS=%d",
          if (SMOKE) "SMOKE" else "PRODUCTION", NGENES, NBOOT, NSIMS),
  ""
)

# -----------------------------------------------------------------------
# Helper: run both methods + compare for one dataset
# -----------------------------------------------------------------------
.run_comparison <- function(pilot_data, pilot_times, bio_diff_opts,
                             N_grid, B_val, design_type, label, color,
                             out_prefix) {
  cat(sprintf("\n--- %s ---\n", label))
  cat(sprintf("    Pilot: n=%d subjects | design=%s | B=%d\n",
              ncol(pilot_data), design_type, B_val))

  # Design vector for active designs (evenly-spaced B points)
  # For passive: use observed pilot TOD directly
  design_vec <- if (design_type == "active") {
    seq(0, 24, length.out = B_val + 1)[seq_len(B_val)]
  } else {
    pilot_times  # empirical TOD distribution
  }

  # --- Two-stage ---
  cat("  Running two-stage...\n")
  design_opts_ts <- CircadianDesignOptions(
    sample_sizes = N_grid,
    nsims        = NSIMS,
    design       = design_type,
    cts          = if (design_type == "active") design_vec else pilot_times,
    test_types   = "DR"
  )

  ts_result <- tryCatch(
    runTwoStagePower(
      pilot_data      = pilot_data,
      pilot_times     = pilot_times,
      design.opts     = design_opts_ts,
      analysis.opts   = opts_analysis,
      bio_diff.opts   = bio_diff_opts,
      min_rhythm_pval = RHYTHM_PVAL,
      test_type       = "DR",
      verbose         = FALSE
    ),
    error = function(e) {
      warning(sprintf("  Two-stage failed: %s", e$message)); NULL
    }
  )

  # --- Bootstrap (single B) ---
  cat(sprintf("  Running bootstrap (B=%d, %d draws)...\n", B_val, NBOOT))
  boot_opts <- CircadianBootstrapOptions(
    design_vector = design_vec,
    B_values      = B_val,
    N_values      = N_grid,
    nboot         = NBOOT,
    nsims_inner   = NSIMS_INNER,
    design        = design_type,
    seed          = 42L
  )

  boot_result <- tryCatch(
    runBootstrapDesignGrid(
      pilot_data    = pilot_data,
      pilot_times   = pilot_times,
      boot.opts     = boot_opts,
      analysis.opts = opts_analysis,
      bio_diff.opts = bio_diff_opts,
      verbose       = FALSE
    ),
    error = function(e) {
      warning(sprintf("  Bootstrap failed: %s", e$message)); NULL
    }
  )

  if (is.null(ts_result) || is.null(boot_result)) return(NULL)

  # --- Compare ---
  comparison <- compareDesignApproaches(
    two_stage_result = ts_result,
    bootstrap_result = boot_result,
    test_type        = "DR",
    target_power     = 0.80
  )

  # Save RDS
  saveRDS(list(two_stage = ts_result, boot = boot_result, comparison = comparison),
          paste0(out_prefix, ".rds"))

  # Per-N CI widths
  comp_df   <- comparison$comparison
  ci_widths <- comp_df$boot_ci_hi - comp_df$boot_ci_lo

  cat(sprintf("  Two-stage n80:        %s\n",
              ifelse(is.na(comparison$n80_two_stage), ">max(N)",
                     comparison$n80_two_stage)))
  cat(sprintf("  Bootstrap n80 median: %s  [95%% CI: %s, %s]\n",
              ifelse(is.na(comparison$n80_boot_median), ">max(N)",
                     round(comparison$n80_boot_median)),
              ifelse(is.na(comparison$n80_boot_lo), "NA",
                     comparison$n80_boot_lo),
              ifelse(is.na(comparison$n80_boot_hi), "NA",
                     comparison$n80_boot_hi)))
  cat(sprintf("  Mean CI width: %.3f (%.0f pp)\n",
              mean(ci_widths, na.rm = TRUE),
              100 * mean(ci_widths, na.rm = TRUE)))

  # Figure
  plotDesignComparison(
    comparison,
    target_power = 0.80,
    panels       = "A",
    output_file  = paste0(out_prefix, ".pdf")
  )
  cat(sprintf("  Figure: %s.pdf\n", out_prefix))

  list(
    label         = label,
    color         = color,
    n_pilot       = ncol(pilot_data),
    design_type   = design_type,
    B             = B_val,
    N_grid        = N_grid,
    comparison    = comparison,
    ci_widths     = ci_widths,
    ts_result     = ts_result,
    boot_result   = boot_result
  )
}


# =======================================================================
# SECTION 1: BABOON — LUN pilot (n=12, smallest → widest CI)
# =======================================================================
cat("\n====================================================================\n")
cat("SECTION 1: Baboon LUN vs CER (n=12 pilot)\n")
cat("====================================================================\n\n")

load(DATA_BABOON)
bab_expr <- baboon_withTOD$baboon
bab_tod  <- baboon_withTOD$tod

prep_lun  <- prepCircadianData(bab_expr[["LUN"]], times = bab_tod[["LUN"]], input_type = "cpm")
prep_cerb <- prepCircadianData(bab_expr[["CER"]], times = bab_tod[["CER"]], input_type = "cpm")
mat_lun   <- prep_lun$data[rowSums(prep_lun$data > 0) >= 6, , drop = FALSE]
mat_cerb  <- prep_cerb$data[rowSums(prep_cerb$data > 0) >= 6, , drop = FALSE]
tod_lun   <- prep_lun$times
common_b  <- intersect(rownames(mat_lun), rownames(mat_cerb))
set.seed(2)
g_idx_b   <- sample(common_b, min(NGENES, length(common_b)))
mat_lun_s <- mat_lun[g_idx_b, , drop = FALSE]
mat_cerb_s <- mat_cerb[g_idx_b, , drop = FALSE]

fit_lun   <- fitCosinorAll(mat_lun_s, times = tod_lun, period = 24)
fit_cerb  <- fitCosinorAll(mat_cerb_s, times = prep_cerb$times, period = 24)
rhy_lun   <- !is.na(fit_lun$pvalue) & fit_lun$pvalue < RHYTHM_PVAL
rhy_cerb  <- !is.na(fit_cerb$pvalue) & fit_cerb$pvalue < RHYTHM_PVAL
prop_DR_b <- mean(xor(rhy_lun, rhy_cerb))

bio_bab <- estCircadianParam(
  data = mat_lun_s, times = tod_lun, period = 24,
  prop_DR = prop_DR_b, prop_DP = 0,
  min_rhythm_pval = RHYTHM_PVAL, verbose = FALSE
)
cat(sprintf("Baboon: prop_DR=%.1f%%  r_med=%.2f\n",
            100 * prop_DR_b,
            median(as.numeric(fit_lun$A[rhy_lun]) /
                   as.numeric(fit_lun$sigma[rhy_lun]), na.rm = TRUE)))

s1 <- .run_comparison(
  pilot_data    = mat_lun_s,
  pilot_times   = tod_lun,
  bio_diff_opts = bio_bab,
  N_grid        = N_GRID_BABOON,
  B_val         = B_BABOON,
  design_type   = "active",
  label         = sprintf("Baboon LUN (n=%d)", ncol(mat_lun_s)),
  color         = "darkorange",
  out_prefix    = file.path(out_dir, "s1_baboon_comparison")
)
if (!is.null(s1)) {
  results_all[["baboon"]] <- s1
  summary_lines <- c(summary_lines,
    sprintf("BABOON  (n_pilot=%d, B=%d): two-stage n80=%s  boot_median n80=%s  mean_CI_width=%.0f pp",
            ncol(mat_lun_s), B_BABOON,
            ifelse(is.na(s1$comparison$n80_two_stage), ">max", s1$comparison$n80_two_stage),
            ifelse(is.na(s1$comparison$n80_boot_median), ">max", round(s1$comparison$n80_boot_median)),
            100 * mean(s1$ci_widths, na.rm = TRUE)))
}


# =======================================================================
# SECTION 2: MOUSE D1D2 — D1 pilot (n=45, medium)
# =======================================================================
cat("\n====================================================================\n")
cat("SECTION 2: Mouse D1 vs D2 (n=45 pilot)\n")
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
cat(sprintf("D1D2: prop_DR=%.1f%%  r_med=%.2f\n",
            100 * prop_DR_d,
            median(as.numeric(fit_d1$A[rhy_d1]) /
                   as.numeric(fit_d1$sigma[rhy_d1]), na.rm = TRUE)))

s2 <- .run_comparison(
  pilot_data    = mat_d1_s,
  pilot_times   = tod_d1,
  bio_diff_opts = bio_d1d2,
  N_grid        = N_GRID_D1D2,
  B_val         = B_D1D2,
  design_type   = "active",
  label         = sprintf("D1D2 D1 (n=%d)", ncol(mat_d1_s)),
  color         = "forestgreen",
  out_prefix    = file.path(out_dir, "s2_d1d2_comparison")
)
if (!is.null(s2)) {
  results_all[["d1d2"]] <- s2
  summary_lines <- c(summary_lines,
    sprintf("D1D2    (n_pilot=%d, B=%d): two-stage n80=%s  boot_median n80=%s  mean_CI_width=%.0f pp",
            ncol(mat_d1_s), B_D1D2,
            ifelse(is.na(s2$comparison$n80_two_stage), ">max", s2$comparison$n80_two_stage),
            ifelse(is.na(s2$comparison$n80_boot_median), ">max", round(s2$comparison$n80_boot_median)),
            100 * mean(s2$ci_widths, na.rm = TRUE)))
}


# =======================================================================
# SECTION 3: SENEY CONTROL — CTL pilot (n=60, passive, largest)
# =======================================================================
cat("\n====================================================================\n")
cat("SECTION 3: Seney MDD vs Control ACC (n=60 CTL pilot, passive)\n")
cat("====================================================================\n\n")

meta_seney <- read_excel(DATA_SENEY_META)
tod_seney  <- read_excel(DATA_SENEY_TOD)
expr_seney <- as.matrix(read.csv(DATA_SENEY_EXPR, row.names = 1, check.names = FALSE))

col_ids_s  <- gsub("[A-Za-z]+$", "", colnames(expr_seney))
meta_idx_s <- match(col_ids_s, as.character(meta_seney$HU_NUM))
tod_idx_s  <- match(col_ids_s, as.character(tod_seney$HU_NUM))
tod_hour_s <- as.numeric(format(tod_seney$OFFC_TIME[tod_idx_s], "%H")) +
              as.numeric(format(tod_seney$OFFC_TIME[tod_idx_s], "%M")) / 60
disease_s  <- meta_seney$Disease[meta_idx_s]

ok_s <- !is.na(disease_s) & !is.na(tod_hour_s)
prep_seney <- prepCircadianData(expr_seney[, ok_s], times = tod_hour_s[ok_s],
                                input_type = "log2")
expr_seney_f <- prep_seney$data
tod_hour_f   <- prep_seney$times
disease_f    <- disease_s[ok_s]

ctrl_idx <- disease_f == 1
mdd_idx  <- disease_f == 2
tod_ctrl <- tod_hour_f[ctrl_idx]
tod_mdd  <- tod_hour_f[mdd_idx]

set.seed(3)
g_idx_s  <- sample(nrow(expr_seney_f), min(NGENES, nrow(expr_seney_f)))
mat_ctrl <- expr_seney_f[g_idx_s, ctrl_idx, drop = FALSE]
mat_mdd  <- expr_seney_f[g_idx_s, mdd_idx,  drop = FALSE]

# Use two-group estimation to get prop_DR from both groups
bio_seney <- estCircadianParamTwoGroup(
  data_1          = mat_ctrl,
  data_2          = mat_mdd,
  times_1         = tod_ctrl,
  times_2         = tod_mdd,
  period          = 24,
  min_rhythm_pval = RHYTHM_PVAL,
  verbose         = FALSE
)
cat(sprintf("Seney: prop_DR=%.1f%%  r_med(CTL)=%.2f\n",
            100 * bio_seney$prop_DR,
            bio_seney$diagnostics$r1_snr["median"]))

# For the comparison, pilot = CTL group; prop_DR already in bio_seney
s3 <- .run_comparison(
  pilot_data    = mat_ctrl,
  pilot_times   = tod_ctrl,
  bio_diff_opts = bio_seney,
  N_grid        = N_GRID_SENEY,
  B_val         = B_SENEY,
  design_type   = "passive",
  label         = sprintf("Seney CTL (n=%d, passive)", ncol(mat_ctrl)),
  color         = "firebrick",
  out_prefix    = file.path(out_dir, "s3_seney_comparison")
)
if (!is.null(s3)) {
  results_all[["seney"]] <- s3
  summary_lines <- c(summary_lines,
    sprintf("SENEY   (n_pilot=%d, B=%d, passive): two-stage n80=%s  boot_median n80=%s  mean_CI_width=%.0f pp",
            ncol(mat_ctrl), B_SENEY,
            ifelse(is.na(s3$comparison$n80_two_stage), ">max", s3$comparison$n80_two_stage),
            ifelse(is.na(s3$comparison$n80_boot_median), ">max", round(s3$comparison$n80_boot_median)),
            100 * mean(s3$ci_widths, na.rm = TRUE)))
}


# =======================================================================
# SECTION 4: SUMMARY — CI width vs pilot size
# =======================================================================
# Key message: bootstrap CI width ∝ 1/√n_pilot.
# Two-stage gives a single number regardless of pilot size → false precision
# when pilot is small (Baboon n=12).

cat("\n====================================================================\n")
cat("SECTION 4: CI width summary across pilot sizes\n")
cat("====================================================================\n\n")

if (length(results_all) >= 2) {

  fig_s4 <- file.path(out_dir, "s4_ci_width_summary.pdf")
  pdf(fig_s4, width = 12, height = 5)
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.5, 1.5))

  # Panel A: power curves for all 3 datasets on same axes
  # Two-stage = solid line, bootstrap = dashed ± shaded CI
  all_N <- sort(unique(unlist(lapply(results_all, `[[`, "N_grid"))))
  plot(NA, xlim = range(all_N), ylim = c(0, 1),
       xlab = "N per group", ylab = "Power",
       main = "Two-stage (solid) vs Bootstrap (dashed ± CI)\nby dataset",
       las = 1)
  abline(h = 0.80, lty = 3, col = "gray50")

  for (nm in names(results_all)) {
    sc  <- results_all[[nm]]
    cmp <- sc$comparison$comparison
    col <- sc$color

    # Bootstrap CI ribbon first (so two-stage line is drawn on top)
    polygon(c(cmp$n, rev(cmp$n)),
            c(cmp$boot_ci_lo, rev(cmp$boot_ci_hi)),
            col = adjustcolor(col, 0.12), border = NA)
    lines(cmp$n, cmp$boot_power_mean, col = col, lwd = 2, lty = 2)

    # Two-stage: solid line drawn last (on top) + SE bars
    lines(cmp$n, cmp$two_stage_power, col = col, lwd = 2.5, lty = 1)
    points(cmp$n, cmp$two_stage_power, col = col, pch = 16, cex = 0.8)
    arrows(cmp$n, cmp$two_stage_power - cmp$two_stage_se,
           cmp$n, cmp$two_stage_power + cmp$two_stage_se,
           length = 0.05, angle = 90, code = 3, col = col, lwd = 1.2)
  }

  legend("bottomright",
         legend = sapply(results_all, `[[`, "label"),
         col    = sapply(results_all, `[[`, "color"),
         lwd    = 2, bty = "n", cex = 0.78)

  # Panel B: mean CI width vs pilot size — the money plot
  pilot_ns   <- sapply(results_all, `[[`, "n_pilot")
  mean_widths <- sapply(results_all, function(sc)
                          100 * mean(sc$ci_widths, na.rm = TRUE))
  cols_bar   <- sapply(results_all, `[[`, "color")
  labels_bar <- sapply(results_all, `[[`, "label")

  ord <- order(pilot_ns)
  bp  <- barplot(mean_widths[ord],
                 names.arg = labels_bar[ord],
                 col       = cols_bar[ord],
                 ylim      = c(0, max(mean_widths) * 1.25),
                 ylab      = "Mean bootstrap CI width (pp)",
                 main      = "CI width vs pilot size\n(wider = more uncertainty)",
                 las       = 2, cex.names = 0.75)

  # Annotate with pilot n
  text(bp, mean_widths[ord] + max(mean_widths) * 0.04,
       sprintf("n=%d", pilot_ns[ord]),
       cex = 0.85, col = "gray20")

  dev.off()
  cat(sprintf("Figure: %s\n", fig_s4))

} else {
  cat("Fewer than 2 datasets succeeded; skipping summary figure.\n")
}


# =======================================================================
# WRAP-UP
# =======================================================================
summary_lines <- c(summary_lines, "",
  "Key interpretation:",
  "  - Two-stage n80 and bootstrap median n80 should be close (same point estimate)",
  "  - Bootstrap CI width reveals estimation uncertainty from finite pilot",
  "  - CI width should decrease as pilot size increases (∝ 1/sqrt(n_pilot))",
  "  - When pilot is small (Baboon n=12), two-stage gives false precision;",
  "    bootstrap CI is the only honest guide to n80 uncertainty."
)

writeLines(summary_lines, file.path(out_dir, "comparison_summary.txt"))

cat("\n====================================================================\n")
cat("08_two_stage_vs_bootstrap_realdata COMPLETE\n")
cat("====================================================================\n\n")
cat(paste(summary_lines, collapse = "\n"), "\n")
cat(sprintf("\nOutput: %s/\n", out_dir))
cat("Done.\n")
