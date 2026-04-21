#' =======================================================================
#' 08c_bootstrap_seney.R — Two-Stage vs Bootstrap: Seney CTL ACC (n=60)
#' =======================================================================
#' Split from 08_two_stage_vs_bootstrap_realdata.R. Section 3 only.
#' See 08a_bootstrap_baboon.R header for parallel launch instructions.

POWERSIM_ROOT <- {
  env <- Sys.getenv("POWERSIM_ROOT", unset = "")
  if (nzchar(env)) env else
    "/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim"
}
setwd(POWERSIM_ROOT)

SMOKE       <- nzchar(Sys.getenv("POWERSIM_SMOKE"))
NGENES      <- if (SMOKE) 300L else 3000L
NBOOT       <- if (SMOKE) 5L  else 50L
NSIMS       <- if (SMOKE) 5L  else 30L
NSIMS_INNER <- if (SMOKE) 5L  else 20L
NCORES      <- if (SMOKE) 2L  else 20L
RHYTHM_PVAL <- 0.05

B_SENEY       <- 4L
N_GRID_SENEY  <- if (SMOKE) c(40L, 100L, 200L) else c(40L, 80L, 120L, 160L, 200L, 300L)
DATA_SENEY_META <- "data/MD5_MetaData_1-15-25.xlsx"
DATA_SENEY_TOD  <- "data/TOD.xlsx"
DATA_SENEY_EXPR <- "data/ACC_RNA_filtered_normalized.csv"

RUN_TAG <- Sys.getenv("RUN_TAG", format(Sys.Date(), "%Y%m%d"))
out_dir <- file.path("output", paste0("08_two_stage_vs_bootstrap_", RUN_TAG))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Mode: %s | NGENES=%d | NBOOT=%d | NSIMS=%d | NSIMS_INNER=%d\n",
            if (SMOKE) "SMOKE" else "PRODUCTION", NGENES, NBOOT, NSIMS, NSIMS_INNER))
cat(sprintf("Output -> %s/\n\n", out_dir))

suppressPackageStartupMessages(library(readxl))
source_dir <- file.path(getwd(), "code")
old_wd <- setwd(source_dir); source("setup.R"); setwd(old_wd)

opts_analysis <- CircadianAnalysisOptions(
  alpha = 0.05, p.adjust.method = "BH",
  fdr_thresholds = c(0.05), reference_n = 60L
)

.run_comparison <- function(pilot_data, pilot_times, bio_diff_opts,
                             N_grid, B_val, design_type, label, color,
                             out_prefix) {
  cat(sprintf("\n--- %s ---\n", label))
  cat(sprintf("    Pilot: n=%d | design=%s | B=%d\n", ncol(pilot_data), design_type, B_val))

  design_vec <- if (design_type == "active") {
    seq(0, 24, length.out = B_val + 1)[seq_len(B_val)]
  } else {
    pilot_times
  }

  cat("  Running two-stage...\n")
  design_opts_ts <- CircadianDesignOptions(
    sample_sizes = N_grid, nsims = NSIMS, design = design_type,
    cts = if (design_type == "active") design_vec else pilot_times,
    test_types = "DR"
  )
  ts_result <- tryCatch(
    runTwoStagePower(pilot_data = pilot_data, pilot_times = pilot_times,
                     design.opts = design_opts_ts, analysis.opts = opts_analysis,
                     bio_diff.opts = bio_diff_opts, min_rhythm_pval = RHYTHM_PVAL,
                     test_type = "DR", verbose = FALSE),
    error = function(e) { warning(sprintf("Two-stage failed: %s", e$message)); NULL }
  )

  cat(sprintf("  Running bootstrap (B=%d, %d draws)...\n", B_val, NBOOT))
  boot_opts <- CircadianBootstrapOptions(
    design_vector = design_vec, B_values = B_val, N_values = N_grid,
    nboot = NBOOT, nsims_inner = NSIMS_INNER, design = design_type, seed = 42L
  )
  boot_result <- tryCatch(
    runBootstrapDesignGrid(pilot_data = pilot_data, pilot_times = pilot_times,
                           boot.opts = boot_opts, analysis.opts = opts_analysis,
                           bio_diff.opts = bio_diff_opts, mode = "differential",
                           verbose = FALSE, mc.cores = NCORES),
    error = function(e) { warning(sprintf("Bootstrap failed: %s", e$message)); NULL }
  )

  if (is.null(ts_result) || is.null(boot_result)) return(NULL)

  comparison <- compareDesignApproaches(
    two_stage_result = ts_result, bootstrap_result = boot_result,
    test_type = "DR", target_power = 0.80
  )

  saveRDS(list(two_stage = ts_result, boot = boot_result, comparison = comparison),
          paste0(out_prefix, ".rds"))

  comp_df   <- comparison$comparison
  ci_widths <- comp_df$boot_ci_hi - comp_df$boot_ci_lo

  cat(sprintf("  Two-stage n80:        %s\n",
              ifelse(is.na(comparison$n80_two_stage), ">max(N)", comparison$n80_two_stage)))
  cat(sprintf("  Bootstrap n80 median: %s  [95%% CI: %s, %s]\n",
              ifelse(is.na(comparison$n80_boot_median), ">max(N)", round(comparison$n80_boot_median)),
              ifelse(is.na(comparison$n80_boot_lo), "NA", comparison$n80_boot_lo),
              ifelse(is.na(comparison$n80_boot_hi), "NA", comparison$n80_boot_hi)))
  cat(sprintf("  Mean CI width: %.0f pp\n", 100 * mean(ci_widths, na.rm = TRUE)))

  plotDesignComparison(comparison, target_power = 0.80, panels = "A",
                       output_file = paste0(out_prefix, ".pdf"))
  cat(sprintf("  Figure: %s.pdf\n", out_prefix))

  list(label = label, color = color, n_pilot = ncol(pilot_data),
       design_type = design_type, B = B_val, N_grid = N_grid,
       comparison = comparison, ci_widths = ci_widths,
       ts_result = ts_result, boot_result = boot_result)
}

# =======================================================================
# SECTION 3: SENEY CONTROL — CTL pilot (n=60, passive)
# =======================================================================
cat("====================================================================\n")
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
prep_seney   <- prepCircadianData(expr_seney[, ok_s], times = tod_hour_s[ok_s], input_type = "log2")
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

bio_seney <- estCircadianParamTwoGroup(
  data_1 = mat_ctrl, data_2 = mat_mdd,
  times_1 = tod_ctrl, times_2 = tod_mdd,
  period = 24, min_rhythm_pval = RHYTHM_PVAL, verbose = FALSE
)
cat(sprintf("Seney: prop_DR=%.1f%%  r_med(CTL)=%.2f\n",
            100 * bio_seney$prop_DR,
            bio_seney$diagnostics$r1_snr["median"]))

s3 <- .run_comparison(
  pilot_data = mat_ctrl, pilot_times = tod_ctrl, bio_diff_opts = bio_seney,
  N_grid = N_GRID_SENEY, B_val = B_SENEY, design_type = "passive",
  label = sprintf("Seney CTL (n=%d, passive)", ncol(mat_ctrl)), color = "firebrick",
  out_prefix = file.path(out_dir, "s3_seney_comparison")
)

cat(sprintf("\nDone. Output: %s/\n", out_dir))
