#' =======================================================================
#' run_pipeline2.R — Real-Pilot Bootstrap & Sensitivity Extensions
#' =======================================================================
#'
#' PURPOSE
#'   Real-pilot application of the new bootstrap and sensitivity modules.
#'   Extends run_pipeline.R; does NOT duplicate it.
#'
#'   Clean division:
#'     code/example.R      = methodological validation (synthetic data,
#'                           known-truth calibration, justification)
#'     examples/run_pipeline.R   = DR/DP/DA power from real pilot
#'     examples/run_pipeline2.R  = bootstrap design grid + uncertainty
#'                                 + optional sensitivity on real pilot
#'
#' WHAT THIS FILE ANSWERS
#'   "Given our real pilot data, what is the design recommendation
#'    and how uncertain is that recommendation?"
#'
#'   Section 3: Bootstrap design grid
#'     "For fixed total post-mortem subjects, how should I split
#'      across TOD bins (B) and subjects per bin (m)?"
#'     → passive design, real TOD distribution
#'     → primary deliverable
#'
#'   Section 4 [OPTIONAL]: Two-stage vs bootstrap on real pilot
#'     "How much is the n80 recommendation inflated or deflated
#'      by using a point estimate vs propagating pilot uncertainty?"
#'     → uncertainty illustration only; not validation
#'     → validation of the method lives in code/example.R Scenario 7-8
#'
#'   Section 5 [OPTIONAL]: Fourier misspecification sensitivity
#'     "If true waveforms have harmonics the cosinor misses,
#'      does that change the design recommendation?"
#'     → sensitivity analysis on the real-pilot conclusion
#'     → method justification lives in code/example.R Scenario 6
#'
#' WHAT NOT TO CONFLATE
#'   - known-truth validation  ≠  real-data application
#'   - bootstrap uncertainty   ≠  model misspecification
#'   - active B vs m tradeoff  ≠  passive planning under empirical TOD
#'
#' USAGE
#'   # Standalone (re-estimates parameters from raw data):
#'   Rscript examples/run_pipeline2.R
#'
#'   # Reuse already-estimated opts_bio from run_pipeline.R:
#'   Set PIPELINE1_RDS below to the path of a saved run_pipeline opts_bio.
#'
#' RUN SETTINGS (edit here)
#'   RUN_OPTIONAL_S4   <- TRUE    # two-stage vs bootstrap comparison
#'   RUN_OPTIONAL_S5   <- TRUE    # Fourier misspecification sensitivity
#'   S3_NBOOT          <- 50L     # bootstrap draws (≥50 for production)
#'   S3_NSIMS_INNER    <- 10L     # sims per bootstrap draw
#'
#' @author Thien Quy Pham


# =====================================================================
# RUN SETTINGS — edit these before launching
# =====================================================================

# Smoke test: POWERSIM_SMOKE=1 Rscript examples/publication/04_power_design.R
SMOKE <- nzchar(Sys.getenv("POWERSIM_SMOKE"))

RUN_OPTIONAL_S4 <- TRUE
RUN_OPTIONAL_S5 <- TRUE

S3_NBOOT       <- if (SMOKE) 3L   else 50L
S3_NSIMS_INNER <- if (SMOKE) 3L   else 10L
S4_NSIMS       <- if (SMOKE) 5L   else 20L
S5_NSIMS       <- if (SMOKE) 5L   else 10L
S_NGENES       <- if (SMOKE) 500L else 5000L

# If run_pipeline.R already ran and saved opts_bio, point here to skip re-estimation.
# Set to NULL to re-estimate from raw data.
PIPELINE1_RDS  <- NULL     # e.g. "output/run_20240101_1200/opts_bio.rds"


# =====================================================================
# SECTION 1: SETUP
# =====================================================================

# Set POWERSIM_ROOT as env var for portability:
#   Linux/Mac shell: export POWERSIM_ROOT=/path/to/PowerSim
POWERSIM_ROOT <- {
  env <- Sys.getenv("POWERSIM_ROOT", unset = "")
  if (nzchar(env)) env else
    "/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim"
}

# DATA_HUMAN path (outside PowerSim root; override via env var on server):
#   export DATA_HUMAN=/path/to/combined_data.rds
DATA_HUMAN <- {
  env <- Sys.getenv("DATA_HUMAN", unset = "")
  if (nzchar(env)) env else
    "/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Collaborative/Paper/Congruence/PNAS_aging/data/combined_data.rds"
}

setwd(POWERSIM_ROOT)
old_wd <- setwd(file.path(getwd(), "code"))
source("setup.R")
setwd(old_wd)

run_tag <- format(Sys.time(), "%Y%m%d_%H%M")
base_out <- file.path("output", paste0("run_pipeline2_", run_tag))
fig_dir  <- file.path(base_out, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("Output: %s/\n\n", base_out))

t_start <- proc.time()

opts_analysis <- CircadianAnalysisOptions(
  fdr_thresholds = c(0.01, 0.05, 0.10, 0.20),
  reference_n    = 60
)


# =====================================================================
# SECTION 2: PILOT DATA & PARAMETER ESTIMATION
# =====================================================================
# Either load pre-saved opts_bio from run_pipeline.R, or re-estimate
# from the raw data. The two arms give identical results; reusing avoids
# running estCircadianParam() twice.

cat("====================================================================\n")
cat("SECTION 2: PILOT DATA\n")
cat("====================================================================\n\n")

if (!is.null(PIPELINE1_RDS) && file.exists(PIPELINE1_RDS)) {

  cat(sprintf("Loading pre-estimated opts_bio from:\n  %s\n\n", PIPELINE1_RDS))
  saved <- readRDS(PIPELINE1_RDS)
  opts_bio_real  <- saved$opts_bio
  expr_younger   <- saved$expr_younger
  times_young    <- saved$times_young
  rm(saved)

} else {

  cat("Loading BA11/BA47 younger brain expression data...\n")

  COMBINED <- readRDS(DATA_HUMAN)

  expr_sample_names <- colnames(COMBINED$expr)
  pheno_order       <- match(expr_sample_names, COMBINED$pheno$sample_name)
  valid_samples     <- !is.na(pheno_order)
  COMBINED$expr     <- COMBINED$expr[, valid_samples]
  pheno_order       <- pheno_order[valid_samples]
  pheno_data        <- COMBINED$pheno[pheno_order, ]
  pheno_data$tod    <- if ("TOD.x" %in% colnames(pheno_data)) pheno_data$TOD.x
                       else pheno_data$TOD.y
  pheno_data$age_group_final <- if ("AgeGroup" %in% colnames(pheno_data)) pheno_data$AgeGroup
                                else pheno_data$age_group
  complete_samples  <- !is.na(pheno_data$age_group_final) & !is.na(pheno_data$tod) &
                       pheno_data$age_group_final %in% c("younger", "older")
  pheno_clean       <- pheno_data[complete_samples, ]
  younger_idx       <- pheno_clean$age_group_final == "younger"
  expr_raw_young    <- COMBINED$expr[, complete_samples][, younger_idx]
  times_raw         <- pheno_clean$tod[younger_idx]

  rm(COMBINED, expr_sample_names, pheno_order, valid_samples,
     pheno_data, complete_samples, pheno_clean, younger_idx)

  # Use prepCircadianData for standard validation (log2 already; coerces to matrix)
  prep_young  <- prepCircadianData(expr_raw_young, times = times_raw, input_type = "log2")
  expr_younger <- prep_young$data
  times_young  <- prep_young$times
  rm(expr_raw_young, times_raw)

  cat(sprintf("  Expression matrix: %d genes x %d samples\n", nrow(expr_younger), ncol(expr_younger)))
  cat(sprintf("  Pilot subjects (younger): n = %d\n", length(times_young)))
  cat(sprintf("  TOD range: %.1f – %.1f h  (unique: %d)\n\n",
              min(times_young), max(times_young), length(unique(times_young))))

  cat("Estimating circadian parameters...\n\n")
  opts_bio_real <- estCircadianParam(
    data       = expr_younger,
    times      = times_young,
    period     = 24,
    prop_DR    = 0.15,
    prop_DP    = 0.10,
    phase_diff = c(-6, 6),
    amp_diff   = c(0.5, 2),
    verbose    = TRUE
  )
  opts_bio_real <- updateBioOptions(opts_bio_real, ngenes = S_NGENES)
}

cat("\nEstimated bio options:\n")
print(opts_bio_real)

# DR-only bio for the design grid (clean single-test comparison)
opts_bio_DR <- updateBioOptions(opts_bio_real,
  prop_DR = 0.15, prop_DP = 0.00,
  phase_diff = c(0, 0), amp_diff = c(1, 1)
)


# =====================================================================
# SECTION 3: BOOTSTRAP DESIGN GRID  [PRIMARY DELIVERABLE]
# =====================================================================
# Scientific question (post-mortem context):
#   "For fixed total subjects, how should I split across TOD bins (B)
#    and subjects per bin (m) to maximise DR detection power?"
#
# Design: passive — subjects drawn from the real pilot TOD distribution.
# This reflects the post-mortem recruitment reality (you cannot place
# subjects at fixed time points; you recruit whoever died at each TOD).
# Each simulation draws N subjects from the empirical TOD distribution,
# producing realistic irregular spacing across the 24h cycle.
#
# The bootstrap outer loop propagates uncertainty from finite pilot size:
# each draw resamples gene rows from the pilot fit, giving a different
# (M, A, phi, sigma) distribution and therefore a different power curve.
# The resulting CI answers: "How stable is this design recommendation
# given that our parameter estimates came from a finite pilot?"
#
# B_values and m_values must satisfy N = B * m exactly (no rounding).
# Choose B_values and m_values so LCM(B_values) divides all N = B * m.

cat("\n====================================================================\n")
cat("SECTION 3: BOOTSTRAP DESIGN GRID  [primary deliverable]\n")
cat("====================================================================\n\n")

# Passive design: subject times are drawn from the real TOD distribution.
# B is not identifiable in passive mode (sampling is independent of B),
# so we sweep over N directly using a single placeholder B = 4.
# The grid covers biologically relevant total sample sizes.
S3_B_VALUE <- 4L
S3_N_VALUES <- c(24L, 36L, 48L, 60L, 72L, 96L, 120L)

cat(sprintf("Design:       passive (real TOD distribution)\n"))
cat(sprintf("N_values:     %s\n", paste(S3_N_VALUES, collapse = ", ")))
cat(sprintf("nboot:        %d\n", S3_NBOOT))
cat(sprintf("nsims_inner:  %d\n\n", S3_NSIMS_INNER))

boot.opts_s3 <- CircadianBootstrapOptions(
  design_vector = times_young,     # real pilot TOD distribution (passive)
  B_values      = S3_B_VALUE,      # single B — passive sweeps N only
  N_values      = S3_N_VALUES,
  nboot         = S3_NBOOT,
  nsims_inner   = S3_NSIMS_INNER,
  design        = "passive",
  seed          = 42L
)
print(boot.opts_s3)

boot_grid <- runBootstrapDesignGrid(
  pilot_data    = expr_younger,
  pilot_times   = times_young,
  boot.opts     = boot.opts_s3,
  analysis.opts = opts_analysis,
  bio_diff.opts = opts_bio_DR,
  verbose       = TRUE
)

saveRDS(boot_grid, file.path(base_out, "s3_boot_grid.rds"))

cat("\n--- Bootstrap Design Grid Summary (DR) ---\n")
summaryBootstrapDesignGrid(boot_grid, test_type = "DR")

plotBootstrapDesignGrid(
  boot_grid,
  test_type   = "DR",
  output_file = file.path(fig_dir, "s3_bootstrap_design_grid.pdf")
)
cat(sprintf("Figure: %s\n", file.path(fig_dir, "s3_bootstrap_design_grid.pdf")))


# =====================================================================
# SECTION 4 [OPTIONAL]: TWO-STAGE VS BOOTSTRAP ON REAL PILOT
# =====================================================================
# Purpose: uncertainty illustration, NOT validation.
#   Two-stage treats the pilot parameter estimates as exact.
#   Bootstrap propagates the finite-pilot uncertainty.
#   The difference between them shows how much the point-estimate
#   recommendation could shift if the pilot were resampled.
#
# This is practically useful for communicating to collaborators:
#   "Our n80 recommendation is X, but given the pilot size, the true
#    recommendation could reasonably be anywhere in [lo, hi]."
#
# Method validation (does bootstrap actually cover oracle truth?) lives
# in code/example.R Scenario 8 — NOT here.
#
# Implementation: single fixed B = modal optimal from Section 3.
# Single-B ensures the CI reflects only parameter uncertainty,
# not variation across incompatible (B, m) design alternatives.

if (RUN_OPTIONAL_S4) {

  cat("\n====================================================================\n")
  cat("SECTION 4 [OPTIONAL]: TWO-STAGE VS BOOTSTRAP COMPARISON\n")
  cat("====================================================================\n\n")
  cat("Interpretation: uncertainty illustration only.\n")
  cat("Validation of the method lives in code/example.R Scenarios 7-8.\n\n")

  # Passive design: use same N grid from Section 3, single B=4 (placeholder)
  S4_B <- 4L
  all_N <- sort(unique(boot_grid$N_values))
  cat(sprintf("Fixed B = %d (passive design — sweeps N only)\n", S4_B))
  cat(sprintf("N grid: %s\n\n", paste(all_N, collapse = ", ")))

  opts_design_s4 <- CircadianDesignOptions(
    sample_sizes = all_N,
    nsims        = S4_NSIMS,
    design       = "passive",
    cts          = times_young
  )

  # Two-stage: estimate once → single power curve
  cat("Running two-stage (point estimate)...\n")
  two_stage_s4 <- runTwoStagePower(
    pilot_data    = expr_younger,
    pilot_times   = times_young,
    design.opts   = opts_design_s4,
    analysis.opts = opts_analysis,
    bio_diff.opts = opts_bio_DR,
    test_type     = "DR",
    verbose       = TRUE
  )

  # Bootstrap: single fixed B → power ± CI (parameter uncertainty only)
  cat(sprintf("\nRunning bootstrap (fixed B=%d, parameter uncertainty)...\n", S4_B))
  boot.opts_s4 <- CircadianBootstrapOptions(
    design_vector = times_young,
    B_values      = S4_B,
    N_values      = all_N,
    nboot         = S3_NBOOT,
    nsims_inner   = S3_NSIMS_INNER,
    design        = "passive",
    seed          = 43L
  )
  boot_s4 <- runBootstrapDesignGrid(
    pilot_data    = expr_younger,
    pilot_times   = times_young,
    boot.opts     = boot.opts_s4,
    analysis.opts = opts_analysis,
    bio_diff.opts = opts_bio_DR,
    verbose       = TRUE
  )

  comparison_s4 <- compareDesignApproaches(
    two_stage_result = two_stage_s4,
    bootstrap_result = boot_s4,
    test_type        = "DR",
    target_power     = 0.80
  )

  saveRDS(list(two_stage  = two_stage_s4,
               boot_s4    = boot_s4,
               comparison = comparison_s4),
          file.path(base_out, "s4_comparison.rds"))

  cat("\n--- Uncertainty Summary ---\n")
  cat(sprintf("  Two-stage n80 (point estimate):  %s\n",
              ifelse(is.na(comparison_s4$n80_two_stage), ">max(N)", comparison_s4$n80_two_stage)))
  cat(sprintf("  Bootstrap n80 median:            %s\n",
              ifelse(is.na(comparison_s4$n80_boot_median), ">max(N)", comparison_s4$n80_boot_median)))
  cat(sprintf("  Bootstrap n80 95%% CI:            [%s, %s]\n",
              ifelse(is.na(comparison_s4$n80_boot_lo), "NA", comparison_s4$n80_boot_lo),
              ifelse(is.na(comparison_s4$n80_boot_hi), "NA", comparison_s4$n80_boot_hi)))
  cat("  (CI width = pilot-estimation uncertainty; method validity → code/example.R S8)\n")

  plotDesignComparison(
    comparison_s4,
    target_power = 0.80,
    output_file  = file.path(fig_dir, "s4_design_comparison.pdf")
  )
  cat(sprintf("Figure: %s\n", file.path(fig_dir, "s4_design_comparison.pdf")))

} else {
  cat("\nSection 4 skipped (RUN_OPTIONAL_S4 = FALSE).\n")
}


# =====================================================================
# SECTION 5 [OPTIONAL]: FOURIER MISSPECIFICATION SENSITIVITY
# =====================================================================
# Purpose: sensitivity analysis on the real-pilot design recommendation.
#   The pilot data informs the base parameter distribution (M, A, phi, sigma).
#   This section asks: if the TRUE waveform also has 2nd/3rd harmonics that
#   the cosinor model ignores, how much would power shift?
#
# Interpretation:
#   - This is NOT a property learned from the pilot data itself.
#     The pilot fit is always a pure cosinor; harmonics are added
#     synthetically at simulation time.
#   - Use this to answer: "Is our design recommendation (from Section 3)
#     still valid if the biology deviates from a pure cosinor?"
#   - Method justification (why this matters) → code/example.R Scenario 6.
#
# The N grid and design (passive, real TOD) match Section 3/4 so the
# sensitivity result is directly comparable to the main recommendation.

if (RUN_OPTIONAL_S5) {

  cat("\n====================================================================\n")
  cat("SECTION 5 [OPTIONAL]: FOURIER MISSPECIFICATION SENSITIVITY\n")
  cat("====================================================================\n\n")
  cat("Interpretation: sensitivity of the Section 3 recommendation to\n")
  cat("waveform misspecification. Method justification → code/example.R S6.\n\n")

  # Use the recommended N grid from Section 3
  N_for_s5 <- sort(unique(boot_grid$N_values))

  opts_design_s5 <- CircadianDesignOptions(
    sample_sizes = N_for_s5,
    nsims        = S5_NSIMS,
    design       = "passive",
    cts          = times_young
  )

  harmonic_grid_s5 <- expand.grid(
    alpha2 = c(0, 0.25, 0.50),
    alpha3 = c(0, 0.25)
  )

  cat(sprintf("Harmonic grid: %d combinations\n", nrow(harmonic_grid_s5)))
  cat(sprintf("N grid: %s\n\n", paste(N_for_s5, collapse = ", ")))

  fourier_s5 <- runFourierDeviationPower(
    bio.opts      = opts_bio_DR,
    design.opts   = opts_design_s5,
    analysis.opts = opts_analysis,
    harmonic_grid = harmonic_grid_s5,
    test_type     = "DR",
    verbose       = TRUE
  )

  saveRDS(fourier_s5, file.path(base_out, "s5_fourier_sensitivity.rds"))

  # Power loss summary at reference N
  ref_col  <- which.min(abs(fourier_s5$sample_sizes - opts_analysis$reference_n))
  pure_idx <- which(fourier_s5$harmonic_grid$alpha2 == 0 &
                    fourier_s5$harmonic_grid$alpha3 == 0)
  harm_idx <- which(fourier_s5$harmonic_grid$alpha2 >= 0.5)

  if (length(pure_idx) > 0 && length(harm_idx) > 0) {
    pure_pwr <- 100 * mean(fourier_s5$power_mean[pure_idx, ref_col], na.rm = TRUE)
    harm_pwr <- 100 * mean(fourier_s5$power_mean[harm_idx, ref_col], na.rm = TRUE)
    cat(sprintf("\n--- Sensitivity at N=%d ---\n", fourier_s5$sample_sizes[ref_col]))
    cat(sprintf("  Pure cosinor DR power:         %.1f%%\n", pure_pwr))
    cat(sprintf("  High-harmonic DR power (a2≥.5): %.1f%%\n", harm_pwr))
    cat(sprintf("  Absolute power loss:           %.1f pp\n", pure_pwr - harm_pwr))
    cat(sprintf("  Recommendation change: %s\n",
                ifelse(abs(pure_pwr - harm_pwr) < 5,
                       "negligible (<5pp) — design recommendation robust",
                       "notable (≥5pp) — consider larger N under harmonic misspecification")))
  }

  plotFourierDeviation(
    fourier_s5,
    test_type   = "DR",
    output_file = file.path(fig_dir, "s5_fourier_sensitivity.pdf")
  )
  cat(sprintf("Figure: %s\n", file.path(fig_dir, "s5_fourier_sensitivity.pdf")))

} else {
  cat("\nSection 5 skipped (RUN_OPTIONAL_S5 = FALSE).\n")
}


# =====================================================================
# SECTION 6: WRAP-UP
# =====================================================================

t_elapsed <- proc.time() - t_start
cat("\n====================================================================\n")
cat("PIPELINE 2 COMPLETE\n")
cat("====================================================================\n\n")
cat(sprintf("Runtime: %.1f min\n\n", t_elapsed[3] / 60))

cat("Saved outputs:\n")
for (f in list.files(base_out, recursive = TRUE, full.names = TRUE)) {
  cat(sprintf("  %s\n", f))
}

cat("\n--- Design recommendation summary ---\n")
cat("Section 3 (bootstrap design grid):\n")
for (ni in seq_along(boot_grid$N_values)) {
  N   <- boot_grid$N_values[ni]
  opt <- boot_grid$optimal_B[ni]
  lo  <- boot_grid$optimal_B_ci_lo[ni]
  hi  <- boot_grid$optimal_B_ci_hi[ni]
  cat(sprintf("  N=%d: optimal B = %d  (bootstrap range [%d, %d])\n", N, opt, lo, hi))
}

if (RUN_OPTIONAL_S4 && exists("comparison_s4")) {
  cat("\nSection 4 (pilot-estimation uncertainty):\n")
  cat(sprintf("  n80 two-stage:  %s\n",
              ifelse(is.na(comparison_s4$n80_two_stage), ">max(N)", comparison_s4$n80_two_stage)))
  cat(sprintf("  n80 bootstrap:  %s  [%s, %s]\n",
              ifelse(is.na(comparison_s4$n80_boot_median), ">max(N)", comparison_s4$n80_boot_median),
              ifelse(is.na(comparison_s4$n80_boot_lo), "NA", comparison_s4$n80_boot_lo),
              ifelse(is.na(comparison_s4$n80_boot_hi), "NA", comparison_s4$n80_boot_hi)))
}

cat("\n")
