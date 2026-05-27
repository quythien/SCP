#' =======================================================================
#' Single-Cohort Circadian Power Analysis — Figure 1
#' =======================================================================
#'
#' WHAT THIS DOES:
#'   Simulation-based power analysis for transcriptome-wide circadian
#'   rhythmicity detection in a single cohort (one-group setting).
#'   Uses the SCP framework: pilot-calibrated semi-parametric simulation
#'   with r = A/sigma stratification and BH correction.
#'
#' SCIENTIFIC QUESTION:
#'   "How many post-mortem brain samples are needed to detect
#'    circadian gene expression with adequate power, given the
#'    empirical effect-size distribution from the Seney human PFC data?"
#'
#' FRAMEWORK:
#'   - Parametric cosinor model: y = M + A*cos(omega*t - phi) + epsilon
#'   - Pilot: Seney Human PFC (ACC), Control group, n=60, 14,455 genes
#'   - Passive design: sampling times from pilot TOD distribution
#'   - BH-FDR correction across all genes; nominal alpha = 0.05
#'   - r-stratified power curves (r = A/sigma SNR strata)
#'   - CircaPower closed-form grid initialisation
#'
#' OUTPUTS:
#'   figures/single_cohort_power.pdf — 3 panels (marginal, discoveries, stratified)
#'   results/single_cohort_power.rds — full simulation results
#'
#' PILOT DATA:
#'   data/gse160521_nac_ctrl_pilot.rds — pre-estimated CircadianBioOptions
#'   Source: GSE160521, Human NAc (nucleus accumbens), Control group
#'           n=59, 15,330 genes, r range 0.43-2.01, r_med=0.55
#'   Metadata: Kyle_multiBrainRegion/NAc_clinical_*_matchIndex34.csv
#'
#' USAGE:
#'   cd /path/to/PowerSim
#'   Rscript examples/publication/10_single_cohort_power.R
#'
#' @author Thien Quy Pham

# =====================================================================
# 0. Setup
# =====================================================================
SMOKE_TEST <- identical(Sys.getenv("SMOKE_TEST"), "true")
GLOBAL_SEED <- 2025L
set.seed(GLOBAL_SEED)

old_wd <- setwd("code")
source("setup.R")
setwd(old_wd)

out_dir_fig <- "output/single_cohort/figures"
out_dir_res <- "output/single_cohort/results"
dir.create(out_dir_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_res, recursive = TRUE, showWarnings = FALSE)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
cat(sprintf("\n=== Single-Cohort Power Analysis  [%s] ===\n", timestamp))

# =====================================================================
# 1. Simulation grid
# =====================================================================
if (SMOKE_TEST) {
  sample_sizes <- c(20, 40, 60, 80, 100)
  nsims        <- 5L
  n_cores      <- 1L
} else {
  sample_sizes <- c(20, 30, 40, 50, 60, 70, 80, 100, 120, 150, 200)
  nsims        <- 200L
  n_cores      <- as.integer(Sys.getenv("MC_CORES", unset = "4"))
}

cat(sprintf("Global seed  : %d\n", GLOBAL_SEED))
cat(sprintf("Sample sizes : %s\n", paste(sample_sizes, collapse = ", ")))
cat(sprintf("nsims        : %d\n", nsims))
cat(sprintf("mc.cores     : %d\n", n_cores))
cat(sprintf("alpha (BH)   : 0.05\n"))
cat(sprintf("design       : passive (pilot TOD distribution)\n"))
cat(sprintf("bin_width    : 0.25\n\n"))

# =====================================================================
# 2. Loop over brain regions
# =====================================================================
regions_all <- list(
  NAc     = "data/gse160521_nac_ctrl_pilot.rds",
  Caudate = "data/gse160521_caudate_ctrl_pilot.rds",
  Putamen = "data/gse160521_putamen_ctrl_pilot.rds"
)
region_filter <- Sys.getenv("REGION", unset = "")
regions <- if (nchar(region_filter) > 0 && region_filter %in% names(regions_all))
             regions_all[region_filter] else regions_all

for (region_name in names(regions)) {

  pilot_path <- regions[[region_name]]
  if (!file.exists(pilot_path)) {
    warning("Pilot not found, skipping: ", pilot_path); next
  }

  cat(sprintf("\n====================================================\n"))
  cat(sprintf("  Region: %s\n", region_name))
  cat(sprintf("====================================================\n"))

  bio <- readRDS(pilot_path)
  cat(sprintf("Pilot: n_genes=%d  prop_rhythmic=%.1f%%  r_med=%.3f  n_pilot=%d\n",
              bio$ngenes,
              100 * bio$prop_rhythmic,
              median(bio$amplitude / bio$sigma_rhythmic, na.rm = TRUE),
              length(bio$cts)))

  design <- CircadianDesignOptions(
    sample_sizes = sample_sizes,
    nsims        = nsims,
    design       = "passive",
    cts          = bio$cts
  )

  analysis <- CircadianAnalysisOptions(
    alpha           = 0.05,
    p.adjust.method = "BH",
    r_strata        = makeAdaptiveRStrata(bio, bin_width = 0.25)
  )

  cat(sprintf("r_strata: %s\n",
              paste(round(analysis$r_strata[is.finite(analysis$r_strata)], 2),
                    collapse = ", "), "Inf"))

  # --- Run simulation ---
  fig_path <- file.path(out_dir_fig,
    sprintf("single_cohort_power_GSE160521_%s_Control.pdf", region_name))
  rds_path <- file.path(out_dir_res,
    sprintf("single_cohort_power_GSE160521_%s_Control_%s.rds",
            region_name, timestamp))

  set.seed(GLOBAL_SEED)
  res <- runSingleCohortPower(bio, design, analysis,
                               methods  = "DCP",
                               plot     = FALSE,
                               verbose  = TRUE,
                               mc.cores = n_cores)

  saveRDS(res, rds_path)
  cat(sprintf("Results saved → %s\n", rds_path))

  plotSingleCohortPower(
    res     = res,
    out_pdf = fig_path,
    title   = sprintf("GSE160521 %s Control — Single-Cohort Power", region_name)
  )
  cat(sprintf("Figure  saved → %s\n", fig_path))
}

cat("\n=== All regions done ===\n")
