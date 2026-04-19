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
#'   data/seney_ctrl_pilot.rds  — pre-estimated CircadianBioOptions
#'   Source: Seney et al., human ACC RNA-seq, Disease=1 (control)
#'           n=60, TOD via official time-of-death (OFFC_TIME)
#'
#' USAGE:
#'   cd /path/to/PowerSim
#'   Rscript examples/publication/10_single_cohort_power.R
#'
#' @author Thien Pham

# =====================================================================
# 0. Setup
# =====================================================================
SMOKE_TEST <- identical(Sys.getenv("SMOKE_TEST"), "true")

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
# 1. Load pilot data
# =====================================================================
pilot_path <- "data/seney_ctrl_pilot.rds"
if (!file.exists(pilot_path)) {
  stop("Pilot file not found: ", pilot_path,
       "\nRun data/estimate_pilot_params.R first, or see session notes.")
}

bio <- readRDS(pilot_path)
cat(sprintf("Pilot: n_genes=%d  prop_rhythmic=%.1f%%  r_med=%.3f\n",
            bio$ngenes,
            100 * bio$prop_rhythmic,
            median(bio$amplitude / bio$sigma_rhythmic, na.rm = TRUE)))

# =====================================================================
# 2. Design & analysis options
# =====================================================================
if (SMOKE_TEST) {
  sample_sizes <- c(20, 40, 60)
  nsims        <- 5L
  n_cores      <- 1L
} else {
  sample_sizes <- c(20, 30, 40, 50, 60, 80, 100, 120, 150, 200)
  nsims        <- 200L
  n_cores      <- as.integer(Sys.getenv("MC_CORES", unset = "4"))
}

design <- CircadianDesignOptions(
  sample_sizes = sample_sizes,
  nsims        = nsims,
  design       = "passive",
  cts          = bio$cts    # pilot TOD distribution for passive sampling
)

analysis <- CircadianAnalysisOptions(
  alpha          = 0.05,
  p.adjust.method = "BH",
  r_strata       = c(0, 0.25, 0.5, 0.75, 1.0, 1.5, Inf)
)

cat(sprintf("Sample sizes: %s\n", paste(sample_sizes, collapse = ", ")))
cat(sprintf("nsims = %d   mc.cores = %d\n", nsims, n_cores))

# =====================================================================
# 3. Run simulation
# =====================================================================
cat("\n--- Running runSimsSingleCohort ---\n")
t_start <- proc.time()

res <- runSimsSingleCohort(bio, design, analysis,
                           verbose = TRUE, mc.cores = n_cores)

elapsed <- (proc.time() - t_start)[["elapsed"]]
cat(sprintf("\nDone in %.1f seconds.\n", elapsed))
cat(sprintf("CircaPower n80 estimate (median r=%.3f): n0=%s\n",
            median(bio$amplitude / bio$sigma_rhythmic, na.rm = TRUE),
            ifelse(is.na(res$n0_circapower), "NA", as.character(res$n0_circapower))))

# =====================================================================
# 4. Console summary
# =====================================================================
cat("\n--- Marginal power summary ---\n")
cat(sprintf("%-6s  %-8s  %-8s  %-8s\n", "n", "Power", "FDR", "Avg TD"))
for (j in seq_along(sample_sizes)) {
  cat(sprintf("%-6d  %6.1f%%  %6.4f  %7.1f\n",
              sample_sizes[j],
              100 * mean(res$marginal_power[j, ], na.rm = TRUE),
              mean(res$marginal_FDR[j, ], na.rm = TRUE),
              mean(res$marginal_TD[j, ], na.rm = TRUE)))
}

# =====================================================================
# 5. Figure 1 — 3 panels
# =====================================================================
fig_path <- file.path(out_dir_fig, "single_cohort_power.pdf")
cat(sprintf("\nGenerating Figure 1 → %s\n", fig_path))

plotSingleCohortFig1(
  res          = res,
  bio.opts     = bio,
  out_pdf      = fig_path,
  title        = "Single-cohort rhythmicity power (Seney PFC Control, n0=60)",
  alpha        = 0.05,
  strata_to_show = which(res$strata_labels %in%
                           c("(0,0.25]", "(0.25,0.5]", "(0.5,0.75]",
                             "(0.75,1]", "(1,1.5]")),
  width  = 12,
  height = 4.5
)

# =====================================================================
# 6. Save results
# =====================================================================
rds_path <- file.path(out_dir_res,
                      sprintf("single_cohort_power_%s.rds", timestamp))
saveRDS(res, rds_path)
cat(sprintf("Results saved → %s\n", rds_path))

cat("\n=== Done ===\n")
