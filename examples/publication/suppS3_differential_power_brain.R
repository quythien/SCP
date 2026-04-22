#' =======================================================================
#' Differential Circadian Power Analysis — Figure 2
#' =======================================================================
#'
#' WHAT THIS DOES:
#'   Simulation-based power analysis for transcriptome-wide differential
#'   circadian detection (DR, DP, DM) across two comparisons:
#'     (A) NAc vs Putamen (Ctrl-Ctrl): strong DR + DP + DM, balanced n=59/59
#'     (B) Putamen Ctrl vs SCZ: disease context, DR-dominant, unbalanced n=59/28
#'
#' FIGURE STRUCTURE:
#'   18 panels: 6 rows x 3 columns
#'   Row 1: Comp A — DR (marginal power | TD by r | stratified power)
#'   Row 2: Comp B — DR
#'   Row 3: Comp A — DP
#'   Row 4: Comp B — DP
#'   Row 5: Comp A — DM
#'   Row 6: Comp B — DM
#'
#' INTERMEDIATE OUTPUTS (saved before figure generation):
#'   output/differential/results/diff_power_NAc_vs_Putamen_Ctrl_<ts>.rds
#'   output/differential/results/diff_power_Putamen_Ctrl_vs_SCZ_<ts>.rds
#'
#' FINAL OUTPUT:
#'   output/differential/figures/diff_power_fig2_<ts>.pdf
#'
#' USAGE:
#'   cd /path/to/PowerSim
#'   Rscript examples/publication/11_differential_power.R
#'
#'   # To regenerate figure only from saved results:
#'   RESULTS_RDS_A=output/differential/results/diff_power_NAc_vs_Putamen_Ctrl_<ts>.rds \
#'   RESULTS_RDS_B=output/differential/results/diff_power_Putamen_Ctrl_vs_SCZ_<ts>.rds \
#'   Rscript examples/publication/11_differential_power.R
#'
#' @author Thien Pham

# =====================================================================
# 0. Setup
# =====================================================================
SMOKE_TEST  <- identical(Sys.getenv("SMOKE_TEST"), "true")
GLOBAL_SEED <- 2025L
set.seed(GLOBAL_SEED)

old_wd <- setwd("code")
source("setup.R")
setwd(old_wd)

out_dir_fig <- "output/differential/figures"
out_dir_res <- "output/differential/results"
dir.create(out_dir_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_res, recursive = TRUE, showWarnings = FALSE)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
cat(sprintf("\n=== Differential Circadian Power Analysis  [%s] ===\n", timestamp))

# =====================================================================
# 1. Simulation grid
# =====================================================================
if (SMOKE_TEST) {
  sample_sizes <- c(30, 60, 100)
  nsims        <- 5L
  n_cores      <- as.integer(Sys.getenv("MC_CORES", unset = "1"))
} else {
  sample_sizes <- c(20, 30, 40, 50, 60, 80, 100, 120, 150, 200)
  nsims        <- 200L
  n_cores      <- as.integer(Sys.getenv("MC_CORES", unset = "4"))
}

cat(sprintf("Global seed  : %d\n", GLOBAL_SEED))
cat(sprintf("Sample sizes : %s\n", paste(sample_sizes, collapse = ", ")))
cat(sprintf("nsims        : %d\n", nsims))
cat(sprintf("mc.cores     : %d\n", n_cores))

# =====================================================================
# 2. Comparisons
# =====================================================================
comparisons <- list(
  NAc_vs_Putamen_Ctrl = list(
    pilot_path  = "data/gse160521_nac_vs_putamen_ctrl_pilot.rds",
    label       = "NAc vs Putamen (Ctrl-Ctrl)",
    label_short = "NAc_vs_Putamen_Ctrl"
  ),
  Putamen_Ctrl_vs_SCZ = list(
    pilot_path  = "data/gse160521_putamen_ctrl_vs_scz_pilot.rds",
    label       = "Putamen Ctrl vs SCZ",
    label_short = "Putamen_Ctrl_vs_SCZ"
  )
)

# =====================================================================
# 3. Run simulations (or load existing results)
# =====================================================================

# Allow pre-computed results to be passed via env vars for figure-only re-runs
env_rds <- list(
  NAc_vs_Putamen_Ctrl = Sys.getenv("RESULTS_RDS_A", unset = ""),
  Putamen_Ctrl_vs_SCZ = Sys.getenv("RESULTS_RDS_B", unset = "")
)

results <- list()

comp_filter    <- Sys.getenv("COMP", unset = "")
comp_names_run <- if (nchar(comp_filter) > 0 && comp_filter %in% names(comparisons))
                    comp_filter else names(comparisons)

for (comp_name in comp_names_run) {
  comp <- comparisons[[comp_name]]

  # --- Check for pre-existing results ---
  rds_override <- env_rds[[comp_name]]
  if (nchar(rds_override) > 0 && file.exists(rds_override)) {
    cat(sprintf("\n[%s] Loading pre-computed results: %s\n", comp$label, rds_override))
    results[[comp_name]] <- readRDS(rds_override)
    next
  }

  if (!file.exists(comp$pilot_path)) {
    warning(sprintf("Pilot not found, skipping: %s", comp$pilot_path)); next
  }

  bio <- readRDS(comp$pilot_path)
  cat(sprintf("\n====================================================\n"))
  cat(sprintf("  %s\n", comp$label))
  cat(sprintf("====================================================\n"))
  cat(sprintf("Pilot: n1=%d  n2=%d  DR=%.3f  DP=%.3f  DM=%.3f\n",
              length(bio$cts), length(bio$cts2),
              bio$prop_DR, bio$prop_DP, bio$prop_DM))

  design <- CircadianDesignOptions(
    sample_sizes = sample_sizes,
    nsims        = nsims,
    design       = "passive",
    cts          = bio$cts,
    test_types   = c("DR", "DP", "DM")
  )

  analysis <- CircadianAnalysisOptions(
    alpha           = 0.05,
    p.adjust.method = "BH"
  )

  # --- Save intermediate results immediately after run ---
  rds_path <- file.path(out_dir_res,
    sprintf("diff_power_%s_%s.rds", comp$label_short, timestamp))

  set.seed(GLOBAL_SEED)
  res <- runDifferentialPower(bio, design, analysis,
                               methods    = "DCP",
                               test_types = c("DR", "DP", "DM"),
                               plot       = FALSE,
                               verbose    = TRUE,
                               mc.cores   = n_cores)

  saveRDS(res, rds_path)
  cat(sprintf("Results saved -> %s\n", rds_path))

  results[[comp_name]] <- res
}

# =====================================================================
# 4. Generate Figure 2
# =====================================================================
if (length(results) < 2) {
  cat(sprintf("\nOnly %d comparison(s) run — skipping Figure 2 (need both).\n",
              length(results)))
  cat(sprintf("Results saved to: %s\n", out_dir_res))
  quit(save = "no", status = 0)
}

fig_path <- file.path(out_dir_fig, sprintf("diff_power_fig2_%s.pdf", timestamp))
cat(sprintf("\nGenerating Figure 2 -> %s\n", fig_path))

plotDiffPower(
  res_list   = results,
  comp_labels = c(
    comparisons$NAc_vs_Putamen_Ctrl$label,
    comparisons$Putamen_Ctrl_vs_SCZ$label
  ),
  out_pdf    = fig_path,
  width      = 15,
  height     = 30
)

cat("\n=== Done ===\n")
cat(sprintf("Results: %s\n", out_dir_res))
cat(sprintf("Figures: %s\n", out_dir_fig))
