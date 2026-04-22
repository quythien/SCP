#' =======================================================================
#' replot_fig1_fig2.R — Regenerate Figure 1 and Figure 2 from saved RDS
#' =======================================================================
#'
#' Loads the most-recent saved simulation RDS files and replots using
#' plotSingleCohortPower() and plotDiffPower() — both of which recompute
#' power at all FDR thresholds directly from the raw p-values stored in
#' the RDS. The old plotAllStratifiedPower() used fabricated scalar
#' multipliers (x0.7, x1.15, x1.25) for non-5% FDR curves; these
#' functions do not.
#'
#' FIGURE 1 (GTEx Adrenal Gland + Liver, single cohort):
#'   output/single_cohort/figures/fig1_AdrenalGland_fixed.pdf
#'   output/single_cohort/figures/fig1_Liver_fixed.pdf
#'
#' FIGURE 2 (NAc vs Putamen, Putamen Ctrl vs SCZ, differential):
#'   output/differential/figures/fig2_fixed.pdf
#'
#' USAGE:
#'   cd /path/to/PowerSim
#'   Rscript examples/publication/replot_fig1_fig2.R
#'
#'   # To override which RDS files are used:
#'   RDS_ADR=output/.../...rds RDS_LIV=... RDS_DIFF_A=... RDS_DIFF_B=... \
#'   Rscript examples/publication/replot_fig1_fig2.R

old_wd <- setwd("code"); source("setup.R"); setwd(old_wd)

out_fig_sc   <- "output/single_cohort/figures"
out_fig_diff <- "output/differential/figures"
dir.create(out_fig_sc,   recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig_diff, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------
# Helper: pick most-recent RDS matching a glob, with env-var override
# -----------------------------------------------------------------------
pick_rds <- function(env_var, glob_dir, glob_pattern) {
  override <- Sys.getenv(env_var, unset = "")
  if (nchar(override) > 0) {
    if (!file.exists(override)) stop(sprintf("File from %s not found: %s", env_var, override))
    return(override)
  }
  candidates <- list.files(glob_dir, pattern = glob_pattern, full.names = TRUE)
  if (length(candidates) == 0)
    stop(sprintf("No RDS found matching '%s' in %s", glob_pattern, glob_dir))
  candidates[order(file.mtime(candidates), decreasing = TRUE)[1]]
}

# -----------------------------------------------------------------------
# Figure 1A: GTEx Adrenal Gland
# -----------------------------------------------------------------------
cat("\n=== Figure 1A: GTEx Adrenal Gland ===\n")
rds_adr <- pick_rds(
  "RDS_ADR",
  "output/single_cohort/results",
  "single_cohort_power_GTEx_AdrenalGland.*\\.rds$"
)
cat(sprintf("  RDS: %s\n", basename(rds_adr)))
res_adr <- readRDS(rds_adr)

out_adr <- file.path(out_fig_sc, "fig1_AdrenalGland_fixed.pdf")
plotSingleCohortPower(
  res             = res_adr,
  out_pdf         = out_adr,
  title           = "GTEx Adrenal Gland - Single-Cohort Power",
  fdr_thresholds  = c(0.01, 0.05, 0.10, 0.20),
  panel_fdr       = 0.05,
  vline_power     = 0.80,
  vline_fdr       = 0.05,
  reference_n     = NULL,
  display_sizes   = NULL,
  width           = 15,
  height          = 5.5
)
cat(sprintf("  Saved: %s\n", out_adr))

# -----------------------------------------------------------------------
# Figure 1B: GTEx Liver
# -----------------------------------------------------------------------
cat("\n=== Figure 1B: GTEx Liver ===\n")
rds_liv <- pick_rds(
  "RDS_LIV",
  "output/single_cohort/results",
  "single_cohort_power_GTEx_Liver.*\\.rds$"
)
cat(sprintf("  RDS: %s\n", basename(rds_liv)))
res_liv <- readRDS(rds_liv)

out_liv <- file.path(out_fig_sc, "fig1_Liver_fixed.pdf")
plotSingleCohortPower(
  res             = res_liv,
  out_pdf         = out_liv,
  title           = "GTEx Liver - Single-Cohort Power",
  fdr_thresholds  = c(0.01, 0.05, 0.10, 0.20),
  panel_fdr       = 0.05,
  vline_power     = 0.80,
  vline_fdr       = 0.05,
  reference_n     = NULL,
  display_sizes   = NULL,
  width           = 15,
  height          = 5.5
)
cat(sprintf("  Saved: %s\n", out_liv))

# -----------------------------------------------------------------------
# Figure 2: Differential power
# -----------------------------------------------------------------------
cat("\n=== Figure 2: Differential power ===\n")
rds_diff_a <- pick_rds(
  "RDS_DIFF_A",
  "output/differential/results",
  "diff_power_NAc_vs_Putamen_Ctrl.*\\.rds$"
)
rds_diff_b <- pick_rds(
  "RDS_DIFF_B",
  "output/differential/results",
  "diff_power_Putamen_Ctrl_vs_SCZ.*\\.rds$"
)
cat(sprintf("  RDS A: %s\n", basename(rds_diff_a)))
cat(sprintf("  RDS B: %s\n", basename(rds_diff_b)))

res_diff_a <- readRDS(rds_diff_a)
res_diff_b <- readRDS(rds_diff_b)

out_fig2 <- file.path(out_fig_diff, "fig2_fixed.pdf")
plotDiffPower(
  res_list       = list(res_diff_a, res_diff_b),
  comp_labels    = c("NAc vs Putamen (Ctrl-Ctrl)", "Putamen Ctrl vs SCZ"),
  endpoints      = c("DR", "DP", "DM"),
  fdr_thresholds = c(0.01, 0.05, 0.10, 0.20),
  panel_fdr      = 0.05,
  vline_power    = 0.80,
  vline_fdr      = 0.05,
  out_pdf        = out_fig2,
  width          = 15,
  height         = 30
)
cat(sprintf("  Saved: %s\n", out_fig2))

cat("\n=== Done ===\n")
cat(sprintf("Fig 1A : %s\n", out_adr))
cat(sprintf("Fig 1B : %s\n", out_liv))
cat(sprintf("Fig 2  : %s\n", out_fig2))
