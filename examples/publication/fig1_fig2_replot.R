#' =======================================================================
#' fig1_fig2_replot.R - Regenerate Fig 1A, Fig 1B, Fig 2 from cached RDS
#' =======================================================================
#'
#' Default behavior: Fig 1 marginal-power Panel A is computed from the
#' alternate single-cohort RDS, while Panels B and C continue to be
#' computed from the primary cached RDS. Fig 2 uses the trimmed
#' 1 x n_endpoints layout (Panel A only) by default.
#'
#' USAGE: Rscript examples/publication/fig1_fig2_replot.R

setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
old_wd <- setwd("code"); source("setup.R"); setwd(old_wd)

pick_latest <- function(dir, pattern) {
  cands <- list.files(dir, pattern, full.names = TRUE)
  if (length(cands) == 0) stop("No match for ", pattern, " in ", dir)
  cands[order(file.mtime(cands), decreasing = TRUE)[1]]
}

mirror_pdf <- function(src, dest_paths) {
  for (d in dest_paths) {
    dir.create(dirname(d), recursive = TRUE, showWarnings = FALSE)
    file.copy(src, d, overwrite = TRUE)
    cat(sprintf("Saved: %s\n", d))
  }
}

DISPLAY_N <- c(20, 40, 60, 80, 100, 120, 150, 200)
R_MAX     <- 3

# ---------------------------------------------------------------------------
# Fig 1A: GTEx Adrenal Gland
# ---------------------------------------------------------------------------
cat("\n=== Fig 1A: GTEx Adrenal Gland ===\n")
rds_adr_primary <- pick_latest("output/single_cohort/results",
                       "^single_cohort_power_GTEx_AdrenalGland_2026.*\\.rds$")
rds_adr_panelA  <- tryCatch(
  pick_latest("output/single_cohort/results",
              "^single_cohort_power_GTEx_AdrenalGland_paired_.*\\.rds$"),
  error = function(e) NULL
)
cat(sprintf("  Panel B/C: %s\n", basename(rds_adr_primary)))
if (!is.null(rds_adr_panelA))
  cat(sprintf("  Panel A  : %s\n", basename(rds_adr_panelA)))

res_adr_primary <- readRDS(rds_adr_primary)
res_adr_panelA  <- if (!is.null(rds_adr_panelA)) readRDS(rds_adr_panelA) else NULL

tmp_adr <- tempfile(fileext = ".pdf")
plotSingleCohortPower(
  res            = res_adr_primary,
  panel_a_res    = res_adr_panelA,
  out_pdf        = tmp_adr,
  title          = "Single-Cohort Power Analysis (GTEx Adrenal Gland)",
  fdr_thresholds = c(0.01, 0.05, 0.10, 0.20),
  panel_fdr      = 0.05,
  vline_power    = 0.80,
  vline_fdr      = 0.05,
  display_sizes  = DISPLAY_N,
  r_max          = R_MAX,
  width          = 13.5, height = 4.8
)
mirror_pdf(tmp_adr,
           c("output/main_figures/Fig1A_single_cohort_AdrenalGland.pdf",
             "submission/figures/Fig1A_single_cohort_AdrenalGland.pdf"))

# ---------------------------------------------------------------------------
# Fig 1B: GTEx Liver
# ---------------------------------------------------------------------------
cat("\n=== Fig 1B: GTEx Liver ===\n")
rds_liv_primary <- pick_latest("output/single_cohort/results",
                       "^single_cohort_power_GTEx_Liver_2026.*\\.rds$")
rds_liv_panelA  <- tryCatch(
  pick_latest("output/single_cohort/results",
              "^single_cohort_power_GTEx_Liver_paired_.*\\.rds$"),
  error = function(e) NULL
)
cat(sprintf("  Panel B/C: %s\n", basename(rds_liv_primary)))
if (!is.null(rds_liv_panelA))
  cat(sprintf("  Panel A  : %s\n", basename(rds_liv_panelA)))

res_liv_primary <- readRDS(rds_liv_primary)
res_liv_panelA  <- if (!is.null(rds_liv_panelA)) readRDS(rds_liv_panelA) else NULL

tmp_liv <- tempfile(fileext = ".pdf")
plotSingleCohortPower(
  res            = res_liv_primary,
  panel_a_res    = res_liv_panelA,
  out_pdf        = tmp_liv,
  title          = "Single-Cohort Power Analysis (GTEx Liver)",
  fdr_thresholds = c(0.01, 0.05, 0.10, 0.20),
  panel_fdr      = 0.05,
  vline_power    = 0.80,
  vline_fdr      = 0.05,
  display_sizes  = DISPLAY_N,
  r_max          = R_MAX,
  width          = 13.5, height = 4.8
)
mirror_pdf(tmp_liv,
           c("output/main_figures/Fig1B_single_cohort_Liver.pdf",
             "submission/figures/Fig1B_single_cohort_Liver.pdf"))

# ---------------------------------------------------------------------------
# Fig 2: ADR vs LIV differential (Panel A only, 1 x 3 layout)
# ---------------------------------------------------------------------------
cat("\n=== Fig 2: ADR vs LIV differential ===\n")
# Prefer alternate Fig 2 RDS if present, otherwise fall back to primary cache
rds_diff <- tryCatch(
  pick_latest("output/differential/results",
              "^diff_power_ADR_vs_LIV_paired_.*\\.rds$"),
  error = function(e) pick_latest("output/differential/results",
                                  "^diff_power_ADR_vs_LIV_2026.*\\.rds$")
)
cat(sprintf("  RDS: %s\n", basename(rds_diff)))
res_diff <- readRDS(rds_diff)

tmp_diff <- tempfile(fileext = ".pdf")
plotDiffPower(
  res_list       = list(res_diff),
  comp_labels    = c("GTEx Adrenal Gland vs Liver"),
  endpoints      = c("DR", "DP", "DM"),
  fdr_thresholds = c(0.01, 0.05, 0.10, 0.20),
  panel_fdr      = 0.05,
  vline_power    = 0.80,
  # N80 vline anchored at FDR = 0.20 because DR/DM do not reach 80% power
  # at FDR = 0.05 within the displayed N grid.
  vline_fdr      = 0.20,
  display_sizes  = DISPLAY_N,
  main_title     = "Differential Circadian Power Analysis (GTEx Adrenal Gland vs Liver)",
  out_pdf        = tmp_diff
)
mirror_pdf(tmp_diff,
           c("output/main_figures/Fig2_differential_ADR_vs_LIV.pdf",
             "submission/figures/Fig2_differential_ADR_vs_LIV.pdf"))

cat("\n=== Done ===\n")
