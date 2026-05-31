# =====================================================================
# Fig 1 (1A AdrenalGland, 1B Liver) + Fig 2 (ADR vs LIV) replot from cache.
#
# Encodes the advisor's sizing preferences (May 2026 review):
#  - Designed narrow (Fig1 13x5.05", Fig2 12x4.7") so the manuscript scales
#    them down LESS -> bigger on-page elements; Fig1 panels share Fig2's
#    near-square plot-box ratio.
#  - Titles sized for CONSISTENT on-page height at \textwidth: title cex/width
#    ~ 0.12 (Fig1 title_cex 1.60 @ 13"; Fig2 outer-title cex 1.50 @ 12").
#  - Title sits close to the panels (small top white space); no clipping.
#  - Legends bottom-right, compact (FDR key + 2-column n-key) so they clear
#    the curves; longer SE caps (se_cap 0.55); extra gap from r~ axis title to
#    the rotated stratum ticks.
# Run: Rscript examples/publication/fig1_fig2_replot.R
# =====================================================================
setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
suppressMessages(library(SCP))
source("R/plot_single_cohort.R")
source("R/plot_diff.R")

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
R_MAX <- 3

# Shared sizing for the enlarged 3-panel single-cohort figures (Fig 1A/1B).
SC_ARGS <- list(
  fdr_thresholds = c(0.01, 0.05, 0.10, 0.20),
  panel_fdr = 0.05, vline_power = 0.80, vline_fdr = 0.05,
  display_sizes = DISPLAY_N, r_max = R_MAX,
  cex_main = 1.80, cex_lab = 1.70, cex_axis = 1.55,
  line_lwd = 2.8, pt_cex = 1.35, legend_cex = 1.05,
  title_cex = 1.60, oma_top = 2.1,
  legend_pos = "bottomright", se_cap = 0.55,
  width = 13.0, height = 5.05
)

draw_sc <- function(primary_rx, panelA_rx, title, dests) {
  rds_p <- pick_latest("output/single_cohort/results", primary_rx)
  rds_a <- pick_latest("output/single_cohort/results", panelA_rx)
  cat(sprintf("  B/C: %s\n  A  : %s\n", basename(rds_p), basename(rds_a)))
  tmp <- tempfile(fileext = ".pdf")
  do.call(plotSingleCohortPower,
          c(list(res = readRDS(rds_p), panel_a_res = readRDS(rds_a),
                 out_pdf = tmp, title = title), SC_ARGS))
  mirror_pdf(tmp, dests)
}

cat("\n=== Fig 1A: GTEx Adrenal Gland ===\n")
draw_sc("^single_cohort_power_GTEx_AdrenalGland_2026.*\\.rds$",
        "^single_cohort_power_GTEx_AdrenalGland_paired_.*\\.rds$",
        "Single-Cohort Power Analysis (GTEx Adrenal Gland)",
        c("figures/Fig1A_single_cohort_AdrenalGland.pdf",
          "submission/figures/Fig1A_single_cohort_AdrenalGland.pdf"))

cat("\n=== Fig 1B: GTEx Liver ===\n")
draw_sc("^single_cohort_power_GTEx_Liver_2026.*\\.rds$",
        "^single_cohort_power_GTEx_Liver_paired_.*\\.rds$",
        "Single-Cohort Power Analysis (GTEx Liver)",
        c("figures/Fig1B_single_cohort_Liver.pdf",
          "submission/figures/Fig1B_single_cohort_Liver.pdf"))

cat("\n=== Fig 2: ADR vs LIV differential ===\n")
rds_diff <- tryCatch(
  pick_latest("output/differential/results", "^diff_power_ADR_vs_LIV_paired_.*\\.rds$"),
  error = function(e) pick_latest("output/differential/results",
                                  "^diff_power_ADR_vs_LIV_2026.*\\.rds$"))
cat(sprintf("  RDS: %s\n", basename(rds_diff)))
tmp_diff <- tempfile(fileext = ".pdf")
plotDiffPower(
  res_list       = list(readRDS(rds_diff)),
  comp_labels    = c("GTEx Adrenal Gland vs Liver"),
  endpoints      = c("DR", "DP", "DM"),
  fdr_thresholds = c(0.01, 0.05, 0.10, 0.20),
  panel_fdr      = 0.05, vline_power = 0.80, vline_fdr = 0.20,
  display_sizes  = DISPLAY_N,
  main_title     = "Differential Circadian Power Analysis (GTEx Adrenal Gland vs Liver)",
  out_pdf        = tmp_diff, width = 12.0, height = 4.7)
mirror_pdf(tmp_diff,
           c("figures/Fig2_differential_ADR_vs_LIV.pdf",
             "submission/figures/Fig2_differential_ADR_vs_LIV.pdf"))

cat("\nDONE\n")
