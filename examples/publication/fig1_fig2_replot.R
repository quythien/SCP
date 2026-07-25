# Replot Fig 1 (Adrenal, Liver) and Fig 2 (Adrenal vs Liver) from cache.
# Run from the package root: Rscript examples/publication/fig1_fig2_replot.R
# (set SCP_ROOT to run from elsewhere).
setwd(Sys.getenv("SCP_ROOT", unset = "."))
suppressMessages(library(SCP))
# Use the installed package's plotSingleCohortPower()/plotDiffPower() directly
# (rather than source()-ing R/plot_single_cohort.R / R/plot_diff.R) so their
# internal helpers (e.g. add_se_bars(), .prepDiffStratified()) resolve from
# within the SCP namespace instead of requiring every helper file to be
# sourced by hand.

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

draw_sc <- function(primary_rx, panelA_rx, title, dests,
                    panel_labels = c("A", "B", "C")) {
  rds_p <- pick_latest("output/single_cohort/results", primary_rx)
  rds_a <- pick_latest("output/single_cohort/results", panelA_rx)
  cat(sprintf("  B/C: %s\n  A  : %s\n", basename(rds_p), basename(rds_a)))
  tmp <- tempfile(fileext = ".pdf")
  do.call(plotSingleCohortPower,
          c(list(res = readRDS(rds_p), panel_a_res = readRDS(rds_a),
                 out_pdf = tmp, title = title, panel_labels = panel_labels),
            SC_ARGS))
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
          "submission/figures/Fig1B_single_cohort_Liver.pdf"),
        panel_labels = c("D", "E", "F"))

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
