#' Replot Fig 2 with r_display_max = 3 (clip stratified r-axis at 3).
setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
old_wd <- setwd("code"); source("setup.R"); setwd(old_wd)

res <- readRDS("output/differential/results/diff_power_ADR_vs_LIV_20260423_085737.rds")

fig_path  <- "output/differential/figures/diff_power_fig_ADR_vs_LIV_rmax3.pdf"
fig_main  <- "output/main_figures/Fig2_differential_ADR_vs_LIV.pdf"
fig_sub   <- "submission/figures/Fig2_differential_ADR_vs_LIV.pdf"

plotDiffPower(
  res_list      = list(res),
  comp_labels   = "GTEx Adrenal Gland vs Liver",
  endpoints     = c("DR", "DP", "DM"),
  display_sizes = c(20, 40, 60, 80, 100, 120, 150, 200),
  vline_fdr     = 0.20,
  vline_power   = 0.80,
  r_display_max = 3,
  out_pdf       = fig_path,
  width         = 15,
  height        = 18
)

file.copy(fig_path, fig_main, overwrite = TRUE)
file.copy(fig_path, fig_sub,  overwrite = TRUE)
cat(sprintf("Saved: %s\n", fig_path))
cat(sprintf("Saved: %s\n", fig_main))
cat(sprintf("Saved: %s\n", fig_sub))
