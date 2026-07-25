#' Driver: regenerate every SCP manuscript figure and the supplementary table.
#'
#' Runs each current generator (one per figure; see FIGURES.md) as an isolated
#' Rscript process, so a setwd or par() in one script cannot affect another. It
#' does not stop on the first failure; a per-figure OK/FAILED summary is printed
#' at the end.
#'
#' Data requirements differ by script:
#'   * Replot scripts (fig1_fig2_replot, fig5_twoharm_framework) redraw from
#'     cached results under output/ and run against the current tree.
#'   * Sim-calibrated scripts (fig3_bootstrap_subject, fig4_twoharm_demo) need the
#'     controlled-access pilot matrices (GSE160521 striatal CSVs, GTEx CPM.all.norm
#'     / v10) and the fast C++ path, so they run on the calibration server only;
#'     they cache to output/ and skip the heavy stages on a re-run.
#'   * Fig 6 (GTEx recalibration, 2026-07) is now a PLOT-ONLY render
#'     (fig6_active_design_5panel_gtex.R) that reads two caches: A/B from
#'     two_harmonic/fig6AB_liver_active.R (GTEx Liver, K=1/K=2 active B-sweep), and
#'     C-E from archive/fig6_differential_Bsweep.R (GTEx Adrenal-vs-Liver DR/DP/DM).
#'     Those two upstream sim scripts use the bundled/cached GTEx Liver + Adr-vs-Liv
#'     pilots and write output/two_harmonic/results/ + output/diagnostics/.
#' Outputs are written under submission/figures/ (Supp Table 1 to submission/).

# Run from the package root (set SCP_ROOT to run from elsewhere).
ROOT <- Sys.getenv("SCP_ROOT", unset = ".")
setwd(ROOT)

GENERATORS <- c(
  "Fig 1A/1B + Fig 2 (single-cohort + differential, replot)" =
    "examples/publication/fig1_fig2_replot.R",
  "Fig 3 (bootstrap uncertainty; needs pilot data)" =
    "examples/publication/fig3_bootstrap_subject.R",
  "Fig 4 (two-harmonic demo; GTEx Liver)" =
    "examples/publication/two_harmonic/fig4_twoharm_demo.R",
  "Fig 5 (two-harmonic framework; replot)" =
    "examples/publication/two_harmonic/fig5_twoharm_framework.R",
  "Fig 6 (active design, 5-panel GTEx; plot-only render from caches)" =
    "examples/publication/fig6_active_design_5panel_gtex.R",
  "Supp Table 1 (pilot signal summary)" =
    "examples/publication/suppT1_pilot_dataset_summary.R"
)

results <- character(0)
for (nm in names(GENERATORS)) {
  script <- GENERATORS[[nm]]
  cat(sprintf("\n==== %s\n     %s\n", nm, script))
  code   <- system2("Rscript", script, stdout = "", stderr = "")
  status <- if (identical(code, 0L)) "OK" else sprintf("FAILED (exit %s)", code)
  results <- c(results, sprintf("  [%s] %s", status, nm))
}

cat("\n===== summary =====\n", paste(results, collapse = "\n"), "\n", sep = "")
