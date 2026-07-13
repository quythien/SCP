# Figure provenance for the SCP manuscript

Which script produces each `submission/figures/*.pdf`. Script *filenames* carry
stale figure numbers from earlier revisions, so the map below is by **content**,
not filename. All generators live under `examples/publication/` (plus the
`two_harmonic/` subfolder).

| Paper figure (submission/figures) | Generator script | Notes |
|---|---|---|
| `Fig1A_single_cohort_AdrenalGland.pdf`, `Fig1B_single_cohort_Liver.pdf` | `fig1_fig2_replot.R` (draws from cache) | uses `library(SCP)` |
| `Fig2_differential_ADR_vs_LIV.pdf` | `fig1_fig2_replot.R` | draws from cache |
| `Fig3_bootstrap_singlecohort.pdf` | **`fig3_bootstrap_subject.R`** (current, Efron subject bootstrap; Caudate control n=59 + GTEx Muscle n=748) | supersedes the old `fig5_bootstrap_sc.R` (Putamen SCZ + Thyroid, per-gene) |
| `Fig4_twoharm_demo.pdf` | `two_harmonic/fig4_twoharm_demo.R` | GTEx v10 Liver; no `source()` calls, runs against the current tree |
| `Fig5_twoharm_framework.pdf` | `two_harmonic/fig5_twoharm_framework.R` (**BROKEN**, see caveats) | GTEx v10 Liver |
| `Fig6_active_design_5panel.pdf` | **`fig6_active_design_5panel.R`** (current, 5-panel) | A/B biomarker (K=1/K=2) on Putamen control n=59; C/D/E differential (DR/DP/DM) on Putamen-vs-Caudate control n=59/group. Supersedes the 2-panel `Fig6_active_BvsM.pdf` from `two_harmonic/fig6_v9_putamen_paired.R` (deleted in `a1ac962`; its A/B cache `output/two_harmonic/results/fig6_v9_data.rds` is retained and reused here). |
| `SuppFig_clock_gene_profiles.pdf`, `SuppFig_FMM_diagnostic.pdf` | **no generator on disk** (not even under `archive/`) | undocumented gap; regenerate or drop from the supplement |

## Fig 3 (bootstrap) reproduction (current)

`examples/publication/fig3_bootstrap_subject.R` regenerates Fig 3 end to end:
Efron **subject** bootstrap (resample pilot subjects with replacement, carrying
their collection times, refit per draw) on **Caudate control** (GSE160521, n=59)
and **GTEx skeletal muscle** (n=748); reports the projected power
curve + the 2.5-97.5 percentile 95% bootstrap CI; `B_out=50`, `N_sim_inner=25`. Requires the raw pilot matrices
(controlled-access GSE160521 CSVs + GTEx `CPM.all.norm.RData`), so it runs on
this server only. The corresponding package entry point is
`runBootstrapDesignGrid(..., resample = "subject")` (the default) in
`R/bootstrap_sim.R`.

## Reproducibility caveats

- The *sim* generators for Figs 1-2 and 4-6 historically opened with
  `setwd("code"); source("setup.R")`, but `code/setup.R` was removed in the
  v0.4.0 release-prep commit. Only the *replot* scripts (`fig1_fig2_replot.R`,
  `fig3_bootstrap_replot.R`), `two_harmonic/fig4_twoharm_demo.R`, and the new
  `fig3_bootstrap_subject.R` run against the current tree.
- **`fig5_twoharm_framework.pdf`'s generator is broken independently of the
  `setup.R` issue.** `two_harmonic/fig5_twoharm_framework.R` directly sources
  `code/utils.R`, `code/simulation.R`, `code/bootstrap_sim.R`, and
  `code/detection.R`; all four were removed along with `code/setup.R`, and only
  `code/options.R`/`code/estimation.R`/`code/pilot_database.R`/
  `code/plot_diff.R`/`code/plot_single_cohort.R`/`code/runner.R` survive in
  `code/`. Swap these `source()` calls for `library(SCP)` (and confirm the
  script only calls exported functions) before Fig 5 can be regenerated.
- **`fig3_bootstrap_replot.R` is a landmine, not a harmless duplicate.** It
  writes to the *same* `Fig3_bootstrap_singlecohort.pdf` path as
  `fig3_bootstrap_subject.R`, but plots the old, superseded Putamen-SCZ /
  Thyroid per-gene panels from a stale cache. Running it after the real
  generator silently clobbers the current Fig 3 with obsolete data. Move it to
  `archive/` (or delete) once Fig 3's provenance above is trusted.
- **Orphan generators** under `examples/publication/` that are not referenced
  anywhere in this file and are also broken via the `setup.R`/`code/*` pattern:
  `fig1_single_cohort_power.R`, `fig2_paired_sims.R`,
  `fig3_bvsm_method_comparison.R` (writes an unrelated `fig3_method_comparison.pdf`),
  `fig4_sensitivity.R` (writes `fig4_sensitivity.pdf`, not the current Fig 4),
  `fig6_phase_mse.R`. `suppT1_pilot_dataset_summary.R` (writes
  `submission/tissue_signal_summary.csv`, Supp Table 1) is also unreferenced
  but not broken. `fig5_bootstrap_sc.R` and `fig6_cosinor_rebuild.R` are
  intentionally-deprecated/non-adopted alternates, not orphans.
- Fast C++ path (about 195x on cosinor fits): as of the current session,
  `NAMESPACE` now has `useDynLib(SCP, .registration = TRUE)` and
  `importFrom(Rcpp, sourceCpp)`, so the compiled `.Call()` routines resolve
  correctly for the first time. The path is still **not** engaged
  automatically, though: every call site (`R/simulation.R`, `R/bootstrap_sim.R`,
  `R/runner.R`) gates on a `.CPP_LOADED` flag that is checked but never set
  anywhere in `R/`. Opt in manually with
  `library(Rcpp); .CPP_LOADED <- TRUE` after `library(SCP)`.
- Cached intermediate results and diagnostics live under gitignored `output/`.

## Minimal consolidation recommendation

There is no single script that regenerates every figure; a maintainer has to
know which of ~15 scripts under `examples/publication/` (several dead) is
current. Recommended: one `examples/publication/run_all_figures.R` (or
Makefile) that (1) fails fast if anything it calls still references
`code/setup.R` or a deleted `code/*.R` file, (2) runs exactly the six
generators in the table above plus `suppT1_pilot_dataset_summary.R` and
whatever regenerates the two missing SuppFigs, (3) writes only into
`submission/figures/`, and (4) has `fig3_bootstrap_replot.R`,
`fig5_bootstrap_sc.R`, `fig6_cosinor_rebuild.R`, and the orphan scripts above
moved into `examples/publication/archive/` so they stop shadowing the current
generator for the same output file.

## Fig 6 (active design, 5-panel) reproduction

`fig6_active_design_5panel.R` regenerates the adopted Fig 6 end to end in three
stages:
1. Two-group pilot calibration: Putamen control (group 1) vs Caudate control
   (group 2) from GSE160521 via `estCircadianParamTwoGroup(paired_sigma=TRUE)`,
   seed 2025 (cache `output/diagnostics/pilot_put_vs_cau_control.rds`).
2. Differential active B-sweep for DR, DP, DM: B in {4,6,8,12,24}, N in
   {24,48,96,144,192,240}, NSIMS=60, NGENES=3000, BH-FDR 0.05, seed 2025+B
   (cache `output/diagnostics/Fig6_diff_Bsweep_PutCau.rds`).
3. Render: panels A/B read from the biomarker cache
   `output/two_harmonic/results/fig6_v9_data.rds` (Putamen control, K=1/K=2);
   panels C/D/E from the sweep cache. Output
   `submission/figures/Fig6_active_design_5panel.pdf`.

Stages 1-2 need the controlled-access GSE160521 striatal CPM matrices (server
only) and the fast C++ path; they cache to RDS, so a re-run skips straight to
the render. Verified against the cached output: A/B power is B-invariant to
within about 3 (K=1) and 2 (K=2) percentage points and reaches 80% near N=50;
DR and DM reach 80% near N=165 and N=140, DP stays below 80% over the grid, and
the DR/DP/DM curves collapse across B to within about 2 percentage points.

The earlier `fig6_differential_Bsweep.R` was the exploratory precursor on the
GTEx Adrenal-vs-Liver pair (output `output/diagnostics/Fig6_differential_Bsweep.pdf`);
it is kept for reference but is not the adopted figure.
