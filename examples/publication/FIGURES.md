# Figure provenance — SCP manuscript

Which script produces each `submission/figures/*.pdf`. Script *filenames* carry
stale figure numbers from earlier revisions, so the map below is by **content**,
not filename. All generators live under `examples/publication/` (plus the
`two_harmonic/` subfolder).

| Paper figure (submission/figures) | Generator script | Notes |
|---|---|---|
| `Fig1A_single_cohort_AdrenalGland.pdf`, `Fig1B_single_cohort_Liver.pdf` | `fig1_fig2_replot.R` (draws from cache) | uses `library(SCP)` |
| `Fig2_differential_ADR_vs_LIV.pdf` | `fig1_fig2_replot.R` | draws from cache |
| `Fig3_bootstrap_singlecohort.pdf` | **`fig3_bootstrap_subject.R`** (current, Efron subject bootstrap; Caudate control n=59 + GTEx Muscle n=748) | supersedes the old `fig5_bootstrap_sc.R` (Putamen SCZ + Thyroid, per-gene) |
| `Fig4_twoharm_demo.pdf` | `two_harmonic/fig4_twoharm_demo.R` | GTEx v10 Liver |
| `Fig5_twoharm_framework.pdf` | `two_harmonic/fig5_twoharm_framework.R` | GTEx v10 Liver |
| `Fig6_active_BvsM.pdf` | `two_harmonic/fig6_v9_putamen_paired.R` (DELETED in cleanup commit `a1ac962`; recover via `git show a1ac962^:...`; writes directly to submission/figures, no rename) | Putamen control n=59, active B-sweep, K=1 + K=2 panels. Shipped PDF is a light hand-edit (added r-tilde, panel relabels). `fig6_cosinor_rebuild.R` is a different, non-adopted rebuild. |

## Fig 3 (bootstrap) reproduction — current

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

- The *sim* generators for Figs 1–2 and 4–6 historically opened with
  `setwd("code"); source("setup.R")`, but `code/setup.R` was removed in the
  v0.4.0 release-prep commit. Only the *replot* scripts (`fig1_fig2_replot.R`,
  `fig3_bootstrap_replot.R`) and the new `fig3_bootstrap_subject.R` run against
  the current tree; the others need `source("code/setup.R")` swapped to
  `library(SCP)` before they will run.
- Fast C++ path (≈195× on cosinor fits) is off by default because the package
  `NAMESPACE` lacks `useDynLib`; scripts activate it with
  `library(Rcpp); dyn.load(system.file("libs","SCP.so",package="SCP")); .CPP_LOADED <- TRUE`.
- Cached intermediate results and diagnostics live under gitignored `output/`.

## Fig 6 differential extension (advisor addition, 2026-07-07)

`fig6_differential_Bsweep.R` tests whether B-invariance also holds for the
differential endpoints (DR, differential rhythmicity; DP, differential phase;
DM, differential mesor) under a balanced active design. Pilot: GTEx
Adrenal-vs-Liver paired (same as Fig 2); B in {4,6,8,12,24}, N in
{24,48,96,144,192,240}, NSIMS=80, BH-FDR 0.05. Result: all three endpoints are
B-invariant (spread <=0.3 percentage points across B). Output:
`output/diagnostics/Fig6_differential_Bsweep.pdf`.
