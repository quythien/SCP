# SCP 0.4.60

## Plotting

* `runBootstrapDesignGrid()` now also computes the plug-in (point-estimate) power
  from the full pilot without resampling and returns it as `power_plugin`.
  `plotBootstrapDesignGrid()` overlays it as the plug-in line, so the single-panel
  default plot reproduces the manuscript bootstrap figure: the plug-in estimate
  against the bootstrap mean and its 95 percent percentile interval. The per-cell
  simulation was refactored into one shared routine used by both the plug-in and
  the bootstrap draws. The README bootstrap section is now drawn by this code on
  the bundled caudate control pilot.

# SCP 0.4.59

## README and bundled data

* The README walkthrough now runs on control-region pilots (nucleus accumbens,
  caudate, putamen) from the public GSE160521 human striatal study, bundled with
  the package. The single-cohort and differential examples run end to end on that
  data, and the two-harmonic section shows the manuscript GTEx Liver figures.
  Each section carries citations, with a References list at the end.

## Plotting

* `plotSingleCohortPower()` now accepts an optional pilot (`bio.opts`,
  `design.opts`, `analysis.opts`); given those, it draws the effect-size
  stratified panels directly from the pilot so every band is populated. Passing a
  simulation result as before behaves exactly as it did.

## Documentation

* `lOD`/`lOD2` are documented as the log residual noise (`sigma = exp(lOD)`)
  rather than "over-dispersion", with a matching glossary entry in the README.
* Maintainer email updated to qtp1@pitt.edu.

# SCP 0.4.58

## Public API

* **Cosinor-only public surface.** The frequency-modulated Mobius (FMM)
  parameters have been removed from the exported interface: `CircadianDesignOptions()`
  no longer takes `omega`/`beta`, `CircadianBioOptions()` no longer takes
  `omega_rhythmic`/`alpha_rhythmic`/`omega_dist`/`alpha_dist`/`paired_omega`/`paired_alpha`,
  and `runSimsSingleCohort()` no longer accepts `method = "FMM"` or the
  `harmonics` argument. The public detectors are the single-harmonic (K = 1)
  and two-harmonic (K = 2) cosinor F-tests. The internal K-harmonic engine is
  unchanged, so `detect_cosinor(K = 2)` and the two-harmonic simulator are
  unaffected.
* **Runnable examples on every exported function.** Most examples run directly
  on the bundled demo pilot; the heavier two-group and bootstrap examples are
  wrapped in `\donttest{}`, and only `launchShiny()` remains `\dontrun{}`.
* **All contributors credited** in `DESCRIPTION` and `CITATION.cff`.

## Documentation

* **Reference manual rewritten and slimmed** (35 pages): a grouped overview of
  the main functions, real per-function descriptions in place of title-repeats,
  plain-language argument text, and no dedicated pages for the trivial S3
  `print`/`plot` methods. User-facing text now says "cosinor" rather than the
  legacy internal name "DCP".
* Dropped the unused `reshape2` dependency.

# SCP 0.4.57

## Release hardening

* **Package now passes `R CMD check` with no errors.** Added the missing
  `useDynLib` and namespace imports so the compiled cosinor path is reachable
  via `.Call()`; declared the vignette builder and its suggested packages;
  moved `stats`/`graphics`/`utils` to Imports; excluded stray top-level files
  from the build. Non-ASCII characters in console output and comments were
  converted to ASCII.
* **Subject bootstrap is the default** for `runBootstrapDesignGrid()`
  (`resample = "subject"`); the reported band is the 2.5-97.5 percentile
  bootstrap interval. See the README "Bootstrap uncertainty" section.
* **More parameters exposed.** `makeAdaptiveRStrata(r_pctile_cap=)` is now an
  argument rather than a hardcoded constant, and `estCircadianParam*()` accept a
  user-settable `top_k`.

## Bug fixes

* `SeqModelSel()` no longer returns a function object on an unrecognized
  `method` string (a shadowed `stop` variable); it now errors clearly.
* The legacy flat-argument path of `runSimsDiff()` no longer errors on missing
  waveform/effect-size defaults.
* Differential simulation now errors on an unimplemented `DCmethod` instead of
  silently reporting near-zero power.
* Removed documented-but-unimplemented options (`DCP_DiffR2()`
  "permutation"/"bootstrap"; differential "LimoRhyde"/"DODR") so the documented
  interface matches the implementation.
* Residual-noise estimation guards against zero degrees of freedom at exactly
  three sampled times.
* Bootstrap and ground-truth grids now warn on a failed cell instead of
  silently leaving it `NA`.
* Loop bounds use `seq_len()`/`seq_along()` and result-file loading uses an
  isolated environment, removing 0-length and scope-leak edge cases.

# SCP 0.4.56

## Shiny app

* **Opens ready to use.** The app now launches on the GTEx Adrenal (passive)
  pilot, so a power curve and the gene explorer are visible before any clicks.
* **Biomarker explorer upgrades.** The rhythmic-gene table is sortable with
  server-side paging and an adjustable row count; the gene-lookup list can be
  sorted and reversed; enrichment controls sit on a single row. Reactome's
  annotation database installs on demand rather than at launch.
* **Uploads are fuller-featured.** The upload cap is raised to 5 GB, and an
  uploaded pilot is fit for both the single-harmonic (K = 1) and two-harmonic
  (K = 2) detectors, so switching the detector changes both the gene table and
  the power curve. Enrichment uses the uploaded matrix's genes as background.

## Bug fixes

* The core clock-gene panel no longer errors ("figure margins too large") when
  no clock gene passes the chosen threshold; it shows a short note instead.

## Internal changes

* Verbose console output and code comments were tidied (plain one-line
  diagnostics in place of decorative banners); no change to results.

# SCP 0.4.48

## Major changes

* **Pilot database enriched and expanded.** The bundled database now holds 161
  pilots across human, mouse, rat, and baboon. Every pilot's per-gene
  `rhythm_fit` table now carries `gene` (native ID), `symbol` (mapped gene
  symbol, original ID retained when none maps), and `mesor`, in addition to
  `pvalue`, `A`, `phi`, and `sigma`. The pilot rds are xz-compressed.

* **Shiny app is now a full study-design and biomarker tool.** Two panels:
  "Circadian Power Study" (pilot selection, threshold, sampling design, sample
  size grid, power curve, recommended N) and "Circadian Biomarker Detection"
  (a core clock-gene cosinor panel, a ranked rhythmic-gene table with BH q-values
  and effect sizes, a per-gene cosinor view, and Metascape-style pathway
  enrichment via Enrichr for KEGG / Reactome / GO). Both work on bundled pilots
  and on user-uploaded data.

* **Upload path.** Uploaded data is fit on the fly; the app detects the species
  from the gene IDs and maps them to symbols for nine model organisms (human,
  mouse, rat, zebrafish, Drosophila, C. elegans, S. cerevisiae, P. falciparum,
  Arabidopsis) when the matching annotation package is installed, keeping the
  original ID otherwise. Uploaded cosinor plots overlay the actual sample points.

* **Runtime-selectable rhythmicity threshold.** `alpha_pilot` (raw p or
  BH-FDR) is chosen at load time and reslices each pilot's `rhythm_fit`, so the
  threshold is no longer frozen at build time.

* **`launchShiny()`** checks for and installs the app's optional helper packages
  (enrichR, writexl, org.Hs/Mm.eg.db) on first launch.

# SCP 0.4.0

## Major changes

* **Unified detector API.** A single entry point, `detect_cosinor(expr, times, K)`,
  replaces the prior `detect_DCP` (single-harmonic) and `detect_FMM` (K-harmonic)
  surface. `K = 1` is the single-harmonic cosinor F-test; `K = 2` is the
  two-harmonic extension used throughout the SCP manuscript. The legacy names
  remain as internal aliases for backward compatibility.

* **Two-harmonic framework.** New generator `simCircadianSingleCohort2H()` and
  pilot estimator `estCircadianParam2H()` support the two-harmonic case
  end-to-end. The frequency-modulated Mobius (FMM) family is no longer part of
  the recommended workflow; legacy FMM generators are retained but undocumented.

* **Pilot dataset database.** Bundled pilots are now organised under
  `inst/extdata/pilots/<species>/<dataset>_<tissue>_<condition>.rds` with a
  central `inst/extdata/pilots/manifest.csv`. Three new accessor functions
  expose the database:
  - `scp_pilots()` lists the manifest.
  - `scp_load_pilot()` returns a single pilot as a `CircadianBioOptions`.
  - `scp_pilot_search()` filters by species, tissue, sample type, or text.

* **Tissue summary table.** `output/supp_tissue_summary/tissue_signal_summary.csv`
  (also at `submission/tissue_signal_summary.csv`) records per-pilot effect
  size and rhythmicity counts. Four new columns describe the sampling pattern:
  `tod_n_unique`, `tod_cycles`, `tod_step_hr`, `tod_phases`.

## Documentation

* `README.md`, `doc/README.md`, `doc/TUTORIAL.md`, and the vignette have been
  realigned to the cosinor-only framework and the `detect_cosinor(K)` API.
* A new runnable Rmd tutorial is available at `vignettes/SCP_tutorial.Rmd`.
* `CITATION.cff` is now provided.

## Internal changes

* Runtime output that was previously emitted via `cat()` is now emitted via
  `message()` for progress and diagnostics, so it can be suppressed via
  `suppressMessages()`. `cat()` is reserved for `print` and `summary` methods
  whose explicit job is to display formatted output.
* All R/ files mirror `code/` exactly. Roxygen blocks have been tightened
  across the public-facing functions.
* `R/runner.R` `print.SCPDiffResult` no longer crashes when called: the
  `intersect()` call that previously passed an unsupported argument has been
  fixed.

## Bug fixes

* `prepCircadianData` example is now wrapped in `\dontrun{}` so it no longer
  triggers an R CMD check ERROR.
* Mouse GSE54651 pilot fits previously returned all-NA p-values because
  `count_clean` was a data.frame and per-gene rows broke `mean()`. The
  generator now coerces to a numeric matrix before fitting.
