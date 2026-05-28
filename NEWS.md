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
