# SCP: Simulation-Based Circadian Power Analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-%E2%89%A54.4.0-blue.svg)](https://www.r-project.org/)

Circadian omics studies are frequently constrained by limited sample
availability, uncontrolled collection times, and the multiple-testing burden
imposed by genome-wide analyses. Existing power analysis tools do not support
pilot-calibrated planning for circadian biomarker detection or differential
rhythmicity analysis.

**SCP** is a simulation-based circadian power framework for sample-size
planning in transcriptomic studies. SCP leverages pilot-derived distributions
of gene-level amplitude, noise, phase, and time-of-day sampling to capture
transcriptome-wide signal heterogeneity, and estimates FDR-controlled power
for single-cohort rhythmic-biomarker detection as well as two-group
differential analyses of rhythmicity (DR), phase (DP), and MESOR (DM). A
bootstrap layer propagates pilot-estimation uncertainty into confidence
intervals on every power estimate. The framework is extensible to any
model-based circadian approach, cosinor or non-cosinor.

To make planning practical, SCP ships a curated **pilot database** of public
circadian transcriptomic datasets across human, baboon, mouse, and rat,
spanning 100+ tissue contexts. Users pick a pilot matching their planned
study and get a recommended sample size in seconds.

## Try it: launch the Shiny app locally

The fastest way to use SCP is the bundled Shiny GUI. Three lines:

```r
install.packages("remotes")
remotes::install_github("quythien/SCP")
SCP::launchShiny()
```

A browser window opens with cascading dropdowns (species → dataset → tissue
→ condition), sliders for target power and FDR, a sampling-design panel,
and a one-click "Run simulation" button. Power curve and recommended sample
size render in 2-5 seconds per click. All 127 bundled pilots are
immediately available.

The Shiny app runs entirely on your machine, no account or hosting needed.

## Install (R package only)

If you want the R API without the GUI:

```r
remotes::install_github("quythien/SCP")
library(SCP)
```

Or from a local clone:

```r
remotes::install_local("path/to/SCP")
```

System requirements: R >= 4.4.0, a C++ compiler. Rcpp + RcppArmadillo are
used for the cosinor fit hot path. On Linux: `apt install r-base-dev
libblas-dev liblapack-dev`. On macOS: `xcode-select --install`. Suggested
packages for the Shiny GUI: `shiny`, `rmarkdown`.

## Quick start (R API)

```r
library(SCP)

# 1. Discover what's bundled
scp_pilots()                                  # full manifest tibble
scp_pilot_search(species = "human",
                 tissue_pattern = "adrenal")  # filter

# 2. Load a pilot
bio <- scp_load_pilot("human", "GTEx", "Adrenal", "All")

# 3. Set the design and run a small power sim
design   <- CircadianDesignOptions(sample_sizes = c(20, 40, 60, 80),
                                    nsims  = 20,
                                    design = "active",        # or "passive"
                                    cts    = seq(0, 22, by = 4))
analysis <- CircadianAnalysisOptions(alpha = 0.05)
res      <- runSimsSingleCohort(bio.opts      = bio,
                                design.opts   = design,
                                analysis.opts = analysis,
                                K             = 1,    # K = 2 for two-harmonic
                                mc.cores      = 4)

# 4. Recommended N at 80% power, FDR 5%
np <- npower(res, target_power = 0.80, fdr = 0.05)
np$n
#> [1] 36

# 5. Plot
plotSingleCohortPower(res, fdr = 0.05)
```

## Pilot database

The bundled pilot database lives in `inst/extdata/pilots/` and is indexed by
`manifest.csv` (species, dataset, tissue, condition, design, sample size,
effect size r-tilde, status). Three accessors:

| Function | Use |
|---|---|
| `scp_pilots()` | Return the manifest as a data.frame |
| `scp_load_pilot(species, dataset, tissue, condition)` | Load one pilot's `CircadianBioOptions` |
| `scp_pilot_search(query, species, tissue_pattern)` | Filtered search |

A per-pilot effect-size summary (median r-tilde, IQR, rhythmic-gene counts at
FDR thresholds, time-of-day sampling pattern) lives at
`summary/tissue_signal_summary.csv`.

## Detection API

A single entry point for the cosinor F-test:

```r
detect_cosinor(expr, times, K = 1)   # single-harmonic
detect_cosinor(expr, times, K = 2)   # two-harmonic (ultradian extension)
```

Identifiability requires at least `2K + 1` distinct sampling phases per
period (`B >= 3` for `K = 1`, `B >= 5` for `K = 2`).

## Documentation

- Runnable vignette: `vignettes/SCP_tutorial.Rmd`
- Function reference: `?SCP` or `help(package = "SCP")`
- News and release notes: `NEWS.md`

## Citation

If you use SCP in your research, please cite:

```
Pham, T. Q. (2026). SCP: Simulation-Based Circadian Power Analysis.
  R package version 0.4.3. https://github.com/quythien/SCP
```

A `CITATION.cff` is provided.

## License

MIT. See `LICENSE`.

## Contact

Issues and feature requests: <https://github.com/quythien/SCP/issues>
