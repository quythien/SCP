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

The fastest way to use SCP is the bundled Shiny GUI:

```r
install.packages(c("remotes", "BiocManager"))
# SCP depends on the Bioconductor package 'limma'; point remotes at the
# Bioconductor repositories so it is resolved automatically:
options(repos = BiocManager::repositories())
# SCP bundles ~100 MB of pilot data, so raise the download timeout from the
# 60s default to avoid "download ... failed" on slower connections:
options(timeout = 600)
remotes::install_github("quythien/SCP", upgrade = "never")
SCP::launchShiny()
```

If the GitHub download still times out, clone and install locally instead
(`git` handles the large download more reliably):

```r
# in a shell:  git clone https://github.com/quythien/SCP.git
options(repos = BiocManager::repositories())
remotes::install_local("SCP", upgrade = "never")
```

The app has two linked halves. **Circadian Power Study** (left) drives the
sample-size calculation: cascading dropdowns (species → dataset → tissue →
condition), a pilot rhythmicity threshold, the sampling design (active
timecourse or passive time-of-death), a sample-size grid, and sliders for
target power and FDR. One click on "Run simulation" draws the FDR-controlled
power curve and reads off the recommended sample size in a few seconds.
**Circadian Biomarker Detection** (right) explores the pilot itself: a core
clock-gene cosinor panel, a ranked table of rhythmic genes with BH q-values and
effect sizes, a per-gene cosinor view, and pathway enrichment (KEGG / Reactome /
GO). On launch the app opens on the GTEx Adrenal (passive) pilot so there is
something to look at right away; all 161 bundled pilots are available from the
dropdowns.

![Power curve for the GTEx Adrenal (passive) pilot, the app's default view](man/figures/app_power_curve_adrenal.png)

The app runs entirely on your machine, no account or hosting needed.

You can also calibrate from **your own pilot** via "Upload my own pilot": an
expression-matrix CSV (genes in rows, samples in columns, up to 5 GB) and a
one-column time-of-day CSV. The upload is fit on the fly for both the
single-harmonic (K = 1) and two-harmonic (K = 2) detectors, gene IDs are mapped
to symbols where an annotation package is available, and the gene table and
enrichment work on your data exactly as on a bundled pilot. A ready-made example
pilot (**simulated to resemble a GTEx Adrenal pilot** so it is freely shareable)
ships with the package:

```r
system.file("extdata/example/example_expression.csv", package = "SCP")
system.file("extdata/example/example_tod.csv",        package = "SCP")
```

## Install (R package only)

If you want the R API without the GUI:

```r
install.packages("BiocManager")
options(repos = BiocManager::repositories())   # resolves the Bioconductor dep 'limma'
options(timeout = 600)                          # bundled pilot data is ~100 MB
remotes::install_github("quythien/SCP", upgrade = "never")
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

**macOS note (Fortran):** RcppArmadillo links against a Fortran/BLAS toolchain,
and CRAN R for macOS expects the official **gfortran** build. If the install
fails to compile with an error mentioning `gfortran`, `libgfortran`, or
`-lgfortran`, install the matching compiler from
<https://mac.r-project.org/tools/> (for R 4.4 on Apple Silicon/Intel this is
**gfortran 12.2-universal** -- `gfortran-12.2-universal.pkg`), then retry the
install. This is a one-time setup.

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

## Bootstrap uncertainty

Quantifying uncertainty in a projected power curve is a two-step call:

```r
# 1. Compute: resample pilot SUBJECTS with replacement (Efron subject bootstrap;
#    resample = "subject", the default), refit the pilot summary on each draw,
#    and re-estimate power across the (N, B) grid.
boot <- runBootstrapDesignGrid(pilot_data, pilot_times,
                               boot.opts, analysis.opts, bio_diff.opts)

# 2. Plot: the bootstrap-mean curve plus the pointwise 95% bootstrap CI.
plotBootstrapDesignGrid(boot)
```

`runBootstrapDesignGrid()` returns the arrays (`power_mean`, and `power_ci_lo` /
`power_ci_hi`, the 2.5th / 97.5th percentiles of the outer resamples);
`plotBootstrapDesignGrid()` draws the mean line and that 95% bootstrap CI band.
The Shiny app wires the two together, so the figure appears automatically there.

## Documentation

- Runnable vignette: `vignettes/SCP_tutorial.Rmd`
- Function reference: `?SCP` or `help(package = "SCP")`
- News and release notes: `NEWS.md`

## Citation

If you use SCP in your research, please cite:

```
Pham, T. Q. (2026). SCP: Simulation-Based Circadian Power Analysis.
  R package, latest version at https://github.com/quythien/SCP
```

A `CITATION.cff` is provided.

## License

MIT. See `LICENSE`.

## Contact

Issues and feature requests: <https://github.com/quythien/SCP/issues>
