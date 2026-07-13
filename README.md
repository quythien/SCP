# SCP: Simulation-Based Circadian Power Analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-%E2%89%A54.4.0-blue.svg)](https://www.r-project.org/)

*How many samples, collected at what times of day, in which tissue, do you need
to find the genes that keep 24 hour time?*

## What SCP does

Circadian omics studies are usually short on samples, and the genes that follow
a daily rhythm are only a fraction of the transcriptome, so it is easy to plan a
study that is too small to see them. SCP answers the sample-size question by
simulation. You give it a small *pilot* dataset (a preliminary experiment in a
comparable tissue, or one of the pilots bundled with the package), and SCP
learns from that pilot how strong the rhythms are, how noisy the genes are, and
when during the day samples were collected. It then simulates studies of many
different sizes and reports how often each size would succeed. The result is a
recommended number of samples for detecting daily rhythms in one group, or for
comparing rhythms between two groups.

A few terms used throughout, in plain language. **Power** is the chance that a
study of a given size detects a gene that truly is rhythmic. **False discovery
rate (FDR)** is the share of the genes a study flags as rhythmic that are not
really rhythmic; you pick the level you are willing to tolerate (5 percent is
common). **Effect size** here is a gene's rhythm amplitude divided by its noise,
written `r_tilde = A / sigma`; larger means easier to detect. **Cosinor** is the
standard model that fits a cosine wave to a gene's expression over the day.

## Install

SCP is an R package. It depends on the Bioconductor package `limma`, so point
the installer at the Bioconductor repositories first.

```r
install.packages(c("remotes", "BiocManager"))
options(repos = BiocManager::repositories())   # resolves the Bioconductor dep 'limma'
options(timeout = 600)                          # bundled pilot data is about 100 MB
remotes::install_github("quythien/SCP", upgrade = "never", build_vignettes = TRUE)
library(SCP)
```

`build_vignettes = TRUE` builds the runnable tutorial so you can open it with
`browseVignettes("SCP")` after installing. If the GitHub download times out,
clone the repository and install from the local copy instead (`git` handles the
large download more reliably):

```r
# in a shell:  git clone https://github.com/quythien/SCP.git
options(repos = BiocManager::repositories())
remotes::install_local("SCP", upgrade = "never", build_vignettes = TRUE)
```

System requirements: R version 4.4.0 or newer and a C++ compiler (SCP uses Rcpp
and RcppArmadillo for the cosinor fit). On Linux: `apt install r-base-dev
libblas-dev liblapack-dev`. On macOS: `xcode-select --install`.

**macOS note (Fortran):** RcppArmadillo links against a Fortran and BLAS
toolchain, and CRAN R for macOS expects the official **gfortran** build. If the
install fails to compile with an error mentioning `gfortran`, `libgfortran`, or
`-lgfortran`, install the matching compiler from
<https://mac.r-project.org/tools/> (for R 4.4 on Apple Silicon or Intel this is
**gfortran 12.2-universal**, the file `gfortran-12.2-universal.pkg`), then retry
the install. This is a one time setup.

## Two ways to use SCP

**Point and click.** The bundled Shiny app runs the whole workflow in a browser
with no scripting. Install as above, then:

```r
SCP::launchShiny()
```

It runs entirely on your own machine. See the [Shiny app](#the-shiny-app)
section below for what the two halves of the app do.

**From R.** The rest of this page is the scripted API, which is what you want
for reproducible analyses.

## Quick start

Load a bundled pilot, set the study design, simulate, and read off the
recommended sample size.

```r
library(SCP)

# 1. Load a bundled pilot as a ready-to-use options object.
bio <- scp_load_pilot("human", "GTEx", "Adrenal", "All")

# Power is a proportion, so it barely depends on how many genes we simulate; we
# cap the gene count here only to keep this example fast. Drop this line for
# real runs.
bio$ngenes <- 1200L

# 2. Describe the planned study: a range of total sample sizes, an active
#    design that collects samples every 4 hours over 24 hours (6 time points).
design   <- CircadianDesignOptions(sample_sizes = c(24, 36, 48, 72, 96),
                                   nsims  = 30,
                                   design = "active",          # or "passive"
                                   cts    = seq(0, 20, by = 4))
analysis <- CircadianAnalysisOptions(alpha = 0.05)

# 3. Simulate and score with the single-harmonic cosinor test (K = 1).
res <- runSimsSingleCohort(bio.opts      = bio,
                           design.opts   = design,
                           analysis.opts = analysis,
                           K             = 1,
                           mc.cores      = 4)

# 4. Smallest N reaching 80 percent power at FDR 5 percent.
npower(res, target_power = 0.80, fdr = 0.05)$n
#> [1] 33

# 5. Plot the power curve.
plotSingleCohortPower(res, fdr = 0.05)
```

## Setting up a run: the three options

Every power run combines three small settings objects. Keeping them separate lets
you change one thing at a time: the biology, the study, or the analysis.

**1. `bio` (the pilot).** What the biology looks like. Load a bundled pilot or
estimate one from your own data, then tune a field if you want.

```r
bio <- scp_load_pilot("human", "GTEx", "Adrenal", "All")   # or estCircadianParam(expr, tod)
bio$ngenes <- 1200L                                        # fewer genes = faster
```

**2. `design` (the study).** Which sample sizes to try, how many simulations, and
the sampling design. This is where the active versus passive choice lives.

```r
# active: an animal study where you control the collection times
design <- CircadianDesignOptions(sample_sizes = c(24, 48, 72, 96), nsims = 100,
                                 design = "active", cts = seq(0, 20, by = 4))

# passive: a human post-mortem study where collection times cannot be controlled,
# so SCP draws them from the pilot's own time-of-day distribution (no cts needed)
design <- CircadianDesignOptions(sample_sizes = c(50, 100, 150, 200), nsims = 100,
                                 design = "passive")
```

**3. `analysis` (the inference).** The significance level, the multiple-testing
correction, and the FDR levels to report.

```r
analysis <- CircadianAnalysisOptions(alpha = 0.05, p.adjust.method = "BH",
                                     fdr_thresholds = c(0.05, 0.10, 0.20))
```

Then pass all three to a run function:

```r
res <- runSimsSingleCohort(bio, design, analysis, K = 1, mc.cores = 4)
```

For a bootstrap run you add a fourth object, `CircadianBootstrapOptions`; for a
two-group differential run, name the endpoints with `test_types = c("DR", "DP",
"DM")` in the design.

## Walkthrough with example outputs

Every block below runs against the installed package on data that ships with it,
so you can reproduce each figure. The runs are kept small (few simulations,
coarse sample-size grids) so they finish in seconds; use more simulations and a
finer grid for study-grade numbers.

### a. Pick a pilot

The pilots differ in how strong their rhythms are, and that is what decides how
many samples you will need. Loading a pilot and looking at its effect-size
distribution (`r_tilde = A / sigma`, amplitude over noise, per gene) tells you
what you are working with before you simulate anything.

```r
bio <- scp_load_pilot("human", "GTEx", "Adrenal", "All")
bio$ngenes <- 1200L

r_tilde <- bio$amplitude / bio$sigma_rhythmic
hist(r_tilde, breaks = 40, col = "#4C79A6", border = "white",
     main = "GTEx Adrenal pilot: per-gene effect size",
     xlab = expression("effect size  " * tilde(r) == A / sigma))
abline(v = median(r_tilde), col = "#C0392B", lwd = 2)
```

![Effect-size distribution of the GTEx Adrenal pilot](man/figures/pilot_effect_size.png)

### b. Single-cohort power

The core question for a one-group study: how many samples reach your target
power at your chosen FDR? `runSimsSingleCohort()` simulates the study across a
grid of sample sizes, `plotSingleCohortPower()` draws the curves, and `npower()`
reads off the recommended N.

```r
design   <- CircadianDesignOptions(sample_sizes = c(24, 36, 48, 72, 96),
                                   nsims = 30, design = "active",
                                   cts = seq(0, 20, by = 4))
analysis <- CircadianAnalysisOptions(alpha = 0.05,
                                     fdr_thresholds = c(0.05, 0.10, 0.20))

res <- runSimsSingleCohort(bio, design, analysis, K = 1, mc.cores = 4)

npower(res, target_power = 0.80, fdr = 0.05)$n
#> [1] 33

plotSingleCohortPower(res, title = "GTEx Adrenal, active q4h x 24h",
                      fdr = 0.05, reference_n = 48)
```

![Single-cohort power curve, effect-size strata, and discovery counts](man/figures/single_cohort_power.png)

Panel A is power against sample size at several FDR levels. Panel B breaks power
down by effect size, showing which genes are within reach at a given N. Panel C
is the number of true rhythmic genes recovered per effect-size band.

### c. Differential power (DR, DP, DM)

For a two-group study you often care not about rhythms in one group but about
how rhythms *differ* between groups. SCP scores three kinds of difference:
**DR** (differential rhythmicity, a gene oscillates in one group but not the
other), **DP** (differential phase, it oscillates in both but its peak time
shifts), and **DM** (differential mesor, where the **mesor** is the rhythm
adjusted average expression level, so DM is a shift in that baseline).
`runDifferentialPower()` returns power for each, and `plotDiffPower()` draws
them side by side.

```r
ex <- readRDS(system.file("extdata", "example_pilot_raw.rds", package = "SCP"))

# Estimate a two-group pilot from the two arms of the shipped demo matrix.
bio_diff <- estCircadianParamTwoGroup(
  data_1  = ex$adrenal$expr, data_2  = ex$liver$expr,
  times_1 = ex$adrenal$times, times_2 = ex$liver$times,
  paired_sigma = TRUE)

design_d <- CircadianDesignOptions(sample_sizes = c(40, 80, 120), nsims = 15,
                                   design = "passive", cts = bio_diff$cts,
                                   test_types = c("DR", "DP", "DM"))

res_diff <- runDifferentialPower(bio_diff, design_d,
                                 CircadianAnalysisOptions(alpha = 0.05),
                                 methods = "DCP", test_types = c("DR", "DP", "DM"),
                                 mc.cores = 4, plot = FALSE)

# Recommended N per endpoint at 80 percent power, FDR 20 percent.
sapply(c("DR", "DP", "DM"),
       function(ep) npower(res_diff, 0.80, 0.20, endpoint = ep)$n)
#>  DR  DP  DM
#>  74  NA  NA

plotDiffPower(list(res_diff), comp_labels = "Adrenal vs Liver (demo)",
              endpoints = c("DR", "DP", "DM"), panel_fdr = 0.05, vline_fdr = 0.20)
```

![Differential power for DR, DP, and DM between two demo groups](man/figures/differential_power.png)

This example is fully reproducible from shipped data. The two demo arms differ
mostly in rhythmicity, so DR reaches 80 percent power near N = 74 while DP and DM
have little to detect here and never cross 80 percent in this tiny run (hence the
`NA`). To plan a study around a hypothesized effect instead of a real second
group, set the fractions of differing genes directly on the pilot object
(`bio_diff$prop_DR`, `$prop_DP`, `$prop_DM`) before simulating; the user guide
vignette walks through that.

### d. Bootstrap uncertainty

When the pilot itself is small, the recommended N is uncertain. A **bootstrap**
resamples the pilot's subjects with replacement many times and re-runs the whole
estimate, so you can see how much the power curve wobbles.
`runBootstrapDesignGrid()` needs raw pilot expression (not a pre-summarized
pilot), so this example uses the shipped example expression matrix.

```r
expr <- as.matrix(read.csv(
  system.file("extdata/example/example_expression.csv", package = "SCP"),
  row.names = 1, check.names = FALSE))
tod <- scan(system.file("extdata/example/example_tod.csv", package = "SCP"),
            what = double())

bio_ex <- estCircadianParam(expr, tod, period = 24)

boot_opts <- CircadianBootstrapOptions(design_vector = tod, B_values = 6,
                                       N_values = c(24, 48, 96),
                                       nboot = 8, nsims_inner = 8,
                                       design = "passive", seed = 42)

boot <- runBootstrapDesignGrid(pilot_data = expr, pilot_times = tod,
                               boot.opts = boot_opts,
                               analysis.opts = CircadianAnalysisOptions(alpha = 0.05),
                               bio_diff.opts = bio_ex, mode = "single",
                               methods = "DCP", mc.cores = 4)

plotBootstrapDesignGrid(boot, panels = "A")
```

![Bootstrap power curve with a 95 percent confidence band](man/figures/bootstrap_ci.png)

The line is the bootstrap-mean power and the band is the 95 percent bootstrap
confidence interval. This example pilot has strong rhythms, so power is already
high by N = 24 and the band is narrow; a weaker or smaller pilot produces a wider
band, a signal that you should treat the recommended N cautiously. Use a larger
`nboot` (50 or more) and `nsims_inner` (20 or more) for real runs.

### e. Two-harmonic detection (K = 2)

Some genes peak in a shape that is not a clean cosine (asymmetric peaks, or two
peaks a day). The two-harmonic detector (`K = 2`) adds a 12 hour component and
picks those up better, at the cost of a little power on genes that are simple
cosines. You select it with a single argument.

```r
p_k2 <- detect_cosinor(expr, tod, K = 2)   # two-harmonic cosinor p-values

res_K2 <- runSimsSingleCohort(bio, design, analysis, K = 2, mc.cores = 4)

data.frame(N        = design$sample_sizes,
           power_K1 = round(rowMeans(res$marginal_power,    na.rm = TRUE), 3),
           power_K2 = round(rowMeans(res_K2$marginal_power, na.rm = TRUE), 3))
```

![Power of the single-harmonic (K = 1) versus two-harmonic (K = 2) detector](man/figures/two_harmonic_compare.png)

On this sinusoidal pilot K = 1 is slightly ahead, which is expected; on a pilot
with real 12 hour structure K = 2 pulls ahead. The two-harmonic model needs at
least 5 distinct sampling times per day to be identifiable (K = 1 needs 3).

### f. Use your own pilot

In practice you start from your own data: a gene by sample expression matrix and
a vector of collection times in hours. `prepCircadianData()` cleans and aligns
the inputs, and `estCircadianParam()` learns the pilot's rhythm parameters. Two
of those parameters are shown below: the per-gene effect size, and the
**acrophase** (the time of day at which each gene peaks). For a two-group pilot,
use `estCircadianParamTwoGroup()` instead.

```r
expr <- as.matrix(read.csv(
  system.file("extdata/example/example_expression.csv", package = "SCP"),
  row.names = 1, check.names = FALSE))
tod <- scan(system.file("extdata/example/example_tod.csv", package = "SCP"),
            what = double())

prep    <- prepCircadianData(expr, tod, input_type = "log2")
bio_own <- estCircadianParam(expr, tod, period = 24)

r_own <- bio_own$amplitude / bio_own$sigma_rhythmic
phase <- (bio_own$phase %% (2 * pi)) / (2 * pi) * 24   # radians to hours

par(mfrow = c(1, 2))
hist(r_own, breaks = 35, col = "#4C79A6", border = "white",
     main = "Estimated effect size", xlab = expression(tilde(r) == A / sigma))
hist(phase, breaks = seq(0, 24, by = 1.5), col = "#6FA46F", border = "white",
     main = "Estimated acrophase (peak time)", xlab = "time of day (h)")
```

![Effect-size and acrophase distributions estimated from your own pilot](man/figures/own_pilot_fit.png)

`bio_own` now plugs into `runSimsSingleCohort()` exactly like a bundled pilot.
The shipped example is simulated data (designed to resemble a GTEx Adrenal pilot
so it is freely shareable); swap in your own `expr` and `tod` and nothing else
changes.

## The Shiny app

`launchShiny()` opens a browser app with two linked halves.

**Circadian Power Study** (left) is the sample-size calculator. Cascading
dropdowns pick a pilot (species, then dataset, then tissue, then condition), and
you set a rhythmicity threshold, the sampling design (active timecourse or
passive time-of-death), a grid of sample sizes, and sliders for target power and
FDR. Clicking *Run simulation* draws the power curve and the recommended sample
size. On launch the app opens on the GTEx Adrenal (passive) pilot.

![Power curve for the GTEx Adrenal (passive) pilot, the app's default view](man/figures/app_power_curve_adrenal.png)

**Circadian Biomarker Detection** (right) explores the pilot itself: a core
clock-gene cosinor panel, a ranked table of rhythmic genes with their FDR
adjusted p-values and effect sizes, a per-gene cosinor fit, and pathway
enrichment (KEGG, Reactome, GO). All 161 bundled pilots are available, and you
can also upload your own pilot (an expression-matrix CSV and a one-column
time-of-day CSV) and get the same analyses on your data. An example pair of
upload files ships with the package:

```r
system.file("extdata/example/example_expression.csv", package = "SCP")
system.file("extdata/example/example_tod.csv",        package = "SCP")
```

## Pilot database

SCP ships 161 public circadian transcriptomic pilots across human, baboon,
mouse, and rat, spanning more than 100 tissue contexts. Three accessors browse
and load them.

| Function | Use |
|---|---|
| `scp_pilots()` | Return the manifest as a data frame (species, dataset, tissue, condition, design, sample size, effect size, status) |
| `scp_pilot_search(query, species, tissue_pattern)` | Filtered search |
| `scp_load_pilot(species, dataset, tissue, condition)` | Load one pilot as a `CircadianBioOptions` object |

```r
scp_pilots()                                             # full manifest
scp_pilot_search(species = "human", tissue_pattern = "adrenal")
bio <- scp_load_pilot("human", "GTEx", "Adrenal", "All")
```

A per-pilot effect-size summary (median effect size, its spread, rhythmic-gene
counts at several FDR thresholds, and the time-of-day sampling pattern) lives at
`summary/tissue_signal_summary.csv`.

## Documentation

- Runnable tutorial: `browseVignettes("SCP")`, then open "Getting started with
  SCP", or the fuller "User Guide".
- Function reference: `?SCP`, or `help(package = "SCP")`, or the bundled
  `SCP-manual.pdf`.
- News and release notes: `NEWS.md`.

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
