# SCP — Simulation-Based Power Calculation for Circadian Rhythmic Analysis

SCP is an R framework for sample size planning and study design in circadian
transcriptomic experiments. Given a small pilot RNA-seq dataset, it answers:
*How many samples do I need, and how should I space my time points?*

It supports both active (animal) and passive (human post-mortem, clinical) designs,
single-cohort rhythmicity detection, and two-group differential comparisons
(differential rhythmicity DR, differential phase DP, differential mesor DM).

---

## Quick start

```r
source("code/setup.R")   # loads all framework functions

# Estimate pilot parameters
bio <- estCircadianParam(expr_matrix, times_vector, period = 24)

# Specify study design and analysis options
des  <- CircadianDesignOptions(sample_sizes = c(24, 48, 72, 96), nsims = 100L)
aopt <- CircadianAnalysisOptions(alpha = 0.05)

# Run single-cohort power simulation
res  <- runSingleCohortPower(bio, des, aopt, methods = "DCP", plot = FALSE)

# Find N for 80% genome-wide power (with linear interpolation)
npower(res, target_power = 0.80, fdr = 0.05)

# Plot power summary
plotSingleCohortPower(res, out_pdf = "power_summary.pdf")
```

## Differential power (two groups)

```r
# Estimate two-group pilot parameters
bio_diff <- estCircadianParamTwoGroup(expr1, times1, expr2, times2)

# Run differential simulation
res_diff <- runDifferentialPower(bio_diff, des, aopt, methods = "DCP", plot = FALSE)

# Find N for 80% DR power
npower(res_diff, endpoint = "DR")

# Plot full differential power figure
plotDiffPower(list(comp1 = res_diff), comp_labels = "Control vs Disease")
```

## B vs m trade-off (DCP B-invariance)

DCP power under the cosinor truth is invariant to the number of distinct
time bins B — a study using B=6 bins with m=8 replicates achieves the same
power as B=24 with m=2 at the same total N. Demonstrated by `runSingleCohortGrid()`
and plotted by `plotBvsMPower()`.

## Waveform robustness (FMM)

Real circadian waveforms are often non-sinusoidal. SCP tests power under
waveform misspecification via the Frequency Modulated Möbius (FMM) model
(Rueda, Rodríguez-Collado & Peddada 2019, *Sci Rep*,
doi:10.1038/s41598-019-54569-1).

```r
# Simulate data under FMM waveform (omega=1: cosinor; omega->0: sharp peak)
cts <- seq(0, 22, by = 2)
sim <- simCircadianFMM(bio, cts, omega = 0.7)

# Plot FMM robustness figure
plotFMMViolation(fmm_df, nsims = 30, omega_fixed = 0.5)
```

## Key functions

| Function | Purpose |
|---|---|
| `estCircadianParam()` | Fit pilot parameters from real data (single group) |
| `estCircadianParamTwoGroup()` | Fit pilot parameters for two-group comparison |
| `runSingleCohortPower()` | Single-cohort rhythmicity power |
| `runDifferentialPower()` | Differential (DR/DP/DM) power |
| `runSingleCohortGrid()` | B vs m grid sweep |
| `runBootstrapDesignGrid()` | Bootstrap uncertainty quantification |
| `npower()` | Find N for target power (with interpolation) |
| `simCircadianFMM()` | Simulate single-cohort FMM waveform |
| `simCircadianDiffFMM()` | Simulate differential FMM waveform |
| `plotSingleCohortPower()` | 3-panel single-cohort figure |
| `plotDiffPower()` | Multi-panel differential power figure |
| `plotBvsMPower()` | B vs m trade-off figure |
| `plotFMMViolation()` | FMM robustness figure (single-cohort) |
| `plotFMMDifferential()` | FMM robustness figure (differential) |

## Documentation

- [`doc/TUTORIAL.md`](doc/TUTORIAL.md) — scenario-based usage guide
- [`vignettes/SCP_vignette.md`](vignettes/SCP_vignette.md) — full user guide with examples

## Reference

FMM model:
> Rueda C, Rodríguez-Collado A, Peddada SD (2019). *Sci Rep* 9, 17982.
> doi:10.1038/s41598-019-54569-1

Manuscript in preparation.
