# SCP: Semiparametric Circadian Power — User Guide

**Author:** Thien Pham  
**Package:** SCP (Semiparametric Circadian Power)  
**Repository:** [DiffCircaPower](https://github.com/quythien/DiffCircaPower)

---

## Table of Contents

1. [Overview](#1-overview)
2. [Quick Start](#2-quick-start)
3. [Single-Cohort Power](#3-single-cohort-power)
4. [Differential Power](#4-differential-power)
5. [B vs m Trade-off and Method Recommendation](#5-b-vs-m-trade-off-and-method-recommendation)
6. [Bootstrap Uncertainty](#6-bootstrap-uncertainty)
7. [Working with Result Objects](#7-working-with-result-objects)
8. [Complete Publication Example](#8-complete-publication-example)
9. [Function Reference](#9-function-reference)

---

## 1. Overview

**SCP** is an R framework for transcriptome-wide power analysis of circadian rhythm studies. It addresses four questions:

| Question | Function |
|----------|----------|
| How many samples to detect rhythmic genes in one group? | `runSingleCohortPower()` |
| How many samples per group to detect DR / DP / DM? | `runDifferentialPower()` |
| Given fixed N, how to spread across time points vs replicates? | `recommendDesign()` |
| How wide are power estimates when the pilot is small? | `runBootstrapDesignGrid()` |

The framework is **semiparametric**: amplitude, noise, and proportion-rhythmic distributions are estimated from real pilot data; the cosinor waveform structure governs simulation.

### Detection endpoints

| Endpoint | Definition | Biological example |
|----------|------------|--------------------|
| **DR** — differential rhythmicity | Gene oscillates in one condition but not the other | Clock gene disrupted in neurodegeneration |
| **DP** — differential phase | Gene oscillates in both but peak timing shifts | Circadian misalignment, tissue-specific timing |
| **DM** — differential mesor | Gene oscillates in both but baseline shifts | Inflammatory state superimposed on intact oscillation; circadian analogue of DE |

### Detection methods

| Method | Single-cohort | Differential | B-sensitive? |
|--------|:---:|:---:|--------------|
| **DCP** | ✓ | DR / DP / DM | No — NCP = N·r²/2 is B-invariant |
| **JTK** | ✓ | — | Favors low B (mean-collapse artifact) |
| **RAIN** | ✓ | — | Yes — umbrella test benefits from distinct ZTs |
| **MH** | ✓ | — | Yes — adaptive K = ⌊(B−1)/2⌋ captures harmonics |
| **CircaCompare** | — | DP / DM | — |
| **LimoRhyde** | — | DR | — |
| **DODR** | — | DR | — |

---

## 2. Quick Start

### Load the framework

```r
old_wd <- setwd("code")
source("setup.R")
setwd(old_wd)
```

### Build option objects

Three objects drive every analysis:

```r
# Estimate biological parameters from pilot data
bio <- estCircadianParam(
  data    = pilot_expr,   # gene × sample expression matrix (log2-CPM)
  times   = pilot_times,  # time of day in hours (0–24)
  period  = 24,
  verbose = TRUE
)

# Study design
design <- CircadianDesignOptions(
  sample_sizes = c(20, 30, 40, 50, 60, 80, 100),
  nsims        = 200L,
  design       = "passive",   # "passive" = pilot TOD; "active" = equispaced
  cts          = bio$cts
)

# Multiple-testing settings
analysis <- CircadianAnalysisOptions(
  alpha           = 0.05,
  p.adjust.method = "BH",
  r_strata        = makeAdaptiveRStrata(bio, bin_width = 0.25)
)
```

---

## 3. Single-Cohort Power

`runSingleCohortPower()` sweeps N × B × α₂ × method and returns an `SCPSingleResult` object.

### Minimal example

```r
bio <- readRDS("data/gse160521_nac_ctrl_pilot.rds")

design <- CircadianDesignOptions(
  sample_sizes = c(20, 40, 60, 80, 100),
  nsims        = 200L,
  design       = "passive",
  cts          = bio$cts
)

analysis <- CircadianAnalysisOptions(
  alpha           = 0.05,
  p.adjust.method = "BH",
  r_strata        = makeAdaptiveRStrata(bio, bin_width = 0.25)
)

res <- runSingleCohortPower(
  bio.opts      = bio,
  design.opts   = design,
  analysis.opts = analysis,
  methods       = "DCP",
  mc.cores      = 4L,
  plot          = TRUE,
  output_file   = "figures/single_cohort_power.pdf"
)

print(res)        # n80 summary table
plot(res)         # power vs N, faceted by B × method
npower(res)       # interpolated N for 80% power
```

### Multi-method B vs m sweep

```r
design_bvsm <- CircadianDesignOptions(
  sample_sizes = seq(12, 96, by = 12),
  nsims        = 100L,
  design       = "active",
  cts          = seq(0, 20, by = 4),     # placeholder — overridden by B_values
  B_values     = c(3L, 4L, 6L, 8L, 12L)
)

# Print method guidance before running
printMethodGuidance(methods = c("DCP", "JTK", "RAIN", "MH"), verbose = TRUE)

res_multi <- runSingleCohortPower(
  bio.opts      = bio,
  design.opts   = design_bvsm,
  analysis.opts = analysis,
  methods       = c("DCP", "JTK", "MH"),
  alpha2        = c(0, 0.5, 1.0),        # cosinor violation sweep
  mc.cores      = 8L,
  plot          = FALSE
)

head(res_multi$power_df)
#   N  B alpha2 alpha3 method    power power_se
#  12  3      0      0    DCP   0.143     ...
```

### Cosinor violation parameter (α₂)

`alpha2` adds a second harmonic to the simulated signal:

$$y \sim M + A\bigl[\cos(\omega t - \phi) + \alpha_2 \cos(2\omega t - \phi)\bigr] + \varepsilon$$

| α₂ | Waveform | Effect on DCP | Effect on MH |
|----|----------|---------------|--------------|
| 0 | Pure sinusoid | Nominal power | Nominal power |
| 0.5 | Moderate non-sinusoidal | Unchanged | Slight gain |
| 1.0 | Equal 1st + 2nd harmonic | Unchanged | Clear gain at B ≥ 6 |

DCP is omnibus and B-invariant; MH benefits from both higher B (more harmonics fit) and higher α₂ (more signal in harmonics).

---

## 4. Differential Power

`runDifferentialPower()` sweeps N × α₂ × method × test_type.  
Pilot parameters come from `estCircadianParamTwoGroup()`.

### Method × endpoint support matrix

| Method | DR | DP | DM | Notes |
|--------|:--:|:--:|:--:|-------|
| DCP | ✓ | ✓ | ✓ | Full hierarchical pipeline |
| CircaCompare | — | ✓ | ✓ | Parametric fit |
| LimoRhyde | ✓ | — | — | limma interaction model |
| DODR | ✓ | — | — | Differential oscillation |

Unsupported combinations return `NA` silently.

### Example

```r
bio_diff <- readRDS("data/gse160521_nac_vs_putamen_ctrl_pilot.rds")

design <- CircadianDesignOptions(
  sample_sizes = c(20, 40, 60, 80, 100, 120),
  nsims        = 200L,
  design       = "passive",
  cts          = bio_diff$cts,
  test_types   = c("DR", "DP", "DM")
)

analysis <- CircadianAnalysisOptions(alpha = 0.05, p.adjust.method = "BH")

res_diff <- runDifferentialPower(
  bio.opts      = bio_diff,
  design.opts   = design,
  analysis.opts = analysis,
  methods       = "DCP",
  test_types    = c("DR", "DP", "DM"),
  mc.cores      = 4L,
  plot          = TRUE,
  output_file   = "figures/diff_power.pdf"
)

print(res_diff)
```

### Multi-method differential

```r
res_diff_multi <- runDifferentialPower(
  bio.opts      = bio_diff,
  design.opts   = design,
  analysis.opts = analysis,
  methods       = c("DCP", "LimoRhyde"),
  test_types    = "DR",
  mc.cores      = 8L,
  plot          = FALSE
)
```

---

## 5. B vs m Trade-off and Method Recommendation

`recommendDesign()` is the full B vs m orchestrator. It runs three steps:

1. **Guidance** — prints method B-sensitivity table
2. **Analytical** — DCP closed-form CircaPower estimate (fast baseline, B-invariant)
3. **Simulation** — calls `runSingleCohortPower()` or `runDifferentialPower()` over B × N grid

```r
bio <- readRDS("data/gse160521_nac_ctrl_pilot.rds")

design_bvsm <- CircadianDesignOptions(
  sample_sizes = seq(12, 96, by = 12),
  nsims        = 100L,
  design       = "active",
  cts          = seq(0, 20, by = 4),
  B_values     = c(3L, 4L, 6L, 8L, 12L)
)

analysis <- CircadianAnalysisOptions(alpha = 0.05, p.adjust.method = "BH")

rec <- recommendDesign(
  bio.opts       = bio,
  design.opts    = design_bvsm,
  analysis.opts  = analysis,
  methods        = c("DCP", "MH"),
  target_power   = 0.80,
  mode           = "single",
  run_simulation = TRUE,
  mc.cores       = 8L
)

print(rec)    # CircaPower estimate + simulation n80 per method × B
plot(rec)
```

### Reuse a previous result

Pass a previous `SCPSingleResult` or `SCPDiffResult` as `prior_result` to skip re-running:

```r
rec2 <- recommendDesign(
  bio.opts      = bio,
  design.opts   = design_bvsm,
  analysis.opts = analysis,
  methods       = c("DCP", "MH"),
  prior_result  = res_multi,    # SCPSingleResult from a previous call
  target_power  = 0.80,
  mode          = "single"
)
```

### Standalone method guidance

```r
printMethodGuidance(methods = c("DCP", "JTK", "RAIN", "MH"), verbose = TRUE)
```

```
Method    B-sensitive?   Why
------    ------------   ---
DCP       No             NCP = N·r²/2 for equispaced B≥3; B cancels analytically
JTK       Favors low B   MetaCycle averages replicates before ranking — always prefer m
RAIN      Yes (genuine)  Umbrella test uses individual obs; more distinct ZTs help
MH        Yes            Adaptive K = floor((B-1)/2); more harmonics captured at B≥6
```

**Design decision rule:**

- Planning to use **DCP**: B doesn't matter — minimize B (e.g., B = 4–6) and maximize m.
- Planning to use **JTK**: same conclusion but for a different reason (method artifact, not biology).
- Planning to use **RAIN** or **MH**: prefer B ≥ 6 when total N is fixed and SNR is low.

---

## 6. Bootstrap Uncertainty

`runBootstrapDesignGrid()` quantifies pilot uncertainty by re-sampling the pilot and re-estimating power at each draw, returning pointwise 95% CIs.

```r
boot_opts <- CircadianBootstrapOptions(
  design_vector = seq(0, 22, by = 2),   # 12 equispaced ZTs
  B_values      = 12L,
  N_values      = c(20, 40, 60, 80, 100),
  nboot         = 50L,
  nsims_inner   = 20L,
  design        = "active",
  seed          = 42L
)

boot_res <- runBootstrapDesignGrid(
  pilot_data    = pilot_expr,
  pilot_times   = pilot_times,
  boot.opts     = boot_opts,
  analysis.opts = analysis,
  bio_diff.opts = bio_diff,
  mode          = "differential",    # "single" or "differential"
  methods       = "DCP",
  test_types    = "DR",
  mc.cores      = 8L
)
```

Bootstrap CIs are informatively wide when:
- Pilot n < 20
- Pilot signal-to-noise r̃ < 0.5

In those cases, treat the median power estimate conservatively and consider additional pilot data before designing the full study.

---

## 7. Working with Result Objects

All runners return S3 objects with consistent `print`, `plot`, and `npower` methods.

```r
# SCPSingleResult
res                          # prints n80 summary table
plot(res)                    # power curves faceted by B × method × alpha2
npower(res, target = 0.80)  # interpolated N for 80% power

# SCPDiffResult
res_diff
plot(res_diff)
npower(res_diff, target = 0.80)

# SCPRecommendResult
rec
plot(rec)
```

### Access the raw power table

All result objects store a tidy data frame at `$power_df`:

```r
df <- res$power_df
# Columns: N, B, alpha2, alpha3, method, power, power_se

# Filter to MH at B=6
df[df$method == "MH" & df$B == 6, c("N", "power")]

# n80 for each method × B combination
res$n80_df
```

---

## 8. Complete Publication Example

```r
# 0. Load framework
old_wd <- setwd("code"); source("setup.R"); setwd(old_wd)
GLOBAL_SEED <- 2025L

# 1. Load pilot (pre-estimated with estCircadianParam)
bio <- readRDS("data/gse160521_nac_ctrl_pilot.rds")

# 2. Options
design <- CircadianDesignOptions(
  sample_sizes = c(20, 30, 40, 50, 60, 80, 100, 120, 150, 200),
  nsims        = 200L,
  design       = "passive",
  cts          = bio$cts
)
analysis <- CircadianAnalysisOptions(
  alpha           = 0.05,
  p.adjust.method = "BH",
  r_strata        = makeAdaptiveRStrata(bio, bin_width = 0.25)
)

# 3. Run
set.seed(GLOBAL_SEED)
res <- runSingleCohortPower(
  bio, design, analysis,
  methods     = "DCP",
  mc.cores    = as.integer(Sys.getenv("MC_CORES", "4")),
  plot        = TRUE,
  output_file = "output/single_cohort/figures/fig1_nac.pdf"
)

# 4. Save
saveRDS(res, sprintf("output/single_cohort/results/fig1_nac_%s.rds",
                     format(Sys.time(), "%Y%m%d_%H%M%S")))
```

### Environment variables

| Variable | Default | Purpose |
|----------|---------|---------|
| `MC_CORES` | `4` | Parallel cores |
| `SMOKE_TEST` | `""` | Set to `true` for fast debug run |
| `TISSUE` / `REGION` / `COMP` | `""` | Run a single tissue/region/comparison |

---

## 9. Function Reference

### Core runners

| Function | Returns | Purpose |
|----------|---------|---------|
| `runSingleCohortPower()` | `SCPSingleResult` | Single-group power: N × B × α₂ × method |
| `runDifferentialPower()` | `SCPDiffResult` | Two-group power: N × α₂ × method × test_type |
| `recommendDesign()` | `SCPRecommendResult` | Full B vs m orchestrator |
| `runBootstrapDesignGrid()` | list with CIs | Pilot uncertainty quantification |

### Pilot estimation

| Function | Purpose |
|----------|---------|
| `estCircadianParam()` | Single-group: A, σ, proportion rhythmic, pilot TOD |
| `estCircadianParamTwoGroup()` | Two-group: extends above with DR/DP/DM proportions |

### Options builders

| Function | Purpose |
|----------|---------|
| `CircadianDesignOptions()` | Sample sizes, nsims, design type, B_values |
| `CircadianAnalysisOptions()` | FDR threshold, adjustment method, r-strata |
| `CircadianBootstrapOptions()` | Bootstrap-specific: nboot, nsims_inner, seed |
| `makeAdaptiveRStrata()` | Adaptive r-strata breaks from pilot distribution |

### Utilities

| Function | Purpose |
|----------|---------|
| `printMethodGuidance()` | Print B-sensitivity table for chosen methods |
| `npower()` | Interpolate N for a target power from a result object |
| `prepCircadianData()` | Preprocess expression matrix (filter, log-transform) |
| `fitCosinorAll()` | Fit cosinor model gene-by-gene (returns A, σ, phase, p-value) |

### S3 methods

| Method | Purpose |
|--------|---------|
| `print.SCPSingleResult` | n80 summary per method × B |
| `plot.SCPSingleResult` | Power curves faceted by B × method × α₂ |
| `print.SCPDiffResult` | n80 per method × test_type |
| `plot.SCPDiffResult` | Power curves faceted by test_type × method |
| `print.SCPRecommendResult` | CircaPower + simulation n80 comparison |
| `plot.SCPRecommendResult` | B vs m heatmap + power curves |

---

## Publication Scripts

The `examples/publication/` directory contains ready-to-run scripts for each paper figure:

| Script | Figure | Description |
|--------|--------|-------------|
| `10_single_cohort_power.R` | Fig 1 (Seney) | Single-cohort, NAc/Caudate/Putamen |
| `14_single_cohort_gtex_ADR_LIV.R` | Fig 1 (GTEx) | Single-cohort, Adrenal Gland + Liver |
| `11_differential_power.R` | Fig 2 | Differential DR/DP/DM, NAc vs Putamen |
| `12_differential_power_gtex_ADR_LIV.R` | Fig 2 (GTEx) | Differential, Adrenal vs Liver |
| `15_bvsm_method_comparison.R` | Fig 3 | DCP/JTK/MH × 3 datasets × α₂=0 |
| `15b_bvsm_rain.R` | Fig 3 (RAIN) | RAIN B vs m, N≤48 |
| `15c_bvsm_rain_extended.R` | Fig 4 | RAIN + α₂ sweep, N≤72 |
| `08a_bootstrap_baboon.R` | Fig 6 | Bootstrap CI, Baboon LUN (n=12) |
| `08c_bootstrap_seney.R` | Fig 6 | Bootstrap CI, Seney CTL (n=60) |

All scripts respect the `SMOKE_TEST=true` environment variable for fast debug runs and `MC_CORES` for parallelism.
