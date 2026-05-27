# SCP Framework — Usage Tutorial

> **Note (2026-05): cosinor-only framework.** SCP now centers on a single
> detector entry point, `detect_cosinor(K)`: `K = 1` is the single-harmonic
> cosinor F-test and `K = 2` adds the 12-hour second harmonic. The
> frequency-modulated Mobius (FMM) model has been removed from the
> recommended framework; JTK_CYCLE, RAIN, and the multi-harmonic detector
> remain only for benchmarking. Any section below that describes FMM or
> `omega`-based waveform switching is legacy and does not reflect the
> current method.

## Date: 2026-04-22

This tutorial walks through the three main use cases with complete working R code.
All examples assume the working directory is the project root (where `code/setup.R` lives).

```r
source("code/setup.R")   # loads all framework functions
```

---

## Scenario A — Single-Cohort Rhythmicity Power

**When to use:** You have a single group (one tissue/condition) and want to know
how many samples are needed to detect rhythmic genes with adequate power.

`runSingleCohortPower()` estimates power from pilot data, sweeping N under the
empirical distribution of A, σ, and phase. Use `methods=` to choose the
detection algorithm; DCP (cosinor F-test + BH) is the default.

```r
# Load or estimate pilot parameters
bio <- estCircadianParam(
  data    = pilot_expr,   # genes × samples matrix (log2-CPM or similar)
  times   = pilot_times,  # collection time of day (hours, 0–24)
  period  = 24,
  verbose = TRUE
)

design <- CircadianDesignOptions(
  sample_sizes = c(20L, 30L, 40L, 60L, 80L, 100L),
  nsims        = 200L,
  design       = "passive",   # use pilot TOD distribution
  cts          = bio$cts
)

analysis <- CircadianAnalysisOptions(
  alpha           = 0.05,
  p.adjust.method = "BH",
  r_strata        = makeAdaptiveRStrata(bio, bin_width = 0.25)
)

set.seed(2025L)
res <- runSingleCohortPower(
  bio, design, analysis,
  methods     = "DCP",    # also: "JTK", "RAIN", "MH"
  plot        = FALSE,
  verbose     = TRUE,
  mc.cores    = 4L
)

# Save and plot separately
saveRDS(res, "output/single_cohort_power.rds")
plotSingleCohortPower(res, out_pdf = "output/single_cohort_power.pdf",
                      title = "My Cohort — Single-Cohort Power")
```

**Output:** a list with `$marginal_power` [N × nsims], `$strat_power` [N × r_strata × nsims],
`$pvalues` [N × genes × nsims], and related arrays. Pass directly to `plotSingleCohortPower()`.

### Cosinor violation sensitivity

To test robustness against non-sinusoidal waveforms, use the `alpha2` / `alpha3`
arguments (2nd and 3rd harmonic amplitudes relative to the fundamental):

```r
res_viol <- runSingleCohortPower(bio, design, analysis,
                                  methods = "DCP", alpha2 = 0.5, plot = FALSE)
```

`alpha2 = 0` (default) is a pure cosinor. `alpha2 = 0.5` adds a second harmonic
of half the amplitude — the kind of asymmetric waveform common in real tissues.
DCP power is insensitive to this violation; MH power increases with `alpha2`.

---

## Scenario B — Two-Group Differential Power (DR / DP / DM)

**When to use:** You have two groups (tissue A vs tissue B, disease vs control)
and want to detect genes that differ in rhythmicity (DR), peak timing (DP), or
baseline expression (DM) between groups.

### Pilot estimation

**Two separate pilots** (one pilot per group): use `estCircadianParamTwoGroup()`.
It derives `prop_DR`, `prop_DP`, `phase_diff`, and the group-2 mesor distribution
directly from the empirical between-group differences.

**One shared pilot**: use `estCircadianParam()` and specify differential
proportions manually with `updateBioOptions()`.

```r
# Two-pilot path (recommended when data exist for both groups)
bio_diff <- estCircadianParamTwoGroup(
  data_1  = expr_g1, times_1 = times_g1,
  data_2  = expr_g2, times_2 = times_g2,
  period  = 24,
  verbose = TRUE
)

design <- CircadianDesignOptions(
  sample_sizes = c(20L, 30L, 40L, 60L, 80L, 100L, 120L),
  nsims        = 200L,
  design       = "passive",
  cts          = bio_diff$cts,
  test_types   = c("DR", "DP", "DM")
)

analysis <- CircadianAnalysisOptions(alpha = 0.05, p.adjust.method = "BH")

set.seed(2025L)
res <- runDifferentialPower(
  bio_diff, design, analysis,
  methods    = "DCP",               # also: "CircaCompare", "LimoRhyde", "DODR"
  test_types = c("DR", "DP", "DM"),
  plot       = FALSE,
  verbose    = TRUE,
  mc.cores   = 4L
)

saveRDS(res, "output/diff_power.rds")
plotDiffPower(
  res_list    = list(res),
  comp_labels = "Group A vs Group B",
  endpoints   = c("DR", "DP", "DM"),
  out_pdf     = "output/diff_power.pdf"
)
```

**Output:** a list with `$fdr_DR` [genes × N × nsims], `$fdr_DP`, `$fdr_DM`,
`$diff_type` (list), `$effectsize` (list). Pass directly to `plotDiffPower()`.

### Biological interpretation of endpoints

- **DR** (differential rhythmicity): gene oscillates in one condition but not the other —
  e.g., clock-regulated genes disrupted in neurodegeneration.
- **DP** (differential phase): gene oscillates in both conditions but peak shifts —
  e.g., circadian misalignment, tissue-specific timing differences.
- **DM** (differential mesor): gene oscillates in both but baseline shifts —
  e.g., inflammatory state superimposed on an intact oscillation. Analogous to
  differential expression in RNA-seq.

### Comparing two cohorts side-by-side

Run `runDifferentialPower()` once per comparison and pass both results to `plotDiffPower()`:

```r
res_a <- runDifferentialPower(bio_a, design, analysis, plot = FALSE, mc.cores = 4L)
res_b <- runDifferentialPower(bio_b, design, analysis, plot = FALSE, mc.cores = 4L)

plotDiffPower(
  res_list    = list(res_a, res_b),
  comp_labels = c("Comparison A", "Comparison B"),
  endpoints   = c("DR", "DP", "DM"),
  out_pdf     = "output/fig2.pdf",
  width = 15, height = 30
)
```

---

## Scenario C — B vs m Trade-off and Method Recommendation

**When to use:** You have a fixed sample budget N and want to know:
*how should I spread samples across time points (B) vs replicates per time point (m = N/B)?*
And: *does my choice of detection method (DCP, JTK, RAIN, MH) matter?*

`recommendDesign()` runs in three steps automatically:
1. **Guidance** — prints method × B sensitivity table
2. **Analytical** — DCP closed-form CircaPower estimate
3. **Simulation** — power sweep across N × B × method

```r
bio <- readRDS("data/my_pilot.rds")

design_bvsm <- CircadianDesignOptions(
  sample_sizes = seq(12L, 96L, by = 12L),
  nsims        = 100L,
  design       = "active",
  cts          = seq(0, 20, by = 4),   # placeholder — overridden by B_values
  B_values     = c(3L, 4L, 6L, 8L, 12L)
)

analysis <- CircadianAnalysisOptions(alpha = 0.05, p.adjust.method = "BH")

rec <- recommendDesign(
  bio.opts      = bio,
  design.opts   = design_bvsm,
  analysis.opts = analysis,
  methods       = c("DCP", "JTK", "MH"),
  target_power  = 0.80,
  mode          = "single",
  run_simulation = TRUE,
  mc.cores      = 8L
)

print(rec)   # optimal B per method, n80 summary
plot(rec)    # power vs N faceted by B × method
```

**Method guidance (printed automatically):**

```
Method    B-sensitive?  Why
DCP       No            NCP = N·r²/2 for equispaced B≥3; B cancels analytically
JTK       Favors low B  MetaCycle averages replicates — always prefer m
RAIN      Yes (genuine) Umbrella test uses individual obs; more distinct ZTs help
MH        Yes           Adaptive K = floor((B-1)/2); more harmonics at B≥6
```

### Reusing previous simulation results

Pass a previous result via `prior_result` to skip re-running the simulation:

```r
rec2 <- recommendDesign(bio.opts = bio, design.opts = design_bvsm,
                         analysis.opts = analysis, methods = c("DCP", "MH"),
                         prior_result = rec$simulation, target_power = 0.80)
```

### Cosinor violation sweep

Test how method power changes when the true waveform is non-sinusoidal:

```r
rec_viol <- recommendDesign(
  bio.opts = bio, design.opts = design_bvsm, analysis.opts = analysis,
  methods  = c("DCP", "MH"),
  alpha2   = c(0, 0.5, 1.0),   # sweep second-harmonic amplitude
  mode     = "single", run_simulation = TRUE, mc.cores = 8L
)
```

---

## Scenario D — Bootstrap Uncertainty

**When to use:** Your pilot is small (n < 20) or you want honest confidence
intervals on power estimates, accounting for uncertainty in the pilot parameter
estimates themselves.

`runBootstrapDesignGrid()` resamples the pilot genes with replacement, re-estimates
parameters, and re-runs the power simulation at each draw.

```r
boot_opts <- CircadianBootstrapOptions(
  N_values      = c(20L, 40L, 60L, 80L, 100L),
  B_values      = 12L,
  nboot         = 50L,
  nsims_inner   = 20L,
  design        = "active",
  design_vector = seq(0, 22, by = 2),
  seed          = 42L
)

boot_res <- runBootstrapDesignGrid(
  pilot_data    = pilot_expr,
  pilot_times   = pilot_times,
  boot.opts     = boot_opts,
  analysis.opts = analysis,
  mc.cores      = 8L
)

summaryBootstrapDesignGrid(boot_res, test_type = "DR")
plotBootstrapDesignGrid(boot_res, test_type = "DR")
```

**Interpretation:** Bootstrap CIs are wide when pilot n < 20 or signal-to-noise
r̃ < 0.5. In those cases, treat the median power estimate conservatively and
consider collecting additional pilot data before committing to a study design.

---

## Gene Type Reference

| Type | Label | Rhythmic G1 | Rhythmic G2 | Differential feature |
|------|-------|-------------|-------------|----------------------|
| 0    | Arrhythmic | No | No | — |
| 1    | Rhythmic (same) | Yes | Yes | — |
| 2    | DR (G1 only) | Yes | No | Loss of rhythmicity in G2 |
| 3    | DR (G2 only) | No | Yes | Gain of rhythmicity in G2 |
| 4    | DP (phase shift) | Yes | Yes | Peak timing differs |
| 5    | DM (mesor shift) | Yes | Yes | Baseline differs |

**Proportion budget constraint:**
```
prop_DR + prop_DP + prop_DM  <=  prop_rhythmic
```

---

## Key Function Reference

| Task | Function | Notes |
|------|----------|-------|
| Estimate single-group pilot params | `estCircadianParam()` | |
| Estimate two-group pilot params | `estCircadianParamTwoGroup()` | |
| Build design options | `CircadianDesignOptions()` | `B_values` for B vs m |
| Build analysis options | `CircadianAnalysisOptions()` | `r_strata` for stratified curves |
| Adaptive r-strata | `makeAdaptiveRStrata(bio)` | |
| **Single-cohort power** | **`runSingleCohortPower()`** | methods=, alpha2=, plot= |
| **Differential power** | **`runDifferentialPower()`** | methods=, test_types=, alpha2= |
| **B vs m + method recommendation** | **`recommendDesign()`** | calls internal grid engine |
| Bootstrap uncertainty | `runBootstrapDesignGrid()` | for small pilots |
| Plot single-cohort results | `plotSingleCohortPower(res)` | 3-panel figure |
| Plot differential results | `plotDiffPower(list(res), ...)` | 6-panel figure |
| Adaptive r-strata breaks | `makeAdaptiveRStrata()` | |
| Method guidance table | `printMethodGuidance()` | |

### Internal engines (advanced use)

These are called by the user-facing functions above. Direct use is only needed
for custom post-processing of raw per-gene FDR arrays.

| Engine | Called by | Returns |
|--------|-----------|---------|
| `runSimsSingleCohort()` | `runSingleCohortPower()` | rich list: pvalues, strat_power, r_values_list |
| `runSimsDiff()` | `runDifferentialPower()` | rich list: fdr_DR, fdr_DP, fdr_DM, diff_type |

---

## Parallelization

All three user-facing runners support `mc.cores` (uses `parallel::mclapply`
internally — works on Linux/macOS; set `mc.cores = 1L` on Windows):

```r
n_cores <- max(1L, parallel::detectCores() - 1L)

res <- runSingleCohortPower(bio, design, analysis,
                             mc.cores = n_cores, plot = FALSE)

res_diff <- runDifferentialPower(bio_diff, design, analysis,
                                  mc.cores = n_cores, plot = FALSE)

rec <- recommendDesign(bio, design_bvsm, analysis,
                        methods = c("DCP", "MH"), mc.cores = n_cores)
```

On a SLURM cluster:
```r
mc.cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
```
