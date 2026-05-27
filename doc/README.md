# SCP: Simulation-Based Power Analysis for Circadian Rhythm Studies

> **Note (2026-05): cosinor-only framework.** SCP now centers on a single
> detector entry point, `detect_cosinor(K)`: `K = 1` is the single-harmonic
> cosinor F-test and `K = 2` adds the 12-hour second harmonic. The
> frequency-modulated Mobius (FMM) model has been removed from the
> recommended framework; JTK_CYCLE, RAIN, and the multi-harmonic detector
> remain only for benchmarking. Any section below that describes FMM or
> `omega`-based waveform switching is legacy and does not reflect the
> current method.


SCP (Simulation-based Circadian Power) is an R framework for sample size planning and study design in **circadian transcriptomic** experiments. Given a small pilot RNA-seq dataset, it answers: *How many samples do I need, and how should I space my time points?*

It supports both a fast two-stage (plug-in) approach and a bootstrap approach that additionally quantifies how much uncertainty in the sample size estimate comes from having a noisy pilot.

---

## Background

Circadian studies measure gene expression at multiple times of day (Zeitgeber time, ZT) to detect 24-hour rhythms. Study design involves two interdependent choices:

- **N** — number of subjects per group
- **B** — number of distinct time points sampled

SCP helps answer:

| Question | Method |
|----------|--------|
| What N achieves 80% power? | Bootstrap median from `runBootstrapDesignGrid()` |
| How uncertain is that estimate given my pilot size? | Bootstrap 95% CI on n80 |
| Is denser time sampling (↑B) or more replicates per ZT (↑m) better? | Sweep `B_values` in one call |
| Is the cosinor model robust if real waveforms are non-sinusoidal? | `simCircadianSingleCohort2H()` (two-harmonic generator) + `detect_cosinor(K = 2)` |

### Detectable difference types

SCP uses the **DiffCircaPipeline (DCP)** likelihood ratio framework to detect three gene-level differential patterns between two groups (e.g., disease vs control, tissue A vs tissue B):

| Type | Meaning | Example |
|------|---------|---------|
| **DR** — Differential Rhythmicity | Gene oscillates in one group but is flat in the other | Clock-controlled gene lost in neurodegeneration |
| **DP** — Differential Phase | Same amplitude, but peak time shifts between groups | Circadian misalignment or disease-altered phasing |
| **DM** — Differential Mesor | Same rhythm, but baseline expression level differs | Inflammatory state elevating average gene expression |

### Effect size: r = A/σ

The key parameter is the amplitude-to-noise ratio **r = A/σ**:

| r | Rhythmicity | Typical n80 (active design, B=12) |
|---|-------------|----------------------------------|
| < 0.5 | Weak | > 200 |
| 0.5–1.0 | Moderate | 80–200 |
| 1.0–2.0 | Strong | 24–80 |
| ≥ 2.0 | Very strong | < 24 |

---

## Installation

No package installation needed. Source the framework directly:

```r
POWERSIM_ROOT <- "/path/to/SCP"   # set to your local clone path
setwd(POWERSIM_ROOT)
source_dir <- file.path(getwd(), "code")
old_wd <- setwd(source_dir); source("setup.R"); setwd(old_wd)
```

---

## Quick Start: End-to-End Example

This example uses **synthetic data** so it runs without any real dataset.
For a real RNA-seq pilot, replace Step 1 with `prepCircadianData()` on your expression matrix.

### Step 1 — Generate a synthetic pilot dataset

```r
set.seed(42)
pilot <- simulate_two_group(
  G         = 500,                      # 500 genes
  n_A       = 8, n_B = 8,              # 8 samples per group
  times     = seq(0, 22, by = 2),      # 12 ZT points (active design)
  diff_type = "rhythmicity",            # genes rhythmic in A, flat in B
  prop_diff = 0.20,                     # 20% of genes are true DR
  period    = 24
)
# Returns: pilot$data_A [500 x 8], pilot$data_B [500 x 8], pilot$times [12]
```

### Step 2 — Estimate parameters from pilot

```r
bio_diff <- estCircadianParamTwoGroup(
  data_1 = pilot$data_A,  data_2 = pilot$data_B,
  times_1 = pilot$times,  times_2 = pilot$times,
  period = 24, min_rhythm_pval = 0.01, verbose = TRUE
)
# Output (example):
#   prop_DR:   0.19    (19% genes estimated as rhythmically different)
#   r_median:  1.02    (amplitude-to-noise ratio in group A)
#   r_median:  0.98    (amplitude-to-noise ratio in group B)
```

### Step 3 — Bootstrap design grid (recommended)

Answers N, CI, and B vs m in a single call:

```r
boot_opts <- CircadianBootstrapOptions(
  B_values    = c(6L, 12L),              # compare 6 vs 12 ZT points
  N_values    = c(20L, 40L, 60L, 80L),  # sample sizes to evaluate
  nboot       = 50L,                     # bootstrap resamples of pilot
  nsims_inner = 20L,                     # power sims per resample
  design      = "active",
  seed        = 42L
)

boot_result <- runBootstrapDesignGrid(
  pilot_data    = pilot$data_A,
  pilot_times   = pilot$times,
  boot.opts     = boot_opts,
  bio_diff.opts = bio_diff,
  analysis.opts = CircadianAnalysisOptions(alpha = 0.05, p.adjust.method = "BH"),
  verbose       = TRUE
)
# Output (example):
#   B=6,  N=40:  power = 62% [55%, 69%]
#   B=6,  N=60:  power = 81% [74%, 87%]  ← n80 with B=6
#   B=12, N=40:  power = 74% [68%, 80%]
#   B=12, N=40:  power = 74% [68%, 80%]  ← n80 with B=12 = fewer subjects needed
#
#   n80 (B=6):  ~60  [95% CI: 52, 68]
#   n80 (B=12): ~40  [95% CI: 35, 46]   → denser time sampling saves ~20 subjects
```

### Step 4 — Two-stage comparison (optional baseline)

```r
design_opts <- CircadianDesignOptions(
  sample_sizes = c(20L, 40L, 60L, 80L),
  nsims        = 50L,
  design       = "active",
  cts          = seq(0, 22, by = 2),    # 12 ZT points
  test_types   = "DR"
)

ts_result <- runTwoStagePower(
  pilot_data    = pilot$data_A,
  pilot_times   = pilot$times,
  design.opts   = design_opts,
  bio_diff.opts = bio_diff,
  analysis.opts = CircadianAnalysisOptions(alpha = 0.05, p.adjust.method = "BH"),
  test_type     = "DR", verbose = TRUE
)
# Two-stage gives a single power curve with no CI.
# It typically overestimates power by 5–10 pp vs bootstrap at small N.
```

### Step 5 — Compare and visualize

```r
comparison <- compareDesignApproaches(
  two_stage_result = ts_result,
  bootstrap_result = boot_result,
  test_type        = "DR",
  target_power     = 0.80
)

cat("Two-stage n80:", comparison$n80_two_stage, "\n")
cat("Bootstrap n80:", comparison$n80_boot_median,
    "[", comparison$n80_boot_lo, ",", comparison$n80_boot_hi, "]\n")

plotDesignComparison(comparison, target_power = 0.80, panels = "A",
                     output_file = "my_design_comparison.pdf")
```

---

## Recommended Workflow

```
Your pilot RNA-seq data
        ↓
prepCircadianData()          prep: normalize, filter, align times
        ↓
estCircadianParamTwoGroup()  estimate: prop_DR, r, A, σ per gene
        ↓
runBootstrapDesignGrid()     design: power vs N, across B values, with CI
        ↓                    (answers Q1 + Q2 + Q3 simultaneously)
compareDesignApproaches()    compare: bootstrap median vs two-stage n80
        ↓
plotDesignComparison()       visualize: power curves + CI ribbon
```

**Why bootstrap over two-stage?**

The two-stage method estimates pilot parameters once and treats them as exact — this systematically overestimates power (by ~5–10 pp) and underestimates n80. The bootstrap resamples the pilot many times, propagating parameter uncertainty into the power estimate:

- Bootstrap median ≈ two-stage point estimate
- Bootstrap 95% CI reveals how reliable that estimate is
- With a small pilot (e.g., n=12), the CI can be wide — meaning the single two-stage number has false precision
- Sweeping `B_values` in one call directly answers the B vs m tradeoff

---

## Using Real Pilot Data

Replace the synthetic pilot in Step 1 with your own data:

```r
# expr_matrix: genes x samples (raw counts, CPM, or log2)
# times: numeric vector of ZT/TOD values, one per sample

pilot_log  <- prepCircadianData(expr_matrix, times = sample_times, input_type = "counts")
data_grp1  <- pilot_log$data[, group == "A"]
data_grp2  <- pilot_log$data[, group == "B"]
times_grp1 <- pilot_log$times[group == "A"]
times_grp2 <- pilot_log$times[group == "B"]

# Then continue from Step 2 above using data_grp1, data_grp2, times_grp1, times_grp2
```

**Active vs passive design:**

| Design | Description | Set in code |
|--------|-------------|-------------|
| Active | Controlled ZT — samples collected at fixed, evenly-spaced times | `design = "active"` |
| Passive | Uncontrolled TOD — post-mortem or convenience sampling | `design = "passive"` |

For passive designs, set `B_values` to a single value (the effective number of distinct TOD bins) and `cts` to the observed TOD distribution from your pilot.

---

## Key Functions

### Data preparation
```r
prepCircadianData(data, times, input_type = c("log2", "counts", "cpm"), ...)
# Normalizes and filters; returns $data (log2 matrix) and $times
```

### Parameter estimation
```r
# Two groups — use this for differential expression power
estCircadianParamTwoGroup(data_1, data_2, times_1, times_2,
                           period = 24, min_rhythm_pval = 0.01, ...)

# Single group — use this for rhythmicity detection power only
estCircadianParam(data, times, period = 24, min_rhythm_pval = 0.01,
                  prop_DR = 0.15, prop_DP = 0.10, prop_DM = 0.00, ...)
```
Both return a `CircadianBioOptions` object with per-gene empirical distributions of amplitude (A), noise (σ), effect size (r = A/σ), phase (φ), and the proportion of DR/DP/DM genes estimated from the pilot. A and σ are stored as paired gene-level tuples (`sigma_rhythmic`) so that downstream simulation preserves the empirical r = A/σ distribution via joint sampling.

### Options objects
```r
CircadianBootstrapOptions(B_values, N_values, nboot, nsims_inner, design, seed)
CircadianDesignOptions(sample_sizes, nsims, design, cts, B_values, test_types)
CircadianAnalysisOptions(alpha, p.adjust.method, fdr_thresholds, reference_n)
```

### Power analysis — main runners

```r
# Single-cohort rhythmicity power (Figs 1, 3, 4)
# One method per call. For multi-method B vs m sweeps use recommendDesign().
runSingleCohortPower(
  bio.opts, design.opts, analysis.opts,
  methods  = "DCP",          # one of "DCP" | "JTK" | "RAIN" | "MH"
  alpha2   = 0,              # 2nd-harmonic waveform deviation
  alpha3   = 0,
  mc.cores = 1L,
  plot     = FALSE,          # set TRUE to auto-plot, or call plotSingleCohortPower(res) post hoc
  output_file = NULL         # PDF path for auto-plot
)
# Returns a rich list: $marginal_power [N x nsims], $strat_power [N x r_strata x nsims],
# $pvalues [N x genes x nsims], $sample_sizes, $r_values_list.
# Plot with: plotSingleCohortPower(res, out_pdf = "fig.pdf")

# Two-group differential power (Fig 2)
# One method per call; test_types controls which of DR/DP/DM are evaluated.
runDifferentialPower(
  bio.opts, design.opts, analysis.opts,
  methods    = "DCP",        # one of "DCP" | "CircaCompare" | "LimoRhyde" | "DODR"
  test_types = c("DR","DP","DM"),  # NA returned silently for unsupported combinations
  alpha2     = 0,
  alpha3     = 0,
  mc.cores   = 1L,
  plot       = FALSE,
  output_file = NULL
)
# Returns a rich list: $fdr_DR [genes x N x nsims], $fdr_DP, $fdr_DM,
# $diff_type [list[nsims]], $effectsize [list[nsims]], $sample_sizes.
# Plot with: plotDiffPower(list(res), comp_labels = "A vs B", endpoints = c("DR","DP","DM"))

# Bootstrap uncertainty wrapper — works for BOTH single-cohort and differential (Fig 5/6)
# mode auto-detected: "single" if pilot_data_2 is NULL, "differential" otherwise.
runBootstrapDesignGrid(
  pilot_data, pilot_times,                   # group 1 (or only group)
  boot.opts, bio_diff.opts, analysis.opts,
  pilot_data_2  = NULL, pilot_times_2 = NULL,  # group 2 for differential mode
  mode       = NULL,     # "single" | "differential" (auto-detected)
  methods    = "DCP",
  test_types = c("DR","DP","DM"),
  alpha2     = 0,
  alpha3     = 0
)
# Returns SCPBootstrapResult with bootstrap CI on power curve and n80.
```

### Detection method support matrix

| Method | Single-cohort | DR | DP | DM | DA |
|--------|--------------|----|----|----|----|
| `"DCP"` | ✓ (detect_DCP) | ✓ | ✓ | — | ✓ |
| `"JTK"` | ✓ (detect_JTK) | — | — | — | — |
| `"RAIN"` | ✓ (detect_RAIN) | — | — | — | — |
| `"MH"` | ✓ (detect_MH) | — | — | — | — |
| `"CircaCompare"` | — | — | ✓ | ✓ | ✓ |
| `"LimoRhyde"` | — | ✓ | — | — | — |
| `"DODR"` | — | ✓ | — | — | — |

### B vs m design study

`recommendDesign()` is the full B vs m study orchestrator. It runs three steps:
1. Prints method guidance (statistical reason per method)
2. Runs CircaPower analytical estimate for DCP (fast, B-invariant, no simulation)
3. Optionally runs `runSingleCohortPower()` / `runDifferentialPower()` for all methods

```r
# Full study (analytical + simulation)
rec <- recommendDesign(
  bio.opts, design.opts, analysis.opts,
  methods        = c("DCP","JTK","RAIN","MH"),
  target_power   = 0.80,
  mode           = "single",       # "single" | "differential"
  run_simulation = TRUE,           # FALSE = CircaPower only (fast)
  prior_result   = NULL,           # pass rec$simulation from a previous recommendDesign() to skip re-running
  alpha2         = 0,
  mc.cores       = 60L
)
# Returns SCPRecommendResult:
#   $guidance       — method guidance table
#   $analytical_df  — CircaPower power at each N (DCP)
#   $simulation     — runSingleCohortGrid() result (SCPSingleResult with $power_df) or runDifferentialPower() rich list
#   $recommendation — optimal B and n_target per method

# Reuse previous simulation (avoid re-running):
prev <- recommendDesign(bio.opts, design.opts, analysis.opts,
                        methods = c("DCP","RAIN","MH"), run_simulation = TRUE)
rec  <- recommendDesign(bio.opts, design.opts, analysis.opts,
                         methods = c("DCP","RAIN","MH"),
                         prior_result = prev$simulation,  # simulation step skipped
                         run_simulation = FALSE)

# Analytical only (instant):
rec_fast <- recommendDesign(bio.opts, design.opts, analysis.opts,
                             run_simulation = FALSE)
```

**Key insight printed by `recommendDesign()`:**
If you plan to use DCP/JTK/LimoRhyde/DODR for your final analysis → B doesn't matter, invest in N.
If you plan to use RAIN or multi-harmonic regression → B=6–8 gives meaningful gains at fixed N.

`runSingleCohortPower()` also prints the method guidance table automatically when `verbose=TRUE`
via `printMethodGuidance()`, which is also available standalone.

### Post-hoc calls
```r
res <- runSingleCohortPower(..., plot = FALSE)
plotSingleCohortPower(res, out_pdf = "fig.pdf", title = "My Cohort")
saveRDS(res, "output/single_cohort_power.rds")   # save for replotting at any FDR

res_diff <- runDifferentialPower(..., plot = FALSE)
plotDiffPower(list(res_diff), comp_labels = "A vs B", endpoints = c("DR","DP","DM"))
saveRDS(res_diff, "output/diff_power.rds")

rec <- recommendDesign(...)
print(rec)                            # B vs m guidance + n80 summary
plot(rec)                             # power vs N figure
printMethodGuidance(methods = "RAIN") # guidance for a specific method only
```

### Comparison and visualization
```r
compareDesignApproaches(two_stage_result, bootstrap_result, test_type, target_power)
plotDesignComparison(comparison, target_power, panels, output_file)
```

### Fourier robustness (waveform sensitivity analysis)
```r
# Tests how power changes if the true waveform is not a pure sinusoid
simCircadianFMM(
  bio_opts, design_opts, analysis_opts,
  harmonic_grid = expand.grid(alpha2 = c(0, 0.25, 0.5), alpha3 = 0),
  test_type = "DR"
)
# alpha2: 2nd harmonic amplitude / fundamental (0 = pure sinusoid, 1.0 = near square-wave)
# For bulk RNA-seq, alpha2 rarely exceeds 0.3 — cosinor-DCP is robust in this range
```

---

## File Structure

```
SCP/
├── doc/
│   ├── README.md                  # This file — API reference and usage guide
│   ├── TUTORIAL.md                # Worked examples for all four use cases
│   └── DATASET_SURVEY.md          # Pilot dataset inventory
├── code/
│   ├── setup.R             # Entry point — sources all code below
│   ├── options.R           # Options constructors (CircadianBioOptions, etc.)
│   ├── runner.R            # runSingleCohortPower, runDifferentialPower, recommendDesign, runBootstrapDesignGrid
│   ├── bootstrap_sim.R     # Bootstrap logic, compareDesignApproaches, fitCosinorAll
│   ├── estimation.R        # estCircadianParam, estCircadianParamTwoGroup
│   ├── simulation.R        # simCircadianDiff — joint A-sigma sampling
│   ├── detection.R         # DCP pipeline + detect_JTK, detect_RAIN
│   ├── fourier_sim.R       # simCircadianFMM
│   ├── design_comparison.R # plotDesignComparison
│   ├── plot_single_cohort.R # plotSingleCohortPower (3-panel)
│   ├── plot_diff.R         # plotDiffPower (18-panel differential figure)
│   ├── npower.R            # npower(): interpolated n for target power at given FDR
│   ├── plot_dr_power.R     # Stratified DR power plotting
│   ├── plot_with_se.R      # SE-bar helpers + DR/DP 6-panel figures
│   ├── summarize_dcp_pairs.R # Cross-pair DCP summary table (supplementary)
│   └── utils.R             # prepCircadianData, helpers
├── data/                  # Pilot datasets (not committed — see ANALYSIS_PLAN.md)
├── examples/
│   ├── publication/       # Reproducible scripts for all paper figures
│   │   ├── fig1_bootstrap_comparison.R   # Fig 1: two-stage vs bootstrap (3-panel)
│   │   ├── fig3_multi_dataset_dr_power.R # Fig 3: DR power across 4 datasets
│   │   ├── 03b_power_core_active.R       # Baboon B vs m bootstrap grid
│   │   ├── 03c_power_core_mouse.R         # D1D2 B vs m bootstrap grid
│   │   ├── 03d_power_core_mouse_gse.R    # Mouse GSE B vs m + bm_tradeoff figure
│   │   ├── 08a_bootstrap_baboon.R        # Bootstrap: Baboon LUN vs CER
│   │   ├── 08b_bootstrap_d1d2.R          # Bootstrap: Mouse D1 vs D2
│   │   ├── 08c_bootstrap_seney.R         # Bootstrap: Human PFC (Seney)
│   │   ├── 08d_bootstrap_summary.R       # Summary: CI width across datasets
│   │   ├── 08e_bootstrap_mouse_gse.R     # Bootstrap: Mouse LIV vs CER (GSE)
│   │   ├── 09_dm_singlecohort_smoke.R    # Smoke: DM endpoint + single-cohort sim
│   │   ├── 10_single_cohort_power.R      # Single-cohort rhythmicity power (3 datasets)
│   │   ├── 11_differential_power.R       # Figure 2: differential power (DR/DP/DM)
│   │   ├── 12_differential_power_gtex_ADR_LIV.R  # Figure 2 variant: GTEx ADR vs LIV
│   │   ├── 13_supp_tissue_summary.R      # Supplementary tissue-level signal table
│   │   ├── 14_single_cohort_gtex_ADR_LIV.R       # Single-cohort power: GTEx ADR & LIV
│   │   ├── 15_bvsm_method_comparison.R   # B vs m: DCP / JTK / multi-harmonic (3 datasets)
│   │   ├── 15b_bvsm_rain.R              # B vs m: RAIN only (N ≤ 48, parallel companion)
│   │   └── bm_tradeoff_twostage.R       # Regenerates bm_tradeoff.pdf via two-stage DR
│   └── exploratory/       # Dataset-specific and sensitivity analyses
├── paper/
│   └── PowerSim/
│       ├── SCP.tex               # Main manuscript (working file)
│       ├── supplementary.tex     # Supplementary material
│       ├── references.bib        # Bibliography
│       └── figures/              # All publication figures (PDF)
└── output/                # All results (not committed)
```

---

## References

1. Wu H, Wang C, Wu Z. Comprehensive power evaluation for differential expression using RNA-seq. *Bioinformatics*. 2015. (PROPER framework — inspiration for bootstrap design)
2. Zong X, et al. Experimental design and power calculation for circadian rhythm detection. *Bioinformatics*. 2023. (CircaPower — single-condition baseline)
3. Cornelissen G. Cosinor-based rhythmometry. *Theor Biol Med Model*. 2014.
4. Thaben PF, Westermark PO. Detecting rhythms in time series with RAIN. *J Biol Rhythms*. 2014.
