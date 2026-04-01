# PowerSim: Bootstrap-Based Power Analysis for Circadian Rhythm Studies

A simulation framework for sample size and study design in circadian differential expression experiments. Supports both two-stage (plug-in) and bootstrap-based design, using real pilot RNA-seq data to estimate parameters empirically.

## Overview

**Target Publication**: Statistics in Medicine

**Core Questions Answered**:
1. What N achieves 80% power? (`runBootstrapDesignGrid` bootstrap median)
2. How uncertain is that n80 given noisy pilot data? (bootstrap 95% CI)
3. Is more temporal coverage (B↑) or more replicates per ZT (m↑) better? (sweep `B_values`)
4. How robust is the framework when the true waveform deviates from cosinor? (Fourier simulation)

**Difference types detected**:
- **DR** — Differential Rhythmicity: rhythmic in one group, flat in the other
- **DP** — Differential Phase: same amplitude, shifted peak time
- **DA** — Differential Amplitude: same phase, different peak-to-trough swing

## Installation

```r
# Source all code files directly (no package install needed)
source_dir <- file.path(POWERSIM_ROOT, "code")
old_wd <- setwd(source_dir); source("setup.R"); setwd(old_wd)
```

## Quick Start

### Prepare pilot data

```r
# Input can be raw counts, CPM, or log2-normalized
pilot_log <- prepCircadianData(expr_matrix, times = pilot_times, input_type = "log2")
data_1 <- pilot_log$data[, group1_idx]
data_2 <- pilot_log$data[, group2_idx]
times_1 <- pilot_log$times[group1_idx]
times_2 <- pilot_log$times[group2_idx]
```

### Estimate parameters from pilot (two-group)

```r
bio_diff <- estCircadianParamTwoGroup(
  data_1 = data_1, data_2 = data_2,
  times_1 = times_1, times_2 = times_2,
  period = 24, min_rhythm_pval = 0.1, verbose = TRUE
)
# Reports: prop_DR, r_median per group, prop_rhythmic
```

### Bootstrap design grid (recommended)

Answers Q1, Q2, and Q3 in a single call:

```r
boot_opts <- CircadianBootstrapOptions(
  B_values     = c(4L, 8L, 12L),      # time point counts to sweep
  N_values     = c(24L, 48L, 72L, 96L, 120L),
  nboot        = 50L,                  # bootstrap draws
  nsims_inner  = 20L,                  # sims per draw
  design       = "active",             # "active" or "passive"
  seed         = 42L
)

boot_result <- runBootstrapDesignGrid(
  pilot_data   = data_1,
  pilot_times  = times_1,
  boot.opts    = boot_opts,
  bio_diff.opts = bio_diff,
  analysis.opts = CircadianAnalysisOptions(alpha = 0.05, p.adjust.method = "BH"),
  verbose      = TRUE
)
# Outputs: power median + 95% CI at each (N, B) combination; n80 with CI
```

### Two-stage approach (alternative, no CI)

```r
design_opts <- CircadianDesignOptions(
  sample_sizes = c(24L, 48L, 72L, 96L, 120L),
  nsims        = 50L,
  design       = "active",
  cts          = seq(0, 22, by = 2),   # 12 ZT points
  test_types   = "DR"
)

ts_result <- runTwoStagePower(
  pilot_data    = data_1,
  pilot_times   = times_1,
  design.opts   = design_opts,
  bio_diff.opts = bio_diff,
  analysis.opts = CircadianAnalysisOptions(alpha = 0.05, p.adjust.method = "BH"),
  test_type     = "DR", verbose = TRUE
)
```

### Compare approaches

```r
comparison <- compareDesignApproaches(
  two_stage_result  = ts_result,
  bootstrap_result  = boot_result,
  test_type         = "DR",
  target_power      = 0.80
)
# comparison$n80_two_stage, $n80_boot_median, $n80_boot_lo, $n80_boot_hi
plotDesignComparison(comparison, target_power = 0.80, panels = "A",
                     output_file = "comparison.pdf")
```

---

## Recommended Workflow

```
prepCircadianData()
    ↓
estCircadianParamTwoGroup()          ← empirical parameter estimation
    ↓
runBootstrapDesignGrid()             ← Q1 + Q2 + Q3 (recommended)
    or
runTwoStagePower()                   ← Q1 only, no CI (two-stage baseline)
    ↓
compareDesignApproaches()            ← side-by-side n80 + CI width
    ↓
plotDesignComparison()
```

**Why bootstrap over two-stage?**
- Same point estimate (bootstrap median ≈ two-stage)
- Adds honest uncertainty quantification via 95% CI
- Small pilots (e.g., baboon n=12) → two-stage gives a single number with false precision; bootstrap reveals how wide the CI really is
- Directly sweeps B in one call → answers B vs m tradeoff without extra scripts

---

## File Structure

```
PowerSim/
├── doc/
│   ├── README.md            # This file
│   ├── ANALYSIS_PLAN.md     # Publication plan, dataset table, script map
│   └── CHANGES.md           # Session-by-session log
├── code/
│   ├── setup.R              # Loads all framework code
│   ├── options.R            # CircadianBioOptions, CircadianDesignOptions, etc.
│   ├── runner.R             # runPowerAnalysis, runBootstrapDesignGrid, runTwoStagePower
│   ├── bootstrap_sim.R      # Bootstrap draw logic, compareDesignApproaches
│   ├── estimation.R         # estCircadianParam, estCircadianParamTwoGroup
│   ├── simulation.R         # simulate_circadian_data, simulate_two_group
│   ├── detection.R          # DCP pipeline, runSimsDiff
│   ├── fourier_sim.R        # runFourierDeviationPower
│   ├── design_comparison.R  # plotDesignComparison
│   └── utils.R              # Helpers, prepCircadianData, fitCosinorAll
├── data/                    # Real pilot datasets (not committed)
├── examples/
│   ├── publication/         # Final publication scripts (01–08)
│   └── exploratory/         # Layer 1 (05/06/07) and utility scripts
└── output/                  # All results (not committed)
```

---

## Key Functions

### Data preparation
```r
prepCircadianData(data, times, input_type = c("log2", "counts", "cpm"), ...)
```

### Parameter estimation
```r
# Single group (rhythmicity power)
estCircadianParam(data, times, period = 24, prop_DR = 0.1, min_rhythm_pval = 0.1, ...)

# Two groups (DR/DP/DA power) — use this for real pilot data
estCircadianParamTwoGroup(data_1, data_2, times_1, times_2, period = 24, min_rhythm_pval = 0.1, ...)
```

Both return a `CircadianBioOptions` object with empirical distributions of A, σ, r=A/σ, φ, mesor, and the proportions prop_DR/prop_DP/prop_DA estimated from the pilot.

### Design specification
```r
CircadianBootstrapOptions(B_values, N_values, nboot, nsims_inner, design, seed)
CircadianDesignOptions(sample_sizes, nsims, design, cts, test_types)
CircadianAnalysisOptions(alpha, p.adjust.method, fdr_thresholds, reference_n)
```

### Power analysis
```r
runBootstrapDesignGrid(pilot_data, pilot_times, boot.opts, bio_diff.opts, analysis.opts, ...)
runTwoStagePower(pilot_data, pilot_times, design.opts, bio_diff.opts, analysis.opts, ...)
runPowerAnalysis(bio.opts, design.opts, analysis.opts, ...)   # single-group only
```

### Comparison and visualization
```r
compareDesignApproaches(two_stage_result, bootstrap_result, test_type, target_power)
plotDesignComparison(comparison, target_power, panels, output_file)
```

### Fourier robustness (supplement)
```r
runFourierDeviationPower(bio_opts, design_opts, analysis_opts,
                          harmonic_grid = expand.grid(alpha2 = c(0, 0.25, 0.5), alpha3 = 0),
                          test_type = "DR", verbose = FALSE)
# harmonic_grid: alpha2 = 2nd harmonic amplitude / fundamental amplitude
#                alpha3 = 3rd harmonic amplitude / fundamental amplitude
# alpha2 = 0 → pure cosinor; alpha2 = 1.0 → near square-wave
```

---

## Effect Size

We use **r = A/σ** (amplitude-to-noise ratio):

| r value | Interpretation | Typical n80 (active B=12) |
|---------|----------------|--------------------------|
| r < 0.5 | Weak rhythm | > 200 |
| 0.5 ≤ r < 1 | Moderate rhythm | 80–200 |
| 1 ≤ r < 2 | Strong rhythm | 24–80 |
| r ≥ 2 | Very strong rhythm | < 24 |

---

## Difference Types

| Type | Null hypothesis | DCP test |
|------|----------------|----------|
| DR (Differential Rhythmicity) | R²_A = R²_B (same rhythmicity fraction) | `LR_deltaR2()` — LR test on R² |
| DP (Differential Phase) | φ_A = φ_B | `LRTest_diff_phase()` — LR test constraining equal phase |
| DA (Differential Amplitude) | A_A = A_B | `LRTest_diff_amp()` — LR test constraining equal amplitude |

All three use likelihood ratio tests from DiffCircaPipeline (DCP) with BH FDR correction.
Priority ordering: DR is tested first; DP and DA are applied only to genes not classified DR.

---

## Publication Scripts

| Script | Purpose | Status |
|--------|---------|--------|
| `examples/publication/01_validation.R` | Framework validation | Complete |
| `examples/publication/02_calibration.R` | Two-stage vs bootstrap (synthetic) | Complete |
| `examples/publication/03_power_core.R` | Human aging passive design | Complete |
| `examples/publication/03b_power_core_active.R` | Baboon LUN vs CER | Complete |
| `examples/publication/03c_power_core_mouse.R` | Mouse D1 vs D2 | Complete |
| `examples/publication/04_power_design.R` | Design grid (S3 complete) | Partial |
| `examples/publication/05_method_comparison.R` | DCP vs CircaCompare | Complete |
| `examples/publication/08a_bootstrap_baboon.R` | Two-stage vs bootstrap: Baboon | Run in parallel with 08b/c |
| `examples/publication/08b_bootstrap_d1d2.R` | Two-stage vs bootstrap: D1D2 | Run in parallel with 08a/c |
| `examples/publication/08c_bootstrap_seney.R` | Two-stage vs bootstrap: Seney | Run in parallel with 08a/b |
| `examples/publication/08d_bootstrap_summary.R` | Summary figure (run after 08a/b/c) | Run after 08a/b/c complete |
| `examples/exploratory/07a_fourier_mouse.R` | Fourier robustness: GSE54651 | Run in parallel with 07b/c |
| `examples/exploratory/07b_fourier_baboon.R` | Fourier robustness: Baboon | Run in parallel with 07a/c |
| `examples/exploratory/07c_fourier_d1d2.R` | Fourier robustness: D1D2 | Run in parallel with 07a/b |
| `examples/exploratory/07d_fourier_summary.R` | Fourier summary figure | Run after 07a/b/c complete |
| `examples/exploratory/07e_a_fourier_extreme_mouse.R` | Extreme harmonics: GSE54651 | Run in parallel with 07e_b/c |
| `examples/exploratory/07e_b_fourier_extreme_baboon.R` | Extreme harmonics: Baboon | Run in parallel with 07e_a/c |
| `examples/exploratory/07e_c_fourier_extreme_d1d2.R` | Extreme harmonics: D1D2 | Run in parallel with 07e_a/b |
| `examples/exploratory/07e_d_fourier_extreme_summary.R` | Extreme harmonics summary | Run after 07e_a/b/c complete |

---

## References

1. Wu H, Wang C, Wu Z. Comprehensive power evaluation for differential expression using RNA-seq. *Bioinformatics*. 2015. (PROPER framework)
2. Zong X, et al. Experimental design and power calculation for circadian rhythm detection. *Bioinformatics*. 2023. (CircaPower)
3. Cornelissen G. Cosinor-based rhythmometry. *Theor Biol Med Model*. 2014.
4. Thaben PF, Westermark PO. Detecting rhythms in time series with RAIN. *J Biol Rhythms*. 2014.
