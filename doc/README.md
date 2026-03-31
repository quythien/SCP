# PowerSim: Simulation-Based Power Analysis for Circadian Rhythm Studies

A comprehensive R package for power analysis in circadian rhythm detection and comparison, combining CircaPower's circadian-specific approach with PROPER's simulation-based framework.

## Overview

**Target Publication**: Statistics in Medicine (improvement over CircaPower)

**Two Scenarios**:
1. **Single Condition**: Power to detect rhythmic genes
2. **Two-Group Comparison**: Power to detect rhythm differences between conditions (NOVEL)

## Key Features

- Semi-parametric simulation (like PROPER)
- Multiple detection methods (cosinor, JTK, RAIN, RJMCMC)
- Flexible noise models (Gaussian, t, negative binomial)
- Stratified power by effect size
- False discovery cost (FDC) metrics
- Two-group comparison for amplitude, phase, and rhythmicity differences

## Installation

```r
# From local directory
devtools::install_local("PowerSim")
```

## Quick Start

### Scenario 1: Single Condition

```r
# Estimate parameters from pilot data
params = estimate_circadian_params(pilot_data, times)

# Run power analysis
results = power_single_condition(
  params = params,
  times = seq(0, 48, by = 4),
  n_range = 1:6,
  methods = c("cosinor"),
  n_sim = 100
)

# Visualize
plot_power_curves(results)
```

### Scenario 2: Two-Group Comparison

```r
# Estimate from two pilot datasets
params_A = estimate_circadian_params(pilot_A, times)
params_B = estimate_circadian_params(pilot_B, times)

# Power for detecting amplitude differences
results = power_two_group(
  params_A = params_A,
  params_B = params_B,
  times = seq(0, 48, by = 4),
  diff_type = "amplitude",
  delta_A = 0.3,
  n_sim = 100
)

# Visualize
plot_two_group_heatmap(results$power_matrix, n_A_range, n_B_range)
```

## File Structure

```
PowerSim/
├── README.md              # This file
├── PLAN.md                # Overall project plan
├── MAPPING.md             # PROPER → PowerSim mapping
├── SIMULATION.md          # Semi-parametric simulation details
├── R/
│   ├── simulate_circadian.R   # Data simulation functions
│   ├── power_single.R         # Scenario 1: Single condition
│   ├── power_two_group.R      # Scenario 2: Two-group comparison
│   ├── estimate_params.R      # Parameter estimation from pilot data
│   └── utils.R                # Helper functions and visualization
└── vignettes/                 # Tutorials (to be created)
```

## What We Borrow from PROPER

1. **Semi-parametric simulation**: Resample parameters from pilot data
2. **Stratified power**: Power by effect size strata
3. **False discovery cost**: FP/TP ratio
4. **Targeted power**: Power for "meaningful" effects only
5. **Method comparison**: Evaluate multiple detection methods

## What We Add Beyond CircaPower

1. **Scenario 2**: Two-group comparison (NOVEL)
2. **Method comparison**: Not limited to cosinor
3. **Flexible noise**: Beyond Gaussian
4. **Waveform robustness**: Beyond pure sinusoid
5. **FDR context**: Multiple testing correction
6. **Empirical estimation**: Data-driven parameter priors

## Key Functions

### Parameter Estimation

```r
estimate_circadian_params(data, times, period = 24)
```

Estimates distributions of:
- M (mesor)
- A (amplitude)
- σ (noise)
- r = A/σ (effect size)
- φ (phase)

### Simulation

```r
simulate_circadian_data(G, n, times, params, noise_type, waveform)
simulate_two_group(G, n_A, n_B, times, diff_type, delta_A, delta_phi)
```

### Power Analysis

```r
power_single_condition(params, times, n_range, methods, n_sim)
power_two_group(params_A, params_B, times, n_A_range, n_B_range,
                diff_type, test_type, n_sim)
```

### Power Assessment

```r
stratified_power(pvals, ground_truth, effect_sizes, breaks)
targeted_power(pvals, ground_truth, effect_sizes, min_effect)
false_discovery_cost(TP, FP)
```

## Effect Size Definition

We use **r = A/σ** (amplitude-to-noise ratio) as the effect size:

| r value | Interpretation |
|---------|----------------|
| r < 0.5 | Weak rhythm (hard to detect) |
| 0.5 ≤ r < 1 | Moderate rhythm |
| 1 ≤ r < 2 | Strong rhythm |
| r ≥ 2 | Very strong rhythm |

## Two-Group Difference Types

### 1. Amplitude Difference (ΔA)
- H0: A_A = A_B
- Test: Wald test or likelihood ratio test

### 2. Phase Difference (Δφ)
- H0: φ_A = φ_B
- Test: Watson-Williams test (circular ANOVA)

### 3. Rhythmicity Change
- H0: Gene rhythmic in both or neither
- Test: Compare individual rhythm p-values

### 4. Joint Test
- H0: A_A = A_B AND φ_A = φ_B
- Test: F-test comparing full vs reduced model

## Validation

- Match CircaPower for Gaussian + sinusoid + single condition
- Apply to real Young vs Old comparisons
- Sensitivity to model misspecification

## References

1. PROPER: Wu H, Wang C, Wu Z. Comprehensive power evaluation for differential expression using RNA-seq. *Bioinformatics*. 2015.

2. CircaPower: Zong X, et al. Experimental design and power calculation for circadian rhythm detection. *Bioinformatics*. 2023.

3. Cosinor: Cornelissen G. Cosinor-based rhythmometry. *Theor Biol Med Model*. 2014.

## License

MIT

## Contact

[Your contact information]
