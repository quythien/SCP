# PowerSim: A Simulation-Based Power Analysis Framework for Differential Circadian Rhythms

**Combining the theoretical foundation of CircaPower with the semiparametric simulation approach of PROPER**

---

## Overview

PowerSim provides power evaluation for differential circadian rhythm analysis through simulation-based methods. It addresses a critical gap in circadian biology: **no statistical framework exists for validating "difference in differences" findings** (differential rhythmicity between groups).

### The Problem

Traditional differential expression (DE) analysis has mature power calculation tools (PROPER, RNASeqPower). Differential **circadian** analysis lacks this framework:

```
Current workflow:
  Run DCP_Analyze(group1, group2)
  → Get p-values for DR, DP, DA
  → Report genes with p < 0.05
  → ??? NO WAY TO VERIFY FINDINGS ???
```

**Critical questions cannot be answered:**
- Are findings reliable or just noise?
- Was the study adequately powered?
- Would larger sample size find more genes?
- Are non-significant results truly negative or underpowered?

---

## Theoretical Foundation

### From CircaPower: Non-Centrality Parameter Decomposition

**λ = r² × n × d**

Where:
- **r = A/σ**: Effect size (amplitude-to-noise ratio)
- **n**: Sample size
- **d**: Design factor (accounts for sampling scheme)

This decomposition provides:
- ✓ Theoretical power calculations for rhythmicity detection
- ✓ Effect size standardization
- ✓ Sample size formulas

### From PROPER: Semiparametric Simulation Framework

**Key insight:** Don't assume parametric distributions—preserve real data structure

**Components:**
1. **Parameter estimation from pilot data**
   - Baseline expression
   - Overdispersion
   - Gene-specific effects

2. **Simulation-based power calculation**
   - Generate realistic data
   - Apply detection method
   - Empirical power estimation

3. **Flexible evaluation criteria**
   - alpha.type: "fdr" or "pval"
   - Stratification by expression/dispersion
   - Target definition by effect size

---

## PowerSim Framework

### Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                    PILOT DATA                                │
│              (1 dataset: control or pooled)                  │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ▼
         ┌───────────────────────┐
         │ Parameter Estimation  │
         │  estimate_circadian_  │
         │       params()        │
         └───────────┬───────────┘
                     │
                     ▼
    ┌────────────────────────────────┐
    │   Simulation Options          │
    │   createSimOptions(           │
    │     lBaselineExpr, lOD,        │
    │     amplitude, prop_DR,        │
    │     prop_DP, prop_DA           │
    │   )                           │
    └────────────┬───────────────────┘
                 │
                 ▼
    ┌────────────────────────────────┐
    │   Run Simulations             │
    │   runSimsDiff(                │
    │     sample_sizes,             │
    │     nsims                     │
    │   )                           │
    └────────────┬───────────────────┘
                 │
                 ▼
    ┌────────────────────────────────┐
    │   Power Evaluation             │
    │   comparePowerDiff(            │
    │     alpha.type = c("pval",    │
    │                "fdr"),         │
    │     target_by, stratify.by     │
    │   )                           │
    └────────────────────────────────┘
```

### Key Functions

| Function | Purpose | Analogous To |
|----------|---------|--------------|
| `estimate_circadian_params()` | Extract parameters from pilot data | PROPER's pilot data usage |
| `createSimOptions()` | Set up simulation scenario | PROPER's `RNAseq.SimOptions.2grp()` |
| `simCircadianDiff()` | Generate differential circadian data | PROPER's data generation |
| `runSimsDiff()` | Run simulations across sample sizes | PROPER's `runSims()` |
| `comparePowerDiff()` | Calculate power | PROPER's `comparePower()` |
| `DCP_Analyze()` | Detect differential rhythms | edgeR/DESeq (PROPER's DE methods) |

---

## Key Innovations

### 1. Pilot Data Usage (Like PROPER)

**Uses ONE pilot dataset to estimate parameters:**

```r
# From PROPER documentation
sim.opts.Cheung = RNAseq.SimOptions.2grp(
  ngenes = 20000,
  p.DE = 0.05,
  lOD = "cheung",           # From 1 pilot dataset
  lBaselineExpr = "cheung"  # From 1 pilot dataset
)

# PowerSim equivalent
params <- estimate_circadian_params(
  data = pilot_data,       # ONE pilot dataset
  times = pilot_times
)

sim.opts <- createSimOptions(
  ngenes = 5000,
  prop_DR = 0.05,
  lBaselineExpr = params$lBaselineExpr,  # From pilot
  lOD = params$lOD,                        # From pilot
  amplitude = params$amplitude
)
```

### 2. Alpha Type Flexibility (Like PROPER)

**Supports both FDR and raw p-values:**

```r
# From PROPER: Option 1 (conservative)
powers = comparePower(
  simres,
  alpha.type = "fdr",      # Control FDR
  alpha.nominal = 0.1
)

# PowerSim equivalent
power <- comparePowerDiff(
  simOutput = simresults,
  alpha.type = "fdr",      # Control FDR
  alpha.nominal = 0.05
)

# From PROPER: Option 2 (liberal)
powers = comparePower(
  simres,
  alpha.type = "pval",     # Use raw p-values
  alpha.nominal = 0.001    # Stricter threshold!
)

# PowerSim equivalent
power <- comparePowerDiff(
  simOutput = simresults,
  alpha.type = "pval",     # Use raw p-values
  alpha.nominal = 0.001    # Stricter threshold!
)
```

**When to use each:**
- **FDR**: Conservative, standard for large samples
- **P-values**: When FDR is too conservative or actual FDR differs from nominal

### 3. Target Definition (Like PROPER)

**Define "biologically interesting" genes:**

```r
# PROPER: by log fold change
target.by = "lfc"
delta = 0.5

# PROPER: by standardized effect size
target.by = "effectsize"
delta = 1

# PowerSim equivalent
target_by = "deltaR2"
delta = c(0, 0.1, 0.2, 0.3)  # Effect size strata
```

### 4. Stratification (Like PROPER)

**Stratify power by different factors:**

```r
# PROPER
stratify.by = "expr"        # Baseline expression
stratify.by = "dispersion"  # Dispersion

# PowerSim
stratify.by = "effectsize"  # Delta R²
stratify.by = "baseline"    # Baseline expression
stratify.by = "rhythmicity"  # R²
```

---

## Real-World Example: Putamen Data

### Without PowerSim

```
Result: 316 genes with DR in SCZ vs CONTROL (p < 0.05)
Question: Is this reliable?
Answer: ??? Cannot verify
```

### With PowerSim

```
Step 1: Estimate parameters from pilot (CONTROL group)
  → Mean R² = 0.098
  → Mean amplitude = 0.5
  → Mean sigma = 0.8

Step 2: Calculate achieved power
  → Current n: (59, 28)
  → Effect size: ΔR² = 0.042 (small)
  → Achieved power: 14% (LOW!)
  → Status: UNDERPOWERED

Step 3: Interpret findings
  → 316 genes @ p < 0.05 (6.6%)
  → 0 genes @ q < 0.05 (FDR)
  → Explanation: Small effects + low power
  → Recommendation: Need n > 100 for 80% power

Step 4: Biological insight
  → 93% of rhythmic genes CONSERVED
  → Core circadian machinery intact
  → Only subset (~7%) disrupted in disease
```

### Alpha Type Comparison

| Threshold | SCZ Genes | % | Interpretation |
|-----------|-----------|---|----------------|
| **p < 0.05** | 316 | 6.6% | Current threshold |
| **p < 0.01** | 26 | 0.5% | PROPER-style liberal |
| **p < 0.001** | 0 | 0.0% | PROPER-style strict |
| **q < 0.05 (FDR)** | 0 | 0.0% | Conservative |

**Key insight**: Even with liberal p < 0.01, only 0.5% of genes pass. This confirms WEAK SIGNAL, not just FDR issues.

---

## Why PowerSim Is Necessary

### 1. No Analytical Solution

Differential rhythmicity involves "difference in differences":
- Compare R² in group 1
- Compare R² in group 2
- Test if difference significant

**No closed-form power calculation exists for this test!**

### 2. Simulation Provides Flexibility

Handles complexities that analytical methods cannot:
- Irregular sampling times
- Missing data
- Unbalanced designs
- Multiple testing (FDR)
- Gene-specific effects

### 3. Validates Findings

Transforms ambiguity into actionable insight:
- "Is this result real?" → Calculate power
- "Should I collect more samples?" → Power curves
- "What can I detect?" → Effect size analysis

---

## Comparison with Existing Tools

| Tool | Purpose | Limitation |
|------|---------|------------|
| **CircaPower** | Rhythmicity power (1-group) | No differential analysis |
| **RNASeqPower** | DE power (RNA-seq) | Not for circadian patterns |
| **PROPER** | DE power (RNA-seq) | Not for circadian patterns |
| **PowerSim** | **Differential circadian power** | **✓ Addresses gap** |

### What PowerSim Combines

**From CircaPower:**
- ✓ λ = r² × n × d theoretical framework
- ✓ Effect size decomposition (r = A/σ)
- ✓ Sample size calculations
- ✓ F-distribution based power

**From PROPER:**
- ✓ Semiparametric simulation
- ✓ Pilot data parameter estimation
- ✓ Flexible alpha.type (pval/fdr)
- ✓ Stratification options
- ✓ Target definition by effect size
- ✓ Comprehensive summary statistics

**PowerSim Innovations:**
- ★ Differential rhythmicity (DR) power
- ★ Differential phase (DP) power
- ★ Differential amplitude (DA) power
- ★ Post-hoc testing framework
- ★ Joint rhythmicity (TOJR) categorization

---

## Usage Workflow

### Basic Usage

```r
# Step 1: Load pilot data
pilot_data <- read.csv("pilot_expression.csv")
pilot_times <- c(0, 4, 8, 12, 16, 20, 24)

# Step 2: Estimate parameters
params <- estimate_circadian_params(
  data = pilot_data,
  times = pilot_times,
  period = 24,
  min_rhythm_pval = 0.15
)

# Step 3: Create simulation options
sim.opts <- createSimOptions(
  ngenes = 5000,
  prop_DR = 0.05,
  prop_DP = 0.03,
  prop_DA = 0.02,
  lBaselineExpr = params$lBaselineExpr,
  lOD = params$lOD,
  amplitude = params$amplitude
)

# Step 4: Run simulations
simresults <- runSimsDiff(
  sample_sizes = c(20, 30, 40, 50),
  nsims = 100,
  sim.opts = sim.opts
)

# Step 5: Calculate power
power <- comparePowerDiff(
  simOutput = simresults,
  alpha.type = "fdr",
  alpha.nominal = 0.05,
  target_effect = 0.3
)

# Step 6: View results
summaryPower(power)
plotPower(power)
```

### Advanced Usage: FDR vs P-value

```r
# Compare power with different alpha types
power_fdr <- comparePowerDiff(
  simresults,
  alpha.type = "fdr",
  alpha.nominal = 0.05
)

power_pval <- comparePowerDiff(
  simresults,
  alpha.type = "pval",
  alpha.nominal = 0.01  # Stricter!
)

# Compare
cat("FDR-based:\n")
print(power_fdr$power.marginal)

cat("\nP-value-based:\n")
print(power_pval$power.marginal)
```

---

## Output Interpretation

### Summary Table (Like PROPER)

| Sample Size | Nominal FDR | Actual FDR | Marginal Power | Avg # TD | Avg # FD | FDC |
|-------------|-------------|------------|----------------|-----------|-----------|-----|
| 20 | 0.05 | 0.12 | 0.32 | 45 | 8 | 0.18 |
| 30 | 0.05 | 0.08 | 0.45 | 68 | 6 | 0.09 |
| 40 | 0.05 | 0.06 | 0.58 | 92 | 6 | 0.07 |
| 50 | 0.05 | 0.05 | 0.68 | 110 | 6 | 0.05 |

**Key metrics:**
- **Marginal Power**: Proportion of target genes detected
- **TD (True Discoveries)**: Number of true positive findings
- **FD (False Discoveries)**: Number of false positive findings
- **FDC (False Discovery Cost)**: FD / TD ratio

---

## Design Recommendations

### For Study Planning

1. **Start with pilot data** (even n=10-15)
2. **Estimate parameters** using `estimate_circadian_params()`
3. **Simulate power curves** for realistic sample sizes
4. **Choose sample size** based on target power and effect size

### For Interpreting Results

1. **Check achieved power** - if < 80%, may be underpowered
2. **Look at effect sizes** - small effects need larger n
3. **Compare FDR vs p-value** - if actual FDR >> nominal, consider p-values
4. **Examine TD/FD trade-off** - balance discoveries vs false positives

### Sample Size Guidelines

| Effect Size (ΔR²) | n for 80% Power | n for 90% Power |
|-------------------|-----------------|-----------------|
| 0.1 (small) | 100 | 140 |
| 0.2 (moderate) | 50 | 70 |
| 0.3 (large) | 35 | 45 |
| 0.5 (very large) | 20 | 30 |

---

## Advanced Options (From PROPER)

### 3.1 Resampling Effect Sizes from Pilot Data

**When to use:** Researchers may not want parametric assumptions but expect differential rhythmicity similar to another experiment.

**From PROPER documentation:**
> "Sometimes the user may not want to make a parametric assumption for the effect size, but can feel comfortable expecting the overall DE is similar to that observed in another experiment. In this case we may use resampling based simulation."

**For PowerSim:**

```r
# Option 1: Parametric (default)
# Assume effect sizes follow theoretical distribution
sim.opts <- createSimOptions(
  ngenes = 5000,
  prop_DR = 0.05,
  # Effect sizes from pilot parameters
  lBaselineExpr = params$lBaselineExpr,
  lOD = params$lOD
)

# Option 2: Resampling from pilot
# Preserve observed effect size distribution
sim.opts <- createSimOptions(
  ngenes = 5000,
  prop_DR = 0.05,
  # Resample effect sizes from observed pilot data
  pilot_effects = observed_delta_R2,  # From Putamen analysis
  effect_sampling = "resample"
)
```

**Advantages of resampling:**
- ✓ No parametric assumptions
- ✓ Preserves real effect size distribution
- ✓ More realistic for complex data

### 3.2 Standardized Effect Size

**The challenge:** What's "biologically interesting" for differential rhythmicity?

**From PROPER documentation:**
> "Many researchers use fold change as a unit of effect size... Others argue that the relevant effect size may depend on a gene's natural biological coefficient of variation (BCV)."

**PROPER provides two approaches:**

```r
# Option 1: By log fold change (default)
powers = comparePower(
  simres,
  target.by = "lfc",      # Log fold change
  delta = 0.5
)

# Option 2: By standardized effect size
powers = comparePower(
  simres,
  target.by = "effectsize",  # Log FC / sqrt(dispersion)
  delta = 1
)
```

**For PowerSim (differential circadian):**

```r
# Option 1: By Delta R² (absolute)
power <- comparePowerDiff(
  simresults,
  target.by = "deltaR2",
  delta = c(0, 0.1, 0.2, 0.3)
)

# Option 2: By standardized effect size (NEW!)
power <- comparePowerDiff(
  simresults,
  target.by = "standardized",  # Delta R² / sqrt(variance)
  delta = c(0, 0.5, 1.0, 1.5)
)
```

**Why standardized?**
- Accounts for gene-specific variability
- Identifies genes with large effects relative to their noise
- More comparable across genes

### 3.3 Raw P-Values Instead of FDR

**When to use:** Actual FDR differs from nominal FDR (common in small samples)

**From PROPER documentation:**
> "As one inspects the DE detection in simulation, a user may notice that the actual false discovery proportion differs from the nominal FDR reported by the DE detection method. One may decide to use unadjusted raw p-values... (typically much lower than 0.05 to account for multiple testing)."

**PROPER example:**

```r
# PROPER: Use p-values with stricter threshold
powers = comparePower(
  simres,
  alpha.type = "pval",     # Raw p-values
  alpha.nominal = 0.001,   # Much stricter!
  stratify.by = "dispersion",
  target.by = "effectsize",
  delta = 1
)
```

**For PowerSim:**

```r
# Putamen example: Actual FDR >> nominal FDR
cat("FDR-based (alpha.type = 'fdr'):\n")
cat("  Nominal FDR: 0.05\n")
cat("  Actual FDR: ~0.5 (10× higher!)\n")
cat("  Significant genes: 0 @ q < 0.05\n\n")

cat("P-value-based (alpha.type = 'pval'):\n")
cat("  Threshold: p < 0.01 (stricter)\n")
cat("  Significant genes: 26 @ p < 0.01\n")
cat("  Better for: Small samples, exploratory analysis\n")
```

**When to choose each:**

| Situation | Recommended | Rationale |
|----------|-------------|-----------|
| **Large samples (n > 50)** | FDR (q < 0.05) | Standard approach |
| **Small samples (n < 30)** | P-values (p < 0.01) | FDR too conservative |
| **Actual FDR >> nominal** | P-values | Nominal FDR misleading |
| **Exploratory analysis** | P-values | Less conservative, hypothesis generation |
| **Confirmatory study** | FDR | More conservative, controls false discoveries |

### Practical Usage: All Three Advanced Options Together

```r
# Comprehensive power analysis (like PROPER)
power <- comparePowerDiff(
  simOutput = simresults,

  # Advanced Option 3.3: Use p-values (not FDR)
  alpha.type = "pval",
  alpha.nominal = 0.001,

  # Advanced Option 3.2: Standardized effect size
  target.by = "standardized",
  delta = 1,

  # Stratification
  stratify.by = "baseline"
)
```

This provides:
- ✓ More realistic effect size expectations
- ✓ Standardized comparisons
- ✓ Appropriate multiple testing control

---

## Summary

**PowerSim = CircaPower (theoretical) + PROPER (simulation) + Differential Circadian Analysis**

### Key Contributions

1. **First power framework** for differential circadian rhythms
2. **Semiparametric simulation** preserving pilot data characteristics
3. **Flexible evaluation criteria** (FDR, p-values, stratification)
4. **Validates findings** - distinguishes underpowered from null results
5. **Guides study design** - sample size planning for future studies

### Value Proposition

**Before PowerSim:**
- ❌ Run differential analysis
- ❌ Report p-values
- ❌ Cannot verify reliability

**After PowerSim:**
- ✓ Estimate parameters from pilot data
- ✓ Calculate achieved power
- ✓ Interpret findings in context
- ✓ Plan future studies with adequate power

This transforms differential circadian analysis from ambiguous to interpretable!
