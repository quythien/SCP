# PowerSim: Motivation and Framework

## The Problem: No Statistical Framework for "Difference in Differences"

### Current Situation in Differential Circadian Analysis

When researchers perform differential circadian analysis (e.g., SCZ vs CONTROL):

```
Step 1: Run DCP_Analyze(expr_SCZ, expr_CONTROL)
Step 2: Get p-values for DR, DP, DA
Step 3: Report genes with p < 0.05
Step 4: ??? ← NO STATISTICAL VERIFICATION!
```

**Critical questions cannot be answered:**
- Are my findings reliable or just noise?
- Was my study adequately powered?
- Would a larger sample size find more genes?
- Are non-significant results truly negative or just underpowered?

### The "Difference in Differences" Problem

Differential rhythmicity is a **second-order comparison**:

```
Level 1: Is gene rhythmic in CONTROL?   (R² > 0)
Level 2: Is gene rhythmic in SCZ?        (R² > 0)
Level 3: Does R² differ between groups?  (R²_CTRL ≠ R²_SCZ)
```

This is a "difference in differences":
- Compare rhythm strength in group 1
- Compare rhythm strength in group 2
- Test if the difference is significant

**No analytical framework exists to:**
- Calculate power for this test
- Determine required sample sizes
- Verify if findings are reliable
- Distinguish underpowered from null results

---

## The Solution: Simulation-Based Power Analysis

### Why Simulation?

**Simulation provides flexibility to:**

1. **Model the complex data structure**
   - Circadian patterns (cosinor model)
   - Sampling times (irregular, cross-sectional)
   - Mean-variance relationships
   - Gene-specific effects

2. **Estimate parameters from pilot data**
   - Preserve real data characteristics
   - Gene-specific baseline expression
   - Gene-specific rhythmicity parameters
   - Realistic effect size distributions

3. **Calculate power empirically**
   - Simulate 100+ datasets
   - Apply detection method
   - Count true positives / false positives
   - Estimate actual power

4. **Validate findings**
   - "If effect is truly ΔR² = 0.04, what's my power?"
   - "Is my significant result reliable?"
   - "Would larger n find more genes?"

---

## PowerSim Framework: CircaPower + PROPER

### Theoretical Foundation (from CircaPower)

**Non-centrality parameter decomposition:**

```
λ = r² × n × d

where:
  r = A/σ (effect size: amplitude/signal-noise ratio)
  n = sample size
  d = design factor (accounts for sampling scheme)
```

**Power calculation:**

```
F ~ F(df1, df2, ncp=λ)
Power = P(F > F_crit | λ)
```

This provides:
- ✓ Theoretical power for detecting rhythmicity
- ✓ Effect size decomposition
- ✓ Sample size calculations

### Semiparametric Simulation (from PROPER)

**Key insight:** Don't assume normal distribution, preserve real data structure

**Parameter estimation from pilot data:**

```r
params <- estimate_circadian_params(
  data = pilot_expression,
  times = pilot_times
)

# Returns:
# - lBaselineExpr: Gene-specific baseline (from pilot)
# - lOD: Gene-specific overdispersion (from pilot)
# - amplitude: Gene-specific amplitudes (from pilot)
# - sigma: Gene-specific noise (from pilot)
# - prop_rhythmic: Proportion rhythmic (from pilot)
```

**Simulate realistic data:**

```r
sim_data <- simCircadianDiff(
  ngenes = 5000,
  n1 = 30, n2 = 30,
  lBaselineExpr = params$lBaselineExpr,
  lOD = params$lOD,
  amplitude = params$amplitude,
  prop_DR = 0.10,
  prop_DP = 0.05,
  prop_DA = 0.05
)
```

**Power analysis:**

```r
results <- runSimsDiff(
  sample_sizes = c(20, 30, 40, 50, 60),
  nsims = 100,
  sim.opts = params
)

power <- comparePowerDiff(
  simOutput = results,
  target_by = "effectsize",
  delta = c(0, 0.1, 0.2, 0.3, 0.4)
)
```

---

## What Makes PowerSim Unique

### Comparison with Existing Tools

| Tool | Purpose | Limitation |
|------|---------|------------|
| **CircaPower** | Rhythmicity power | Single group only, no differential |
| **RNASeqPower** | DE power | Not for circadian data |
| **PROPER** | DE power (RNA-seq) | Not for circadian patterns |
| **PowerSim** | Differential circadian power | ✓ Addresses gap |

### PowerSim Innovations

1. **Differential Rhythmicity Power**
   - DR (Differential Rhythmicity): R²₁ ≠ R²₂
   - DP (Differential Phase): φ₁ ≠ φ₂
   - DA (Differential Amplitude): A₁ ≠ A₂

2. **Semiparametric Framework**
   - Preserves pilot data structure
   - Accounts for irregular sampling
   - Gene-specific effects
   - Mean-variance relationships

3. **Flexible Sampling Designs**
   - Longitudinal (repeated measures)
   - Cross-sectional (different individuals)
   - Mixed designs
   - Irregular time points

4. **Theoretical + Empirical**
   - Theoretical power (λ = r² × n × d)
   - Empirical validation via simulation
   - Both approaches confirm findings

---

## Real-World Example: Putamen Data

### Without PowerSim:

```
Result: 316 genes with DR in SCZ vs CONTROL (p < 0.05)
Question: Is this reliable?
Answer: ???
```

### With PowerSim:

```
Step 1: Estimate parameters from pilot
  → Mean R²_CTRL = 0.098
  → Mean R²_SCZ = 0.056
  → ΔR² = 0.042 (small effect)

Step 2: Calculate power
  → Current n: (59, 28)
  → Achieved power: 14%
  → Status: UNDERPOWERED

Step 3: Interpret findings
  → 316 genes @ p < 0.05 (6.6%)
  → Zero genes @ q < 0.05 (FDR)
  → Explanation: Small effects + low power
  → Recommendation: Need n > 100 for 80% power

Step 4: Biological insight
  → 93% of rhythmic genes CONSERVED
  → Core circadian machinery intact
  → Only subset disrupted in disease
```

**PowerSim transforms ambiguity into actionable insight!**

---

## Summary

### The Problem
- No statistical framework for differential circadian analysis
- Cannot verify if findings are reliable
- Cannot determine required sample sizes
- Cannot distinguish underpowered from null results

### The Solution
**PowerSim = CircaPower + PROPER**

- **CircaPower**: Theoretical framework (λ = r² × n × d)
- **PROPER**: Semiparametric simulation framework
- **PowerSim**: Extended to differential circadian rhythms

### The Value
Simulation-based approach provides:
- ✓ Power calculations for complex tests
- ✓ Sample size planning
- ✓ Validation of findings
- ✓ Interpretation of negative results
- ✓ Flexible scenario testing

**This is what makes PowerSim a necessary tool for differential circadian analysis!**
