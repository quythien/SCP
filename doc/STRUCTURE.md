# PowerSim Structure Reference

## Quick Reference

### Code Modules (7 files)

| File | Purpose | Key Functions |
|------|---------|---------------|
| **simulation.R** | Generate simulated data | `simCircadian()`, `simCircadianDiff()`, `simulate_circadian_data()` |
| **detection.R** | Detect rhythms & differential patterns | `fitCosinor()`, `DCP_DiffR2()`, `DCP_DiffPar()`, `DCP_Analyze()` |
| **power.R** | Calculate power | `comparePower()`, `comparePowerDiff()`, `power_single_condition()`, `power_two_group()` |
| **estimation.R** | Estimate parameters from pilot data | `estimate_circadian_params()` |
| **options.R** | Configure simulation options | `createSimOptions()`, `updateSimOptions()` |
| **runner.R** | Run multiple simulations | `runSims()`, `runSimsDiff()` |
| **utils.R** | Helper functions | Various utilities |

### Test Scripts (4 files)

| File | Purpose |
|------|---------|
| **test_simple.R** | Basic validation - single group |
| **test_differential.R** | Differential pattern tests (DR/DP/DA) |
| **test_diffcircadian_integration.R** | Validate DiffCircadian integration |
| **test_powersim_proper.R** | PROPER-style comprehensive validation |

### Analysis Scripts (6 files)

| File | Purpose |
|------|---------|
| **CircaPower_Theoretical_Foundation.R** | Explains λ = r² × n × d decomposition |
| **Concept_Explanation_Phase_ThreeHypotheses.R** | Explains DR/DP/DA patterns |
| **PowerSim_Critical_Challenges.R** | Identifies 7 unique challenges |
| **PROPER_vs_PowerSim_Parameters.R** | Parameter mapping documentation |
| **create_proper_style_figures.R** | Generate PROPER-style figures |
| **create_reference_figures.R** | Generate reference comparison figures |

## Usage Example

```r
# Set working directory
setwd("/path/to/PowerSim")

# Load core modules
source("code/options.R")
source("code/simulation.R")
source("code/detection.R")
source("code/power.R")
source("code/runner.R")

# Create simulation options
sim.opts <- createSimOptions(
  ngenes = 5000,
  lBaselineExpr = rnorm(5000, 5, 2),
  lOD = rnorm(5000, -1.0, 0.3),
  prop_rhythmic = 0.25,
  amplitude = function(n) pmax(rlnorm(n, log(0.5), 0.5), 0.05),
  sim.seed = 42
)

# Run simulations
results <- runSims(
  sample_sizes = c(12, 24, 36, 48),
  nsims = 50,
  sim.opts = sim.opts,
  design = "active"
)

# Calculate power
power <- comparePower(
  simOutput = results,
  stratify.by = "effectsize",
  strata = c(0, 0.5, 1.0, 1.5, 2.0, Inf)
)

# View results
summaryPower(power)
```

## Differential Analysis Example

```r
# Load differential functions (if not already loaded)
source("code/simulation.R")  # Contains simCircadianDiff
source("code/detection.R")   # Contains DCP_Analyze
source("code/power.R")       # Contains comparePowerDiff

# Simulate differential data
sim.diff <- simCircadianDiff(
  ngenes = 5000,
  n1 = 24, n2 = 24,
  lBaselineExpr = rnorm(5000, 5, 2),
  lOD = rnorm(5000, -1.0, 0.3),
  prop_DR = 0.10,
  prop_DP = 0.05,
  prop_DA = 0.05,
  design = "active"
)

# Detect differential patterns
dcp_result <- DCP_Analyze(
  expr1 = sim.diff$expr1,
  expr2 = sim.diff$expr2,
  times1 = sim.diff$times1,
  times2 = sim.diff$times2,
  alpha = 0.05
)

# Calculate power for each test type
power_DR <- sum(dcp_result$classification == "DR" &
                sim.diff$ground_truth$category == "DR") /
              sum(sim.diff$ground_truth$category == "DR")
```

## Running Tests

```bash
# Basic test
Rscript scripts/test_simple.R

# Differential patterns test
Rscript scripts/test_differential.R

# DiffCircadian integration test
Rscript scripts/test_diffcircadian_integration.R

# PROPER-style validation
Rscript scripts/test_powersim_proper.R
```

## Archive Locations

- **scripts_archive/** - Dataset-specific analyses (baboon, mouse, etc.)
- **analysis_archive/** - Exploratory scripts
- **code_backup/** - Original code files before merging

## Key Improvements from Restructuring

1. **No duplicate functions** - All implementations consolidated
2. **Clear module boundaries** - Each file has single responsibility
3. **Easier navigation** - 7 files instead of 19 in code/
4. **Better maintainability** - Changes isolated to specific modules
5. **Professional structure** - Ready for package development
