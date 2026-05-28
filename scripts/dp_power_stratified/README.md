# DP Power Analysis - Stratified by r = A/σ (15 bins)
## Purpose
PROPER-style stratified power analysis by signal-to-noise ratio (r = A/sigma).

## Scripts
- `dp_power_signal_stratified.R` - Main simulation (6 sizes × 15 strata × 50 sims)
- `replot_signal_stratified.R` - Replot with CI and FDC

## Parameters
- Sample sizes: 10, 20, 40, 60, 80, 100
- r bins: (0,0.25], (0.25,0.5], (0.5,0.75], (0.75,1], (1,1.25], ..., >5
- Simulations: 50 per scenario
- Target effect: 0.1

## Output
- `dp_power_signal_stratified.pdf` - 4 plots with CI and FDC
- `dp_power_signal_stratified_results.rds` - Saved results for replotting

## Key Finding
Power increases dramatically with r. At n=60:
- r < 0.5: < 25% power
- r = 0.75-1: > 90% power
- r > 1.25: ~100% power

## Your Study (n=62/74, ~9hr phase shift)
Expected marginal power: ~86%
