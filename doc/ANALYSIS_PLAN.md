# PowerSim — Analysis Plan for Publication
## Date: 2026-03-19

---

## Overview

The paper demonstrates a bootstrap-simulation framework for circadian study design, covering
three complementary questions (Q1–Q3) across three real datasets and two simulation layers.
The framework supports both two-stage (plug-in) and bootstrap-based design; the paper
emphasizes bootstrap while providing two-stage as an alternative.

---

## Research Questions

| Q | Question | Framework function |
|---|---|---|
| Q1 | What N achieves 80% power? | `runPowerAnalysis()` or bootstrap median from `runBootstrapDesignGrid()` |
| Q2 | How uncertain is that n80 given noisy pilot data? | `runBootstrapDesignGrid()` — 95% CI on n80 |
| Q3 | For fixed N, is more time-point coverage (B↑) or more replicates (m↑) better? | `runBootstrapDesignGrid()` sweeping `B_values` |
| Q4 | How robust is the framework when the true waveform deviates from cosinor? | Fourier simulation (see Layer 2) |
| Q5 | How does two-stage compare to bootstrap-based design? | Side-by-side: `runPowerAnalysis()` vs `runBootstrapDesignGrid()` |

**Note:** Q2 subsumes Q1 — the bootstrap median IS the two-stage point estimate, with
added uncertainty quantification. Q3 is answered by sweeping `B_values` in the design grid.

---

## Four Datasets — Two Crossed Comparisons

Four datasets are selected to cleanly separate two distinct effects:
- **GSE54651 vs Baboon**: same design (active, m=1), different r → isolates effect of signal strength
- **Baboon vs D1D2**: similar r range, different m → isolates effect of replicates per ZT
- **Human aging**: passive baseline — shows cost of uncontrolled TOD

| Dataset | Design | m/ZT | r | prop_DR | n80 | Unique lesson |
|---|---|---|---|---|---|---|
| GSE54651 LIV vs CER | Active | 1 | ~2.9 | ~25% | <24 | Strong r → tiny n80 despite modest prop_DR |
| Baboon LUN vs CER | Active | 1 | ~1.72 | ~41% | ~24–36 | Same design as GSE54651; lower r → higher n80 |
| Mouse D1D2 | Active | ~8 | ~0.66 | ~20% | ~80–150 | Multi-replicate → B vs m tradeoff; tight CI |
| Seney MDD vs Control ACC | Passive | N/A | ~0.55 | ~19% | ~280–300 | Passive baseline — sex heterogeneity adds design lesson |
| ↳ Seney male-only | Passive | N/A | ~0.82 | ~19% | ~200–220 | Sex stratification reduces n80 by ~30% vs combined |

### Dataset 3: Seney et al. — MDD vs Control ACC
**Story:** Passive design — uncontrolled post-mortem TOD. Shows passive design cost AND introduces
sex stratification as a partial remedy. Sex heterogeneity in circadian rhythmicity is documented
in ACC and in MDD biology; stratifying by sex reduces r penalty and cuts n80 by ~30%.

| Property | Value |
|---|---|
| Files | `data/ACC_RNA_filtered_normalized.csv` + `data/MD5_MetaData_1-15-25.xlsx` + `data/TOD.xlsx` |
| Load | `read_excel()` + `read.csv()` |
| Groups | Control (n=60) vs MDD (n=60); or male-only (n=30/group) |
| Design | **Passive** — irregular post-mortem TOD, range [0.5, 23.5]h |
| Signal (combined) | r_median ≈ 0.55, prop_DR ≈ 19% (p<0.05) |
| Signal (male-only) | r_median ≈ 0.82, prop_DR ≈ 19% (p<0.05) |
| Expected n80 (combined) | **~280–300** |
| Expected n80 (male-only) | **~200–220** |
| B_values in bootstrap | Fixed at 4L (passive design, B not identifiable) |
| Framework lesson | Passive design is costly regardless of biology; BUT sex stratification restores within-sex phase coherence, raising r from 0.55→0.82 and cutting n80 by ~30%. Practical recommendation: control known confounders even in passive designs. |
| Script | `examples/exploratory/05_multi_dataset_DR.R` (Section 3) + `07_seney_sex_DR.R` |
| Alt dataset | Human aging young vs old — see commented-out block in `05_multi_dataset_DR.R` |

### Dataset 2: Baboon Cross-Tissue — LUN vs CER
**Story:** Active design, single replicate per ZT. Controlled timing helps enormously.

| Property | Value |
|---|---|
| File | `data/CAMO_PRC_hmb.RData` (copied to local data/) |
| Load | `load()` → `baboon_withTOD$baboon`, `baboon_withTOD$tod` |
| Groups | LUN (lung, n=12) vs CER (cerebellum, n=12) |
| Design | **Active** — controlled ZT 0,2,...,22h; n=1 per ZT per tissue |
| Transform | `data.matrix(log2(mat + 1))` — raw CPM, data.frame |
| Signal | r_median(LUN) ≈ 1.72, prop_DR ≈ 41% (p<0.05) |
| Expected n80 | **~24–36** |
| B_values in bootstrap | c(4L, 6L, 8L, 12L) |
| Framework lesson | Active design with B=12 ZT achieves n80~24; shows value of temporal coverage |

### Dataset 3: Mouse GSE54651 — LIV vs CER
**Story:** Active design, m=1 per ZT, very strong signal. Contrasts with baboon (same design, lower r).

| Property | Value |
|---|---|
| File | `data/mice_GSE54651_CPM.RData` (copied to local data/) |
| Load | `readRDS()` — non-standard magic number, load() fails |
| Groups | LIV (liver, n=8) vs CER (cerebellum, n=8) |
| Design | **Active** — CT 22,28,...,64h (8 points, 2 days); use all 8 CT, period=24 |
| Transform | Already log2-scaled; wrap `data.matrix()` |
| Signal | r_median(LIV) ≈ 2.9, prop_DR ≈ 20–25% (p<0.05) |
| Expected n80 | **< 24** |
| B_values in bootstrap | c(4L, 8L) |
| Framework lesson | r dominates: strong signal → n80<24 despite prop_DR only ~25%; shows r² scaling |
| Gotcha | CT 22,28,...,64 folds to ZT {4,10,16,22}×2 — do NOT split by day (df=1 → r inflates to ~60) |

### Dataset 4: Mouse D1 vs D2 Cell Types
**Story:** Active design, multiple replicates per ZT. Adds B vs m tradeoff question.

| Property | Value |
|---|---|
| Files | `data/mouse_clinicalinfo_03082021_rmOutliers.csv` + `data/mouse_D1D2_logCPMfiltered_counts.csv` |
| Load | `read.csv()` for both; expression values are **raw counts** (not log despite filename) |
| Groups | D1 cell type (n=45) vs D2 cell type (n=46) |
| Design | **Active** — controlled ZT 2,6,10,14,18,22h; m≈7–8 replicates per ZT |
| Transform | `log2(sweep(expr_raw, 2, colSums(expr_raw), "/") * 1e6 + 1)` |
| Signal | r_median ≈ 0.66, prop_DR ≈ 18–22% (p<0.05) |
| Expected n80 | **~80–150** |
| B_values in bootstrap | c(3L, 6L) — every 8h vs every 4h |
| Framework lesson | Multi-replicate pilot gives narrow bootstrap CI; B=6 beats B=3 (more coverage > more replicates at r~0.66); direct B vs m tradeoff answer |

### Pairwise comparison structure

| Comparison | What it isolates |
|---|---|
| GSE54651 vs Baboon | Effect of **r** (same active m=1 design, r=2.9 vs 1.72) |
| Baboon vs D1D2 | Effect of **m** (similar r~0.6–1.7, m=1 vs m=8) |
| Any active vs Human | Effect of **passive vs active** design |

---

## Analysis Layer 1: Bootstrap Design Grid (main analysis)

**Script:** `examples/exploratory/05_multi_dataset_DR.R` (4 datasets: GSE54651 + baboon + D1D2 + human)

For each of the 4 datasets:
- `runBootstrapDesignGrid()` sweeping N and B
- Summary: power vs N (bootstrap median ± 95% CI ribbon)
- n80 ± CI horizontal bar chart
- B vs m tradeoff plot (active designs: baboon and D1D2 only)

**Three-way comparison figure:**
- Panel A: power curves for all 4 datasets on one plot (log-x optional)
- Panel B: n80 bar chart with CI whiskers
- Colors: human=firebrick, baboon=darkorange, GSE54651=steelblue, D1D2=forestgreen

---

## Analysis Layer 2: Fourier Robustness Simulation

**Purpose:** The cosinor model assumes a perfect sinusoidal waveform. Real circadian
expression often has non-sinusoidal features (harmonics, asymmetric peaks, sharp transitions).
This layer tests: if the truth deviates from cosinor, how does that affect power recommendations?
Does the optimal B and m choice change under waveform misspecification?

**Approach:**
- Generate synthetic data using Fourier synthesis (`fourier_sim.R` already in codebase)
  adding harmonics: y(t) = A₁cos(2πt/24 - φ₁) + A₂cos(4πt/24 - φ₂) + ε
  where A₂/A₁ = deviation parameter (0 = pure cosinor, 0.5 = strong harmonic)
- Fix N (e.g., N=60), sweep B ∈ {3, 4, 6, 8, 12} and m = N/B
- Run DCP test (cosinor-based) on Fourier-generated data
- Measure: does power degrade with harmonic deviation? Does optimal B shift?

**Expected findings:**
- At low deviation: cosinor framework works well regardless of B
- At high deviation: more time points (higher B) become more important because
  sparse B misses the harmonic peaks; more replicates (high m) don't help
- This motivates recommending higher B when waveform shape is uncertain

**Key output:** 2D heatmap or line plot: power vs deviation level, stratified by B

---

## Analysis Layer 3: Two-Stage vs Bootstrap Comparison

**Purpose:** Demonstrate that bootstrap subsumes two-stage (same point estimate, adds CI),
and show when they diverge (small pilot → two-stage overconfident; bootstrap honest).

**Two-stage approach:**
- Estimate (A, σ, prop_DR) once from pilot (plug-in)
- Run `runPowerAnalysis()` at each N → single power curve, no CI

**Bootstrap approach:**
- Resample pilot NBOOT times; estimate parameters each time
- Run power simulation for each bootstrap draw → distribution of power curves
- Median = point estimate (≈ two-stage); 95% CI = uncertainty quantification

**Comparison design:**
- Use same 3 real datasets
- Compare two-stage n80 vs bootstrap median n80 (should be similar)
- Show bootstrap CI width as function of pilot size:
  - Human aging: n_pilot=62 → moderate CI
  - Baboon: n_pilot=12 → wide CI (small pilot → two-stage most misleading here)
  - Mouse D1D2: n_pilot=45 → narrow CI (well-estimated pilot)
- Key message: when pilot is small (baboon n=12), two-stage gives a single number
  with false precision; bootstrap honestly reflects that uncertainty

**Script:** `examples/publication/02_calibration.R` already does this comparison on
synthetic data; extend to real datasets or add as a new section.

---

## Paper Focus and Framework Flexibility

The framework supports both approaches:

```r
# Two-stage (plug-in):
bio_opts  <- estCircadianParam(pilot_data, pilot_times, ...)
power_res <- runPowerAnalysis(bio.opts=bio_opts, design.opts=design_opts, ...)

# Bootstrap (recommended):
boot_res  <- runBootstrapDesignGrid(pilot_data, pilot_times, boot.opts=boot_opts, ...)
```

**Paper emphasis:** Bootstrap-based design (Q1+Q2+Q3 in one call). Recommended because:
1. Honest about pilot uncertainty — CI reflects how much n80 could shift with different pilot
2. Directly answers B vs m tradeoff (sweep B_values in single call)
3. Small-pilot setting (baboon n=12) — bootstrap CI is the only way to know the estimate is unreliable

**Two-stage role in paper:** Presented as the baseline approach; bootstrap shown to be strictly
better (same point estimate, adds CI at modest computational cost).

---

## Script Map

| Layer | Script | Status |
|---|---|---|
| Validation | `examples/publication/01_validation.R` | **Production running** (2026-03-30) |
| Calibration | `examples/publication/02_calibration.R` | **Production running** (2026-03-30) |
| Core pipeline — passive (human aging) | `examples/publication/03_power_core.R` | **Production running** (2026-03-30) |
| Core pipeline — active (Baboon LUN vs CER) | `examples/publication/03b_power_core_active.R` | **Production running** (2026-03-30) |
| Design grid | `examples/publication/04_power_design.R` | **Production running** (2026-03-30) |
| Method comparison (DCP vs CircaCompare) | `examples/publication/05_method_comparison.R` | **Production running** (2026-03-30); est. 14–42h for CircaCompare |
| Layer 1 — 4 datasets | `examples/exploratory/05_multi_dataset_DR.R` | Smoke-tested; all 4 sections complete |
| Layer 1 — sex-stratified | `examples/exploratory/07_seney_sex_DR.R` | Smoke-tested; combined n80~300, male-only n80~200 |
| Layer 1 — D1D2 only | `examples/exploratory/06_mouse_D1D2_DR.R` | Smoke-tested |
| Layer 2 — Fourier robustness | `examples/exploratory/07_fourier_robustness.R` | **Production running** (2026-03-30) |
| Layer 3 — 2-stage vs bootstrap (synthetic) | `examples/publication/02_calibration.R` | Smoke-tested (synthetic pilot n=30) |
| Layer 3 — 2-stage vs bootstrap (real data) | `examples/publication/08_two_stage_vs_bootstrap_realdata.R` | **Production running** (2026-03-30) |

---

## Immediate Next Steps

1. ✓ `05_multi_dataset_DR.R`: 4-dataset script complete (GSE54651, Baboon, Seney, D1D2)
2. ✓ `07_seney_sex_DR.R`: sex-stratified Seney analysis (combined vs male vs female)
3. ✓ `prepCircadianData()` integrated in all scripts; `POWERSIM_ROOT` / `DATA_HUMAN` env vars added
4. ✓ `07_fourier_robustness.R`: written and running
5. ✓ `08_two_stage_vs_bootstrap_realdata.R`: written and running
6. ✓ `03b_power_core_active.R`: written, bug-fixed, running
7. ✓ `05_method_comparison.R`: written, running (CircaCompare expected to finish in ~14–42h)
8. ✓ **Production run launched** — all 8 scripts running in parallel screen sessions (2026-03-30)
9. TODO: Collect and review output from all 8 production runs
10. TODO: Decide whether `05_method_comparison.R` results warrant inclusion in paper or supplement
