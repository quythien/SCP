# PowerSim — Change Log

---

## Session: 2026-03-30

### New Scripts

**`examples/publication/03b_power_core_active.R`**
Companion to `03_power_core.R` for the active Baboon LUN vs CER dataset.
- Loads `data/CAMO_PRC_hmb.RData`; uses `prepCircadianData` with `input_type="cpm"`
- Uses `estCircadianParamTwoGroup(LUN, CER)` so that `prop_rhythmic` is computed from the
  union of both tissues — correctly handles DR genes rhythmic in CER but not LUN
- Sections: DR power (B=12), DP power (B=12), B vs m tradeoff (B=4/6/12)
- Env var: `POWERSIM_ROOT` only (no `DATA_HUMAN` — baboon data is inside the repo)
- Smoke: `POWERSIM_SMOKE=1 Rscript examples/publication/03b_power_core_active.R`

**`examples/publication/05_method_comparison.R`**
Standalone DCP vs CircaCompare power and FDR comparison, extracted from the commented-out
Section 10 of `03_power_core.R`.
- Compares: DP power (DCP LRT vs CircaCompare Wald NLS) and DM type I error
- Uses human PFC younger pilot (`DATA_HUMAN` env var or fallback path)
- `NGENES_CMP = 500`, `NSIMS_CMP = 20` — reduced from 5000 because CircaCompare is ~1–3 sec/gene
- Outputs: `dp_dcp.rds`, `dp_cc.rds`, `comparison_results.rds`, `.xlsx` tables, 3 PDF figures
- Smoke: `POWERSIM_SMOKE=1 Rscript examples/publication/05_method_comparison.R`

---

### Bug Fixes

**`examples/publication/03b_power_core_active.R` — `estCircadianParamTwoGroup` call**
- Original script passed `prop_DP`, `prop_DA`, `phase_diff`, `amp_diff` to
  `estCircadianParamTwoGroup()`, which does not accept these arguments (error: "unused arguments").
- Fix: removed those args from the `estCircadianParamTwoGroup()` call; they are applied
  per-section via `updateBioOptions()` (DR section sets `prop_DP=0`; DP section sets `prop_DR=0`).
- Base `opts_bio` now only calls `updateBioOptions(opts_bio, ngenes = NGENES_CORE)`.

**`examples/publication/03b_power_core_active.R` — prop_rhythmic constraint**
Earlier version computed `prop_DR_b <- mean(xor(rhy_lun, rhy_cer))` (XOR = 40.7%) then
passed it to `estCircadianParam()`, which estimated `prop_rhythmic` from LUN only (37.8%).
`updateBioOptions()` enforces `prop_rhythmic >= prop_DR`, causing a hard error.
- Root cause: XOR counts genes rhythmic in CER-only, which are not in LUN's rhythmic budget.
- Fix: replaced `estCircadianParam(mat_lun, ...)` with `estCircadianParamTwoGroup(mat_lun, mat_cer, ...)`
  so `prop_rhythmic` is computed from the union of both tissues (42.8%), satisfying the constraint.

---

### Server Production Run (2026-03-30)

All 8 publication scripts launched on the server in parallel `screen` sessions.
Pattern used (POWERSIM_ROOT inline; DATA_HUMAN for human-data scripts):

```bash
ROOT=/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim
DATA=$ROOT/data/combined_data.rds

screen -S validation    -dm bash -c "POWERSIM_ROOT=$ROOT Rscript examples/publication/01_validation.R > $ROOT/output/01_validation.log 2>&1"
screen -S calibration   -dm bash -c "POWERSIM_ROOT=$ROOT Rscript examples/publication/02_calibration.R > $ROOT/output/02_calibration.log 2>&1"
screen -S power_core    -dm bash -c "POWERSIM_ROOT=$ROOT DATA_HUMAN=$DATA Rscript examples/publication/03_power_core.R > $ROOT/output/03_power_core.log 2>&1"
screen -S active_core   -dm bash -c "POWERSIM_ROOT=$ROOT Rscript examples/publication/03b_power_core_active.R > $ROOT/output/03b_active_core.log 2>&1"
screen -S power_design  -dm bash -c "POWERSIM_ROOT=$ROOT DATA_HUMAN=$DATA Rscript examples/publication/04_power_design.R > $ROOT/output/04_power_design.log 2>&1"
screen -S method_cmp    -dm bash -c "POWERSIM_ROOT=$ROOT DATA_HUMAN=$DATA Rscript examples/publication/05_method_comparison.R > $ROOT/output/05_method_comparison.log 2>&1"
screen -S fourier       -dm bash -c "POWERSIM_ROOT=$ROOT Rscript examples/exploratory/07_fourier_robustness.R > $ROOT/output/07_fourier.log 2>&1"
screen -S bootstrap     -dm bash -c "POWERSIM_ROOT=$ROOT Rscript examples/publication/08_two_stage_vs_bootstrap_realdata.R > $ROOT/output/08_bootstrap.log 2>&1"
```

Note: `combined_data.rds` is present at `data/combined_data.rds` within the PowerSim root
(no need for an external path). Scripts 03 and 04 still expose `DATA_HUMAN` for portability.

| Screen | Script | Started | Notes |
|---|---|---|---|
| `validation` | 01 | 11:42 AM | Silent during null sim (5000 genes × 50 sims × 5 N) |
| `calibration` | 02 | 11:42 AM | Running |
| `power_core` | 03 | 11:41 AM | Running — passive, human aging |
| `active_core` | 03b | 11:42 AM | Running — active, Baboon |
| `power_design` | 04 | 11:40 AM | Running — bootstrap inner loops |
| `method_cmp` | 05 | 11:42 AM | Running — DCP fast; CircaCompare est. 14–42h |
| `fourier` | 07 | 11:42 AM | Running |
| `bootstrap` | 08 | 11:41 AM | Running |

---

## Session: 2026-03-19

---

## 1. Folder Reorganization: `examples/`

Scripts in `examples/` were reorganized into subfolders:

```
examples/
  publication/        ← 4 numbered scripts for paper (see below)
  exploratory/        ← informal analysis scripts
  figures/            ← standalone figure-generation scripts
  panelC160/          ← panel C sensitivity at n=160
```

Publication scripts were renamed with numeric prefixes:

| Old name              | New name                    | Purpose |
|-----------------------|-----------------------------|---------|
| `run_validation.R`    | `01_validation.R`           | Type I error, QQ plots, monotonicity, FDR control |
| `run_calibration.R`   | `02_calibration.R`          | Two-stage vs bootstrap calibration on synthetic data |
| `run_pipeline.R`      | `03_power_core.R`           | Core DR/DP/DA power from real pilot data |
| `run_pipeline2.R`     | `04_power_design.R`         | Bootstrap design grid + sensitivity on real pilot |

---

## 2. Code Fixes

### `examples/publication/01_validation.R`
- **Made self-contained**: removed hardcoded `load("output/run_20260219.../...")` path.
  Sections V3–V6 now run `runPowerAnalysis()` and `runPhaseShiftAnalysis()` internally.
- **Timestamped output**: `val_dir` now uses `format(Sys.time(), "%Y%m%d_%H%M")`.
- **QQ plots** instead of histograms: Section V2 now produces `-log10(p)` QQ plots with:
  - Genomic inflation factor λ = `median(χ²_obs) / qchisq(0.5, df=1)`; acceptable: 0.95–1.05
  - 95% pointwise confidence band from beta distribution
  - Red diagonal reference line (expected under null)
- **Fixed `j_60`**: changed hardcoded `j_60 <- 3` → `which.min(abs(val_design$sample_sizes - 60))`

### `examples/publication/02_calibration.R`
- **Removed Panel B**: `plotDesignComparison(...)` now called with `panels = "A"` to suppress
  the n80 bar chart (Panel B is redundant for calibration output).

### `examples/publication/04_power_design.R`
- **Fixed `cat()` format bug**: `cat("\nRunning bootstrap (fixed B=%d...)\n", S4_B)` →
  `cat(sprintf(..., S4_B))`.
- **Passive + multi-B guard**: `CircadianBootstrapOptions` now uses `B_values = 4L` (single B)
  for passive design. In passive mode, B does not affect which times are sampled, so sweeping
  B has no meaning. Section 3 sweeps N only.
- **Section 4 B**: changed from `modal optimal B` (which no longer applies with single B) to
  `S4_B <- 4L`.

### `code/design_comparison.R` — `plotDesignComparison()`
- Added `panels = c("A", "B")` argument (default = both panels).
- `panels = "A"` renders only the power-curve panel; skips the n80 bar chart.
- PDF width auto-adjusts: 7 inches for 1 panel, 12 inches for 2 panels.

### `code/bootstrap_sim.R`
- **Dimension-drop fix in `runBootstrapDesignGrid()`** (lines 342–345):
  `apply(boot_power_arr, c(2,3,4), FUN)` drops singleton dimensions when `n_B=1` or
  `n_tests=1`. Fixed by wrapping with `array(..., dim = c(n_N, n_B, n_tests))`:
  ```r
  power_mean  <- array(apply(..., mean),     dim = c(n_N, n_B, n_tests))
  power_se    <- array(apply(..., sd),        dim = c(n_N, n_B, n_tests))
  power_ci_lo <- array(apply(..., quantile), dim = c(n_N, n_B, n_tests))
  power_ci_hi <- array(apply(..., quantile), dim = c(n_N, n_B, n_tests))
  ```
- **Dimension-drop fix in `plotBootstrapDesignGrid()`** (lines 484–486):
  `result$power_mean[, , t_idx]` on an `[n_N × 1 × n_tests]` array drops the n_B dimension.
  Fixed by wrapping with `matrix(..., nrow = n_N, ncol = n_B)`.

### `code/estimation.R` — `estCircadianParam()`
- **Graceful rhythmic budget handling**: if `prop_DR + prop_DP + prop_DA > prop_rhythmic`
  (estimated from pilot data), the differential proportions are now silently scaled down
  proportionally to fit the budget, with an informative `message()`.
  Previously this caused a hard error from `CircadianBioOptions()`.
  This matters when: pilot data has a low fraction of rhythmic genes (e.g., ~20%) but the
  user specifies 15% DR + 10% DP = 25% — common in practice.

---

## 3. Smoke Test Status (2026-03-19)

| Script              | Status  | Notes |
|---------------------|---------|-------|
| `01_validation.R`   | PASSED  | QQ plots, monotonicity, FDR checks all green |
| `02_calibration.R`  | PASSED  | Panel A only (panels="A") |
| `03_power_core.R`   | PASSED  | DR/DP only; Section 10 (CircaCompare) removed (2026-03-19) |
| `04_power_design.R` | S3 PASSED; S4/S5 verified structurally | Long runtime (~30+ min at production settings) |

Quick smoke settings for `04_power_design.R`:
```r
S3_NBOOT       <- 3L
S3_NSIMS_INNER <- 3L
```

---

## 4. Key Design Decisions

### Passive design + single B
In passive mode, subjects are drawn from the empirical TOD distribution regardless of B.
B is not identifiable in passive mode. `04_power_design.R` uses `B_values = 4L` as a fixed
placeholder and sweeps only N. This is documented in comments in the script.

### `prop_rhythmic` budget
`CircadianBioOptions` enforces `prop_DR + prop_DP + prop_DA <= prop_rhythmic` because all
differential genes must be rhythmic in at least one group. When estimating from real pilot
data, the empirical rhythmic fraction may be less than the user-specified differential
proportions. `estCircadianParam` now handles this gracefully.

### Bootstrap CI interpretation
- The bootstrap CI in `runBootstrapDesignGrid` reflects **parameter-estimation uncertainty**
  (i.e., how much the power recommendation would shift if the pilot were resampled).
- It does NOT model within-subject variability or random effects.
- The model is a population-average cosinor; `sigma_g` captures pooled residual spread.

---

## 5. Files Changed (this session)

```
code/
  bootstrap_sim.R        — dimension-drop fix (array + matrix wraps)
  design_comparison.R    — panels argument added to plotDesignComparison
  estimation.R           — rhythmic budget auto-scaling in estCircadianParam

examples/publication/
  01_validation.R        — self-contained, QQ plots, timestamped output
  02_calibration.R       — panels="A" in plotDesignComparison call
  03_power_core.R        — Section 10 (DCP vs CircaCompare) removed; DA removed; DR+DP only
  04_power_design.R      — cat/sprintf fix, single B for passive design
```

---

## Session addendum: 2026-03-19 (dataset survey + new scripts)

### New exploratory scripts created

| Script | Status | Purpose |
|---|---|---|
| `examples/exploratory/baboon_tissue_scan.R` | Smoke-tested | Scan all 61 baboon tissues for circadian signal; select top DR pairs; run power |
| `examples/exploratory/dr_camo_power.R` | Smoke-tested | Baboon vs GTEx LUN (cross-species); Baboon LUN vs ILE (cross-tissue) DR power |
| `examples/exploratory/05_multi_dataset_DR.R` | Smoke-tested | Unified 3-dataset DR comparison with bootstrap design grid |
| `examples/exploratory/06_mouse_D1D2_DR.R` | Smoke-tested | Mouse D1 vs D2 cell type DR power; B vs m tradeoff |

### New docs created
- `doc/DATASET_SURVEY.md` — comprehensive dataset comparison, data paths, loading gotchas
- `doc/ANALYSIS_PLAN.md` — full publication analysis plan (3 datasets × 3 layers)

### Final dataset selection for publication (4 datasets)

Four datasets enabling two clean pairwise contrasts:

| Dataset | Design | Replicates/ZT | r_median | prop_DR | n80 |
|---|---|---|---|---|---|
| GSE54651 LIV vs CER | **Active** CT | 1 | ~2.9 | ~25% | <24 |
| Baboon LUN vs CER | **Active** ZT | 1 | ~1.72 | ~41% | ~24–36 |
| Mouse D1 vs D2 cell types | **Active** ZT | ~8 | ~0.66 | ~18–22% | ~80–150 |
| Human aging young vs old | **Passive** post-mortem | N/A | ~0.50 | ~14% | >300 |

Pairwise contrasts:
- GSE54651 vs Baboon → isolates **r** (same active m=1 design, different signal strength)
- Baboon vs D1D2 → isolates **m** (similar r range, single vs multi-replicate per ZT)
- Any active vs Human → isolates **passive vs active** design cost

Mouse D1D2 added (not replacing GSE54651) because it adds the multi-replicate dimension
and directly motivates the B vs m tradeoff question with real data.

### Three analysis layers planned

1. **Bootstrap design grid** — Q1+Q2+Q3 across all 3 datasets (`05_multi_dataset_DR.R`)
2. **Fourier robustness** — deviation from cosinor × B/m tradeoff (`07_fourier_robustness.R`, planned)
3. **Two-stage vs bootstrap** — same point estimate, bootstrap adds CI; small pilot → bootstrap honest (`02_calibration.R` extension)

### Key finding: CT day-by-day analysis (GSE54651)

Comparing all 8 CT vs day1-only vs day2-only for LIV:
- Day1/Day2 (n=4, df=1): r_median inflates to ~60 (overfitting artifact) → misleadingly high apparent power + FDR inflated to ~13%
- All 8 CT (df=5): r_median=2.89 (reliable) → proper FDR control (~4%)
- **Decision:** always use all 8 CT. Documented in Section 0 of `05_multi_dataset_DR.R`.

### 03_power_core.R correction
Section 10 was previously listed as "removed" but was actually **commented out** (not deleted).
Each line prefixed with `#`; header block added with "TO RUN: uncomment this section" note.

### Section 10 removal (`03_power_core.R`)
Section 10 (DCP vs CircaCompare method comparison) was removed entirely. Reason: CircaCompare
is ~1–3 sec/gene (NLS), making full runs infeasible. The comparison is not needed for the
core publication pipeline. Removed: `opts_design_compare`, `opts_bio_compare`,
`opts_analysis_DCP/CC`, `dp_dcp`, `dp_cc`, all helper functions (`.compute_comparison`,
`.compute_dm_typeI`, `.compute_stratified`), all ggplot2/openxlsx output, and the 6 PDF figures.
Wrap-up section cleaned of all `compare_base`/`compare_dir`/`xlsx_file` references.

---

## Session addendum: 2026-03-19 (prepCircadianData integration + bug fixes)

### New utility: `prepCircadianData()` in `code/utils.R`
Standard entry-point for loading and normalizing expression data before passing to any
framework function. Handles three input types:
- `"counts"` — raw integer counts → `log2(CPM + 1)`
- `"cpm"`    — raw CPM values → `log2(CPM + 1)`
- `"log2"`   — already log2-scale → `data.matrix()` only (no transformation)

Supports file-path loading (CSV/TSV), pheno-based time alignment via `sample_col`, and
automatic removal of samples with NA times. Returns `list(data, times, n_genes, n_samples)`.

### `prepCircadianData()` integrated into all analysis scripts

| Script | Dataset | `input_type` |
|---|---|---|
| `05_multi_dataset_DR.R` (Mouse) | GSE54651 log2-scale | `"log2"` |
| `05_multi_dataset_DR.R` (Baboon) | raw CPM data.frame | `"cpm"` |
| `05_multi_dataset_DR.R` (Seney) | log2-scale ACC | `"log2"` (after metadata join) |
| `05_multi_dataset_DR.R` (D1D2) | raw counts CSV | `"counts"` + `sample_col="sample"` |
| `06_mouse_D1D2_DR.R` | raw counts CSV | `"counts"` + `sample_col="sample"` |
| `07_seney_sex_DR.R` | log2-scale ACC | `"log2"` (after metadata join) |
| `03_power_core.R` | human aging log2 | `"log2"` (after pheno alignment) |
| `04_power_design.R` | human aging log2 | `"log2"` (after pheno alignment) |

For Seney (multi-file Excel + HU_NUM matching) and human aging (COMBINED object + pheno join),
`prepCircadianData` is called after the metadata preprocessing is done, on the already-subset
matrix, to provide validation (NA time removal, range check, coercion to numeric matrix).

### Server portability: `POWERSIM_ROOT` + `DATA_HUMAN` env vars
All five scripts now use the `POWERSIM_ROOT` env var pattern (was hardcoded in 06/03/04).
`03_power_core.R` and `04_power_design.R` also expose `DATA_HUMAN` for the path to
`combined_data.rds` (which lives outside the PowerSim root).

```bash
export POWERSIM_ROOT=/path/to/PowerSim
export DATA_HUMAN=/path/to/PNAS_aging/data/combined_data.rds
Rscript examples/exploratory/05_multi_dataset_DR.R
```

### Bug fixes

**`code/runner.R` — suppress NaN r-strata in verbose output**
Empty r-strata (no pilot genes in that r range) printed `NaN%` in the per-N progress line.
Fixed: `if (!is.nan(mean_p))` guard added before `cat(sprintf(...))`.
Output now stops at the highest occupied stratum rather than printing empty high-r buckets.

**`examples/exploratory/07_seney_sex_DR.R` — `utils.R` missing from source list**
`prepCircadianData` (and other utilities) live in `code/utils.R`. This file was not in the
manual `src_files` list, causing "could not find function" errors.
Fixed: `"code/utils.R"` added between `options.R` and `estimation.R` in `src_files`.

---

## Session addendum: 2026-03-19 (Layer 2 + Layer 3 real-data scripts)

### New scripts written

**`examples/exploratory/07_fourier_robustness.R`** — Layer 2: Fourier waveform robustness

Applies `runFourierDeviationPower()` (already in `code/fourier_sim.R`) to 3 real active pilots.
For each dataset (GSE54651, Baboon, D1D2), sweeps harmonic grid (α₂, α₃) at two B values.

Sections:
- S1–S3: per-dataset Fourier run (B_low vs B_high); 3-panel PDF: heatmap B_low, heatmap B_high, power vs α₂ comparison line
- S4: cross-dataset summary — power vs α₂ at fixed N for all 3 datasets and both B levels

Answers:
- Q1: How much does DR power degrade as harmonics increase? (α₂=0 → 0.75)
- Q2: Does higher B (denser time coverage) protect against harmonic misspecification?

Note: Seney (passive) excluded — B is not identifiable in passive designs.

Smoke: `POWERSIM_SMOKE=1 Rscript examples/exploratory/07_fourier_robustness.R`

---

**`examples/publication/08_two_stage_vs_bootstrap_realdata.R`** — Layer 3: real-data CI width comparison

Extends `02_calibration.R` (synthetic only) to three real pilot datasets with different n:

| Dataset | Pilot n | Design | Expected CI width |
|---|---|---|---|
| Baboon LUN | 12 | Active, B=12 | Wide (small pilot → high uncertainty) |
| Mouse D1D2 D1 | 45 | Active, B=6 | Medium |
| Seney CTL ACC | 60 | Passive, B=4 | Narrower |

Note: n=62 (human aging) was the original plan; replaced by Seney CTL n=60 since Seney is
now the primary passive dataset in the publication plan (human aging retained as alternative).

For each dataset: `runTwoStagePower()` + `runBootstrapDesignGrid()` (single B) → `compareDesignApproaches()`.

Sections:
- S1–S3: per-dataset comparison PDF (two-stage solid vs bootstrap dashed ± CI)
- S4: summary figure — Panel A: all 3 power curves overlaid; Panel B: mean CI width vs pilot n bar chart

Key message: two-stage n80 ≈ bootstrap median n80 (same point estimate), but bootstrap CI
width ∝ 1/√n_pilot. Baboon (n=12) CI is wide — two-stage gives false precision there.
Bootstrap is the honest quantification of pilot-induced uncertainty.

Smoke: `POWERSIM_SMOKE=1 Rscript examples/publication/08_two_stage_vs_bootstrap_realdata.R`

### Smoke test status (2026-03-19, this addendum)

| Script | Status |
|---|---|
| `07_seney_sex_DR.R` | PASSED — all 3 scenarios (combined / male / female), clean exit |
| `05_multi_dataset_DR.R` | Previously smoke-tested; prepCircadianData changes are non-transforming for log2/cpm (cosmetic coercion only); counts path identical to manual normalization |
| `06_mouse_D1D2_DR.R` | prepCircadianData replaces 4 lines; equivalent output |
