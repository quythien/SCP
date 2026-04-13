# SCP Framework — Usage Tutorial
## Date: 2026-04-13

This tutorial walks through six common use cases with complete working R code.
All examples assume the working directory is the project root (where `code/setup.R` lives).

```r
source("code/setup.R")   # loads all framework functions
```

---

## Scenario A — Single-Cohort Rhythmicity: Closed-Form Power

**When to use:** You have a single group (one tissue/condition), want to detect
rhythmic genes, and your pilot has clean cosinor parameter estimates (r, n, B, period).

The framework wraps the `CircaPower()` closed-form formula gene-by-gene and
averages across the transcriptome-wide r distribution estimated from a pilot.

```r
# --- Closed-form: runPowerAnalysis() with synthetic parameters ---

set.seed(1L)

bio_opts <- CircadianBioOptions(
  ngenes        = 1000L,
  prop_rhythmic = 0.30,                       # 30% of genes are rhythmic
  lBaselineExpr = rnorm(1000, 5, 1),          # log-scale baseline expression
  lOD           = rnorm(1000, -0.5, 0.3),     # log residual SD
  amplitude     = abs(rnorm(300, 0.6, 0.2)) + 0.05,
  sigma_rhythmic = NULL,                      # use lOD for both groups
  phase         = "uniform",
  prop_DR       = 0,                          # no differential (single cohort)
  prop_DP       = 0,
  prop_DM       = 0,
  period        = 24L,
  sim.seed      = 1L
)

# Active design: 6 ZT points, balanced
design_opts <- CircadianDesignOptions(
  sample_sizes = c(12L, 24L, 36L, 48L),
  nsims        = 50L,
  design       = "active",
  cts          = rep(seq(0, 20, by = 4), each = 8),  # 48 samples, 6 ZTs
  test_types   = "DR"
)

analysis_opts <- CircadianAnalysisOptions(alpha = 0.05, p.adjust.method = "BH")

# runPowerAnalysis uses the closed-form CircaPower formula per gene
power_res <- runPowerAnalysis(bio.opts = bio_opts,
                               design.opts = design_opts,
                               analysis.opts = analysis_opts)

# Summary: mean power at each N
cat(sprintf("Power at N=12: %.1f%%\n", 100 * power_res$summary$power_mean[1]))
cat(sprintf("Power at N=48: %.1f%%\n", 100 * power_res$summary$power_mean[4]))
```

**Output:** a list with `summary` (data frame of N vs power_mean/power_se) and
`marginal_power` (matrix of per-simulation power values for uncertainty).

---

## Scenario B — Single-Cohort Rhythmicity: Simulation-Based Power

**When to use:** You have real pilot data from a single group and want
transcriptome-wide power estimates that account for the full empirical
distribution of A, σ, and phase — not just a summary statistic.

`runSimsSingleCohort()` simulates gene-by-gene from the pilot distributions,
applies the cosinor F-test + BH FDR correction, and reports power across N.

```r
# --- Simulation-based single-cohort power from real pilot ---

# Load a pilot dataset (genes x samples matrix)
# pilot_data <- readRDS("data/my_pilot.rds")  # your data here
# pilot_times <- c(0,2,4,6,8,10,12,14,16,18,20,22)  # ZT times

# For this example, use synthetic pilot
set.seed(42L)
pilot_data  <- matrix(rnorm(500 * 24, 5, 1), nrow = 500)
pilot_times <- rep(seq(0, 22, by = 4), each = 4)

# Estimate parameters from the pilot
bio_pilot <- estCircadianParam(
  data            = pilot_data,
  times           = pilot_times,
  period          = 24,
  prop_DR         = 0,   # single cohort — no differential
  prop_DP         = 0,
  prop_DM         = 0,
  verbose         = TRUE
)

design_single <- CircadianDesignOptions(
  sample_sizes = c(12L, 24L, 48L, 96L),
  nsims        = 30L,
  design       = "active",
  cts          = pilot_times   # reuse pilot time grid
)

res_single <- runSimsSingleCohort(bio_pilot, design_single, analysis_opts, verbose = TRUE)

# Power and FDR at each N
mean_power <- rowMeans(res_single$marginal_power, na.rm = TRUE)
mean_fdr   <- rowMeans(res_single$marginal_FDR,   na.rm = TRUE)

for (i in seq_along(design_single$sample_sizes)) {
  cat(sprintf("N = %3d  power = %.1f%%  FDR = %.1f%%\n",
              design_single$sample_sizes[i],
              100 * mean_power[i],
              100 * mean_fdr[i]))
}
```

**Output:** `list(marginal_power, marginal_FDR, marginal_TD, marginal_FD, sample_sizes, nsims)`
where each matrix is [n_sizes × nsims].

---

## Scenario C — Two-Group Differential Power (DR / DP / DM)

**When to use:** You have two groups (tissue A vs tissue B, disease vs control)
and want to detect genes that differ in rhythmicity (DR), peak timing (DP), or
baseline expression (DM) between groups.

### Pilot estimation: one pilot vs two pilots

All differential endpoints use the same two estimation paths — the choice
depends on what pilot data you have:

**One shared pilot** (one dataset that spans both groups, or a single-group
dataset used as a proxy for both): use `estCircadianParam()`. It estimates
the base parameter distribution (A, σ, φ, mesor) from that single dataset
and you specify the differential proportions (prop_DR, prop_DP, prop_DM) by
hand based on prior knowledge.

**Two separate pilots** (one pilot per group, e.g., pilot tissue A + pilot
tissue B): use `estCircadianParamTwoGroup()`. It fits cosinor models
independently for each group and derives differential simulation inputs
directly from the empirical between-group differences — `prop_DR`, `prop_DP`,
`phase_diff`, and the group-2 mesor distribution (`lBaselineExpr2`) are all
estimated from data rather than guessed. This applies equally to DR, DP,
and DM scenarios.

```r
# --- Two-group differential power from two separate pilots ---

set.seed(7L)
n_pilot     <- 12L
times_pilot <- seq(0, 22, by = 2)   # 12 ZT points, n=1 each
data_g1     <- matrix(rnorm(500 * n_pilot, 5.0, 1), nrow = 500)
data_g2     <- matrix(rnorm(500 * n_pilot, 5.3, 1), nrow = 500)  # shifted mesor

# estCircadianParamTwoGroup estimates prop_DR, prop_DP, phase_diff, amp_diff,
# and lBaselineExpr2 (group-2 mesor distribution) from the two pilots.
# This is the recommended path when pilot data exists for both groups.
bio_twog <- estCircadianParamTwoGroup(
  data_1  = data_g1, data_2 = data_g2,
  times_1 = times_pilot, times_2 = times_pilot,
  period  = 24,
  verbose = TRUE
)

# Adjust differential proportions if prior knowledge warrants it.
# prop_DM can also be set here — it applies to the same two-group simulation.
bio_twog <- updateBioOptions(bio_twog,
  prop_DR = 0.20,   # 20% of genes differentially rhythmic
  prop_DP = 0.05,   # 5% phase-shifted
  prop_DM = 0.08    # 8% have shifted baseline (DM); lBaselineExpr2 from pilot
)

# Design: sweep N per group, test DR + DP + DM simultaneously
design_dr <- CircadianDesignOptions(
  sample_sizes = c(12L, 24L, 36L, 48L, 60L),
  nsims        = 30L,
  design       = "active",
  cts          = rep(seq(0, 22, by = 2), 2),   # 24 samples, 12 ZTs
  test_types   = c("DR", "DP", "DM")
)

analysis_opts <- CircadianAnalysisOptions(alpha = 0.05, p.adjust.method = "BH")
sims_dr <- runSimsDiff(bio_twog, design_dr, analysis_opts)

# Power at each N — helper to avoid repeating the loop
power_by_N <- function(fdr_arr, target_fn, nsims) {
  sapply(seq_len(dim(fdr_arr)[2]), function(j) {
    mean(sapply(seq_len(nsims), function(s) {
      fdr_vec <- fdr_arr[, j, s]
      idx <- target_fn(sims_dr$diff_type[[s]])
      if (sum(idx) == 0) return(NA_real_)
      sum(fdr_vec[idx] <= 0.05, na.rm = TRUE) / sum(idx)
    }), na.rm = TRUE)
  })
}

dr_power <- power_by_N(sims_dr$fdr_DR, function(dt) dt %in% c(2, 3), design_dr$nsims)
dp_power <- power_by_N(sims_dr$fdr_DP, function(dt) dt == 4,          design_dr$nsims)
dm_power <- power_by_N(sims_dr$fdr_DM, function(dt) dt == 5,          design_dr$nsims)

cat("Power by N:\n")
for (i in seq_along(design_dr$sample_sizes)) {
  cat(sprintf("  N = %3d  DR=%.1f%%  DP=%.1f%%  DM=%.1f%%\n",
              design_dr$sample_sizes[i],
              100 * dr_power[i], 100 * dp_power[i], 100 * dm_power[i]))
}
```

### When only one pilot is available

If you only have data from one group (or want to specify priors by hand),
use `estCircadianParam()` and set the differential proportions manually:

```r
bio_onepilot <- estCircadianParam(
  data   = data_g1, times = times_pilot,
  period = 24,
  prop_DR = 0.20, prop_DP = 0.05, prop_DM = 0.08,
  mesor_diff = c(0.5, 2.0),   # assumed DM shift range
  verbose = TRUE
)
# lBaselineExpr2 defaults to NULL → group-2 mesor drawn from same
# distribution as group-1 mesor, plus the DM shift for DM genes.
```

---

## Scenario D — Differential Mesor (DM) as a Primary Endpoint

**When to use:** You want to focus specifically on DM power — e.g., you expect
an inflammatory or metabolic condition to shift baseline expression of rhythmic
genes without disrupting their timing or amplitude. DM genes are rhythmic in
both groups at the same A and φ, but at different mesor levels.

This scenario shows DM-only power with synthetic parameters. For real-pilot
DM estimation alongside DR/DP, use `estCircadianParamTwoGroup()` as shown in
Scenario C.

```r
# --- DM-focused power sweep: synthetic parameters ---

set.seed(42L)

bio_dm <- CircadianBioOptions(
  ngenes        = 1000L,
  prop_rhythmic = 0.35,
  lBaselineExpr = rnorm(1000, 5, 1),
  lOD           = rnorm(1000, -0.5, 0.3),
  amplitude     = abs(rnorm(350, 0.5, 0.2)) + 0.05,
  sigma_rhythmic = NULL,
  phase         = "uniform",
  prop_DR       = 0.10,
  prop_DP       = 0.00,
  prop_DM       = 0.15,   # 15% of genes: both rhythmic, mesor shifted
  mesor_diff    = c(0.8, 2.5),   # shift magnitude ~ Uniform[0.8, 2.5]
  period        = 24L,
  sim.seed      = 42L
)

# Active design: 6 ZT points, 8 replicates → N=48 per group
cts_active <- rep(seq(0, 20, by = 4), each = 8)

design_dm <- CircadianDesignOptions(
  sample_sizes = c(24L, 48L, 72L),
  nsims        = 20L,
  design       = "active",
  cts          = cts_active,
  test_types   = c("DR", "DM")
)

analysis_opts <- CircadianAnalysisOptions(alpha = 0.05, p.adjust.method = "BH")
sims_dm <- runSimsDiff(bio_dm, design_dm, analysis_opts)

# DM power at each N
dm_power <- sapply(seq_along(design_dm$sample_sizes), function(j) {
  mean(sapply(seq_len(design_dm$nsims), function(sim_i) {
    fdr_vec <- sims_dm$fdr_DM[, j, sim_i]
    is_dm   <- sims_dm$diff_type[[sim_i]] == 5
    if (sum(is_dm) == 0) return(NA_real_)
    sum(fdr_vec[is_dm] <= 0.05, na.rm = TRUE) / sum(is_dm)
  }), na.rm = TRUE)
})

cat("DM power by N:\n")
for (i in seq_along(design_dm$sample_sizes)) {
  cat(sprintf("  N = %3d: %.1f%%\n", design_dm$sample_sizes[i], 100 * dm_power[i]))
}
```

**Note on `mesor_diff`:** Each DM gene receives a shift drawn from
`Uniform[mesor_diff[1], mesor_diff[2]]` with a random sign (up or down).
The shift is applied on top of whatever `lBaselineExpr2` provides as the
group-2 background mesor — so when `estCircadianParamTwoGroup()` is used,
the simulation already captures any population-level mesor difference between
groups, and `mesor_diff` specifies the *additional* per-gene DM effect.

---

## Scenario E — B vs m Tradeoff (Bootstrap Design Grid)

**When to use:** You have a real pilot dataset and want to answer:
*given a fixed total budget N, is it better to collect more time points (larger B)
or more replicates per time point (larger m)?*

`runBootstrapDesignGrid()` sweeps (N, B) combinations and propagates pilot
uncertainty via bootstrap resampling of the gene-parameter distribution.

```r
# --- B vs m tradeoff with bootstrap uncertainty ---

# Using real pilot data (e.g., baboon lung pilot, n=12 samples)
# pilot_data  <- baboon_mat   # genes x 12
# pilot_times <- baboon_tod   # length 12

# Synthetic stand-in:
set.seed(5L)
pilot_data  <- matrix(rnorm(800 * 12, 5, 1), nrow = 800)
pilot_times <- seq(0, 22, by = 2)   # 12 ZT points

bio_diff <- CircadianBioOptions(
  ngenes        = 800L,
  prop_rhythmic = 0.35,
  lBaselineExpr = rnorm(800, 5, 1),
  lOD           = rnorm(800, -0.5, 0.3),
  amplitude     = abs(rnorm(280, 0.6, 0.25)) + 0.05,
  phase         = "uniform",
  prop_DR       = 0.25,
  prop_DP       = 0.05,
  prop_DM       = 0.00,
  period        = 24L,
  sim.seed      = 5L
)

boot_opts <- CircadianBootstrapOptions(
  N_values      = c(12L, 24L, 36L, 48L),
  B_values      = c(4L, 6L, 8L, 12L),      # candidate # of time points
  nboot         = 20L,                       # outer bootstrap draws
  nsims_inner   = 10L,                       # inner simulations per bootstrap
  design        = "active",
  design_vector = seq(0, 22, by = 2),        # full possible time grid
  seed          = 42L
)

boot_res <- runBootstrapDesignGrid(
  pilot_data    = pilot_data,
  pilot_times   = pilot_times,
  boot.opts     = boot_opts,
  analysis.opts = analysis_opts,
  bio_diff.opts = bio_diff,
  verbose       = TRUE,
  mc.cores      = 1L        # increase to use parallel cores (see Parallelization below)
)

# Summarize: power ± 95% CI for DR test
summaryBootstrapDesignGrid(boot_res, test_type = "DR")

# Plot power curves (Panel A) and heatmap (Panel B)
plotBootstrapDesignGrid(boot_res, test_type = "DR", panels = c("A", "B"))
```

**Reading the output:**
- `boot_res$power_mean[n_idx, B_idx, test_idx]` — mean power over bootstrap draws
- `boot_res$power_ci_lo / ci_hi` — 95% bootstrap CI on power
- `boot_res$optimal_B[n_idx]` — B that maximises power at each N

---

## Scenario F — Two-Stage vs Bootstrap Comparison

**When to use:** You want to compare the point-estimate (two-stage) approach to
the bootstrap approach, and quantify how much uncertainty the pilot size adds.

Two-stage uses pilot parameters once (no CI); bootstrap resamples to produce
a distribution of power curves. When the pilot is small (n<20), bootstrap CIs
are wide — this is *honest* about uncertainty rather than falsely precise.

```r
# --- Side-by-side comparison ---

# Two-stage (point-estimate, no CI)
two_stage_res <- runTwoStagePower(
  pilot_data    = pilot_data,
  pilot_times   = pilot_times,
  design.opts   = CircadianDesignOptions(
    sample_sizes = c(12L, 24L, 36L, 48L),
    nsims        = 20L,
    design       = "active",
    cts          = rep(seq(0, 22, by = 2), 2),
    test_types   = "DR"
  ),
  analysis.opts = analysis_opts,
  bio_diff.opts = bio_diff,
  test_type     = "DR",
  verbose       = TRUE,
  mc.cores      = 1L   # parallelizes over sample_sizes
)

# Bootstrap (single fixed B for apples-to-apples CI comparison)
boot_res_fixed <- runBootstrapDesignGrid(
  pilot_data    = pilot_data,
  pilot_times   = pilot_times,
  boot.opts     = CircadianBootstrapOptions(
    N_values      = c(12L, 24L, 36L, 48L),
    B_values      = 12L,                    # fix B=12 for this comparison
    nboot         = 30L,
    nsims_inner   = 10L,
    design        = "active",
    design_vector = seq(0, 22, by = 2),
    seed          = 42L
  ),
  analysis.opts = analysis_opts,
  bio_diff.opts = bio_diff,
  verbose       = TRUE
)

# Compare: pull DR mean ± CI from bootstrap
t_idx <- 1  # DR is first test
for (n_idx in seq_along(c(12, 24, 36, 48))) {
  N      <- c(12, 24, 36, 48)[n_idx]
  ts_pow <- two_stage_res$power_mean[n_idx]
  bt_pow <- boot_res_fixed$power_mean[n_idx, 1, t_idx]
  bt_lo  <- boot_res_fixed$power_ci_lo[n_idx, 1, t_idx]
  bt_hi  <- boot_res_fixed$power_ci_hi[n_idx, 1, t_idx]
  cat(sprintf("N=%2d  two-stage=%.1f%%  bootstrap=%.1f%% [%.1f%%, %.1f%%]\n",
              N, 100*ts_pow, 100*bt_pow, 100*bt_lo, 100*bt_hi))
}

# compareDesignApproaches() merges both into a comparison data frame
cmp <- compareDesignApproaches(two_stage_res, boot_res_fixed, test_type = "DR")
```

**Interpretation guide:**
- If CI is narrow: pilot is large enough that two-stage and bootstrap agree — either is reliable.
- If CI is wide: pilot is small (e.g., baboon n=12) — two-stage's single number is misleadingly precise; trust the bootstrap CI for study planning.
- Practical rule of thumb: if the 95% CI width on n80 exceeds ±50%, consider collecting more pilot samples before committing to a study design.

---

## Parallelization

Two functions support parallel execution via `mc.cores` (uses `parallel::mclapply`
internally — works on Linux/macOS; on Windows use `mc.cores = 1L`):

| Function | What is parallelized | Typical speedup |
|----------|---------------------|-----------------|
| `runBootstrapDesignGrid()` | Outer bootstrap draws (each `nboot` draw is independent) | ~linear up to `nboot` cores |
| `runTwoStagePower()` | Per-`sample_sizes` iterations | ~linear up to `length(sample_sizes)` cores |

```r
# Detect available cores (leave one free for the OS)
n_cores <- max(1L, parallel::detectCores() - 1L)

# Bootstrap grid — parallelized over nboot=50 draws
boot_res <- runBootstrapDesignGrid(
  pilot_data    = pilot_data,
  pilot_times   = pilot_times,
  boot.opts     = boot_opts,
  analysis.opts = analysis_opts,
  bio_diff.opts = bio_diff,
  mc.cores      = n_cores
)

# Two-stage — parallelized over sample_sizes
ts_res <- runTwoStagePower(
  pilot_data    = pilot_data,
  pilot_times   = pilot_times,
  design.opts   = design_opts,
  analysis.opts = analysis_opts,
  bio_diff.opts = bio_diff,
  test_type     = "DR",
  mc.cores      = n_cores
)
```

**Practical guidance:**
- `nboot = 50`, `mc.cores = 8` → bootstrap grid runs in ~1/8 the time
- `nsims_inner` is NOT parallelized — keep it moderate (10–30) and increase `nboot` for precision
- On a cluster with SLURM, set `mc.cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))`
- Seeds: each bootstrap draw gets its own independent seed derived from `boot.opts$seed`; results are reproducible regardless of `mc.cores`

---

## Gene Type Reference

| Type | Label | Rhythmic G1 | Rhythmic G2 | Differential feature |
|------|-------|-------------|-------------|----------------------|
| 0    | Arrhythmic | No | No | — |
| 1    | Rhythmic (same) | Yes | Yes | — |
| 2    | DR (G1 only) | Yes | No | Loss of rhythmicity in G2 |
| 3    | DR (G2 only) | No | Yes | Gain of rhythmicity in G2 |
| 4    | DP (phase shift) | Yes | Yes | Phase differs by ≥ phase_diff |
| 5    | DM (mesor shift) | Yes | Yes | Baseline differs by mesor_diff |
| 6    | DA (amplitude diff) | Yes | Yes | Amplitude ratio ≥ amp_diff |

**Proportion budget constraint:**
```
prop_DR + prop_DP + prop_DM + prop_DA  <=  prop_rhythmic
```
All differential genes must be rhythmic in at least one group.

---

## Key Function Quick Reference

| Task | Function | File |
|------|----------|------|
| Build simulation parameters | `CircadianBioOptions()` | `code/options.R` |
| Update simulation parameters | `updateBioOptions()` | `code/options.R` |
| Build design options | `CircadianDesignOptions()` | `code/options.R` |
| Estimate from single-group pilot | `estCircadianParam()` | `code/estimation.R` |
| Estimate from two-group pilot | `estCircadianParamTwoGroup()` | `code/estimation.R` |
| Closed-form single-cohort power | `runPowerAnalysis()` | `code/runner.R` |
| Simulation-based single-cohort power | `runSimsSingleCohort()` | `code/runner.R` |
| Two-group differential simulation | `runSimsDiff()` | `code/runner.R` |
| Bootstrap design grid (N × B) | `runBootstrapDesignGrid()` | `code/bootstrap_sim.R` |
| Two-stage point-estimate power | `runTwoStagePower()` | `code/design_comparison.R` |
| Compare two-stage vs bootstrap | `compareDesignApproaches()` | `code/design_comparison.R` |
| Summarize bootstrap grid | `summaryBootstrapDesignGrid()` | `code/bootstrap_sim.R` |
| Plot bootstrap grid | `plotBootstrapDesignGrid()` | `code/bootstrap_sim.R` |
