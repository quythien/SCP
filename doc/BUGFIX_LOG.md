# SCP Framework — Bug Fix Log

All confirmed bugs found during development and code review, with root cause,
affected files, and the fix applied. Ordered from most recent session backward.

---

## Session 2026-04-20: RAIN replicated-design bug

### BUG-RAIN-01: `detection.R` — `detect_RAIN()` broken for replicated designs

**Root cause:** Two errors in the original `detect_RAIN()` wrapper:

1. `deltat <- median(diff(times_sorted))` computes the median gap across all
   consecutive sorted times. With m replicates per time point, consecutive
   entries are at the same ZT → gap = 0 → `deltat = 0`. RAIN then errors or
   produces degenerate results.

2. `nr.series` was never passed to `rain::rain()`. RAIN's umbrella test assumes
   a single time series when `nr.series=1`; without it, each replicate is treated
   as an independent time point rather than a replicate within a series.

**Effect:** For any replicated active design (B < N), RAIN either errored or
treated the data as if m=1 with B×m distinct time points. Power estimates were
wrong for all B < N cells.

**Fix (`code/detection.R` line ~1730):**
```r
# Before
deltat <- median(diff(times_sorted))
rain::rain(x = t(expr_sorted), deltat = deltat, period = period)

# After
unique_times <- sort(unique(times_sorted))
B_pts  <- length(unique_times)
deltat <- if (B_pts > 1) median(diff(unique_times)) else period
counts <- tabulate(match(times_sorted, unique_times))
nr     <- if (length(unique(counts)) == 1L && counts[1L] > 1L) counts[1L] else 1L
rain::rain(x = t(expr_sorted), deltat = deltat, period = period, nr.series = nr)
```

**Affected analyses:** Any `detect_RAIN()` call on a replicated design.
`15b_bvsm_rain.R` is the first production script using the corrected version.
Prior exploratory RAIN results (D1D2 smoke tests) were also wrong; the
corrected results show RAIN favors higher B (not lower B as was initially
suggested by the broken wrapper).

---

## Session 2026-04-13 (Round 2): Post-DM type-numbering bugs

After adding DM (Differential Mesor) as gene type 5 and reassigning DA
(Differential Amplitude) to type 6, several dispatch blocks and helper
functions still used the old numbering.

---

### BUG-01: `bootstrap_sim.R` — DA dispatch targeted DM genes

**Root cause:** When DM was added as type 5, DA was moved to type 6. The inner
dispatch block in `runBootstrapDesignGrid()` still used `diff_type_vec == 5`
for "DA", so requesting a "DA" test measured power on DM genes instead.

**Affected file:** `code/bootstrap_sim.R` (inner `powers <- sapply(...)` block)

**Before:**
```r
} else if (tt == "DA") {
  diff_type_vec == 5
```

**After:**
```r
} else if (tt == "DM") {
  diff_type_vec == 5
} else if (tt == "DA") {
  diff_type_vec == 6
```

---

### BUG-02: `bootstrap_sim.R` — DM never included in `all_tests`

**Root cause:** `runBootstrapDesignGrid()` builds `all_tests` by checking
`prop_DR`, `prop_DP`, `prop_DA` from `bio_diff.opts` but never checked
`prop_DM`. Users setting `prop_DM > 0` got no DM power estimates from the
bootstrap grid — silently dropped.

**Affected file:** `code/bootstrap_sim.R` (`all_tests` construction block)

**Before:**
```r
if (bio_diff.opts$prop_DA > 0) all_tests <- c(all_tests, "DA")
```

**After:**
```r
if (!is.null(bio_diff.opts$prop_DM) && bio_diff.opts$prop_DM > 0) all_tests <- c(all_tests, "DM")
if (!is.null(bio_diff.opts$prop_DA) && bio_diff.opts$prop_DA > 0) all_tests <- c(all_tests, "DA")
```

---

### BUG-03: `bootstrap_sim.R` `.buildBioFromBoot()` — DM params dropped on bootstrap resample

**Root cause:** `.buildBioFromBoot()` builds a `CircadianBioOptions` from each
bootstrap draw. It copied `prop_DR`, `prop_DP`, `prop_DA` from `bio_diff.opts`
but not `prop_DM`, `mesor_diff`, or `lBaselineExpr2`. Every bootstrap draw
therefore silently zeroed out DM genes, making BUG-02 moot even when fixed.

**Affected file:** `code/bootstrap_sim.R` (`.buildBioFromBoot()`, `opts <- list(...)`)

**After:** Added to the list:
```r
prop_DM        = bio_diff.opts$prop_DM   %||% 0,
mesor_diff     = bio_diff.opts$mesor_diff %||% c(0.5, 2.0),
lBaselineExpr2 = bio_diff.opts$lBaselineExpr2,
```

---

### BUG-04: `design_comparison.R` `runTwoStagePower()` — DA dispatch targeted DM genes

**Root cause:** Same type-numbering error as BUG-01, in a separate file.
`runTwoStagePower()` used `diff_type == 5` for its "DA" branch.

**Affected file:** `code/design_comparison.R` (inner sim loop)

**Before:**
```r
} else if (test_type == "DA") {
  target_idx <- diff_type == 5
```

**After:**
```r
} else if (test_type == "DM") {
  target_idx <- diff_type == 5
} else if (test_type == "DA") {
  target_idx <- diff_type == 6
```

---

### BUG-05: `design_comparison.R` `runTwoStagePower()` — `prop_DM` not passed to `estCircadianParam()`

**Root cause:** `runTwoStagePower()` calls `estCircadianParam()` to estimate
pilot parameters, but did not pass `prop_DM` or `mesor_diff`. The two-stage
approach could not run a DM power analysis.

**Affected file:** `code/design_comparison.R` (`runTwoStagePower()`, Step 1)

**After:** Added:
```r
prop_DM    = bio_diff.opts$prop_DM   %||% 0,
mesor_diff = bio_diff.opts$mesor_diff %||% c(0.5, 2.0),
```

---

### BUG-06: `options.R` `updateBioOptions()` — `prop_DM`, `mesor_diff`, `lBaselineExpr2` silently dropped

**Root cause:** `updateBioOptions()` reconstructs a `CircadianBioOptions` by
pulling a fixed list of fields from the existing opts object, then merging user
overrides. `prop_DM`, `mesor_diff`, and `lBaselineExpr2` were not in that
fixed list, so `updateBioOptions(opts, prop_DM = 0.1)` would silently lose all
three fields in the returned object.

**Affected file:** `code/options.R` (`updateBioOptions()`)

**Before:**
```r
current_args <- opts[c("ngenes", "prop_rhythmic", "period",
                       "prop_DR", "prop_DP", "prop_DA",
                       "phase_diff", ...)]
current_args$lBaselineExpr <- opts$lBaselineExpr_spec
```

**After:**
```r
current_args <- opts[c("ngenes", "prop_rhythmic", "period",
                       "prop_DR", "prop_DP", "prop_DA", "prop_DM",
                       "phase_diff", ..., "mesor_diff", "sim.seed")]
current_args$lBaselineExpr  <- opts$lBaselineExpr_spec
current_args$lBaselineExpr2 <- opts$lBaselineExpr2
```

---

### BUG-07: `options.R` `print.CircadianBioOptions()` — `prop_DM` not displayed

**Root cause:** The S3 print method listed `prop_DR`, `prop_DP`, `prop_DA` but
not `prop_DM`. Interactive inspection of options with DM set would not show the
value, making debugging harder.

**Affected file:** `code/options.R` (`print.CircadianBioOptions()`)

**After:** Added one line:
```r
cat(sprintf("  prop_DM:        %.0f%%\n", 100 * (x$prop_DM %||% 0)))
```

---

### BUG-08: `estimation.R` `estCircadianParam()` — no `prop_DM` in signature or budget check

**Root cause:** `estCircadianParam()` accepted `prop_DR`, `prop_DP`, `prop_DA`
but not `prop_DM`. The proportion budget check (`total_diff > prop_rhythmic`)
only summed the three old endpoints, so a caller specifying DM would have it
silently ignored.

**Affected file:** `code/estimation.R` (`estCircadianParam()`)

**After:**
- Added `prop_DM = 0.00` and `mesor_diff = c(0.5, 2.0)` to signature
- `total_diff <- prop_DR + prop_DP + prop_DA + prop_DM`
- `prop_DM` scaled with the others if budget exceeded
- Passed through to `CircadianBioOptions()`

---

## Session 2026-04-13 (Round 1): DM / single-cohort implementation bugs

---

### BUG-09: `simulation.R` `simCircadianDiff()` — hard stop on `cts` length mismatch

**Root cause:** The function required `length(cts) == n1` exactly. When called
from `runSimsDiff()` with a compact ZT grid (e.g., 12 elements) and a target
`n1 = 60`, it threw a hard error. The correct behavior is to tile the ZT grid
to match the sample count.

**Affected file:** `code/simulation.R`

**Before:**
```r
stop("Active design: length(cts) must equal n1")
```

**After:**
```r
cts <- sort(rep_len(cts, n1))   # tile grid to fill n1 samples
```
Same fix applied to `cts2`.

---

### BUG-10: `runner.R` — `cts2` length mismatch when sample size differs from pilot

**Root cause:** `runSimsDiff()` stored `bio.opts$cts2` from the pilot (e.g.,
24 elements) and passed it directly to `simCircadianDiff()` at each candidate
sample size (e.g., N=48 or N=80). Since `simCircadianDiff()` then checked for
`length(cts2) == n2`, any N ≠ pilot size crashed.

**Affected file:** `code/runner.R` (`runSimsDiff()`)

**Fix:** Added expansion logic before calling sim_args:
```r
cts2_n <- if (design == "active" && length(cts2_raw) != n) {
  sort(rep_len(cts2_raw, n))
} else { cts2_raw }
```

---

### BUG-11: `09_dm_singlecohort_smoke.R` — `sys.frame(1)` crash when run via Rscript

**Root cause:** `sys.frame(1)$ofile` is used to locate the script's directory
for `POWERSIM_ROOT` resolution. When the script is invoked via `Rscript -e`
(inline evaluation context), there is no frame 1 on the stack, causing:
`Error in sys.frame(1): not that many frames on the stack`.

**Affected file:** `examples/publication/09_dm_singlecohort_smoke.R`

**Fix:** Wrapped in `tryCatch` with fallback to `getwd()`:
```r
POWERSIM_ROOT <- tryCatch({
  script_path <- sys.frame(1)$ofile
  normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
}, error = function(e) {
  if (file.exists(file.path(getwd(), "code", "setup.R"))) return(getwd())
  ...
})
```

---

## Cross-cutting statistical notes

### Gene type numbering (definitive reference)

| Type | Label | Meaning |
|------|-------|---------|
| 0 | Arrhythmic | Not rhythmic in either group |
| 1 | Rhythmic (same) | Rhythmic in both; no differential feature |
| 2 | DR (G1 only) | Rhythmic in G1, arrhythmic in G2 |
| 3 | DR (G2 only) | Arrhythmic in G1, rhythmic in G2 |
| 4 | DP | Phase shift ≥ phase_diff between groups |
| 5 | DM | Differential mesor; both rhythmic, same A/φ |
| 6 | DA | Differential amplitude; both rhythmic, same φ |

DA was originally type 5 before DM was added. BUG-01 through BUG-05 all stem
from this renumbering not being propagated to all dispatch blocks.

### Proportion budget constraint

All differential types (DR, DP, DM, DA) draw from the rhythmic gene pool:

```
prop_DR + prop_DP + prop_DM + prop_DA  <=  prop_rhythmic
```

`estCircadianParam()` now scales all four props together when the budget is
exceeded (BUG-08 fix).
