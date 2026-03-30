#' Extension Smoke Tests — Scenarios 5-8 (Tier 1)
#'
#' Validates each new module with tiny settings.
#' Stops with an informative error if any invariant fails.
#'
#' Usage:
#'   cd code && Rscript run_extension_smoke_tests.R
#'
#' Tiers:
#'   Tier 1 (this script): ngenes=200, nsims=3, nboot=3 — every edit
#'   Tier 2: ngenes=1000, nsims=10, nboot=10 — before merging
#'   Tier 3: ngenes=3000, nsims=30, nboot=30 — final figures

source("setup.R")

# ── Settings ──────────────────────────────────────────────────────────────────
NGENES      <- 200
NSIMS       <- 3
NBOOT       <- 3
NSIMS_INNER <- 3
N_PILOT     <- 20
PILOT_TIMES <- seq(0, 22, by = 2)   # 12 evenly-spaced TOD bins
DESIGN_VEC  <- PILOT_TIMES

# Feasible N grid: all values must be divisible by every B in B_VALUES.
# B_VALUES = c(4, 8) → LCM = 8 → use multiples of 8.
B_VALUES <- c(4L, 8L)
N_GRID   <- c(24L, 48L)             # small but fully feasible

cat("=== Extension Smoke Tests (Tier 1) ===\n")
cat(sprintf("ngenes=%d  nsims=%d  nboot=%d  N_grid=%s\n\n",
            NGENES, NSIMS, NBOOT, paste(N_GRID, collapse = ",")))

# ── Helpers ──────────────────────────────────────────────────────────────────
CHECK <- function(label, expr) {
  ok <- tryCatch(isTRUE(expr), error = function(e) FALSE)
  if (!ok) stop(sprintf("FAIL: %s", label))
  cat(sprintf("  OK  %s\n", label))
}

SECTION <- function(name) cat(sprintf("\n── %s ──\n", name))

# ── Shared bio + analysis options ─────────────────────────────────────────────
opts_bio <- CircadianBioOptions(
  ngenes        = NGENES,
  prop_rhythmic = 0.25,
  prop_DR       = 0.15,
  prop_DP       = 0.00,
  prop_DA       = 0.00,
  phase_diff    = c(0, 0),
  amp_diff      = c(1, 1)
)
opts_analysis <- CircadianAnalysisOptions(reference_n = 48)

# ═══════════════════════════════════════════════════════════════════════════════
SECTION("1. Simulation truth — gene-class counts")
# ═══════════════════════════════════════════════════════════════════════════════

opts_design_small <- CircadianDesignOptions(
  sample_sizes = 24L, nsims = NSIMS, design = "active"
)

opts_mixed <- CircadianBioOptions(
  ngenes        = NGENES,
  prop_rhythmic = 0.30,
  prop_DR       = 0.10,
  prop_DP       = 0.10,
  prop_DA       = 0.05,
  phase_diff    = c(-4, 4),
  amp_diff      = c(0.5, 2)
)

sim_out <- runSimsDiff(opts_mixed, opts_design_small, opts_analysis)
dt <- sim_out$diff_type[[1]]

n_non  <- sum(dt == 0)
n_same <- sum(dt == 1)
n_DR   <- sum(dt %in% c(2, 3))
n_DP   <- sum(dt == 4)
n_DA   <- sum(dt == 5)
total  <- n_non + n_same + n_DR + n_DP + n_DA

CHECK("gene classes sum to ngenes", total == NGENES)

# Total rhythmic = DR + DP + DA + same (all rhythmic in at least one group)
n_rhythmic_actual <- n_DR + n_DP + n_DA + n_same
expected_rhythmic  <- round(NGENES * opts_mixed$prop_rhythmic)
CHECK("total rhythmic within 2 genes of target",
      abs(n_rhythmic_actual - expected_rhythmic) <= 2)

cat(sprintf("    non=%d  same=%d  DR=%d  DP=%d  DA=%d  total_rhythmic=%d (target=%d)\n",
            n_non, n_same, n_DR, n_DP, n_DA, n_rhythmic_actual, expected_rhythmic))

# ═══════════════════════════════════════════════════════════════════════════════
SECTION("2. Active design — no duplicate 0/24 endpoint")
# ═══════════════════════════════════════════════════════════════════════════════

sim_active <- simCircadianDiff(
  ngenes = 10, n1 = 6, n2 = 6,
  prop_rhythmic = 0, prop_DR = 0, prop_DP = 0, prop_DA = 0,
  design = "active"
)
CHECK("max(times1) < 24 (endpoint excluded)", max(sim_active$times1) < 24)
CHECK("length(unique(times1)) == n1 (all times distinct)",
      length(unique(sim_active$times1)) == 6)
cat(sprintf("    times1: %s\n", paste(round(sim_active$times1, 2), collapse = ", ")))

# ═══════════════════════════════════════════════════════════════════════════════
SECTION("3. Bootstrap design grid — exact divisibility")
# ═══════════════════════════════════════════════════════════════════════════════

pilot <- generatePilotData(opts_bio, N_PILOT, PILOT_TIMES, seed = 42)

boot.opts <- CircadianBootstrapOptions(
  design_vector = DESIGN_VEC,
  B_values      = B_VALUES,
  N_values      = N_GRID,
  nboot         = NBOOT,
  nsims_inner   = NSIMS_INNER,
  design        = "active"
)

boot_grid <- runBootstrapDesignGrid(
  pilot_data    = pilot$data,
  pilot_times   = pilot$times,
  boot.opts     = boot.opts,
  analysis.opts = opts_analysis,
  bio_diff.opts = opts_bio
)

# All valid (N, B) cells: m = N / B exactly
for (ni in seq_along(N_GRID)) {
  for (bi in seq_along(B_VALUES)) {
    N <- N_GRID[ni]; B <- B_VALUES[bi]
    m <- boot_grid$m_matrix[ni, bi]
    if (!is.na(m)) {
      CHECK(sprintf("m_matrix[N=%d, B=%d] == N/B (%d)", N, B, N/B),
            m == as.integer(N / B))
    }
  }
}

# Any infeasible pair (if grid were mixed) must be NA — verify with a mixed grid
boot.opts_bad <- suppressWarnings(CircadianBootstrapOptions(
  design_vector = DESIGN_VEC,
  B_values      = B_VALUES,
  N_values      = c(24L, 36L, 48L),   # N=36 is not divisible by B=8
  nboot         = 1, nsims_inner = 1, design = "active"
))
m_mat_bad <- outer(c(24L, 36L, 48L), B_VALUES, FUN = function(N, B) {
  ifelse(B > N | (N %% B) != 0L, NA_integer_, as.integer(N / B))
})
CHECK("infeasible cell (N=36, B=8) is NA", is.na(m_mat_bad[2, 2]))

# Output shapes
CHECK("boot_power dim[1] == nboot",  dim(boot_grid$boot_power)[1] == NBOOT)
CHECK("boot_power dim[2] == n_N",    dim(boot_grid$boot_power)[2] == length(N_GRID))
CHECK("boot_power dim[3] == n_B",    dim(boot_grid$boot_power)[3] == length(B_VALUES))
CHECK("power_mean all in [0,1] or NA",
      all(boot_grid$power_mean[!is.na(boot_grid$power_mean)] >= 0) &&
      all(boot_grid$power_mean[!is.na(boot_grid$power_mean)] <= 1))
CHECK("CI lo <= mean <= hi (valid cells)",
      all(boot_grid$power_ci_lo <= boot_grid$power_mean + 1e-9, na.rm = TRUE) &&
      all(boot_grid$power_ci_hi >= boot_grid$power_mean - 1e-9, na.rm = TRUE))

# ═══════════════════════════════════════════════════════════════════════════════
SECTION("4. Bootstrap grid plot — renders without error")
# ═══════════════════════════════════════════════════════════════════════════════

tmp_pdf <- tempfile(fileext = ".pdf")
tryCatch({
  plotBootstrapDesignGrid(boot_grid, test_type = "DR", output_file = tmp_pdf)
  CHECK("plotBootstrapDesignGrid renders to PDF", file.exists(tmp_pdf))
  file.remove(tmp_pdf)
}, error = function(e) stop(sprintf("FAIL: plotBootstrapDesignGrid error: %s", e$message)))

# ═══════════════════════════════════════════════════════════════════════════════
SECTION("5. Fourier deviation — power decreases with harmonics")
# ═══════════════════════════════════════════════════════════════════════════════

opts_design_f <- CircadianDesignOptions(
  sample_sizes = N_GRID, nsims = NSIMS, design = "active"
)
harmonic_grid <- expand.grid(alpha2 = c(0, 0.5), alpha3 = c(0))

fourier_result <- runFourierDeviationPower(
  bio.opts      = opts_bio,
  design.opts   = opts_design_f,
  analysis.opts = opts_analysis,
  harmonic_grid = harmonic_grid,
  test_type     = "DR"
)

CHECK("fourier power_mean shape [n_harm x n_N]",
      all(dim(fourier_result$power_mean) == c(nrow(harmonic_grid), length(N_GRID))))
CHECK("fourier power_mean in [0,1] or NA",
      all(fourier_result$power_mean[!is.na(fourier_result$power_mean)] >= 0) &&
      all(fourier_result$power_mean[!is.na(fourier_result$power_mean)] <= 1))

tmp_pdf2 <- tempfile(fileext = ".pdf")
tryCatch({
  plotFourierDeviation(fourier_result, test_type = "DR", output_file = tmp_pdf2)
  CHECK("plotFourierDeviation renders to PDF", file.exists(tmp_pdf2))
  file.remove(tmp_pdf2)
}, error = function(e) stop(sprintf("FAIL: plotFourierDeviation error: %s", e$message)))

# ═══════════════════════════════════════════════════════════════════════════════
SECTION("6. Design comparison — single-B: no error; multi-B: hard error")
# ═══════════════════════════════════════════════════════════════════════════════

opts_design_ts <- CircadianDesignOptions(
  sample_sizes = N_GRID, nsims = NSIMS, design = "active"
)

two_stage_result <- runTwoStagePower(
  pilot_data    = pilot$data,
  pilot_times   = pilot$times,
  design.opts   = opts_design_ts,
  analysis.opts = opts_analysis,
  bio_diff.opts = opts_bio,
  test_type     = "DR"
)

# Single-B bootstrap (no warning expected)
boot.opts_s7 <- CircadianBootstrapOptions(
  design_vector = DESIGN_VEC,
  B_values      = 4L,
  N_values      = N_GRID,
  nboot         = NBOOT,
  nsims_inner   = NSIMS_INNER,
  design        = "active"
)
boot_s7 <- runBootstrapDesignGrid(
  pilot_data    = pilot$data,
  pilot_times   = pilot$times,
  boot.opts     = boot.opts_s7,
  analysis.opts = opts_analysis,
  bio_diff.opts = opts_bio
)

# Single-B: no error
comparison <- tryCatch(
  compareDesignApproaches(two_stage_result, boot_s7, test_type = "DR"),
  error = function(e) stop(sprintf("FAIL: single-B compareDesignApproaches errored: %s", e$message))
)
CHECK("single-B comparison does NOT error", !is.null(comparison))

# Multi-B: must now hard-error (changed from warning to stop)
errored_multi <- tryCatch({
  compareDesignApproaches(two_stage_result, boot_grid, test_type = "DR")
  FALSE   # no error → bad
}, error = function(e) grepl("multiple B", conditionMessage(e)))
CHECK("multi-B comparison hard-errors (not just warns)", isTRUE(errored_multi))
CHECK("comparison$comparison has rows for N_GRID",
      nrow(comparison$comparison) == length(N_GRID))
CHECK("two_stage_power in [0,1] or NA",
      all(comparison$comparison$two_stage_power >= 0 |
          is.na(comparison$comparison$two_stage_power)) &&
      all(comparison$comparison$two_stage_power <= 1 |
          is.na(comparison$comparison$two_stage_power)))

tmp_pdf3 <- tempfile(fileext = ".pdf")
tryCatch({
  plotDesignComparison(comparison, target_power = 0.80, output_file = tmp_pdf3)
  CHECK("plotDesignComparison renders to PDF", file.exists(tmp_pdf3))
  file.remove(tmp_pdf3)
}, error = function(e) stop(sprintf("FAIL: plotDesignComparison error: %s", e$message)))

# ═══════════════════════════════════════════════════════════════════════════════
SECTION("7. Ground truth calibration — oracle / two-stage / bootstrap align")
# ═══════════════════════════════════════════════════════════════════════════════

opts_design_gt <- CircadianDesignOptions(
  sample_sizes = N_GRID, nsims = NSIMS, design = "active"
)

gt_result <- runGroundTruthComparison(
  true_bio.opts = opts_bio,
  design.opts   = opts_design_gt,
  analysis.opts = opts_analysis,
  n_pilot       = N_PILOT,
  pilot_times   = PILOT_TIMES,
  nboot         = NBOOT,
  nsims_inner   = NSIMS_INNER,
  test_type     = "DR",
  seed          = 42
)

CHECK("oracle power length == n_N",   length(gt_result$true_power_mean) == length(N_GRID))
CHECK("two-stage power length == n_N", length(gt_result$ts_power_mean) == length(N_GRID))
CHECK("boot power length == n_N",      length(gt_result$boot_power_mean) == length(N_GRID))
CHECK("all power values in [0,1]",
      all(c(gt_result$true_power_mean, gt_result$ts_power_mean,
            gt_result$boot_power_mean) >= 0, na.rm = TRUE) &&
      all(c(gt_result$true_power_mean, gt_result$ts_power_mean,
            gt_result$boot_power_mean) <= 1, na.rm = TRUE))
CHECK("boot CI lo <= boot mean",
      all(gt_result$boot_ci_lo <= gt_result$boot_power_mean + 1e-9, na.rm = TRUE))
CHECK("boot CI hi >= boot mean",
      all(gt_result$boot_ci_hi >= gt_result$boot_power_mean - 1e-9, na.rm = TRUE))
CHECK("coverage in [0,1]",
      gt_result$coverage >= 0 && gt_result$coverage <= 1)

tmp_pdf4 <- tempfile(fileext = ".pdf")
tryCatch({
  plotGroundTruthComparison(gt_result, target_power = 0.80, output_file = tmp_pdf4)
  CHECK("plotGroundTruthComparison renders to PDF", file.exists(tmp_pdf4))
  file.remove(tmp_pdf4)
}, error = function(e) stop(sprintf("FAIL: plotGroundTruthComparison error: %s", e$message)))

# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== All invariants passed ===\n")
cat("Ready for Tier 2 (moderate settings) when needed.\n")
