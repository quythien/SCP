## =============================================================================
## test_all_functions.R
## Comprehensive smoke + invariant tests for the SCP / CircadianPower framework.
##
## Run from project root (PowerSim/):
##   Rscript tests/test_all_functions.R
## or inside an R session:
##   source("tests/test_all_functions.R")
## =============================================================================

options(warn = 1)   # surface warnings immediately

# ---------------------------------------------------------------------------
# 0.  Test harness
# ---------------------------------------------------------------------------

.PASS <- 0L
.FAIL <- 0L
.RESULTS <- list()

.record <- function(name, ok, msg = "") {
  status <- if (isTRUE(ok)) "PASS" else "FAIL"
  cat(sprintf("  [%s] %s%s\n", status, name,
              if (nchar(msg) > 0) paste0("  -- ", msg) else ""))
  if (isTRUE(ok)) .PASS <<- .PASS + 1L else .FAIL <<- .FAIL + 1L
  .RESULTS[[length(.RESULTS) + 1L]] <<- list(name = name, status = status, msg = msg)
  invisible(ok)
}

.expect_true <- function(name, expr, ...) {
  ok <- tryCatch(isTRUE(expr), error = function(e) FALSE)
  .record(name, ok, ...)
}

.expect_no_error <- function(name, expr_fn, ...) {
  ok <- tryCatch({ expr_fn(); TRUE }, error = function(e) { .record(name, FALSE, conditionMessage(e)); FALSE })
  if (ok) .record(name, TRUE, ...)
  invisible(ok)
}

.expect_class <- function(name, obj, cls) {
  ok <- inherits(obj, cls)
  .record(name, ok, if (!ok) sprintf("got class '%s', expected '%s'", paste(class(obj), collapse="/"), cls) else "")
}

.expect_range <- function(name, x, lo = 0, hi = 1) {
  ok <- all(!is.na(x) & x >= lo & x <= hi)
  .record(name, ok, if (!ok) sprintf("values outside [%g, %g]: %s", lo, hi,
                                      paste(x[!is.na(x) & (x < lo | x > hi)][seq_len(min(3, sum(!is.na(x) & (x<lo|x>hi))))], collapse=", ")) else "")
}

# ---------------------------------------------------------------------------
# 1.  Load all functions via setup.R
# ---------------------------------------------------------------------------
cat("\n--- Section 1: Loading setup.R ---\n")

old_wd <- setwd("code")
tryCatch({
  source("setup.R")
  setwd(old_wd)
  .record("setup.R loads without error", TRUE)
}, error = function(e) {
  setwd(old_wd)
  .record("setup.R loads without error", FALSE, conditionMessage(e))
  cat("\nFATAL: setup.R failed to load. Aborting tests.\n")
  quit(status = 1)
})

# ---------------------------------------------------------------------------
# 2.  Helper: small synthetic pilot matrix
# ---------------------------------------------------------------------------

set.seed(1)
N_PILOT  <- 48
G_PILOT  <- 200
pilot_times <- seq(0, 24, length.out = N_PILOT + 1L)[seq_len(N_PILOT)]
pilot_data  <- {
  omega <- 2 * pi / 24
  mat   <- matrix(NA_real_, nrow = G_PILOT, ncol = N_PILOT)
  for (g in seq_len(G_PILOT)) {
    A   <- if (g <= 60) rlnorm(1, log(0.5), 0.3) else 0
    phi <- runif(1, 0, 24)
    M   <- rnorm(1, 5, 1)
    sig <- exp(rnorm(1, -1, 0.3))
    mat[g, ] <- rnorm(N_PILOT, M + A * cos(omega * pilot_times - omega * phi), sig)
  }
  mat
}
rownames(pilot_data) <- paste0("Gene", seq_len(G_PILOT))

# ---------------------------------------------------------------------------
# 3.  Options constructors
# ---------------------------------------------------------------------------
cat("\n--- Section 3: Options constructors ---\n")

# estCircadianParam
bio <- tryCatch(
  suppressWarnings(suppressMessages(
    estCircadianParam(pilot_data, pilot_times, period = 24, verbose = FALSE)
  )),
  error = function(e) NULL
)
.expect_true("estCircadianParam returns non-NULL", !is.null(bio))
.expect_class("estCircadianParam returns CircadianBioOptions", bio, "CircadianBioOptions")
.expect_true("bio$prop_rhythmic in [0,1]", !is.null(bio$prop_rhythmic) &&
             bio$prop_rhythmic >= 0 && bio$prop_rhythmic <= 1)
.expect_true("bio$lBaselineExpr is numeric", is.numeric(bio$lBaselineExpr) && length(bio$lBaselineExpr) > 0)
.expect_true("bio$lOD is numeric", is.numeric(bio$lOD) && length(bio$lOD) > 0)
.expect_true("bio$amplitude is numeric", is.numeric(bio$amplitude) && length(bio$amplitude) > 0)
.expect_true("bio$cts stored (same times as input)", !is.null(bio$cts) && length(bio$cts) == N_PILOT)

# CircadianBioOptions directly
bio_direct <- tryCatch(
  CircadianBioOptions(
    ngenes        = 200,
    lBaselineExpr = rnorm(200, 5, 1),
    lOD           = rnorm(200, -1, 0.3),
    amplitude     = abs(rnorm(80, 0.4, 0.2)),
    prop_rhythmic = 0.30,
    prop_DR       = 0.10,
    prop_DP       = 0.10,
    sim.seed      = 99
  ),
  error = function(e) NULL
)
.expect_class("CircadianBioOptions direct construction", bio_direct, "CircadianBioOptions")

# CircadianDesignOptions
design <- tryCatch(
  CircadianDesignOptions(
    sample_sizes = c(10, 20, 30),
    nsims        = 3L,
    design       = "active"
  ),
  error = function(e) NULL
)
.expect_class("CircadianDesignOptions construction", design, "CircadianDesignOptions")

# CircadianAnalysisOptions
analysis <- tryCatch(
  CircadianAnalysisOptions(alpha = 0.05, p.adjust.method = "BH"),
  error = function(e) NULL
)
.expect_class("CircadianAnalysisOptions construction", analysis, "CircadianAnalysisOptions")

# CircadianBootstrapOptions
boot_opts <- tryCatch(
  CircadianBootstrapOptions(
    design_vector = seq(0, 23, by = 2),
    B_values      = c(4, 6),
    N_values      = c(12, 24),
    nboot         = 2L,
    nsims_inner   = 2L,
    design        = "active",
    seed          = 42L
  ),
  error = function(e) NULL
)
.expect_class("CircadianBootstrapOptions construction", boot_opts, "CircadianBootstrapOptions")

# ---------------------------------------------------------------------------
# 4.  Simulation functions
# ---------------------------------------------------------------------------
cat("\n--- Section 4: Simulation functions ---\n")

# Active design CTS for testing
cts_active <- seq(0, 24, length.out = 21)[seq_len(20)]

# simCircadianSingleCohort
sc <- tryCatch(
  simCircadianSingleCohort(bio, cts_active, seed = 1L),
  error = function(e) { cat("  ERROR simCircadianSingleCohort:", conditionMessage(e), "\n"); NULL }
)
.expect_true("simCircadianSingleCohort returns non-NULL", !is.null(sc))
.expect_true("sc$expr is ngenes x N matrix",
             !is.null(sc) && is.matrix(sc$expr) &&
             nrow(sc$expr) == bio$ngenes && ncol(sc$expr) == length(cts_active))
.expect_true("sc$is_rhythmic is logical vector of length ngenes",
             !is.null(sc) && is.logical(sc$is_rhythmic) && length(sc$is_rhythmic) == bio$ngenes)
.expect_true("sc$r_values is numeric of length ngenes",
             !is.null(sc) && is.numeric(sc$r_values) && length(sc$r_values) == bio$ngenes)
.expect_true("sc$r_values non-negative", !is.null(sc) && all(sc$r_values >= 0, na.rm = TRUE))

# simCircadianDiff
sd_sim <- tryCatch(
  suppressWarnings(simCircadianDiff(
    ngenes = 200, n1 = 20, n2 = 20,
    lBaselineExpr = bio$lBaselineExpr,
    lOD           = bio$lOD,
    amplitude     = bio$amplitude,
    prop_rhythmic = 0.30,
    prop_DR       = 0.10, prop_DP = 0.10,
    period        = 24, design = "active",
    sim.seed      = 7L
  )),
  error = function(e) { cat("  ERROR simCircadianDiff:", conditionMessage(e), "\n"); NULL }
)
.expect_true("simCircadianDiff returns non-NULL", !is.null(sd_sim))
.expect_true("sd_sim$expr1 is 200 x 20",
             !is.null(sd_sim) && nrow(sd_sim$expr1) == 200 && ncol(sd_sim$expr1) == 20)
.expect_true("sd_sim ground_truth has 200 rows",
             !is.null(sd_sim) && nrow(sd_sim$ground_truth) == 200)
.expect_true("sd_sim diff_type values in 0:5",
             !is.null(sd_sim) && all(sd_sim$ground_truth$diff_type %in% 0:5))

# simCircadianFMM (requires FMM package)
if (requireNamespace("FMM", quietly = TRUE)) {
  sc_fmm1 <- tryCatch(
    simCircadianFMM(bio, cts_active, omega = 1.0, seed = 5L),
    error = function(e) { cat("  ERROR simCircadianFMM omega=1:", conditionMessage(e), "\n"); NULL }
  )
  .expect_true("simCircadianFMM(omega=1) returns non-NULL", !is.null(sc_fmm1))
  .expect_true("simCircadianFMM(omega=1) has same structure as simCircadianSingleCohort",
               !is.null(sc_fmm1) && is.matrix(sc_fmm1$expr) &&
               nrow(sc_fmm1$expr) == bio$ngenes && ncol(sc_fmm1$expr) == length(cts_active))
  .expect_true("simCircadianFMM stores omega field",
               !is.null(sc_fmm1) && !is.null(sc_fmm1$omega) && sc_fmm1$omega == 1.0)

  sc_fmm0 <- tryCatch(
    simCircadianFMM(bio, cts_active, omega = 0.0, seed = 5L),
    error = function(e) NULL
  )
  .expect_true("simCircadianFMM(omega=0) returns non-NULL", !is.null(sc_fmm0))

  sc_fmm5 <- tryCatch(
    simCircadianFMM(bio, cts_active, omega = 0.5, seed = 5L),
    error = function(e) NULL
  )
  .expect_true("simCircadianFMM(omega=0.5) returns non-NULL", !is.null(sc_fmm5))

  # FMM with omega outside [0,1] must error
  ok_err <- tryCatch({ simCircadianFMM(bio, cts_active, omega = 1.5); FALSE },
                     error = function(e) TRUE)
  .expect_true("simCircadianFMM(omega=1.5) errors as expected", ok_err)
} else {
  cat("  [SKIP] simCircadianFMM tests — FMM package not installed\n")
}

# simCircadianDiffFMM (requires FMM package)
if (requireNamespace("FMM", quietly = TRUE)) {
  sdiff_fmm <- tryCatch(
    suppressWarnings(simCircadianDiffFMM(
      ngenes = 200, n1 = 12, n2 = 12,
      lBaselineExpr = bio$lBaselineExpr,
      lOD           = bio$lOD,
      amplitude     = bio$amplitude,
      prop_rhythmic = 0.30, prop_DR = 0.10, prop_DP = 0.10,
      design = "active", sim.seed = 3L,
      omega = 0.5, beta = pi
    )),
    error = function(e) { cat("  ERROR simCircadianDiffFMM:", conditionMessage(e), "\n"); NULL }
  )
  .expect_true("simCircadianDiffFMM(omega=0.5) returns non-NULL", !is.null(sdiff_fmm))
  .expect_true("simCircadianDiffFMM has expr1/expr2 fields",
               !is.null(sdiff_fmm) && !is.null(sdiff_fmm$expr1) && !is.null(sdiff_fmm$expr2))

  # omega=1.0 must return same structure as simCircadianDiff
  sdiff_fmm1 <- tryCatch(
    suppressWarnings(simCircadianDiffFMM(
      ngenes = 200, n1 = 12, n2 = 12,
      lBaselineExpr = bio$lBaselineExpr,
      lOD           = bio$lOD,
      amplitude     = bio$amplitude,
      prop_rhythmic = 0.30, prop_DR = 0.10, prop_DP = 0.10,
      design = "active", sim.seed = 3L,
      omega = 1.0
    )),
    error = function(e) NULL
  )
  .expect_true("simCircadianDiffFMM(omega=1.0) returns same structure as simCircadianDiff",
               !is.null(sdiff_fmm1) && !is.null(sdiff_fmm1$ground_truth) &&
               nrow(sdiff_fmm1$ground_truth) == 200)
} else {
  cat("  [SKIP] simCircadianDiffFMM tests — FMM package not installed\n")
}

# ---------------------------------------------------------------------------
# 5.  Estimation functions
# ---------------------------------------------------------------------------
cat("\n--- Section 5: Estimation functions ---\n")

# one_cosinor_OLS basic correctness
ols_t   <- seq(0, 23, by = 2)
ols_y   <- 5 + 1.5 * cos(2 * pi / 24 * ols_t - 2 * pi / 24 * 6) + rnorm(12, 0, 0.2)
ols_res <- tryCatch(one_cosinor_OLS(ols_t, ols_y, period = 24),
                    error = function(e) NULL)
.expect_true("one_cosinor_OLS returns list with M, A, phi, pvalue", !is.null(ols_res) &&
             all(c("M","A","phi","pvalue") %in% names(ols_res)))
.expect_true("one_cosinor_OLS amplitude near truth (1.5 +/-0.5)",
             !is.null(ols_res) && !is.na(ols_res$A) && abs(ols_res$A - 1.5) < 0.5)
.expect_true("one_cosinor_OLS pvalue in [0,1]",
             !is.null(ols_res) && !is.na(ols_res$pvalue) &&
             ols_res$pvalue >= 0 && ols_res$pvalue <= 1)

# one_cosinor_OLS sigma denominator: n-3
n_obs   <- 12L
ols_yhat <- ols_res$M + ols_res$A * cos(2 * pi / 24 * ols_t - 2 * pi / 24 * ols_res$phi)
sigma_manual <- sqrt(sum((ols_y - ols_yhat)^2) / (n_obs - 3L))
# The denominator inside fitCosinorAll uses (n-3), verify sigma_vals matches
.expect_true("fitCosinorAll sigma uses n-3 denominator (manually verified)", {
  df_test <- fitCosinorAll(
    matrix(ols_y, nrow = 1), ols_t, period = 24
  )
  abs(df_test$sigma[1] - sigma_manual) < 1e-10
})

# one_cosinor_OLS single-observation edge case (n <= 3 must return NA pvalue)
ols_edge <- tryCatch(one_cosinor_OLS(c(0, 8, 16), c(1, 2, 3), period = 24),
                     error = function(e) NULL)
.expect_true("one_cosinor_OLS n<=3 returns NA pvalue", !is.null(ols_edge) && is.na(ols_edge$pvalue))

# estCircadianParamTwoGroup
bio2 <- tryCatch(
  suppressWarnings(suppressMessages(
    estCircadianParamTwoGroup(
      data_1 = pilot_data,  times_1 = pilot_times,
      data_2 = pilot_data,  times_2 = pilot_times,
      verbose = FALSE
    )
  )),
  error = function(e) { cat("  ERROR estCircadianParamTwoGroup:", conditionMessage(e), "\n"); NULL }
)
.expect_class("estCircadianParamTwoGroup returns CircadianBioOptions", bio2, "CircadianBioOptions")
.expect_true("bio2$diagnostics is a list", !is.null(bio2) && is.list(bio2$diagnostics))

# ---------------------------------------------------------------------------
# 6.  Detection functions
# ---------------------------------------------------------------------------
cat("\n--- Section 6: Detection functions ---\n")

small_expr  <- sc$expr[seq_len(50), seq_len(20)]
small_times <- cts_active[seq_len(20)]

# detect_DCP
pv_dcp <- tryCatch(detect_DCP(small_expr, small_times, period = 24),
                   error = function(e) { cat("  ERROR detect_DCP:", conditionMessage(e), "\n"); NULL })
.expect_true("detect_DCP returns numeric vector", !is.null(pv_dcp) && is.numeric(pv_dcp))
.expect_true("detect_DCP length equals nrow(expr)", !is.null(pv_dcp) && length(pv_dcp) == 50)
.expect_range("detect_DCP p-values in [0,1]", pv_dcp[!is.na(pv_dcp)])

# detect_JTK (optional — MetaCycle)
if (requireNamespace("MetaCycle", quietly = TRUE)) {
  pv_jtk <- tryCatch(detect_JTK(small_expr, small_times, period = 24),
                     error = function(e) NULL)
  .expect_true("detect_JTK returns length-50 numeric", !is.null(pv_jtk) && length(pv_jtk) == 50)
  .expect_range("detect_JTK p-values in [0,1]", pv_jtk[!is.na(pv_jtk)])
} else {
  cat("  [SKIP] detect_JTK tests — MetaCycle not installed\n")
}

# detect_RAIN (optional)
if (requireNamespace("rain", quietly = TRUE)) {
  pv_rain <- tryCatch(detect_RAIN(small_expr, small_times, period = 24),
                      error = function(e) NULL)
  .expect_true("detect_RAIN returns length-50 numeric", !is.null(pv_rain) && length(pv_rain) == 50)
  .expect_range("detect_RAIN p-values in [0,1]", pv_rain[!is.na(pv_rain)])
} else {
  cat("  [SKIP] detect_RAIN tests — rain not installed\n")
}

# ---------------------------------------------------------------------------
# 7.  fitCosinorAll
# ---------------------------------------------------------------------------
cat("\n--- Section 7: fitCosinorAll ---\n")

fca <- tryCatch(
  fitCosinorAll(pilot_data[seq_len(50), ], pilot_times, period = 24),
  error = function(e) { cat("  ERROR fitCosinorAll:", conditionMessage(e), "\n"); NULL }
)
.expect_true("fitCosinorAll returns data.frame", !is.null(fca) && is.data.frame(fca))
.expect_true("fitCosinorAll has 50 rows", !is.null(fca) && nrow(fca) == 50)
.expect_true("fitCosinorAll has is_rhythmic column", !is.null(fca) && "is_rhythmic" %in% names(fca))
.expect_true("fitCosinorAll sigma > 0 where not NA",
             !is.null(fca) && all(fca$sigma[!is.na(fca$sigma)] > 0))
# Verify denominator is n-3 (not n)
.expect_true("fitCosinorAll sigma denominator is n-3", {
  g <- 1L
  y <- pilot_data[g, ]
  fit <- one_cosinor_OLS(pilot_times, y, period = 24)
  yhat <- fit$M + fit$A * cos(2*pi/24 * pilot_times - 2*pi/24 * fit$phi)
  sigma_n3 <- sqrt(sum((y - yhat)^2) / (length(y) - 3L))
  abs(fca$sigma[g] - sigma_n3) < 1e-10
})

# ---------------------------------------------------------------------------
# 8.  npower function
# ---------------------------------------------------------------------------
cat("\n--- Section 8: npower ---\n")

# Build a tiny runSimsSingleCohort result to feed npower
design_np <- CircadianDesignOptions(sample_sizes = c(10, 20, 30), nsims = 3L, design = "active")
analysis_np <- CircadianAnalysisOptions(alpha = 0.05)

sc_res <- tryCatch(
  suppressWarnings(suppressMessages(
    runSimsSingleCohort(bio, design_np, analysis_np,
                        method = "DCP", verbose = FALSE)
  )),
  error = function(e) { cat("  ERROR runSimsSingleCohort for npower:", conditionMessage(e), "\n"); NULL }
)

if (!is.null(sc_res)) {
  np <- tryCatch(npower(sc_res, target_power = 0.80, fdr = 0.05),
                 error = function(e) { cat("  ERROR npower:", conditionMessage(e), "\n"); NULL })
  .expect_class("npower returns class 'npower'", np, "npower")
  .expect_true("npower has $n, $n_grid, $n_interp, $power",
               !is.null(np) && all(c("n","n_grid","n_interp","power") %in% names(np)))
  .expect_true("npower$power values in [0,1]",
               !is.null(np) && all(np$power >= 0 & np$power <= 1, na.rm = TRUE))
  .expect_true("npower interpolate=TRUE returns n <= n_grid when n_grid not NA", {
    if (!is.null(np) && !is.na(np$n_grid)) {
      is.na(np$n_interp) || np$n_interp <= np$n_grid
    } else {
      TRUE  # can't test if target never reached
    }
  })
  .expect_true("npower interpolate=FALSE returns n_grid", {
    np_nogrid <- tryCatch(npower(sc_res, target_power = 0.80, fdr = 0.05, interpolate = FALSE),
                          error = function(e) NULL)
    !is.null(np_nogrid) && identical(np_nogrid$n, np_nogrid$n_grid)
  })
} else {
  cat("  [SKIP] npower tests — runSimsSingleCohort failed\n")
}

# ---------------------------------------------------------------------------
# 9.  runSingleCohortGrid (core invariants)
# ---------------------------------------------------------------------------
cat("\n--- Section 9: runSingleCohortGrid ---\n")

design_grid <- CircadianDesignOptions(
  sample_sizes = c(12, 24),
  nsims        = 3L,
  design       = "active",
  B_values     = c(4, 6)
)

sc_grid <- tryCatch(
  suppressWarnings(suppressMessages(
    runSingleCohortGrid(bio, design_grid, analysis, methods = "DCP",
                        alpha2 = 0, mc.cores = 1L, verbose = FALSE)
  )),
  error = function(e) { cat("  ERROR runSingleCohortGrid:", conditionMessage(e), "\n"); NULL }
)
.expect_class("runSingleCohortGrid returns SCPSingleResult", sc_grid, "SCPSingleResult")
.expect_true("SCPSingleResult has $power_df",
             !is.null(sc_grid) && !is.null(sc_grid$power_df) && is.data.frame(sc_grid$power_df))
.expect_true("SCPSingleResult has $n80_df",
             !is.null(sc_grid) && !is.null(sc_grid$n80_df) && is.data.frame(sc_grid$n80_df))
.expect_range("power_df$power values in [0,1]",
              if (!is.null(sc_grid)) sc_grid$power_df$power else NA_real_)
.expect_true("power_df columns include N, B, power, power_se",
             !is.null(sc_grid) &&
             all(c("N","B","power","power_se") %in% names(sc_grid$power_df)))

# ---------------------------------------------------------------------------
# 10.  runSimsSingleCohort full structure
# ---------------------------------------------------------------------------
cat("\n--- Section 10: runSimsSingleCohort full structure ---\n")

design_sc <- CircadianDesignOptions(sample_sizes = c(10, 20, 30), nsims = 3L, design = "active")

sc_full <- tryCatch(
  suppressWarnings(suppressMessages(
    runSimsSingleCohort(bio, design_sc, analysis, method = "DCP", verbose = FALSE)
  )),
  error = function(e) { cat("  ERROR runSimsSingleCohort full:", conditionMessage(e), "\n"); NULL }
)
.expect_true("runSimsSingleCohort returns non-NULL", !is.null(sc_full))
.expect_true("$marginal_power is 3x3 matrix",
             !is.null(sc_full) && is.matrix(sc_full$marginal_power) &&
             nrow(sc_full$marginal_power) == 3 && ncol(sc_full$marginal_power) == 3)
.expect_range("$marginal_power in [0,1]",
              if (!is.null(sc_full)) sc_full$marginal_power[!is.na(sc_full$marginal_power)] else NA_real_)
.expect_true("$pvalues is array [3 x ngenes x 3]",
             !is.null(sc_full) && is.array(sc_full$pvalues) &&
             dim(sc_full$pvalues)[1] == 3 && dim(sc_full$pvalues)[2] == bio$ngenes)
.expect_true("$r_values_list is a list of length 3",
             !is.null(sc_full) && is.list(sc_full$r_values_list) &&
             length(sc_full$r_values_list) == 3)

# ---------------------------------------------------------------------------
# 11.  runSimsDiff
# ---------------------------------------------------------------------------
cat("\n--- Section 11: runSimsDiff ---\n")

design_diff <- CircadianDesignOptions(
  sample_sizes = c(10, 20),
  nsims        = 3L,
  design       = "active",
  test_types   = c("DR", "DP")
)
bio_diff <- updateBioOptions(bio,
  prop_DR = 0.10, prop_DP = 0.10, prop_DM = 0.00
)

sd_full <- tryCatch(
  suppressWarnings(suppressMessages(
    runSimsDiff(bio_diff, design_diff, analysis, verbose = FALSE)
  )),
  error = function(e) { cat("  ERROR runSimsDiff:", conditionMessage(e), "\n"); NULL }
)
.expect_true("runSimsDiff returns non-NULL", !is.null(sd_full))
.expect_true("$fdr_DR is array [ngenes x 2 x 3]",
             !is.null(sd_full) && is.array(sd_full$fdr_DR) &&
             dim(sd_full$fdr_DR)[1] == bio$ngenes && dim(sd_full$fdr_DR)[2] == 2)
.expect_true("$diff_type is a list of length 3",
             !is.null(sd_full) && is.list(sd_full$diff_type) && length(sd_full$diff_type) == 3)
.expect_range("$fdr_DR values in [0,1]",
              if (!is.null(sd_full)) sd_full$fdr_DR[!is.na(sd_full$fdr_DR)] else NA_real_)

# ---------------------------------------------------------------------------
# 12.  npower on differential result
# ---------------------------------------------------------------------------
cat("\n--- Section 12: npower on differential result ---\n")

if (!is.null(sd_full)) {
  np_dr <- tryCatch(npower(sd_full, target_power = 0.80, fdr = 0.05, endpoint = "DR"),
                    error = function(e) { cat("  ERROR npower(DR):", conditionMessage(e), "\n"); NULL })
  .expect_class("npower(endpoint='DR') returns class 'npower'", np_dr, "npower")
  .expect_true("npower(DR)$endpoint == 'DR'", !is.null(np_dr) && identical(np_dr$endpoint, "DR"))
  .expect_range("npower(DR)$power in [0,1]",
                if (!is.null(np_dr)) np_dr$power else NA_real_)

  # Must error when endpoint not specified on differential result
  err_ok <- tryCatch({ npower(sd_full); FALSE }, error = function(e) TRUE)
  .expect_true("npower on diff result errors without endpoint arg", err_ok)
}

# ---------------------------------------------------------------------------
# 13.  runBootstrapDesignGrid
# ---------------------------------------------------------------------------
cat("\n--- Section 13: runBootstrapDesignGrid ---\n")

boot_opts_small <- tryCatch(
  CircadianBootstrapOptions(
    design_vector = seq(0, 22, by = 2),
    B_values      = c(4, 6),
    N_values      = c(12, 24),
    nboot         = 2L,
    nsims_inner   = 2L,
    design        = "active",
    seed          = 7L
  ),
  error = function(e) NULL
)

br <- tryCatch(
  suppressWarnings(suppressMessages(
    runBootstrapDesignGrid(
      pilot_data    = pilot_data[seq_len(100), ],
      pilot_times   = pilot_times,
      boot.opts     = boot_opts_small,
      analysis.opts = analysis,
      bio_diff.opts = bio,
      mode          = "single",
      verbose       = FALSE,
      mc.cores      = 1L
    )
  )),
  error = function(e) { cat("  ERROR runBootstrapDesignGrid:", conditionMessage(e), "\n"); NULL }
)
.expect_true("runBootstrapDesignGrid returns non-NULL", !is.null(br))
.expect_true("$power_mean exists and is array",
             !is.null(br) && !is.null(br$power_mean) && is.array(br$power_mean))
.expect_true("$boot_power 4D array [nboot x n_N x n_B x n_tests]",
             !is.null(br) && is.array(br$boot_power) &&
             length(dim(br$boot_power)) == 4 && dim(br$boot_power)[1] == 2)
.expect_range("$power_mean in [0,1]",
              if (!is.null(br)) br$power_mean[!is.na(br$power_mean)] else NA_real_)

# ---------------------------------------------------------------------------
# 14.  circaPowerApproxN80
# ---------------------------------------------------------------------------
cat("\n--- Section 14: circaPowerApproxN80 ---\n")

n80_approx <- tryCatch(circaPowerApproxN80(bio, alpha = 0.05, target_power = 0.80),
                        error = function(e) NULL)
.expect_true("circaPowerApproxN80 returns single numeric or NA",
             !is.null(n80_approx) && (is.na(n80_approx) || (is.numeric(n80_approx) && length(n80_approx) == 1)))
.expect_true("circaPowerApproxN80 > 0 when not NA",
             !is.null(n80_approx) && (is.na(n80_approx) || n80_approx > 0))

# ---------------------------------------------------------------------------
# 15.  updateBioOptions
# ---------------------------------------------------------------------------
cat("\n--- Section 15: updateBioOptions ---\n")

bio_upd <- tryCatch(
  updateBioOptions(bio, prop_DR = 0.20, prop_DP = 0.05),
  error = function(e) NULL
)
.expect_class("updateBioOptions returns CircadianBioOptions", bio_upd, "CircadianBioOptions")
.expect_true("updateBioOptions$prop_DR updated to 0.20",
             !is.null(bio_upd) && abs(bio_upd$prop_DR - 0.20) < 1e-9)

# ---------------------------------------------------------------------------
# 16.  Invariants
# ---------------------------------------------------------------------------
cat("\n--- Section 16: Key invariants ---\n")

# I1: simCircadianFMM(omega=1) has same dimensions as simCircadianSingleCohort
if (requireNamespace("FMM", quietly = TRUE) && !is.null(sc)) {
  sc1_dims <- dim(sc$expr)
  fmm1 <- tryCatch(simCircadianFMM(bio, cts_active, omega = 1.0, seed = 1L), error = function(e) NULL)
  .expect_true("INVARIANT: simCircadianFMM(omega=1) expr dims == simCircadianSingleCohort",
               !is.null(fmm1) && identical(dim(fmm1$expr), sc1_dims))
}

# I2: npower interpolate=TRUE  =>  $n <= $n_grid always (when target reached)
if (!is.null(sc_full)) {
  np_interp <- tryCatch(npower(sc_full, target_power = 0.40, fdr = 0.05, interpolate = TRUE),
                        error = function(e) NULL)
  np_grid   <- tryCatch(npower(sc_full, target_power = 0.40, fdr = 0.05, interpolate = FALSE),
                        error = function(e) NULL)
  .expect_true("INVARIANT: npower(interpolate=TRUE)$n <= npower(interpolate=FALSE)$n_grid",
               !is.null(np_interp) && !is.null(np_grid) &&
               (is.na(np_interp$n) || is.na(np_grid$n_grid) ||
                np_interp$n <= np_grid$n_grid))
}

# I3: Power is in [0,1] for all grid cells in runSingleCohortGrid
if (!is.null(sc_grid)) {
  pwr <- sc_grid$power_df$power
  .expect_true("INVARIANT: all power_df$power in [0,1]",
               all(!is.na(pwr) & pwr >= 0 & pwr <= 1, na.rm = TRUE))
}

# I4: simCircadianDiff gene assignments are reproducible across sample sizes
set.seed(1)
d_n20 <- suppressWarnings(simCircadianDiff(
  ngenes=200, n1=20, n2=20,
  lBaselineExpr=bio$lBaselineExpr, lOD=bio$lOD, amplitude=bio$amplitude,
  prop_rhythmic=0.3, prop_DR=0.1, prop_DP=0.1,
  design="active", sim.seed=77L))
set.seed(1)
d_n40 <- suppressWarnings(simCircadianDiff(
  ngenes=200, n1=40, n2=40,
  lBaselineExpr=bio$lBaselineExpr, lOD=bio$lOD, amplitude=bio$amplitude,
  prop_rhythmic=0.3, prop_DR=0.1, prop_DP=0.1,
  design="active", sim.seed=77L))
.expect_true("INVARIANT: diff_type identical across sample sizes (n-independent gene assignment)",
             identical(d_n20$ground_truth$diff_type, d_n40$ground_truth$diff_type))

# I5: prop_rhythmic constraints — prop_DR + prop_DP > prop_rhythmic should error
err5 <- tryCatch(
  CircadianBioOptions(
    ngenes=100, lBaselineExpr=rnorm(100,5,1), lOD=rnorm(100,-1,0.3),
    amplitude=abs(rnorm(40,0.4,0.2)), prop_rhythmic=0.10,
    prop_DR=0.20, prop_DP=0.20),
  error = function(e) "ERROR"
)
.expect_true("INVARIANT: CircadianBioOptions errors when prop_DR+prop_DP > prop_rhythmic",
             identical(err5, "ERROR"))

# I6: circular_difference is in (-period/2, period/2)
phi_test <- runif(100, 0, 24)
delta    <- circular_difference(phi_test, rev(phi_test), period = 24)
.expect_true("INVARIANT: circular_difference in (-12, 12)",
             all(delta >= -12 & delta <= 12, na.rm = TRUE))

# I7: r_values = 0 for non-rhythmic genes in simCircadianSingleCohort
if (!is.null(sc)) {
  .expect_true("INVARIANT: r_values == 0 for non-rhythmic genes",
               all(sc$r_values[!sc$is_rhythmic] == 0))
}

# ---------------------------------------------------------------------------
# 17.  Print methods (should not error)
# ---------------------------------------------------------------------------
cat("\n--- Section 17: Print methods ---\n")

for (obj_name in c("bio", "design", "analysis", "boot_opts")) {
  obj <- tryCatch(get(obj_name), error = function(e) NULL)
  if (!is.null(obj)) {
    ok <- tryCatch({ capture.output(print(obj)); TRUE }, error = function(e) FALSE)
    .record(sprintf("print(%s) does not error", obj_name), ok)
  }
}
if (!is.null(np)) {
  ok <- tryCatch({ capture.output(print(np)); TRUE }, error = function(e) FALSE)
  .record("print.npower does not error", ok)
}

# ---------------------------------------------------------------------------
# 18.  Summary
# ---------------------------------------------------------------------------
total <- .PASS + .FAIL
cat(sprintf("\n=================================================\n"))
cat(sprintf("  %d/%d tests passed\n", .PASS, total))
if (.FAIL > 0) {
  cat(sprintf("  FAILED tests:\n"))
  for (r in .RESULTS) {
    if (r$status == "FAIL") {
      cat(sprintf("    - %s%s\n", r$name,
                  if (nchar(r$msg) > 0) paste0(": ", r$msg) else ""))
    }
  }
}
cat(sprintf("=================================================\n\n"))

invisible(.RESULTS)
