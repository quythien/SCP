#' Validation Suite for SCP Framework
#'
#' Call runSCPValidation() to execute all checks. Nothing runs at source() time.
#'
#' Validations:
#'   1. Type I error control under null
#'   2. Analytical vs simulation agreement (single-cohort)
#'   3. Phase effect size predicts power monotonically
#'   4. Power curve increases with sample size
#'   5. Differential power sanity (large > small effect; DR > subtle DP)
#'   6. CircaPower analytical formula consistency

# ==============================================================================
# Internal helper: pass/fail printer
# ==============================================================================
.vcat <- function(label, pass, detail = NULL) {
  sym <- if (pass) "✓ PASS" else "✗ FAIL"
  cat(sprintf("%-40s %s\n", label, sym))
  if (!is.null(detail)) cat("   ", detail, "\n")
}

# ==============================================================================
# Validation 1: Type I error control
# ==============================================================================
.val1_type_I_error <- function(nsims, n, ngenes, alpha) {
  cat("\n--- Validation 1: Type I Error Control ---\n")

  bio   <- CircadianBioOptions(ngenes = ngenes, prop_rhythmic = 0.25,
                               lBaselineExpr = rnorm(ngenes, 5, 2),
                               lOD = rnorm(ngenes, -1, 0.3),
                               amplitude = pmax(rlnorm(round(ngenes * 0.25), log(0.4), 0.5), 0.05),
                               prop_DR = 0, prop_DP = 0)
  des   <- CircadianDesignOptions(sample_sizes = n, nsims = nsims,
                                  design = "active")
  aopt  <- CircadianAnalysisOptions(alpha = alpha)

  out   <- runSimsDiff(bio, des, aopt)

  em_DR <- mean(out$pval_DR[, 1, ] < alpha, na.rm = TRUE)
  em_DP <- mean(out$pval_DP[, 1, ] < alpha, na.rm = TRUE)

  lo <- alpha * 0.6
  hi <- alpha * 1.4

  res <- data.frame(
    Test             = c("DR", "DP"),
    Empirical_alpha  = c(em_DR, em_DP),
    Pass             = c(em_DR >= lo & em_DR <= hi,
                         em_DP >= lo & em_DP <= hi)
  )
  print(res, row.names = FALSE)
  cat(sprintf("Acceptable range: [%.3f, %.3f]\n", lo, hi))
  .vcat("Type I error control", all(res$Pass))
  invisible(res)
}

# ==============================================================================
# Validation 2: Analytical vs simulation agreement
# ==============================================================================
.val2_analytical_vs_sim <- function(r_values, n_test, nsims, alpha) {
  cat("\n--- Validation 2: Analytical vs Simulation Agreement ---\n")

  res <- data.frame(
    r              = r_values,
    Analytical     = NA_real_,
    Simulation     = NA_real_,
    AbsDiff        = NA_real_,
    Agreement      = NA
  )

  for (i in seq_along(r_values)) {
    r      <- r_values[i]
    an     <- CircaPower(n = n_test, power = NULL, r = r, alpha = alpha,
                         design_factor = 0.5)

    sigma  <- 0.3
    A      <- r * sigma
    bio    <- CircadianBioOptions(ngenes = 500, prop_rhythmic = 1,
                                  lOD = log(sigma), amplitude = A,
                                  phase = "uniform")
    des    <- CircadianDesignOptions(sample_sizes = n_test, nsims = nsims,
                                    design = "active")
    aopt   <- CircadianAnalysisOptions(alpha = alpha)

    sim    <- runSimsSingleCohort(bio, des, aopt)
    sim_pw <- mean(sim$power[1, ], na.rm = TRUE)

    res$Analytical[i] <- an$power
    res$Simulation[i] <- sim_pw
    res$AbsDiff[i]    <- abs(an$power - sim_pw)
    res$Agreement[i]  <- res$AbsDiff[i] < 0.10
  }

  print(res, row.names = FALSE)
  .vcat("Analytical vs simulation (within 10%)", all(res$Agreement))
  invisible(res)
}

# ==============================================================================
# Validation 3: Phase effect size monotonicity
# ==============================================================================
.val3_phase_formula <- function(phase_shifts, n_test, nsims, alpha) {
  cat("\n--- Validation 3: Phase Effect Size Monotonicity ---\n")

  powers <- numeric(length(phase_shifts))

  for (i in seq_along(phase_shifts)) {
    shift <- phase_shifts[i]
    bio   <- CircadianBioOptions(ngenes = 300, prop_rhythmic = 0.5,
                                 lBaselineExpr = rnorm(300, 5, 2),
                                 lOD = rnorm(300, -1, 0.3),
                                 amplitude = pmax(rlnorm(round(300 * 0.5), log(0.4), 0.5), 0.05),
                                 prop_DR = 0, prop_DP = 1,
                                 phase_diff = c(shift, shift))
    des   <- CircadianDesignOptions(sample_sizes = n_test, nsims = nsims,
                                    design = "active")
    aopt  <- CircadianAnalysisOptions(alpha = alpha)

    sim   <- runSimsDiff(bio, des, aopt)
    powers[i] <- mean(sim$pval_DP[, 1, ] < alpha, na.rm = TRUE)
  }

  res <- data.frame(Phase_shift_h = phase_shifts, Power = round(powers, 3))
  print(res, row.names = FALSE)

  is_mono <- all(diff(powers) >= -0.05)  # allow small sampling noise
  .vcat("Power increases with phase shift", is_mono)
  invisible(res)
}

# ==============================================================================
# Validation 4: Monotonicity with sample size
# ==============================================================================
.val4_monotonicity <- function(sample_sizes, nsims, alpha) {
  cat("\n--- Validation 4: Power Curve Monotonicity ---\n")

  bio  <- CircadianBioOptions(ngenes = 300, prop_rhythmic = 0.3,
                              lBaselineExpr = rnorm(300, 5, 2),
                              lOD = rnorm(300, -1, 0.3),
                              amplitude = pmax(rlnorm(round(300 * 0.3), log(0.4), 0.5), 0.05))
  des  <- CircadianDesignOptions(sample_sizes = sample_sizes, nsims = nsims,
                                 design = "active")
  aopt <- CircadianAnalysisOptions(alpha = alpha)

  sim  <- runSimsSingleCohort(bio, des, aopt)
  pw   <- rowMeans(sim$power, na.rm = TRUE)

  res <- data.frame(N = sample_sizes, Power = round(pw, 3))
  print(res, row.names = FALSE)

  is_mono <- all(diff(pw) >= -0.03)
  if (!is_mono) {
    drops <- which(diff(pw) < -0.03)
    cat(sprintf("  Non-monotone at N: %s\n",
                paste(sample_sizes[drops], collapse = ", ")))
  }
  .vcat("Power increases monotonically with N", is_mono)
  invisible(res)
}

# ==============================================================================
# Validation 5: Differential power sanity
# ==============================================================================
.val5_differential_sanity <- function(n_test, nsims, alpha) {
  cat("\n--- Validation 5: Differential Power Sanity ---\n")

  run_dp <- function(phase_shift) {
    bio  <- CircadianBioOptions(ngenes = 300, prop_rhythmic = 0.5,
                                lBaselineExpr = rnorm(300, 5, 2),
                                lOD = rnorm(300, -1, 0.3),
                                amplitude = pmax(rlnorm(round(300 * 0.5), log(0.4), 0.5), 0.05),
                                prop_DR = 0, prop_DP = 0.5,
                                phase_diff = c(phase_shift, phase_shift))
    des  <- CircadianDesignOptions(sample_sizes = n_test, nsims = nsims,
                                   design = "active")
    aopt <- CircadianAnalysisOptions(alpha = alpha)
    sim  <- runSimsDiff(bio, des, aopt)
    mean(sim$pval_DP[, 1, ] < alpha, na.rm = TRUE)
  }

  run_dr <- function() {
    bio  <- CircadianBioOptions(ngenes = 300, prop_rhythmic = 0.5,
                                lBaselineExpr = rnorm(300, 5, 2),
                                lOD = rnorm(300, -1, 0.3),
                                amplitude = pmax(rlnorm(round(300 * 0.5), log(0.4), 0.5), 0.05),
                                prop_DR = 0.5, prop_DP = 0)
    des  <- CircadianDesignOptions(sample_sizes = n_test, nsims = nsims,
                                   design = "active")
    aopt <- CircadianAnalysisOptions(alpha = alpha)
    sim  <- runSimsDiff(bio, des, aopt)
    mean(sim$pval_DR[, 1, ] < alpha, na.rm = TRUE)
  }

  pw_small <- run_dp(2)
  pw_large <- run_dp(8)
  pw_DR    <- run_dr()

  cat(sprintf("  Power DP 2h shift: %.3f\n", pw_small))
  cat(sprintf("  Power DP 8h shift: %.3f\n", pw_large))
  cat(sprintf("  Power DR (lose rhythm): %.3f\n", pw_DR))

  t1 <- pw_large > pw_small
  .vcat("Large phase shift > small phase shift", t1)
  .vcat("DR >= DP (2h) [note: failure is a warning]",
        pw_DR >= pw_small,
        if (pw_DR < pw_small) "DR unexpectedly weaker than 2h DP" else NULL)

  invisible(list(pw_small = pw_small, pw_large = pw_large, pw_DR = pw_DR,
                 pass_t1 = t1))
}

# ==============================================================================
# Validation 6: CircaPower analytical consistency
# ==============================================================================
.val6_circapower <- function(r_values, alpha) {
  cat("\n--- Validation 6: CircaPower Analytical Consistency ---\n")

  # Check that n_for_power=0.8 decreases as r increases
  ns <- vapply(r_values, function(r)
    CircaPower(n = NULL, power = 0.8, r = r, alpha = alpha)$n,
    numeric(1))

  res <- data.frame(r = r_values, n_for_80pct_power = ns)
  print(res, row.names = FALSE)

  is_decreasing <- all(diff(ns) <= 0)
  .vcat("Required N decreases as r increases", is_decreasing)
  invisible(res)
}

# ==============================================================================
# Main entry point
# ==============================================================================

#' Run the SCP validation suite
#'
#' @param nsims     Simulation replicates per validation (default 100).
#' @param n_test    Sample size for single-n validations (default 24).
#' @param alpha     Nominal significance level (default 0.05).
#' @param seed      RNG seed for reproducibility (default 42).
#' @export
runSCPValidation <- function(nsims    = 100L,
                             n_test   = 24L,
                             alpha    = 0.05,
                             seed     = 42L) {

  set.seed(seed)

  cat("SCP validation suite\n")
  cat(sprintf("  nsims=%d  n_test=%d  alpha=%.2f  seed=%d\n\n",
              nsims, n_test, alpha, seed))

  v1 <- .val1_type_I_error(nsims = nsims, n = n_test,
                            ngenes = 500L, alpha = alpha)

  v2 <- .val2_analytical_vs_sim(r_values = c(0.5, 1.0, 1.5, 2.0),
                                 n_test = n_test, nsims = nsims, alpha = alpha)

  v3 <- .val3_phase_formula(phase_shifts = c(2, 4, 6, 8, 12),
                             n_test = n_test, nsims = nsims, alpha = alpha)

  v4 <- .val4_monotonicity(sample_sizes = c(6, 12, 18, 24, 36, 48),
                            nsims = nsims, alpha = alpha)

  v5 <- .val5_differential_sanity(n_test = n_test, nsims = nsims, alpha = alpha)

  v6 <- .val6_circapower(r_values = c(0.5, 1.0, 1.5, 2.0), alpha = alpha)

  # ---- Summary ----
  cat("\nSummary\n")

  .vcat("Val 1 (Type I error)",           all(v1$Pass))
  .vcat("Val 2 (Analytical match)",        all(v2$Agreement))
  .vcat("Val 3 (Phase monotonicity)",      all(diff(v3$Power) >= -0.05))
  .vcat("Val 4 (N monotonicity)",          all(diff(v4$Power) >= -0.03))
  .vcat("Val 5 (Diff power sanity)",       v5$pass_t1)
  .vcat("Val 6 (CircaPower consistency)",  all(diff(v6$n_for_80pct_power) <= 0))
  cat("\n")

  invisible(list(v1 = v1, v2 = v2, v3 = v3, v4 = v4, v5 = v5, v6 = v6))
}
