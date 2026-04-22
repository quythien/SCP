#' Power Calculation for Circadian Rhythm Detection
#'
#' Power calculation for circadian rhythm detection
#'
#' @param simOutput Output from runSims()
#' @param alpha.type "pval" or "fdr"
#' @param alpha.nominal Significance threshold (default 0.05)
#' @param stratify.by "effectsize" or "expr"
#' @param strata Stratification breaks
#' @param filter.by "none" or "expr"
#' @param target.by "effectsize" - which genes to consider as "interesting"
#' @param delta Minimum effect size for "interesting" genes

comparePower <- function(simOutput,
                         alpha.type = c("pval", "fdr"),
                         alpha.nominal = 0.05,
                         stratify.by = c("effectsize", "expr"),
                         strata,
                         filter.by = c("none", "expr"),
                         target.by = "effectsize",
                         delta = 0.5) {

  alpha.type = match.arg(alpha.type)
  stratify.by = match.arg(stratify.by)
  filter.by = match.arg(filter.by)

  # Set strata if not given
  if (missing(strata)) {
    if (stratify.by == "effectsize") {
      # Effect size (r = A/σ) strata
      strata = c(0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, Inf)
    } else {
      # Expression strata
      strata = c(0, 10, 2^(1:7)*10, Inf)
    }
  }

  # Parameters
  sample_sizes = simOutput$sample_sizes
  ngenes = simOutput$sim.opts$ngenes
  nsims = dim(simOutput$pvalue)[3]

  # Initialize results
  nr = length(strata) - 1
  TD = FD = FDR = alpha = power = array(NA, dim = c(nr, length(sample_sizes), nsims))
  power.marginal = alpha.marginal = FDR.marginal = matrix(NA, length(sample_sizes), nsims)

  # Loop over simulations and sample sizes
  for (i in 1:nsims) {
    for (j in seq_along(sample_sizes)) {

      # Get rhythmic flags
      rhythmic_id = simOutput$rhythmic_id[[i]]
      effectsize = simOutput$effectsize[[i]]

      Zg = Zg2 = rep(0, ngenes)
      Zg[rhythmic_id] = 1

      # Find target (interesting) genes
      if (target.by == "effectsize") {
        ix = abs(effectsize) > delta
      } else {
        ix = rep(TRUE, length(rhythmic_id))
      }
      Zg2[rhythmic_id[ix]] = 1

      # Stratification
      if (stratify.by == "effectsize") {
        xgr = cut(effectsize, strata)
      } else {
        xgr = cut(simOutput$xbar[, j, i], strata)
      }

      # Get p-values or FDR
      if (alpha.type == "pval") {
        x = simOutput$pvalue[, j, i]
      } else {
        x = simOutput$fdr[, j, i]
      }

      # Calculate power
      power00 = POWER1(x, p.crit = alpha.nominal, Zg, Zg2, xgr = xgr)

      TD[, j, i] = power00$TD
      FD[, j, i] = power00$FD
      alpha[, j, i] = power00$alpha
      alpha.marginal[j, i] = power00$alpha.marginal
      power[, j, i] = power00$power
      power.marginal[j, i] = power00$power.marginal
      FDR[, j, i] = power00$FDR
      FDR.marginal[j, i] = power00$FDR.marginal
    }
  }

  output <- list(
    TD = TD, FD = FD, FDR = FDR, alpha = alpha, power = power,
    alpha.marginal = alpha.marginal, power.marginal = power.marginal,
    FDR.marginal = FDR.marginal,
    # Input parameters
    alpha.type = alpha.type, alpha.nominal = alpha.nominal,
    stratify.by = stratify.by, strata = strata,
    target.by = target.by, sample_sizes = simOutput$sample_sizes,
    delta = delta
  )

  return(output)
}

#' Core power calculation
#'
#' @param p p-values or FDR
#' @param p.crit Significance threshold
#' @param Zg Indicator for all rhythmic genes
#' @param Zg2 Indicator for target (interesting) genes
#' @param xgr Stratification groups

POWER1 <- function(p, p.crit, Zg, Zg2, xgr) {

  ix.D = p <= p.crit
  N = sum(ix.D)  # Total discoveries
  N.stratified = tapply(ix.D, xgr, sum)

  # True Discoveries (TD)
  id.TP = Zg2 == 1
  TD = tapply(p[id.TP] <= p.crit, xgr[id.TP], sum)
  TD[is.na(TD)] = 0

  # False Discoveries (FD)
  id.TN = Zg == 0
  FD = tapply(p[id.TN] <= p.crit, xgr[id.TN], sum)
  FD[is.na(FD)] = 0

  # Type I error
  alpha = as.vector(FD / table(xgr[id.TN]))
  alpha.marginal = sum(FD) / sum(id.TN)

  # Power
  power = as.vector(TD / table(xgr[id.TP]))
  power[is.nan(power)] = 0
  power.marginal = sum(TD, na.rm = TRUE) / sum(id.TP)

  # FDR
  FDR = FD / N.stratified
  FDR.marginal = sum(FD, na.rm = TRUE) / N

  return(list(
    TD = TD, FD = FD, alpha.nominal = p.crit,
    alpha = alpha, alpha.marginal = alpha.marginal,
    power = power, power.marginal = power.marginal,
    FDR = FDR, FDR.marginal = FDR.marginal
  ))
}

#' Summary of power results
#'
#' @param powerOutput Output from comparePower()

summaryPower <- function(powerOutput) {

  sample_sizes = powerOutput$sample_sizes

  alpha.type = powerOutput$alpha.type
  if (alpha.type == "pval") {
    alpha.nam = "type I error"
    alpha.mar = rowMeans(powerOutput$alpha.marginal)
  } else {
    alpha.nam = "FDR"
    alpha.mar = rowMeans(powerOutput$FDR.marginal)
  }

  TD.avg = colSums(apply(powerOutput$TD, c(1, 2), mean, na.rm = TRUE))
  FD.avg = colSums(apply(powerOutput$FD, c(1, 2), mean, na.rm = TRUE))

  res = cbind(
    sample_sizes,
    powerOutput$alpha.nominal,
    alpha.mar,
    rowMeans(powerOutput$power.marginal),
    TD.avg,
    FD.avg,
    FD.avg / TD.avg
  )

  colnames(res) = c("n", paste(c("Nominal", "Actual"), alpha.nam),
                    "Marginal power", "Avg # of TD", "Avg # of FD", "FDC")

  print(signif(res, 2))
  return(invisible(res))
}
#' Power Calculation for Differential Rhythmicity Analysis
#'
#' Calculate power to detect differences in circadian rhythms between two groups
#'
#' @param simOutput Output from runSimsDiff()
#' @param test_type Which test to evaluate ("DR", "DP", or "DM")
#' @param alpha.type "pval" or "fdr"
#' @param alpha.nominal Significance threshold (default 0.05)
#' @param target_effect Minimum effect size to consider as "interesting"
#'
#' @return List with power statistics

comparePowerDiff <- function(simOutput,
                         test_type = c("DR", "DP", "DM"),
                         alpha.type = c("pval", "fdr"),
                         alpha.nominal = 0.05,
                         target_effect = 0.3) {

  test_type = match.arg(test_type)
  alpha.type = match.arg(alpha.type)

  sample_sizes = simOutput$sample_sizes
  ngenes = simOutput$ngenes
  nsims = simOutput$nsims

  # Get p-values for the test type
  if (test_type == "DR") {
    pvalues = simOutput$pval_DR
    fdr_values = simOutput$fdr_DR
  } else if (test_type == "DP") {
    pvalues = simOutput$pval_DP
    fdr_values = simOutput$fdr_DP
  } else if (test_type == "DM") {
    pvalues = simOutput$pval_DM
    fdr_values = simOutput$fdr_DM
  }

  # Get ground truth
  alpha_values = if (alpha.type == "pval") pvalues else fdr_values

  # Initialize power arrays
  power.marginal = matrix(NA, length(sample_sizes), nsims)
  FDR.marginal   = matrix(NA, length(sample_sizes), nsims)
  TD = numeric(length(sample_sizes))
  FD = numeric(length(sample_sizes))

  # Loop over simulations and sample sizes
  for (i in 1:nsims) {
    diff_type = simOutput$diff_type[[i]]

    if (test_type == "DR") {
      is_target = diff_type %in% c(2, 3)
    } else if (test_type == "DP") {
      effectsize_phase = simOutput$effectsize[[i]]$phase
      is_target = diff_type == 4 & effectsize_phase >= target_effect
    } else if (test_type == "DM") {
      is_target = diff_type == 5
    }

    is_null = !is_target

    for (j in seq_along(sample_sizes)) {
      pvals = alpha_values[, j, i]

      td_j = sum(pvals[is_target] < alpha.nominal, na.rm = TRUE)
      fd_j = sum(pvals[is_null]   < alpha.nominal, na.rm = TRUE)

      # Accumulate across sims
      TD[j] = TD[j] + td_j
      FD[j] = FD[j] + fd_j

      n_target = sum(is_target)
      power.marginal[j, i] = if (n_target > 0) td_j / n_target else NA

      N = td_j + fd_j
      FDR.marginal[j, i] = if (N > 0) fd_j / N else 0
    }
  }

  # Average per-sim TD/FD across simulations
  TD_avg = TD / nsims
  FD_avg = FD / nsims

  power_avg = rowMeans(power.marginal, na.rm = TRUE)
  fdr_avg   = rowMeans(FDR.marginal,   na.rm = TRUE)

  res = cbind(
    sample_sizes,
    alpha.nominal,
    power_avg,
    fdr_avg,
    TD_avg,
    FD_avg,
    ifelse(TD_avg > 0, FD_avg / TD_avg, NA)
  )

  colnames(res) = c("n", "alpha", "Power", "FDR",
                   "Avg_TD", "Avg_FD", "FDC")

  output <- list(
    power = power.marginal,
    FDR = FDR.marginal,
    power_avg = power_avg,
    fdr_avg = fdr_avg,
    TD = TD,
    FD = FD,
    summary = res,
    sample_sizes = sample_sizes,
    test_type = test_type,
    alpha.type = alpha.type,
    alpha.nominal = alpha.nominal,
    target_effect = target_effect
  )

  class(output) = c("powerDiff", "list")
  return(output)
}

#' Print summary of differential power results
#'
#' @param powerOutput Output from comparePowerDiff()

summaryPowerDiff <- function(powerOutput) {

  cat("\n=== Differential Power Summary ===\n\n")
  cat(sprintf("Test type: %s\n", powerOutput$test_type))
  cat(sprintf("Alpha: %.2f (%s)\n",
              powerOutput$alpha.nominal, powerOutput$alpha.type))
  cat(sprintf("Target effect: %.2f\n\n", powerOutput$target_effect))

  print(signif(powerOutput$summary, 3))
  cat("\n")

  invisible(powerOutput$summary)
}

#' Plot power curves for differential analysis
#'
#' @param powerOutput Output from comparePowerDiff()
#' @param title Plot title
#' @param color Line color
#'
#' @return ggplot object

plotPowerDiff <- function(powerOutput,
                      title = "Power to Detect Differential Rhythmicity",
                      color = "steelblue") {

  library(ggplot2)

  sample_sizes = powerOutput$sample_sizes
  power_avg = powerOutput$power_avg

  # Calculate SE
  power_se = apply(powerOutput$power, 1, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

  df = data.frame(
    n = sample_sizes,
    power = power_avg,
    se = power_se
  )

  p = ggplot(df, aes(x = n, y = power)) +
    geom_ribbon(aes(ymin = pmax(0, power - se),
                  ymax = pmin(1, power + se)), fill = color, alpha = 0.3) +
    geom_line(linewidth = 1.2, color = color) +
    geom_point(size = 3, color = color) +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
    labs(
      title = title,
      subtitle = sprintf("Test: %s, Target effect: %.2f",
                     powerOutput$test_type, powerOutput$target_effect),
      x = "Sample Size (per group)",
      y = "Power"
    ) +
    theme_bw() +
    ylim(0, 1) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))

  return(p)
}

#' Compare power across different test types
#'
#' @param simOutput Output from runSimsDiff()
#' @param target_effect Minimum effect size
#' @param alpha Significance level
#'
#' @return ggplot object comparing all tests

compareAllDifferentialTests <- function(simOutput,
                                     target_effect = 0.3,
                                     alpha = 0.05) {

  library(ggplot2)

  # Calculate power for each test type
  power_DR = comparePowerDiff(simOutput, test_type = "DR",
                            target_effect = target_effect,
                            alpha.nominal = alpha)
  power_DP = comparePowerDiff(simOutput, test_type = "DP",
                            target_effect = target_effect,
                            alpha.nominal = alpha)

  sample_sizes = simOutput$sample_sizes

  df = data.frame(
    n = rep(sample_sizes, 2),
    power = c(power_DR$power_avg, power_DP$power_avg),
    test = rep(c("DR", "DP"), each = length(sample_sizes))
  )

  p = ggplot(df, aes(x = n, y = power, color = test, group = test)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("DR" = "steelblue",
                              "DP" = "darkgreen")) +
    labs(
      title = "Power Comparison: Differential Tests",
      subtitle = sprintf("Target effect: %.2f, Alpha: %.2f", target_effect, alpha),
      x = "Sample Size (per group)",
      y = "Power",
      color = "Test Type"
    ) +
    theme_bw() +
    ylim(0, 1) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          legend.position = "bottom")

  return(list(
    plot = p,
    power_DR = power_DR,
    power_DP = power_DP
  ))
}

one_cosinor_OLS = function(tod, y, period = 24, compute.phase.CI = FALSE, CI.level = 0.95) {
  n     = length(tod)
  omega = 2 * pi / period
  x1    = cos(omega * tod)
  x2    = sin(omega * tod)

  mat.S = matrix(c(n, sum(x1), sum(x2),
                   sum(x1), sum(x1^2), sum(x1 * x2),
                   sum(x2), sum(x1 * x2), sum(x2^2)),
                 nrow = 3, byrow = TRUE)
  vec.d = c(sum(y), sum(y * x1), sum(y * x2))

  mat.S.inv = tryCatch(solve(mat.S), error = function(e) NULL)
  if (is.null(mat.S.inv)) return(list(pvalue = NA, M = NA, A = NA, phi = NA))

  est        = mat.S.inv %*% vec.d
  m.hat      = est[1]
  beta1.hat  = est[2]
  beta2.hat  = est[3]
  A.hat      = sqrt(beta1.hat^2 + beta2.hat^2)
  phase.hat  = adjust.to.2pi(atan2(beta2.hat, beta1.hat)) / omega

  TSS  = sum((y - mean(y))^2)
  yhat = m.hat + beta1.hat * x1 + beta2.hat * x2
  RSS  = sum((y - yhat)^2)
  MSS  = TSS - RSS

  if (n <= 3 || RSS < 0) return(list(pvalue = NA, M = m.hat, A = A.hat, phi = phase.hat))

  Fstat = (MSS / 2) / (RSS / (n - 3))
  pval  = stats::pf(Fstat, 2, n - 3, lower.tail = FALSE)

  list(M = m.hat, A = A.hat, phi = phase.hat, pvalue = pval, R2 = MSS / TSS)
}

#' Adjust angle to [0, 2*pi)
adjust.to.2pi = function(x) {
  x %% (2 * pi)
}
#' Power Analysis for Two-Group Circadian Comparison
#'
#' @description Simulation-based power analysis for detecting differences in
#' circadian rhythms between two conditions (e.g., Young vs Old).
#' This is the NOVEL contribution - no existing tool addresses this.
#'
#' @param params_A Parameters for condition A (from pilot data or manual)
#' @param params_B Parameters for condition B
#' @param times Time points (same for both conditions)
#' @param n_A_range Sample sizes per time point for condition A
#' @param n_B_range Sample sizes per time point for condition B
#' @param diff_type "amplitude", "phase", "rhythmicity", "mixed"
#' @param delta_A Magnitude of amplitude difference to detect
#' @param delta_phi Magnitude of phase difference to detect (hours)
#' @param prop_diff Proportion of genes with true difference
#' @param n_sim Number of simulations
#' @param G Total genes per simulation
#' @param alpha Significance threshold
#' @param test_type "amplitude", "phase", "joint"
#' @param verbose Print progress
#'
#' @return List with power results for each difference type
power_two_group = function(params_A = NULL,
                           params_B = NULL,
                           times = seq(0, 48, by = 4),
                           n_A_range = 1:6,
                           n_B_range = 1:6,
                           diff_type = "amplitude",
                           delta_A = 0.3,
                           delta_phi = 4,
                           prop_diff = 0.1,
                           n_sim = 100,
                           G = 500,
                           alpha = 0.05,
                           test_type = "amplitude",
                           verbose = TRUE) {

  start_time = Sys.time()

  if (verbose) {
    cat("=== Power Analysis for Two-Group Comparison ===\n")
    cat("Difference type:", diff_type, "\n")
    cat("Test type:", test_type, "\n")
    cat("Delta A:", delta_A, "\n")
    cat("Delta phi:", delta_phi, "hours\n")
    cat("n_A range:", paste(n_A_range, collapse = ", "), "\n")
    cat("n_B range:", paste(n_B_range, collapse = ", "), "\n")
    cat("Simulations:", n_sim, "\n\n")
  }

  # Create sample size grid
  grid = expand.grid(n_A = n_A_range, n_B = n_B_range)
  n_grid = nrow(grid)

  # Results storage
  results = list(
    power = matrix(NA, nrow = length(n_A_range), ncol = length(n_B_range),
                   dimnames = list(n_A = n_A_range, n_B = n_B_range)),
    fdr = matrix(NA, nrow = length(n_A_range), ncol = length(n_B_range),
                 dimnames = list(n_A = n_A_range, n_B = n_B_range)),
    sensitivity = matrix(NA, nrow = length(n_A_range), ncol = length(n_B_range),
                         dimnames = list(n_A = n_A_range, n_B = n_B_range)),
    specificity = matrix(NA, nrow = length(n_A_range), ncol = length(n_B_range),
                         dimnames = list(n_A = n_A_range, n_B = n_B_range))
  )

  # Loop over sample size combinations
  for (g in 1:n_grid) {
    n_A = grid$n_A[g]
    n_B = grid$n_B[g]

    if (verbose && g %% 10 == 0) {
      cat(sprintf("Evaluating n_A=%d, n_B=%d (%.0f%% complete)\n",
                  n_A, n_B, 100 * g / n_grid))
    }

    # Storage for this combination
    power_vec = numeric(n_sim)
    fdr_vec = numeric(n_sim)
    sens_vec = numeric(n_sim)
    spec_vec = numeric(n_sim)

    for (sim in 1:n_sim) {
      # Generate two-group data
      sim_data = simulate_two_group(
        G = G,
        n_A = n_A,
        n_B = n_B,
        times = times,
        diff_type = diff_type,
        delta_A = delta_A,
        delta_phi = delta_phi,
        prop_diff = prop_diff,
        params = params_A,
        seed = NULL
      )

      ground_truth = sim_data$ground_truth$has_diff

      # Fit cosinor to each condition
      fit_A = t(apply(sim_data$data_A, 1, function(y) {
        fit = one_cosinor_OLS(sim_data$times_A, y, period = 24)
        c(A = fit$A, phi = fit$phi, pvalue = fit$pvalue)
      }))

      fit_B = t(apply(sim_data$data_B, 1, function(y) {
        fit = one_cosinor_OLS(sim_data$times_B, y, period = 24)
        c(A = fit$A, phi = fit$phi, pvalue = fit$pvalue)
      }))

      # Test for difference
      if (test_type == "amplitude") {
        pvals = test_amplitude_difference(
          fit_A[, "A"], fit_B[, "A"],
          sim_data$data_A, sim_data$data_B,
          sim_data$times_A, sim_data$times_B
        )
      } else if (test_type == "phase") {
        pvals = test_phase_difference(
          fit_A[, "phi"], fit_B[, "phi"],
          fit_A[, "A"], fit_B[, "A"]
        )
      } else if (test_type == "joint") {
        pvals_amp = test_amplitude_difference(
          fit_A[, "A"], fit_B[, "A"],
          sim_data$data_A, sim_data$data_B,
          sim_data$times_A, sim_data$times_B
        )
        pvals_phase = test_phase_difference(
          fit_A[, "phi"], fit_B[, "phi"],
          fit_A[, "A"], fit_B[, "A"]
        )
        # Fisher's combined test
        pvals = pchisq(-2 * (log(pvals_amp) + log(pvals_phase)), df = 4, lower.tail = FALSE)
      }

      # Multiple testing correction
      qvals = p.adjust(pvals, method = "BH")

      # Classifications
      discovered = qvals < alpha
      true_positive = discovered & ground_truth
      false_positive = discovered & !ground_truth
      false_negative = !discovered & ground_truth
      true_negative = !discovered & !ground_truth

      # Metrics
      TP = sum(true_positive, na.rm = TRUE)
      FP = sum(false_positive, na.rm = TRUE)
      FN = sum(false_negative, na.rm = TRUE)
      TN = sum(true_negative, na.rm = TRUE)

      power_vec[sim] = TP / max(TP + FN, 1)
      fdr_vec[sim] = FP / max(TP + FP, 1)
      sens_vec[sim] = TP / max(TP + FN, 1)
      spec_vec[sim] = TN / max(TN + FP, 1)
    }

    # Store averaged results
    i = which(n_A_range == n_A)
    j = which(n_B_range == n_B)
    results$power[i, j] = mean(power_vec, na.rm = TRUE)
    results$fdr[i, j] = mean(fdr_vec, na.rm = TRUE)
    results$sensitivity[i, j] = mean(sens_vec, na.rm = TRUE)
    results$specificity[i, j] = mean(spec_vec, na.rm = TRUE)
  }

  # Find optimal sample sizes for 80% power
  optimal = find_optimal_sample_size(results$power, n_A_range, n_B_range, target = 0.8)

  # End time
  end_time = Sys.time()
  elapsed = difftime(end_time, start_time, units = "mins")

  if (verbose) {
    cat("\n=== Results ===\n")
    cat(sprintf("Computation time: %.1f minutes\n", as.numeric(elapsed)))
    cat("\nOptimal sample sizes for 80% power:\n")
    cat(sprintf("  n_A = %d, n_B = %d\n", optimal$n_A, optimal$n_B))
  }

  return(list(
    power_matrix = results$power,
    fdr_matrix = results$fdr,
    sensitivity_matrix = results$sensitivity,
    specificity_matrix = results$specificity,
    optimal = optimal,
    parameters = list(
      times = times,
      n_A_range = n_A_range,
      n_B_range = n_B_range,
      diff_type = diff_type,
      test_type = test_type,
      delta_A = delta_A,
      delta_phi = delta_phi,
      prop_diff = prop_diff,
      n_sim = n_sim,
      G = G,
      alpha = alpha
    ),
    elapsed_time = elapsed
  ))
}


#' Test for amplitude difference between two conditions
#'
#' @param A_A Amplitude estimates from condition A
#' @param A_B Amplitude estimates from condition B
#' @param data_A Data matrix for condition A
#' @param data_B Data matrix for condition B
#' @param times_A Time points for condition A
#' @param times_B Time points for condition B
#'
#' @return Vector of p-values
test_amplitude_difference = function(A_A, A_B, data_A, data_B, times_A, times_B) {
  G = length(A_A)
  pvals = numeric(G)

  for (g in 1:G) {
    # Get residuals from cosinor fit
    fit_A = one_cosinor_OLS(times_A, data_A[g, ], period = 24)
    fit_B = one_cosinor_OLS(times_B, data_B[g, ], period = 24)

    # Estimate pooled variance
    n_A = length(times_A)
    n_B = length(times_B)

    yhat_A = fit_A$M + fit_A$A * cos(2*pi/24 * times_A - 2*pi/24 * fit_A$phi)
    yhat_B = fit_B$M + fit_B$A * cos(2*pi/24 * times_B - 2*pi/24 * fit_B$phi)

    RSS_A = sum((data_A[g, ] - yhat_A)^2)
    RSS_B = sum((data_B[g, ] - yhat_B)^2)

    sigma2_pooled = (RSS_A + RSS_B) / (n_A + n_B - 6)

    # Wald test: H0: A_A = A_B
    # SE of amplitude difference (approximate)
    SE_diff = sqrt(sigma2_pooled * (1/n_A + 1/n_B) * 2)  # Approximate

    z_stat = abs(A_A[g] - A_B[g]) / max(SE_diff, 0.01)
    pvals[g] = 2 * pnorm(-z_stat)
  }

  return(pvals)
}


#' Test for phase difference between two conditions
#'
#' @param phi_A Phase estimates from condition A
#' @param phi_B Phase estimates from condition B
#' @param A_A Amplitude estimates from condition A (for weighting)
#' @param A_B Amplitude estimates for condition B
#'
#' @return Vector of p-values
test_phase_difference = function(phi_A, phi_B, A_A, A_B) {
  G = length(phi_A)
  pvals = numeric(G)

  for (g in 1:G) {
    # Circular difference
    delta_phi = ((phi_B[g] - phi_A[g] + 12) %% 24) - 12

    # Watson-Williams test (approximation)
    # Weight by amplitude (more confident in phase when amplitude is high)
    weight = sqrt(A_A[g]^2 + A_B[g]^2)

    # Simple z-test approximation
    # Phase SE is roughly proportional to sigma/A (smaller for larger A)
    SE_phi = 2 / max(weight, 0.1)  # Approximate

    z_stat = abs(delta_phi) / SE_phi
    pvals[g] = 2 * pnorm(-z_stat)
  }

  return(pvals)
}


#' Find optimal sample size combination for target power
#'
#' @param power_matrix Power matrix (n_A x n_B)
#' @param n_A_range n_A values
#' @param n_B_range n_B values
#' @param target Target power (default 0.8)
#'
#' @return List with optimal n_A, n_B, and achieved power
find_optimal_sample_size = function(power_matrix, n_A_range, n_B_range, target = 0.8) {

  # Find all combinations achieving target power
  idx = which(power_matrix >= target, arr.ind = TRUE)

  if (nrow(idx) == 0) {
    return(list(n_A = NA, n_B = NA, power = NA, achieved = FALSE))
  }

  # Find minimum total sample size
  total_n = outer(n_A_range, n_B_range, "+")
  achievable = power_matrix >= target

  if (!any(achievable)) {
    return(list(n_A = NA, n_B = NA, power = NA, achieved = FALSE))
  }

  # Set non-achievable to Inf
  total_n[!achievable] = Inf

  # Find minimum
  min_idx = which(total_n == min(total_n), arr.ind = TRUE)[1, ]

  return(list(
    n_A = n_A_range[min_idx[1]],
    n_B = n_B_range[min_idx[2]],
    power = power_matrix[min_idx[1], min_idx[2]],
    total_n = n_A_range[min_idx[1]] + n_B_range[min_idx[2]],
    achieved = TRUE
  ))
}


#'=============================================================================
#' STRATIFIED POWER ANALYSIS
#'=============================================================================
#' Functions for analyzing power stratified by signal-to-noise ratio (r)
#' Stratified power analysis for circadian rhythms

#' Calculate Stratified Power with TD, FD, and FDC
#'
#' @param strat_power Array [sample_size, r_stratum, sim] of power values
#' @param strat_TD Array [sample_size, r_stratum, sim] of true discoveries
#' @param strat_FD Array [sample_size, r_stratum, sim] of false discoveries
#' @param strat_n_targets Array [sample_size, r_stratum, sim] of target counts
#' @param sample_sizes Vector of sample sizes
#' @param strata_labels Labels for r strata
#' @param idx_n Index of reference sample size (e.g., n=60)
#'
#' @return List with stratified power statistics
calculateStratifiedPower <- function(strat_power, strat_TD, strat_FD,
                                     strat_n_targets,
                                     sample_sizes, strata_labels,
                                     idx_n = NULL) {

  n_size <- length(sample_sizes)
  n_r_strata <- length(strata_labels)

  # Use middle sample size if not specified
  if (is.null(idx_n)) {
    idx_n <- ceiling(n_size / 2)
  }

  #---------------------------------------------------------------------------
  # Calculate mean power by stratum at reference sample size
  #---------------------------------------------------------------------------
  mean_power_n60 <- apply(strat_power[idx_n, , ], 1, mean, na.rm = TRUE)
  se_power_n60 <- apply(strat_power[idx_n, , ], 1, sd, na.rm = TRUE) /
    apply(strat_power[idx_n, , ], 1, function(x) sum(!is.na(x)))
  ci_power_n60 <- 1.96 * se_power_n60

  #---------------------------------------------------------------------------
  # Calculate True Discoveries by stratum at reference sample size
  #---------------------------------------------------------------------------
  mean_TD_n60 <- apply(strat_TD[idx_n, , ], 1, mean, na.rm = TRUE)
  se_TD_n60 <- apply(strat_TD[idx_n, , ], 1, sd, na.rm = TRUE) /
    apply(strat_TD[idx_n, , ], 1, function(x) sum(!is.na(x)))
  ci_TD_n60 <- 1.96 * se_TD_n60

  #---------------------------------------------------------------------------
  # Calculate False Discoveries by stratum at reference sample size
  #---------------------------------------------------------------------------
  mean_FD_n60 <- apply(strat_FD[idx_n, , ], 1, mean, na.rm = TRUE)
  se_FD_n60 <- apply(strat_FD[idx_n, , ], 1, sd, na.rm = TRUE) /
    apply(strat_FD[idx_n, , ], 1, function(x) sum(!is.na(x)))

  #---------------------------------------------------------------------------
  # Calculate FDC by stratum at reference sample size
  #---------------------------------------------------------------------------
  # FDC = FD / TD (false discoveries per true discovery)
  fdc_n60 <- mean_FD_n60 / mean_TD_n60
  fdc_n60[is.infinite(fdc_n60) | is.nan(fdc_n60)] <- NA

  #---------------------------------------------------------------------------
  # Marginal power across all sample sizes
  #---------------------------------------------------------------------------
  marginal_power_by_n <- numeric(n_size)
  marginal_TD_by_n <- numeric(n_size)
  marginal_FD_by_n <- numeric(n_size)
  marginal_FDC_by_n <- numeric(n_size)

  for (j in seq_along(sample_sizes)) {
    total_TD <- sum(strat_TD[j, , ], na.rm = TRUE)
    total_FD <- sum(strat_FD[j, , ], na.rm = TRUE)
    total_targets <- sum(strat_n_targets[j, , ], na.rm = TRUE)

    marginal_power_by_n[j] <- if (total_targets > 0) total_TD / total_targets else NA
    marginal_TD_by_n[j] <- total_TD
    marginal_FD_by_n[j] <- total_FD
    marginal_FDC_by_n[j] <- if (total_TD > 0) total_FD / total_TD else NA
  }

  return(list(
    mean_power = mean_power_n60,
    se_power = se_power_n60,
    ci_power = ci_power_n60,
    mean_TD = mean_TD_n60,
    se_TD = se_TD_n60,
    ci_TD = ci_TD_n60,
    mean_FD = mean_FD_n60,
    se_FD = se_FD_n60,
    fdc = fdc_n60,
    marginal_power = marginal_power_by_n,
    marginal_TD = marginal_TD_by_n,
    marginal_FD = marginal_FD_by_n,
    marginal_FDC = marginal_FDC_by_n,
    idx_n = idx_n,
    n_size = n_size,
    n_r_strata = n_r_strata
  ))
}

#' Plot Stratified True Discoveries by r
#'
#' @param power_results Output from calculateStratifiedPower()
#' @param strata_labels Labels for r strata
#' @param title Plot title
#' @param color Bar color
plotStratifiedTD <- function(power_results, strata_labels,
                             title = "True Discoveries by Signal-to-Noise",
                             color = "darkgreen") {

  n_r_strata <- power_results$n_r_strata
  mean_TD <- power_results$mean_TD
  ci_TD <- power_results$ci_TD

  plot(1:n_r_strata, mean_TD,
       type = "b", pch = 19, col = color, lwd = 2,
       xlim = c(0.5, n_r_strata + 0.5),
       ylim = c(0, max(mean_TD + ci_TD, na.rm = TRUE) * 1.1),
       xlab = "r = A/sigma (Signal-to-Noise Ratio)",
       ylab = "Mean # True Discoveries",
       main = title,
       xaxt = "n")
  axis(1, at = 1:n_r_strata, labels = strata_labels, las = 2, cex.axis = 0.6)
  grid()

  # Add error bars
  arrows(1:n_r_strata, mean_TD - ci_TD,
         1:n_r_strata, mean_TD + ci_TD,
         angle = 90, code = 3, length = 0.02, col = color, lwd = 1.5)
}

#' Plot Stratified False Discovery Cost by r
#'
#' @param power_results Output from calculateStratifiedPower()
#' @param strata_labels Labels for r strata
#' @param title Plot title
#' @param color Bar color
plotStratifiedFDC <- function(power_results, strata_labels,
                              title = "False Discovery Cost by Signal-to-Noise",
                              color = "darkred") {

  n_r_strata <- power_results$n_r_strata
  fdc <- power_results$fdc

  plot(1:n_r_strata, fdc,
       type = "b", pch = 19, col = color, lwd = 2,
       xlim = c(0.5, n_r_strata + 0.5),
       ylim = c(0, min(max(fdc, na.rm = TRUE) * 1.2, 5)),  # Cap at 5 for readability
       xlab = "r = A/sigma (Signal-to-Noise Ratio)",
       ylab = "False Discovery Cost (FD/TD)",
       main = title,
       xaxt = "n")
  axis(1, at = 1:n_r_strata, labels = strata_labels, las = 2, cex.axis = 0.6)

  # Add reference lines
  abline(h = 0.1, lty = 2, col = "green")  # Excellent: 1 FP per 10 TP
  abline(h = 0.5, lty = 2, col = "orange") # Borderline: 1 FP per 2 TP
  abline(h = 1.0, lty = 2, col = "red")    # Bad: equal FP and TP

  grid()
  legend("topright",
         legend = c("FDC = 0.1 (Excellent)", "FDC = 0.5 (Borderline)", "FDC = 1.0 (Bad)"),
         col = c("green", "orange", "red"), lty = 2, cex = 0.7)
}

# plotAllStratifiedPower removed — use plotSingleCohortPower() (plot_single_cohort.R)
# or plotDiffPower() (plot_diff.R), both of which compute FDR curves from raw p-values.
if (FALSE) {  # kept for reference only — never called
plotAllStratifiedPower <- function(strat_power, strat_TD, strat_FD, strat_n_targets,
                                   strata_labels, sample_sizes,
                                   test_name = "DR",
                                   output_file = NULL) {

  pdf(output_file, width = 14, height = 10)
  par(mfrow = c(2, 2), mai = c(1.0, 1.0, 0.6, 0.3), mgp = c(3, 0.5, 0))

  colors <- list(
    DR = "darkgreen",
    DP = "darkblue"
  )
  col <- colors[[test_name]]

  n_size <- length(sample_sizes)
  n_r_strata <- length(strata_labels)

  # Use rainbow colors for different sample sizes
  size_colors <- rainbow(n_size, s = 0.6, v = 0.8)

  #---------------------------------------------------------------------------
  # Panel A: Stratified Power by r (multiple lines for each n)
  #---------------------------------------------------------------------------
  mean_power <- apply(strat_power, c(1, 2), mean, na.rm = TRUE)
  matplot(1:n_r_strata, t(100 * mean_power),
          type = "l", lwd = 2,
          col = size_colors, lty = 1,
          xlim = c(0.5, n_r_strata + 0.5), ylim = c(0, 100),
          xlab = "r = A/sigma (Signal-to-Noise Ratio)",
          ylab = "Targeted Power (%)",
          main = paste0(test_name, " Power by r"),
          xaxt = "n")
  axis(1, at = 1:n_r_strata, labels = strata_labels, las = 2, cex.axis = 0.6)
  abline(h = 80, lty = 2, col = "gray")
  grid()
  legend("bottomright", paste0("n=", sample_sizes),
         col = size_colors, lty = 1, lwd = 2, cex = 0.7)
  mtext("A", side = 3, at = -0.02, font = 2)

  #---------------------------------------------------------------------------
  # Panel B: True Discoveries by r (multiple lines for each n)
  #---------------------------------------------------------------------------
  mean_TD <- apply(strat_TD, c(1, 2), mean, na.rm = TRUE)
  matplot(1:n_r_strata, t(mean_TD),
          type = "l", lwd = 2,
          col = size_colors, lty = 1,
          xlim = c(0.5, n_r_strata + 0.5),
          ylim = c(0, max(mean_TD, na.rm = TRUE) * 1.1),
          xlab = "r = A/sigma (Signal-to-Noise Ratio)",
          ylab = "Mean # True Discoveries",
          main = paste0(test_name, " True Discoveries by r"),
          xaxt = "n")
  axis(1, at = 1:n_r_strata, labels = strata_labels, las = 2, cex.axis = 0.6)
  grid()
  mtext("B", side = 3, at = -0.02, font = 2)

  #---------------------------------------------------------------------------
  # Panel C: Marginal Power vs Sample Size at different FDR thresholds
  #---------------------------------------------------------------------------
  # FDR thresholds to test
  fdr_thresholds <- c(0.01, 0.05, 0.10, 0.20)
  fdr_colors <- c("darkgreen", "steelblue", "orange", "red")
  fdr_labels <- c("FDR 1%", "FDR 5%", "FDR 10%", "FDR 20%")

  # Calculate marginal power for each threshold
  # Note: strat_TD was calculated at FDR < 0.05
  # For other thresholds, we approximate (power scales with threshold)
  marginal_power_05 <- numeric(n_size)
  for (j in 1:n_size) {
    total_TD <- sum(strat_TD[j, , ], na.rm = TRUE)
    total_targets <- sum(strat_n_targets[j, , ], na.rm = TRUE)
    marginal_power_05[j] <- if (total_targets > 0) total_TD / total_targets else NA
  }

  # Approximate power at other thresholds
  marginal_power_by_threshold <- matrix(NA, nrow = n_size, ncol = 4)
  marginal_power_by_threshold[, 2] <- marginal_power_05  # FDR 5%
  marginal_power_by_threshold[, 1] <- pmin(1, marginal_power_05 * 0.7)   # FDR 1%
  marginal_power_by_threshold[, 3] <- pmin(1, marginal_power_05 * 1.15) # FDR 10%
  marginal_power_by_threshold[, 4] <- pmin(1, marginal_power_05 * 1.25) # FDR 20%

  matplot(sample_sizes, 100 * marginal_power_by_threshold,
          type = "b", pch = 19, lwd = 2,
          col = fdr_colors, lty = 1,
          xlim = c(0, max(sample_sizes) * 1.1), ylim = c(0, 100),
          xlab = "Sample Size (per group)",
          ylab = "Marginal Power (%)",
          main = paste0(test_name, " Power vs Sample Size by FDR"))
  abline(h = 80, lty = 2, col = "gray")
  grid()
  legend("bottomright", fdr_labels, col = fdr_colors, lty = 1, pch = 19, lwd = 2, cex = 0.7)
  mtext("C", side = 3, at = -0.02, font = 2)

  # Note: Thresholds 1%, 10%, 20% are approximations based on 5% FDR results

  #---------------------------------------------------------------------------
  # Panel D: Total Discoveries vs Sample Size by FDR threshold
  #---------------------------------------------------------------------------
  total_TD_05 <- numeric(n_size)
  for (j in 1:n_size) {
    total_TD_05[j] <- sum(strat_TD[j, , ], na.rm = TRUE)
  }

  # Approximate TD at other thresholds (scales with power)
  total_TD_by_threshold <- matrix(NA, nrow = n_size, ncol = 4)
  total_TD_by_threshold[, 2] <- total_TD_05  # FDR 5%
  total_TD_by_threshold[, 1] <- total_TD_05 * 0.7   # FDR 1%
  total_TD_by_threshold[, 3] <- total_TD_05 * 1.15 # FDR 10%
  total_TD_by_threshold[, 4] <- total_TD_05 * 1.25 # FDR 20%

  matplot(sample_sizes, total_TD_by_threshold,
          type = "b", pch = 19, lwd = 2,
          col = fdr_colors, lty = 1,
          xlim = c(0, max(sample_sizes) * 1.1),
          ylim = c(0, max(total_TD_by_threshold, na.rm = TRUE) * 1.1),
          xlab = "Sample Size (per group)",
          ylab = "Total True Discoveries",
          main = paste0(test_name, " Discoveries vs Sample Size by FDR"))
  grid()
  mtext("D", side = 3, at = -0.02, font = 2)

  dev.off()
  cat("Figure saved:", output_file, "\n")
}
}  # end if(FALSE) — plotAllStratifiedPower reference block

#'=============================================================================
#' SAMPLING DENSITY VS SAMPLE SIZE ANALYSIS
#'=============================================================================
#' Tradeoff: more time points per subject vs more subjects

#' Run Sampling Density vs Sample Size Simulation
#'
#' @param n_subjects Sample sizes for subjects per group to test
#' @param n_time_points Number of time points per subject to test
#' @param total_samples Fixed total number of measurements (optional)
#' @param test_type "DR" or "DP"
#' @param nsims Number of simulations per scenario
#' @param prop_DR Proportion of DR genes (for DR test)
#' @param prop_DP Proportion of DP genes (for DP test)
#' @param phase_diff Phase shift range (for DP test)
#' @param amp_diff Unused; retained for interface compatibility
#' @param times_base Base time distribution (will be subsampled)
#'
#' @return List with simulation results
runSimsDensity <- function(n_subjects = c(10, 20, 30, 60),
                           n_time_points = c(2, 4, 6, 12),
                           total_samples = NULL,
                           test_type = "DR",
                           nsims = 50,
                           prop_DR = 0.15,
                           prop_DP = 0.15,
                           phase_diff = c(0, 0),
                           amp_diff = c(0.5, 2),
                           times_base = NULL) {

  cat("====================================================================\n")
  cat("SAMPLING DENSITY VS SAMPLE SIZE ANALYSIS\n")
  cat("====================================================================\n\n")

  cat(sprintf("Test type: %s\n", test_type))
  cat(sprintf("Subject sizes: %s\n", paste(n_subjects, collapse = ", ")))
  cat(sprintf("Time points: %s\n", paste(n_time_points, collapse = ", ")))
  cat(sprintf("Simulations: %d\n\n", nsims))

  if (is.null(times_base)) {
    stop("'times_base' must be provided: a numeric vector of sample collection times (hours) ",
         "from your pilot dataset. Example: times_base = c(2.1, 6.4, 14.0, ...)")
  }

  n_n <- length(n_subjects)
  n_nt <- length(n_time_points)

  # Storage: [n_subjects, n_time_points, nsims]
  power <- array(NA, dim = c(n_n, n_nt, nsims))
  TD <- array(NA, dim = c(n_n, n_nt, nsims))
  FD <- array(NA, dim = c(n_n, n_nt, nsims))

  for (i in seq_along(n_subjects)) {
    n <- n_subjects[i]

    for (j in seq_along(n_time_points)) {
      nt <- n_time_points[j]

      cat(sprintf(">>> n=%d subjects, %d time points\n", n, nt))

      # Sample nt time points from base times (with replacement if needed)
      times_idx <- sample(1:length(times_base), nt, replace = nt > length(times_base))
      times_sub <- times_base[times_idx]

      # Run simulation with these times
      sim_out <- runSimsDiff(
        sample_sizes = c(n),
        nsims = nsims,
        ngenes = 5000,
        prop_rhythmic = 0.30,
        prop_DR = if(test_type == "DR") prop_DR else 0.00,
        prop_DP = if(test_type == "DP") prop_DP else 0.00,
        phase_diff = phase_diff,
        amp_diff = amp_diff,
        cts = times_sub,
        design = "passive",
        test_types = c(test_type),
        verbose = FALSE
      )

      # Get results
      if (test_type == "DR") {
        fdr <- sim_out$fdr_DR[, 1, ]
      } else if (test_type == "DP") {
        fdr <- sim_out$fdr_DP[, 1, ]
      } else {
        fdr <- sim_out$fdr_DM[, 1, ]
      }

      for (s in 1:nsims) {
        diff_type <- sim_out$diff_type[[s]]

        # Define targets based on test type
        if (test_type == "DR") {
          is_target <- diff_type %in% c(2, 3)
          effectsize_DR1 <- sim_out$effectsize[[s]]$DR1
          effectsize_DR2 <- sim_out$effectsize[[s]]$DR2
          effectsize_R2 <- abs(effectsize_DR1 - effectsize_DR2)
          is_target <- is_target & effectsize_R2 >= 0.1
        } else if (test_type == "DP") {
          is_target <- diff_type == 4
          effectsize_phase <- sim_out$effectsize[[s]]$phase
          is_target <- is_target & effectsize_phase >= 6  # 6 hour phase shift
        } else if (test_type == "DM") {
          is_target <- diff_type == 5
        }

        discoveries <- fdr[, s] <= 0.05

        TD[i, j, s] <- sum(discoveries & is_target, na.rm = TRUE)
        FD[i, j, s] <- sum(discoveries & !is_target, na.rm = TRUE)
        power[i, j, s] <- if (sum(is_target) > 0) TD[i, j, s] / sum(is_target) else NA
      }

      mean_p <- mean(power[i, j, ], na.rm = TRUE)
      cat(sprintf("  Power: %.1f%%\n", 100 * mean_p))
    }
    cat("\n")
  }

  # Calculate summary statistics
  power_avg <- apply(power, c(1, 2), mean, na.rm = TRUE)
  TD_avg <- apply(TD, c(1, 2), mean, na.rm = TRUE)
  FD_avg <- apply(FD, c(1, 2), mean, na.rm = TRUE)
  FDC_avg <- FD_avg / TD_avg

  results <- list(
    n_subjects = n_subjects,
    n_time_points = n_time_points,
    power = power,
    TD = TD,
    FD = FD,
    power_avg = power_avg,
    TD_avg = TD_avg,
    FD_avg = FD_avg,
    FDC_avg = FDC_avg,
    test_type = test_type,
    nsims = nsims
  )

  return(results)
}

#' Plot Sampling Density Results
#'
#' @param density_results Output from runSimsDensity()
#' @param output_file Path to save PDF
plotDensityResults <- function(density_results, output_file = NULL) {

  n_subjects <- density_results$n_subjects
  n_time_points <- density_results$n_time_points
  power_avg <- density_results$power_avg
  TD_avg <- density_results$TD_avg
  FDC_avg <- density_results$FDC_avg
  test_type <- density_results$test_type

  pdf(output_file, width = 14, height = 10)
  par(mfrow = c(2, 2), mai = c(0.9, 1.0, 0.6, 0.3), mgp = c(3, 0.5, 0))

  colors <- rainbow(length(n_time_points), s = 0.6, v = 0.8)

  #---------------------------------------------------------------------------
  # Panel A: Power vs Sample Size for different time point densities
  #---------------------------------------------------------------------------
  matplot(n_subjects, 100 * power_avg,
          type = "b", pch = 19, lwd = 2,
          col = colors, lty = 1,
          xlim = c(0, max(n_subjects) * 1.1), ylim = c(0, 100),
          xlab = "Number of Subjects (per group)",
          ylab = "Power (%)",
          main = paste0(test_type, " Power: Subjects vs Time Points"))
  abline(h = 80, lty = 2, col = "gray")
  grid()
  legend("bottomright", paste0(n_time_points, " time points"),
         col = colors, lty = 1, pch = 19, lwd = 2, cex = 0.7)
  mtext("A", side = 3, at = -0.02, font = 2)

  #---------------------------------------------------------------------------
  # Panel B: True Discoveries vs Sample Size
  #---------------------------------------------------------------------------
  matplot(n_subjects, TD_avg,
          type = "b", pch = 19, lwd = 2,
          col = colors, lty = 1,
          xlim = c(0, max(n_subjects) * 1.1),
          ylim = c(0, max(TD_avg) * 1.1),
          xlab = "Number of Subjects (per group)",
          ylab = "Mean # True Discoveries",
          main = paste0(test_type, " Discoveries: Subjects vs Time Points"))
  grid()
  mtext("B", side = 3, at = -0.02, font = 2)

  #---------------------------------------------------------------------------
  # Panel C: Power Heatmap
  #---------------------------------------------------------------------------
  par(mai = c(0.9, 1.0, 0.6, 0.5))
  z_matrix <- 100 * t(power_avg)
  image(z_matrix,
        xlab = "Number of Subjects", ylab = "Time Points",
        main = paste0(test_type, " Power Heatmap (%)"),
        col = colorRampPalette(c("white", "yellow", "orange", "red", "darkred"))(100),
        xaxt = "n", yaxt = "n")
  axis(1, at = seq(0, 1, length.out = length(n_subjects)), labels = n_subjects)
  axis(2, at = seq(0, 1, length.out = length(n_time_points)), labels = n_time_points)
  box()
  mtext("C", side = 3, at = -0.02, font = 2)

  #---------------------------------------------------------------------------
  # Panel D: FDC Heatmap
  #---------------------------------------------------------------------------
  par(mai = c(0.9, 1.0, 0.6, 0.5))
  z_matrix2 <- t(FDC_avg)
  z_matrix2[is.na(z_matrix2) | is.infinite(z_matrix2)] <- 5
  image(z_matrix2,
        xlab = "Number of Subjects", ylab = "Time Points",
        main = paste0(test_type, " FDC Heatmap"),
        col = colorRampPalette(c("white", "yellow", "orange", "red", "darkred"))(100),
        xaxt = "n", yaxt = "n")
  axis(1, at = seq(0, 1, length.out = length(n_subjects)), labels = n_subjects)
  axis(2, at = seq(0, 1, length.out = length(n_time_points)), labels = n_time_points)
  box()
  mtext("D", side = 3, at = -0.02, font = 2)

  dev.off()
  cat("Figure saved:", output_file, "\n")
}

#' Print Density Analysis Summary Table
#'
#' @param density_results Output from runSimsDensity()
printDensitySummary <- function(density_results) {

  n_subjects <- density_results$n_subjects
  n_time_points <- density_results$n_time_points
  power_avg <- density_results$power_avg
  TD_avg <- density_results$TD_avg
  FD_avg <- density_results$FD_avg
  FDC_avg <- density_results$FDC_avg

  cat("\n====================================================================\n")
  cat(sprintf("%s POWER: SUBJECTS VS TIME POINTS\n", density_results$test_type))
  cat("====================================================================\n\n")

  for (j in seq_along(n_time_points)) {
    nt <- n_time_points[j]
    cat(sprintf("Time points: %d\n", nt))
    cat(sprintf("%-10s | %-8s | %-8s | %-8s | %-8s\n",
                "Subjects", "Power", "TD", "FD", "FDC"))
    cat(paste0(rep("-", 50), collapse = ""), "\n")

    for (i in seq_along(n_subjects)) {
      cat(sprintf("%-10d | %-8.1f%% | %-8.1f | %-8.1f | %-8.3f\n",
                  n_subjects[i],
                  100 * power_avg[i, j],
                  TD_avg[i, j],
                  FD_avg[i, j],
                  FDC_avg[i, j]))
    }
    cat("\n")
  }
}

#' Compare multiple test statistics for two-group comparison
#'
#' @description Compare different approaches to testing rhythm differences
#'
#' @param sim_data Simulated two-group data
#' @param methods Test methods to compare
#'
#' @return Comparison results
compare_two_group_methods = function(sim_data, methods = c("amplitude", "phase", "joint")) {

  G = nrow(sim_data$data_A)
  ground_truth = sim_data$ground_truth$has_diff

  results = list()

  # Fit cosinor to each condition
  fit_A = t(apply(sim_data$data_A, 1, function(y) {
    fit = one_cosinor_OLS(sim_data$times_A, y, period = 24)
    c(A = fit$A, phi = fit$phi, pvalue = fit$pvalue)
  }))

  fit_B = t(apply(sim_data$data_B, 1, function(y) {
    fit = one_cosinor_OLS(sim_data$times_B, y, period = 24)
    c(A = fit$A, phi = fit$phi, pvalue = fit$pvalue)
  }))

  for (method in methods) {
    if (method == "amplitude") {
      pvals = test_amplitude_difference(
        fit_A[, "A"], fit_B[, "A"],
        sim_data$data_A, sim_data$data_B,
        sim_data$times_A, sim_data$times_B
      )
    } else if (method == "phase") {
      pvals = test_phase_difference(
        fit_A[, "phi"], fit_B[, "phi"],
        fit_A[, "A"], fit_B[, "A"]
      )
    } else if (method == "joint") {
      pvals_amp = test_amplitude_difference(
        fit_A[, "A"], fit_B[, "A"],
        sim_data$data_A, sim_data$data_B,
        sim_data$times_A, sim_data$times_B
      )
      pvals_phase = test_phase_difference(
        fit_A[, "phi"], fit_B[, "phi"],
        fit_A[, "A"], fit_B[, "A"]
      )
      pvals = pchisq(-2 * (log(pvals_amp) + log(pvals_phase)), df = 4, lower.tail = FALSE)
    }

    qvals = p.adjust(pvals, method = "BH")
    discovered = qvals < 0.05

    TP = sum(discovered & ground_truth)
    FP = sum(discovered & !ground_truth)
    FN = sum(!discovered & ground_truth)
    TN = sum(!discovered & !ground_truth)

    results[[method]] = list(
      power = TP / max(TP + FN, 1),
      fdr = FP / max(TP + FP, 1),
      sensitivity = TP / max(TP + FN, 1),
      specificity = TN / max(TN + FP, 1)
    )
  }

  return(results)
}

#'=============================================================================
#' POWER CALCULATION AT MULTIPLE THRESHOLDS
#'=============================================================================
#' Functions for calculating exact power at different FDR/p-value thresholds
#' from saved q-values and p-values

#' Calculate Power at Multiple Thresholds from Raw q-values
#'
#' @param qvalues Array [genes, nsims] or [sample_size, genes, nsims] of FDR values
#' @param is_target Logical vector or list indicating which genes are targets
#' @param is_null Logical vector or list indicating which genes are null
#' @param thresholds Vector of FDR thresholds to test
#' @param sample_sizes Vector of sample sizes (if qvalues is 3D array)
#'
#' @return List with power, TD, FD at each threshold
calculatePowerByThreshold <- function(qvalues, is_target, is_null,
                                      thresholds = c(0.01, 0.05, 0.10, 0.20),
                                      sample_sizes = NULL) {

  n_thresholds <- length(thresholds)

  # Check if qvalues is 3D [sample_size, genes, nsims] or 2D [genes, nsims]
  if (is.null(dim(qvalues)) || length(dim(qvalues)) == 2) {
    # 2D case: single sample size
    n_sims <- ncol(qvalues)
    n_genes <- nrow(qvalues)

    # Handle is_target and is_null as lists (one per simulation) or vectors
    if (is.list(is_target)) {
      # Calculate per simulation
      power <- TD <- FD <- matrix(NA, nrow = n_thresholds, ncol = n_sims)

      for (s in 1:n_sims) {
        target_vec <- if(is.matrix(is_target[[s]])) is_target[[s]] else is_target[[s]]
        null_vec <- if(is.matrix(is_null[[s]])) is_null[[s]] else is_null[[s]]

        for (t in 1:n_thresholds) {
          discoveries <- qvalues[, s] <= thresholds[t]
          TD[t, s] <- sum(discoveries & target_vec, na.rm = TRUE)
          FD[t, s] <- sum(discoveries & null_vec, na.rm = TRUE)
          n_targets <- sum(target_vec, na.rm = TRUE)
          power[t, s] <- if (n_targets > 0) TD[t, s] / n_targets else NA
        }
      }
    } else {
      # Single vector for all simulations
      power <- TD <- FD <- numeric(n_thresholds)
      for (t in 1:n_thresholds) {
        discoveries <- qvalues <= thresholds[t]
        TD[t] <- sum(discoveries & is_target, na.rm = TRUE)
        FD[t] <- sum(discoveries & is_null, na.rm = TRUE)
        n_targets <- sum(is_target, na.rm = TRUE)
        power[t] <- if (n_targets > 0) TD[t] / n_targets else NA
      }
    }

  } else if (length(dim(qvalues)) == 3) {
    # 3D case: multiple sample sizes
    n_sizes <- dim(qvalues)[1]
    n_genes <- dim(qvalues)[2]
    n_sims <- dim(qvalues)[3]

    power <- TD <- FD <- array(NA, dim = c(n_sizes, n_thresholds, n_sims))

    for (j in 1:n_sizes) {
      for (s in 1:n_sims) {
        # Get target/null for this simulation
        target_vec <- if(is.matrix(is_target[[s]])) is_target[[s]] else is_target[[s]]
        null_vec <- if(is.matrix(is_null[[s]])) is_null[[s]] else is_null[[s]]

        for (t in 1:n_thresholds) {
          discoveries <- qvalues[j, , s] <= thresholds[t]
          TD[j, t, s] <- sum(discoveries & target_vec, na.rm = TRUE)
          FD[j, t, s] <- sum(discoveries & null_vec, na.rm = TRUE)
          n_targets <- sum(target_vec, na.rm = TRUE)
          power[j, t, s] <- if (n_targets > 0) TD[j, t, s] / n_targets else NA
        }
      }
    }
  }

  return(list(
    power = power,
    TD = TD,
    FD = FD,
    thresholds = thresholds
  ))
}

#' Calculate Stratified Power at Multiple Thresholds
#'
#' @param values Array of p-values or q-values
#' @param r_values Vector of r values per gene (signal-to-noise)
#' @param is_target Logical indicator for target genes
#' @param is_null Logical indicator for null genes
#' @param r_strata Stratification breaks for r
#' @param thresholds Vector of FDR thresholds
#' @param sample_sizes Vector of sample sizes (if 3D array)
#' @param apply_fdr Apply FDR correction (BH) to p-values before thresholding
#'
#' @return List with stratified power by r and threshold
calculateStratifiedPowerByThreshold <- function(values, r_values,
                                                is_target_list, is_null_list,
                                                r_strata, strata_labels,
                                                thresholds = c(0.01, 0.05, 0.10, 0.20),
                                                sample_sizes = NULL,
                                                apply_fdr = TRUE) {

  n_strata <- length(r_strata) - 1
  n_thresholds <- length(thresholds)

  # Check dimensions
  if (is.null(dim(values)) || length(dim(values)) == 2) {
    # Single sample size
    n_sims <- ncol(values)
    n_genes <- nrow(values)

    # Storage: [stratum, threshold, sim]
    power <- array(NA, dim = c(n_strata, n_thresholds, n_sims))
    TD <- array(NA, dim = c(n_strata, n_thresholds, n_sims))
    FD <- array(NA, dim = c(n_strata, n_thresholds, n_sims))

    for (s in 1:n_sims) {
      # Get r values for this simulation
      r_vec <- r_values[[s]]
      target_vec <- if(is.matrix(is_target_list[[s]])) is_target_list[[s]] else is_target_list[[s]]
      null_vec <- if(is.matrix(is_null_list[[s]])) is_null_list[[s]] else is_null_list[[s]]

      # Stratify by r
      xgr <- cut(r_vec, breaks = r_strata, include.lowest = TRUE, labels = FALSE)

      for (t in 1:n_thresholds) {
        # Apply FDR correction if requested (values are p-values)
        test_values <- if (apply_fdr) p.adjust(values[, s], method = "BH") else values[, s]
        discoveries <- test_values <= thresholds[t]

        for (k in 1:n_strata) {
          in_stratum <- xgr == k
          TD[k, t, s] <- sum(discoveries & target_vec & in_stratum, na.rm = TRUE)
          FD[k, t, s] <- sum(discoveries & null_vec & in_stratum, na.rm = TRUE)
          n_targets <- sum(target_vec & in_stratum, na.rm = TRUE)
          power[k, t, s] <- if (n_targets > 0) TD[k, t, s] / n_targets else NA
        }
      }
    }

  } else if (length(dim(values)) == 3) {
    # Multiple sample sizes
    n_sizes <- dim(values)[1]
    n_sims <- dim(values)[3]

    # Storage: [sample_size, stratum, threshold, sim]
    power <- array(NA, dim = c(n_sizes, n_strata, n_thresholds, n_sims))
    TD <- array(NA, dim = c(n_sizes, n_strata, n_thresholds, n_sims))
    FD <- array(NA, dim = c(n_sizes, n_strata, n_thresholds, n_sims))

    for (j in 1:n_sizes) {
      for (s in 1:n_sims) {
        r_vec <- r_values[[s]]
        target_vec <- if(is.matrix(is_target_list[[s]])) is_target_list[[s]] else is_target_list[[s]]
        null_vec <- if(is.matrix(is_null_list[[s]])) is_null_list[[s]] else is_null_list[[s]]

        xgr <- cut(r_vec, breaks = r_strata, include.lowest = TRUE, labels = FALSE)

        for (t in 1:n_thresholds) {
          # Apply FDR correction if requested (values are p-values)
          test_values <- if (apply_fdr) p.adjust(values[j, , s], method = "BH") else values[j, , s]
          discoveries <- test_values <= thresholds[t]

          for (k in 1:n_strata) {
            in_stratum <- xgr == k
            TD[j, k, t, s] <- sum(discoveries & target_vec & in_stratum, na.rm = TRUE)
            FD[j, k, t, s] <- sum(discoveries & null_vec & in_stratum, na.rm = TRUE)
            n_targets <- sum(target_vec & in_stratum, na.rm = TRUE)
            power[j, k, t, s] <- if (n_targets > 0) TD[j, k, t, s] / n_targets else NA
          }
        }
      }
    }
  }

  return(list(
    power = power,
    TD = TD,
    FD = FD,
    r_strata = r_strata,
    strata_labels = strata_labels,
    thresholds = thresholds,
    sample_sizes = sample_sizes
  ))
}

#' Plot Stratified Power by r for Different FDR Thresholds
#'
#' @param power_by_threshold Output from calculateStratifiedPowerByThreshold
#' @param strata_labels Labels for r strata
#' @param sample_sizes Sample sizes (if multiple)
#' @param threshold_labels Labels for thresholds
#' @param ref_n Reference sample size to plot (if multiple)
#' @param test_name "DR" or "DP"
#' @param output_file Path to save PDF
plotPowerByThreshold <- function(power_by_threshold, strata_labels,
                                  sample_sizes = NULL,
                                  threshold_labels = c("FDR 1%", "FDR 5%", "FDR 10%", "FDR 20%"),
                                  ref_n = 60,
                                  test_name = "DR",
                                  output_file = NULL) {

  pdf(output_file, width = 10, height = 8)
  par(mfrow = c(1, 1), mai = c(1.2, 1.2, 0.8, 0.4), mgp = c(3, 0.5, 0))

  # Extract power at reference sample size if multiple
  if (!is.null(dim(power_by_threshold$power)) && length(dim(power_by_threshold$power)) == 4) {
    # [sample_size, stratum, threshold, sim]
    idx_n <- which.min(abs(sample_sizes - ref_n))
    power_mat <- apply(power_by_threshold$power[idx_n, , , ], c(1, 2), mean, na.rm = TRUE)
  } else {
    # [stratum, threshold, sim]
    power_mat <- apply(power_by_threshold$power, c(1, 2), mean, na.rm = TRUE)
  }

  n_strata <- nrow(power_mat)
  n_thresholds <- ncol(power_mat)

  # Colors for thresholds
  threshold_colors <- c("darkgreen", "steelblue", "orange", "red")[1:n_thresholds]

  matplot(1:n_strata, 100 * power_mat,
          type = "b", pch = 19, lwd = 2,
          col = threshold_colors, lty = 1,
          xlim = c(0.5, n_strata + 0.5), ylim = c(0, 100),
          xlab = "r = A/sigma (Signal-to-Noise Ratio)",
          ylab = "Targeted Power (%)",
          main = paste0(test_name, " Power by r at Different FDR Thresholds (n=", ref_n, ")"),
          xaxt = "n")
  axis(1, at = 1:n_strata, labels = strata_labels, las = 2, cex.axis = 0.8)
  abline(h = 80, lty = 2, col = "gray")
  grid()
  legend("bottomright", threshold_labels[1:n_thresholds],
         col = threshold_colors, lty = 1, pch = 19, lwd = 2, cex = 0.8)

  dev.off()
  cat("Figure saved:", output_file, "\n")
}

#' Plot Marginal Power vs Sample Size for Different Thresholds
#'
#' @param power_by_threshold Output from calculateStratifiedPowerByThreshold
#' @param sample_sizes Vector of sample sizes
#' @param threshold_labels Labels for thresholds
#' @param test_name "DR" or "DP"
#' @param output_file Path to save PDF
plotMarginalPowerByThreshold <- function(power_by_threshold, sample_sizes,
                                           threshold_labels = c("FDR 1%", "FDR 5%", "FDR 10%", "FDR 20%"),
                                           test_name = "DR",
                                           output_file = NULL) {

  pdf(output_file, width = 10, height = 8)
  par(mfrow = c(1, 1), mai = c(1.2, 1.2, 0.8, 0.4), mgp = c(3, 0.5, 0))

  # power is [sample_size, stratum, threshold, sim]
  # Marginal power = total TD / total targets across all strata
  n_sizes <- length(sample_sizes)
  n_thresholds <- length(power_by_threshold$thresholds)

  marginal_power <- matrix(NA, nrow = n_sizes, ncol = n_thresholds)

  for (j in 1:n_sizes) {
    for (t in 1:n_thresholds) {
      total_TD <- sum(power_by_threshold$TD[j, , t, ], na.rm = TRUE)
      # Get total targets (independent of threshold)
      total_targets <- sum(power_by_threshold$TD[j, , t, ] / power_by_threshold$power[j, , t, ], na.rm = TRUE)
      # Alternative: use first threshold to get target count
      total_targets <- sum(power_by_threshold$TD[j, , 1, ], na.rm = TRUE) /
        mean(power_by_threshold$power[j, , 1, ], na.rm = TRUE)
      marginal_power[j, t] <- total_TD / total_targets
    }
  }

  threshold_colors <- c("darkgreen", "steelblue", "orange", "red")[1:n_thresholds]

  matplot(sample_sizes, 100 * marginal_power,
          type = "b", pch = 19, lwd = 2,
          col = threshold_colors, lty = 1,
          xlim = c(0, max(sample_sizes) * 1.1), ylim = c(0, 100),
          xlab = "Sample Size (per group)",
          ylab = "Marginal Power (%)",
          main = paste0(test_name, " Power vs Sample Size by FDR Threshold"),
          xaxt = "n")
  axis(1, at = sample_sizes, labels = sample_sizes)
  abline(h = 80, lty = 2, col = "gray")
  grid()
  legend("bottomright", threshold_labels[1:n_thresholds],
         col = threshold_colors, lty = 1, pch = 19, lwd = 2, cex = 0.8)

  dev.off()
  cat("Figure saved:", output_file, "\n")
}
