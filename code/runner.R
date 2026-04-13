#' Run Multiple Simulations for Differential Power Analysis
#'
#' Test power to detect differential rhythmicity, phase, or amplitude
#' between two groups using rigorous DCP pipeline.
#'
#' @param sample_sizes Vector of sample sizes (per group)
#' @param nsims Number of simulations per sample size
#' @param ngenes Number of genes
#' @param prop_rhythmic Proportion of rhythmic genes
#' @param prop_DR Proportion with differential rhythmicity
#' @param prop_DP Proportion with differential phase
#' @param prop_DA Proportion with differential amplitude
#' @param phase_diff Range of phase shift for DP genes (c(min, max))
#' @param amp_diff Range of amplitude ratio for DA genes (c(min, max))
#' @param design "active" or "passive"
#' @param cts TOD distribution for passive design
#' @param test_types Which tests to run ("DR", "DP", "DA", "all")
#' @param verbose Print progress
#'
#' @return List with p-values, FDR, and ground truth for each simulation
#'
#' @details This function uses the full DiffCircadian pipeline (DCP_Rhythmicity,
#' DCP_DiffR2, DCP_DiffPar) to ensure rigorous likelihood ratio tests with
#' proper hierarchical testing and Sidak adjustments.

# Detection functions (DCP pipeline) should be sourced before this file
# e.g., source("code/detection.R") in calling script

runSimsDiff <- function(sample_sizes = c(12, 24, 36),
                        nsims = 50,
                        ngenes = 5000,
                        prop_rhythmic = 0.25,
                        prop_DR = 0.1,
                        prop_DP = 0.1,
                        prop_DA = 0.1,
                        phase_diff = NULL,
                        amp_diff = NULL,
                        design = c("active", "passive"),
                        cts = NULL,
                        test_types = c("DR", "DP", "DM"),
                        verbose = TRUE,
                        harmonics = NULL,
                        # New config-object arguments (used when first arg is CircadianBioOptions)
                        design.opts = NULL,
                        analysis.opts = NULL) {

  # --- Dual-signature: detect config objects vs legacy flat args ---
  if (inherits(sample_sizes, "CircadianBioOptions")) {
    bio.opts <- sample_sizes
    if (is.null(design.opts) || !inherits(design.opts, "CircadianDesignOptions")) {
      # Second positional arg might be design.opts
      if (inherits(nsims, "CircadianDesignOptions")) {
        design.opts <- nsims
      } else {
        stop("When using CircadianBioOptions, must also provide CircadianDesignOptions")
      }
    }
    if (is.null(analysis.opts)) {
      if (inherits(ngenes, "CircadianAnalysisOptions")) {
        analysis.opts <- ngenes
      } else {
        analysis.opts <- CircadianAnalysisOptions()
      }
    }

    # Extract from config objects
    sample_sizes  <- design.opts$sample_sizes
    nsims         <- design.opts$nsims
    ngenes        <- bio.opts$ngenes
    prop_rhythmic <- bio.opts$prop_rhythmic
    prop_DR       <- bio.opts$prop_DR
    prop_DP       <- bio.opts$prop_DP
    prop_DA       <- bio.opts$prop_DA
    prop_DM       <- bio.opts$prop_DM %||% 0
    phase_diff    <- bio.opts$phase_diff
    amp_diff      <- bio.opts$amp_diff
    mesor_diff    <- bio.opts$mesor_diff %||% c(0.5, 2.0)
    dp_shift_mode <- bio.opts$dp_shift_mode %||% "fixed"
    dr_amp_scale  <- bio.opts$dr_amp_scale %||% 1.0
    dr_sigma_scale <- bio.opts$dr_sigma_scale %||% 1.0
    design        <- design.opts$design
    cts           <- design.opts$cts
    test_types    <- design.opts$test_types
    period        <- bio.opts$period
    alpha         <- analysis.opts$alpha
    p.adjust.method <- analysis.opts$p.adjust.method
    parallel.ncores <- analysis.opts$parallel.ncores
    amp.cutoff    <- analysis.opts$amp.cutoff
    DCmethod      <- analysis.opts$DCmethod %||% "DCP"
  } else {
    # Legacy path: use hardcoded defaults for DCP params
    design         <- match.arg(design)
    period         <- 24
    alpha          <- 0.05
    p.adjust.method <- "BH"
    parallel.ncores <- 1
    amp.cutoff     <- 0
    DCmethod       <- "DCP"
    dp_shift_mode  <- "fixed"
    dr_amp_scale   <- 1.0
    dr_sigma_scale <- 1.0
  }

  # For passive design, use TOD distribution
  if (design == "passive" && is.null(cts)) {
    # Create default TOD distribution
    set.seed(123)
    cts = c(
      rnorm(30, 6, 2),    # Morning peak
      rnorm(30, 14, 3),   # Afternoon
      runif(20, 0, 24)    # Uniform background
    )
    cts = cts %% 24
  }

  # Initialize results
  pval_DR = fdr_DR = array(NA, dim = c(ngenes, length(sample_sizes), nsims))
  pval_DP = fdr_DP = array(NA, dim = c(ngenes, length(sample_sizes), nsims))
  pval_DA = fdr_DA = array(NA, dim = c(ngenes, length(sample_sizes), nsims))
  pval_DM = fdr_DM = array(NA, dim = c(ngenes, length(sample_sizes), nsims))

  # Ground truth storage
  diff_type_list = vector("list", nsims)
  effectsize_list = vector("list", nsims)

  # Gene names
  gene_names = paste0("Gene", 1:ngenes)

  # Loop over simulations
  for (i in 1:nsims) {
    if (verbose && i %% 10 == 0) {
      cat("Simulation", i, "of", nsims, "\n")
    }

    # Loop over sample sizes
    for (j in seq_along(sample_sizes)) {
      n = sample_sizes[j]

      # For active design, cts is the ZT template (length B). Expand to length n
      # by repeating the template so each ZT gets floor(n/B) or ceil(n/B) replicates.
      cts_n <- if (design == "active" && !is.null(cts) && length(cts) != n) {
        sort(rep_len(cts, n))
      } else {
        cts
      }

      # Simulate two-group data
      # Build call args; include empirical params only when using config objects
      sim_args <- list(
        ngenes = ngenes,
        n1 = n,
        n2 = n,
        prop_rhythmic = prop_rhythmic,
        prop_DR = prop_DR,
        prop_DP = prop_DP,
        prop_DA = prop_DA,
        prop_DM = prop_DM,
        phase_diff = phase_diff,
        dp_shift_mode = dp_shift_mode,
        amp_diff = amp_diff,
        mesor_diff = mesor_diff,
        design = design,
        cts = cts_n,
        sim.seed = i * 1000,
        harmonics = harmonics
      )
      # Pass empirical distributions when available (config-object path)
      if (exists("bio.opts", inherits = FALSE) && !is.null(bio.opts)) {
    if (!is.null(bio.opts$lBaselineExpr)) sim_args$lBaselineExpr <- bio.opts$lBaselineExpr
    if (!is.null(bio.opts$lBaselineExpr2)) sim_args$lBaselineExpr2 <- bio.opts$lBaselineExpr2
    if (!is.null(bio.opts$lOD)) {
      sim_args$lOD <- bio.opts$lOD + log(dr_sigma_scale)
    }
    if (!is.null(bio.opts$lOD2)) {
      sim_args$lOD2 <- bio.opts$lOD2 + log(dr_sigma_scale)
    }
    if (!is.null(bio.opts$amplitude)) {
      sim_args$amplitude <- bio.opts$amplitude * dr_amp_scale
    }
    if (!is.null(bio.opts$amplitude2)) {
      sim_args$amplitude2 <- bio.opts$amplitude2 * dr_amp_scale
    }
    if (!is.null(bio.opts$sigma_rhythmic)) {
      sim_args$sigma_rhythmic <- bio.opts$sigma_rhythmic * dr_sigma_scale
    }
    if (!is.null(bio.opts$cts2)) {
      # For active design, cts2 from the pilot may differ in length from the
      # target n (the pilot had a different number of samples). Expand the
      # ZT template (or, if length matches cts_n, use directly).
      cts2_raw <- bio.opts$cts2
      cts2_n <- if (design == "active" && length(cts2_raw) != n) {
        sort(rep_len(cts2_raw, n))
      } else {
        cts2_raw
      }
      sim_args$cts2 <- cts2_n
    }
      }
      sim_data = do.call(simCircadianDiff, sim_args)

      # Store ground truth (only need to do this once per simulation)
      if (j == 1) {
        diff_type_list[[i]] = sim_data$ground_truth$diff_type
        effectsize_list[[i]] = list(
          DR1 = sim_data$effectsize_DR1,
          DR2 = sim_data$effectsize_DR2,
          phase = sim_data$effectsize_phase,
          amp = sim_data$effectsize_amp,
          mesor = sim_data$effectsize_mesor
        )
      }

      # Initialize p-values (default to 1 for untested genes)
      pval_DR_g = rep(1, ngenes)
      pval_DP_g = rep(1, ngenes)
      pval_DA_g = rep(1, ngenes)
      pval_DM_g = rep(1, ngenes)

      # =================================================================
      # DIFFERENTIAL TESTING (method dispatcher)
      # =================================================================

      if (DCmethod == "DCP") {
        # ---------------------------------------------------------------
        # DCP Pipeline: TOJR + LR tests (rigorous hierarchical approach)
        # ---------------------------------------------------------------
        tryCatch({
          x1_dcp = format_for_DCP(sim_data$expr1, sim_data$times1, gene_names)
          x2_dcp = format_for_DCP(sim_data$expr2, sim_data$times2, gene_names)

          # Step 1: TOJR classification with hierarchical Sidak testing
          rhythm_res = DCP_Rhythmicity(
            x1 = x1_dcp, x2 = x2_dcp,
            method = "Sidak_FS",
            period = period, amp.cutoff = amp.cutoff,
            alpha = alpha, CI = FALSE,
            p.adjust.method = p.adjust.method,
            parallel.ncores = parallel.ncores
          )

          # Step 2: DR test (genes not "arrhy")
          if ("DR" %in% test_types) {
            n_testable_DR = sum(rhythm_res$rhythm.joint$TOJR != "arrhy")
            if (n_testable_DR > 0) {
              dr_results = DCP_DiffR2(
                rhythm_res, method = "LR", TOJR = NULL,
                alpha = alpha, p.adjust.method = p.adjust.method,
                parallel.ncores = parallel.ncores
              )
              match_idx = match(gene_names, dr_results$gname)
              pval_DR_g[!is.na(match_idx)] = dr_results$p.R2[match_idx[!is.na(match_idx)]]
            }
          }

          # Step 3: DP, DA, DM tests (only for "both" genes)
          needs_par <- any(c("DP", "DA", "DM") %in% test_types)
          if (needs_par) {
            n_testable = sum(rhythm_res$rhythm.joint$TOJR == "both")
            if (n_testable > 0) {
              # Use A&phase&M to get all three parameter tests
              Par_mode <- if ("DM" %in% test_types) "A&phase&M" else "A&phase"
              dp_da_results = DCP_DiffPar(
                rhythm_res, Par = Par_mode, TOJR = NULL,
                alpha = alpha, p.adjust.method = p.adjust.method,
                parallel.ncores = parallel.ncores
              )
              match_idx = match(gene_names, dp_da_results$gname)

              if ("DP" %in% test_types) {
                pval_DP_g[!is.na(match_idx)] = dp_da_results$p.delta.peak[match_idx[!is.na(match_idx)]]
              }
              if ("DA" %in% test_types) {
                pval_DA_g[!is.na(match_idx)] = dp_da_results$p.delta.A[match_idx[!is.na(match_idx)]]
              }
              if ("DM" %in% test_types && "p.delta.M" %in% colnames(dp_da_results)) {
                pval_DM_g[!is.na(match_idx)] = dp_da_results$p.delta.M[match_idx[!is.na(match_idx)]]
              }
            }
          }

        }, error = function(e) {
          if (verbose) {
            warning(sprintf("DCP pipeline failed for sim %d, n=%d: %s",
                           i, n, e$message))
          }
        })

      } else if (DCmethod == "CircaCompare") {
        # ---------------------------------------------------------------
        # CircaCompare: per-gene NLS with Wald t-tests for DA, DP, DM
        # No DR test available (CircaCompare has no delta-R2 equivalent)
        # ---------------------------------------------------------------
        tryCatch({
          cc_result = detect_CircaCompare(
            expr1 = sim_data$expr1, times1 = sim_data$times1,
            expr2 = sim_data$expr2, times2 = sim_data$times2,
            gene_names = gene_names, period = period
          )
          pval_DP_g = cc_result$pval_DP
          pval_DA_g = cc_result$pval_DA
          pval_DM_g = cc_result$pval_DM
          # pval_DR_g stays at 1 (CircaCompare has no DR test)
        }, error = function(e) {
          if (verbose) {
            warning(sprintf("CircaCompare failed for sim %d, n=%d: %s",
                           i, n, e$message))
          }
        })
      }

      # Adjust for multiple testing (BH correction)
      fdr_DR_g = pval_DR_g
      fdr_DP_g = pval_DP_g
      fdr_DA_g = pval_DA_g
      fdr_DM_g = pval_DM_g

      # Only adjust p-values that were actually tested (< 1)
      ix.tested = pval_DR_g < 1
      if (sum(ix.tested) > 0) {
        fdr_DR_g[ix.tested] = p.adjust(pval_DR_g[ix.tested], method = p.adjust.method)
      }

      ix.tested = pval_DP_g < 1
      if (sum(ix.tested) > 0) {
        fdr_DP_g[ix.tested] = p.adjust(pval_DP_g[ix.tested], method = p.adjust.method)
      }

      ix.tested = pval_DA_g < 1
      if (sum(ix.tested) > 0) {
        fdr_DA_g[ix.tested] = p.adjust(pval_DA_g[ix.tested], method = p.adjust.method)
      }

      ix.tested = pval_DM_g < 1
      if (sum(ix.tested) > 0) {
        fdr_DM_g[ix.tested] = p.adjust(pval_DM_g[ix.tested], method = p.adjust.method)
      }

      # Store results
      pval_DR[, j, i] = pval_DR_g
      fdr_DR[, j, i] = fdr_DR_g
      pval_DP[, j, i] = pval_DP_g
      fdr_DP[, j, i] = fdr_DP_g
      pval_DA[, j, i] = pval_DA_g
      fdr_DA[, j, i] = fdr_DA_g
      pval_DM[, j, i] = pval_DM_g
      fdr_DM[, j, i] = fdr_DM_g
    }
  }

  return(list(
    pval_DR = pval_DR,
    fdr_DR = fdr_DR,
    pval_DP = pval_DP,
    fdr_DP = fdr_DP,
    pval_DA = pval_DA,
    fdr_DA = fdr_DA,
    pval_DM = pval_DM,
    fdr_DM = fdr_DM,
    diff_type = diff_type_list,
    effectsize = effectsize_list,
    sample_sizes = sample_sizes,
    nsims = nsims,
    ngenes = ngenes,
    sim_params = list(
      prop_rhythmic = prop_rhythmic,
      prop_DR = prop_DR,
      prop_DP = prop_DP,
      prop_DA = prop_DA,
      design = design,
      DCmethod = DCmethod
    )
  ))
}


# =====================================================================
# Stratified Power Analysis (DR or DP)
# =====================================================================
#' Run power analysis with r = A/sigma stratification.
#'
#' Loops over sample sizes, calls runSimsDiff() for each, and computes
#' stratified power / TD / FD at FDR 5%.
#'
#' @param bio.opts      CircadianBioOptions
#' @param design.opts   CircadianDesignOptions
#' @param analysis.opts CircadianAnalysisOptions
#' @param test_type     "DR" or "DP"
#'
#' @return List compatible with plotWithSE()
runPowerAnalysis <- function(bio.opts, design.opts, analysis.opts,
                             test_type = "DR", verbose = TRUE) {

  sample_sizes  <- design.opts$sample_sizes
  nsims         <- design.opts$nsims
  ngenes        <- bio.opts$ngenes
  target_effect <- analysis.opts$target_effect
  r_strata      <- analysis.opts$r_strata
  strata_labels <- analysis.opts$strata_labels
  n_r_strata    <- length(r_strata) - 1

  # Storage: raw p-values [sample_size, gene, sim]
  pvalues <- array(NA, dim = c(length(sample_sizes), ngenes, nsims))

  # Nested ground truth [[sample_size]][[sim]]
  r_values_list  <- vector("list", length(sample_sizes))
  is_target_list <- vector("list", length(sample_sizes))
  is_null_list   <- vector("list", length(sample_sizes))
  for (.j in seq_along(sample_sizes)) {
    r_values_list[[.j]]  <- vector("list", nsims)
    is_target_list[[.j]] <- vector("list", nsims)
    is_null_list[[.j]]   <- vector("list", nsims)
  }

  # Stratified power at FDR 5% [sample_size, r_stratum, sim]
  strat_power     <- array(NA, dim = c(length(sample_sizes), n_r_strata, nsims))
  strat_TD        <- array(NA, dim = c(length(sample_sizes), n_r_strata, nsims))
  strat_FD        <- array(NA, dim = c(length(sample_sizes), n_r_strata, nsims))
  strat_n_targets <- array(NA, dim = c(length(sample_sizes), n_r_strata, nsims))
  strat_n_nulls   <- array(NA, dim = c(length(sample_sizes), n_r_strata, nsims))
  strat_FDR       <- array(NA, dim = c(length(sample_sizes), n_r_strata, nsims))
  strat_alpha     <- array(NA, dim = c(length(sample_sizes), n_r_strata, nsims))
  strat_FDC       <- array(NA, dim = c(length(sample_sizes), n_r_strata, nsims))

  # Marginal (unstratified) quantities [sample_size, sim]
  marginal_power <- matrix(NA, nrow = length(sample_sizes), ncol = nsims)
  marginal_FDR   <- matrix(NA, nrow = length(sample_sizes), ncol = nsims)
  marginal_alpha <- matrix(NA, nrow = length(sample_sizes), ncol = nsims)
  marginal_TD    <- matrix(NA, nrow = length(sample_sizes), ncol = nsims)
  marginal_FD    <- matrix(NA, nrow = length(sample_sizes), ncol = nsims)

  for (j in seq_along(sample_sizes)) {
    n <- sample_sizes[j]
    if (verbose) cat(sprintf("  >>> n = %d\n", n))

    # Create single-n design for this iteration
    iter_design <- CircadianDesignOptions(
      sample_sizes = c(n),
      nsims        = nsims,
      design       = design.opts$design,
      cts          = design.opts$cts,
      test_types   = c(test_type)
    )

    sim_out <- runSimsDiff(bio.opts, iter_design, analysis.opts)

    # Extract p-values (genes x nsims, sample_size dim = 1)
    pval_key <- paste0("pval_", test_type)
    fdr_key  <- paste0("fdr_", test_type)
    pval <- sim_out[[pval_key]][, 1, ]
    fdr  <- sim_out[[fdr_key]][, 1, ]

    pvalues[j, , ] <- pval

    for (i in 1:nsims) {
      diff_type       <- sim_out$diff_type[[i]]
      effectsize_DR1  <- sim_out$effectsize[[i]]$DR1
      effectsize_DR2  <- sim_out$effectsize[[i]]$DR2

      if (test_type == "DR") {
        r_values <- pmax(effectsize_DR1, effectsize_DR2)
        effectsize_R2 <- abs(effectsize_DR1 - effectsize_DR2)
        is_diff <- diff_type %in% c(2, 3)
        is_target <- is_diff & (effectsize_R2 >= target_effect)
      } else if (test_type == "DP") {
        r_values <- pmin(effectsize_DR1, effectsize_DR2)
        effectsize_phase <- sim_out$effectsize[[i]]$phase
        is_diff <- diff_type == 4
        is_target <- is_diff & (effectsize_phase >= target_effect)
      } else if (test_type == "DA") {
        r_values <- pmin(effectsize_DR1, effectsize_DR2)
        effectsize_amp <- sim_out$effectsize[[i]]$amp
        is_diff <- diff_type == 6   # DA is now type 6 (legacy)
        is_target <- is_diff & (effectsize_amp >= target_effect)
      } else if (test_type == "DM") {
        r_values <- pmin(effectsize_DR1, effectsize_DR2)
        effectsize_m <- sim_out$effectsize[[i]]$mesor
        is_diff  <- diff_type == 5   # DM is type 5
        is_target <- is_diff & (effectsize_m >= target_effect)
      }
      is_null <- !is_diff

      r_values_list[[j]][[i]]  <- r_values
      is_target_list[[j]][[i]] <- is_target
      is_null_list[[j]][[i]]   <- is_null

      Zg  <- ifelse(is_diff, 1, 0)
      Zg2 <- ifelse(is_target, 1, 0)

      r_for_strat <- r_values
      r_for_strat[!is_diff] <- 0
      xgr <- cut(r_for_strat, breaks = r_strata, include.lowest = TRUE, labels = FALSE)
      xgr[!is_diff] <- NA

      discoveries <- fdr[, i] <= 0.05

      # Stratified quantities (per r-stratum)
      for (k in 1:n_r_strata) {
        in_stratum <- xgr == k
        TD    <- sum(discoveries & Zg2 == 1 & in_stratum, na.rm = TRUE)
        FD    <- sum(discoveries & Zg == 0 & in_stratum, na.rm = TRUE)
        n_tgt <- sum(Zg2 == 1 & in_stratum, na.rm = TRUE)
        n_nul <- sum(Zg == 0 & !is.na(xgr) & in_stratum, na.rm = TRUE)
        N_disc <- TD + FD  # total discoveries in stratum

        strat_power[j, k, i]     <- if (n_tgt > 0) TD / n_tgt else NA
        strat_TD[j, k, i]        <- TD
        strat_FD[j, k, i]        <- FD
        strat_n_targets[j, k, i] <- n_tgt
        strat_n_nulls[j, k, i]   <- n_nul
        strat_FDR[j, k, i]       <- if (N_disc > 0) FD / N_disc else NA
        strat_alpha[j, k, i]     <- if (n_nul > 0) FD / n_nul else NA
        strat_FDC[j, k, i]       <- if (TD > 0) FD / TD else NA
      }

      # Marginal quantities (pooled across ALL genes, not strata)
      total_TD_marginal   <- sum(discoveries & Zg2 == 1, na.rm = TRUE)
      total_FD_marginal   <- sum(discoveries & Zg == 0, na.rm = TRUE)
      total_tgt_marginal  <- sum(Zg2 == 1, na.rm = TRUE)
      total_null_marginal <- sum(Zg == 0, na.rm = TRUE)
      total_disc_marginal <- total_TD_marginal + total_FD_marginal
      marginal_power[j, i] <- if (total_tgt_marginal > 0) total_TD_marginal / total_tgt_marginal else NA
      marginal_FDR[j, i]   <- if (total_disc_marginal > 0) total_FD_marginal / total_disc_marginal else NA
      marginal_alpha[j, i] <- if (total_null_marginal > 0) total_FD_marginal / total_null_marginal else NA
      marginal_TD[j, i]    <- total_TD_marginal
      marginal_FD[j, i]    <- total_FD_marginal
    }

    if (verbose) {
      cat(sprintf("    n=%d: ", n))
      for (k in 1:n_r_strata) {
        mean_p <- mean(strat_power[j, k, ], na.rm = TRUE)
        if (!is.nan(mean_p)) cat(sprintf("%s=%.0f%% ", strata_labels[k], 100 * mean_p))
      }
      cat("\n")
    }
  }

  list(
    sample_sizes    = sample_sizes,
    nsims           = nsims,
    target_effect   = target_effect,
    r_strata        = r_strata,
    strata_labels   = strata_labels,
    pvalues         = pvalues,
    r_values_list   = r_values_list,
    is_target_list  = is_target_list,
    is_null_list    = is_null_list,
    # Stratified quantities [sample_size, r_stratum, sim]
    strat_power     = strat_power,
    strat_TD        = strat_TD,
    strat_FD        = strat_FD,
    strat_n_targets = strat_n_targets,
    strat_n_nulls   = strat_n_nulls,
    strat_FDR       = strat_FDR,
    strat_alpha     = strat_alpha,
    strat_FDC       = strat_FDC,
    # Marginal quantities [sample_size, sim]
    marginal_power  = marginal_power,
    marginal_FDR    = marginal_FDR,
    marginal_alpha  = marginal_alpha,
    marginal_TD     = marginal_TD,
    marginal_FD     = marginal_FD
  )
}


# =====================================================================
# Phase Shift Sensitivity Analysis
# =====================================================================
#' Sweep over phase shift magnitudes AND sample sizes.
#'
#' @param bio.opts      CircadianBioOptions (base; prop_DP/phase_diff will be overridden)
#' @param design.opts   CircadianDesignOptions
#' @param analysis.opts CircadianAnalysisOptions (phase_shifts read from here)
#' @param prop_DP       Proportion of DP genes for each scenario (default 0.15)
#' @param amp_diff      Amplitude ratio range for DP genes (default c(0.5, 2))
#'
#' @return List with 4D arrays [phase, size, stratum, sim] compatible with plotPhaseShiftWithSE()
runPhaseShiftAnalysis <- function(bio.opts, design.opts, analysis.opts,
                                  prop_DP  = 0.15,
                                  amp_diff = c(0.5, 2),
                                  mc.cores = 1L) {

  sample_sizes  <- design.opts$sample_sizes
  nsims         <- design.opts$nsims
  phase_shifts  <- analysis.opts$phase_shifts
  target_effect <- analysis.opts$target_effect
  r_strata      <- analysis.opts$r_strata
  strata_labels <- analysis.opts$strata_labels
  n_r_strata    <- length(r_strata) - 1
  n_phase       <- length(phase_shifts)
  n_size        <- length(sample_sizes)

  cat(sprintf("  runPhaseShiftAnalysis: %d phases x %d N x %d sims | mc.cores=%d\n",
              n_phase, n_size, nsims, mc.cores))

  # Flatten (phase, N) into a combo grid for parallel dispatch
  combo_grid <- expand.grid(p = seq_along(phase_shifts), j = seq_along(sample_sizes))

  run_one_phasexN <- function(ci) {
    p           <- combo_grid$p[ci]
    j           <- combo_grid$j[ci]
    phase_shift <- phase_shifts[p]
    n           <- sample_sizes[j]

    opts_bio_ps <- updateBioOptions(bio.opts,
      prop_DR    = 0.00,
      prop_DP    = prop_DP,
      prop_DA    = 0.00,
      phase_diff = c(-phase_shift, phase_shift),
      amp_diff   = amp_diff
    )
    iter_design <- CircadianDesignOptions(
      sample_sizes = c(n),
      nsims        = nsims,
      design       = design.opts$design,
      cts          = design.opts$cts,
      test_types   = c("DP")
    )

    sim_out <- tryCatch(
      runSimsDiff(opts_bio_ps, iter_design, analysis.opts),
      error = function(e) { warning(sprintf("phase=%.1f N=%d failed: %s", phase_shift, n, e$message)); NULL }
    )
    if (is.null(sim_out)) return(NULL)

    fdr_DP <- sim_out$fdr_DP[, 1, ]

    # Per-sim accumulators
    strat_power     <- array(NA, dim = c(n_r_strata, nsims))
    strat_TD        <- array(NA, dim = c(n_r_strata, nsims))
    strat_FD        <- array(NA, dim = c(n_r_strata, nsims))
    strat_n_targets <- array(NA, dim = c(n_r_strata, nsims))
    strat_n_nulls   <- array(NA, dim = c(n_r_strata, nsims))
    strat_FDR       <- array(NA, dim = c(n_r_strata, nsims))
    strat_alpha     <- array(NA, dim = c(n_r_strata, nsims))
    strat_FDC       <- array(NA, dim = c(n_r_strata, nsims))
    marginal_power  <- rep(NA, nsims)
    marginal_FDR    <- rep(NA, nsims)
    marginal_alpha  <- rep(NA, nsims)
    marginal_TD     <- rep(NA, nsims)
    marginal_FD     <- rep(NA, nsims)

    for (i in seq_len(nsims)) {
      diff_type        <- sim_out$diff_type[[i]]
      effectsize_phase <- sim_out$effectsize[[i]]$phase
      effectsize_DR1   <- sim_out$effectsize[[i]]$DR1
      effectsize_DR2   <- sim_out$effectsize[[i]]$DR2

      r_values    <- pmin(effectsize_DR1, effectsize_DR2)
      is_DP       <- diff_type == 4
      Zg          <- ifelse(is_DP, 1, 0)
      Zg2         <- ifelse(is_DP & effectsize_phase >= target_effect, 1, 0)
      r_for_strat <- r_values
      r_for_strat[!is_DP] <- 0
      xgr         <- cut(r_for_strat, breaks = r_strata, include.lowest = TRUE, labels = FALSE)
      xgr[!is_DP] <- NA
      discoveries <- fdr_DP[, i] <= 0.05

      for (k in seq_len(n_r_strata)) {
        in_stratum <- xgr == k
        TD     <- sum(discoveries & Zg2 == 1 & in_stratum, na.rm = TRUE)
        FD     <- sum(discoveries & Zg  == 0 & in_stratum, na.rm = TRUE)
        n_tgt  <- sum(Zg2 == 1 & in_stratum, na.rm = TRUE)
        n_nul  <- sum(Zg  == 0 & !is.na(xgr) & in_stratum, na.rm = TRUE)
        N_disc <- TD + FD
        strat_power[k, i]     <- if (n_tgt  > 0) TD / n_tgt  else NA
        strat_TD[k, i]        <- TD
        strat_FD[k, i]        <- FD
        strat_n_targets[k, i] <- n_tgt
        strat_n_nulls[k, i]   <- n_nul
        strat_FDR[k, i]       <- if (N_disc > 0) FD / N_disc else NA
        strat_alpha[k, i]     <- if (n_nul  > 0) FD / n_nul  else NA
        strat_FDC[k, i]       <- if (TD     > 0) FD / TD     else NA
      }

      total_TD_m   <- sum(discoveries & Zg2 == 1, na.rm = TRUE)
      total_FD_m   <- sum(discoveries & Zg  == 0, na.rm = TRUE)
      total_tgt_m  <- sum(Zg2 == 1, na.rm = TRUE)
      total_nul_m  <- sum(Zg  == 0, na.rm = TRUE)
      total_disc_m <- total_TD_m + total_FD_m
      marginal_power[i] <- if (total_tgt_m  > 0) total_TD_m  / total_tgt_m  else NA
      marginal_FDR[i]   <- if (total_disc_m > 0) total_FD_m  / total_disc_m else NA
      marginal_alpha[i] <- if (total_nul_m  > 0) total_FD_m  / total_nul_m  else NA
      marginal_TD[i]    <- total_TD_m
      marginal_FD[i]    <- total_FD_m
    }

    list(p = p, j = j,
         strat_power = strat_power, strat_TD = strat_TD, strat_FD = strat_FD,
         strat_n_targets = strat_n_targets, strat_n_nulls = strat_n_nulls,
         strat_FDR = strat_FDR, strat_alpha = strat_alpha, strat_FDC = strat_FDC,
         marginal_power = marginal_power, marginal_FDR = marginal_FDR,
         marginal_alpha = marginal_alpha, marginal_TD = marginal_TD, marginal_FD = marginal_FD)
  }

  combo_results <- parallel::mclapply(
    seq_len(nrow(combo_grid)), run_one_phasexN,
    mc.cores = mc.cores, mc.set.seed = TRUE
  )

  # Reassemble into original 4D/3D arrays
  ps_strat_power     <- array(NA, dim = c(n_phase, n_size, n_r_strata, nsims))
  ps_strat_TD        <- array(NA, dim = c(n_phase, n_size, n_r_strata, nsims))
  ps_strat_FD        <- array(NA, dim = c(n_phase, n_size, n_r_strata, nsims))
  ps_strat_n_targets <- array(NA, dim = c(n_phase, n_size, n_r_strata, nsims))
  ps_strat_n_nulls   <- array(NA, dim = c(n_phase, n_size, n_r_strata, nsims))
  ps_strat_FDR       <- array(NA, dim = c(n_phase, n_size, n_r_strata, nsims))
  ps_strat_alpha     <- array(NA, dim = c(n_phase, n_size, n_r_strata, nsims))
  ps_strat_FDC       <- array(NA, dim = c(n_phase, n_size, n_r_strata, nsims))
  ps_marginal_power  <- array(NA, dim = c(n_phase, n_size, nsims))
  ps_marginal_FDR    <- array(NA, dim = c(n_phase, n_size, nsims))
  ps_marginal_alpha  <- array(NA, dim = c(n_phase, n_size, nsims))
  ps_marginal_TD     <- array(NA, dim = c(n_phase, n_size, nsims))
  ps_marginal_FD     <- array(NA, dim = c(n_phase, n_size, nsims))

  for (ci in seq_along(combo_results)) {
    r <- combo_results[[ci]]
    if (is.null(r)) next
    p <- r$p; j <- r$j
    ps_strat_power[p, j, , ]     <- r$strat_power
    ps_strat_TD[p, j, , ]        <- r$strat_TD
    ps_strat_FD[p, j, , ]        <- r$strat_FD
    ps_strat_n_targets[p, j, , ] <- r$strat_n_targets
    ps_strat_n_nulls[p, j, , ]   <- r$strat_n_nulls
    ps_strat_FDR[p, j, , ]       <- r$strat_FDR
    ps_strat_alpha[p, j, , ]     <- r$strat_alpha
    ps_strat_FDC[p, j, , ]       <- r$strat_FDC
    ps_marginal_power[p, j, ]    <- r$marginal_power
    ps_marginal_FDR[p, j, ]      <- r$marginal_FDR
    ps_marginal_alpha[p, j, ]    <- r$marginal_alpha
    ps_marginal_TD[p, j, ]       <- r$marginal_TD
    ps_marginal_FD[p, j, ]       <- r$marginal_FD
  }

  list(
    sample_sizes    = sample_sizes,
    phase_shifts    = phase_shifts,
    nsims           = nsims,
    target_effect   = target_effect,
    r_strata        = r_strata,
    strata_labels   = strata_labels,
    # Stratified quantities [phase, size, stratum, sim]
    strat_power     = ps_strat_power,
    strat_TD        = ps_strat_TD,
    strat_FD        = ps_strat_FD,
    strat_n_targets = ps_strat_n_targets,
    strat_n_nulls   = ps_strat_n_nulls,
    strat_FDR       = ps_strat_FDR,
    strat_alpha     = ps_strat_alpha,
    strat_FDC       = ps_strat_FDC,
    # Marginal quantities [phase, size, sim]
    marginal_power  = ps_marginal_power,
    marginal_FDR    = ps_marginal_FDR,
    marginal_alpha  = ps_marginal_alpha,
    marginal_TD     = ps_marginal_TD,
    marginal_FD     = ps_marginal_FD
  )
}


# =====================================================================
# Summary Function (analogous to PROPER's summaryPower)
# =====================================================================
#' Summarize power analysis results in a compact table.
#'
#' Computes mean marginal power, type I error, FDR, TD, FD, and FDC
#' across simulation replicates for each sample size.
#'
#' @param powerOutput Output from runPowerAnalysis()
#' @param verbose     Print the table (default TRUE)
#'
#' @return Data frame with one row per sample size (invisibly)
summaryRunPower <- function(powerOutput, verbose = TRUE) {
  ss <- powerOutput$sample_sizes

  avg_power <- rowMeans(powerOutput$marginal_power, na.rm = TRUE)
  avg_FDR   <- rowMeans(powerOutput$marginal_FDR, na.rm = TRUE)
  avg_alpha <- rowMeans(powerOutput$marginal_alpha, na.rm = TRUE)
  avg_TD    <- rowMeans(powerOutput$marginal_TD, na.rm = TRUE)
  avg_FD    <- rowMeans(powerOutput$marginal_FD, na.rm = TRUE)
  avg_FDC   <- avg_FD / avg_TD
  avg_FDC[!is.finite(avg_FDC)] <- NA

  res <- data.frame(
    n           = ss,
    Power       = avg_power,
    TypeI_Error = avg_alpha,
    FDR         = avg_FDR,
    Avg_TD      = avg_TD,
    Avg_FD      = avg_FD,
    FDC         = avg_FDC,
    stringsAsFactors = FALSE
  )

  if (verbose) {
    cat(sprintf("%-6s | %-8s | %-11s | %-8s | %-8s | %-8s | %-8s\n",
                "n", "Power", "Type I Err", "FDR", "Avg TD", "Avg FD", "FDC"))
    cat(paste0(rep("-", 72), collapse = ""), "\n")
    for (i in seq_len(nrow(res))) {
      cat(sprintf("%-6d | %6.1f%% | %9.4f | %6.4f | %7.1f | %7.1f | %6.2f\n",
                  res$n[i],
                  100 * res$Power[i],
                  res$TypeI_Error[i],
                  res$FDR[i],
                  res$Avg_TD[i],
                  res$Avg_FD[i],
                  ifelse(is.na(res$FDC[i]), 0, res$FDC[i])))
    }
  }

  invisible(res)
}


# =====================================================================
# Simulation-Based Single-Cohort Power
# =====================================================================
#' Simulation-Based Power for Single-Cohort Rhythmicity Detection
#'
#' Extends the closed-form CircaPower approach by running full simulations:
#' for each sample size and simulation replicate, generates data from the
#' empirical pilot parameter distributions, applies the cosinor F-test per
#' gene, applies BH correction, and aggregates empirical power and FDR.
#'
#' This is the simulation counterpart to \code{runPowerAnalysis()} for the
#' single-cohort (one-group) scenario.  It does not require a closed-form
#' solution and therefore works for any design (active or passive) and any
#' pilot parameter distribution.
#'
#' @param bio.opts      \code{CircadianBioOptions} — pilot parameter distributions.
#'   Use \code{estCircadianParam()} to build from real data.
#' @param design.opts   \code{CircadianDesignOptions} — sample sizes, nsims, design.
#' @param analysis.opts \code{CircadianAnalysisOptions} — alpha, p.adjust.method.
#' @param verbose       Print progress (default TRUE).
#'
#' @return List with:
#'   \item{marginal_power}{Matrix [sample_sizes x nsims] — per-simulation power}
#'   \item{marginal_FDR}{Matrix [sample_sizes x nsims] — per-simulation FDR}
#'   \item{sample_sizes}{Vector of sample sizes}
#'   \item{nsims}{Number of simulation replicates}
#'
#' @examples
#' \dontrun{
#' bio  <- estCircadianParam(expr, times)
#' dopt <- CircadianDesignOptions(sample_sizes = c(20, 40, 80), nsims = 50)
#' aopt <- CircadianAnalysisOptions(alpha = 0.05)
#' res  <- runSimsSingleCohort(bio, dopt, aopt)
#' rowMeans(res$marginal_power)  # mean power at each N
#' }
#'
#' @export
runSimsSingleCohort <- function(bio.opts, design.opts, analysis.opts,
                                verbose = TRUE) {

  stopifnot(inherits(bio.opts, "CircadianBioOptions"))
  stopifnot(inherits(design.opts, "CircadianDesignOptions"))
  if (missing(analysis.opts)) analysis.opts <- CircadianAnalysisOptions()

  sample_sizes    <- design.opts$sample_sizes
  nsims           <- design.opts$nsims
  design          <- design.opts$design
  cts             <- design.opts$cts
  alpha           <- analysis.opts$alpha
  p.adjust.method <- analysis.opts$p.adjust.method

  ngenes        <- bio.opts$ngenes
  prop_rhythmic <- bio.opts$prop_rhythmic
  period        <- bio.opts$period

  # Storage
  marginal_power <- matrix(NA_real_, nrow = length(sample_sizes), ncol = nsims)
  marginal_FDR   <- matrix(NA_real_, nrow = length(sample_sizes), ncol = nsims)
  marginal_TD    <- matrix(NA_real_, nrow = length(sample_sizes), ncol = nsims)
  marginal_FD    <- matrix(NA_real_, nrow = length(sample_sizes), ncol = nsims)

  # Build simOptions structure expected by simCircadian()
  sim_base <- list(
    ngenes        = ngenes,
    prop_rhythmic = prop_rhythmic,
    period        = period,
    lBaselineExpr = bio.opts$lBaselineExpr,
    lOD           = bio.opts$lOD,
    amplitude     = bio.opts$amplitude,
    phase         = bio.opts$phase,
    sim.seed      = bio.opts$sim.seed %||% 42L
  )

  for (j in seq_along(sample_sizes)) {
    n <- sample_sizes[j]
    if (verbose) cat(sprintf("  >>> n = %d\n", n))

    cts_n <- if (design == "active" && !is.null(cts) && length(cts) != n) {
      sort(rep_len(cts, n))
    } else {
      cts
    }

    for (i in seq_len(nsims)) {
      # Re-seed per (j, i) so results are reproducible and independent
      sim_opts_i <- sim_base
      sim_opts_i$sim.seed <- (j * 1000L + i) * 7919L
      set.seed(sim_opts_i$sim.seed)

      # Assign rhythmic gene indices and draw parameters
      n_rhythmic <- round(ngenes * prop_rhythmic)
      rhythmic_id <- sample(ngenes, n_rhythmic)
      is_rhythmic <- rep(FALSE, ngenes)
      is_rhythmic[rhythmic_id] <- TRUE

      mesor_g  <- bio.opts$lBaselineExpr
      sigma_g  <- exp(bio.opts$lOD)

      amp_g   <- rep(0, ngenes)
      phase_g <- rep(0, ngenes)
      if (n_rhythmic > 0) {
        # Joint sampling of A and sigma if sigma_rhythmic is available
        if (!is.null(bio.opts$sigma_rhythmic) &&
            length(bio.opts$sigma_rhythmic) == length(bio.opts$amplitude)) {
          joint_idx <- sample(length(bio.opts$amplitude), n_rhythmic, replace = TRUE)
          amp_g[rhythmic_id]   <- pmax(bio.opts$amplitude[joint_idx], 0.05)
          sigma_g[rhythmic_id] <- pmax(bio.opts$sigma_rhythmic[joint_idx], 1e-6)
        } else {
          amp_g[rhythmic_id] <- pmax(
            sample(bio.opts$amplitude, n_rhythmic, replace = TRUE), 0.05)
        }
        phase_g[rhythmic_id] <- sample(bio.opts$phase, n_rhythmic, replace = TRUE)
      }

      # Generate time points
      if (design == "active") {
        if (!is.null(cts_n)) {
          times_i <- cts_n
        } else {
          times_i <- seq(0, period, length.out = n + 1L)[seq_len(n)]
        }
      } else {
        times_i <- sampleTimesFromDist(n, cts_n)
      }

      # Simulate expression matrix [ngenes x n]
      omega <- 2 * pi / period
      expr  <- matrix(NA_real_, nrow = ngenes, ncol = n)
      for (g in seq_len(ngenes)) {
        mu <- mesor_g[g] + amp_g[g] * cos(omega * times_i - omega * phase_g[g])
        expr[g, ] <- rnorm(n, mu, sigma_g[g])
      }

      # Fit cosinor F-test per gene and compute p-values
      pvals <- tryCatch({
        fitCosinorAll(expr, times_i, period = period)$pvalue
      }, error = function(e) {
        rep(NA_real_, ngenes)
      })

      # BH correction
      pvals[is.na(pvals)] <- 1
      fdr_g <- p.adjust(pvals, method = p.adjust.method)

      discoveries  <- fdr_g <= alpha
      TD <- sum(discoveries &  is_rhythmic, na.rm = TRUE)
      FD <- sum(discoveries & !is_rhythmic, na.rm = TRUE)
      n_tgt <- sum(is_rhythmic)
      n_disc <- TD + FD

      marginal_power[j, i] <- if (n_tgt  > 0) TD / n_tgt else NA_real_
      marginal_FDR[j, i]   <- if (n_disc > 0) FD / n_disc else NA_real_
      marginal_TD[j, i]    <- TD
      marginal_FD[j, i]    <- FD
    }

    if (verbose) {
      cat(sprintf("    n=%d: mean power = %.1f%%  mean FDR = %.3f\n",
                  n,
                  100 * mean(marginal_power[j, ], na.rm = TRUE),
                  mean(marginal_FDR[j, ], na.rm = TRUE)))
    }
  }

  list(
    marginal_power = marginal_power,
    marginal_FDR   = marginal_FDR,
    marginal_TD    = marginal_TD,
    marginal_FD    = marginal_FD,
    sample_sizes   = sample_sizes,
    nsims          = nsims
  )
}
