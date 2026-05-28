#' Low-level differential power simulation engine
#'
#' Preferred interface: \code{runDifferentialPower(bio.opts, design.opts, analysis.opts)}.
#' This function supports both the config-object API (primary) and a legacy
#' flat-argument API (retained for back-compatibility).
#'
#' @param sample_sizes \code{CircadianBioOptions} (primary) or numeric vector of
#'   per-group sample sizes (legacy flat API).
#' @param nsims \code{CircadianDesignOptions} (primary) or number of simulations (legacy).
#' @param ngenes \code{CircadianAnalysisOptions} (primary) or number of genes (legacy).
#' @param prop_rhythmic Overall proportion of rhythmic genes (legacy only).
#' @param prop_DR Proportion with differential rhythmicity (legacy only).
#' @param prop_DP Proportion with differential phase (legacy only).
#' @param phase_diff Phase-shift range in hours for DP genes (legacy only).
#' @param amp_diff Unused; retained for interface compatibility.
#' @param design \code{"active"} or \code{"passive"} (legacy only).
#' @param cts TOD distribution for passive design (legacy only).
#' @param test_types Which endpoints to evaluate: any of \code{"DR"}, \code{"DP"}, \code{"DM"}.
#' @param verbose Print progress.
#' @param mc.cores Parallel cores.
#' @param design.opts \code{CircadianDesignOptions} (used with legacy positional first arg).
#' @param analysis.opts \code{CircadianAnalysisOptions}.
#'
#' @return List (class \code{SCPDiffResult}) with fields:
#'   \code{pval_DR}, \code{fdr_DR}, \code{pval_DP}, \code{fdr_DP},
#'   \code{pval_DM}, \code{fdr_DM} -- arrays \code{[ngenes x n_sizes x nsims]};
#'   \code{diff_type}, \code{effectsize}, \code{sample_sizes}, \code{nsims}.
#'   Pass directly to \code{plotDiffPower()} or \code{npower()}.
#'
#' @seealso \code{\link{runDifferentialPower}} (user-facing wrapper)

# Detection functions (DCP pipeline) should be sourced before this file
# e.g., source("code/detection.R") in calling script

#' @export
runSimsDiff <- function(sample_sizes = c(12, 24, 36),
                        nsims = 50,
                        ngenes = 5000,
                        prop_rhythmic = 0.25,
                        prop_DR = 0.1,
                        prop_DP = 0.1,
                        phase_diff = NULL,
                        amp_diff = NULL,
                        design = c("active", "passive"),
                        cts = NULL,
                        test_types = c("DR", "DP", "DM"),
                        verbose = TRUE,
                        harmonics = NULL,
                        mc.cores = 1L,
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
    harmonics     <- design.opts$harmonics %||% NULL
    fmm_omega     <- design.opts$omega %||% 1.0   # 1 = cosinor; <1 = FMM
    fmm_beta      <- design.opts$beta  %||% pi
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
    prop_DM        <- 0
    mesor_diff     <- c(0.5, 2.0)
  }

  # For passive design, use TOD distribution
  if (design == "passive" && is.null(cts)) {
    warning("passive design requested but no 'cts' provided; using a synthetic TOD distribution. ",
            "Supply cts from your pilot data via estCircadianParam()$cts for accurate power estimates.")
    cts = c(
      rnorm(30, 6, 2),
      rnorm(30, 14, 3),
      runif(20, 0, 24)
    )
    cts = cts %% 24
  }

  # Gene names
  gene_names = paste0("Gene", 1:ngenes)

  # Capture bio.opts for use inside mclapply closure
  bio_opts_inner <- if (exists("bio.opts", inherits = FALSE)) bio.opts else NULL

  # Run simulations (parallel over sims if mc.cores > 1)
  sim_results <- parallel::mclapply(seq_len(nsims), function(i) {
    pval_DR_i <- matrix(1,   nrow = ngenes, ncol = length(sample_sizes))
    pval_DP_i <- matrix(1,   nrow = ngenes, ncol = length(sample_sizes))
    pval_DM_i <- matrix(1,   nrow = ngenes, ncol = length(sample_sizes))
    fdr_DR_i  <- matrix(NA,  nrow = ngenes, ncol = length(sample_sizes))
    fdr_DP_i  <- matrix(NA,  nrow = ngenes, ncol = length(sample_sizes))
    fdr_DM_i  <- matrix(NA,  nrow = ngenes, ncol = length(sample_sizes))
    diff_type_i  <- NULL
    effectsize_i <- NULL

    for (j in seq_along(sample_sizes)) {
      n = sample_sizes[j]

      cts_n <- if (design == "active" && !is.null(cts) && length(cts) != n) {
        sort(rep_len(cts, n))
      } else {
        cts
      }

      sim_args <- list(
        ngenes = ngenes,
        n1 = n,
        n2 = n,
        prop_rhythmic = prop_rhythmic,
        prop_DR = prop_DR,
        prop_DP = prop_DP,
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

      if (!is.null(bio_opts_inner)) {
        bo <- bio_opts_inner
        if (!is.null(bo$lBaselineExpr))  sim_args$lBaselineExpr  <- bo$lBaselineExpr
        if (!is.null(bo$lBaselineExpr2)) sim_args$lBaselineExpr2 <- bo$lBaselineExpr2
        if (!is.null(bo$lOD))  sim_args$lOD  <- bo$lOD  + log(dr_sigma_scale)
        if (!is.null(bo$lOD2)) sim_args$lOD2 <- bo$lOD2 + log(dr_sigma_scale)
        if (!is.null(bo$amplitude))      sim_args$amplitude      <- bo$amplitude  * dr_amp_scale
        if (!is.null(bo$amplitude2))     sim_args$amplitude2     <- bo$amplitude2 * dr_amp_scale
        if (!is.null(bo$sigma_rhythmic)) sim_args$sigma_rhythmic <- bo$sigma_rhythmic * dr_sigma_scale
        if (!is.null(bo$cts2)) {
          cts2_raw <- bo$cts2
          sim_args$cts2 <- if (design == "active" && length(cts2_raw) != n)
                             sort(rep_len(cts2_raw, n)) else cts2_raw
        }
      }

      # FMM path: non-sinusoidal waveform when fmm_omega < 1
      if (fmm_omega < 1.0) {
        sim_data <- do.call(simCircadianDiffFMM,
                            c(sim_args, list(omega = fmm_omega, beta = fmm_beta)))
      } else {
        sim_data <- do.call(simCircadianDiff, sim_args)
      }

      if (j == 1) {
        diff_type_i  <- sim_data$ground_truth$diff_type
        effectsize_i <- list(
          DR1   = sim_data$effectsize_DR1,
          DR2   = sim_data$effectsize_DR2,
          phase = sim_data$effectsize_phase,
          amp   = sim_data$effectsize_amp,
          mesor = sim_data$effectsize_mesor
        )
      }

      pval_DR_g = rep(1, ngenes)
      pval_DP_g = rep(1, ngenes)
      pval_DM_g = rep(1, ngenes)

      if (DCmethod == "DCP") {
        tryCatch({
          x1_dcp = format_for_DCP(sim_data$expr1, sim_data$times1, gene_names)
          x2_dcp = format_for_DCP(sim_data$expr2, sim_data$times2, gene_names)

          rhythm_res = DCP_Rhythmicity(
            x1 = x1_dcp, x2 = x2_dcp,
            method = "Sidak_FS",
            period = period, amp.cutoff = amp.cutoff,
            alpha = alpha, CI = FALSE,
            p.adjust.method = p.adjust.method,
            parallel.ncores = parallel.ncores
          )

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

          needs_par <- any(c("DP", "DM") %in% test_types)
          if (needs_par) {
            n_testable = sum(rhythm_res$rhythm.joint$TOJR == "both")
            if (n_testable > 0) {
              Par_mode <- if ("DM" %in% test_types) "A&phase&M" else "A&phase"
              dp_da_results = DCP_DiffPar(
                rhythm_res, Par = Par_mode, TOJR = NULL,
                alpha = alpha, p.adjust.method = p.adjust.method,
                parallel.ncores = parallel.ncores
              )
              match_idx = match(gene_names, dp_da_results$gname)
              if ("DP" %in% test_types)
                pval_DP_g[!is.na(match_idx)] = dp_da_results$p.delta.peak[match_idx[!is.na(match_idx)]]
              if ("DM" %in% test_types && "p.delta.M" %in% colnames(dp_da_results))
                pval_DM_g[!is.na(match_idx)] = dp_da_results$p.delta.M[match_idx[!is.na(match_idx)]]
            }
          }
        }, error = function(e) {
          warning(sprintf("DCP pipeline failed for sim %d, n=%d: %s", i, n, e$message))
        })

      } else if (DCmethod == "CircaCompare") {
        tryCatch({
          cc_result = detect_CircaCompare(
            expr1 = sim_data$expr1, times1 = sim_data$times1,
            expr2 = sim_data$expr2, times2 = sim_data$times2,
            gene_names = gene_names, period = period
          )
          pval_DP_g = cc_result$pval_DP
          pval_DM_g = cc_result$pval_DM
        }, error = function(e) {
          warning(sprintf("CircaCompare failed for sim %d, n=%d: %s", i, n, e$message))
        })
      }

      # BH adjustment over all ngenes (including untested genes with pval=1 sentinel).
      # Subsetting to pval<1 before adjusting uses a smaller denominator than the true
      # hypothesis pool, inflating significance for testable genes.
      fdr_DR_g <- p.adjust(pval_DR_g, method = p.adjust.method)
      fdr_DP_g <- p.adjust(pval_DP_g, method = p.adjust.method)
      fdr_DM_g <- p.adjust(pval_DM_g, method = p.adjust.method)

      pval_DR_i[, j] = pval_DR_g; fdr_DR_i[, j] = fdr_DR_g
      pval_DP_i[, j] = pval_DP_g; fdr_DP_i[, j] = fdr_DP_g
      pval_DM_i[, j] = pval_DM_g; fdr_DM_i[, j] = fdr_DM_g
    }

    list(pval_DR = pval_DR_i, fdr_DR = fdr_DR_i,
         pval_DP = pval_DP_i, fdr_DP = fdr_DP_i,
         pval_DM = pval_DM_i, fdr_DM = fdr_DM_i,
         diff_type = diff_type_i, effectsize = effectsize_i)
  }, mc.cores = mc.cores, mc.set.seed = TRUE)

  if (verbose) cat(sprintf("Completed %d simulations.\n", nsims))

  # Combine parallel results into arrays
  pval_DR = fdr_DR = array(NA, dim = c(ngenes, length(sample_sizes), nsims))
  pval_DP = fdr_DP = array(NA, dim = c(ngenes, length(sample_sizes), nsims))
  pval_DM = fdr_DM = array(NA, dim = c(ngenes, length(sample_sizes), nsims))
  diff_type_list  <- vector("list", nsims)
  effectsize_list <- vector("list", nsims)

  for (i in seq_len(nsims)) {
    r <- sim_results[[i]]
    pval_DR[,,i] <- r$pval_DR; fdr_DR[,,i] <- r$fdr_DR
    pval_DP[,,i] <- r$pval_DP; fdr_DP[,,i] <- r$fdr_DP
    pval_DM[,,i] <- r$pval_DM; fdr_DM[,,i] <- r$fdr_DM
    diff_type_list[[i]]  <- r$diff_type
    effectsize_list[[i]] <- r$effectsize
  }

  return(list(
    pval_DR = pval_DR,
    fdr_DR = fdr_DR,
    pval_DP = pval_DP,
    fdr_DP = fdr_DP,
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

  .Deprecated("runDifferentialPower",
              msg = "runPowerAnalysis() is deprecated. Use runDifferentialPower() instead.")

  sample_sizes  <- design.opts$sample_sizes
  nsims         <- design.opts$nsims
  ngenes        <- bio.opts$ngenes
  alpha         <- analysis.opts$alpha
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

      discoveries <- fdr[, i] <= alpha

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
  alpha         <- analysis.opts$alpha
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
      discoveries <- fdr_DP[, i] <= alpha

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
#' Simulation-based alternative to the closed-form \code{CircaPower} formula:
#' works for any design (active or passive) and any pilot parameter distribution.
#'
#' @param bio.opts      \code{CircadianBioOptions} -- pilot parameter distributions.
#'   Build with \code{estCircadianParam()}.
#' @param design.opts   \code{CircadianDesignOptions} -- sample sizes, nsims, design.
#' @param analysis.opts \code{CircadianAnalysisOptions} -- alpha, p.adjust.method.
#' @param method        Detection method: \code{"DCP"} (1-harmonic cosinor F-test),
#'   \code{"JTK"}, \code{"RAIN"}, \code{"MH"} (multi-harmonic), or \code{"FMM"}
#'   (K-harmonic F-test motivated by the FMM Fourier expansion; default detector
#'   for non-cosinor signals).
#' @param K             Integer. Number of harmonics for \code{method = "FMM"}
#'   (default 2). K=1 is equivalent to DCP; K=2 captures the 1st and 2nd
#'   Fourier harmonics; K=3+ for very sharp peaks.
#' @param harmonics     Numeric length-2: \code{c(alpha2, alpha3)} harmonic coefficients (default \code{c(0,0)}).
#' @param verbose       Print progress (default TRUE).
#' @param mc.cores      Parallel cores (default 1).
#'
#' @return List with:
#'   \item{marginal_power}{Matrix [sample_sizes x nsims]}
#'   \item{marginal_FDR}{Matrix [sample_sizes x nsims]}
#'   \item{marginal_TD}{Matrix [sample_sizes x nsims]}
#'   \item{marginal_FD}{Matrix [sample_sizes x nsims]}
#'   \item{marginal_alpha}{Matrix [sample_sizes x nsims] -- empirical type-I error}
#'   \item{pvalues}{Array [sample_sizes x ngenes x nsims] -- raw p-values}
#'   \item{r_values_list}{List[[sample_size]][[sim]] -- per-gene A/sigma}
#'   \item{strat_power}{Array [sample_sizes x n_strata x nsims]}
#'   \item{strat_TD}{Array [sample_sizes x n_strata x nsims]}
#'   \item{strat_FD}{Array [sample_sizes x n_strata x nsims]}
#'   \item{strat_n_targets}{Array [sample_sizes x n_strata x nsims]}
#'   \item{strat_n_nulls}{Array [sample_sizes x n_strata x nsims]}
#'   \item{strat_FDR}{Array [sample_sizes x n_strata x nsims]}
#'   \item{strat_alpha}{Array [sample_sizes x n_strata x nsims]}
#'   \item{r_strata}{Stratum breakpoints (from analysis.opts)}
#'   \item{strata_labels}{Stratum labels}
#'   \item{n0_circapower}{CircaPower n80 estimate from pilot median r}
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
                                method = "cosinor",
                                K = 1L,
                                harmonics = c(0, 0),
                                verbose = TRUE, mc.cores = 1L) {
  # New default: method = "cosinor" + K = 1. Users wanting the two-harmonic
  # test set K = 2 only; method = "DCP" / "FMM" / "Kharmonic" are retained
  # as legacy aliases.

  stopifnot(inherits(bio.opts, "CircadianBioOptions"))
  stopifnot(inherits(design.opts, "CircadianDesignOptions"))
  if (missing(analysis.opts)) analysis.opts <- CircadianAnalysisOptions()

  method    <- match.arg(method, c("DCP", "cosinor", "Kharmonic", "JTK", "RAIN", "MH", "FMM"))
  K         <- as.integer(K)
  harmonics <- rep_len(as.numeric(harmonics), 2L)

  sample_sizes    <- design.opts$sample_sizes
  nsims           <- design.opts$nsims
  design          <- design.opts$design
  cts             <- design.opts$cts

  # Nyquist guard for the K-harmonic detector. K-harmonic regression has
  # 2K+1 free parameters and requires at least 2K+1 distinct sampling
  # phases per period to be identifiable. For active designs this is
  # decided by design.opts$cts; for passive designs the effective B is
  # the number of unique TODs (mod period). If the design fails this
  # check, warn loudly at dispatch time (before mclapply swallows it)
  # and suggest the largest K the design supports.
  if (method == "FMM" && !is.null(cts)) {
    period_check <- bio.opts$period %||% 24
    n_unique     <- length(unique(cts %% period_check))
    K_max        <- max(1L, (n_unique - 1L) %/% 2L)
    if (n_unique < 2L * K + 1L) {
      warning(sprintf(
        "runSimsSingleCohort: method=\"FMM\" with K=%d requires at least %d distinct sampling phases per period (Nyquist condition), but the design has only %d unique time(s) modulo period=%g. The K=%d test is not identifiable on this design; expect zero power. Use K = %d (the largest identifiable K for this design) or switch to method=\"DCP\".",
        K, 2L * K + 1L, n_unique, period_check, K, K_max),
        call. = FALSE)
    }
  } else if (method == "FMM" && is.null(cts)) {
    # Passive design without supplied TOD vector: sampling phases are
    # drawn from F_TOD at runtime, so the exact identifiability check is
    # deferred. We can still flag the unambiguous case where even the
    # smallest requested N is below the 2K+1 threshold; at that N the
    # K-harmonic test is rank-deficient regardless of how phases land.
    min_N <- min(sample_sizes)
    if (min_N < 2L * K + 1L) {
      warning(sprintf(
        "runSimsSingleCohort: method=\"FMM\" with K=%d, passive design (cts=NULL): the smallest requested N=%d is below the identifiability threshold 2K+1=%d, so the K=%d test will be conservative or undefined at that N. Per-replicate Nyquist identifiability is otherwise determined by the runtime TOD draws from F_TOD.",
        K, min_N, 2L * K + 1L, K),
        call. = FALSE)
    }
  }

  fmm_omega       <- design.opts$omega %||% 1.0   # 1 = cosinor; <1 = FMM
  fmm_beta        <- design.opts$beta  %||% pi
  # Per-gene FMM activation: bio.opts carries omega_dist / alpha_dist /
  # omega_rhythmic / alpha_rhythmic (set by estCircadianParamFMM or by Fig 4
  # sweep code). When any of these is present, force the FMM simulation path
  # so the per-gene parameters are honoured even though design.opts$omega == 1.
  has_fmm_per_gene <- !is.null(bio.opts$omega_dist) ||
                       !is.null(bio.opts$alpha_dist) ||
                       (!is.null(bio.opts$omega_rhythmic) &&
                        length(bio.opts$omega_rhythmic) > 0) ||
                       (!is.null(bio.opts$alpha_rhythmic) &&
                        length(bio.opts$alpha_rhythmic) > 0)
  # Two-harmonic cosinor truth: bio.opts carries per-gene (A2, phi2) vectors
  # paired with (amplitude, phase, sigma_rhythmic) by gene index. Set by
  # estCircadianParam2H(). Takes precedence over the existing alpha2_dist /
  # alpha3_dist scalar-harmonic mechanism (which we keep for Fig 4 sensitivity).
  has_2h_per_gene <- isTRUE(bio.opts$paired_2h) &&
                      !is.null(bio.opts$amplitude2) &&
                      length(bio.opts$amplitude2) > 0L &&
                      !is.null(bio.opts$phase2) &&
                      length(bio.opts$phase2) > 0L
  alpha           <- analysis.opts$alpha
  p.adjust.method <- analysis.opts$p.adjust.method
  r_strata        <- analysis.opts$r_strata
  strata_labels   <- analysis.opts$strata_labels
  n_r_strata      <- length(r_strata) - 1

  ngenes        <- bio.opts$ngenes
  prop_rhythmic <- bio.opts$prop_rhythmic
  period        <- bio.opts$period

  has_joint <- !is.null(bio.opts$sigma_rhythmic) &&
               length(bio.opts$sigma_rhythmic) == length(bio.opts$amplitude)

  # method labels: "cosinor" (canonical, K from K argument) is the recommended
  # public name; "DCP" forces single-harmonic; "Kharmonic" is an alias of
  # "cosinor"; "FMM" is retained for backward compatibility with figure scripts.
  detect_fn <- switch(method,
    DCP       = function(e, t) detect_cosinor(e, t, K = 1L, period = period),
    cosinor   = function(e, t) detect_cosinor(e, t, K = K,  period = period),
    Kharmonic = function(e, t) detect_cosinor(e, t, K = K,  period = period),
    FMM       = function(e, t) detect_cosinor(e, t, K = K,  period = period),
    JTK       = function(e, t) detect_JTK(e, t, period = period),
    RAIN      = function(e, t) detect_RAIN(e, t, period = period),
    MH        = function(e, t) detect_MH(e, t, period = period)
  )

  # CircaPower n80 estimate from median r
  n0_circapower <- circaPowerApproxN80(bio.opts, alpha = alpha)

  # Storage
  pvalues         <- array(NA_real_, dim = c(length(sample_sizes), ngenes, nsims))
  r_values_list   <- vector("list", length(sample_sizes))
  for (.j in seq_along(sample_sizes)) r_values_list[[.j]] <- vector("list", nsims)

  marginal_power  <- matrix(NA_real_, nrow = length(sample_sizes), ncol = nsims)
  marginal_FDR    <- matrix(NA_real_, nrow = length(sample_sizes), ncol = nsims)
  marginal_alpha  <- matrix(NA_real_, nrow = length(sample_sizes), ncol = nsims)
  marginal_TD     <- matrix(NA_real_, nrow = length(sample_sizes), ncol = nsims)
  marginal_FD     <- matrix(NA_real_, nrow = length(sample_sizes), ncol = nsims)

  strat_power     <- array(NA_real_, dim = c(length(sample_sizes), n_r_strata, nsims))
  strat_TD        <- array(NA_real_, dim = c(length(sample_sizes), n_r_strata, nsims))
  strat_FD        <- array(NA_real_, dim = c(length(sample_sizes), n_r_strata, nsims))
  strat_n_targets <- array(NA_real_, dim = c(length(sample_sizes), n_r_strata, nsims))
  strat_n_nulls   <- array(NA_real_, dim = c(length(sample_sizes), n_r_strata, nsims))
  strat_FDR       <- array(NA_real_, dim = c(length(sample_sizes), n_r_strata, nsims))
  strat_alpha     <- array(NA_real_, dim = c(length(sample_sizes), n_r_strata, nsims))

  for (j in seq_along(sample_sizes)) {
    n <- sample_sizes[j]
    if (verbose) cat(sprintf("  >>> n = %d\n", n))

    cts_n <- if (design == "active" && !is.null(cts) && length(cts) != n) {
      sort(rep_len(cts, n))
    } else {
      cts
    }

    # Run nsims replicates (parallel if mc.cores > 1).
    # mc.set.seed=TRUE activates L'Ecuyer-CMRG streams per worker,
    # making results reproducible across mc.cores values.
    sim_results <- parallel::mclapply(seq_len(nsims), function(i) {
      set.seed((j * 1000L + i) * 7919L)

      # Draw gene parameters
      n_rhythmic  <- round(ngenes * prop_rhythmic)
      rhythmic_id <- sample(ngenes, n_rhythmic)
      is_rhythmic <- logical(ngenes)
      is_rhythmic[rhythmic_id] <- TRUE

      mesor_g <- bio.opts$lBaselineExpr
      sigma_g <- exp(bio.opts$lOD)
      # Resample to match ngenes if pilot vectors are shorter (bio$ngenes override)
      if (length(mesor_g) != ngenes) {
        idx     <- sample(length(mesor_g), ngenes, replace = TRUE)
        mesor_g <- mesor_g[idx]
        sigma_g <- sigma_g[idx]
      }
      amp_g   <- numeric(ngenes)
      phase_g <- numeric(ngenes)

      if (n_rhythmic > 0) {
        if (has_joint) {
          ji <- sample(length(bio.opts$amplitude), n_rhythmic, replace = TRUE)
          amp_g[rhythmic_id]   <- pmax(bio.opts$amplitude[ji], 0.05)
          sigma_g[rhythmic_id] <- pmax(bio.opts$sigma_rhythmic[ji], 1e-6)
        } else {
          amp_g[rhythmic_id] <- pmax(
            sample(bio.opts$amplitude, n_rhythmic, replace = TRUE), 0.05)
        }
        phase_g[rhythmic_id] <- sample(bio.opts$phase, n_rhythmic, replace = TRUE)
      }

      # Per-gene true r = A/sigma
      r_values <- amp_g / sigma_g

      # Time points
      times_i <- if (design == "active") {
        if (!is.null(cts_n)) cts_n
        else seq(0, period, length.out = n + 1L)[seq_len(n)]
      } else {
        sampleTimesFromDist(n, cts_n)
      }

      # Simulate expression [ngenes x n]
      # Path selection (mutually exclusive, in priority order):
      #   1. FMM:  fmm_omega < 1 OR per-gene FMM params on bio.opts
      #   2. 2H:   bio.opts carries paired (A2, phi2) vectors (paired_2h = TRUE)
      #   3. cosinor (Cpp or R fallback) with optional alpha2_dist / alpha3_dist
      if (fmm_omega < 1.0 || has_fmm_per_gene) {
        fmm_out     <- simCircadianFMM(bio.opts, times_i, omega = fmm_omega,
                                       beta = fmm_beta)
        expr        <- fmm_out$expr
        is_rhythmic <- fmm_out$is_rhythmic
        r_values    <- fmm_out$r_values
      } else if (has_2h_per_gene) {
        # Two-harmonic cosinor truth.
        # Re-draw per-replicate joint (A1, phi1, A2, phi2, sigma) tuples from
        # the pilot using a SHARED index so the empirical cross-parameter
        # dependence is preserved exactly as in the pilot.
        omega_circ <- 2 * pi / period
        amp1_src   <- bio.opts$amplitude
        amp2_src   <- bio.opts$amplitude2
        phi1_src   <- bio.opts$phase
        phi2_src   <- bio.opts$phase2
        sig_src    <- bio.opts$sigma_rhythmic
        L_pilot    <- length(amp1_src)

        # Overwrite the rhythmic-gene slots with a paired draw. Non-rhythmic
        # genes already have amp_g[g] = 0 from the assignment block above.
        if (n_rhythmic > 0L) {
          ji2 <- sample.int(L_pilot, n_rhythmic, replace = TRUE)
          A1_r   <- pmax(amp1_src[ji2], 0.05)
          A2_r   <- pmax(amp2_src[ji2], 0)          # A2 can be ~0 (noise spike)
          phi1_r <- phi1_src[ji2]
          phi2_r <- phi2_src[ji2]
          sig_r  <- pmax(sig_src[ji2], 1e-6)

          amp_g[rhythmic_id]   <- A1_r
          phase_g[rhythmic_id] <- phi1_r
          sigma_g[rhythmic_id] <- sig_r
          # Update r_values with the refreshed (A1, sigma) pair.
          r_values <- amp_g / sigma_g
        }

        # Per-gene 2nd-harmonic vectors (zero for non-rhythmic genes).
        A2_g   <- numeric(ngenes)
        phi2_g <- numeric(ngenes)
        if (n_rhythmic > 0L) {
          A2_g[rhythmic_id]   <- A2_r
          phi2_g[rhythmic_id] <- phi2_r
        }

        expr <- matrix(NA_real_, nrow = ngenes, ncol = n)
        for (g in seq_len(ngenes)) {
          if (is_rhythmic[g]) {
            mu <- mesor_g[g] +
                  amp_g[g] * cos(omega_circ * (times_i - phase_g[g])) +
                  A2_g[g]  * cos(2 * omega_circ * (times_i - phi2_g[g]))
          } else {
            mu <- rep(mesor_g[g], n)
          }
          expr[g, ] <- rnorm(n, mu, sigma_g[g])
        }
      } else {
        # Traditional cosinor path (alpha2/alpha3 Fourier harmonics).
        # Per-gene alpha2/alpha3 distributions (for MH sensitivity sweeps):
        # bio.opts$alpha2_dist = list(family = "beta", a = 1, b = eta_a2) or
        #                       list(family = "fixed", value = a2)
        # When either is set, fall back to the R loop (per-gene math).
        omega_circ <- 2 * pi / period
        has_a2_dist <- !is.null(bio.opts$alpha2_dist)
        has_a3_dist <- !is.null(bio.opts$alpha3_dist)
        if (has_a2_dist) {
          spec2 <- bio.opts$alpha2_dist
          alpha2_g <- switch(spec2$family,
            beta  = rbeta(ngenes, spec2$a, spec2$b),
            fixed = rep(spec2$value, ngenes),
            stop("Unknown alpha2_dist$family: ", spec2$family))
        } else {
          alpha2_g <- rep(harmonics[1], ngenes)
        }
        if (has_a3_dist) {
          spec3 <- bio.opts$alpha3_dist
          alpha3_g <- switch(spec3$family,
            beta  = rbeta(ngenes, spec3$a, spec3$b),
            fixed = rep(spec3$value, ngenes),
            stop("Unknown alpha3_dist$family: ", spec3$family))
        } else {
          alpha3_g <- rep(harmonics[2], ngenes)
        }
        if (!has_a2_dist && !has_a3_dist &&
            exists(".CPP_LOADED", inherits = TRUE) &&
            isTRUE(get(".CPP_LOADED", inherits = TRUE)) &&
            exists("sim_cosinor_expr_fast", mode = "function")) {
          expr <- sim_cosinor_expr_fast(mesor_g, amp_g, phase_g, sigma_g, times_i,
                                        period, harmonics[1], harmonics[2])
        } else {
          expr <- matrix(NA_real_, nrow = ngenes, ncol = n)
          for (g in seq_len(ngenes)) {
            mu <- mesor_g[g] + amp_g[g] * (
              cos(omega_circ * times_i - omega_circ * phase_g[g]) +
              alpha2_g[g] * cos(2 * omega_circ * times_i - 2 * omega_circ * phase_g[g]) +
              alpha3_g[g] * cos(3 * omega_circ * times_i - 3 * omega_circ * phase_g[g])
            )
            expr[g, ] <- rnorm(n, mu, sigma_g[g])
          }
        }
      }

      pvals <- tryCatch(
        detect_fn(expr, times_i),
        error = function(e) rep(NA_real_, ngenes)
      )
      pvals[is.na(pvals)] <- 1
      fdr_g <- p.adjust(pvals, method = p.adjust.method)

      list(pvals = pvals, fdr_g = fdr_g, is_rhythmic = is_rhythmic, r_values = r_values)
    }, mc.cores = mc.cores, mc.set.seed = TRUE)

    for (i in seq_len(nsims)) {
      res_i       <- sim_results[[i]]
      pvals       <- res_i$pvals
      fdr_g       <- res_i$fdr_g
      is_rhythmic <- res_i$is_rhythmic
      r_values    <- res_i$r_values

      pvalues[j, , i]        <- pvals
      r_values_list[[j]][[i]] <- r_values

      discoveries <- fdr_g <= alpha
      n_rhythmic  <- sum(is_rhythmic)

      # r-stratum assignment for rhythmic genes only
      r_for_strat <- ifelse(is_rhythmic, r_values, 0)
      xgr <- cut(r_for_strat, breaks = r_strata, include.lowest = TRUE, labels = FALSE)
      xgr[!is_rhythmic] <- NA

      # Stratified quantities
      for (k in seq_len(n_r_strata)) {
        in_k  <- !is.na(xgr) & xgr == k
        TD_k  <- sum(discoveries &  is_rhythmic & in_k, na.rm = TRUE)
        FD_k  <- sum(discoveries & !is_rhythmic & in_k, na.rm = TRUE)  # 0 by construction
        n_tgt_k  <- sum(is_rhythmic & in_k)
        n_nul_k  <- sum(!is_rhythmic & in_k)

        strat_power[j, k, i]     <- if (n_tgt_k > 0) TD_k / n_tgt_k else NA_real_
        strat_TD[j, k, i]        <- TD_k
        strat_FD[j, k, i]        <- FD_k
        strat_n_targets[j, k, i] <- n_tgt_k
        strat_n_nulls[j, k, i]   <- n_nul_k
        disc_k <- TD_k + FD_k
        strat_FDR[j, k, i]   <- if (disc_k > 0) FD_k / disc_k else NA_real_
        strat_alpha[j, k, i] <- if (n_nul_k > 0) FD_k / n_nul_k else NA_real_
      }

      # Marginal
      TD <- sum(discoveries &  is_rhythmic, na.rm = TRUE)
      FD <- sum(discoveries & !is_rhythmic, na.rm = TRUE)
      n_null  <- ngenes - n_rhythmic
      n_disc  <- TD + FD
      marginal_power[j, i]  <- if (n_rhythmic > 0) TD / n_rhythmic else NA_real_
      marginal_FDR[j, i]    <- if (n_disc > 0) FD / n_disc else NA_real_
      marginal_alpha[j, i]  <- if (n_null > 0) FD / n_null else NA_real_
      marginal_TD[j, i]     <- TD
      marginal_FD[j, i]     <- FD
    }

    if (verbose) {
      cat(sprintf("    n=%d: mean power = %.1f%%  mean FDR = %.3f\n",
                  n,
                  100 * mean(marginal_power[j, ], na.rm = TRUE),
                  mean(marginal_FDR[j, ], na.rm = TRUE)))
    }
  }

  list(
    sample_sizes    = sample_sizes,
    nsims           = nsims,
    n0_circapower   = n0_circapower,
    r_strata        = r_strata,
    strata_labels   = strata_labels,
    pvalues         = pvalues,
    r_values_list   = r_values_list,
    marginal_power  = marginal_power,
    marginal_FDR    = marginal_FDR,
    marginal_alpha  = marginal_alpha,
    marginal_TD     = marginal_TD,
    marginal_FD     = marginal_FD,
    strat_power     = strat_power,
    strat_TD        = strat_TD,
    strat_FD        = strat_FD,
    strat_n_targets = strat_n_targets,
    strat_n_nulls   = strat_n_nulls,
    strat_FDR       = strat_FDR,
    strat_alpha     = strat_alpha
  )
}


# =====================================================================
# CircaPower grid initialisation helper
# =====================================================================
#' Estimate n for 80 percent power using CircaPower at median pilot r
#'
#' @param bio.opts \code{CircadianBioOptions} with amplitude and sigma_rhythmic.
#' @param alpha    Significance level (default 0.05).
#' @param target_power Target power (default 0.80).
#' @param n_search Integer vector of candidate sample sizes.
#' @return Single integer -- smallest n achieving target_power, or NA if not found.
#' @export
circaPowerApproxN80 <- function(bio.opts, alpha = 0.05, target_power = 0.80,
                                n_search = seq(5, 500, by = 5)) {
  if (!is.null(bio.opts$sigma_rhythmic) &&
      length(bio.opts$sigma_rhythmic) == length(bio.opts$amplitude)) {
    r_vec <- bio.opts$amplitude / bio.opts$sigma_rhythmic
  } else {
    # Fallback: assume sigma ≈ exp(median lOD)
    sigma_med <- exp(median(bio.opts$lOD, na.rm = TRUE))
    r_vec     <- bio.opts$amplitude / sigma_med
  }
  r_med <- median(r_vec, na.rm = TRUE)
  power_n <- vapply(n_search, function(n) {
    lam <- r_med^2 * n / 2
    q   <- qf(1 - alpha, 2, n - 3)
    1 - pf(q, 2, n - 3, ncp = lam)
  }, numeric(1))
  found <- which(power_n >= target_power)
  if (length(found) == 0L) return(NA_integer_)
  n_search[found[1L]]
}


# =====================================================================
# Adaptive r-strata from pilot data
# =====================================================================
#' Compute adaptive r-strata breakpoints with fixed bin width
#'
#' Generates breakpoints from 0 to ceiling(r_max / bin_width) * bin_width,
#' stepping by bin_width, then appends Inf. Every bin width is identical,
#' but the number of bins adapts to the empirical r range of the pilot data.
#'
#' @param bio.opts  \code{CircadianBioOptions} with amplitude and sigma_rhythmic.
#' @param bin_width Width of each r bin (default 0.25).
#' @param r_min_pct Lower percentile of pilot r used as range floor (default 0
#'   = always start from 0).
#' @return Numeric vector of breakpoints starting at 0 and ending at Inf.
#' @export
makeAdaptiveRStrata <- function(bio.opts, bin_width = 0.25, r_min_pct = 0) {
  if (!is.null(bio.opts$sigma_rhythmic) &&
      length(bio.opts$sigma_rhythmic) == length(bio.opts$amplitude)) {
    r_vec <- bio.opts$amplitude / bio.opts$sigma_rhythmic
  } else {
    sigma_med <- exp(median(bio.opts$lOD, na.rm = TRUE))
    r_vec     <- bio.opts$amplitude / sigma_med
  }
  r_vec  <- r_vec[is.finite(r_vec) & r_vec > 0]
  r_max  <- quantile(r_vec, 0.99, na.rm = TRUE)
  breaks <- seq(0, ceiling(r_max / bin_width) * bin_width, by = bin_width)
  c(breaks, Inf)
}


# =====================================================================
# B vs m design guidance
# =====================================================================

#' Print B vs m method guidance table
#'
#' Prints a table explaining how each detection method responds to B (time
#' points) vs m (replicates per ZT). Called automatically by recommendDesign()
#' and runSingleCohortPower() when verbose=TRUE.
#'
#' @param methods Character vector of methods to include
#' @param verbose Print to console (default TRUE)
#' @return Invisible data.frame[method, recommended_B, B_vs_m, reason]
#' @export
printMethodGuidance <- function(methods = c("DCP","JTK","RAIN","MH"),
                                 verbose = TRUE) {
  tbl <- data.frame(
    method        = c("DCP",      "JTK",   "RAIN",  "MH",
                      "LimoRhyde","DODR"),
    recommended_B = c("≥4, any",  "4–6",   "6–8",   "6 (or 3–4 if sinusoidal)",
                      "≥4, any",  "≥4, any"),
    B_vs_m        = c("N-driven", "↑m",    "↑B",    "↑B (if α₂≥0.5)",
                      "N-driven", "N-driven"),
    reason        = c(
      "NCP = N·r²/2 is B-invariant for equispaced B≥3; time points do not help",
      "Collapses replicates to per-ZT means before ranking; more m = stronger test",
      "Umbrella test uses individual observations; more ZTs improve rank resolution",
      "Adaptive K=⌊(B-1)/2⌋: K=2 at B=6 captures harmonic signal; K≥5 over-fits",
      "limma interaction model: power scales with N, not B",
      "Regression-based; power scales with N; B does not change sensitivity"
    ),
    stringsAsFactors = FALSE
  )
  tbl <- tbl[tbl$method %in% methods, ]

  if (verbose) {
    cat("\n=== B vs m Method Guidance ===\n")
    cat(sprintf("%-14s  %-26s  %-16s  %s\n",
                "Method", "Recommended B", "Preference", "Statistical reason"))
    cat(strrep("-", 105), "\n")
    for (i in seq_len(nrow(tbl))) {
      cat(sprintf("%-14s  %-26s  %-16s  %s\n",
                  tbl$method[i], tbl$recommended_B[i],
                  tbl$B_vs_m[i], tbl$reason[i]))
    }
    cat("\nKey: N-driven = B irrelevant, invest in N.\n")
    cat("     ↑B = denser time coverage helps.  ↑m = more replicates help.\n\n")
  }
  invisible(tbl)
}


#' B vs m design study -- analytical + simulation
#'
#' Full B vs m study orchestrator. Three sequential steps:
#'
#' 1. Print method guidance table (printMethodGuidance)
#' 2. Analytical: CircaPower closed-form power at each N (DCP; B-invariant).
#'    For JTK/RAIN/MH no closed form exists -- simulation is required.
#' 3. Simulation: run runSingleCohortPower() or runDifferentialPower(), or
#'    absorb a prior_result from a previous call to skip re-running.
#'
#' @param bio.opts      CircadianBioOptions from estCircadianParam() or
#'                      estCircadianParamTwoGroup()
#' @param design.opts   CircadianDesignOptions with sample_sizes and B_values
#' @param analysis.opts CircadianAnalysisOptions
#' @param methods       Detection methods to compare. Single-cohort:
#'                      c("DCP","JTK","RAIN","MH"). Differential:
#'                      c("DCP","CircaCompare","LimoRhyde","DODR").
#' @param target_power  Target power for recommendation (default 0.80)
#' @param mode          "single" or "differential" (default "single")
#' @param run_simulation TRUE = run full simulation after analytical step.
#'                       FALSE = analytical (CircaPower) only.
#' @param prior_result  A previous runSingleCohortGrid()/runDifferentialPower() result to reuse
#'                      instead of re-running simulation.
#' @param alpha2        2nd-harmonic deviation (scalar or vector)
#' @param alpha3        3rd-harmonic deviation (scalar or vector)
#' @param mc.cores      Parallel cores for simulation step
#' @param plot          Auto-plot final recommendation
#' @param output_file   PDF path (NULL = screen)
#' @param verbose       Print progress
#'
#' @return SCPRecommendResult with:
#'   $guidance      -- data.frame: method guidance table
#'   $analytical_df -- data.frame[method, N, power_analytical] (DCP only)
#'   $simulation    -- runSingleCohortGrid() or runDifferentialPower() result (NULL if not run)
#'   $recommendation -- data.frame[method, optimal_B, n_target, note]
#' @export
recommendDesign <- function(bio.opts,
                             design.opts,
                             analysis.opts,
                             methods       = c("DCP","JTK","RAIN","MH"),
                             target_power  = 0.80,
                             mode          = c("single","differential"),
                             run_simulation = TRUE,
                             prior_result  = NULL,
                             alpha2        = 0,
                             alpha3        = 0,
                             mc.cores      = 1L,
                             plot          = TRUE,
                             output_file   = NULL,
                             verbose       = TRUE) {

  stopifnot(inherits(bio.opts,      "CircadianBioOptions"))
  stopifnot(inherits(design.opts,   "CircadianDesignOptions"))
  stopifnot(inherits(analysis.opts, "CircadianAnalysisOptions"))
  mode <- match.arg(mode)

  N_vals <- design.opts$sample_sizes
  B_vals <- design.opts$B_values %||% length(unique(design.opts$cts %||% integer(0)))
  if (length(B_vals) == 0) B_vals <- 1L
  alpha  <- analysis.opts$alpha
  period <- bio.opts$period %||% 24

  # ------------------------------------------------------------------
  # Step 1: Print method guidance
  # ------------------------------------------------------------------
  guidance <- printMethodGuidance(methods = methods, verbose = verbose)

  # ------------------------------------------------------------------
  # Step 2: Analytical estimates (CircaPower, DCP only, B-invariant)
  # ------------------------------------------------------------------
  if (verbose) cat("=== Step 2: Analytical estimates (CircaPower / DCP) ===\n")

  # Median r from bio.opts
  if (!is.null(bio.opts$sigma_rhythmic) &&
      length(bio.opts$sigma_rhythmic) == length(bio.opts$amplitude)) {
    r_vec <- bio.opts$amplitude / bio.opts$sigma_rhythmic
  } else {
    r_vec <- bio.opts$amplitude / exp(median(bio.opts$lOD, na.rm = TRUE))
  }
  r_med <- median(r_vec[is.finite(r_vec) & r_vec > 0], na.rm = TRUE)

  analytical_df <- data.frame(
    method = "DCP",
    N      = N_vals,
    power  = vapply(N_vals, function(n) {
      tryCatch(CircaPower(n = n, power = NULL, r = r_med, alpha = alpha)$power,
               error = function(e) NA_real_)
    }, numeric(1)),
    note   = sprintf("B-invariant, r_med=%.2f", r_med),
    stringsAsFactors = FALSE
  )

  n80_analytical <- tryCatch(
    CircaPower(n = NULL, power = target_power, r = r_med, alpha = alpha)$n,
    error = function(e) NA_real_)

  if (verbose) {
    cat(sprintf("  Median pilot r = %.2f\n", r_med))
    cat(sprintf("  CircaPower n%d  = %s (DCP, any B)\n",
                round(target_power * 100),
                if (is.na(n80_analytical)) "not reached" else as.character(n80_analytical)))
    cat("  Note: JTK, RAIN, MH have no closed-form — simulation required.\n\n")
  }

  # ------------------------------------------------------------------
  # Step 3: Simulation (use prior_result or run fresh)
  # ------------------------------------------------------------------
  sim_result <- NULL

  if (!is.null(prior_result)) {
    if (!is.list(prior_result))
      stop("prior_result must be a list returned by runSingleCohortPower() or runDifferentialPower()")
    sim_result <- prior_result
    if (verbose) cat("=== Step 3: Using prior simulation result ===\n\n")

  } else if (run_simulation) {
    if (verbose) cat("=== Step 3: Running simulation ===\n")
    if (mode == "single") {
      sim_result <- runSingleCohortGrid(
        bio.opts, design.opts, analysis.opts,
        methods  = methods, alpha2 = alpha2, alpha3 = alpha3,
        mc.cores = mc.cores, verbose = verbose
      )
    } else {
      sim_result <- runDifferentialPower(
        bio.opts, design.opts, analysis.opts,
        methods  = methods,
        mc.cores = mc.cores, plot = FALSE, verbose = verbose
      )
    }
  } else {
    if (verbose) cat("=== Step 3: Skipped (run_simulation=FALSE) ===\n\n")
  }

  # ------------------------------------------------------------------
  # Step 4: Synthesise recommendation — optimal B per method
  # ------------------------------------------------------------------
  if (verbose) cat(sprintf("=== Step 4: Recommendation at target power=%.0f%% ===\n",
                            target_power * 100))

  # Start with analytical recommendation for DCP
  rec_rows <- list(data.frame(
    method    = "DCP",
    optimal_B = "any (B-invariant)",
    n_target  = n80_analytical %||% NA_real_,
    source    = "analytical",
    note      = sprintf("B does not affect power; invest in N (n%d≈%s)",
                        round(target_power * 100),
                        if (is.na(n80_analytical)) "not reached" else n80_analytical),
    stringsAsFactors = FALSE
  ))

  # Simulation-based recommendations for other methods
  if (!is.null(sim_result) && !is.null(sim_result$power_df)) {
    df <- sim_result$power_df
    sim_methods <- setdiff(unique(df$method), "DCP")

    for (mth in sim_methods) {
      sub <- df[df$method == mth & df$alpha2 == min(df$alpha2), ]
      if (nrow(sub) == 0) next

      # For each B: find smallest N reaching target_power
      best_B <- NA_integer_; best_n <- Inf
      for (b in sort(unique(sub$B))) {
        sb  <- sub[sub$B == b, ]
        sb  <- sb[order(sb$N), ]
        hit <- sb$N[sb$power >= target_power]
        if (length(hit) > 0 && hit[1] < best_n) {
          best_n <- hit[1]; best_B <- b
        }
      }
      rec_rows[[mth]] <- data.frame(
        method    = mth,
        optimal_B = if (is.na(best_B)) "not reached" else as.character(best_B),
        n_target  = if (is.infinite(best_n)) NA_real_ else best_n,
        source    = "simulation",
        note      = guidance$B_vs_m[guidance$method == mth][1] %||% "",
        stringsAsFactors = FALSE
      )
    }
  }

  recommendation <- do.call(rbind, rec_rows)

  if (verbose) {
    cat(sprintf("  %-14s  %-18s  %-10s  %s\n",
                "Method", "Optimal B", sprintf("n%d", round(target_power*100)), "Source"))
    cat(strrep("-", 70), "\n")
    for (i in seq_len(nrow(recommendation))) {
      cat(sprintf("  %-14s  %-18s  %-10s  %s\n",
                  recommendation$method[i],
                  recommendation$optimal_B[i],
                  recommendation$n_target[i] %||% "NA",
                  recommendation$source[i]))
    }
    cat("\n")
  }

  out <- structure(
    list(
      guidance       = guidance,
      analytical_df  = analytical_df,
      simulation     = sim_result,
      recommendation = recommendation,
      target_power   = target_power,
      mode           = mode
    ),
    class = "SCPRecommendResult"
  )

  if (plot) plot(out, output_file = output_file)
  out
}

#' @export
print.SCPRecommendResult <- function(x, ...) {
  cat("SCPRecommendResult\n")
  cat(sprintf("  mode: %s  |  target power: %.0f%%\n",
              x$mode, x$target_power * 100))
  cat("\nRecommendation:\n")
  print(x$recommendation[, c("method","optimal_B","n_target","note")])
  invisible(x)
}

#' @export
plot.SCPRecommendResult <- function(x, output_file = NULL, ...) {
  if (!is.null(x$simulation)) {
    plot(x$simulation, output_file = output_file, ...)
  } else if (requireNamespace("ggplot2", quietly = TRUE)) {
    df <- x$analytical_df
    p  <- ggplot2::ggplot(df, ggplot2::aes(x = N, y = power)) +
      ggplot2::geom_line(colour = "steelblue", linewidth = 1) +
      ggplot2::geom_point(colour = "steelblue") +
      ggplot2::geom_hline(yintercept = x$target_power,
                          linetype = "dashed", colour = "grey40") +
      ggplot2::labs(x = "N", y = "Power",
                    title = "CircaPower analytical (DCP, B-invariant)") +
      ggplot2::theme_bw()
    if (!is.null(output_file)) ggplot2::ggsave(output_file, p, width=6, height=4)
    else print(p)
  }
  invisible(x)
}


# =====================================================================
# Internal grid engine for B vs m recommendation
# =====================================================================

# Not exported. Called by recommendDesign() to sweep N x B x method x alpha2.
# Returns SCPSingleResult with $power_df (N, B, method, alpha2, power, power_se)
# and $n80_df — the compact format needed for the recommendation table.
#' Single-cohort B-versus-m power grid (active designs)
#'
#' Sweeps sample size, number of time points, and detection method to produce
#' the power grid used by \code{\link{recommendDesign}}. Active (balanced)
#' designs only.
#'
#' @param bio.opts Pilot summary (\code{CircadianBioOptions}).
#' @param design.opts Design options (\code{CircadianDesignOptions}).
#' @param analysis.opts Analysis options (\code{CircadianAnalysisOptions}).
#' @param methods Detection method(s); default "DCP".
#' @param alpha2,alpha3 Second- and third-harmonic generative coefficients (default 0).
#' @param mc.cores Number of cores.
#' @param plot Whether to draw the grid.
#' @param output_file Optional output path for the plot.
#' @param verbose Print progress.
#' @return An \code{SCPSingleResult} with the power grid and N80 table.
#' @export
runSingleCohortGrid <- function(bio.opts, design.opts, analysis.opts,
                                methods     = "DCP",
                                alpha2      = 0,
                                alpha3      = 0,
                                mc.cores    = 1L,
                                plot        = FALSE,
                                output_file = NULL,
                                verbose     = TRUE) {

  stopifnot(inherits(bio.opts,      "CircadianBioOptions"))
  stopifnot(inherits(design.opts,   "CircadianDesignOptions"))
  stopifnot(inherits(analysis.opts, "CircadianAnalysisOptions"))

  methods <- match.arg(methods, c("DCP", "cosinor", "Kharmonic", "JTK", "RAIN", "MH", "FMM"), several.ok = TRUE)

  # Passive designs cannot identify a B effect: cts is resampled from the
  # pilot TOD distribution regardless of B, so any cross-B variation is
  # simulation noise. Refuse to run a multi-B grid in passive mode to
  # avoid users misinterpreting the noise as signal.
  if (identical(design.opts$design, "passive") &&
      length(design.opts$B_values %||% integer(0)) > 1L)
    stop("runSingleCohortGrid: passive design cannot identify a B effect. ",
         "Use runSimsSingleCohort() with a single design for passive pilots, ",
         "or use design='active' with multiple B values.")

  N_vals      <- design.opts$sample_sizes
  B_vals      <- design.opts$B_values %||% length(unique(design.opts$cts))
  period      <- bio.opts$period %||% 24
  fdr         <- min(analysis.opts$fdr_thresholds)
  nsims       <- design.opts$nsims
  GLOBAL_SEED <- bio.opts$sim.seed %||% 2025L

  detect_fn <- list(
    DCP  = function(expr, cts) detect_cosinor(expr, cts, K = 1L, period = period),
    JTK  = function(expr, cts) detect_JTK(expr,  cts, period = period),
    RAIN = function(expr, cts) detect_RAIN(expr, cts, period = period),
    MH   = function(expr, cts) detect_MH(expr,   cts, period = period)
  )

  grid <- expand.grid(N = N_vals, B = B_vals, alpha2 = alpha2, alpha3 = alpha3,
                      method = methods, stringsAsFactors = FALSE)
  grid <- grid[grid$N %% grid$B == 0, ]

  if (verbose) {
    printMethodGuidance(methods = methods, verbose = TRUE)
    cat(sprintf("runSingleCohortGrid: %d cells x %d sims = %d runs\n",
                nrow(grid), nsims, nrow(grid) * nsims))
    cat(sprintf("  methods: %s\n", paste(methods, collapse = ", ")))
    cat(sprintf("  N: %s\n", paste(N_vals, collapse = ", ")))
    cat(sprintf("  B: %s\n", paste(B_vals, collapse = ", ")))
    cat(sprintf("  alpha2: %s\n", paste(alpha2, collapse = ", ")))
  }

  design_type <- design.opts$design %||% "active"
  pilot_cts   <- design.opts$cts   # used for passive KDE sampling

  run_cell <- function(i) {
    N   <- grid$N[i];  B  <- grid$B[i]
    a2  <- grid$alpha2[i]; a3 <- grid$alpha3[i]
    mth <- grid$method[i]
    fn  <- detect_fn[[mth]]

    # Respect design: active = equispaced; passive = sample from pilot TOD KDE
    cts_active <- rep(seq(0, period * (1 - 1/B), length.out = B), each = N / B)

    fmm_omega_i <- design.opts$omega %||% 1.0
    fmm_beta_i  <- design.opts$beta  %||% pi
    sims <- vapply(seq_len(nsims), function(s) {
      cts <- if (design_type == "passive" && !is.null(pilot_cts))
               sampleTimesFromDist(N, pilot_cts)
             else
               cts_active
      dat <- simCircadianSingleCohort(bio.opts, cts, alpha2 = a2, alpha3 = a3,
                                      omega = fmm_omega_i, beta = fmm_beta_i,
                                      seed = GLOBAL_SEED + i * 1000L + s)
      pv  <- fn(dat$expr, cts)
      adj <- p.adjust(pv, method = analysis.opts$p.adjust.method)
      sum(adj[dat$is_rhythmic] <= fdr, na.rm = TRUE) / sum(dat$is_rhythmic)
    }, numeric(1))

    c(power = mean(sims, na.rm = TRUE), power_se = sd(sims, na.rm = TRUE) / sqrt(nsims))
  }

  t0 <- proc.time()[["elapsed"]]
  results <- parallel::mclapply(seq_len(nrow(grid)), run_cell, mc.cores = mc.cores)
  if (verbose) cat(sprintf("  Done in %.1f min\n", (proc.time()[["elapsed"]] - t0) / 60))

  safe_get <- function(r, field) {
    if (is.null(r) || inherits(r, "try-error")) NA_real_ else r[[field]]
  }
  n_failed <- sum(vapply(results, function(r) is.null(r) || inherits(r, "try-error"), logical(1)))
  if (n_failed > 0)
    warning(sprintf("runSingleCohortGrid: %d/%d cells failed — set to NA",
                    n_failed, length(results)))
  grid$power    <- vapply(results, safe_get, 0, "power")
  grid$power_se <- vapply(results, safe_get, 0, "power_se")

  n80_df <- do.call(rbind, lapply(
    split(grid, interaction(grid$method, grid$B, grid$alpha2, grid$alpha3)),
    function(sub) {
      n80 <- tryCatch({
        sub  <- sub[order(sub$N), ]
        idx  <- which(sub$power >= 0.80)
        if (length(idx) == 0L) {
          NA_real_
        } else if (min(idx) == 1L) {
          sub$N[1L]
        } else {
          j2 <- min(idx); j1 <- j2 - 1L
          p1 <- sub$power[j1]; p2 <- sub$power[j2]
          n1 <- sub$N[j1];     n2 <- sub$N[j2]
          ceiling(n1 + (0.80 - p1) / (p2 - p1) * (n2 - n1))
        }
      }, error = function(e) NA_real_)
      data.frame(method = sub$method[1], B = sub$B[1],
                 alpha2 = sub$alpha2[1], alpha3 = sub$alpha3[1],
                 n80 = n80, stringsAsFactors = FALSE)
    }
  ))

  result <- structure(
    list(power_df = grid, n80_df = n80_df,
         bio.opts = bio.opts, design.opts = design.opts, analysis.opts = analysis.opts),
    class = "SCPSingleResult"
  )

  if (isTRUE(plot))
    plotBvsMPower(result, nsims = design.opts$nsims %||% 30L,
                  output_file = output_file)

  result
}


# =====================================================================
# Unified single-cohort power runner
# =====================================================================

#' Run single-cohort rhythmicity power analysis
#'
#' Sweeps N x B x alpha2 x alpha3 x method and returns a tidy power table.
#' B values come from design.opts$B_values; if NULL, B is inferred from cts.
#' Automatically calls plot() and npower() unless plot=FALSE.
#'
#' @param bio.opts      CircadianBioOptions from estCircadianParam()
#' @param design.opts   CircadianDesignOptions (sample_sizes, B_values, cts, design)
#' @param analysis.opts CircadianAnalysisOptions
#' @param methods       Character vector: any of "DCP","JTK","RAIN","MH"
#' @param alpha2        Scalar or vector of 2nd-harmonic amplitudes to sweep
#' @param alpha3        Scalar or vector of 3rd-harmonic amplitudes to sweep
#' @param mc.cores      Parallel cores
#' @param plot          Auto-plot on completion
#' @param output_file   PDF path for auto-plot (NULL = screen)
#' @param verbose       Print progress
#'
#' @return Rich list (from runSimsSingleCohort): $marginal_power [N x nsims],
#'   $strat_power [N x r_strata x nsims], $pvalues [N x genes x nsims],
#'   $sample_sizes, $r_values_list. Pass directly to plotSingleCohortPower().
#' @export
runSingleCohortPower <- function(bio.opts,
                                  design.opts,
                                  analysis.opts,
                                  methods     = "DCP",
                                  alpha2      = 0,
                                  alpha3      = 0,
                                  mc.cores    = 1L,
                                  plot        = TRUE,
                                  output_file = NULL,
                                  verbose     = TRUE) {

  stopifnot(inherits(bio.opts,      "CircadianBioOptions"))
  stopifnot(inherits(design.opts,   "CircadianDesignOptions"))
  stopifnot(inherits(analysis.opts, "CircadianAnalysisOptions"))

  methods <- match.arg(methods, c("DCP", "cosinor", "Kharmonic", "JTK", "RAIN", "MH", "FMM"), several.ok = TRUE)
  if (length(methods) > 1L) {
    warning("runSingleCohortPower: multiple methods specified; using only '", methods[1L],
            "'. For multi-method B vs m comparison use recommendDesign().")
    methods <- methods[1L]
  }

  res <- runSimsSingleCohort(bio.opts, design.opts, analysis.opts,
                              method    = methods,
                              harmonics = c(alpha2[1], alpha3[1]),
                              verbose   = verbose,
                              mc.cores  = mc.cores)

  if (isTRUE(plot))
    plotSingleCohortPower(res, out_pdf = output_file)

  res
}

#' @export
print.SCPSingleResult <- function(x, ...) {
  cat("SCPSingleResult\n")
  cat(sprintf("  methods: %s\n", paste(unique(x$power_df$method), collapse = ", ")))
  cat(sprintf("  N range: %d – %d\n", min(x$power_df$N), max(x$power_df$N)))
  cat(sprintf("  B values: %s\n", paste(sort(unique(x$power_df$B)), collapse = ", ")))
  cat(sprintf("  alpha2: %s\n", paste(sort(unique(x$power_df$alpha2)), collapse = ", ")))
  cat("\nn80 summary:\n")
  print(x$n80_df)
  invisible(x)
}

#' @export
plot.SCPSingleResult <- function(x, output_file = NULL, ...) {
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    df <- x$power_df
    p  <- ggplot2::ggplot(df, ggplot2::aes(x = B, y = power, colour = factor(N),
                                            group = factor(N))) +
      ggplot2::geom_line() + ggplot2::geom_point() +
      ggplot2::facet_grid(alpha2 ~ method, labeller = ggplot2::label_both) +
      ggplot2::geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey40") +
      ggplot2::labs(x = "B (time points)", y = "Power", colour = "N",
                    title = "Single-cohort rhythmicity power") +
      ggplot2::theme_bw()
    if (!is.null(output_file)) {
      ggplot2::ggsave(output_file, p, width = 3 * length(unique(df$method)), height = 4)
    } else {
      print(p)
    }
  }
  invisible(x)
}


# =====================================================================
# Unified two-group differential power runner
# =====================================================================

#' Run two-group differential circadian power analysis
#'
#' Sweeps N x alpha2 x alpha3 x method x test_type and returns a tidy power table.
#' Automatically calls plot() and npower() unless plot=FALSE.
#'
#' @param bio.opts      CircadianBioOptions from estCircadianParamTwoGroup()
#' @param design.opts   CircadianDesignOptions (sample_sizes, design, cts)
#' @param analysis.opts CircadianAnalysisOptions
#' @param methods       Any of "DCP","CircaCompare","LimoRhyde","DODR"
#' @param test_types    Any of "DR","DP","DM" (silently NA if method lacks support)
#' @param alpha2        Scalar or vector swept for both groups
#' @param alpha3        Scalar or vector swept for both groups
#' @param mc.cores      Parallel cores
#' @param plot          Auto-plot on completion
#' @param output_file   PDF path (NULL = screen)
#' @param verbose       Print progress
#'
#' @return List (class \code{SCPDiffResult}) from \code{runSimsDiff()}:
#'   \code{pval_DR}, \code{fdr_DR}, \code{pval_DP}, \code{fdr_DP},
#'   \code{pval_DM}, \code{fdr_DM} -- 3-D arrays \code{[ngenes x n_sizes x nsims]};
#'   plus \code{diff_type}, \code{effectsize}, \code{sample_sizes}, \code{nsims}.
#'   Pass to \code{plotDiffPower()} or \code{npower(..., endpoint="DR")}.
#' @export
runDifferentialPower <- function(bio.opts,
                                  design.opts,
                                  analysis.opts,
                                  methods     = "DCP",
                                  test_types  = c("DR", "DP", "DM"),
                                  alpha2      = 0,
                                  alpha3      = 0,
                                  mc.cores    = 1L,
                                  plot        = TRUE,
                                  output_file = NULL,
                                  verbose     = TRUE) {

  stopifnot(inherits(bio.opts,      "CircadianBioOptions"))
  stopifnot(inherits(design.opts,   "CircadianDesignOptions"))
  stopifnot(inherits(analysis.opts, "CircadianAnalysisOptions"))

  methods    <- match.arg(methods,    c("DCP","CircaCompare","LimoRhyde","DODR"), several.ok = TRUE)
  test_types <- match.arg(test_types, c("DR","DP","DM"), several.ok = TRUE)

  if (length(methods) > 1L) {
    warning("runDifferentialPower: multiple methods specified; using only '", methods[1L], "'.")
    methods <- methods[1L]
  }

  analysis.opts$DCmethod  <- methods
  design.opts$test_types  <- test_types
  design.opts$harmonics   <- c(alpha2[1], alpha3[1])

  res <- runSimsDiff(bio.opts, design.opts, analysis.opts,
                     mc.cores = mc.cores,
                     verbose  = verbose)

  # Attach SCPDiffResult class so S3 methods dispatch correctly.
  class(res) <- c("SCPDiffResult", "list")

  if (isTRUE(plot))
    plotDiffPower(list(res),
                  comp_labels = NULL,
                  endpoints   = intersect(test_types, c("DR", "DP", "DM")),
                  out_pdf     = output_file)

  res
}

#' @export
print.SCPDiffResult <- function(x, ...) {
  cat("SCPDiffResult\n")
  # x is a plain list from runSimsDiff(); use its actual fields.
  cat(sprintf("  DCmethod:   %s\n", x$sim_params$DCmethod %||% "DCP"))
  cat(sprintf("  test_types: %s\n",
              paste(intersect(c("DR","DP","DM"),
                              sub("^fdr_", "",
                                  names(x)[startsWith(names(x), "fdr_")])),
                    collapse = ", ")))
  cat(sprintf("  N range:    %d – %d\n",
              min(x$sample_sizes), max(x$sample_sizes)))
  cat(sprintf("  nsims:      %d\n", x$nsims))
  cat(sprintf("  ngenes:     %d\n", x$ngenes))
  invisible(x)
}

#' @export
plot.SCPDiffResult <- function(x, output_file = NULL, ...) {
  # x is a plain runSimsDiff() list — delegate to plotDiffPower() which
  # already knows how to consume this structure.
  endpoints <- intersect(c("DR", "DP", "DM"),
                         sub("^fdr_", "", grep("^fdr_", names(x), value = TRUE)))
  plotDiffPower(list(x),
                comp_labels = NULL,
                endpoints   = if (length(endpoints) > 0) endpoints else c("DR", "DP", "DM"),
                out_pdf     = output_file)
  invisible(x)
}
