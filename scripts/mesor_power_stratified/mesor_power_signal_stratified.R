#' MESOR Power Analysis - Save Raw p-values for Multi-Threshold Analysis
#' Differential Basal/MESOR (M difference) - Stratified by r = A/σ
#' Saves raw p-values (more flexible than q-values) for exact power calculation
#' at any p-value or FDR threshold

setwd("/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")

source("code/options.R")
source("code/simulation.R")
source("code/detection.R")
source("code/runner.R")
source("code/power.R")

cat("====================================================================\n")
cat("MESOR POWER ANALYSIS - SAVE RAW p-VALUES FOR MULTI-THRESHOLD\n")
cat("Stratifying by amplitude (signal strength)\n")
cat("====================================================================\n\n")

# Load real data for time distribution
COMBINED <- readRDS("/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Collaborative/Paper/Congruence/PNAS_aging/data/combined_data.rds")

expr_sample_names <- colnames(COMBINED$expr)
pheno_order <- match(expr_sample_names, COMBINED$pheno$sample_name)
valid_samples <- !is.na(pheno_order)
expr_sample_names <- expr_sample_names[valid_samples]
COMBINED$expr <- COMBINED$expr[, valid_samples]
pheno_order <- pheno_order[valid_samples]
pheno_data <- COMBINED$pheno[pheno_order, ]
pheno_data$tod <- if("TOD.x" %in% colnames(pheno_data)) pheno_data$TOD.x else pheno_data$TOD.y
pheno_data$age_group_final <- if("AgeGroup" %in% colnames(pheno_data)) pheno_data$AgeGroup else pheno_data$age_group
complete_samples <- !is.na(pheno_data$age_group_final) & !is.na(pheno_data$tod) &
  pheno_data$age_group_final %in% c("younger", "older")
pheno_clean <- pheno_data[complete_samples, ]
younger_idx <- pheno_clean$age_group_final == "younger"
times_young <- pheno_clean$tod[younger_idx]

cat(sprintf("Reference data: n=%d\n\n", length(times_young)))

# =====================================================================
# SIMULATION PARAMETERS
# =====================================================================

sample_sizes <- c(10, 20, 40, 60, 80, 100)
nsims <- 50
target_effect <- 0.5  # Threshold for "interesting" MESOR genes (0.5 expression units)

cat(sprintf("Sample sizes: %s\n", paste(sample_sizes, collapse=", ")))
cat(sprintf("Simulations: %d per scenario\n", nsims))
cat(sprintf("Target effect threshold: %.2f\n\n", target_effect))

# r = A/σ strata
r_strata <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.5, 3, 3.5, 4, 4.5, 5, Inf)
strata_labels <- c("(0,0.25]", "(0.25,0.5]", "(0.5,0.75]", "(0.75,1]",
                  "(1,1.25]", "(1.25,1.5]", "(1.5,1.75]", "(1.75,2]",
                  "(2,2.5]", "(2.5,3]", "(3,3.5]", "(3.5,4]",
                  "(4,4.5]", "(4.5,5]", ">5")

# =====================================================================
# RUN SIMULATIONS - MESOR (Differential Basal)
# =====================================================================

cat("====================================================================\n")
cat("RUNNING SIMULATIONS...\n")
cat("====================================================================\n\n")

# Storage for raw p-values: [sample_size, genes, nsims]
pvalues_MESOR <- array(NA, dim = c(length(sample_sizes), 5000, nsims))

# Storage for ground truth per simulation
r_values_list <- vector("list", nsims)
is_target_list <- vector("list", nsims)
is_null_list <- vector("list", nsims)

# Storage for stratified power at FDR 0.05 (for quick summary)
strat_power <- array(NA, dim = c(length(sample_sizes), length(r_strata)-1, nsims))
strat_TD <- array(NA, dim = c(length(sample_sizes), length(r_strata)-1, nsims))
strat_FD <- array(NA, dim = c(length(sample_sizes), length(r_strata)-1, nsims))
strat_n_targets <- array(NA, dim = c(length(sample_sizes), length(r_strata)-1, nsims))

for (j in seq_along(sample_sizes)) {
  n <- sample_sizes[j]
  cat(sprintf(">>> n = %d\n", n))

  sim_out <- runSimsDiff(
    sample_sizes = c(n),
    nsims = nsims,
    ngenes = 5000,
    prop_rhythmic = 0.30,
    prop_DR = 0.00,
    prop_DP = 0.00,
    prop_DA = 0.15,
    phase_diff = c(0, 0),  # No phase shift
    amp_diff = c(1, 1),     # No amplitude change
    cts = times_young,
    design = "passive",
    test_types = c("DA"),  # DA test actually includes basal (MESOR) in our implementation
    verbose = FALSE
  )

  pval_MESOR <- sim_out$pval_DA[, 1, ]  # genes x nsims
  fdr_MESOR <- sim_out$fdr_DA[, 1, ]

  pvalues_MESOR[j, , ] <- pval_MESOR

  for (i in 1:nsims) {
    diff_type <- sim_out$diff_type[[i]]
    effectsize_amp <- sim_out$effectsize[[i]]$amp
    effectsize_DR1 <- sim_out$effectsize[[i]]$DR1
    effectsize_DR2 <- sim_out$effectsize[[i]]$DR2
    r_values <- pmin(effectsize_DR1, effectsize_DR2)

    # DA genes (diff_type = 5) include both amplitude and basal changes
    # For MESOR analysis, we focus on genes with significant basal/mesor difference
    is_DA <- diff_type == 5

    if (j == 1) {
      r_values_list[[i]] <- r_values
      # Target: DA genes with large effect size (both amp and basal contribute)
      is_target_list[[i]] <- is_DA & (effectsize_amp >= target_effect)
      is_null_list[[i]] <- !is_DA
    }

    # Calculate stratified power at FDR 0.05 for summary
    Zg <- ifelse(is_DA, 1, 0)
    Zg2 <- ifelse(is_DA & effectsize_amp >= target_effect, 1, 0)

    r_for_strat <- r_values
    r_for_strat[!is_DA] <- 0
    n_strata <- length(r_strata) - 1
    xgr <- cut(r_for_strat, breaks = r_strata, include.lowest = TRUE, labels = FALSE)
    xgr[!is_DA] <- NA

    discoveries <- fdr_MESOR[, i] <= 0.05

    for (k in 1:n_strata) {
      in_stratum <- xgr == k
      TD <- sum(discoveries & Zg2 == 1 & in_stratum, na.rm = TRUE)
      FD <- sum(discoveries & Zg == 0 & in_stratum, na.rm = TRUE)
      n_targets <- sum(Zg2 == 1 & in_stratum, na.rm = TRUE)

      if (n_targets > 0) {
        strat_power[j, k, i] <- TD / n_targets
      } else {
        strat_power[j, k, i] <- NA
      }
      strat_TD[j, k, i] <- TD
      strat_FD[j, k, i] <- FD
      strat_n_targets[j, k, i] <- n_targets
    }
  }

  cat(sprintf("  n=%d: ", n))
  for (k in 1:(length(r_strata)-1)) {
    mean_p <- mean(strat_power[j, k, ], na.rm = TRUE)
    cat(sprintf("%s=%.0f%% ", strata_labels[k], 100 * mean_p))
  }
  cat("\n\n")
}

# =====================================================================
# SAVE RESULTS
# =====================================================================

cat("====================================================================\n")
cat("SAVING RESULTS\n")
cat("====================================================================\n\n")

mesor_power_raw <- list(
  sample_sizes = sample_sizes,
  nsims = nsims,
  target_effect = target_effect,
  r_strata = r_strata,
  strata_labels = strata_labels,
  pvalues = pvalues_MESOR,
  r_values_list = r_values_list,
  is_target_list = is_target_list,
  is_null_list = is_null_list,
  strat_power = strat_power,
  strat_TD = strat_TD,
  strat_FD = strat_FD,
  strat_n_targets = strat_n_targets
)

# Create output directory if needed
if (!dir.exists("output/mesor_power_stratified")) {
  dir.create("output/mesor_power_stratified", recursive = TRUE)
}

save(mesor_power_raw, file = "output/mesor_power_stratified/mesor_power_raw_pvalues.rds")
cat("Results saved: output/mesor_power_stratified/mesor_power_raw_pvalues.rds\n")
cat("  (Contains raw p-values for exact multi-threshold analysis)\n")
cat("  Can calculate power at any p-value or FDR threshold\n\n")

cat("\nAll simulations complete!\n")
