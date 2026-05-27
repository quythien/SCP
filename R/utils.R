#' Utility Functions for Circadian Power Analysis
#'



# QUICK ANALYTICAL POWER 
# ==============================================================================

#' Quick analytical power for single-group rhythmicity
#'
#' @description Fast closed-form power calculation based on F-test.
#' Useful for quick estimates before running full simulations.
#' Based on CircaPower's analytical approach.
#'
#' @param n Sample size (NULL to solve for n)
#' @param power Target power (NULL to solve for power)
#' @param r Effect size (A/σ)
#' @param alpha Significance level
#' @param design_factor Sampling design factor (0.5 for even spacing)
#'
#' @return List with n, power, r, alpha
#' 
#' @examples
#' # How many samples for 80% power?
#' CircaPower(n = NULL, power = 0.8, r = 1.5)
#' 
#' # What power with n=24?
#' CircaPower(n = 24, power = NULL, r = 1.2)
#' 
#' @export
CircaPower = function(n = NULL, power = NULL, r = NULL, 
                      alpha = 0.05, design_factor = 0.5) {
  
  # Exactly one parameter must be NULL
  if (sum(sapply(list(n, r, power, alpha), is.null)) != 1) {
    stop("Exactly one of n, r, power, and alpha must be NULL")
  }
  
  # Based on F-test with ncp = r² × n × d
  k = 2  # 2 parameters (cos, sin)
  df1 = k
  
  if (is.null(power)) {
    # Calculate power
    df2 = n - k - 1
    ncp = r^2 * n * design_factor
    F_crit = qf(1 - alpha, df1, df2, ncp = 0)
    power = 1 - pf(F_crit, df1, df2, ncp = ncp)
    
  } else if (is.null(n)) {
    # Solve for n
    n = uniroot(function(n) {
      df2 = n - k - 1
      ncp = r^2 * n * design_factor
      F_crit = qf(1 - alpha, df1, df2, ncp = 0)
      (1 - pf(F_crit, df1, df2, ncp = ncp)) - power
    }, c(4, 1e9))$root
    n = ceiling(n)
    
  } else if (is.null(r)) {
    # Solve for r
    r = uniroot(function(r) {
      df2 = n - k - 1
      ncp = r^2 * n * design_factor
      F_crit = qf(1 - alpha, df1, df2, ncp = 0)
      (1 - pf(F_crit, df1, df2, ncp = ncp)) - power
    }, c(1e-10, 1e5))$root
    
  } else if (is.null(alpha)) {
    # Solve for alpha
    alpha = uniroot(function(alpha) {
      df2 = n - k - 1
      ncp = r^2 * n * design_factor
      F_crit = qf(1 - alpha, df1, df2, ncp = 0)
      (1 - pf(F_crit, df1, df2, ncp = ncp)) - power
    }, c(1e-10, 1 - 1e-10))$root
  }
  
  return(list(n = n, power = power, r = r, alpha = alpha, 
              design_factor = design_factor))
}

# ==============================================================================
# STRATIFIED POWER ANALYSIS
# ==============================================================================

#' Stratified Power Analysis
# ... (rest of utility_functions.R continues)


#' Stratified Power Analysis
#'
#' @description Calculate power stratified by effect size (r = A/sigma)
#' stratified power by expression level.
#'
#' @param pvals Vector of p-values
#' @param ground_truth Logical vector of true positives
#' @param effect_sizes Vector of effect sizes (r = A/σ)
#' @param breaks Break points for effect size strata
#' @param alpha Significance threshold
#'
#' @return Data frame with power by effect size stratum
stratified_power = function(pvals, ground_truth, effect_sizes,
                            breaks = c(0, 0.5, 1, 2, Inf),
                            alpha = 0.05) {

  qvals = p.adjust(pvals, method = "BH")
  discovered = qvals < alpha

  # Stratify by effect size
  strata = cut(effect_sizes, breaks = breaks,
               labels = paste0("r_", head(breaks, -1), "-", tail(breaks, -1)),
               include.lowest = TRUE)

  results = data.frame(
    stratum = levels(strata),
    n_total = as.numeric(table(strata)),
    n_true = as.numeric(tapply(ground_truth, strata, sum)),
    n_discovered = as.numeric(tapply(discovered, strata, sum)),
    stringsAsFactors = FALSE
  )

  results$TP = as.numeric(tapply(discovered & ground_truth, strata, sum))
  results$FP = results$n_discovered - results$TP
  results$FN = results$n_true - results$TP
  results$TN = results$n_total - results$n_discovered - results$FN

  results$power = results$TP / pmax(results$n_true, 1)
  results$FDR = results$FP / pmax(results$n_discovered, 1)
  results$precision = results$TP / pmax(results$n_discovered, 1)
  results$recall = results$power

  return(results)
}


#' False Discovery Cost (FDC)
#'
#' @description Calculate the cost of each true discovery, measured as the
#' number of false discoveries per true discovery.
#'
#' @param TP True positives
#' @param FP False positives
#'
#' @return FDC (FP / TP)
false_discovery_cost = function(TP, FP) {
  if (TP == 0) return(Inf)
  return(FP / TP)
}


#' Targeted Power (for biologically meaningful effects)
#'
#' @description Calculate power only for genes with effect sizes above
#' a meaningful threshold.
#'
#' @param pvals P-values
#' @param ground_truth Logical vector of true positives
#' @param effect_sizes Effect sizes (r = A/σ)
#' @param min_effect Minimum effect size to be considered "meaningful"
#' @param alpha Significance threshold
#'
#' @return Targeted power (proportion of meaningful effects discovered)
targeted_power = function(pvals, ground_truth, effect_sizes,
                          min_effect = 0.5, alpha = 0.05) {

  meaningful = ground_truth & effect_sizes >= min_effect

  qvals = p.adjust(pvals, method = "BH")
  discovered = qvals < alpha

  TP_meaningful = sum(discovered & meaningful)
  n_meaningful = sum(meaningful)

  return(TP_meaningful / max(n_meaningful, 1))
}


#' Power Assessment Summary
#'
#' @description Comprehensive power assessment with stratified and marginal metrics
#'
#' @param pvals P-values
#' @param ground_truth True positive indicator
#' @param effect_sizes Effect sizes
#' @param alpha Significance threshold
#'
#' @return List with comprehensive power metrics
power_assessment = function(pvals, ground_truth, effect_sizes, alpha = 0.05) {

  qvals = p.adjust(pvals, method = "BH")
  discovered = qvals < alpha

  # Overall metrics
  TP = sum(discovered & ground_truth)
  FP = sum(discovered & !ground_truth)
  FN = sum(!discovered & ground_truth)
  TN = sum(!discovered & !ground_truth)

  overall = list(
    power = TP / max(TP + FN, 1),
    FDR = FP / max(TP + FP, 1),
    precision = TP / max(TP + FP, 1),
    recall = TP / max(TP + FN, 1),
    specificity = TN / max(TN + FP, 1),
    FDC = false_discovery_cost(TP, FP)
  )

  # Stratified power
  stratified = stratified_power(pvals, ground_truth, effect_sizes,
                                breaks = c(0, 0.5, 1, 2, Inf), alpha = alpha)

  # Targeted power (meaningful effects only)
  targeted = list(
    r05 = targeted_power(pvals, ground_truth, effect_sizes, min_effect = 0.5, alpha = alpha),
    r10 = targeted_power(pvals, ground_truth, effect_sizes, min_effect = 1.0, alpha = alpha),
    r15 = targeted_power(pvals, ground_truth, effect_sizes, min_effect = 1.5, alpha = alpha)
  )

  return(list(
    overall = overall,
    stratified = stratified,
    targeted = targeted
  ))
}


#' Plot Power Curves
#'
#' @description Visualize power vs sample size for different methods
#'
#' @param power_results Results from power_single_condition
#' @param target_power Target power level (default 0.8)
#'
#' @return ggplot object
plot_power_curves = function(power_results, target_power = 0.8) {

  df = power_results$power_curves
  n_methods = length(grep("^power_", names(df)))
  methods = gsub("^power_", "", grep("^power_", names(df), value = TRUE))

  # Reshape for plotting
  plot_df = data.frame()
  for (method in methods) {
    temp = data.frame(
      n = df$n,
      power = df[[paste0("power_", method)]],
      fdr = df[[paste0("fdr_", method)]],
      method = method
    )
    plot_df = rbind(plot_df, temp)
  }

  p = ggplot(plot_df, aes(x = n, y = power, color = method)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_hline(yintercept = target_power, linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = 0.5, linetype = "dotted", color = "gray30") +
    labs(
      title = "Power vs Sample Size",
      subtitle = paste("Target power:", target_power),
      x = "Replicates per time point (n)",
      y = "Power",
      color = "Method"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "bottom"
    ) +
    ylim(0, 1)

  return(p)
}


#' Plot Two-Group Power Heatmap
#'
#' @description Heatmap of power for different (n_A, n_B) combinations
#'
#' @param power_matrix Power matrix from power_two_group
#' @param n_A_range n_A values
#' @param n_B_range n_B values
#'
#' @return ggplot object
plot_two_group_heatmap = function(power_matrix, n_A_range, n_B_range) {
  if (!requireNamespace("reshape2", quietly = TRUE)) {
    stop("Package 'reshape2' is required for this function. Please install it.")
  }
  
  
  df = reshape2::melt(power_matrix)
  names(df) = c("n_A", "n_B", "power")

  p = ggplot(df, aes(x = n_A, y = n_B, fill = power)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f", power)), color = "white", size = 3) +
    scale_fill_gradient2(
      low = "red", mid = "yellow", high = "green",
      midpoint = 0.5, limits = c(0, 1),
      name = "Power"
    ) +
    labs(
      title = "Two-Group Comparison Power",
      x = "Replicates per time point (Condition A)",
      y = "Replicates per time point (Condition B)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )

  return(p)
}


#' Plot Stratified Power
#'
#' @description Bar plot of power by effect size stratum
#'
#' @param stratified_results Results from stratified_power
#'
#' @return ggplot object
plot_stratified_power = function(stratified_results) {

  p = ggplot(stratified_results, aes(x = stratum, y = power)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_text(aes(label = sprintf("%.2f", power)), vjust = -0.5) +
    geom_errorbar(aes(ymin = 0, ymax = power), width = 0.2) +
    labs(
      title = "Stratified Power by Effect Size",
      x = "Effect Size Stratum (r = A/σ)",
      y = "Power"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    ylim(0, 1)

  return(p)
}


#' Compare Multiple Methods
#'
#' @description Compare power, FDR, and FDC across methods
#'
#' @param power_results_list List of results from multiple methods
#'
#' @return Comparison summary table
compare_methods = function(power_results_list) {

  methods = names(power_results_list)

  comparison = data.frame(
    Method = methods,
    Power_Overall = sapply(power_results_list, function(x) x$overall$power),
    FDR = sapply(power_results_list, function(x) x$overall$FDR),
    Precision = sapply(power_results_list, function(x) x$overall$precision),
    Recall = sapply(power_results_list, function(x) x$overall$recall),
    FDC = sapply(power_results_list, function(x) x$overall$FDC),
    Power_r05 = sapply(power_results_list, function(x) x$targeted$r05),
    Power_r10 = sapply(power_results_list, function(x) x$targeted$r10),
    stringsAsFactors = FALSE
  )

  return(comparison)
}


#' Generate Power Report
#'
#' @description Generate comprehensive power analysis report
#'
#' @param power_results Results from power analysis
#' @param output_file Output file path (NULL = return as text)
#'
#' @return Report text or NULL if written to file
generate_power_report = function(power_results, output_file = NULL) {

  report = c()
  report = c(report, "=== CIRCADIAN POWER ANALYSIS REPORT ===\n")
  report = c(report, sprintf("Generated: %s\n\n", Sys.time()))

  # Parameters
  report = c(report, "=== SIMULATION PARAMETERS ===\n")
  params = power_results$parameters
  report = c(report, sprintf("Time points: %s\n", paste(range(params$times), collapse = " to ")))
  report = c(report, sprintf("Sample sizes tested: %s\n", paste(params$n_range, collapse = ", ")))
  report = c(report, sprintf("Simulations per sample size: %d\n", params$n_sim))
  report = c(report, sprintf("Genes per simulation: %d\n", params$G))
  report = c(report, sprintf("Proportion rhythmic: %.1f%%\n", 100 * params$prop_rhythmic))
  report = c(report, sprintf("Alpha level: %.3f\n\n", params$alpha))

  # Optimal sample size
  report = c(report, "=== OPTIMAL SAMPLE SIZE ===\n")
  for (method in names(power_results$optimal_n)) {
    n_opt = power_results$optimal_n[[method]]
    report = c(report, sprintf("%s: n = %s\n", method,
                               ifelse(is.na(n_opt), "Not achieved", n_opt)))
  }
  report = c(report, "\n")

  # Power at different sample sizes
  report = c(report, "=== POWER CURVES ===\n")
  df = power_results$power_curves
  report = c(report, capture.output(print(df, row.names = FALSE)))
  report = c(report, "\n")

  # Computation time
  report = c(report, sprintf("Computation time: %.1f minutes\n", as.numeric(power_results$elapsed_time)))

  report_text = paste(report, collapse = "")

  if (!is.null(output_file)) {
    writeLines(report_text, output_file)
    return(invisible(NULL))
  } else {
    return(report_text)
  }
}


# ==============================================================================
# SIGNAL-TO-NOISE RATIO (r) ESTIMATION FROM REAL DATA
# ==============================================================================

#' Estimate Signal-to-Noise Ratio (r) from DCP Results
#'
#' @description Estimate the median signal-to-noise ratio r = A/σ from the
#' top circadian genes in a DCP (Differential Circadian Pipeline) analysis.
#' Uses the relationship r = 1/sqrt(1/R² - 1) where R² is rhythm strength.
#'
#' @param dcp_results List from DCP analysis containing DR and DP components
#' @param top_n Number of top circadian genes to use (default: 100)
#' @param rank_by Column to rank genes by (default: "p.overall" from DP results)
#' @param group Which group to use for R² ("min" for minimum across groups, "1" or "2" for specific group)
#'
#' @return List with median r, mean r, and individual r values
#'
#' @examples
#' # dcp_results is the output from runDiffCircadian()
#' r_estimate <- estimate_r_from_data(dcp_results, top_n = 100)
#' # Use median r for power simulation
#' target_r <- r_estimate$median_r
#'
#' @export
estimate_r_from_data = function(dcp_results, top_n = 100,
                                rank_by = "q.R2",
                                use_group = 1) {

  # Extract DR (Differential Rhythm) results
  if (is.null(dcp_results$DR)) {
    stop("dcp_results must contain DR component with R² values")
  }

  dr <- dcp_results$DR

  # Get top rhythmic genes by significance (default: DR q-values)
  # This ranks by how rhythmic genes are, NOT by DP significance
  # For power analysis, we want the top circadian genes from control/pilot group
  if (!rank_by %in% names(dr)) {
    stop(paste("rank_by must be a column in DR. Options:", paste(names(dr), collapse=", ")))
  }

  top_genes <- head(dr[order(dr[[rank_by]]), "gname"], top_n)

  # Merge with DR data to get R² values
  dr_top <- merge(dr, data.frame(gname = top_genes), by = "gname")

  # Calculate r from R² for each group: r = 1/sqrt(1/R² - 1)
  # Handle edge cases where R² = 0 or R² = 1
  dr_top$r1 <- with(dr_top, {
    valid <- R2.1 > 0 & R2.1 < 1
    r <- ifelse(valid, 1/sqrt(1/R2.1 - 1), NA)
    r[!valid] <- ifelse(R2.1 <= 0, 0, Inf)
    r
  })

  dr_top$r2 <- with(dr_top, {
    valid <- R2.2 > 0 & R2.2 < 1
    r <- ifelse(valid, 1/sqrt(1/R2.2 - 1), NA)
    r[!valid] <- ifelse(R2.2 <= 0, 0, Inf)
    r
  })

  # Select which group to use for r estimation
  # use_group = 1: control/pilot group (Group 1)
  # use_group = 2: treatment group (Group 2)
  # use_group = "min": minimum across both groups (conservative)
  if (use_group == "min") {
    dr_top$r <- pmin(dr_top$r1, dr_top$r2, na.rm = TRUE)
    group_label <- "min(G1, G2)"
  } else if (use_group == 1) {
    dr_top$r <- dr_top$r1
    group_label <- "Group 1 (control/pilot)"
  } else if (use_group == 2) {
    dr_top$r <- dr_top$r2
    group_label <- "Group 2 (treatment)"
  } else {
    stop("use_group must be 1, 2, or 'min'")
  }

  # Remove infinite/NA values for summary statistics
  r_finite <- dr_top$r[is.finite(dr_top$r)]

  # Return results
  result <- list(
    median_r = median(r_finite, na.rm = TRUE),
    mean_r = mean(r_finite, na.rm = TRUE),
    sd_r = sd(r_finite, na.rm = TRUE),
    r_values = dr_top$r,
    top_n = top_n,
    n_valid = length(r_finite),
    group_used = group_label,
    dr_top = dr_top
  )

  class(result) <- "estimate_r_from_data"
  return(result)
}


#' Print method for estimate_r_from_data
#'
#' @param x Result from estimate_r_from_data
#' @param ... Additional arguments (ignored)
#'
#' @export
print.estimate_r_from_data = function(x, ...) {
  cat("=== Signal-to-Noise Ratio (r) Estimation ===\n")
  cat(sprintf("Based on top %d rhythmic genes (ranked by DR q-value)\n", x$top_n))
  cat(sprintf("Group used: %s\n", x$group_used))
  cat(sprintf("Valid r values: %d\n\n", x$n_valid))
  cat(sprintf("Median r: %.3f\n", x$median_r))
  cat(sprintf("Mean r:   %.3f\n", x$mean_r))
  cat(sprintf("SD r:     %.3f\n", x$sd_r))
  cat("\nSuggested r for simulation: ", round(x$median_r, 1), "\n")
  invisible(x)
}


#' Get Recommended r Stratum for Power Simulation
#'
#' @description Given estimated r, return the appropriate r stratum label
#' for use in power simulation stratification.
#'
#' @param r_estimate Result from estimate_r_from_data or numeric r value
#' @param r_strata Standard r strata breaks (default: c(0, 0.25, 0.5, 0.75, 1, ...))
#'
#' @return List with stratum index, label, and suggested r value
#'
#' @export
get_r_stratum = function(r_estimate, r_strata = c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.5, 3, 3.5, 4, 4.5, 5, Inf)) {

  if (inherits(r_estimate, "estimate_r_from_data")) {
    r_value <- r_estimate$median_r
  } else {
    r_value <- as.numeric(r_estimate)
  }

  # Find which stratum
  stratum_idx <- cut(r_value, breaks = r_strata, include.lowest = TRUE, labels = FALSE)

  # Get label
  strata_labels <- c("(0,0.25]", "(0.25,0.5]", "(0.5,0.75]", "(0.75,1]",
                    "(1,1.25]", "(1.25,1.5]", "(1.5,1.75]", "(1.75,2]",
                    "(2,2.5]", "(2.5,3]", "(3,3.5]", "(3.5,4]",
                    "(4,4.5]", "(4.5,5]", ">5")
  stratum_label <- strata_labels[stratum_idx]

  return(list(
    r_value = r_value,
    stratum_idx = stratum_idx,
    stratum_label = stratum_label,
    suggested_r_for_sim = round(r_value * 2) / 2  # Round to nearest 0.5
  ))
}


# ==============================================================================
# DATA PREPARATION
# ==============================================================================

#' Prepare expression data for circadian power analysis
#'
#' @description
#' Converts raw expression data into the standard (matrix, times) format
#' expected by \code{estCircadianParam()}, \code{estCircadianParamTwoGroup()},
#' and \code{runBootstrapDesignGrid()}.
#'
#' The pipeline always expects:
#'   \itemize{
#'     \item \code{data}  — numeric matrix, \strong{genes × samples}, log2-scale
#'     \item \code{times} — numeric vector, hours in [0, 24), same length as \code{ncol(data)}
#'   }
#'
#' This function handles the three most common raw formats:
#' \describe{
#'   \item{"counts"}{Raw integer counts (e.g. featureCounts, STAR).
#'     Converts to log2(CPM + 1): \code{log2(sweep(x, 2, colSums(x), "/") * 1e6 + 1)}.}
#'   \item{"cpm"}{Raw CPM values (e.g. edgeR \code{cpm()}, baboon data).
#'     Converts to log2(CPM + 1).}
#'   \item{"log2"}{Already log2-scale (e.g. GSE54651, Seney ACC).
#'     Wraps with \code{data.matrix()} only; no numeric transformation.}
#' }
#'
#' @param expr  A matrix, data.frame, or file path (CSV/TSV) of expression values.
#'              Rows = genes, columns = samples. Row names used as gene IDs.
#' @param times Numeric vector of time-of-day values (hours, 0–24) for each sample,
#'              OR a character column name in \code{pheno} to extract times from.
#' @param input_type One of \code{"counts"}, \code{"cpm"}, or \code{"log2"}.
#'              Default \code{"counts"}.
#' @param pheno Optional data.frame of sample metadata. Must have the same column
#'              order as \code{expr} if \code{times} is a column name string.
#' @param sample_col Optional: name of the column in \code{pheno} whose values match
#'              \code{colnames(expr)} for alignment. If NULL, assumes columns are
#'              already in the same order as rows of \code{pheno}.
#' @param time_col  Name of the TOD column in \code{pheno} (used when \code{times}
#'              is a column name string).
#' @param period  Circadian period in hours. Default 24. Used only for validation.
#' @param verbose Print a short summary of what was done. Default TRUE.
#'
#' @return A list with:
#'   \item{data}{Numeric matrix (genes × samples), log2-scale}
#'   \item{times}{Numeric vector of hours (0–24), length == ncol(data)}
#'   \item{n_genes}{Number of genes}
#'   \item{n_samples}{Number of samples}
#'   \item{input_type}{The \code{input_type} used}
#'
#' @examples
#' # From raw counts matrix
#' prep <- prepCircadianData(counts_matrix, times = tod_vector, input_type = "counts")
#' bio  <- estCircadianParam(prep$data, prep$times)
#'
#' # From data.frame with metadata
#' prep <- prepCircadianData(
#'   expr       = expr_df,
#'   times      = "TOD",
#'   input_type = "log2",
#'   pheno      = pheno_df,
#'   sample_col = "SampleID"
#' )
#'
#' # Baboon raw CPM (load gives a data.frame)
#' prep <- prepCircadianData(baboon_df, times = tod_vector, input_type = "cpm")
#' @export
prepCircadianData <- function(expr,
                              times,
                              input_type = c("counts", "cpm", "log2"),
                              pheno      = NULL,
                              sample_col = NULL,
                              time_col   = NULL,
                              period     = 24,
                              verbose    = TRUE) {

  input_type <- match.arg(input_type)

  # --- 1. Load from file if a path was given ---
  if (is.character(expr) && length(expr) == 1 && file.exists(expr)) {
    if (verbose) cat(sprintf("Reading expression from: %s\n", expr))
    sep <- if (grepl("\\.tsv$|\\.txt$", expr, ignore.case=TRUE)) "\t" else ","
    expr <- as.matrix(read.csv(expr, row.names = 1, sep = sep, check.names = FALSE))
  }

  # --- 2. Coerce to numeric matrix ---
  mat <- data.matrix(expr)
  if (!is.numeric(mat))
    stop("Expression data must be numeric after coercion.")

  # --- 3. Align pheno and extract times ---
  if (is.character(times) && length(times) == 1) {
    # times is a column name — extract from pheno
    if (is.null(pheno))
      stop("'times' is a column name but 'pheno' was not provided.")
    time_col_use <- times

    if (!is.null(sample_col)) {
      # Align pheno rows to expression columns
      idx <- match(colnames(mat), as.character(pheno[[sample_col]]))
      if (any(is.na(idx)))
        warning(sprintf("%d expression columns could not be matched in pheno.", sum(is.na(idx))))
      pheno_aligned <- pheno[idx, , drop = FALSE]
    } else {
      pheno_aligned <- pheno
    }

    times_vec <- as.numeric(pheno_aligned[[time_col_use]])
  } else {
    times_vec <- as.numeric(times)
  }

  # --- 4. Validate alignment ---
  if (length(times_vec) != ncol(mat))
    stop(sprintf(
      "Length of times (%d) must equal ncol(expr) (%d).",
      length(times_vec), ncol(mat)))

  # Remove samples with missing times
  ok <- !is.na(times_vec)
  if (any(!ok)) {
    if (verbose)
      message(sprintf("Removing %d samples with NA times.", sum(!ok)))
    mat       <- mat[, ok, drop = FALSE]
    times_vec <- times_vec[ok]
  }

  # Warn if times are outside [0, 24)
  if (any(times_vec < 0 | times_vec >= 24, na.rm = TRUE))
    warning("Some times are outside [0, 24). Verify TOD extraction is correct.")

  # --- 5. Normalize ---
  if (input_type == "counts") {
    lib_size <- colSums(mat)
    if (any(lib_size == 0))
      stop("Some samples have zero library size. Check input.")
    mat <- log2(sweep(mat, 2, lib_size, "/") * 1e6 + 1)
    norm_desc <- "log2(CPM + 1) from raw counts"

  } else if (input_type == "cpm") {
    mat <- log2(mat + 1)
    norm_desc <- "log2(CPM + 1) from raw CPM"

  } else {
    # log2: already normalized — just ensure numeric matrix
    norm_desc <- "as-is (already log2-scale)"
  }

  # --- 6. Summary ---
  if (verbose) {
    cat(sprintf(
      "prepCircadianData: %d genes x %d samples | input=%s | norm=%s\n",
      nrow(mat), ncol(mat), input_type, norm_desc))
    cat(sprintf("  TOD range: [%.1f, %.1f] h  |  unique times: %d\n",
                min(times_vec), max(times_vec), length(unique(times_vec))))
  }

  list(
    data       = mat,
    times      = times_vec,
    n_genes    = nrow(mat),
    n_samples  = ncol(mat),
    input_type = input_type
  )
}
