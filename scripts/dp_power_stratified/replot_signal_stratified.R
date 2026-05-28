#' Replot Signal Stratified Power - Without re-running simulations
#' Stratified by r = A/σ (CircaPower-style) - 15 bins
#' With CI and power vs n plot

setwd("/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")

# Load saved results
load("output/dp_power_stratified/dp_power_signal_stratified_results.rds")
# Results are now in dp_power_results list
sample_sizes <- dp_power_results$sample_sizes
nsims <- dp_power_results$nsims
phase_shift <- dp_power_results$phase_shift
target_effect <- dp_power_results$target_effect
r_strata <- dp_power_results$r_strata
strata_labels <- dp_power_results$strata_labels
strat_power <- dp_power_results$strat_power
strat_TD <- dp_power_results$strat_TD
strat_FD <- dp_power_results$strat_FD
strat_n_targets <- dp_power_results$strat_n_targets
strat_n_null <- dp_power_results$strat_n_null
r_list <- dp_power_results$r_list

pdf("output/figures/dp_power_signal_stratified.pdf", width = 16, height = 12)

par(mfrow = c(2, 2), mai = c(2.2, 1.3, 0.6, 0.3), mgp = c(4, 0.6, 0))

# Colors for 15 strata
colors <- rainbow(15, s = 0.6, v = 0.8)

# =====================================================================
# Plot A: Stratified Power by r (PROPER-style)
# X-axis: r bins, Lines: Sample sizes
# =====================================================================

mean_power <- apply(strat_power, c(1, 2), mean, na.rm = TRUE)
matplot(1:15, t(100 * mean_power),  # Transpose: strata on x-axis, sizes as lines
        type = "l", lwd = 2,
        col = rainbow(6, s = 0.6, v = 0.8), lty = 1,  # 6 colors for 6 sample sizes
        xlim = c(0.5, 15.5), ylim = c(0, 100),
        xlab = "r = A/sigma (Signal-to-Noise Ratio)",
        ylab = "Targeted Power (%)",
        main = "Stratified Power by Signal-to-Noise Ratio",
        xaxt = "n")
axis(1, at = 1:15, labels = strata_labels, las = 2, cex.axis = 0.6)
abline(h = 80, lty = 2, col = "gray")
grid()
legend("bottomright", paste0("n=", sample_sizes),
       col = rainbow(6, s = 0.6, v = 0.8), lty = 1, lwd = 2, cex = 0.7)
mtext("A", side = 3, at = -0.02, font = 2)

# =====================================================================
# Plot B: Power by r Stratum at n=60 with Error Bars (95% CI)
# =====================================================================

idx_n60 <- which.min(abs(sample_sizes - 60))
power_n60_matrix <- strat_power[idx_n60, , ]  # strata x sims
mean_power_n60 <- apply(power_n60_matrix, 1, mean, na.rm = TRUE)
se_power_n60 <- apply(power_n60_matrix, 1, sd, na.rm = TRUE) / sqrt(apply(power_n60_matrix, 1, function(x) sum(!is.na(x))))
ci_power_n60 <- 1.96 * se_power_n60

# Plot with error bars
par(mai = c(2.2, 1.3, 0.6, 0.3), mgp = c(4, 0.6, 0))  # Extra margins for axis labels
plot(1:15, 100 * mean_power_n60,
     type = "b", pch = 19, col = "darkblue", lwd = 2,
     xlim = c(0.5, 15.5), ylim = c(0, 100),
     xlab = "r = A/sigma (Signal-to-Noise Ratio)",
     ylab = "Targeted Power (%)",
     main = "Power at n=60 by r (with 95% CI)",
     xaxt = "n")
axis(1, at = 1:15, labels = strata_labels, las = 2, cex.axis = 0.6)
abline(h = 80, lty = 2, col = "gray")
grid()

# Add error bars
arrows(1:15, 100 * (mean_power_n60 - ci_power_n60),
       1:15, 100 * (mean_power_n60 + ci_power_n60),
       angle = 90, code = 3, length = 0.02, col = "darkblue", lwd = 1.5)

mtext("B", side = 3, at = -0.02, font = 2)

# =====================================================================
# Plot C: Distribution of Differential Phase Genes Across r Strata
# Shows what proportion of phase-shifting genes fall in each signal-to-noise bin
# =====================================================================

par(mai = c(2.2, 1.3, 0.6, 0.3), mgp = c(4, 0.6, 0))  # Extra margins for axis labels
gene_counts <- apply(strat_n_targets[idx_n60, , ], 1, mean, na.rm = TRUE)
total_genes <- sum(gene_counts, na.rm = TRUE)
prop_genes <- 100 * gene_counts / total_genes

barplot(prop_genes, names.arg = strata_labels,
        col = colors, las = 2, cex.names = 0.6,
        ylim = c(0, 25),
        xlab = "r = A/sigma (Signal-to-Noise Ratio)",
        ylab = "Percent of Genes with Phase Shift (%)",
        main = "Distribution of Phase-Shifting Genes by r")
abline(h = 0, lwd = 1)
grid(nx = NA, ny = TRUE)
mtext("C", side = 3, at = -0.02, font = 2)

# =====================================================================
# Plot D: Marginal Power vs Sample Size
# =====================================================================

par(mai = c(2.2, 1.3, 0.6, 0.3), mgp = c(4, 0.6, 0))  # Extra margins for axis labels
marginal_power_by_n <- numeric(length(sample_sizes))
for (j in seq_along(sample_sizes)) {
  marginal_TD <- sum(apply(strat_TD[j, , ], 1, mean, na.rm = TRUE))
  marginal_targets <- sum(apply(strat_n_targets[j, , ], 1, mean, na.rm = TRUE))
  marginal_power_by_n[j] <- if (marginal_targets > 0) marginal_TD / marginal_targets else NA
}

plot(sample_sizes, 100 * marginal_power_by_n,
     type = "b", pch = 19, col = "darkred", lwd = 2,
     xlim = c(0, 110), ylim = c(0, 100),
     xlab = "Sample Size (per group)",
     ylab = "Marginal Power (%)",
     main = "Marginal Power vs Sample Size")
abline(h = 80, lty = 2, col = "gray")
abline(v = 60, lty = 3, col = "darkgreen", lwd = 2)
abline(v = 62, lty = 2, col = "steelblue", lwd = 1.5)
abline(v = 74, lty = 2, col = "darkorange", lwd = 1.5)
grid()
legend("bottomright",
       c("n=60 (simulated)", "Young (n=62)", "Aging (n=74)"),
       col = c("darkgreen", "steelblue", "darkorange"), lty = c(3, 2, 2), lwd = 2, cex = 0.7)
mtext("D", side = 3, at = -0.02, font = 2)

dev.off()

# Print summary table
cat("\n")
cat("====================================================================\n")
cat("STRATIFIED POWER SUMMARY (CircaPower-style) - n=5000 genes, 50 sims\n")
cat("====================================================================\n\n")

cat("Sample Size | ")
for (k in 1:length(strata_labels)) {
  cat(sprintf(" %-4s", gsub("[\\(\\)]", "", strata_labels[k])))
}
cat(" | Marginal\n")
cat(paste0(rep("-", 15 + length(strata_labels)*5), collapse = ""), "\n")

for (j in seq_along(sample_sizes)) {
  cat(sprintf("n = %-7d | ", sample_sizes[j]))

  marginal_TD <- 0
  marginal_targets <- 0

  for (k in 1:length(strata_labels)) {
    mean_p <- mean(strat_power[j, k, ], na.rm = TRUE)
    cat(sprintf(" %-4s", sprintf("%.0f", 100 * mean_p)))

    marginal_TD <- marginal_TD + mean(strat_TD[j, k, ], na.rm = TRUE)
    marginal_targets <- marginal_targets + mean(strat_n_targets[j, k, ], na.rm = TRUE)
  }

  marginal_power <- if (marginal_targets > 0) marginal_TD / marginal_targets else NA
  cat(sprintf(" | %.1f%%\n", 100 * marginal_power))
}

cat("\n")
cat("Figures saved: output/figures/dp_power_signal_stratified.pdf\n")
