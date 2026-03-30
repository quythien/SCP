# DR 4-panel figure (A, B, D, F) with n up to 160
# Usage:
#   Rscript examples/replot_dr_ABDF_panelC160.R

setwd("/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
source("code/plot_with_se.R")

load("output/run_panelC160/dr_power_raw_pvalues.rds")

obj <- dr_power_raw
sample_sizes <- obj$sample_sizes
nsims <- obj$nsims
r_strata <- obj$r_strata
strata_labels <- obj$strata_labels
pvalues <- obj$pvalues

n_sizes <- length(sample_sizes)
n_r_strata <- length(r_strata) - 1
fdr_thresholds <- c(0.01, 0.05, 0.10, 0.20)
threshold_labels <- paste0("FDR ", fdr_thresholds * 100, "%")
threshold_colors <- c("darkgreen", "steelblue", "orange", "red")

# Prepare arrays
power_arr <- array(NA, dim = c(n_sizes, n_r_strata, length(fdr_thresholds), nsims))
TD_arr <- array(NA, dim = c(n_sizes, n_r_strata, length(fdr_thresholds), nsims))
FD_arr <- array(NA, dim = c(n_sizes, n_r_strata, length(fdr_thresholds), nsims))

nested_gt <- is.list(obj$is_target_list[[1]]) && length(obj$is_target_list) == n_sizes

for (j in 1:n_sizes) {
  for (s in 1:nsims) {
    pvals <- pvalues[j, , s]

    if (nested_gt) {
      r_vec <- obj$r_values_list[[j]][[s]]
      is_target <- obj$is_target_list[[j]][[s]]
      is_null <- obj$is_null_list[[j]][[s]]
    } else {
      r_vec <- obj$r_values_list[[s]]
      is_target <- obj$is_target_list[[s]]
      is_null <- obj$is_null_list[[s]]
    }

    xgr <- cut(r_vec, breaks = r_strata, include.lowest = TRUE, labels = FALSE)

    qvals <- rep(1, length(pvals))
    tested <- pvals < 1
    if (sum(tested) > 0) qvals[tested] <- p.adjust(pvals[tested], method = "BH")

    for (t in 1:length(fdr_thresholds)) {
      discoveries <- qvals <= fdr_thresholds[t]
      for (k in 1:n_r_strata) {
        in_stratum <- xgr == k
        td <- sum(discoveries & is_target & in_stratum, na.rm = TRUE)
        fd <- sum(discoveries & is_null & in_stratum, na.rm = TRUE)
        nt <- sum(is_target & in_stratum, na.rm = TRUE)
        TD_arr[j, k, t, s] <- td
        FD_arr[j, k, t, s] <- fd
        power_arr[j, k, t, s] <- if (nt > 0) td / nt else NA
      }
    }
  }
}

# Marginal power by n and threshold
marginal_mean <- matrix(NA, nrow = n_sizes, ncol = length(fdr_thresholds))
marginal_se <- matrix(NA, nrow = n_sizes, ncol = length(fdr_thresholds))
for (j in 1:n_sizes) {
  for (t in 1:length(fdr_thresholds)) {
    td <- apply(TD_arr[j, , t, ], 2, sum, na.rm = TRUE)
    nt <- apply(obj$strat_n_targets[j, , ], 2, sum, na.rm = TRUE)
    power <- ifelse(nt > 0, td / nt, NA_real_)
    marginal_mean[j, t] <- mean(power, na.rm = TRUE)
    marginal_se[j, t] <- sd(power, na.rm = TRUE) / sqrt(sum(!is.na(power)))
  }
}

# Indices
idx_n60 <- which.min(abs(sample_sizes - 60))
idx_fdr5 <- 2

# Panel A: power by r (FDR 5%), lines per n (<=100)
idx_other <- which(sample_sizes <= 100)
size_colors <- rainbow(length(idx_other), s = 0.6, v = 0.8)
mean_power_A <- apply(power_arr[idx_other, , idx_fdr5, ], c(1, 2), mean, na.rm = TRUE)
se_power_A <- apply(power_arr[idx_other, , idx_fdr5, ], c(1, 2),
                    function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

# Panel B: power by r at different FDR thresholds (n=60)
mean_power_B <- apply(power_arr[idx_n60, , , ], c(1, 2), mean, na.rm = TRUE)
se_power_B <- apply(power_arr[idx_n60, , , ], c(1, 2),
                    function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

# Panel F: min r stratum reaching 80% power at n=60
min_r_stratum <- rep(NA, n_r_strata)
for (k in 1:n_r_strata) {
  pwr <- mean_power_B[k, idx_fdr5]
  if (!is.na(pwr) && pwr >= 0.80) {
    min_r_stratum[k] <- k
  }
}

out_file <- "output/run_panelC160/figures/dr_abdf_panelC160.pdf"
pdf(out_file, width = 12, height = 9)
par(mfrow = c(2, 2), mai = c(0.9, 1.0, 0.5, 0.2), mgp = c(3, 0.5, 0), oma = c(0, 0, 2, 0))

# Panel A
matplot(1:n_r_strata, t(mean_power_A),
        type = "l", lwd = 2, col = size_colors, lty = 1,
        xlim = c(0.5, n_r_strata + 0.5), ylim = c(0, 1),
        xlab = expression(r == A/sigma), ylab = "Power",
        main = "DR Power by r (FDR 5%)", xaxt = "n")
axis(1, at = 1:n_r_strata, labels = strata_labels, las = 2, cex.axis = 0.6)
for (j in seq_along(idx_other)) {
  points(1:n_r_strata, mean_power_A[j, ], pch = 19, col = size_colors[j], cex = 0.6)
  add_se_bars(1:n_r_strata, mean_power_A[j, ], se_power_A[j, ], col = size_colors[j])
}
grid()
legend("bottomright", paste0("n=", sample_sizes[idx_other]), col = size_colors, lty = 1, lwd = 2, cex = 0.6)
mtext("A", side = 3, at = -0.02, font = 2, line = 0.5)

# Panel B
matplot(1:n_r_strata, mean_power_B,
        type = "l", lwd = 2, col = threshold_colors, lty = 1,
        xlim = c(0.5, n_r_strata + 0.5), ylim = c(0, 1),
        xlab = expression(r == A/sigma), ylab = "Power",
        main = sprintf("DR Power by r (n=%d)", sample_sizes[idx_n60]), xaxt = "n")
axis(1, at = 1:n_r_strata, labels = strata_labels, las = 2, cex.axis = 0.6)
for (t in 1:length(fdr_thresholds)) {
  points(1:n_r_strata, mean_power_B[, t], pch = 19, col = threshold_colors[t], cex = 0.6)
  add_se_bars(1:n_r_strata, mean_power_B[, t], se_power_B[, t], col = threshold_colors[t])
}
grid()
legend("bottomright", threshold_labels, col = threshold_colors, lty = 1, pch = 19, lwd = 2, cex = 0.7)
mtext("B", side = 3, at = -0.02, font = 2, line = 0.5)

# Panel D (marginal power vs n, FDR 5%, n up to 160)
plot(sample_sizes, 100 * marginal_mean[, idx_fdr5],
     type = "b", pch = 19, lwd = 2, col = "steelblue",
     xlim = c(0, max(sample_sizes) * 1.1), ylim = c(0, 100),
     xlab = "Sample Size (per group)", ylab = "Marginal Power (%)",
     xaxt = "n",
     main = "DR Marginal Power vs Sample Size (FDR 5%)")
axis(1, at = sample_sizes, labels = sample_sizes)
add_se_bars(sample_sizes, 100 * marginal_mean[, idx_fdr5],
            100 * marginal_se[, idx_fdr5], col = "steelblue")
abline(h = 80, lty = 2, col = "gray")
grid()
text(sample_sizes, 100 * marginal_mean[, idx_fdr5] + 4,
     sprintf("%.1f%%", 100 * marginal_mean[, idx_fdr5]), cex = 0.6, col = "steelblue")
mtext("D", side = 3, at = -0.02, font = 2, line = 0.5)

# Panel F: min r stratum reaching 80% power at n=60
bar_vals <- ifelse(is.na(min_r_stratum), NA, min_r_stratum)
bar_vals[is.na(bar_vals)] <- n_r_strata
bar_colors <- ifelse(is.na(min_r_stratum), "gray80", "steelblue")
bar_colors[bar_vals == n_r_strata & is.na(min_r_stratum)] <- "gray80"

bp <- barplot(bar_vals, names.arg = strata_labels,
              col = bar_colors, border = "gray40",
              ylim = c(0, n_r_strata + 1),
              xlab = expression(r == A/sigma), ylab = "Min r Stratum",
              main = sprintf("Min r Stratum for 80%% Power (n=%d)", sample_sizes[idx_n60]),
              las = 2, cex.names = 0.6)
for (i in 1:n_r_strata) {
  if (is.na(min_r_stratum[i])) {
    text(bp[i], bar_vals[i] + 0.3, "N/R", cex = 0.55, font = 3, col = "gray40")
  }
}
grid()
mtext("F", side = 3, at = -0.02, font = 2, line = 0.5)

mtext("Differential Rhythmicity Power Analysis (Panel C up to n=160)", outer = TRUE, font = 2, cex = 1.1)
dev.off()

cat("Saved:", out_file, "\n")
