# Alternate phase-shift figure with revised panels
# - Replaces Panel C with minimum detectable phase shift (80% power) vs n
# - Replaces Panel E with marginal power vs n for fixed phase shifts (4h, 6h)
#
# Usage:
#   Rscript examples/replot_phase_shift_alt_panels.R

setwd("/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
source("code/plot_with_se.R")

load("output/run_20260219_2156_key/dp_power_phase_shift_results.rds")

obj <- dp_phase_results
phase_shifts <- obj$phase_shifts
sample_sizes <- obj$sample_sizes
nsims <- obj$nsims
r_strata <- obj$r_strata
strata_labels <- obj$strata_labels
strat_power <- obj$strat_power
strat_TD <- obj$strat_TD
strat_n_targets <- obj$strat_n_targets

# Drop 0h (null)
keep_phase <- phase_shifts != 0
phase_shifts <- phase_shifts[keep_phase]
strat_power <- strat_power[keep_phase, , , , drop = FALSE]
strat_TD <- strat_TD[keep_phase, , , , drop = FALSE]
strat_n_targets <- strat_n_targets[keep_phase, , , , drop = FALSE]

n_phase <- length(phase_shifts)
n_size <- length(sample_sizes)
n_r_strata <- length(r_strata) - 1

idx_n60 <- which.min(abs(sample_sizes - 60))
show_phase_idx <- which(phase_shifts %in% c(2, 4, 6, 8, 10, 12))
show_phase_labels <- phase_shifts[show_phase_idx]

show_r_idx <- c(4, 6, 8, 10, 12, 14)
show_r_labels <- strata_labels[show_r_idx]
show_r_colors <- rainbow(6, s = 0.6, v = 0.8)

# Panel A stats
mean_power_n60 <- array(NA, dim = c(n_phase, n_r_strata))
se_power_n60 <- array(NA, dim = c(n_phase, n_r_strata))
for (p in 1:n_phase) {
  for (k in 1:n_r_strata) {
    vals <- strat_power[p, idx_n60, k, ]
    mean_power_n60[p, k] <- mean(vals, na.rm = TRUE)
    se_power_n60[p, k] <- sd(vals, na.rm = TRUE) / sqrt(sum(!is.na(vals)))
  }
}

# Panel D stats: marginal power by n, phase
marginal_sim <- array(NA, dim = c(n_size, n_phase, nsims))
for (j in 1:n_size) {
  for (p in 1:n_phase) {
    for (s in 1:nsims) {
      td <- sum(strat_TD[p, j, , s], na.rm = TRUE)
      tgt <- sum(strat_n_targets[p, j, , s], na.rm = TRUE)
      marginal_sim[j, p, s] <- if (tgt > 0) td / tgt else NA
    }
  }
}
marginal_mean <- apply(marginal_sim, c(1, 2), mean, na.rm = TRUE)
marginal_se <- apply(marginal_sim, c(1, 2), function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

# Panel C (new): minimum detectable phase shift at 80% power vs n
target_power <- 0.70
min_shift_80 <- rep(NA_real_, n_size)
for (j in 1:n_size) {
  # marginal power by phase at sample size j
  pow_j <- marginal_mean[j, ]
  idx <- which(pow_j >= target_power)
  if (length(idx) > 0) {
    min_shift_80[j] <- phase_shifts[min(idx)]
  }
}

# Panel E (new): marginal power vs n at fixed shifts
fixed_shifts <- c(4, 6)
fixed_idx <- match(fixed_shifts, phase_shifts)

out_file <- "output/run_20260219_2156_key/figures/phase_shift_alt.pdf"
pdf(out_file, width = 16, height = 10)
par(mfrow = c(2, 3), mai = c(0.9, 1.0, 0.5, 0.2), mgp = c(3, 0.5, 0), oma = c(0, 0, 2, 0))

phase_colors <- rainbow(length(show_phase_idx), s = 0.6, v = 0.8)

# Panel A
matplot(1:n_r_strata, 100 * t(mean_power_n60[show_phase_idx, ]),
        type = "b", pch = 19, lwd = 2,
        col = phase_colors, lty = 1,
        xlim = c(0.5, n_r_strata + 0.5), ylim = c(0, 100),
        xlab = expression(r == A/sigma),
        ylab = "Power (%)",
        main = "Power vs r by Phase Shift (n=60)",
        xaxt = "n")
axis(1, at = 1:n_r_strata, labels = strata_labels, las = 2, cex.axis = 0.6)
for (ii in seq_along(show_phase_idx)) {
  p <- show_phase_idx[ii]
  add_se_bars(1:n_r_strata, 100 * mean_power_n60[p, ], 100 * se_power_n60[p, ], col = phase_colors[ii])
}
abline(h = 80, lty = 2, col = "gray")
grid()
legend("bottomright", paste0(show_phase_labels, "h"),
       col = phase_colors, lty = 1, pch = 19, lwd = 2, cex = 0.7)
mtext("A", side = 3, at = -0.02, font = 2, line = 0.5)

# Panel B (heatmap)
z_matrix <- 100 * t(mean_power_n60)
image(z_matrix,
      xlab = "Phase Shift (hours)", ylab = expression(r == A/sigma),
      main = "Power Heatmap (n=60)",
      col = colorRampPalette(c("white", "yellow", "orange", "red", "darkred"))(100),
      xaxt = "n", yaxt = "n")
axis(1, at = seq(0, 1, length.out = n_phase), labels = phase_shifts)
axis(2, at = seq(0, 1, length.out = n_r_strata), labels = strata_labels, las = 2, cex.axis = 0.6)
box()
mtext("B", side = 3, at = -0.02, font = 2, line = 0.5)

# Panel C (new)
plot(sample_sizes, min_shift_80, type = "n",
     xlab = "Sample Size (per group)",
     ylab = sprintf("Min Phase Shift for %.0f%% Power (h)", 100 * target_power),
     main = "Minimum Detectable Phase Shift",
     ylim = c(0, max(phase_shifts)))
grid()
if (all(is.na(min_shift_80))) {
  text(mean(sample_sizes), max(phase_shifts) * 0.6,
       sprintf("%.0f%% power not reached", 100 * target_power), cex = 0.9)
} else {
  lines(sample_sizes, min_shift_80, type = "b", pch = 19, lwd = 2)
}
mtext("C", side = 3, at = -0.02, font = 2, line = 0.5)

# Panel D (marginal power vs n by phase shift)
matplot(sample_sizes, 100 * marginal_mean[, show_phase_idx],
        type = "b", pch = 19, lwd = 2,
        col = phase_colors, lty = 1,
        xlim = c(0, 110), ylim = c(0, 100),
        xlab = "Sample Size (per group)",
        ylab = "Marginal Power (%)",
        main = "Marginal Power vs n by Phase Shift")
abline(h = 80, lty = 2, col = "gray")
grid()
legend("bottomright", paste0(show_phase_labels, "h"),
       col = phase_colors, lty = 1, pch = 19, lwd = 2, cex = 0.7)
mtext("D", side = 3, at = -0.02, font = 2, line = 0.5)

# Panel E (new): fixed shifts
cols_fixed <- c("steelblue", "darkorange")
matplot(sample_sizes, 100 * marginal_mean[, fixed_idx],
        type = "b", pch = 19, lwd = 2,
        col = cols_fixed, lty = 1,
        xlim = c(0, 110), ylim = c(0, 100),
        xlab = "Sample Size (per group)",
        ylab = "Marginal Power (%)",
        main = "Marginal Power at Fixed Shifts")
abline(h = 80, lty = 2, col = "gray")
grid()
legend("bottomright", paste0(fixed_shifts, "h"),
       col = cols_fixed, lty = 1, pch = 19, lwd = 2, cex = 0.7)
mtext("E", side = 3, at = -0.02, font = 2, line = 0.5)

# Panel F (existing: minimum detectable phase shift by r)
power_by_phase <- matrix(NA, nrow = n_phase, ncol = length(show_r_idx))
se_by_phase <- matrix(NA, nrow = n_phase, ncol = length(show_r_idx))
for (p in 1:n_phase) {
  for (ri in seq_along(show_r_idx)) {
    k <- show_r_idx[ri]
    vals <- strat_power[p, idx_n60, k, ]
    power_by_phase[p, ri] <- mean(vals, na.rm = TRUE)
    se_by_phase[p, ri] <- sd(vals, na.rm = TRUE) / sqrt(sum(!is.na(vals)))
  }
}
matplot(phase_shifts, 100 * power_by_phase,
        type = "b", pch = 19, lwd = 2,
        col = show_r_colors, lty = 1,
        xlim = c(0, 13), ylim = c(0, 100),
        xlab = "Phase Shift (hours)",
        ylab = "Power (%)",
        main = "Power vs Phase Shift by r (n=60)")
for (ri in seq_along(show_r_idx)) {
  add_se_bars(phase_shifts, 100 * power_by_phase[, ri], 100 * se_by_phase[, ri], col = show_r_colors[ri])
}
abline(h = 80, lty = 2, col = "gray")
abline(v = 6, lty = 3, col = "darkgreen", lwd = 1.5)
grid()
legend("bottomright", show_r_labels,
       col = show_r_colors, lty = 1, pch = 19, lwd = 2, cex = 0.6)
mtext("F", side = 3, at = -0.02, font = 2, line = 0.5)

mtext("Phase Shift Sensitivity (Alternate Panels)", outer = TRUE, font = 2, cex = 1.2)
dev.off()

cat("Saved:", out_file, "\n")
