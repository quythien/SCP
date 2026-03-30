# Phase shift figure with panels A, B, D, F only (2x2 layout)
# Usage:
#   Rscript examples/replot_phase_shift_ABDF.R

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

out_file <- "output/run_20260219_2156_key/figures/phase_shift_abdf.pdf"
pdf(out_file, width = 12, height = 9)
par(mfrow = c(2, 2), mai = c(0.9, 1.0, 0.5, 0.2), mgp = c(3, 0.5, 0), oma = c(0, 0, 2, 0))

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

# Panel D
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
mtext("C", side = 3, at = -0.02, font = 2, line = 0.5)

# Panel F (min detectable shift by r at n=60; histogram style)
min_shift <- rep(NA, n_r_strata)
for (k in 1:n_r_strata) {
  for (p in 1:n_phase) {
    pwr <- mean(strat_power[p, idx_n60, k, ], na.rm = TRUE)
    if (!is.na(pwr) && pwr >= 0.80) {
      if (p > 1) {
        pwr_prev <- mean(strat_power[p - 1, idx_n60, k, ], na.rm = TRUE)
        if (!is.na(pwr_prev) && pwr > pwr_prev) {
          frac <- (0.80 - pwr_prev) / (pwr - pwr_prev)
          min_shift[k] <- phase_shifts[p - 1] + frac * (phase_shifts[p] - phase_shifts[p - 1])
        } else {
          min_shift[k] <- phase_shifts[p]
        }
      } else {
        min_shift[k] <- phase_shifts[p]
      }
      break
    }
  }
}
has_data <- sapply(1:n_r_strata, function(k) {
  any(!is.na(strat_power[, idx_n60, k, ]), na.rm = TRUE)
})
bar_vals <- min_shift
bar_vals[is.na(bar_vals)] <- max(phase_shifts)
bar_colors_f <- ifelse(is.na(min_shift), "gray80",
                       ifelse(min_shift <= 2, "steelblue",
                              ifelse(min_shift <= 6, "orange", "darkred")))
bar_colors_f[!has_data] <- "white"
bp <- barplot(bar_vals, names.arg = strata_labels,
              col = bar_colors_f, border = "gray40",
              ylim = c(0, max(phase_shifts) + 1),
              xlab = expression(r == A/sigma), ylab = "Min Phase Shift (hours)",
              main = "Min Detectable Shift for 80% Power (n=60)",
              las = 2, cex.names = 0.6)
for (i in 1:n_r_strata) {
  if (!has_data[i]) next
  if (is.na(min_shift[i])) {
    text(bp[i], bar_vals[i] + 0.3, "N/R", cex = 0.55, font = 3, col = "gray40")
  } else {
    text(bp[i], bar_vals[i] + 0.3, sprintf("%.1fh", min_shift[i]), cex = 0.55, font = 2)
  }
}
abline(h = 6, lty = 2, col = "darkgreen", lwd = 1.5)
text(bp[n_r_strata], 6.3, "6h shift", cex = 0.6, col = "darkgreen", pos = 2)
abline(h = 2, lty = 3, col = "steelblue", lwd = 1.5)
text(bp[n_r_strata], 2.3, "2h shift", cex = 0.6, col = "steelblue", pos = 2)
grid()
legend("topright",
       c("<= 2h", "2-6h", "> 6h", "Not reached"),
       fill = c("steelblue", "orange", "darkred", "gray80"),
       cex = 0.6, border = "gray40")
mtext("D", side = 3, at = -0.02, font = 2, line = 0.5)

mtext("Differential Phase Power Analysis: Effect of Phase Shift Magnitude", outer = TRUE, font = 2, cex = 1.1)
dev.off()

cat("Saved:", out_file, "\n")
