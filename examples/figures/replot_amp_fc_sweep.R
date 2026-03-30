# Amplitude fold-change sensitivity figure (4-panel, analogous to phase_shift_abdf.pdf)
# Usage:
#   Rscript examples/replot_amp_fc_sweep.R

setwd("/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
source("code/plot_with_se.R")

amp_fc_results <- readRDS("output/run_panelC160/da_amp_fc_sweep_results.rds")

amp_fcs       <- c(1.25, 1.5, 2, 3, 4)
fc_labels     <- paste0(amp_fcs, "x")
reference_n   <- 60
fdr_threshold <- 0.05

# Extract metadata from first element
res1         <- amp_fc_results[[1]]
sample_sizes <- res1$sample_sizes
strata_labels <- res1$strata_labels
r_strata     <- res1$r_strata
nsims        <- res1$nsims
n_size       <- length(sample_sizes)
n_r_strata   <- length(strata_labels)
n_fc         <- length(amp_fcs)

idx_n60 <- which.min(abs(sample_sizes - reference_n))
fc_colors <- rainbow(n_fc, s = 0.6, v = 0.8)

# -----------------------------------------------------------------------
# Compute stratified power at n=60 and marginal power for all n
# -----------------------------------------------------------------------

strat_power_n60 <- matrix(NA, nrow = n_fc, ncol = n_r_strata)
strat_se_n60    <- matrix(NA, nrow = n_fc, ncol = n_r_strata)
marginal_mean   <- matrix(NA, nrow = n_size, ncol = n_fc)

for (f in seq_along(amp_fcs)) {
  res <- amp_fc_results[[f]]

  # Stratified power at n=60
  for (k in 1:n_r_strata) {
    vals <- res$strat_power[idx_n60, k, ]
    strat_power_n60[f, k] <- mean(vals, na.rm = TRUE)
    strat_se_n60[f, k]    <- sd(vals, na.rm = TRUE) / sqrt(sum(!is.na(vals)))
  }

  # Marginal power across all n (FDR 5%)
  for (j in 1:n_size) {
    td  <- apply(res$strat_TD[j, , ], 1, sum, na.rm = TRUE)
    tgt <- apply(res$strat_n_targets[j, , ], 1, sum, na.rm = TRUE)
    mp  <- ifelse(tgt > 0, td / tgt, NA_real_)
    marginal_mean[j, f] <- mean(mp, na.rm = TRUE)
  }
}

# Min detectable fold-change by r stratum at n=60
min_fc <- rep(NA, n_r_strata)
for (k in 1:n_r_strata) {
  for (f in 1:n_fc) {
    pwr <- strat_power_n60[f, k]
    if (!is.na(pwr) && pwr >= 0.80) {
      if (f > 1) {
        pwr_prev <- strat_power_n60[f - 1, k]
        if (!is.na(pwr_prev) && pwr > pwr_prev) {
          frac <- (0.80 - pwr_prev) / (pwr - pwr_prev)
          min_fc[k] <- amp_fcs[f - 1] + frac * (amp_fcs[f] - amp_fcs[f - 1])
        } else {
          min_fc[k] <- amp_fcs[f]
        }
      } else {
        min_fc[k] <- amp_fcs[f]
      }
      break
    }
  }
}

# -----------------------------------------------------------------------
# Plot
# -----------------------------------------------------------------------

out_file <- "output/run_panelC160/figures/amp_fc_sweep_abcd.pdf"
pdf(out_file, width = 12, height = 9)
par(mfrow = c(2, 2), mai = c(0.9, 1.0, 0.5, 0.2), mgp = c(3, 0.5, 0), oma = c(0, 0, 2, 0))

# Panel A: Power vs r by fold-change at n=60
matplot(1:n_r_strata, 100 * t(strat_power_n60),
        type = "b", pch = 19, lwd = 2,
        col = fc_colors, lty = 1,
        xlim = c(0.5, n_r_strata + 0.5), ylim = c(0, 100),
        xlab = expression(r == A/sigma),
        ylab = "Power (%)",
        main = sprintf("DA Power vs r by Fold-Change (n=%d)", reference_n),
        xaxt = "n")
axis(1, at = 1:n_r_strata, labels = strata_labels, las = 2, cex.axis = 0.6)
for (f in 1:n_fc) {
  add_se_bars(1:n_r_strata, 100 * strat_power_n60[f, ], 100 * strat_se_n60[f, ], col = fc_colors[f])
}
abline(h = 80, lty = 2, col = "gray")
grid()
legend("bottomright", fc_labels, col = fc_colors, lty = 1, pch = 19, lwd = 2, cex = 0.7)
mtext("A", side = 3, at = -0.02, font = 2, line = 0.5)

# Panel B: Heatmap (fold-change x r) at n=60
z_matrix <- 100 * t(strat_power_n60)
image(z_matrix,
      xlab = expression(r == A/sigma), ylab = "Amplitude Fold-Change",
      main = sprintf("DA Power Heatmap (n=%d)", reference_n),
      col = colorRampPalette(c("white", "yellow", "orange", "red", "darkred"))(100),
      xaxt = "n", yaxt = "n")
axis(1, at = seq(0, 1, length.out = n_r_strata), labels = strata_labels, las = 2, cex.axis = 0.55)
axis(2, at = seq(0, 1, length.out = n_fc), labels = fc_labels, las = 2, cex.axis = 0.7)
box()
mtext("B", side = 3, at = -0.02, font = 2, line = 0.5)

# Panel C: Marginal power vs n by fold-change
matplot(sample_sizes, 100 * marginal_mean,
        type = "b", pch = 19, lwd = 2,
        col = fc_colors, lty = 1,
        xlim = c(0, max(sample_sizes) * 1.1), ylim = c(0, 100),
        xlab = "Sample Size (per group)",
        ylab = "Marginal Power (%)",
        main = "DA Marginal Power vs n by Fold-Change")
abline(h = 80, lty = 2, col = "gray")
grid()
legend("bottomright", fc_labels, col = fc_colors, lty = 1, pch = 19, lwd = 2, cex = 0.7)
mtext("C", side = 3, at = -0.02, font = 2, line = 0.5)

# Panel D: Min detectable fold-change by r at n=60
bar_vals <- min_fc
bar_vals[is.na(bar_vals)] <- max(amp_fcs)
bar_colors_d <- ifelse(is.na(min_fc), "gray80",
                       ifelse(min_fc <= 1.5, "steelblue",
                              ifelse(min_fc <= 2, "orange", "darkred")))
bp <- barplot(bar_vals, names.arg = strata_labels,
              col = bar_colors_d, border = "gray40",
              ylim = c(0, max(amp_fcs) + 0.5),
              xlab = expression(r == A/sigma),
              ylab = "Min Amplitude Fold-Change",
              main = sprintf("Min Detectable Fold-Change for 80%% Power (n=%d)", reference_n),
              las = 2, cex.names = 0.6)
for (i in 1:n_r_strata) {
  if (is.na(min_fc[i])) {
    text(bp[i], bar_vals[i] + 0.1, "N/R", cex = 0.55, font = 3, col = "gray40")
  } else {
    text(bp[i], bar_vals[i] + 0.1, sprintf("%.2fx", min_fc[i]), cex = 0.55, font = 2)
  }
}
abline(h = 2, lty = 2, col = "darkgreen", lwd = 1.5)
text(bp[n_r_strata], 2.1, "2x", cex = 0.6, col = "darkgreen", pos = 2)
abline(h = 1.5, lty = 3, col = "steelblue", lwd = 1.5)
text(bp[n_r_strata], 1.6, "1.5x", cex = 0.6, col = "steelblue", pos = 2)
grid()
legend("topright",
       c("<= 1.5x", "1.5-2x", "> 2x", "Not reached"),
       fill = c("steelblue", "orange", "darkred", "gray80"),
       cex = 0.6, border = "gray40")
mtext("D", side = 3, at = -0.02, font = 2, line = 0.5)

mtext("Differential Amplitude Power Analysis: Effect of Amplitude Fold-Change", outer = TRUE, font = 2, cex = 1.1)
dev.off()

cat("Saved:", out_file, "\n")
