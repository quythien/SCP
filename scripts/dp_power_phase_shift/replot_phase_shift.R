#' Replot Phase Shift Power Results - 6 Panel Figure
#' Shows relationship between phase shift magnitude, sample size, r, and power

setwd("/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")

# Load saved results
load("output/dp_power_phase_shift/dp_power_phase_shift_results.rds")

sample_sizes <- dp_phase_results$sample_sizes
phase_shifts <- dp_phase_results$phase_shifts
nsims <- dp_phase_results$nsims
target_effect <- dp_phase_results$target_effect
r_strata <- dp_phase_results$r_strata
strata_labels <- dp_phase_results$strata_labels
strat_power <- dp_phase_results$strat_power
strat_TD <- dp_phase_results$strat_TD
strat_FD <- dp_phase_results$strat_FD
strat_n_targets <- dp_phase_results$strat_n_targets

n_phase <- length(phase_shifts)
n_size <- length(sample_sizes)
n_r_strata <- length(r_strata) - 1

# Indices
idx_n60 <- which.min(abs(sample_sizes - 60))
show_phase_idx <- which(phase_shifts %in% c(2, 4, 6, 8, 10, 12))
show_phase_labels <- phase_shifts[show_phase_idx]

# Subset of r strata to show (for cleaner plots)
show_r_idx <- c(4, 6, 8, 10, 12, 14)  # (0.75,1], (1.25,1.5], (1.75,2], (2.5,3], (3.5,4], (4.5,5]
show_r_labels <- strata_labels[show_r_idx]
show_r_colors <- rainbow(6, s = 0.6, v = 0.8)

# =====================================================================
# PRE-COMPUTE ALL DATA FOR PLOTS
# =====================================================================

cat("Computing plot data...\n")

# Panel A: Power vs r by Phase Shift (at n=60)
# Dimensions: [phase_shift, r_stratum]
mean_power_n60 <- array(NA, dim = c(n_phase, n_r_strata))
for (p in 1:n_phase) {
  for (k in 1:n_r_strata) {
    mean_power_n60[p, k] <- mean(strat_power[p, idx_n60, k, ], na.rm = TRUE)
  }
}

# Panel B: Heatmap of Power (phase shift x r) at n=60
heatmap_data <- mean_power_n60

# Panel C: Power vs Phase Shift by r stratum (at n=60)
# For each r stratum, power across phase shifts
power_by_phase <- matrix(NA, nrow = n_phase, ncol = length(show_r_idx))
for (p in 1:n_phase) {
  for (j in seq_along(show_r_idx)) {
    k <- show_r_idx[j]
    power_by_phase[p, j] <- mean(strat_power[p, idx_n60, k, ], na.rm = TRUE)
  }
}

# Panel D: Power vs Sample Size by Phase Shift (marginal across r)
marginal_power_size_phase <- matrix(NA, nrow = n_size, ncol = n_phase)
for (j in 1:n_size) {
  for (p in 1:n_phase) {
    total_TD <- sum(strat_TD[p, j, , ], na.rm = TRUE)
    total_targets <- sum(strat_n_targets[p, j, , ], na.rm = TRUE)
    marginal_power_size_phase[j, p] <- if (total_targets > 0) total_TD / total_targets else NA
  }
}

# Panel E: Heatmap of Power (sample size x phase shift)
heatmap_size_phase <- marginal_power_size_phase

# Panel F: Sample size needed for 80% power by phase shift (interpolated + extrapolated)
min_n_for_80 <- numeric(n_phase)
is_extrapolated <- logical(n_phase)
for (p in 1:n_phase) {
  powers <- marginal_power_size_phase[, p]
  found <- FALSE
  is_extrapolated[p] <- FALSE

  # Check if 80% is achievable within our range
  if (any(!is.na(powers) & powers >= 0.80)) {
    # Find first crossing point and interpolate
    for (j in 1:n_size) {
      if (!is.na(powers[j]) && powers[j] >= 0.80) {
        if (j == 1) {
          min_n_for_80[p] <- sample_sizes[1]
        } else {
          # Linear interpolation between j-1 and j
          p_lo <- powers[j - 1]
          p_hi <- powers[j]
          n_lo <- sample_sizes[j - 1]
          n_hi <- sample_sizes[j]
          if (!is.na(p_lo) && p_hi > p_lo) {
            frac <- (0.80 - p_lo) / (p_hi - p_lo)
            min_n_for_80[p] <- n_lo + frac * (n_hi - n_lo)
          } else {
            min_n_for_80[p] <- sample_sizes[j]
          }
        }
        found <- TRUE
        break
      }
    }
  }
  if (!found) {
    # Extrapolate: fit power ~ 1 - exp(-b * n) model on available data
    valid <- !is.na(powers) & powers > 0 & powers < 1
    if (sum(valid) >= 3) {
      nn <- sample_sizes[valid]
      pp <- powers[valid]
      # Transform: log(1 - power) = -b * n => linear in n
      y_trans <- log(1 - pp)
      fit <- tryCatch(lm(y_trans ~ nn - 1), error = function(e) NULL)
      if (!is.null(fit) && coef(fit)[1] < 0) {
        b <- -coef(fit)[1]
        # Solve: 1 - exp(-b * n) = 0.80 => n = -log(0.20) / b
        n_needed <- log(1 / 0.20) / b
        min_n_for_80[p] <- round(n_needed)
        is_extrapolated[p] <- TRUE
      } else {
        min_n_for_80[p] <- NA
      }
    } else {
      min_n_for_80[p] <- NA
    }
  }
}

cat("Plot data computed.\n\n")

# =====================================================================
# CREATE 6-PANEL FIGURE
# =====================================================================

pdf("output/figures/dp_power_phase_shift.pdf", width = 16, height = 10)
par(mfrow = c(2, 3), mai = c(0.9, 1.0, 0.5, 0.2), mgp = c(3, 0.5, 0), oma = c(0, 0, 2, 0))

# =====================================================================
# Panel A: Power vs r by Phase Shift (n=60)
# =====================================================================

matplot(1:n_r_strata, 100 * t(mean_power_n60[show_phase_idx, ]),
        type = "b", pch = 19, lwd = 2,
        col = rainbow(length(show_phase_idx), s = 0.6, v = 0.8), lty = 1,
        xlim = c(0.5, n_r_strata + 0.5), ylim = c(0, 100),
        xlab = "r = A/sigma (Signal-to-Noise Ratio)",
        ylab = "Power (%)",
        main = "Power vs r by Phase Shift (n=60)",
        xaxt = "n")
axis(1, at = 1:n_r_strata, labels = strata_labels, las = 2, cex.axis = 0.6)
abline(h = 80, lty = 2, col = "gray")
grid()
legend("bottomright", paste0(show_phase_labels, "h"),
       col = rainbow(length(show_phase_idx), s = 0.6, v = 0.8),
       lty = 1, pch = 19, lwd = 2, cex = 0.7)
mtext("A", side = 3, at = -0.02, font = 2, line = 0.5)

# =====================================================================
# Panel B: Heatmap of Power (phase shift x r) at n=60
# =====================================================================

par(mai = c(0.9, 1.0, 0.5, 0.5))
# heatmap_data is [phase, r], need [r, phase] for image()
# Use image(z) without x/y args to avoid dimension error, then add custom axes
z_matrix <- 100 * t(heatmap_data)  # [r, phase]
image(z_matrix,
      xlab = "Phase Shift (hours)", ylab = "r = A/sigma",
      main = "Power Heatmap (n=60)",
      col = colorRampPalette(c("white", "yellow", "orange", "red", "darkred"))(100),
      xaxt = "n", yaxt = "n")
# Add custom axes
axis(1, at = seq(0, 1, length.out = n_phase), labels = phase_shifts)
axis(2, at = seq(0, 1, length.out = n_r_strata), labels = strata_labels, las = 2, cex.axis = 0.6)
box()
mtext("B", side = 3, at = -0.02, font = 2, line = 0.5)

# =====================================================================
# Panel C: Power vs Phase Shift by r stratum (n=60)
# =====================================================================

par(mai = c(0.9, 1.0, 0.5, 0.2))
matplot(phase_shifts, 100 * power_by_phase,
        type = "b", pch = 19, lwd = 2,
        col = show_r_colors, lty = 1,
        xlim = c(0, 13), ylim = c(0, 100),
        xlab = "Phase Shift (hours)",
        ylab = "Power (%)",
        main = "Power vs Phase Shift by r (n=60)")
abline(h = 80, lty = 2, col = "gray")
abline(v = 6, lty = 3, col = "darkgreen", lwd = 1.5)
grid()
legend("bottomright", show_r_labels,
       col = show_r_colors, lty = 1, pch = 19, lwd = 2, cex = 0.6)
mtext("C", side = 3, at = -0.02, font = 2, line = 0.5)

# =====================================================================
# Panel D: Power vs Sample Size by Phase Shift (marginal across r)
# =====================================================================

par(mai = c(0.9, 1.0, 0.5, 0.2))
# marginal_power_size_phase is [size, phase], need to subset columns then transpose
matplot(sample_sizes, 100 * marginal_power_size_phase[, show_phase_idx],
        type = "b", pch = 19, lwd = 2,
        col = rainbow(length(show_phase_idx), s = 0.6, v = 0.8), lty = 1,
        xlim = c(0, 110), ylim = c(0, 100),
        xlab = "Sample Size (per group)",
        ylab = "Marginal Power (%)",
        main = "Marginal Power vs Sample Size")
abline(h = 80, lty = 2, col = "gray")
abline(v = 60, lty = 3, col = "darkgreen", lwd = 1.5)
grid()
legend("bottomright", paste0(show_phase_labels, "h"),
       col = rainbow(length(show_phase_idx), s = 0.6, v = 0.8),
       lty = 1, pch = 19, lwd = 2, cex = 0.6)
mtext("D", side = 3, at = -0.02, font = 2, line = 0.5)

# =====================================================================
# Panel E: Heatmap of Power (sample size x phase shift)
# =====================================================================

par(mai = c(0.9, 1.0, 0.5, 0.5))
# heatmap_size_phase is [size, phase], need [phase, size] for image()
# Use image(z) without x/y args to avoid dimension error, then add custom axes
z_matrix2 <- 100 * t(heatmap_size_phase)  # [phase, size]
image(z_matrix2,
      xlab = "Sample Size (per group)", ylab = "Phase Shift (hours)",
      main = "Power Heatmap (Marginal)",
      col = colorRampPalette(c("white", "yellow", "orange", "red", "darkred"))(100),
      xaxt = "n", yaxt = "n")
# Add custom axes
axis(1, at = seq(0, 1, length.out = n_size), labels = sample_sizes)
axis(2, at = seq(0, 1, length.out = n_phase), labels = phase_shifts)
box()
mtext("E", side = 3, at = -0.02, font = 2, line = 0.5)

# =====================================================================
# Panel F: Sample size needed for 80% power by phase shift
# =====================================================================

par(mai = c(0.9, 1.0, 0.5, 0.2))

# Bar chart with log-scale y-axis
non_zero_phase <- phase_shifts[phase_shifts > 0]
n_for_80_nz <- min_n_for_80[phase_shifts > 0]
is_extrap_nz <- is_extrapolated[phase_shifts > 0]

# Color: blue = achievable at study n, darkred = needs more, orange = extrapolated
bar_colors <- ifelse(is.na(n_for_80_nz), "gray80",
                     ifelse(is_extrap_nz, "orange",
                            ifelse(n_for_80_nz <= 60, "steelblue", "darkred")))

# Log-scale bar heights (use log10 for display, plot manually)
log_heights <- ifelse(is.na(n_for_80_nz), 0, log10(n_for_80_nz))

# Y-axis range in log scale
y_max_log <- ceiling(max(log_heights, na.rm = TRUE)) + 0.3

bp <- barplot(log_heights, names.arg = paste0(non_zero_phase, "h"),
              col = bar_colors, border = "gray40",
              ylim = c(0, y_max_log),
              xlab = "Phase Shift (hours)",
              ylab = "",
              main = "Sample Size for 80% Marginal Power",
              yaxt = "n")

# Custom log-scale y-axis
y_ticks <- c(10, 20, 50, 100, 200, 500, 1000, 5000, 10000, 20000)
y_ticks <- y_ticks[log10(y_ticks) <= y_max_log]
axis(2, at = log10(y_ticks), labels = formatC(y_ticks, format = "d", big.mark = ","), las = 1, cex.axis = 0.7)
mtext("Min Sample Size (per group)", side = 2, line = 3.5)

# Add value labels on bars
for (i in seq_along(n_for_80_nz)) {
  if (is.na(n_for_80_nz[i])) {
    text(bp[i], 0.3, "N/A", cex = 0.65, font = 3, col = "gray40")
  } else {
    n_val <- round(n_for_80_nz[i])
    label <- formatC(n_val, format = "d", big.mark = ",")
    if (is_extrap_nz[i]) label <- paste0("~", label, "*")
    else label <- paste0("~", label)
    text(bp[i], log_heights[i] + 0.15, label, cex = 0.6, font = 2)
  }
}

# Reference lines
abline(h = log10(60), lty = 2, col = "steelblue", lwd = 1.5)
text(bp[1], log10(60) + 0.1, "n=60 (study)", cex = 0.6, col = "steelblue", pos = 4)
legend("topright",
       c("Interpolated (n<=60)", "Interpolated (n>60)", "Extrapolated*"),
       fill = c("steelblue", "darkred", "orange"), cex = 0.6, border = "gray40")
mtext("F", side = 3, at = -0.02, font = 2, line = 0.5)

# Add overall title
mtext("Differential Phase Power Analysis: Effect of Phase Shift Magnitude",
      outer = TRUE, font = 2, cex = 1.2)

dev.off()

cat("\nFigure saved: output/figures/dp_power_phase_shift.pdf\n\n")

# =====================================================================
# SUMMARY TABLE
# =====================================================================

cat("====================================================================\n")
cat("POWER SUMMARY BY PHASE SHIFT AND SAMPLE SIZE (Marginal across r)\n")
cat("====================================================================\n\n")

cat("Phase Shift | ")
for (n in sample_sizes) {
  cat(sprintf(" n=%d |", n))
}
cat("\n")
cat(paste0(rep("-", 18 + length(sample_sizes) * 7), collapse = ""), "\n")

for (p in seq_along(phase_shifts)) {
  cat(sprintf("%11s | ", phase_shifts[p]))
  for (j in seq_along(sample_sizes)) {
    cat(sprintf(" %4.1f%% |", 100 * marginal_power_size_phase[j, p]))
  }
  cat("\n")
}

cat("\n")
