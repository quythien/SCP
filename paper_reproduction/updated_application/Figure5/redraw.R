# 1. Read the Liver and Caudate-Putamen power summaries.
source("code/onecol_common.R")
curves <- read.csv("Figure2/data/curves_all_BH.csv")
ab <- readRDS("Figure5/data/liver_active.rds")
endpoints <- c("DR", "DP", "DM")

stopifnot(isTRUE(ab$balanced_B_curves_coincide), isTRUE(ab$nonidentifiable_designs_not_included))
active <- subset(curves, design == "active" & target_population == "all_assigned" & BH_threshold == .2)
B_grid <- sort(unique(active$B))
N_grid <- sort(unique(active$N_per_group))
B_all <- c(3, 4, 6, 8, 12, 24)
B_colors <- setNames(c("#440154", "#3b528b", "#21918c", "#5ec962", "#fde725", "#7d028c"), B_all)

# Restore the original display offsets and zero baselines for unidentified designs.
# These baseline values are a plotting convention, not estimated power.
K2_display <- matrix(0, length(ab$N), length(B_all))
K2_display_se <- K2_display
K2_display[, match(ab$B_K2, B_all)] <- ab$pwr_K2
K2_display_se[, match(ab$B_K2, B_all)] <- ab$se_K2
draw_panel <- function(N, M, SE, B, letter, title, xlab, ylab = TRUE) {
  plot(NA, xlim = c(0, max(N) * 1.04), ylim = c(0, 100), xlab = xlab,
       ylab = if (ylab) "Power (%)" else "", yaxt = if (ylab) "s" else "n")
  if (!ylab) axis(2, labels = FALSE)
  panel_title(letter, title)
  abline(h = 80, lty = 2, col = "grey50", lwd = .9)
  grid(lwd = .5)
  for (j in seq_along(B)) {
    jitter_step <- if (letter %in% c("A", "B")) 1.2 else 2.2
    xj <- N + (j - (length(B) + 1) / 2) * jitter_step
    y <- 100 * M[, j]; s <- 100 * SE[, j]
    ok <- is.finite(s) & s > 1e-6
    if (any(ok)) arrows(xj[ok], pmax(0, y[ok] - s[ok]), xj[ok], pmin(100, y[ok] + s[ok]),
      code = 3, angle = 90, length = .015, col = B_colors[as.character(B[j])], lwd = .7)
    lines(xj, y, type = "o", pch = 19, lwd = 1.1, cex = .35, col = B_colors[as.character(B[j])])
  }
  plot_centre()
}
pdf(file.path(OUT_DIR, "Figure5.pdf"), width = COL_W, height = 3.35, pointsize = PT)
x <- (COL_W - .36) / 3
widths <- c((x + .36)/2, (x + .36)/2, (x - .36)/2, (x + .36)/2, x/2, x/2)
layout(matrix(c(1,1,1,2,2,2, 3,3,4,4,5,5, 6,6,6,6,6,6), nrow = 3, byrow = TRUE),
       widths = widths, heights = c(1.4,1.5,.42))
base_par(c(.34,.44,.24,.05))
centers <- draw_panel(ab$N, ab$pwr_K1, ab$se_K1, ab$B_K1, "A", "Single-harmonic (K=1)", "Sample size (N)")
centers <- c(centers, draw_panel(ab$N, K2_display, K2_display_se, B_all, "B", "Two-harmonic (K=2)", "Sample size (N)"))
for (i in seq_along(endpoints)) {
  endpoint_name <- endpoints[i]
  d <- subset(active, endpoint == endpoint_name)
  mean_power <- sapply(B_grid, function(b) {
    z <- d[d$B == b, ]; z$mean_power[match(N_grid, z$N_per_group)]
  })
  mcse <- sapply(B_grid, function(b) {
    z <- d[d$B == b, ]; z$power_MCSE[match(N_grid, z$N_per_group)]
  })
  base_par(if (i == 1) c(.36,.44,.24,.03) else c(.36,.08,.24, if (i == 3) .13 else .03))
  centers <- c(centers, draw_panel(N_grid, mean_power, mcse, B_grid,
    LETTERS[i + 2], endpoint_name, "N per group", ylab = i == 1))
}
legend_strip(paste0("B=", B_all), B_colors, pch = 19, ncol = length(B_all),
             title = "Time points per cycle (B)", lwd = 1.2, at = mean(centers))
dev.off()
