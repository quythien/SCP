# 1. Read the completed differential experiments.
source("code/onecol_common.R")
data_dir <- Sys.getenv("SCP_DIFFERENTIAL_DATA", "Figure2/data")
curves <- read.csv(file.path(data_dir,"curves_all_BH.csv"))
guides <- read.csv(file.path(data_dir,"N80_all_BH.csv"))

# Figure 2.2 uses all assigned targets and the four original FDR thresholds.
endpoints <- c("DR", "DP", "DM")
thresholds <- c(.01, .05, .10, .20)
colors <- c("darkgreen", "steelblue", "orange", "red")
passive <- subset(curves, design == "passive" & target_population == "all_assigned")
N <- sort(unique(passive$N_per_group))
pdf(file.path(OUT_DIR, "Figure2.pdf"), width = COL_W, height = 1.76, pointsize = PT)
layout(matrix(c(1, 2, 3, 4, 4, 4), nrow = 2, byrow = TRUE),
       widths = c(1.38, 1, 1), heights = c(1.48, .27))
centers <- numeric(0)
for (i in seq_along(endpoints)) {
  endpoint_name <- endpoints[i]
  d <- subset(passive, endpoint == endpoint_name)
  mean_power <- sapply(thresholds, function(q) {
    z <- d[d$BH_threshold == q, ]; z$mean_power[match(N, z$N_per_group)]
  })
  mcse <- sapply(thresholds, function(q) {
    z <- d[d$BH_threshold == q, ]; z$power_MCSE[match(N, z$N_per_group)]
  })
  base_par(if (i == 1) c(.36, .46, .24, .04) else c(.36, .08, .24, .04))
  matplot(N, 100 * mean_power, type = "b", pch = 19, lwd = LWD, cex = PT_CEX,
          col = colors, lty = 1, xlim = c(0, max(N) * 1.05), ylim = c(0, 100),
          xlab = "N per group", ylab = if (i == 1) "Power (%)" else "",
          yaxt = if (i == 1) "s" else "n")
  if (i != 1) axis(2, labels = FALSE)
  panel_title(LETTERS[i], endpoint_name)
  for (j in seq_along(thresholds)) se_bars(N, 100 * mean_power[, j], 100 * mcse[, j], colors[j])
  abline(h = 80, lty = 2, col = "grey50", lwd = .9)
  grid(lwd = .5)
  n80 <- subset(guides, design == "passive" & endpoint == endpoint_name &
                  target_population == "all_assigned" & BH_threshold == .2)$ceil_N80
  stopifnot(length(n80) == 1)
  abline(v = n80, lty = 2, col = adjustcolor("steelblue", .7), lwd = 1)
  label_y <- if (endpoint_name %in% c("DR", "DM")) 60 else 92
  label_x <- if (endpoint_name %in% c("DR", "DM")) n80 + 20 else n80
  label_adj <- if (endpoint_name %in% c("DR", "DM")) c(0, .5) else c(1.08, .5)
  text(label_x, label_y, sprintf("N=%d", n80), col = "steelblue", cex = CEX_LEG,
       adj = label_adj, font = 2)
  centers <- c(centers, plot_centre())
}
legend_strip(paste0("FDR ", thresholds * 100, "%"), colors, pch = 19, ncol = 4,
             at = mean(centers))
dev.off()

