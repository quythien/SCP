# Figure 1: single-cohort power and signal-to-noise strata.

# 1. Set the plotting style ------------------------------------------------
FIGURE_DIR <- "Figure1"
source("code/onecol_common.R")
VLINE_FDR <- 0.05

# 2. Read the supplied panel values ---------------------------------------
plot_data <- readRDS(file.path(FIGURE_DIR, "data", "panel_data.rds"))
adr <- plot_data$Adrenal
liv <- plot_data$Liver

# 3. Define the panels ----------------------------------------------------
draw_A <- function(d, letter) {
  matplot(d$ss_a, 100 * d$marg_mean, type = "b", pch = 19, lwd = LWD, cex = PT_CEX,
          col = d$thresh_cols, lty = 1, xlim = c(0, max(d$ss_a) * 1.05),
          ylim = c(0, 100), xlab = "Sample size (N)", ylab = "", main = "",
          yaxt = "n")
  axis(2, at = seq(0, 100, 20))
  title(ylab = "Power (%)", line = YLAB_LINE)
  panel_title(letter, "Power vs sample size")
  for (t in seq_along(d$thresh_cols))
    se_bars(d$ss_a, 100 * d$marg_mean[, t], 100 * d$marg_se[, t], d$thresh_cols[t])
  abline(h = 80, lty = 2, col = "grey50", lwd = 0.9)
  grid(lwd = 0.5)
  if (!is.na(d$vline_n)) {
    abline(v = d$vline_n, lty = 2, col = adjustcolor("steelblue", 0.7), lwd = 1.0)
    right <- d$vline_n <= 0.6 * max(d$ss_a)      # label on the roomier side of the line
    text(d$vline_n, if (right) 35 else 10, sprintf("N=%d", d$vline_n), col = "steelblue",
         cex = CEX_LEG, adj = c(if (right) -0.1 else 1.1, 0.5), font = 2)
  }
}

strata_axis <- function(d) {
  # ticks at every bin, labels at every other bin (the last bin always
  # labelled). Labelled positions are passed alone: R drops labels it judges
  # overlapping, and with blank alternates it can drop the visible ones.
  k <- length(d$strata)
  axis(1, at = seq_len(k), labels = d$strata, las = 2, cex.axis = CEX_AXIS)
  title(xlab = expression(r == A/sigma), line = 4.0)
}

draw_B <- function(d, letter, ylab = TRUE) {
  k <- length(d$strata)
  matplot(seq_len(k), 100 * t(d$pow), type = "l", lwd = LWD, col = d$size_cols,
          lty = 1, xlim = c(0.5, k + 0.5), ylim = c(0, 100),
          xaxt = "n", yaxt = if (ylab) "s" else "n", xlab = "", ylab = "", main = "")
  if (!ylab) axis(2, labels = FALSE) else title(ylab = "Power (%)", line = YLAB_LINE)
  panel_title(letter, "Stratified power")
  strata_axis(d)
  grid(lwd = 0.5)
  for (j in seq_len(nrow(d$pow))) {
    points(seq_len(k), 100 * d$pow[j, ], pch = 19, col = d$size_cols[j], cex = PT_CEX)
    se_bars(seq_len(k), 100 * d$pow[j, ], 100 * d$pow_se[j, ], d$size_cols[j])
  }
}

draw_C <- function(d, letter) {
  k <- length(d$strata)
  plot(seq_len(k), rep(0, k), type = "n", xaxt = "n",
       xlim = c(0.5, k + 0.5), ylim = c(0, max(d$gene_counts) * 1.15),
       xlab = "", ylab = "", main = "")
  title(ylab = "# True discoveries", line = YLAB_LINE)
  panel_title(letter, "True discoveries")
  strata_axis(d)
  grid(lwd = 0.5)
  step_x <- rep(seq(0.5, k + 0.5, by = 1), each = 2)
  step_y <- c(0, rep(d$gene_counts, each = 2), 0)
  polygon(step_x, step_y, col = "#cccccc55", border = NA)
  lines(step_x, step_y, col = "grey60", lwd = 1.0, lty = 2)
  for (j in seq_len(nrow(d$td))) {
    lines(seq_len(k), d$td[j, ], col = d$size_cols[j], lwd = LWD)
    points(seq_len(k), d$td[j, ], pch = 19, col = d$size_cols[j], cex = PT_CEX)
    se_bars(seq_len(k), d$td[j, ], d$td_se[j, ], d$size_cols[j])
  }
}


# 4. Draw the figure ------------------------------------------------------
OUT <- file.path(OUT_DIR, "Figure1.pdf")
stopifnot(identical(adr$ss_bc, liv$ss_bc), identical(adr$size_cols, liv$size_cols))
# 3 x 2 grid: Adrenal Gland A-C in the left column, Liver D-F in the right
# column; FDR legend under row A/D, target-N legend under the last row.
pdf(OUT, width = COL_W, height = 5.0, pointsize = PT)
layout(matrix(c(1, 2, 3, 3, 4, 5, 6, 7, 8, 8), nrow = 5, byrow = TRUE),   # cells numbered in drawing order
       heights = c(1.24, 0.22, 1.6, 1.6, 0.34))
par(oma = c(0, 0, 1.3, 0))
YLAB_LINE <- 2.5                       # y titles aligned across rows (row C has 4-digit ticks)
LEFT <- 0.54                           # one left margin for every row: equal plot widths
base_par(c(0.34, LEFT, 0.20, 0.05))
draw_A(adr, "A"); cx_adr <- plot_centre()   # plot-box centre of A
draw_A(liv, "D"); cx_liv <- plot_centre()   # plot-box centre of D
cx_mid <- mean(c(cx_adr, cx_liv))           # legends centred on the pair of columns
legend_strip(adr$thresh_labels, adr$thresh_cols, pch = 19, ncol = 4, at = cx_mid)      # A, D
base_par(c(0.66, LEFT, 0.20, 0.05))
draw_B(adr, "B"); draw_B(liv, "E")
draw_C(adr, "C"); draw_C(liv, "F")
legend_strip(paste0("N=", adr$ss_bc), adr$size_cols, pch = NA, ncol = 4, at = cx_mid)  # B, C, E, F
# column headers centred on the plot areas of A and D (device NDC)
mtext("GTEx Adrenal Gland", side = 3, outer = TRUE, at = cx_adr, line = 0.1, font = 2, cex = 1)
mtext("GTEx Liver",         side = 3, outer = TRUE, at = cx_liv, line = 0.1, font = 2, cex = 1)
dev.off()
cat(sprintf("Saved %s (N80 at FDR %g%%: Adrenal %s, Liver %s)\n", OUT, 100 * VLINE_FDR,
            adr$vline_n, liv$vline_n))
