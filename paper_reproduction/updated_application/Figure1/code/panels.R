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


