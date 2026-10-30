draw_boot <- function(d, letter, title, ylab = TRUE) {
  d <- d[order(d$N), ]
  plot(NA, xlim = c(0, max(d$N) * 1.04), ylim = c(0, 100), xlab = "Sample size (N)",
       ylab = if (ylab) "Power (%)" else "", yaxt = if (ylab) "s" else "n", main = "")
  if (!ylab) axis(2, labels = FALSE)
  panel_title(letter, title)
  abline(h = 80, lty = 2, col = "grey50", lwd = 0.9); grid(lwd = 0.5)
  lines(d$N, 100 * d$plugin, lwd = 1.6, col = CP)
  lines(d$N, 100 * d$mean, lwd = 1.0, lty = 2, col = CB)
  ok <- is.finite(d$lo) & is.finite(d$hi) & (d$hi - d$lo) > 1e-6   # zero-length bars are not drawn (as shipped)
  arrows(d$N[ok], 100 * d$lo[ok], d$N[ok], 100 * d$hi[ok], code = 3, angle = 90,
         length = 0.025, lwd = 1.0, col = CB)
  points(d$N, 100 * d$mean, pch = 19, col = CB, cex = PT_CEX + 0.05)
}

