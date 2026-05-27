# =============================================================================
# _pub_style.R
#
# Shared publication-style helpers for SCP figures. Targets Bioinformatics /
# Genome Biology look but designed to read cleanly at higher-bar venues too.
#
# Conventions:
#   - Vector PDF via cairo_pdf
#   - Helvetica-equivalent sans-serif via family = "sans" (default cairo)
#   - Panel letters: bold lowercase ("a", "b", ...) at top-left, no parens
#   - No in-figure figure-level titles; titles go in the manuscript caption
#   - Subtitles per panel carry only a topic word, not parameter values
#   - Detector palette: teal = DCP, orange = K-harmonic F-test (colorblind-safe)
#   - Sequential palette: viridis (perceptually uniform)
# =============================================================================

# Sets par() to publication defaults. Call once after opening the device.
# mfrow:  panel grid c(nrow, ncol)
# outer:  outer margins (typically c(0,0,0,0) since we drop figure titles)
pub_par <- function(mfrow = c(1, 1), mar = c(4.2, 4.2, 2.2, 1.0),
                    oma = c(0, 0, 0, 0)) {
  par(mfrow = mfrow,
      mar   = mar,
      oma   = oma,
      las   = 1,
      tcl   = -0.3,       # tick length inward, modest
      mgp   = c(2.5, 0.6, 0),  # axis title / tick label / line offsets
      family = "Helvetica",    # fontconfig-resolves to Nimbus Sans
                                # (Adobe-clone, metric-compatible with
                                # Helvetica), so PDF can be reprocessed
                                # with the proper font by the publisher.
      cex.lab  = 1.0,
      cex.axis = 0.9,
      cex.main = 1.0,
      font.main = 2)
}

# Adds a bold lowercase panel letter in the top-left corner of the current plot.
# Use immediately after plot() / before any data is drawn so it sits in the
# title area.
panel_label <- function(letter, line = 0.6, cex = 1.25) {
  mtext(letter, side = 3, line = line, adj = -0.18,
        font = 2, cex = cex)
}

# Two-color detector palette. Colorblind-safe (Wong / Okabe-Ito derived):
#   teal/cyan = DCP (cosinor)
#   orange    = K-harmonic F-test (K = 2 default)
pub_palette_detector <- function() {
  c(DCP = "#0072B2", FMM = "#D55E00")  # Wong blue, Wong vermilion
}

# Sequential palette for gradient-coded factors (e.g., B-grid, beta-grid).
# Viridis-derived hex codes; perceptually uniform.
pub_palette_sequential <- function(n) {
  full <- c("#440154", "#3b528b", "#21918c", "#5ec962", "#fde725",
            "#9e2a2b", "#ef8a62")
  if (n <= length(full)) full[seq_len(n)]
  else full[round(seq(1, length(full), length.out = n))]
}

# 80% power reference line. Light grey, dotted, low visual weight.
abline_80pct <- function(power = 0.80) {
  abline(h = power, lty = 3, col = "grey60", lwd = 0.8)
}

# Standard publication legend wrapper: no box, sans-serif, slightly smaller.
pub_legend <- function(x = "bottomright", legend, col, lwd = 1.8, lty = 1,
                       pch = NA, title = NULL, cex = 0.85) {
  legend(x, legend = legend, col = col, lwd = lwd, lty = lty, pch = pch,
         title = title, bty = "n", cex = cex)
}

# Build a legend-ready expression vector from a mixed list of plain-text
# strings and bquote-style math expressions. R's c() drops the expression
# class when mixing types; do.call(expression, list(...)) preserves it.
# Usage:
#   labs <- pub_exprs(bquote(median == .(0.28)),
#                     bquote(Beta(1, hat(beta))))
#   pub_legend("topright", legend = labs, ...)
pub_exprs <- function(...) {
  do.call(expression, list(...))
}
