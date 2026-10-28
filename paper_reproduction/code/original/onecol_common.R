# Shared settings for the one-column (86 mm) figure versions used by the
# Bioinformatics port (examples/publication/onecol/*.R). These scripts are
# plot-only: they read the same cached results as the full-width generators
# and never run simulations. Journal figure spec: single column 86 mm, text
# at least 7 pt, Arial/Helvetica (R's pdf() default family is Helvetica).
COL_W    <- 86 / 25.4   # inches
PT       <- 8           # device point size: axis titles 8 pt
CEX_AXIS <- 0.875       # tick labels 7 pt
CEX_LAB  <- 1.0
CEX_MAIN <- 0.9         # panel titles 7.2 pt bold
CEX_LEG  <- 0.875       # legends 7 pt
LWD      <- 1.3
PT_CEX   <- 0.45
OUT_DIR  <- "figures/onecol"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

pick_latest <- function(dir, pattern) {
  cands <- list.files(dir, pattern, full.names = TRUE)
  if (length(cands) == 0) stop("No match for ", pattern, " in ", dir)
  cands[order(file.mtime(cands), decreasing = TRUE)[1]]
}

# +/- SE bars (same rule as SCP:::add_se_bars: skip NA and near-zero SE)
se_bars <- function(x, y, se, col, len = 0.02, lwd = 0.8) {
  ok <- !is.na(y) & !is.na(se) & se > 1e-6
  if (!any(ok)) return(invisible())
  arrows(x[ok], (y - se)[ok], x[ok], (y + se)[ok], angle = 90, code = 3,
         length = len, col = col, lwd = lwd)
}

panel_title <- function(letter, text, line = 0.3) {
  title(main = paste0(letter, "   ", text), adj = 0.5, font.main = 2,
        cex.main = CEX_MAIN, line = line)
}

# cex = 1 undoes the automatic 0.83 / 0.66 reduction that par(mfrow) and
# layout() apply to multi-panel pages, so sizes above are the printed sizes.
# One legend filling the current layout cell (so it lines up with the panels
# above it), without a border. The entries are shrunk until they fit the cell.
# `at`: device-NDC x on which to centre the legend (use plot_centre() of the
# panel(s) above), so the legend lines up with the plot boxes rather than
# with the cell, whose left part is taken up by the axis margin.
legend_strip <- function(labels, cols, pch = 19, ncol = length(labels),
                         title = NULL, cex = 0.8, lwd = LWD, box = FALSE, at = NULL) {
  par(mai = c(0.03, 0.03, 0.03, 0.03), xpd = NA); plot.new()
  on.exit(par(xpd = FALSE), add = TRUE)   # restore clipping for panels drawn later
  if (box) rect(0, 0, 1, 1, border = "grey60", lwd = 0.5, col = "white")
  pos <- if (is.null(at)) list(x = "center") else
    list(x = grconvertX(at, "ndc", "user"), y = 0.5, xjust = 0.5, yjust = 0.5)
  args <- c(pos, list(legend = labels, col = cols, pch = pch, lty = 1,
               lwd = lwd, ncol = ncol, bty = "n", seg.len = 1.1, x.intersp = 0.5,
               y.intersp = 0.95, pt.cex = PT_CEX, title = title, title.font = 2))
  repeat {
    lg <- do.call(legend, c(args, list(cex = cex, plot = FALSE)))
    if ((lg$rect$w <= 0.97 && lg$rect$h <= 0.97) || cex <= 0.55) break
    cex <- cex - 0.025
  }
  do.call(legend, c(args, list(cex = cex)))
  invisible(cex)
}

# One legend box under a multi-panel figure: FDR entries (line + point) and
# target-N entries (line), two rows that each read left to right.
bottom_legend <- function(fdr_labels, fdr_cols, n_labels, n_cols, cex = 0.78) {
  par(mai = c(0.02, 0.04, 0.02, 0.04)); plot.new()
  lab <- c(fdr_labels, n_labels); col <- c(fdr_cols, n_cols)
  pch <- c(rep(19, length(fdr_labels)), rep(NA, length(n_labels)))
  n <- length(lab); nc <- ceiling(n / 2)
  ord <- as.vector(rbind(seq_len(nc), nc + seq_len(nc)))   # column-major fill
  ord <- ord[ord <= n]
  legend("center", legend = lab[ord], col = col[ord], pch = pch[ord], lty = 1,
         lwd = LWD, ncol = nc, cex = cex, bty = "o", box.col = "grey60",
         box.lwd = 0.5, bg = "white", seg.len = 1.1, x.intersp = 0.5,
         y.intersp = 0.95, pt.cex = PT_CEX)
}

# device-NDC x of the centre of the current panel's plot box
plot_centre <- function() grconvertX(mean(par("usr")[1:2]), "user", "ndc")

base_par <- function(mai) {
  par(mai = mai, mgp = c(1.7, 0.45, 0), tcl = -0.25, cex = 1, las = 1, xpd = FALSE,
      cex.axis = CEX_AXIS, cex.lab = CEX_LAB, font.main = 2)
}
