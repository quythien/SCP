#!/usr/bin/env Rscript
# =============================================================================
# Fig 4 - Two-harmonic motivation on GTEx Liver  (rebuild)
#
# Row 1 (4 panels): exemplar genes ARNTL / RORC / HADHA / TPP1
#   - x-axis: TOD shifted to [-6, 18] convention
#   - main title (centered): gene name
#   - subtitle: detection mode + q-values
#   - vertical dotted lines at numerical K=2 fit peaks; "Peaks at h, h" annotation
#
# Row 2 (2 panels):
#   - Venn at BH-FDR < 0.05 on Liver (K=1 / K=2 labels, no "(DCP)")
#   - KEGG dot plot with separate count + -log10(q) gradient legends
#
# Pilot: GTEx Liver  |  Output: submission/figures/Fig4_twoharm_demo.pdf
# =============================================================================

setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")

GTEX_PATH <- "/home/qtp1/Projects/Collaborative/GTEXdata/data/v10/cpm/v10_CPM_full.RData"
FITS_PATH <- "output/two_harmonic/results/tissue_Liver_fits.rds"
VENN_PATH <- "output/two_harmonic/results/tissue_Liver_venn.rds"
KEGG_PATH <- "output/two_harmonic/results/tissue_Liver_kegg.rds"
FIG_PATH  <- "submission/figures/Fig4_twoharm_demo.pdf"

PERIOD  <- 24
OMEGA0  <- 2 * pi / PERIOD
TOD_MIN <- -6
TOD_MAX <- 18

# ---- Load Liver expression matrix ----
env <- new.env(); load(GTEX_PATH, envir = env)
liv  <- env$CPM.all[["Liver"]]
lcpm <- liv$counts
ann  <- liv$gene_annotation
tod_raw <- as.numeric(liv$tod) %% PERIOD
ok   <- !is.na(tod_raw)
lcpm <- lcpm[, ok]; tod_raw <- tod_raw[ok]
if (!is.null(ann) && "symbol" %in% colnames(ann))
  rownames(lcpm) <- ann$symbol[match(rownames(lcpm), ann$gene_id)]
first <- !duplicated(rownames(lcpm)) & nzchar(rownames(lcpm)) &
         !is.na(rownames(lcpm))
lcpm <- lcpm[first, ]
keep <- rowMeans(lcpm) > 3
lcpm <- lcpm[keep, ]

# Shift TOD to the [-6, 18] convention.
shift_tod <- function(t) ((t + 6) %% PERIOD) - 6
tod_plot  <- shift_tod(tod_raw)
N         <- ncol(lcpm)

cat(sprintf("Liver matrix: %d genes x %d samples\n", nrow(lcpm), N))

# ---- Load cached fits / venn / kegg ----
fits <- readRDS(FITS_PATH)
venn <- readRDS(VENN_PATH)
ek   <- readRDS(KEGG_PATH)

# Design matrices use raw TOD (period-respecting).
X1 <- cbind(1, cos(OMEGA0 * tod_raw), sin(OMEGA0 * tod_raw))
X2 <- cbind(X1, cos(2 * OMEGA0 * tod_raw), sin(2 * OMEGA0 * tod_raw))

# ---- Numerical peak finder for K=2 curve over one period ----
find_peaks_2H <- function(b2, n = 2000) {
  # Evaluate over [0, 24] then wrap to [-6, 18]
  tt <- seq(0, PERIOD, length.out = n + 1)[-(n + 1)]
  y  <- b2[1] +
        b2[2] * cos(OMEGA0 * tt) + b2[3] * sin(OMEGA0 * tt) +
        b2[4] * cos(2 * OMEGA0 * tt) + b2[5] * sin(2 * OMEGA0 * tt)
  # Circular local maxima (use wrap-around neighbours)
  prev <- c(y[n], y[-n])
  nxt  <- c(y[-1], y[1])
  is_pk <- y > prev & y > nxt
  peak_t <- tt[is_pk]; peak_v <- y[is_pk]
  ord <- order(peak_v, decreasing = TRUE)
  list(t = shift_tod(peak_t[ord]), v = peak_v[ord])
}

# ---- Helper: render one gene exemplar panel ----
exemplar_genes <- list(
  c("ARNTL", "both"),
  c("RORC",  "K=2 only"),
  c("IK",    "K=2 only"),    # Hughes 2009 Table S2 verified concordance
  c("GBA",   "K=2 only")     # Lysosome biogenesis KEGG enrichment driver
)

render_gene <- function(gene, status, letter_label) {
  gi <- which(rownames(lcpm) == gene)[1]
  y  <- as.numeric(lcpm[gi, ])
  b1 <- lm.fit(X1, y)$coefficients
  b2 <- lm.fit(X2, y)$coefficients
  row_q <- fits[fits$gene == gene, , drop = FALSE]
  q1 <- row_q$q1H[1]; q2 <- row_q$q2H[1]

  # Continuous curves over the shifted axis
  tt_orig <- seq(0, PERIOD, length.out = 400)
  tt_plot <- shift_tod(tt_orig)
  c1 <- b1[1] + b1[2]*cos(OMEGA0*tt_orig) + b1[3]*sin(OMEGA0*tt_orig)
  c2 <- b2[1] + b2[2]*cos(OMEGA0*tt_orig) + b2[3]*sin(OMEGA0*tt_orig) +
                b2[4]*cos(2*OMEGA0*tt_orig) + b2[5]*sin(2*OMEGA0*tt_orig)

  # Sort by plot-x so polyline draws cleanly across the discontinuity
  ord <- order(tt_plot)
  tt_plot <- tt_plot[ord]; c1 <- c1[ord]; c2 <- c2[ord]

  ylim <- range(c(y, c1, c2), na.rm = TRUE)

  plot(tod_plot, y, pch = 19, col = adjustcolor("grey50", 0.5), cex = 0.55,
       xlim = c(TOD_MIN, TOD_MAX), ylim = ylim,
       xlab = "Time of day (h)", ylab = expression(log[2](CPM+1)),
       xaxt = "n", main = "")
  lines(tt_plot, c1, col = "red",  lwd = 2.2, lty = 2)
  lines(tt_plot, c2, col = "blue", lwd = 2.2, lty = 1)
  axis(1, at = seq(TOD_MIN, TOD_MAX, 6))

  # Peak markers from K=2 fit; in-plot label positioned lower (~85% of plot height)
  pks <- find_peaks_2H(b2)
  in_range <- pks$t >= TOD_MIN & pks$t < TOD_MAX
  pk_t <- pks$t[in_range]; pk_v <- pks$v[in_range]
  usr_y <- par("usr")[3:4]
  lbl_y <- usr_y[1] + (usr_y[2] - usr_y[1]) * 0.85
  for (k in seq_along(pk_t)) {
    abline(v = pk_t[k], lty = 3, col = "blue", lwd = 1.2)
    text(pk_t[k], lbl_y,
         sprintf("%.1fh", pk_t[k]),
         adj = c(-0.1, 0.5), cex = 0.72, col = "blue", xpd = NA)
  }
  box()

  # Centered main title = gene name (bold); subtitle further below for spacing
  title(main = gene, adj = 0.5, font.main = 2, cex.main = 1.30, line = 1.55)
  # Single subtitle line: q-values only (no status label; that's in caption)
  mtext(bquote(q[1*H] == .(formatC(q1, format="g", digits=2)) ~~ "  " ~~
               q[2*H] == .(formatC(q2, format="g", digits=2))),
        side = 3, line = 0.25, adj = 0.5, cex = 0.78, col = "grey20")
  # Optional sub-panel letter (only printed when letter_label != "")
  if (nzchar(letter_label)) {
    mtext(letter_label, side = 3, line = 1.95, at = par("usr")[1],
          adj = 0, font = 2, cex = 1.25, col = "grey20")
  }
}

# ---- Helper: Venn diagram ----
draw_venn_panel <- function(n_k1only, n_both, n_k2only, n_total) {
  par(mar = c(3.2, 1, 4.0, 1))
  plot.new(); plot.window(xlim = c(-2.0, 2.0), ylim = c(-1.3, 1.3), asp = 1)
  theta <- seq(0, 2*pi, length.out = 200)
  r <- 0.95
  c1 <- cbind(-0.55 + r*cos(theta),  r*sin(theta))
  c2 <- cbind( 0.55 + r*cos(theta),  r*sin(theta))
  polygon(c1, col = adjustcolor("#d62728", 0.45), border = "#d62728", lwd = 1.5)
  polygon(c2, col = adjustcolor("#1f77b4", 0.45), border = "#1f77b4", lwd = 1.5)
  text(-1.10, -0.05, format(n_k1only, big.mark = ","), cex = 1.55, font = 2)
  text( 1.10, -0.05, format(n_k2only, big.mark = ","), cex = 1.55, font = 2)
  text( 0.00, -0.05, format(n_both,   big.mark = ","), cex = 1.55, font = 2)
  text(-0.85,  1.05, "K=1", col = "#d62728", cex = 1.25, font = 2)
  text( 0.85,  1.05, "K=2", col = "#1f77b4", cex = 1.25, font = 2)
  text(0, -1.20,
       sprintf("n = %s genes", format(n_total, big.mark = ",")),
       cex = 0.95, col = "grey25")
  title(main = "Rhythmicity Biomarker overlap in GTEx Liver (FDR 5%)",
        adj = 0.5, font.main = 2, cex.main = 1.30, line = 1.4)
  mtext("B", side = 3, line = 2.4, at = par("usr")[1],
        adj = 0, font = 2, cex = 1.05, col = "grey20")
}

# ---- Helper: KEGG dot plot with proper legends ----
draw_kegg_panel <- function(ek, max_terms = 10) {
  par(mar = c(4.5, 15.0, 4.5, 7.0))    # wider left margin for pathway labels; wider right for legends
  df_e <- as.data.frame(ek)
  # Drop KEGG "Human Diseases" category. These pathway gene sets share genes
  # with metabolism/proteostasis/immunity categories already represented in
  # the result, so the underlying biology is preserved without disease labels.
  if ("category" %in% colnames(df_e))
    df_e <- df_e[df_e$category != "Human Diseases", , drop = FALSE]
  df_e <- df_e[order(df_e$p.adjust), ]
  df_e <- head(df_e, max_terms)
  if (nrow(df_e) == 0L) {
    plot.new(); title("KEGG: no terms enriched"); return(invisible())
  }
  gene_ratio_num <- sapply(strsplit(df_e$GeneRatio, "/"),
                            function(z) as.integer(z[1]))
  gene_ratio_den <- sapply(strsplit(df_e$GeneRatio, "/"),
                            function(z) as.integer(z[2]))
  gene_ratio <- gene_ratio_num / gene_ratio_den

  # Reorder so smallest q at top (largest -log10(q) at top)
  ord <- order(df_e$p.adjust, decreasing = TRUE)
  df_e <- df_e[ord, ]; gene_ratio <- gene_ratio[ord]

  qadj <- df_e$p.adjust
  size <- df_e$Count
  z <- -log10(qadj)

  # Color gradient: blue (low sig) -> orange (high sig)
  col_pal <- colorRampPalette(c("#2c7fb8", "#7fcdbb", "#ffeda0", "#feb24c", "#f03b20"))(60)
  z_rng <- range(z); if (diff(z_rng) < 1e-9) z_rng[2] <- z_rng[1] + 1
  col_idx <- pmin(60, pmax(1, ceiling(60 * (z - z_rng[1]) / diff(z_rng))))
  cols <- col_pal[col_idx]

  # Size scaling (capped for legend readability)
  size_rng <- range(size)
  cex_pts <- 0.9 + 3.0 * (size - size_rng[1]) / (diff(size_rng) + 1e-9)

  plot(gene_ratio, seq_len(nrow(df_e)),
       pch = 19, cex = cex_pts, col = cols,
       xlim = c(0, max(gene_ratio) * 1.10),
       ylim = c(0.5, nrow(df_e) + 0.5),
       xlab = "Gene ratio (K=2-only set)",
       ylab = "", yaxt = "n", main = "")
  axis(2, at = seq_len(nrow(df_e)), labels = df_e$Description,
       las = 1, cex.axis = 1.00)
  grid(nx = NULL, ny = NA, lty = "dotted", col = "grey80")
  title(main = "KEGG enrichment (K=2-only)",
        adj = 0.5, font.main = 2, cex.main = 1.35, line = 1.7)
  mtext("C", side = 3, line = 2.7, at = par("usr")[1],
        adj = 0, font = 2, cex = 1.10, col = "grey20")

  # --- Legends in right margin (more space + clear separation) ---
  usr <- par("usr")
  plot_w <- usr[2] - usr[1]
  # Anchor legends in the right margin with explicit horizontal offset
  leg_x <- usr[2] + plot_w * 0.13

  # ---- Size legend (top of right margin) ----
  size_breaks <- pretty(size, 4)
  size_breaks <- size_breaks[size_breaks > 0]
  if (length(size_breaks) > 4)
    size_breaks <- size_breaks[seq(1, length(size_breaks), length.out = 4)]
  size_cex_leg <- 0.55 + 1.5 * (size_breaks - size_rng[1]) /
                  (diff(size_rng) + 1e-9)
  n_sz <- length(size_breaks)
  # Compact count legend; leaves clear visual gap before the gradient bar
  legend_top_y    <- usr[4] - (usr[4]-usr[3]) * 0.05
  legend_bot_y    <- usr[4] - (usr[4]-usr[3]) * 0.28
  sz_y_positions  <- seq(legend_top_y, legend_bot_y, length.out = n_sz + 1)
  text(leg_x, sz_y_positions[1], "Count",
       cex = 0.85, font = 2, xpd = NA, adj = 0.5)
  for (k in seq_along(size_breaks)) {
    yy <- sz_y_positions[k + 1]
    points(leg_x - plot_w * 0.02, yy,
           pch = 19, cex = size_cex_leg[k], col = "grey50", xpd = NA)
    text(leg_x + plot_w * 0.025, yy,
         size_breaks[k], cex = 0.72, adj = 0, xpd = NA)
  }

  # ---- Color gradient bar (bottom of right margin) ----
  bar_x  <- leg_x - plot_w * 0.022
  bar_w  <- plot_w * 0.030
  bar_y_high <- usr[3] + (usr[4]-usr[3]) * 0.45
  bar_y_low  <- usr[3] + (usr[4]-usr[3]) * 0.05
  n_seg <- 60
  ys <- seq(bar_y_low, bar_y_high, length.out = n_seg + 1)
  for (j in seq_len(n_seg))
    rect(bar_x, ys[j], bar_x + bar_w, ys[j + 1],
         col = col_pal[j], border = NA, xpd = NA)
  rect(bar_x, bar_y_low, bar_x + bar_w, bar_y_high,
       col = NA, border = "grey50", xpd = NA)
  z_tick_vals <- pretty(z, 4)
  z_tick_vals <- z_tick_vals[z_tick_vals >= z_rng[1] & z_tick_vals <= z_rng[2]]
  if (length(z_tick_vals) > 4)
    z_tick_vals <- z_tick_vals[seq(1, length(z_tick_vals), length.out = 4)]
  for (zt in z_tick_vals) {
    y_zt <- bar_y_low + (zt - z_rng[1]) / diff(z_rng) * (bar_y_high - bar_y_low)
    segments(bar_x + bar_w, y_zt, bar_x + bar_w + plot_w * 0.005, y_zt, xpd = NA)
    text(bar_x + bar_w + plot_w * 0.012, y_zt,
         formatC(zt, format = "g", digits = 2),
         adj = 0, cex = 0.72, xpd = NA)
  }
  text(bar_x + bar_w/2, bar_y_high + (usr[4]-usr[3]) * 0.05,
       expression(-log[10](q[adj])),
       cex = 0.78, font = 2, xpd = NA, adj = 0.5)
}

# ---- Layout + render ----
# 3 rows now: row 1 = gene exemplars; row 2 = horizontal legend strip; row 3 = Venn + KEGG
pdf(FIG_PATH, width = 12.0, height = 9.6)
lay <- rbind(
  c(1, 2, 3, 4),
  c(5, 5, 5, 5),
  c(6, 6, 7, 7)
)
layout(lay, heights = c(1.00, 0.12, 1.35))
par(mai = c(0.70, 0.78, 0.70, 0.20),
    mgp = c(2.5, 0.6, 0), oma = c(0.4, 0.4, 2.4, 0.4),
    cex.axis = 1.05, cex.lab = 1.15, font.main = 2)

# Row 1: gene exemplars (panels 1-4). Only the first panel carries the "A" label.
panel_letters <- c("A", "", "", "")
for (i in seq_along(exemplar_genes)) {
  g <- exemplar_genes[[i]]
  render_gene(g[1], g[2], panel_letters[i])
}

# Row 2: condensed horizontal legend strip (panel 5), boxed and centered
par(mar = c(0, 0, 0, 0))
plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1))
legend("center",
       legend = c("Raw expression",
                  "Single-harmonic (1H)",
                  "Two-harmonic (2H)",
                  "2H peak times"),
       col    = c(adjustcolor("grey50", 0.6), "red", "blue", "blue"),
       pch    = c(19, NA, NA, NA),
       lty    = c(NA, 2, 1, 3),
       lwd    = c(NA, 2.0, 2.0, 1.0),
       horiz  = TRUE, bty = "o", box.col = "grey70", box.lwd = 0.5,
       bg = "white", cex = 0.95, xpd = NA,
       seg.len = 0.9, text.width = 0.10,
       x.intersp = 0.2)

# Row 3: Venn (panel 6) + KEGG (panel 7)
draw_venn_panel(length(venn$k1_only), length(venn$both),
                length(venn$k2_only), venn$n_total)
draw_kegg_panel(ek, max_terms = 10)

# Outer figure title (positioned closer to the gene panels)
mtext("Single-harmonic (1H) vs. Two-harmonic (2H) cosinor fits on GTEx Liver exemplar genes",
      outer = TRUE, side = 3, line = 0.3, adj = 0.5,
      font = 2, cex = 1.20)

dev.off()
cat(sprintf("Wrote %s (%.1f KB)\n", FIG_PATH, file.size(FIG_PATH)/1024))
cat(sprintf("\nVenn @ BH-FDR < 0.05 (Liver): K1=%d, both=%d, K2-only=%d  (total %d)\n",
            length(venn$k1_only), length(venn$both),
            length(venn$k2_only), venn$n_total))

# Also print peak times for the 4 exemplars
cat("\nNumerical K=2 peak times (h) for Fig 4 Panel A genes:\n")
for (g in exemplar_genes) {
  gi <- which(rownames(lcpm) == g[1])[1]
  y <- as.numeric(lcpm[gi, ])
  b2 <- lm.fit(X2, y)$coefficients
  pks <- find_peaks_2H(b2)
  cat(sprintf("  %-7s  peaks at  %s\n",
              g[1], paste(sprintf("%.1f", pks$t), collapse = ", ")))
}
