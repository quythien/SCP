#!/usr/bin/env Rscript
# =============================================================================
# Fig 4 - Two-harmonic motivation on GTEx Liver  (rebuild)
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

  plot(tod_plot, y, pch = 19, col = adjustcolor("grey50", 0.5), cex = 0.62,
       xlim = c(TOD_MIN, TOD_MAX), ylim = ylim,
       xlab = "Time of day (h)", ylab = expression(log[2](CPM+1)),
       xaxt = "n", main = "")
  lines(tt_plot, c1, col = "red",  lwd = 2.9, lty = 2)
  lines(tt_plot, c2, col = "blue", lwd = 2.9, lty = 1)
  axis(1, at = seq(TOD_MIN, TOD_MAX, 6))

  # Peak markers from K=2 fit
  pks <- find_peaks_2H(b2)
  in_range <- pks$t >= TOD_MIN & pks$t < TOD_MAX
  pk_t <- pks$t[in_range]; pk_v <- pks$v[in_range]
  usr_y <- par("usr")[3:4]
  lbl_y <- usr_y[1] + (usr_y[2] - usr_y[1]) * 0.85
  for (k in seq_along(pk_t)) {
    abline(v = pk_t[k], lty = 3, col = "blue", lwd = 1.2)
    text(pk_t[k], lbl_y,
         sprintf("%.1fh", pk_t[k]),
         adj = c(-0.1, 0.5), cex = 0.92, col = "blue", xpd = NA)
  }
  box()

  title(main = gene, adj = 0.5, font.main = 2, cex.main = 1.65, line = 2.05)
  mtext(bquote(q[1*H] == .(formatC(q1, format="g", digits=2)) ~~ " " ~~
               q[2*H] == .(formatC(q2, format="g", digits=2))),
        side = 3, line = 0.35, adj = 0.5, cex = 0.90, col = "grey20")
  if (nzchar(letter_label)) {
    mtext(letter_label, side = 3, line = 2.05, at = par("usr")[1],
          adj = 0, font = 2, cex = 1.55, col = "grey20")
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
  text(-1.10, -0.05, format(n_k1only, big.mark = ","), cex = 1.95, font = 2)
  text( 1.10, -0.05, format(n_k2only, big.mark = ","), cex = 1.95, font = 2)
  text( 0.00, -0.05, format(n_both,   big.mark = ","), cex = 1.95, font = 2)
  text(-0.85,  1.05, "K=1", col = "#d62728", cex = 1.55, font = 2)
  text( 0.85,  1.05, "K=2", col = "#1f77b4", cex = 1.55, font = 2)
  text(0, -1.20,
       sprintf("n = %s genes", format(n_total, big.mark = ",")),
       cex = 1.20, col = "grey25")
  title(main = "Rhythmicity Biomarker overlap in GTEx Liver (FDR 5%)",
        adj = 0.5, font.main = 2, cex.main = 1.52, line = 1.5)
  mtext("B", side = 3, line = 1.5, at = par("usr")[1],
        adj = 0, font = 2, cex = 1.55, col = "grey20")
}

# ---- Helper: KEGG enrichment as a horizontal -log10(P) bar chart ----
# Metascape convention: bars encode statistical significance as -log10(P)
# (raw enrichment p-value), pathway names sit on the y-axis. This mirrors the
# SCP Shiny app's enrichment view; the dot-plot gene-ratio/count/colour
# encoding is dropped in favour of a single, directly readable significance bar.
draw_kegg_panel <- function(ek, max_terms = 10) {
  par(mar = c(5.2, 17.5, 4.6, 2.4))
  df_e <- as.data.frame(ek)
  # Drop KEGG "Human Diseases" category. These pathway gene sets share genes
  # with metabolism/proteostasis/immunity categories already represented in
  # the result, so the underlying biology is preserved without disease labels.
  if ("category" %in% colnames(df_e))
    df_e <- df_e[df_e$category != "Human Diseases", , drop = FALSE]

  # Significance = -log10(raw enrichment p-value). Fall back to the adjusted
  # p-value only if the raw column is unavailable.
  pcol <- if ("pvalue" %in% colnames(df_e)) df_e$pvalue else df_e$p.adjust
  df_e$.p <- pcol
  df_e <- df_e[order(df_e$.p), ]          # most significant first
  df_e <- head(df_e, max_terms)
  if (nrow(df_e) == 0L) {
    plot.new(); title("KEGG: no terms enriched"); return(invisible())
  }

  # barplot() draws the first entry at the BOTTOM; reverse so the most
  # significant pathway sits at the TOP of the panel.
  df_e <- df_e[rev(seq_len(nrow(df_e))), ]
  val  <- -log10(df_e$.p)
  labs <- df_e$Description

  pal <- colorRampPalette(c("#fee0b6", "#f1a340", "#b35806"))(nrow(df_e))

  bp <- barplot(val, horiz = TRUE, names.arg = labs, las = 1,
                col = pal, border = NA,
                xlab = expression(-log[10]("P-value")),
                cex.names = 1.22, cex.lab = 1.45, cex.axis = 1.22,
                xlim = c(0, max(val) * 1.08), main = "")
  grid(nx = NULL, ny = NA, lty = "dotted", col = "grey80")
  barplot(val, horiz = TRUE, col = pal, border = NA, add = TRUE, axes = FALSE)
  title(main = "KEGG enrichment (K=2-only)",
        adj = 0.5, font.main = 2, cex.main = 1.55, line = 1.7)
  mtext("C", side = 3, line = 1.7, at = par("usr")[1],
        adj = 0, font = 2, cex = 1.55, col = "grey20")
}

# ---- Layout + render ----
pdf(FIG_PATH, width = 13.6, height = 9.6)
lay <- rbind(
  c(1, 2, 3, 4),
  c(5, 5, 5, 5),
  c(6, 6, 7, 7)
)
layout(lay, heights = c(1.00, 0.12, 1.35))
par(mai = c(0.78, 0.84, 0.86, 0.14),
    mgp = c(2.6, 0.6, 0), oma = c(0.4, 0.4, 2.3, 0.4),
    cex.axis = 1.28, cex.lab = 1.45, font.main = 2)

# Row 1: gene exemplars (panels 1-4)
panel_letters <- c("A", "", "", "")
for (i in seq_along(exemplar_genes)) {
  g <- exemplar_genes[[i]]
  render_gene(g[1], g[2], panel_letters[i])
}

# Row 2: legend strip (panel 5)
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
       lwd    = c(NA, 2.8, 2.8, 1.4),
       horiz  = TRUE, bty = "o", box.col = "grey60", box.lwd = 0.8,
       bg = "white", cex = 1.28, xpd = NA,
       seg.len = 1.4, text.width = 0.135,
       x.intersp = 0.3)

# Row 3: Venn (panel 6) + KEGG (panel 7)
draw_venn_panel(length(venn$k1_only), length(venn$both),
                length(venn$k2_only), venn$n_total)
draw_kegg_panel(ek, max_terms = 10)

# Outer figure title
mtext("Single-harmonic (1H) vs. Two-harmonic (2H) cosinor fits on GTEx Liver exemplar genes",
      outer = TRUE, side = 3, line = 0.12, adj = 0.5,
      font = 2, cex = 1.60)

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
