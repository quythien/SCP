#!/usr/bin/env Rscript
# =============================================================================
# Figure 4, one-column (86 mm) pieces for the Bioinformatics port.
#
# Panels B (K=1 / K=2 overlap) and C (KEGG enrichment of the K=2-only genes)
# are redrawn here from the cached results of two_harmonic/fig4_twoharm_demo.R
# (tissue_Liver_venn.rds, tissue_Liver_kegg.rds) at column width with 7-8 pt
# text -> figures/onecol/Fig4BC_1col.pdf.
#
# Panel A (four exemplar genes with raw expression points) needs the GTEx v10
# Liver CPM matrix (controlled access, server only). The manuscript therefore
# clips the four gene panels out of the full-width Fig4_twoharm_demo.pdf with
# graphicx viewports (see port_bioinformatics.py, step 6b) and places this
# B/C strip under them.
# Run from the package root: Rscript examples/publication/onecol/fig4_onecol.R
# =============================================================================
setwd(Sys.getenv("SCP_ROOT", unset = "."))
source("examples/publication/onecol/onecol_common.R")

venn <- readRDS("output/two_harmonic/results/tissue_Liver_venn.rds")
ek   <- readRDS("output/two_harmonic/results/tissue_Liver_kegg.rds")

draw_venn <- function(n_k1only, n_both, n_k2only, n_total) {
  par(mai = c(0.02, 0.02, 0.22, 0.02))
  plot.new(); plot.window(xlim = c(-1.75, 1.75), ylim = c(-1.35, 1.25), asp = 1)
  theta <- seq(0, 2 * pi, length.out = 200); r <- 0.95
  polygon(cbind(-0.55 + r * cos(theta), r * sin(theta)),
          col = adjustcolor("#d62728", 0.45), border = "#d62728", lwd = 1)
  polygon(cbind(0.55 + r * cos(theta), r * sin(theta)),
          col = adjustcolor("#1f77b4", 0.45), border = "#1f77b4", lwd = 1)
  text(-1.10, -0.05, format(n_k1only, big.mark = ","), cex = 1.0, font = 2)
  text( 1.10, -0.05, format(n_k2only, big.mark = ","), cex = 1.0, font = 2)
  text( 0.00, -0.05, format(n_both,   big.mark = ","), cex = 1.0, font = 2)
  text(-0.85,  1.08, "K=1", col = "#d62728", cex = 0.9, font = 2)
  text( 0.85,  1.08, "K=2", col = "#1f77b4", cex = 0.9, font = 2)
  text(0, -1.22, sprintf("n = %s genes", format(n_total, big.mark = ",")),
       cex = CEX_LEG, col = "grey25")
  panel_title("B", "Gene overlap", line = 0.1)
}

draw_kegg <- function(ek, max_terms = 10, wrap = 27) {
  df_e <- as.data.frame(ek)
  if ("category" %in% colnames(df_e))
    df_e <- df_e[df_e$category != "Human Diseases", , drop = FALSE]
  df_e$.p <- if ("pvalue" %in% colnames(df_e)) df_e$pvalue else df_e$p.adjust
  df_e <- head(df_e[order(df_e$.p), ], max_terms)
  df_e <- df_e[rev(seq_len(nrow(df_e))), ]          # most significant on top
  val  <- -log10(df_e$.p)
  labs <- vapply(df_e$Description, function(s)
    paste(strwrap(s, width = wrap), collapse = "\n"), character(1))
  pal  <- colorRampPalette(c("#fee0b6", "#f1a340", "#b35806"))(nrow(df_e))
  par(mai = c(0.34, 1.22, 0.22, 0.06), mgp = c(1.5, 0.4, 0))
  bp <- barplot(val, horiz = TRUE, names.arg = rep("", length(val)), col = pal,
                border = NA, xlab = expression(-log[10]("P-value")),
                cex.lab = CEX_LAB, cex.axis = CEX_AXIS, xlim = c(0, max(val) * 1.08),
                main = "", space = 0.45)
  grid(nx = NULL, ny = NA, lty = "dotted", col = "grey80")
  barplot(val, horiz = TRUE, col = pal, border = NA, add = TRUE, axes = FALSE,
          space = 0.45)
  # pathway names drawn with text(): axis() silently drops labels it judges
  # overlapping, which lost two of the five names in a 1.6 in panel
  text(-0.02 * max(val), bp, labels = labs, adj = c(1, 0.5), xpd = NA, cex = CEX_AXIS)
  title(main = "C   KEGG, K=2-only genes", adj = 1, font.main = 2, cex.main = CEX_MAIN, line = 0.1)
}

OUT <- file.path(OUT_DIR, "Fig4BC_1col.pdf")
pdf(OUT, width = COL_W, height = 1.62, pointsize = PT)
layout(matrix(1:2, nrow = 1), widths = c(0.36, 0.64))
par(cex = 1, font.main = 2, tcl = -0.25)
draw_venn(length(venn$k1_only), length(venn$both), length(venn$k2_only), venn$n_total)
draw_kegg(ek)
dev.off()
cat(sprintf("Saved %s (K1-only=%d, both=%d, K2-only=%d of %d genes)\n", OUT,
            length(venn$k1_only), length(venn$both), length(venn$k2_only), venn$n_total))
