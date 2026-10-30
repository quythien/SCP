# 1. Read the saved overlap and enrichment results.
FIGURE_DIR <- "Figure4"
source("code/onecol_common.R")

venn <- readRDS(file.path(FIGURE_DIR, "data", "tissue_Liver_venn.rds"))
ek <- readRDS(file.path(FIGURE_DIR, "data", "kegg_table.rds"))

# 2. Draw the overlap panel in the original one-column style.
draw_venn <- function(n_k1only, n_both, n_k2only, n_total) {
  par(mai = c(0.02, 0.02, 0.22, 0.02))
  plot.new()
  plot.window(xlim = c(-1.75, 1.75), ylim = c(-1.35, 1.25), asp = 1)
  theta <- seq(0, 2 * pi, length.out = 200)
  r <- 0.95
  polygon(cbind(-0.55 + r * cos(theta), r * sin(theta)),
          col = adjustcolor("#d62728", 0.45), border = "#d62728", lwd = 1)
  polygon(cbind(0.55 + r * cos(theta), r * sin(theta)),
          col = adjustcolor("#1f77b4", 0.45), border = "#1f77b4", lwd = 1)
  text(-1.10, -0.05, format(n_k1only, big.mark = ","), cex = 1.0, font = 2)
  text(1.10, -0.05, format(n_k2only, big.mark = ","), cex = 1.0, font = 2)
  text(0.00, -0.05, format(n_both, big.mark = ","), cex = 1.0, font = 2)
  text(-0.85, 1.08, "K=1", col = "#d62728", cex = 0.9, font = 2)
  text(0.85, 1.08, "K=2", col = "#1f77b4", cex = 0.9, font = 2)
  text(0, -1.22, sprintf("n = %s genes", format(n_total, big.mark = ",")),
       cex = CEX_LEG, col = "grey25")
  panel_title("B", "Gene overlap", line = 0.1)
}

# 3. Draw the KEGG panel from the saved enrichment table.
draw_kegg <- function(ek, max_terms = 10, wrap = 27) {
  df_e <- as.data.frame(ek)
  if ("category" %in% colnames(df_e))
    df_e <- df_e[df_e$category != "Human Diseases", , drop = FALSE]
  df_e$.p <- if ("pvalue" %in% colnames(df_e)) df_e$pvalue else df_e$p.adjust
  df_e <- head(df_e[order(df_e$.p), ], max_terms)
  df_e <- df_e[rev(seq_len(nrow(df_e))), ]
  val <- -log10(df_e$.p)
  labs <- vapply(df_e$Description, function(s)
    paste(strwrap(s, width = wrap), collapse = "\n"), character(1))
  pal <- colorRampPalette(c("#fee0b6", "#f1a340", "#b35806"))(nrow(df_e))
  par(mai = c(0.34, 1.22, 0.22, 0.06), mgp = c(1.5, 0.4, 0))
  bp <- barplot(val, horiz = TRUE, names.arg = rep("", length(val)), col = pal,
                border = NA, xlab = expression(-log[10]("P-value")),
                cex.lab = CEX_LAB, cex.axis = CEX_AXIS,
                xlim = c(0, max(val) * 1.08), main = "", space = 0.45)
  grid(nx = NULL, ny = NA, lty = "dotted", col = "grey80")
  barplot(val, horiz = TRUE, col = pal, border = NA, add = TRUE, axes = FALSE,
          space = 0.45)
  text(-0.02 * max(val), bp, labels = labs, adj = c(1, 0.5),
       xpd = NA, cex = CEX_AXIS)
  title(main = "C   KEGG, K=2-only genes", adj = 1, font.main = 2,
        cex.main = CEX_MAIN, line = 0.1)
}

# 4. Save the redrawn B and C panels.
pdf(file.path(OUT_DIR, "Figure4BC.pdf"), width = COL_W, height = 1.62,
    pointsize = PT)
layout(matrix(1:2, nrow = 1), widths = c(0.36, 0.64))
par(cex = 1, font.main = 2, tcl = -0.25)
draw_venn(length(venn$k1_only), length(venn$both),
          length(venn$k2_only), venn$n_total)
draw_kegg(ek)
dev.off()

# 5. Copy the archived vector figure used as the source for panel A.
stopifnot(file.copy(file.path(FIGURE_DIR, "final", "Fig4_twoharm_demo.pdf"),
                    file.path(OUT_DIR, "Figure4A.pdf"), overwrite = TRUE))
cat("Saved", file.path(OUT_DIR, "Figure4A.pdf"), "and",
    file.path(OUT_DIR, "Figure4BC.pdf"), "\n")

# 6. Retain the assembled observed-data figure used in the manuscript.
stopifnot(file.copy("Figure4/final/Figure4.pdf",
  file.path(OUT_DIR,"Figure4.pdf"),overwrite=TRUE))
