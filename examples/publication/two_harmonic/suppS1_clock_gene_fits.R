#!/usr/bin/env Rscript
# =============================================================================
# Supplementary Fig S1: 1H vs 2H cosinor fits on 9 clock / clock-controlled
# genes from GTEx Liver.
#
# 3x3 grid; rows organised as a "K=2 advantage" gradient:
#   Row 1 -- canonical clock TFs, detected by both K=1 and K=2 (ARNTL, NPAS2,
#            CLOCK)
#   Row 2 -- secondary clock genes, both detected (PER2, NR1D1, CRY1)
#   Row 3 -- K=2-only / borderline: RORC (clock), HADHA (mitochondrial CCG),
#            TPP1 (lysosomal CCG)
#
# For each gene: scatter of raw log2(CPM+1), 1H fit (dashed red), 2H fit
# (solid blue); panel title shows symbol and (q1H, q2H).
#
# Output: submission/figures/SuppS1_clock_gene_fits.pdf
# =============================================================================

setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")

GTEX <- "/home/qtp1/Projects/Collaborative/GTEXdata/data/v10/cpm/v10_CPM_full.RData"
FIG  <- "submission/figures/SuppS1_clock_gene_fits.pdf"
PER  <- 24
OMEGA0 <- 2 * pi / PER

# 3 x 3 gene grid with one-line role label
genes <- list(
  c("ARNTL",  "Core clock TF"),
  c("NPAS2",  "Core clock TF"),
  c("CLOCK",  "Core clock TF"),
  c("PER2",   "Repressor"),
  c("NR1D1",  "REV-ERB"),
  c("CRY1",   "Repressor"),
  c("RORC",   "Clock K=2 only"),
  c("HADHA",  "Mitochondrial CCG"),
  c("TPP1",   "Lysosomal CCG")
)

# ---- Load Liver expression matrix ----
env <- new.env(); load(GTEX, envir = env)
liv  <- env$CPM.all[["Liver"]]
lcpm <- liv$counts
ann  <- liv$gene_annotation
tod  <- as.numeric(liv$tod) %% PER
ok   <- !is.na(tod)
lcpm <- lcpm[, ok]; tod <- tod[ok]
if (!is.null(ann) && "symbol" %in% colnames(ann))
  rownames(lcpm) <- ann$symbol[match(rownames(lcpm), ann$gene_id)]
first <- !duplicated(rownames(lcpm)) & nzchar(rownames(lcpm)) &
         !is.na(rownames(lcpm))
lcpm <- lcpm[first, ]
keep <- rowMeans(lcpm) > 3
lcpm <- lcpm[keep, ]
cat(sprintf("Liver expression matrix: %d genes x %d samples\n",
            nrow(lcpm), ncol(lcpm)))

# ---- Per-gene 1H/2H F-tests on entire matrix (for honest BH-FDR labels) ----
N <- ncol(lcpm)
X1 <- cbind(1, cos(OMEGA0 * tod), sin(OMEGA0 * tod))
X2 <- cbind(X1, cos(2 * OMEGA0 * tod), sin(2 * OMEGA0 * tod))
df1 <- N - 3L; df2 <- N - 5L

fit_all <- function() {
  G <- nrow(lcpm)
  p1 <- p2 <- rep(NA_real_, G)
  for (g in seq_len(G)) {
    y <- as.numeric(lcpm[g, ])
    R0 <- sum((y - mean(y))^2); if (R0 <= 0) next
    f1 <- lm.fit(X1, y); f2 <- lm.fit(X2, y)
    R1 <- sum(f1$residuals^2); R2 <- sum(f2$residuals^2)
    if (R1 <= 0 || R2 <= 0) next
    p1[g] <- pf(((R0 - R1)/2) / (R1/df1), 2, df1, lower.tail = FALSE)
    p2[g] <- pf(((R0 - R2)/4) / (R2/df2), 4, df2, lower.tail = FALSE)
  }
  list(p1 = p1, p2 = p2, q1 = p.adjust(p1, "BH"), q2 = p.adjust(p2, "BH"))
}
fits_all <- fit_all()

# ---- Helper to render one gene panel ----
render_gene <- function(gene, lab) {
  if (!(gene %in% rownames(lcpm))) {
    plot.new(); title(sprintf("%s (missing)", gene)); return(invisible())
  }
  gi <- which(rownames(lcpm) == gene)[1]
  y  <- as.numeric(lcpm[gi, ])
  f1 <- lm.fit(X1, y); f2 <- lm.fit(X2, y)
  b1 <- f1$coefficients; b2 <- f2$coefficients
  q1 <- fits_all$q1[gi]; q2 <- fits_all$q2[gi]

  flag <- if (q1 <= 0.05 && q2 <= 0.05) "both"
          else if (q1 >  0.05 && q2 <= 0.05) "K=2 only"
          else if (q1 <= 0.05 && q2 >  0.05) "K=1 only"
          else "neither"

  tt <- seq(0, 24, length.out = 200)
  fit1_curve <- b1[1] + b1[2]*cos(OMEGA0*tt) + b1[3]*sin(OMEGA0*tt)
  fit2_curve <- b2[1] + b2[2]*cos(OMEGA0*tt) + b2[3]*sin(OMEGA0*tt) +
                       b2[4]*cos(2*OMEGA0*tt) + b2[5]*sin(2*OMEGA0*tt)
  ylim <- range(c(y, fit1_curve, fit2_curve), na.rm = TRUE)

  plot(tod, y, pch = 19, col = adjustcolor("grey50", 0.5), cex = 0.55,
       xlab = "Time of day (h)", ylab = expression(log[2](CPM+1)),
       xlim = c(0, 24), ylim = ylim, main = "")
  lines(tt, fit1_curve, col = "red",  lwd = 2.2, lty = 2)
  lines(tt, fit2_curve, col = "blue", lwd = 2.2, lty = 1)
  axis(1, at = seq(0, 24, 6))
  box()
  title(main = sprintf("%s  (%s)", gene, flag),
        adj = 0, font.main = 2, cex.main = 1.2, line = 0.7)
  mtext(sprintf("q1H=%.2g  q2H=%.2g  [%s]", q1, q2, lab),
        side = 3, line = -0.1, adj = 0.5, cex = 0.75, col = "grey25")
}

# ---- Render figure ----
pdf(FIG, width = 11.5, height = 10.5)
par(mfrow = c(3, 3), mai = c(0.78, 0.85, 0.70, 0.18),
    mgp = c(2.6, 0.7, 0), oma = c(1.4, 0, 2.5, 0),
    cex.axis = 1.15, cex.lab = 1.25, font.main = 2)
for (g in genes) render_gene(g[1], g[2])
mtext("Clock-gene 1H vs 2H cosinor fits (GTEx Liver)",
      outer = TRUE, cex = 1.45, font = 2, line = 0.5)

# Footer legend in outer margin
par(new = TRUE, mfrow = c(1,1), mar = c(0,0,0,0), oma = c(0,0,0,0))
plot.new(); plot.window(0:1, 0:1)
legend("bottom",
       legend = c("Raw data", "1H (DCP) fit", "2H fit"),
       col    = c(adjustcolor("grey50", 0.6), "red", "blue"),
       pch    = c(19, NA, NA),
       lty    = c(NA, 2, 1),
       lwd    = c(NA, 2.2, 2.2),
       horiz  = TRUE, bty = "n", cex = 1.05, xpd = NA, inset = c(0, 0.0))
dev.off()

cat(sprintf("\nWrote %s (%.1f KB)\n", FIG, file.size(FIG)/1024))

cat("\nPer-gene BH-adjusted q-values:\n")
for (g in genes) {
  if (!(g[1] %in% rownames(lcpm))) {
    cat(sprintf("  %-8s  MISSING\n", g[1])); next
  }
  gi <- which(rownames(lcpm) == g[1])[1]
  cat(sprintf("  %-8s  q1H=%.3g   q2H=%.3g\n",
              g[1], fits_all$q1[gi], fits_all$q2[gi]))
}
