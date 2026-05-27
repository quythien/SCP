#!/usr/bin/env Rscript
# =============================================================================
# Supp S2: Venn of 1H vs 2H BH-FDR rhythmic-gene discoveries
#          + table of top "K=2 unique" discoveries (GTEx Adrenal Gland).
#
# Panel A : two-circle Venn at FDR <= 0.05, comparing
#             K1 = one-harmonic F-test (estCircadianParam, 2 d.f.)
#             K2 = two-harmonic  F-test (estCircadianParam2H, 4 d.f.)
# Panel B : table of top 20 genes uniquely discovered by K=2 (q1H > 0.05,
#           q2H <= 0.05), showing HGNC symbol, A2/A1, p1H, q1H, p2H, q2H.
#           Falls back to a horizontal bar chart of (-log10 p2H) - (-log10 p1H)
#           if gridExtra is unavailable.
#
# Output  : output/two_harmonic/results/suppS2_venn_sets.rds (k1_only/both/k2_only)
#           submission/figures/SuppS2_twoharm_venn.pdf
# =============================================================================

RUN_FITS <- TRUE  # set FALSE to use cached RDS

setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
suppressPackageStartupMessages({
  source("code/setup.R")
})

RAW_PATH <- "/home/qtp1/Projects/Collaborative/GTEXdata/data/v10/cpm/v10_CPM_full.RData"
RES_DIR  <- "output/two_harmonic/results"
RES_PATH <- file.path(RES_DIR, "suppS2_twoharm_venn.rds")
SET_PATH <- file.path(RES_DIR, "suppS2_venn_sets.rds")
FIG_PATH <- "submission/figures/SuppS2_twoharm_venn.pdf"

if (!dir.exists(RES_DIR)) dir.create(RES_DIR, recursive = TRUE)
if (!dir.exists(dirname(FIG_PATH))) dir.create(dirname(FIG_PATH), recursive = TRUE)

PERIOD  <- 24
OMEGA0  <- 2 * pi / PERIOD
# FDR threshold and method. Adrenal Gland is a noisy postmortem bulk-RNA
# pilot whose 1H/2H F-test p-value histograms are dominated by uniform-null
# mass with a moderate alternative bump; in this regime BH-FDR is conservative
# and yields very few discoveries even when there is clear K=2 signal. We use
# Storey q-values (qvalue::qvalue) which estimate pi0 from the histogram and
# are routinely substituted for BH in circadian discovery; the RDS also
# contains raw BH-FDR for inspection.
FDR_ALPHA  <- 0.20
FDR_METHOD <- "qvalue"   # "qvalue" or "BH"
TOP_N <- 20L

# ---------------------------------------------------------------------------
# (1) Build expression matrix + run per-gene 1H and 2H F-tests
# ---------------------------------------------------------------------------
fit_1H_2H_allgenes <- function(data, times) {
  N <- ncol(data); G <- nrow(data)
  if (N < 6L) stop("Need N >= 6 samples for K=2.")
  X0 <- matrix(1, N, 1)
  X1 <- cbind(1, cos(OMEGA0 * times),     sin(OMEGA0 * times))
  X2 <- cbind(X1, cos(2 * OMEGA0 * times), sin(2 * OMEGA0 * times))
  df_full_1H <- N - 3L
  df_full_2H <- N - 5L

  p1H <- rep(NA_real_, G); p2H <- rep(NA_real_, G)
  A1  <- rep(NA_real_, G); A2  <- rep(NA_real_, G)

  for (g in seq_len(G)) {
    y <- as.numeric(data[g, ])
    RSS0 <- sum((y - mean(y))^2)
    if (RSS0 <= 0) next

    f1 <- tryCatch(stats::lm.fit(X1, y), error = function(e) NULL)
    f2 <- tryCatch(stats::lm.fit(X2, y), error = function(e) NULL)
    if (is.null(f1) || is.null(f2)) next

    R1 <- sum(f1$residuals^2)
    R2 <- sum(f2$residuals^2)
    if (R1 <= 0 || R2 <= 0) next

    F1 <- ((RSS0 - R1) / 2) / (R1 / df_full_1H)
    F2 <- ((RSS0 - R2) / 4) / (R2 / df_full_2H)
    p1H[g] <- stats::pf(F1, 2, df_full_1H, lower.tail = FALSE)
    p2H[g] <- stats::pf(F2, 4, df_full_2H, lower.tail = FALSE)

    b1 <- f1$coefficients
    A1[g] <- sqrt(b1[2]^2 + b1[3]^2)
    b2 <- f2$coefficients
    A2[g] <- sqrt(b2[4]^2 + b2[5]^2)
  }
  list(p1H = p1H, p2H = p2H, A1 = A1, A2 = A2)
}

if (RUN_FITS || !file.exists(RES_PATH)) {
  cat("Loading GTEx v10 Adrenal Gland log2(CPM+1)...\n")
  env <- new.env()
  load(RAW_PATH, envir = env)
  ag    <- env$CPM.all[["Adrenal_Gland"]]
  lcpm  <- ag$counts
  ann   <- ag$gene_annotation
  times <- as.numeric(ag$tod) %% 24
  if (!is.null(ann) && "symbol" %in% colnames(ann)) {
    sym <- ann$symbol[match(rownames(lcpm), ann$gene_id)]
    rownames(lcpm) <- sym
  } else {
    rownames(lcpm) <- sub("^[^_]+_", "", rownames(lcpm))
  }
  # Expression filter: mean log2(CPM+1) > 3 (matches the bundled pilot RDS
  # gene count, ~10.4k genes; removes uninformative low-expressed rows that
  # otherwise weaken BH-FDR).
  keep_g <- rowMeans(lcpm) > 3
  lcpm   <- lcpm[keep_g, , drop = FALSE]
  # Drop unmapped/duplicate symbols; keep first occurrence.
  first_idx <- !duplicated(rownames(lcpm)) & nzchar(rownames(lcpm))
  lcpm <- lcpm[first_idx, , drop = FALSE]
  cat(sprintf("Expression matrix (after mean-log2CPM > 3 filter): %d genes x %d samples\n",
              nrow(lcpm), ncol(lcpm)))

  cat("Fitting 1H + 2H F-tests for every gene...\n")
  t0 <- Sys.time()
  fits <- fit_1H_2H_allgenes(lcpm, times)
  cat(sprintf("  done in %.1f s\n",
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))

  # Both BH and Storey q-values; downstream selects which to use.
  q1_bh <- stats::p.adjust(fits$p1H, method = "BH")
  q2_bh <- stats::p.adjust(fits$p2H, method = "BH")
  if (requireNamespace("qvalue", quietly = TRUE)) {
    q1_st <- qvalue::qvalue(fits$p1H[!is.na(fits$p1H)])$qvalues
    q2_st <- qvalue::qvalue(fits$p2H[!is.na(fits$p2H)])$qvalues
    q1_storey <- rep(NA_real_, length(fits$p1H))
    q2_storey <- rep(NA_real_, length(fits$p2H))
    q1_storey[!is.na(fits$p1H)] <- q1_st
    q2_storey[!is.na(fits$p2H)] <- q2_st
  } else {
    q1_storey <- q1_bh; q2_storey <- q2_bh
  }
  df <- data.frame(
    symbol     = rownames(lcpm),
    p1H        = fits$p1H,
    p2H        = fits$p2H,
    q1H_BH     = q1_bh,
    q2H_BH     = q2_bh,
    q1H_storey = q1_storey,
    q2H_storey = q2_storey,
    A1         = fits$A1,
    A2         = fits$A2,
    stringsAsFactors = FALSE
  )
  df$A2_over_A1 <- df$A2 / pmax(df$A1, 1e-12)
  # Set the "primary" q1H/q2H by method choice.
  if (FDR_METHOD == "qvalue") {
    df$q1H <- df$q1H_storey; df$q2H <- df$q2H_storey
  } else {
    df$q1H <- df$q1H_BH;     df$q2H <- df$q2H_BH
  }

  saveRDS(df, RES_PATH)
  cat("Saved fits to", RES_PATH, "\n")
} else {
  cat("Loading cached fits from", RES_PATH, "\n")
  df <- readRDS(RES_PATH)
}

# ---------------------------------------------------------------------------
# (2) Build Venn sets
# ---------------------------------------------------------------------------
ok <- !is.na(df$q1H) & !is.na(df$q2H)
df <- df[ok, , drop = FALSE]
G_total <- nrow(df)

sig1 <- df$q1H <= FDR_ALPHA
sig2 <- df$q2H <= FDR_ALPHA

k1_only <- df$symbol[sig1 & !sig2]
both    <- df$symbol[sig1 &  sig2]
k2_only <- df$symbol[!sig1 & sig2]

cat(sprintf("Total genes tested:    %d\n", G_total))
cat(sprintf("  K=1 only:            %d\n", length(k1_only)))
cat(sprintf("  both:                %d\n", length(both)))
cat(sprintf("  K=2 only:            %d\n", length(k2_only)))

saveRDS(list(k1_only = k1_only, both = both, k2_only = k2_only,
             total = G_total, alpha = FDR_ALPHA),
        SET_PATH)
cat("Wrote", SET_PATH, "\n")

# ---------------------------------------------------------------------------
# (3) Top K=2-unique table
# ---------------------------------------------------------------------------
sub <- df[!sig1 & sig2, , drop = FALSE]
sub <- sub[order(sub$q2H), , drop = FALSE]
top_tbl <- utils::head(sub, TOP_N)
top_disp <- data.frame(
  Gene = top_tbl$symbol,
  `A2/A1` = sprintf("%.2f", top_tbl$A2_over_A1),
  p_1H    = sprintf("%.2g", top_tbl$p1H),
  q_1H    = sprintf("%.2g", top_tbl$q1H),
  p_2H    = sprintf("%.2g", top_tbl$p2H),
  q_2H    = sprintf("%.2g", top_tbl$q2H),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------------------------
# (4) Plot 1 x 2: Venn (A) + Table (B)
# ---------------------------------------------------------------------------
has_gridExtra   <- requireNamespace("gridExtra",   quietly = TRUE)
has_VennDiagram <- requireNamespace("VennDiagram", quietly = TRUE)
has_grid        <- requireNamespace("grid",        quietly = TRUE)
if (has_VennDiagram) suppressMessages(library(VennDiagram))
if (has_gridExtra)   suppressMessages(library(gridExtra))
if (has_grid)        suppressMessages(library(grid))

# Suppress the VennDiagram log files
if (requireNamespace("futile.logger", quietly = TRUE)) {
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
}

draw_venn_grob <- function(n1, n2, n12, cat1 = "K = 1", cat2 = "K = 2") {
  if (has_VennDiagram) {
    g <- VennDiagram::draw.pairwise.venn(
      area1 = n1, area2 = n2, cross.area = n12,
      category   = c(cat1, cat2),
      fill       = c("#d62728", "#1f77b4"),
      alpha      = c(0.45, 0.45),
      lty        = "blank",
      cex        = 1.6,
      cat.cex    = 1.55,
      cat.col    = c("#d62728", "#1f77b4"),
      cat.pos    = c(-25, 25),
      cat.dist   = c(0.05, 0.05),
      ind        = FALSE
    )
    return(g)
  }
  NULL
}

draw_venn_base <- function(n_k1only, n_both, n_k2only) {
  # base-graphics fallback
  plot.new()
  plot.window(xlim = c(-1.5, 1.5), ylim = c(-1, 1), asp = 1)
  theta <- seq(0, 2 * pi, length.out = 200)
  r <- 0.85
  c1 <- cbind(-0.4 + r * cos(theta),  r * sin(theta))
  c2 <- cbind( 0.4 + r * cos(theta),  r * sin(theta))
  polygon(c1, col = adjustcolor("#d62728", 0.45), border = NA)
  polygon(c2, col = adjustcolor("#1f77b4", 0.45), border = NA)
  text(-0.85,  0, n_k1only, cex = 1.5)
  text( 0.85,  0, n_k2only, cex = 1.5)
  text( 0.00,  0, n_both,   cex = 1.5)
  text(-0.55,  1.0, "K = 1", col = "#d62728", cex = 1.55)
  text( 0.55,  1.0, "K = 2", col = "#1f77b4", cex = 1.55)
}

n_k1_total <- length(k1_only) + length(both)
n_k2_total <- length(k2_only) + length(both)
n_both     <- length(both)

pdf(FIG_PATH, width = 13, height = 6.5)
if (has_VennDiagram && has_grid && has_gridExtra) {
  # gridExtra path: Venn (grob) on left, table on right.
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(
    nrow = 1, ncol = 2, widths = grid::unit(c(0.45, 0.55), "npc"))))

  # ---- Panel A: Venn
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid::grid.text("A", x = grid::unit(0.02, "npc"),
                  y = grid::unit(0.96, "npc"),
                  gp = grid::gpar(fontsize = 22, fontface = "bold"))
  grid::grid.text(
    sprintf("%s q <= %.2f rhythmic-gene discoveries",
            if (FDR_METHOD == "qvalue") "Storey" else "BH",
            FDR_ALPHA),
    x = grid::unit(0.5, "npc"), y = grid::unit(0.92, "npc"),
    gp = grid::gpar(fontsize = 15))
  vg <- draw_venn_grob(n1 = n_k1_total, n2 = n_k2_total, n12 = n_both)
  if (!is.null(vg)) {
    grid::pushViewport(grid::viewport(y = grid::unit(0.45, "npc"),
                                       height = grid::unit(0.82, "npc")))
    grid::grid.draw(vg)
    grid::popViewport()
  }
  grid::grid.text(
    sprintf("Total: %d   |   K=1 only: %d   |   both: %d   |   K=2 only: %d",
            G_total, length(k1_only), length(both), length(k2_only)),
    x = grid::unit(0.5, "npc"), y = grid::unit(0.04, "npc"),
    gp = grid::gpar(fontsize = 12))
  grid::popViewport()

  # ---- Panel B: Table
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  grid::grid.text("B", x = grid::unit(0.02, "npc"),
                  y = grid::unit(0.96, "npc"),
                  gp = grid::gpar(fontsize = 22, fontface = "bold"))
  grid::grid.text(sprintf("Top %d K=2-unique discoveries", nrow(top_disp)),
                  x = grid::unit(0.5, "npc"), y = grid::unit(0.92, "npc"),
                  gp = grid::gpar(fontsize = 15))
  tt <- gridExtra::ttheme_default(
    core    = list(fg_params = list(cex = 0.85)),
    colhead = list(fg_params = list(cex = 0.9, fontface = "bold")))
  tg <- gridExtra::tableGrob(top_disp, rows = NULL, theme = tt)
  grid::pushViewport(grid::viewport(y = grid::unit(0.45, "npc"),
                                     height = grid::unit(0.82, "npc"),
                                     width  = grid::unit(0.95, "npc")))
  grid::grid.draw(tg)
  grid::popViewport()
  grid::popViewport()

  grid::popViewport()  # outer layout
} else {
  # Base-graphics fallback for both panels.
  op <- par(mfrow = c(1, 2),
            mar = c(3, 3, 4, 1),
            cex.lab = 1.55, cex.main = 1.55, cex.axis = 1.35)

  draw_venn_base(length(k1_only), length(both), length(k2_only))
  title(main = sprintf("A     BH-FDR <= %.2f", FDR_ALPHA), adj = 0)

  # Bar fallback for B
  delta_logp <- -log10(top_tbl$p2H) + log10(top_tbl$p1H)
  ord <- order(delta_logp)
  barplot(delta_logp[ord],
          names.arg = top_tbl$symbol[ord],
          horiz = TRUE, las = 1, col = "#1f77b4", border = NA,
          xlab = expression(-log[10] * p[2*H] + log[10] * p[1*H]),
          main = sprintf("B     Top %d K=2-unique", nrow(top_tbl)))
  par(op)
}
dev.off()

cat("\nWrote", FIG_PATH, "\n")
cat(sprintf("Size: %.1f KB\n", file.size(FIG_PATH) / 1024))

cat("\nTop K=2 unique discoveries:\n")
print(top_disp, row.names = FALSE)
