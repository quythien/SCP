#!/usr/bin/env Rscript
# =============================================================================
# tissue_venn_enrichment.R
#
# End-to-end K=1 (DCP) vs K=2 cosinor rhythmicity comparison on a single GTEx
# tissue, with Venn-set extraction, clock-gene flagging, and pathway enrichment
# (KEGG + GO:BP).
#
# Inputs:
#   - GTEx v10 CPM matrix at GTEX_CPM_PATH (per-tissue list with $counts and $tod)
#   - Tissue name (default: "Liver")
#
# Outputs (under output/two_harmonic/results/):
#   - tissue_<TISSUE>_fits.rds      Full per-gene fits (p1H, q1H, p2H, q2H,
#                                    A1, A2, phi1, phi2, A2/A1)
#   - tissue_<TISSUE>_venn.rds      List with k1_only / both / k2_only symbols
#   - tissue_<TISSUE>_clock.rds     Subset on a curated clock + Hughes-2009
#                                    ultradian gene set with detection status
#   - tissue_<TISSUE>_kegg.rds      enrichKEGG object (clusterProfiler)
#   - tissue_<TISSUE>_gobp.rds      enrichGO  object (clusterProfiler)
#
# Console output: counts, top-N tables, clock-gene table.
#
# Reproducibility: deterministic given the input matrix. No simulation.
#
# Usage:
#   Rscript examples/publication/two_harmonic/tissue_venn_enrichment.R [TISSUE]
#
# Example:
#   Rscript examples/publication/two_harmonic/tissue_venn_enrichment.R Liver
#   Rscript examples/publication/two_harmonic/tissue_venn_enrichment.R Adrenal_Gland
# =============================================================================

# ---- 1. Configuration --------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
TISSUE <- if (length(args) >= 1) args[1] else "Liver"

setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")

GTEX_CPM_PATH <- "/home/qtp1/Projects/Collaborative/GTEXdata/data/v10/cpm/v10_CPM_full.RData"
RES_DIR       <- "output/two_harmonic/results"
PERIOD        <- 24
OMEGA0        <- 2 * pi / PERIOD
MEAN_EXPR_CUT <- 3      # log2(CPM+1) >= 3 to filter low-expressed rows
FDR_ALPHA     <- 0.05   # BH-FDR threshold for Venn

if (!dir.exists(RES_DIR)) dir.create(RES_DIR, recursive = TRUE)

# Curated clock + Hughes 2009 12-h ultradian pathway markers.
CLOCK_AND_CCG <- c(
  # Core clock TFs
  "ARNTL", "BMAL1", "CLOCK", "NPAS2",
  "PER1", "PER2", "PER3", "CRY1", "CRY2",
  "NR1D1", "NR1D2", "RORA", "RORB", "RORC",
  "DBP", "TEF", "HLF", "NFIL3", "BHLHE40", "BHLHE41", "CIART",
  "TIMELESS", "CSNK1D", "CSNK1E", "FBXL3",
  # Hughes 2009 ultradian-enriched pathways (mitochondrial, lysosomal,
  # translation, proteostasis), used as positive-control markers.
  "ACO2", "NDUFS2", "NDUFA8", "ATP5F1B", "SDHA", "HADHA", "HADHB",
  "TPP1", "GBA", "CTSB", "CTSD",
  "EIF3A", "EIF4G1", "PSMD1", "PSMA4",
  # Liver-specific clock-controlled
  "HNF4A", "PPARA", "SREBF1", "INSIG1", "INSIG2",
  "G6PC1", "PCK1", "CYP7A1", "CYP3A4", "CYP2E1",
  "ALAS1", "ELOVL2", "ELOVL3"
)

cat("=========================================================\n")
cat(sprintf("Tissue-level K=1 vs K=2 Venn + enrichment: %s\n", TISSUE))
cat("=========================================================\n")

# ---- 2. Load expression matrix -----------------------------------------------
cat("\n[1/6] Loading GTEx v10 ...\n")
env <- new.env()
load(GTEX_CPM_PATH, envir = env)
if (!TISSUE %in% names(env$CPM.all))
  stop(sprintf("Tissue '%s' not found. Available: %s",
               TISSUE, paste(names(env$CPM.all), collapse = ", ")))

td   <- env$CPM.all[[TISSUE]]
lcpm <- td$counts
ann  <- td$gene_annotation
tod  <- as.numeric(td$tod) %% 24

ok   <- !is.na(tod)
lcpm <- lcpm[, ok, drop = FALSE]
tod  <- tod[ok]

# Filter low-expressed rows.
keep <- rowMeans(lcpm) > MEAN_EXPR_CUT
lcpm <- lcpm[keep, , drop = FALSE]

# Map Ensembl row names to HGNC symbols where possible.
if (!is.null(ann) && "symbol" %in% colnames(ann)) {
  rownames(lcpm) <- ann$symbol[match(rownames(lcpm), ann$gene_id)]
} else {
  rownames(lcpm) <- sub("^[^_]+_", "", rownames(lcpm))
}
first_idx <- !duplicated(rownames(lcpm)) & nzchar(rownames(lcpm)) & !is.na(rownames(lcpm))
lcpm <- lcpm[first_idx, , drop = FALSE]

cat(sprintf("  Filtered expression matrix: %d genes x %d samples (mean log2 CPM > %g)\n",
            nrow(lcpm), ncol(lcpm), MEAN_EXPR_CUT))

# ---- 3. K=1 (DCP) and K=2 F-tests per gene -----------------------------------
cat("\n[2/6] Per-gene 1H (DCP) and 2H F-tests ...\n")
fit_1H_2H <- function(data, times) {
  N <- ncol(data); G <- nrow(data)
  if (N < 6L) stop("Need N >= 6 samples for K=2.")
  X1 <- cbind(1, cos(OMEGA0 * times), sin(OMEGA0 * times))
  X2 <- cbind(X1, cos(2 * OMEGA0 * times), sin(2 * OMEGA0 * times))
  df1 <- N - 3L; df2 <- N - 5L
  p1H <- p2H <- A1 <- A2 <- phi1 <- phi2 <- rep(NA_real_, G)
  for (g in seq_len(G)) {
    y <- as.numeric(data[g, ])
    RSS0 <- sum((y - mean(y))^2); if (RSS0 <= 0) next
    f1 <- tryCatch(stats::lm.fit(X1, y), error = function(e) NULL)
    f2 <- tryCatch(stats::lm.fit(X2, y), error = function(e) NULL)
    if (is.null(f1) || is.null(f2)) next
    R1 <- sum(f1$residuals^2); R2 <- sum(f2$residuals^2)
    if (R1 <= 0 || R2 <= 0) next
    F1 <- ((RSS0 - R1) / 2) / (R1 / df1)
    F2 <- ((RSS0 - R2) / 4) / (R2 / df2)
    p1H[g] <- stats::pf(F1, 2, df1, lower.tail = FALSE)
    p2H[g] <- stats::pf(F2, 4, df2, lower.tail = FALSE)
    b1 <- f1$coefficients
    b2 <- f2$coefficients
    A1[g]   <- sqrt(b1[2]^2 + b1[3]^2)
    phi1[g] <- (atan2(b1[3], b1[2]) / OMEGA0) %% PERIOD
    A2[g]   <- sqrt(b2[4]^2 + b2[5]^2)
    phi2[g] <- (atan2(b2[5], b2[4]) / (2 * OMEGA0)) %% (PERIOD / 2)
  }
  data.frame(
    gene = rownames(data),
    p1H = p1H, p2H = p2H,
    A1 = A1, A2 = A2,
    phi1 = phi1, phi2 = phi2,
    A2_over_A1 = A2 / pmax(A1, 1e-9),
    stringsAsFactors = FALSE
  )
}

t0 <- Sys.time()
df <- fit_1H_2H(lcpm, tod)
cat(sprintf("  done in %.1f s\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# BH-FDR on each F-test, applied across all tested genes.
df$q1H <- stats::p.adjust(df$p1H, method = "BH")
df$q2H <- stats::p.adjust(df$p2H, method = "BH")

# Discard any rows with NA p-values from downstream sets.
df <- df[!is.na(df$q1H) & !is.na(df$q2H), ]

# ---- 4. Venn sets at BH-FDR <= FDR_ALPHA -------------------------------------
cat(sprintf("\n[3/6] Building Venn sets at BH-FDR <= %g ...\n", FDR_ALPHA))
sig1 <- df$q1H <= FDR_ALPHA
sig2 <- df$q2H <= FDR_ALPHA

k1_only <- df$gene[sig1 & !sig2]
both    <- df$gene[sig1 &  sig2]
k2_only <- df$gene[!sig1 & sig2]

cat(sprintf("  Total genes tested:  %d\n", nrow(df)))
cat(sprintf("  K=1 only:            %d\n", length(k1_only)))
cat(sprintf("  Both K=1 and K=2:    %d\n", length(both)))
cat(sprintf("  K=2 only:            %d\n", length(k2_only)))

venn_obj <- list(
  tissue     = TISSUE,
  fdr_alpha  = FDR_ALPHA,
  n_total    = nrow(df),
  k1_only    = k1_only,
  both       = both,
  k2_only    = k2_only
)
saveRDS(df, file.path(RES_DIR, sprintf("tissue_%s_fits.rds", TISSUE)))
saveRDS(venn_obj, file.path(RES_DIR, sprintf("tissue_%s_venn.rds", TISSUE)))

# ---- 5. Clock-gene status table ----------------------------------------------
cat("\n[4/6] Clock + CCG gene detection status ...\n")
clock_present <- intersect(CLOCK_AND_CCG, df$gene)
clock_tbl <- df[df$gene %in% clock_present, ]
clock_tbl$status <- with(clock_tbl,
  ifelse(q1H <= FDR_ALPHA & q2H <= FDR_ALPHA, "both",
  ifelse(q1H >  FDR_ALPHA & q2H <= FDR_ALPHA, "K2_only",
  ifelse(q1H <= FDR_ALPHA & q2H >  FDR_ALPHA, "K1_only", "neither"))))
clock_tbl <- clock_tbl[order(clock_tbl$status, clock_tbl$q2H),
                       c("gene", "status", "q1H", "q2H", "A1", "A2",
                         "A2_over_A1", "phi1", "phi2")]
print(clock_tbl, row.names = FALSE, digits = 3)
saveRDS(clock_tbl, file.path(RES_DIR, sprintf("tissue_%s_clock.rds", TISSUE)))

# ---- 6. Pathway enrichment ---------------------------------------------------
cat("\n[5/6] KEGG + GO:BP pathway enrichment for K=2-only ...\n")
have_cp  <- requireNamespace("clusterProfiler", quietly = TRUE)
have_org <- requireNamespace("org.Hs.eg.db",    quietly = TRUE)

if (!have_cp || !have_org) {
  cat("  clusterProfiler or org.Hs.eg.db unavailable; enrichment skipped.\n")
  cat("  Install with BiocManager::install(c('clusterProfiler','org.Hs.eg.db')).\n")
} else {
  suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Hs.eg.db)
  })
  sym2eid <- function(s) {
    m <- AnnotationDbi::select(org.Hs.eg.db, keys = s,
                                columns = "ENTREZID", keytype = "SYMBOL")
    m <- m[!is.na(m$ENTREZID) & !duplicated(m$SYMBOL), ]
    m$ENTREZID[match(s, m$SYMBOL)]
  }
  bg_eid <- sym2eid(unique(df$gene));    bg_eid <- bg_eid[!is.na(bg_eid)]
  k2_eid <- sym2eid(unique(k2_only));    k2_eid <- k2_eid[!is.na(k2_eid)]

  cat(sprintf("  Background Entrez genes:    %d\n", length(bg_eid)))
  cat(sprintf("  K=2-only Entrez genes:      %d\n", length(k2_eid)))

  # GO Biological Process
  cat("\n  -- GO:BP enrichment --\n")
  ego_bp <- tryCatch(enrichGO(
    gene          = k2_eid,
    universe      = bg_eid,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05,
    readable      = TRUE
  ), error = function(e) { cat("    GO error:", conditionMessage(e), "\n"); NULL })
  if (!is.null(ego_bp) && nrow(as.data.frame(ego_bp)) > 0) {
    print(head(as.data.frame(ego_bp), 15)[, c("Description","GeneRatio",
                                              "p.adjust","Count")],
          row.names = FALSE)
    saveRDS(ego_bp, file.path(RES_DIR, sprintf("tissue_%s_gobp.rds", TISSUE)))
  } else {
    cat("    No GO:BP terms enriched at q < 0.05\n")
  }

  # KEGG (requires internet)
  cat("\n  -- KEGG enrichment --\n")
  ek <- tryCatch(enrichKEGG(
    gene          = k2_eid,
    universe      = bg_eid,
    organism      = "hsa",
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05
  ), error = function(e) { cat("    KEGG error:", conditionMessage(e), "\n"); NULL })
  if (!is.null(ek) && nrow(as.data.frame(ek)) > 0) {
    ek_r <- setReadable(ek, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    print(head(as.data.frame(ek_r), 15)[, c("Description","GeneRatio",
                                            "p.adjust","Count")],
          row.names = FALSE)
    saveRDS(ek_r, file.path(RES_DIR, sprintf("tissue_%s_kegg.rds", TISSUE)))
  } else {
    cat("    No KEGG pathways enriched at q < 0.05\n")
  }
}

# ---- 7. Summary --------------------------------------------------------------
cat("\n[6/6] Files written:\n")
files <- list.files(RES_DIR, pattern = sprintf("^tissue_%s_", TISSUE),
                    full.names = TRUE)
for (f in files) cat("  ", f, "  (", file.size(f), " B)\n", sep = "")

cat("\nDone.\n")
