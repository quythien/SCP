#!/usr/bin/env Rscript
# Fig 6 Option B: Baboon KIM (K=2 cosinor) - matches the OLD Fig 6 style
setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
sink(stderr())
source("code/utils.R"); source("code/options.R"); source("code/simulation.R")
source("code/estimation.R"); source("code/bootstrap_sim.R"); source("code/detection.R")
sink()
suppressPackageStartupMessages(library(parallel))

# Build K=2 pilot from Baboon KIM expression
df <- read.csv("/home/qtp1/Projects/Circadian/Baboon/data/GSE98965_baboon_tissue_expression_FPKM.csv",
               stringsAsFactors = FALSE)
expr_mat <- as.matrix(df[, !colnames(df) %in% c("EnsemblID","Symbol")])
rownames(expr_mat) <- df$Symbol
cn_clean <- gsub("_ZT", ".ZT", colnames(expr_mat))
parts <- strsplit(cn_clean, "\\."); tis <- sapply(parts,`[`,1)
zts <- as.numeric(sub("ZT","",sapply(parts,`[`,2)))
keep_genes <- rowSums(expr_mat > 1, na.rm = TRUE) > 100
expr_mat <- expr_mat[keep_genes, ]
expr_log <- log2(expr_mat + 1)

idx <- which(tis == "KIM"); ts <- zts[idx]
ord <- order(ts); ts <- ts[ord]; idx <- idx[ord]
sub <- expr_log[, idx, drop = FALSE]
first <- !duplicated(rownames(sub)) & nzchar(rownames(sub)) & !is.na(rownames(sub))
sub <- sub[first, , drop = FALSE]
cat(sprintf("Baboon KIM: %d genes x %d samples\n", nrow(sub), ncol(sub)))
psi <- estCircadianParam2H(data = sub, times = ts, period = 24,
                            min_rhythm_pval = 0.05, top_k = 300, verbose = FALSE)
cat(sprintf("KIM K=2 pilot: top_k=%d, A2/A1 med=%.2f, prop_rhythmic=%.3f\n",
            psi$diagnostics$top_k_used,
            psi$diagnostics$A2_over_A1_med, psi$prop_rhythmic))
saveRDS(psi, "output/two_harmonic/results/pilot_2h_BaboonKIM.rds")

PERIOD <- 24; OMEGA0 <- 2*pi/PERIOD
DISPLAY_N    <- c(12L, 16L, 20L, 24L, 32L, 48L, 72L, 96L, 120L, 144L)
B_GRID_K1    <- c(3L, 6L, 8L, 12L, 24L)
B_GRID_K2    <- c(6L, 8L, 12L, 24L)
NSIMS <- 100; NCORES <- 12; SEED_BASE <- 20260526
NGENES_SIM <- 2000

active_times <- function(N, B) {
  if (N < B) return(sort(sample(seq(0, PERIOD - PERIOD/B, length.out = B), N, replace = FALSE)))
  m  <- N %/% B
  ts <- rep(seq(0, PERIOD - PERIOD/B, length.out = B), each = m)
  if (length(ts) < N) ts <- c(ts, sample(seq(0, PERIOD - PERIOD/B, length.out = B),
                                          N - length(ts), replace = TRUE))
  if (length(ts) > N) ts <- ts[seq_len(N)]
  ts
}
sim_one <- function(N, B, seed, K) {
  ts <- active_times(N, B)
  set.seed(seed)
  psi_local <- psi; psi_local$ngenes <- NGENES_SIM
  sim <- simCircadianSingleCohort2H(psi_local, ts, seed = seed)
  if (K == 1L) {
    Xs <- cbind(1, cos(OMEGA0*ts), sin(OMEGA0*ts)); df_full <- length(ts) - 3L
  } else {
    Xs <- cbind(1, cos(OMEGA0*ts), sin(OMEGA0*ts),
                cos(2*OMEGA0*ts), sin(2*OMEGA0*ts)); df_full <- length(ts) - 5L
  }
  if (df_full <= 0) return(NULL)
  G <- nrow(sim$expr); pv <- numeric(G)
  for (g in seq_len(G)) {
    y <- sim$expr[g,]; R0 <- sum((y-mean(y))^2)
    if (R0 <= 0) { pv[g] <- 1; next }
    f <- lm.fit(Xs, y); R <- sum(f$residuals^2)
    if (R <= 0) { pv[g] <- 1; next }
    pv[g] <- pf(((R0-R)/(K*2)) / (R/df_full), K*2, df_full, lower.tail = FALSE)
  }
  list(pvalue = pv, is_rhythmic = sim$is_rhythmic)
}

cat("== K=1 power sweep ==\n")
pwr_K1 <- matrix(NA_real_, nrow=length(DISPLAY_N), ncol=length(B_GRID_K1),
                  dimnames=list(paste0("N=",DISPLAY_N), paste0("B=",B_GRID_K1)))
for (j in seq_along(DISPLAY_N)) {
  for (k in seq_along(B_GRID_K1)) {
    N <- DISPLAY_N[j]; B <- B_GRID_K1[k]
    if (N < B) next
    sims <- mclapply(seq_len(NSIMS), function(i)
      sim_one(N, B, SEED_BASE + (N*1000L+B+i)*7919L, K = 1L),
      mc.cores = NCORES, mc.preschedule = FALSE)
    sims <- sims[vapply(sims, is.list, logical(1))]
    if (length(sims) > 0) {
      pwr_K1[j,k] <- mean(sapply(sims, function(s) {
        q <- p.adjust(s$pvalue, "BH")
        sum(q <= 0.05 & s$is_rhythmic) / max(1, sum(s$is_rhythmic))
      }), na.rm = TRUE)
    }
  }
  cat(sprintf("  N=%d  K=1: %s\n", DISPLAY_N[j],
              paste(sprintf("%.0f%%", 100*pwr_K1[j,]), collapse = " / ")))
}

cat("\n== K=2 power sweep ==\n")
pwr_K2 <- matrix(NA_real_, nrow=length(DISPLAY_N), ncol=length(B_GRID_K2),
                  dimnames=list(paste0("N=",DISPLAY_N), paste0("B=",B_GRID_K2)))
for (j in seq_along(DISPLAY_N)) {
  for (k in seq_along(B_GRID_K2)) {
    N <- DISPLAY_N[j]; B <- B_GRID_K2[k]
    if (N < B) next
    sims <- mclapply(seq_len(NSIMS), function(i)
      sim_one(N, B, SEED_BASE + (N*1000L+B+i)*7919L, K = 2L),
      mc.cores = NCORES, mc.preschedule = FALSE)
    sims <- sims[vapply(sims, is.list, logical(1))]
    if (length(sims) > 0) {
      pwr_K2[j,k] <- mean(sapply(sims, function(s) {
        q <- p.adjust(s$pvalue, "BH")
        sum(q <= 0.05 & s$is_rhythmic) / max(1, sum(s$is_rhythmic))
      }), na.rm = TRUE)
    }
  }
  cat(sprintf("  N=%d  K=2: %s\n", DISPLAY_N[j],
              paste(sprintf("%.0f%%", 100*pwr_K2[j,]), collapse = " / ")))
}

saveRDS(list(pwr_K1 = pwr_K1, pwr_K2 = pwr_K2,
              N = DISPLAY_N, B_K1 = B_GRID_K1, B_K2 = B_GRID_K2,
              source = "Baboon KIM K=2"),
        "output/two_harmonic/results/fig6_optionB_BaboonKIM.rds")
cat("\nSaved: output/two_harmonic/results/fig6_optionB_BaboonKIM.rds\n")
