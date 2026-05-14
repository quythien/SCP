#' =======================================================================
#' suppFig_k_sweep_baboon.R - K-harmonic sweep for supplementary
#' =======================================================================
#'
#' K in {1, 2, 3} on the active baboon KIM pilot (GSE98965). KIM is the
#' most non-cosinor pilot in the diagnostic set (median omega 0.24), so
#' it gives K = 3 the largest possible empirical edge over K = 2. The
#' simulation answers the only K-choice question that the theoretical
#' variance bound does not resolve: does the marginal variance captured
#' by K = 3 (V_3 / V_2 ~ 1.10 at this omega) outweigh the two extra df
#' lost in the F-statistic denominator at finite N?
#'
#' Design: B = 12 (satisfies K = 3 identifiability B >= 7), single tissue,
#' N in {24, 48, 72, 144}, nsims = 30.
#'
#' Output: output/k_sweep/k_sweep_KIM.rds
#'         output/k_sweep/k_sweep_KIM.csv

setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
old_wd <- setwd("code"); source("setup.R"); setwd(old_wd)

GLOBAL_SEED <- 2025L
TISSUE      <- "KIM"
NSIMS       <- 30L
NGENES      <- 2000L
TOP_K_FMM   <- 300L
K_GRID      <- c(1L, 2L, 3L)
N_GRID      <- c(24L, 48L, 72L, 144L)
B_FIXED     <- 12L     # satisfies B >= 2K+1 = 7 for K = 3
N_CORES     <- as.integer(Sys.getenv("MC_CORES", unset = "8"))

out_dir <- "output/k_sweep"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ----- Load cached KIM pilot (from active_vs_passive run) -----
bio_path <- sprintf("output/active_vs_passive/results/bio_baboon_%s_NG%d_K%d.rds",
                    TISSUE, NGENES, TOP_K_FMM)
if (!file.exists(bio_path)) stop("Missing cached pilot: ", bio_path)
bio <- readRDS(bio_path)
cat(sprintf("Loaded KIM pilot: eta=%.3f  omega_med=%.3f  prop_rhythmic=%.3f\n",
            bio$diagnostics$beta_hat,
            median(bio$diagnostics$omega_emp),
            bio$prop_rhythmic))

analysis <- CircadianAnalysisOptions(alpha = 0.05, p.adjust.method = "BH",
                                     fdr_thresholds = 0.05)

# ----- Run K x N grid -----
results <- list()
for (K in K_GRID) {
  pwr <- vapply(N_GRID, function(N) {
    m   <- N %/% B_FIXED
    cts <- rep(seq(0, 24 * (1 - 1/B_FIXED), length.out = B_FIXED), each = m)
    design_i <- CircadianDesignOptions(sample_sizes = N, nsims = NSIMS,
                                       design = "active", cts = cts,
                                       B_values = B_FIXED)
    method <- if (K == 1L) "DCP" else "FMM"
    set.seed(GLOBAL_SEED + 1000L * K + N)
    r <- runSimsSingleCohort(bio, design_i, analysis,
                             method = method, K = K,
                             mc.cores = N_CORES, verbose = FALSE)
    mean(r$marginal_power, na.rm = TRUE)
  }, numeric(1))
  results[[as.character(K)]] <- data.frame(K = K, N = N_GRID, power = pwr)
  cat(sprintf("K=%d  %s\n", K,
              paste(sprintf("N=%d:%.3f", N_GRID, pwr), collapse = "  ")))
}

results_df <- do.call(rbind, results)
rownames(results_df) <- NULL

# ----- Save -----
saveRDS(list(results = results_df,
             bio_eta = bio$diagnostics$beta_hat,
             omega_med = median(bio$diagnostics$omega_emp),
             B = B_FIXED, NSIMS = NSIMS, NGENES = NGENES,
             tissue = TISSUE),
        file.path(out_dir, "k_sweep_KIM.rds"))
write.csv(results_df, file.path(out_dir, "k_sweep_KIM.csv"), row.names = FALSE)

cat("\n=== Results table ===\n")
print(reshape(results_df, idvar = "N", timevar = "K", direction = "wide"),
      row.names = FALSE)
cat("\nSaved: output/k_sweep/k_sweep_KIM.rds + .csv\n")
