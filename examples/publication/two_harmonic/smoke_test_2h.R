# =============================================================================
# Smoke test for the two-harmonic (K=2) cosinor refactor.
#
# Verifies:
#   1. estCircadianParam2H() runs on a real pilot (baboon LUN, ~500 genes).
#   2. The returned CircadianBioOptions carries amplitude2 / phase2 / paired_2h.
#   3. runSimsSingleCohort() dispatches to the new 2-harmonic generator
#      (gated on has_2h_per_gene = TRUE) and method = "FMM", K = 2L runs
#      end-to-end without an explicit `harmonics` argument.
#   4. Empirical power > 0 at N = 48.
#
# Usage (from project root):
#   SMOKE_TEST=true Rscript examples/publication/two_harmonic/smoke_test_2h.R
# =============================================================================

t0_total <- Sys.time()

# --- load package via code/setup.R (uses working directory) -----------------
PROJ <- "/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim"
old_wd <- getwd()
setwd(file.path(PROJ, "code"))
source("setup.R")
setwd(old_wd)

# --- load baboon LUN pilot (4938 genes x 12 samples) ------------------------
load(file.path(PROJ, "data/CAMO_PRC_hmb.RData"))
expr_full <- baboon_withTOD$baboon$LUN
times     <- baboon_withTOD$tod$LUN
cat(sprintf("\nPilot: baboon LUN, %d genes x %d samples\n",
            nrow(expr_full), ncol(expr_full)))

# Sub-sample to ~500 genes for a fast smoke test.
set.seed(2026)
G_smoke <- 500L
keep    <- sort(sample.int(nrow(expr_full), min(G_smoke, nrow(expr_full))))
expr    <- expr_full[keep, , drop = FALSE]
cat(sprintf("Smoke subset: %d genes x %d samples\n\n",
            nrow(expr), ncol(expr)))

# --- (1) fit two-harmonic pilot ---------------------------------------------
cat("--- (1) estCircadianParam2H ---\n")
t0 <- Sys.time()
bio_2h <- estCircadianParam2H(
  data            = expr,
  times           = times,
  period          = 24,
  min_rhythm_pval = 0.05,    # slightly looser to ensure top-K is populated
  top_k           = 200L,
  prop_DR         = 0,
  prop_DP         = 0,
  prop_DM         = 0,
  verbose         = TRUE
)
elapsed_est <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("\nestCircadianParam2H runtime: %.2f s\n\n", elapsed_est))

# --- (2) structure check ----------------------------------------------------
cat("--- (2) structure check ---\n")
stopifnot(inherits(bio_2h, "CircadianBioOptions"))
stopifnot(isTRUE(bio_2h$paired_2h))
stopifnot(!is.null(bio_2h$amplitude2),
          length(bio_2h$amplitude2) == length(bio_2h$amplitude))
stopifnot(!is.null(bio_2h$phase2),
          length(bio_2h$phase2) == length(bio_2h$amplitude))
stopifnot(length(bio_2h$sigma_rhythmic) == length(bio_2h$amplitude))
cat(sprintf("  paired_2h:    %s\n", bio_2h$paired_2h))
cat(sprintf("  L_pilot:      %d\n", length(bio_2h$amplitude)))
cat(sprintf("  range(A1):    [%.3f, %.3f]\n",
            min(bio_2h$amplitude), max(bio_2h$amplitude)))
cat(sprintf("  range(A2):    [%.3f, %.3f]\n",
            min(bio_2h$amplitude2), max(bio_2h$amplitude2)))
cat(sprintf("  range(phi1):  [%.2f, %.2f] h\n",
            min(bio_2h$phase), max(bio_2h$phase)))
cat(sprintf("  range(phi2):  [%.2f, %.2f] h\n",
            min(bio_2h$phase2), max(bio_2h$phase2)))
stopifnot(max(bio_2h$phase2) < 12 + 1e-6)  # phi2 wrapped to [0, 12)
cat("  STRUCTURE: OK\n\n")

# --- (3) runSimsSingleCohort dispatch test ----------------------------------
cat("--- (3) runSimsSingleCohort, method = 'FMM', K = 2 ---\n")
# Configure for a single small replicate at N = 48.
bio_2h$ngenes        <- 500L  # match smoke subset
bio_2h$prop_rhythmic <- max(bio_2h$prop_rhythmic, 0.05)  # ensure rhythmic budget > 0

dopt <- CircadianDesignOptions(
  sample_sizes = 48L,
  nsims        = 5L,   # tiny replicate budget for a smoke test
  design       = "active"
)
aopt <- CircadianAnalysisOptions(alpha = 0.05)

t1 <- Sys.time()
res <- runSimsSingleCohort(
  bio.opts      = bio_2h,
  design.opts   = dopt,
  analysis.opts = aopt,
  method        = "FMM",
  K             = 2L,
  verbose       = TRUE,
  mc.cores      = 1L
)
elapsed_sim <- as.numeric(difftime(Sys.time(), t1, units = "secs"))

power_mean <- mean(res$marginal_power, na.rm = TRUE)
fdr_mean   <- mean(res$marginal_FDR,   na.rm = TRUE)
r_smp      <- unlist(res$r_values_list[[1]])
r_pos      <- r_smp[r_smp > 0]

cat(sprintf("\n  Sim runtime:       %.2f s\n", elapsed_sim))
cat(sprintf("  Mean power (N=48): %.3f\n", power_mean))
cat(sprintf("  Mean FDR (N=48):   %.3f\n", fdr_mean))
cat(sprintf("  r_values summary (rhythmic only, n=%d): median=%.3f  IQR=[%.3f, %.3f]\n",
            length(r_pos),
            stats::median(r_pos),
            stats::quantile(r_pos, 0.25),
            stats::quantile(r_pos, 0.75)))

stopifnot(power_mean > 0)
cat("  POWER > 0: OK\n\n")

elapsed_total <- as.numeric(difftime(Sys.time(), t0_total, units = "secs"))
cat(sprintf("=== SMOKE TEST PASSED (total runtime: %.2f s) ===\n",
            elapsed_total))
