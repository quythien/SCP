#' ============================================================
#' fig1_paired_sims.R - Re-run Fig 1A/1B Panel A with paired sigma
#'
#' Produces paired-sampling RDS files used for the Panel A marginal
#' power curves only. Panels B and C continue to use the cached
#' unpaired RDS files as effect-size-stratified visualizations.
#'
#' Output:
#'   output/single_cohort/results/single_cohort_power_GTEx_AdrenalGland_paired_<ts>.rds
#'   output/single_cohort/results/single_cohort_power_GTEx_Liver_paired_<ts>.rds
#' ============================================================

setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
old_wd <- setwd("code"); source("setup.R"); setwd(old_wd)

GLOBAL_SEED <- 2025L
SAMPLE_SIZES <- c(20, 30, 40, 50, 60, 70, 80, 100, 120, 150, 200)
NSIMS  <- 100L
NCORES <- as.integer(Sys.getenv("MC_CORES", unset = "12"))

out_dir_res <- "output/single_cohort/results"
dir.create(out_dir_res, recursive = TRUE, showWarnings = FALSE)
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

gtex_path <- "/home/qtp1/Projects/Collaborative/GTEXdata/TOD Xiangning/data/CPM/CPM.all.norm.RData"
cat(sprintf("Loading GTEx: %s\n", gtex_path))
load(gtex_path)

extract_tissue <- function(tissue_name) {
  df   <- CPM.all.norm[[tissue_name]]
  ids  <- as.character(colnames(df))
  hhmm <- sapply(strsplit(ids, "\\."), function(x) if (length(x) >= 3) x[3] else NA)
  hrs  <- as.numeric(substr(hhmm, 1, 2)) + as.numeric(substr(hhmm, 3, 4)) / 60
  ok   <- !is.na(hrs)
  list(expr = as.matrix(df[, ok]), times = hrs[ok], n = sum(ok))
}

tissues <- list(
  AdrenalGland = "Adrenal Gland",
  Liver        = "Liver"
)

for (tname in names(tissues)) {
  cat(sprintf("\n=== %s (paired) ===\n", tname))
  t <- extract_tissue(tissues[[tname]])

  set.seed(GLOBAL_SEED)
  bio <- estCircadianParam(t$expr, t$times,
                            paired_sigma = TRUE,
                            verbose = TRUE)
  cat(sprintf("  pilot prop_rhythmic=%.3f, n_pilot=%d, paired_sigma=TRUE\n",
              bio$prop_rhythmic, length(bio$cts)))

  design <- CircadianDesignOptions(
    sample_sizes = SAMPLE_SIZES,
    nsims        = NSIMS,
    design       = "passive",
    cts          = bio$cts
  )

  # Standard bin convention (matches original Fig 1 cache). This RDS
  # contributes Panel A only; Panel B/C continues to come from the
  # cached unpaired RDS.
  analysis <- CircadianAnalysisOptions(
    alpha           = 0.05,
    p.adjust.method = "BH",
    r_strata        = makeAdaptiveRStrata(bio, bin_width = 0.25)
  )

  rds_path <- file.path(out_dir_res,
    sprintf("single_cohort_power_GTEx_%s_paired_%s.rds", tname, timestamp))

  set.seed(GLOBAL_SEED)
  res <- runSingleCohortPower(bio, design, analysis,
                               methods  = "DCP",
                               plot     = FALSE,
                               verbose  = TRUE,
                               mc.cores = NCORES)
  saveRDS(res, rds_path)
  cat(sprintf("  saved: %s\n", rds_path))
}

cat("\n=== Done ===\n")
