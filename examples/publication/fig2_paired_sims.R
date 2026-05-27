#' ============================================================
#' fig2_paired_sims.R - Re-run Fig 2 differential power with paired sigma
#'
#' Re-estimates the two-group ADR-vs-LIV pilot with paired_sigma = TRUE
#' and re-runs runDifferentialPower. Output RDS replaces the unpaired
#' Panel A source for the trimmed 1 x 3 Fig 2 layout.
#'
#' Output:
#'   data/gtex_adr_vs_liv_pilot_paired.rds
#'   output/differential/results/diff_power_ADR_vs_LIV_paired_<ts>.rds
#' ============================================================

setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
old_wd <- setwd("code"); source("setup.R"); setwd(old_wd)

GLOBAL_SEED  <- 2025L
SAMPLE_SIZES <- c(20, 30, 40, 50, 60, 70, 80, 100, 120, 150, 200)
NSIMS  <- 100L
NCORES <- as.integer(Sys.getenv("MC_CORES", unset = "12"))

set.seed(GLOBAL_SEED)
out_dir <- "output/differential/results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
cat(sprintf("=== Paired Fig 2 (ADR vs LIV) [%s] ===\n", timestamp))

# Build paired two-group pilot. The bundled pilot RDS lets this script
# reproduce Figure 2 without the access-controlled raw GTEx matrix; the raw
# rebuild below runs only if the cached pilot is absent.
pilot_rds <- "data/gtex_adr_vs_liv_pilot_paired.rds"
if (file.exists(pilot_rds)) {
  cat(sprintf("Loading cached paired pilot: %s\n", pilot_rds))
  bio <- readRDS(pilot_rds)
} else {
  cat("Cached pilot not found; rebuilding from raw GTEx\n")
  cat("(requires the access-controlled CPM.all.norm.RData; see README data provenance).\n")
  gtex_path <- "/home/qtp1/Projects/Collaborative/GTEXdata/TOD Xiangning/data/CPM/CPM.all.norm.RData"
  load(gtex_path)
  extract_tissue <- function(tissue_name) {
    df   <- CPM.all.norm[[tissue_name]]
    ids  <- as.character(colnames(df))
    hhmm <- sapply(strsplit(ids, "\\."), function(x) if (length(x) >= 3) x[3] else NA)
    hrs  <- as.numeric(substr(hhmm, 1, 2)) + as.numeric(substr(hhmm, 3, 4)) / 60
    ok   <- !is.na(hrs)
    list(expr = as.matrix(df[, ok]), times = hrs[ok], n = sum(ok))
  }
  adr <- extract_tissue("Adrenal Gland")
  liv <- extract_tissue("Liver")
  common_genes <- intersect(rownames(adr$expr), rownames(liv$expr))
  expr_adr <- adr$expr[common_genes, ]
  expr_liv <- liv$expr[common_genes, ]
  set.seed(GLOBAL_SEED)
  bio <- estCircadianParamTwoGroup(
    data_1  = expr_adr, data_2  = expr_liv,
    times_1 = adr$times, times_2 = liv$times,
    paired_sigma = TRUE,
    verbose = TRUE
  )
  saveRDS(bio, pilot_rds)
  cat(sprintf("Paired pilot saved: %s\n", pilot_rds))
}

cat(sprintf("\nPilot: n1=%d n2=%d DR=%.3f DP=%.3f DM=%.3f\n",
            length(bio$cts), length(bio$cts2),
            bio$prop_DR, bio$prop_DP, bio$prop_DM))

design <- CircadianDesignOptions(
  sample_sizes = SAMPLE_SIZES,
  nsims        = NSIMS,
  design       = "passive",
  cts          = bio$cts,
  test_types   = c("DR", "DP", "DM")
)

analysis <- CircadianAnalysisOptions(
  alpha           = 0.05,
  p.adjust.method = "BH"
)

rds_path <- file.path(out_dir,
  sprintf("diff_power_ADR_vs_LIV_paired_%s.rds", timestamp))

cat(sprintf("\nRunning runDifferentialPower (NSIMS=%d, NCORES=%d)...\n",
            NSIMS, NCORES))
t0 <- Sys.time()
set.seed(GLOBAL_SEED)
res <- runDifferentialPower(bio, design, analysis,
                             methods    = "DCP",
                             test_types = c("DR", "DP", "DM"),
                             plot       = FALSE,
                             verbose    = TRUE,
                             mc.cores   = NCORES)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
saveRDS(res, rds_path)
cat(sprintf("Elapsed: %.1f min\n", elapsed))
cat(sprintf("Saved: %s\n", rds_path))
cat("\n=== Done ===\n")
