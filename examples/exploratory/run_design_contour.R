#' =======================================================================
#' Optimal Design Contour Analysis
#' =======================================================================
#'
#' WHAT THIS DOES:
#'   Sweeps a (n, B) grid where n = total sample size per group and
#'   B = number of equally-spaced time bins (quota-sampled active design).
#'   Computes DR and DP power at each cell, then plots iso-power contour
#'   curves in (n, B) space.  A passive reference row (uncontrolled pilot
#'   TOD) is included for comparison.
#'
#' GRID:
#'   n  in {20, 40, 60, 80, 100, 120}
#'   B  in {4, 6, 8, 12, 24}          (cells with B > n are skipped)
#'   + passive reference at each n (pilot TOD distribution)
#'
#' USAGE:
#'   Rscript examples/run_design_contour.R
#'
#' @author Thien Quy Pham

# =====================================================================
# SECTION 1: SETUP
# =====================================================================

setwd("/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
source_dir <- file.path(getwd(), "code")
old_wd <- setwd(source_dir)
source("setup.R")
setwd(old_wd)

source("code/design_contour.R")


# =====================================================================
# SECTION 2: LOAD PILOT DATA
# =====================================================================

cat("Loading pilot expression data (PFC younger: BA11 + BA47)...\n")

COMBINED <- readRDS("/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Collaborative/Paper/Congruence/PNAS_aging/data/combined_data.rds")

expr_sample_names <- colnames(COMBINED$expr)
pheno_order       <- match(expr_sample_names, COMBINED$pheno$sample_name)
valid_samples     <- !is.na(pheno_order)
COMBINED$expr     <- COMBINED$expr[, valid_samples]
pheno_order       <- pheno_order[valid_samples]
pheno_data        <- COMBINED$pheno[pheno_order, ]
pheno_data$tod    <- if ("TOD.x" %in% colnames(pheno_data)) pheno_data$TOD.x else pheno_data$TOD.y
pheno_data$age_group_final <- if ("AgeGroup" %in% colnames(pheno_data)) pheno_data$AgeGroup else pheno_data$age_group

complete_samples <- !is.na(pheno_data$age_group_final) & !is.na(pheno_data$tod) &
  pheno_data$age_group_final %in% c("younger", "older")
pheno_clean  <- pheno_data[complete_samples, ]
younger_idx  <- pheno_clean$age_group_final == "younger"
expr_younger <- COMBINED$expr[, complete_samples][, younger_idx]
times_young  <- pheno_clean$tod[younger_idx]

cat(sprintf("  Expression matrix: %d genes x %d samples\n", nrow(expr_younger), ncol(expr_younger)))
cat(sprintf("  Reference TOD: n=%d younger subjects\n\n", length(times_young)))
rm(COMBINED, expr_sample_names, pheno_order, valid_samples, pheno_data,
   complete_samples, pheno_clean, younger_idx)


# =====================================================================
# SECTION 3: ESTIMATE PARAMETERS FROM PILOT DATA
# =====================================================================

cat("Estimating circadian parameters from pilot data...\n\n")

# DR-isolated scenario: 15% DR genes, no DP/DA (mirrors run_pipeline.R DR analysis)
opts_bio <- estCircadianParam(
  data          = expr_younger,
  times         = times_young,
  period        = 24,
  prop_DR       = 0.15,
  prop_DP       = 0.15,
  phase_diff    = c(-6, 6),
  amp_diff      = c(1, 1)
)
rm(expr_younger)

opts_bio <- updateBioOptions(opts_bio, ngenes = 5000)

# Increase effect size uniformly (scale amplitudes in both groups)
amp_mult <- 1.5
opts_bio$amplitude  <- opts_bio$amplitude  * amp_mult
if (!is.null(opts_bio$amplitude2)) {
  opts_bio$amplitude2 <- opts_bio$amplitude2 * amp_mult
}


# =====================================================================
# SECTION 4: OUTPUT DIRECTORY
# =====================================================================

run_tag  <- format(Sys.time(), "%Y%m%d_%H%M")
base_out <- file.path("output", paste0("design_contour_", run_tag))
fig_dir  <- file.path(base_out, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("Output directory: %s/\n\n", base_out))


# =====================================================================
# SECTION 5: RUN DESIGN CONTOUR GRID
# =====================================================================

cat("====================================================================\n")
cat("DESIGN CONTOUR: (n, B) grid sweep\n")
cat("====================================================================\n\n")

cat("Grid:\n")
cat("  n : 20, 40, 60, 80, 100, 120\n")
cat("  B : 4, 6, 8, 12, 24  (cells with B > n skipped)\n")
cat("  + passive reference row at each n (pilot TOD)\n")
cat("  nsims = 5 per cell (dry run)\n\n")

t_start <- proc.time()

grid_results <- powerDesignGrid(
  bio.opts   = opts_bio,
  n_grid     = c(20, 40, 60, 80, 100, 120),
  B_grid     = c(4, 6, 8, 12, 24),
  nsims      = 5,
  cts        = times_young,          # passive reference
  test_types = "DR",                 # DR only (prop_DP=0, no DP signal simulated)
  alpha      = 0.05,
  verbose    = TRUE
)

t_elapsed <- proc.time() - t_start
cat(sprintf("\nGrid sweep complete: %.1f minutes\n\n", t_elapsed[3] / 60))

# Save raw results
grid_rds <- file.path(base_out, "design_contour_grid.rds")
saveRDS(grid_results, file = grid_rds)
cat(sprintf("Results saved: %s\n\n", grid_rds))

# Print summary
print(grid_results)


# =====================================================================
# SECTION 6: DESIGN EFFICIENCY OF PASSIVE PILOT
# =====================================================================

cat("\nDesign efficiency of pilot TOD distribution:\n")
deff_pilot <- computeDesignEfficiency(times_young, period = 24)
cat(sprintf("  d_mean  = %.3f  (always ~0.5)\n",   deff_pilot$d_mean))
cat(sprintf("  d_min   = %.3f  (worst-case gene)\n", deff_pilot$d_min))
cat(sprintf("  d_sd    = %.3f  (spread across genes)\n", deff_pilot$d_sd))
cat(sprintf("  d_eff   = %.3f  (mean - sd)\n\n",    deff_pilot$d_eff))


# =====================================================================
# SECTION 7: CONTOUR PLOTS
# =====================================================================

cat("Generating contour plots...\n")

# --- DR contour ---
dr_fig <- file.path(fig_dir, "design_contour_DR.pdf")
pdf(dr_fig, width = 8, height = 6)
plotDesignContour(
  grid_results  = grid_results,
  test          = "DR",
  target_powers = c(0.2, 0.4, 0.6, 0.8),
  show_pFDR     = FALSE,
  passive_B     = NULL,   # no horizontal line (passive rows shown separately)
  title         = "DR Power: iso-power curves in (n, B) space\n(active design, fixed times)"
)
dev.off()
cat(sprintf("  DR contour: %s\n", dr_fig))

# --- DP contour ---
dp_fig <- file.path(fig_dir, "design_contour_DP.pdf")
pdf(dp_fig, width = 8, height = 6)
plotDesignContour(
  grid_results  = grid_results,
  test          = "DP",
  target_powers = c(0.2, 0.4, 0.6, 0.8),
  show_pFDR     = FALSE,
  passive_B     = NULL,
  title         = "DP Power: iso-power curves in (n, B) space\n(active design, fixed times)"
)
dev.off()
cat(sprintf("  DP contour: %s\n", dp_fig))

# --- Combined DR+DP side-by-side ---
combined_fig <- file.path(fig_dir, "design_contour_combined.pdf")
pdf(combined_fig, width = 14, height = 6)
par(mfrow = c(1, 2))
plotDesignContour(grid_results, test = "DR",
                  target_powers = c(0.2, 0.4, 0.6, 0.8),
                  title = "DR Power")
plotDesignContour(grid_results, test = "DP",
                  target_powers = c(0.2, 0.4, 0.6, 0.8),
                  title = "DP Power")
dev.off()
cat(sprintf("  Combined:   %s\n\n", combined_fig))

# --- Design factor curves: d(phi) vs gene phase (explains WHY B doesn't matter) ---
dfactor_fig <- file.path(fig_dir, "design_factor_curves.pdf")
pdf(dfactor_fig, width = 8, height = 5)
plotDesignFactor(
  pilot_times = times_young,
  B_active    = c(4, 6, 12, 24),
  period      = 24,
  title       = "Design factor d(\u03c6): active (any B) vs passive pilot TOD"
)
dev.off()
cat(sprintf("  Design factor: %s\n\n", dfactor_fig))

# --- Line plot: active vs passive (clearest comparison when B barely matters) ---
lines_dr_fig <- file.path(fig_dir, "design_lines_DR.pdf")
pdf(lines_dr_fig, width = 7, height = 5)
plotDesignLines(grid_results, test = "DR",
                title = "DR Power: active (equally-spaced) vs passive pilot",
                target_power = 0.8)
dev.off()
cat(sprintf("  DR lines:   %s\n", lines_dr_fig))

lines_dp_fig <- file.path(fig_dir, "design_lines_DP.pdf")
pdf(lines_dp_fig, width = 7, height = 5)
plotDesignLines(grid_results, test = "DP",
                title = "DP Power: active (equally-spaced) vs passive pilot",
                target_power = 0.8)
dev.off()
cat(sprintf("  DP lines:   %s\n\n", lines_dp_fig))

# --- (N, m) contour: sample size vs replicates per time point ---
nm_fig <- file.path(fig_dir, "design_contour_NM.pdf")
pdf(nm_fig, width = 8, height = 6)
plotDesignContourNM(grid_results, test = "DR")
dev.off()
cat(sprintf("  NM contour: %s\n\n", nm_fig))


# =====================================================================
# SECTION 8: PASSIVE vs ACTIVE COMPARISON TABLE
# =====================================================================

cat("====================================================================\n")
cat("PASSIVE vs ACTIVE (B=24) COMPARISON\n")
cat("====================================================================\n\n")

passive_rows <- grid_results[!is.na(grid_results$design_type) &
                               grid_results$design_type == "passive", ]
active_B24   <- grid_results[!is.na(grid_results$B) & grid_results$B == 24, ]

cat(sprintf("%-6s | %-14s %-14s | %-14s %-14s\n",
            "n", "DR_pass(TPow)", "DR_B=24(TPow)", "DP_pass(TPow)", "DP_B=24(TPow)"))
cat(paste0(rep("-", 50), collapse = ""), "\n")
for (nv in c(20, 40, 60, 80, 100, 120)) {
  pr <- passive_rows[passive_rows$n == nv, ]
  ar <- active_B24[active_B24$n == nv, ]
  if (nrow(pr) == 0 || nrow(ar) == 0) next
  cat(sprintf("%-6d | %7.1f%% %7.1f%% | %7.1f%% %7.1f%%\n",
              nv,
              100 * pr$EDR_DR, 100 * ar$EDR_DR,
              100 * pr$EDR_DP, 100 * ar$EDR_DP))
}

cat(sprintf("\nPilot passive d_eff = %.3f vs active d = 0.5\n", deff_pilot$d_eff))
cat(sprintf("Passive pilot needs ~%.1fx more samples than active B=24 for same DR power\n\n",
            0.5 / max(deff_pilot$d_eff, 0.01)))


# =====================================================================
# SECTION 9: OPTIMAL DESIGN SELECTION (utility maximization)
# =====================================================================

cat("====================================================================\n")
cat("OPTIMAL DESIGN SELECTION\n")
cat("====================================================================\n\n")

# DR-only objective (w_DR=1, w_DP=0)
cat("--- DR objective (w_DR=1, w_DP=0) ---\n")
util_dr <- computeDesignUtility(
  grid_results,
  w_DR         = 1,
  w_DP         = 0,
  lambda       = 1,
  alpha        = 0.05,
  target_power = 0.8
)

# Save utility table
util_rds <- file.path(base_out, "design_utility_DR.rds")
saveRDS(util_dr, file = util_rds)
cat(sprintf("\nUtility table saved: %s\n\n", util_rds))


# =====================================================================
# SECTION 10: WRAP-UP
# =====================================================================

t_total <- proc.time() - t_start
cat("====================================================================\n")
cat("DESIGN CONTOUR + UTILITY OPTIMIZATION COMPLETE\n")
cat("====================================================================\n\n")
cat(sprintf("Total runtime: %.1f minutes\n", t_total[3] / 60))
cat(sprintf("Output:\n"))
cat(sprintf("  Results: %s\n", grid_rds))
cat(sprintf("  DR plot: %s\n", dr_fig))
cat(sprintf("  DP plot: %s\n", dp_fig))
cat(sprintf("  Combined: %s\n", combined_fig))
cat("\n")
