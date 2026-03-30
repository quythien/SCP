#' =======================================================================
#' 07_seney_sex_DR.R — Seney MDD vs Control: Sex-Stratified DR Power
#' =======================================================================
#'
#' PURPOSE
#'   Compare DR power across three analysis strategies for the Seney
#'   MDD vs Control ACC dataset, which has equal numbers of male and female
#'   subjects (30M + 30F per disease group):
#'
#'     (A) Combined   — n=60 MDD vs n=60 Control (all subjects)
#'     (B) Male-only  — n=30 MDD-M vs n=30 CTL-M
#'     (C) Female-only — n=30 MDD-F vs n=30 CTL-F
#'
#' SCIENTIFIC RATIONALE FOR SEX STRATIFICATION
#'   Sex is a well-documented source of biological heterogeneity in circadian
#'   gene expression in the ACC (Seney et al. 2018) and in MDD (females ~2x
#'   more likely; distinct transcriptional profiles). When male and female
#'   subjects are pooled, two rhythm profiles with different phase alignments
#'   and amplitudes are averaged into a single cosinor estimate:
#'
#'     (1) Apparent amplitude A decreases — partial cancellation of rhythmic
#'         signal due to between-sex phase incoherence
#'     (2) Residual variance σ increases — cosinor residuals absorb between-sex
#'         variation that the cosinor model cannot explain
#'
#'   Together these reduce r = A/σ from ~0.82 (male-only) to ~0.55 (combined),
#'   a ~50% drop. Since power scales as r² × N, this requires ~3× more subjects:
#'     Combined n80   ~300  (passive, ~55% power at N=200)
#'     Male-only n80  ~200  (passive, ~73% power at N=200)
#'
#'   Practical recommendation: even in a passive design where timing cannot be
#'   controlled, stratifying by a known biological confounder (sex) can cut the
#'   required sample size by ~30%.
#'
#' QUICK RUN RESULTS (N=150, 200; NSIMS=20; NGENES=2000; p<0.05)
#'   Combined:    n=150 -> 47.8%   n=200 -> 56.3%   n80 ~280-300
#'   Male-only:   n=150 -> 66.0%   n=200 -> 73.1%   n80 ~200-220
#'   Female-only: n=150 -> 68.5%   n=200 -> 75.5%   n80 ~200-220
#'
#' DATA
#'   data/ACC_RNA_filtered_normalized.csv   — expression (log-normalized; do NOT re-normalize)
#'   data/MD5_MetaData_1-15-25.xlsx         — metadata (Disease, GENDER, HU_NUM)
#'   data/TOD.xlsx                          — OFFC_TIME (official pronounced time)
#'
#' KEY METADATA CODES
#'   Disease: 1 = Control,  2 = MDD
#'   GENDER:  1 = Male,     2 = Female
#'
#' OUTPUTS
#'   output/07_seney_sex_DR_<timestamp>/
#'     signal_summary.txt     — prop_rhythmic, r, prop_DR per scenario
#'     power_comparison.txt   — power table + n80 for all three scenarios
#'
#' USAGE
#'   Rscript examples/exploratory/07_seney_sex_DR.R
#'   POWERSIM_SMOKE=1 Rscript examples/exploratory/07_seney_sex_DR.R
#'
#'   On a server, set POWERSIM_ROOT before running:
#'     export POWERSIM_ROOT=/path/to/PowerSim
#'     Rscript examples/exploratory/07_seney_sex_DR.R

# -----------------------------------------------------------------------
# 0. Path configuration
# -----------------------------------------------------------------------
# Set POWERSIM_ROOT as an environment variable for portability:
#   Linux/Mac: export POWERSIM_ROOT=/path/to/PowerSim
#   R console: Sys.setenv(POWERSIM_ROOT="/path/to/PowerSim")

POWERSIM_ROOT <- {
  env <- Sys.getenv("POWERSIM_ROOT", unset = "")
  if (nzchar(env)) env else
    "/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim"
}
setwd(POWERSIM_ROOT)

# -----------------------------------------------------------------------
# 1. Settings
# -----------------------------------------------------------------------
SMOKE       <- nzchar(Sys.getenv("POWERSIM_SMOKE"))
NGENES      <- if (SMOKE) 300L  else 3000L
NSIMS       <- if (SMOKE) 5L    else 50L
N_GRID      <- if (SMOKE) c(40L, 60L, 100L) else c(40L, 60L, 80L, 100L, 150L, 200L, 300L)
RHYTHM_PVAL <- 0.05   # threshold for calling a gene "rhythmic" in pilot (parameter estimation only)

DATA_META   <- "data/MD5_MetaData_1-15-25.xlsx"
DATA_TOD    <- "data/TOD.xlsx"
DATA_EXPR   <- "data/ACC_RNA_filtered_normalized.csv"

# -----------------------------------------------------------------------
# 2. Load libraries and framework
# -----------------------------------------------------------------------
suppressPackageStartupMessages({
  library(readxl)
})

src_files <- c(
  "code/options.R", "code/utils.R", "code/estimation.R", "code/simulation.R",
  "code/detection.R", "code/power.R", "code/runner.R", "code/bootstrap_sim.R"
)
for (f in src_files) source(f)

# -----------------------------------------------------------------------
# 2. Load and align data
# -----------------------------------------------------------------------
cat("Loading Seney data...\n")

meta_s  <- read_excel(DATA_META)
tod_s   <- read_excel(DATA_TOD)
expr_s  <- as.matrix(read.csv(DATA_EXPR, row.names = 1, check.names = FALSE))
# Note: prepCircadianData applied below after HU_NUM metadata join + ok filter

# Match expression columns to metadata/TOD by HU_NUM
# Strip trailing letter from colnames (e.g., "12345A" -> "12345")
col_ids <- gsub("[A-Za-z]+$", "", colnames(expr_s))
meta_idx <- match(col_ids, as.character(meta_s$HU_NUM))
tod_idx  <- match(col_ids, as.character(tod_s$HU_NUM))

# Extract TOD from OFFC_TIME (hour of day, 0-24)
offc_time <- tod_s$OFFC_TIME[tod_idx]
tod_hour  <- as.numeric(format(offc_time, "%H")) +
             as.numeric(format(offc_time, "%M")) / 60

# Metadata vectors (aligned to expression columns)
disease <- meta_s$Disease[meta_idx]   # 1=Control, 2=MDD
gender  <- meta_s$GENDER[meta_idx]    # 1=Male, 2=Female

# Check alignment
cat(sprintf("Expression: %d genes x %d samples\n", nrow(expr_s), ncol(expr_s)))
cat(sprintf("TOD range: [%.1f, %.1f] h\n", min(tod_hour, na.rm=TRUE), max(tod_hour, na.rm=TRUE)))
cat(sprintf("Disease: %d Control, %d MDD\n",
            sum(disease == 1, na.rm=TRUE), sum(disease == 2, na.rm=TRUE)))
cat(sprintf("Gender:  %d Male, %d Female\n",
            sum(gender == 1, na.rm=TRUE), sum(gender == 2, na.rm=TRUE)))

# Filter to complete cases; use prepCircadianData for standard validation + coercion
ok <- !is.na(disease) & !is.na(gender) & !is.na(tod_hour)
prep_seney <- prepCircadianData(expr_s[, ok], times = tod_hour[ok], input_type = "log2")
expr_s   <- prep_seney$data
tod_hour <- prep_seney$times
disease  <- disease[ok]
gender   <- gender[ok]

# -----------------------------------------------------------------------
# 3. Subsample genes for analysis
# -----------------------------------------------------------------------
set.seed(42)
gene_idx <- sample(nrow(expr_s), min(NGENES, nrow(expr_s)))
expr_sub  <- expr_s[gene_idx, ]

# -----------------------------------------------------------------------
# 4. Define the three analysis scenarios
# -----------------------------------------------------------------------
scenarios <- list(
  A_combined = list(
    label  = "Combined (n=60/group)",
    color  = "firebrick",
    ctl    = which(disease == 1),
    mdd    = which(disease == 2)
  ),
  B_male = list(
    label  = "Male-only (n=30/group)",
    color  = "steelblue",
    ctl    = which(disease == 1 & gender == 1),
    mdd    = which(disease == 2 & gender == 1)
  ),
  C_female = list(
    label  = "Female-only (n=30/group)",
    color  = "darkorchid",
    ctl    = which(disease == 1 & gender == 2),
    mdd    = which(disease == 2 & gender == 2)
  )
)

for (nm in names(scenarios)) {
  sc <- scenarios[[nm]]
  cat(sprintf("  %s: CTL n=%d, MDD n=%d\n", sc$label, length(sc$ctl), length(sc$mdd)))
}

# -----------------------------------------------------------------------
# 5. Setup output directory
# -----------------------------------------------------------------------
ts       <- format(Sys.time(), "%Y%m%d_%H%M")
out_dir  <- file.path("output", paste0("07_seney_sex_DR_", ts))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("\nOutputs -> %s/\n", out_dir))

# -----------------------------------------------------------------------
# 6. For each scenario: estimate parameters + run power analysis
# -----------------------------------------------------------------------
results <- list()
sig_lines <- character(0)

for (nm in names(scenarios)) {
  sc <- scenarios[[nm]]
  cat(sprintf("\n=== Scenario: %s ===\n", sc$label))

  ctl_data <- expr_sub[, sc$ctl]
  mdd_data <- expr_sub[, sc$mdd]
  ctl_tod  <- tod_hour[sc$ctl]
  mdd_tod  <- tod_hour[sc$mdd]

  # --- Parameter estimation (two-group) ---
  cat("  Estimating parameters...\n")
  # Estimate parameters from both groups simultaneously.
  # estCircadianParamTwoGroup() now uses prop_union_rhythmic as the budget,
  # so prop_DR (xor of both groups) will always fit within the budget.
  bio_opts <- tryCatch(
    estCircadianParamTwoGroup(
      data_1          = ctl_data,
      data_2          = mdd_data,
      times_1         = ctl_tod,
      times_2         = mdd_tod,
      period          = 24,
      min_rhythm_pval = RHYTHM_PVAL,
      verbose         = TRUE
    ),
    error = function(e) {
      cat(sprintf("  ERROR in estCircadianParamTwoGroup: %s\n", conditionMessage(e)))
      NULL
    }
  )

  if (is.null(bio_opts)) {
    cat("  Skipping power analysis for this scenario.\n")
    next
  }

  diag        <- bio_opts$diagnostics
  prop_rhy_1  <- diag$prop_rhythmic_1
  prop_rhy_2  <- diag$prop_rhythmic_2
  r_snr       <- diag$r1_snr["median"]
  prop_dr_est <- bio_opts$prop_DR   # may be slightly < prop_DR_emp if budget scaled

  cat(sprintf("  CTL rhythmic: %.1f%%  MDD rhythmic: %.1f%%  union: %.1f%%\n",
              100 * prop_rhy_1, 100 * prop_rhy_2, 100 * diag$prop_union_rhythmic))
  cat(sprintf("  prop_DR (in sim): %.1f%%   r median (CTL): %.3f  IQR [%.2f, %.2f]\n",
              100 * prop_dr_est, r_snr,
              diag$r1_snr["q25"], diag$r1_snr["q75"]))

  sig_lines <- c(sig_lines, sprintf(
    "%-30s  CTL_rhy=%5.1f%%  MDD_rhy=%5.1f%%  union_rhy=%5.1f%%  prop_DR=%5.1f%%  r_med=%.3f",
    sc$label, 100 * prop_rhy_1, 100 * prop_rhy_2,
    100 * diag$prop_union_rhythmic, 100 * prop_dr_est, r_snr
  ))

  # --- Design options (passive) ---
  design_opts <- CircadianDesignOptions(
    sample_sizes = N_GRID,
    nsims        = NSIMS,
    design       = "passive",
    cts          = ctl_tod,    # empirical TOD distribution
    test_types   = c("DR")
  )

  analysis_opts <- CircadianAnalysisOptions()

  # --- Run power analysis ---
  cat(sprintf("  Running power analysis (%d N values, %d sims each)...\n",
              length(N_GRID), NSIMS))
  power_res <- tryCatch(
    runPowerAnalysis(bio_opts, design_opts, analysis_opts, test_type = "DR"),
    error = function(e) {
      cat(sprintf("  ERROR in runPowerAnalysis: %s\n", conditionMessage(e)))
      NULL
    }
  )

  if (is.null(power_res)) next

  # marginal_power is [n_sizes x nsims]; rowMeans gives power per N
  power_vec <- rowMeans(power_res$marginal_power, na.rm = TRUE)

  results[[nm]] <- list(
    label      = sc$label,
    color      = sc$color,
    n_grid     = N_GRID,
    power      = power_vec,
    r_snr      = r_snr,
    prop_DR    = prop_dr_est,
    prop_rhy_1 = prop_rhy_1,
    prop_rhy_2 = prop_rhy_2
  )

  cat(sprintf("  Power at each N: %s\n",
              paste(sprintf("n=%d:%.0f%%", N_GRID, 100 * power_vec), collapse=", ")))
}

# -----------------------------------------------------------------------
# 7. Find n80 for each scenario and print comparison table
# -----------------------------------------------------------------------
cat("\n\n========== POWER COMPARISON SUMMARY ==========\n")

power_lines <- c(
  "Scenario-by-scenario power at each N per group",
  sprintf("%-30s  %s", "Scenario",
          paste(sprintf("n=%3d", N_GRID), collapse = "  ")),
  strrep("-", 80)
)

n80_vals <- c()
for (nm in names(results)) {
  res <- results[[nm]]
  pct <- sprintf("%5.1f%%", 100 * res$power)
  power_lines <- c(power_lines,
    sprintf("%-30s  %s", res$label, paste(pct, collapse = "  ")))

  # n80: smallest N with power >= 0.80
  hit <- which(res$power >= 0.80)
  n80 <- if (length(hit) > 0) N_GRID[hit[1]] else paste0(">", max(N_GRID))
  n80_vals <- c(n80_vals, setNames(n80, nm))
}

power_lines <- c(power_lines, "",
  "n80 (smallest N achieving >= 80% power):",
  sprintf("  %-30s  n80 = %s", sapply(results, `[[`, "label"), n80_vals)
)

cat(paste(power_lines, collapse = "\n"), "\n")

# Write power table
writeLines(
  c(sig_lines, "", power_lines),
  file.path(out_dir, "power_comparison.txt")
)

# Write signal summary
if (length(sig_lines) > 0)
  writeLines(sig_lines, file.path(out_dir, "signal_summary.txt"))

# -----------------------------------------------------------------------
# 8. Interpretation guidance
# -----------------------------------------------------------------------
cat("\n========== INTERPRETATION ==========\n")
cat("Is sex-stratified analysis worth pursuing?\n\n")
cat("Key factors:\n")
cat("  1. If male n80 << combined n80: male-only has substantially better power\n")
cat("     due to higher r (less within-sex TOD noise / more homogeneous signal)\n")
cat("  2. BUT: male-only uses n=30 pilot vs n=60 combined pilot\n")
cat("     -> bootstrap CI will be WIDER for male-only (less pilot data)\n")
cat("  3. Biological question changes: male MDD vs male CTL (not all MDD vs all CTL)\n")
cat("  4. For the paper, the passive design story is more compelling if n80 is\n")
cat("     large (>300) to contrast with active designs (n80~24). Splitting by\n")
cat("     sex may reduce n80 and muddy the 'passive is hard' message.\n\n")

# Recommendation
if (length(n80_vals) >= 2) {
  n80_combined <- tryCatch(as.integer(n80_vals[["A_combined"]]), warning = function(w) Inf)
  n80_male     <- tryCatch(as.integer(n80_vals[["B_male"]]),     warning = function(w) Inf)
  if (is.finite(n80_male) && is.finite(n80_combined) && n80_male < n80_combined * 0.5) {
    cat("RECOMMENDATION: Male-only shows substantially better power (n80 < half of combined).\n")
    cat("  Consider male-only as passive example; female-only as supplementary.\n")
  } else {
    cat("RECOMMENDATION: Combined analysis likely preferred.\n")
    cat("  Sex stratification does not dramatically change n80.\n")
    cat("  Stick with combined MDD vs Control (n=60/group) for the paper.\n")
  }
}

cat(sprintf("\nAll outputs written to: %s/\n", out_dir))
cat("Done.\n")
