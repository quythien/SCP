#' =======================================================================
#' Supplementary Table — Single-Tissue Circadian Signal Summary
#' =======================================================================
#'
#' DATASETS INCLUDED:
#'   1. GTEx (Human, passive)      — all tissues with n >= 80
#'   2. Baboon GSE98965 (active)   — all 61 tissues
#'   3. Mouse D1D2 (active)        — D1 and D2 striatal cell types
#'   4. Mouse GSE54651 (active)    — all 12 tissues
#'   5. Seney ACC (Human, passive) — controls only
#'   6. BA11/BA47 (Human, passive) — younger + older, split by region
#'   7. GSE160521 (Human, passive) — NAc, Caudate, Putamen, all 3 conditions
#'
#' COLUMNS (left to right):
#'   species, dataset, tissue, design, condition, n, ngenes,
#'   tod_min, tod_max, tod_sd,
#'   r_median_top300, r_q25_top300, r_q75_top300,
#'   phi_median_rhy, phi_q25_rhy, phi_q75_rhy,
#'   rhy_FDR20, rhy_FDR10, rhy_FDR05, rhy_FDR01,
#'   rhy_p05, rhy_p01, rhy_p001
#'
#' Row order: Human -> Mouse -> Baboon, then dataset -> tissue -> condition
#'
#' USAGE:
#'   Rscript examples/publication/13_supp_tissue_summary.R
#'   SMOKE_TEST=true Rscript examples/publication/13_supp_tissue_summary.R
#'
#' @author Thien Quy Pham

# =====================================================================
# 0. Setup
# =====================================================================
SMOKE_TEST <- identical(Sys.getenv("SMOKE_TEST"), "true")
N_CORES    <- as.integer(Sys.getenv("MC_CORES", unset = "4"))

old_wd <- setwd("code")
source("setup.R")
setwd(old_wd)

out_dir   <- "output/supp_tissue_summary"
cache_dir <- file.path(out_dir, "cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
cat(sprintf("\n=== Supplementary Tissue Summary [%s] ===\n", timestamp))
cat(sprintf("Mode     : %s\n", if (SMOKE_TEST) "SMOKE" else "PRODUCTION"))
cat(sprintf("mc.cores : %d\n\n", N_CORES))

# =====================================================================
# Helper: fit or load cached estimate_circadian_params
# input_type: "log2" (no preprocessing) or "cpm" (prepCircadianData first)
# =====================================================================
run_or_cache <- function(tag, expr, times, input_type = "log2") {
  rds <- file.path(cache_dir, sprintf("%s.rds", tag))
  if (file.exists(rds)) {
    cat(sprintf("  [cache] %s\n", tag))
    return(readRDS(rds))
  }
  cat(sprintf("  [fit]   %s  (%d genes x %d samples)\n", tag, nrow(expr), ncol(expr)))
  if (input_type == "cpm") {
    prep  <- prepCircadianData(expr, times = times, input_type = "cpm")
    expr  <- prep$data
    times <- prep$times
  }
  p <- tryCatch(
    estimate_circadian_params(expr, times, verbose = FALSE),
    error = function(e) { message(sprintf("    ERROR: %s", e$message)); NULL }
  )
  if (!is.null(p)) saveRDS(p, rds)
  p
}

# =====================================================================
# Helper: build one summary row
# =====================================================================
make_row <- function(tag, p, times, species, dataset, tissue, design, condition) {
  if (is.null(p)) return(NULL)

  pvals   <- p$raw$pvalue
  qvals   <- p.adjust(pvals, "BH")
  r_all   <- p$raw$r
  phi_all <- p$raw$phi

  # r: top-150 by p-value
  ok  <- !is.na(r_all) & !is.na(pvals) & is.finite(r_all) & r_all > 0
  idx <- order(pvals[ok])[seq_len(min(300L, sum(ok)))]
  rv  <- r_all[ok][idx]
  rv  <- rv[is.finite(rv)]

  # phase: rhythmic genes (p < 0.05)
  rhy_ok  <- !is.na(pvals) & pvals < 0.05 & !is.na(phi_all)
  phi_rhy <- phi_all[rhy_ok] %% 24

  t_mod <- times %% 24

  # TOD design summary: distinct sampling phases (mod 24), cycle count, step
  uniq_phases <- sort(unique(round(t_mod, 2)))
  n_unique    <- length(uniq_phases)
  span_full   <- max(times) - min(times)
  tod_cycles  <- if (span_full > 0) round(span_full / 24, 2) else 0
  step_diffs  <- if (n_unique >= 2) diff(uniq_phases) else numeric(0)
  is_regular  <- length(step_diffs) > 0 && max(step_diffs) - min(step_diffs) < 0.5
  tod_step    <- if (is_regular) round(stats::median(step_diffs), 1) else NA_real_
  tod_phases  <- if (n_unique <= 12)
                   paste(formatC(uniq_phases, format = "g", digits = 3), collapse = ",")
                 else
                   paste0(formatC(uniq_phases[1], format = "g", digits = 3),
                          "..", formatC(uniq_phases[n_unique], format = "g", digits = 3),
                          " (", n_unique, " phases)")

  as.data.frame(list(
    species         = species,
    dataset         = dataset,
    tissue          = tissue,
    design          = design,
    condition       = condition,
    n               = length(times),
    ngenes          = length(pvals),
    tod_min         = round(min(t_mod), 1),
    tod_max         = round(max(t_mod), 1),
    tod_sd          = round(sd(t_mod),  2),
    tod_n_unique    = n_unique,
    tod_cycles      = tod_cycles,
    tod_step_hr     = tod_step,
    tod_phases      = tod_phases,
    r_median_top300 = round(if (length(rv) >= 3) median(rv)                        else NA_real_, 3),
    r_q25_top300    = round(if (length(rv) >= 3) quantile(rv, 0.25, names = FALSE) else NA_real_, 3),
    r_q75_top300    = round(if (length(rv) >= 3) quantile(rv, 0.75, names = FALSE) else NA_real_, 3),
    phi_median_rhy  = round(if (length(phi_rhy) >= 3) median(phi_rhy)                        else NA_real_, 2),
    phi_q25_rhy     = round(if (length(phi_rhy) >= 3) quantile(phi_rhy, 0.25, names = FALSE) else NA_real_, 2),
    phi_q75_rhy     = round(if (length(phi_rhy) >= 3) quantile(phi_rhy, 0.75, names = FALSE) else NA_real_, 2),
    rhy_FDR20       = sum(qvals < 0.20, na.rm = TRUE),
    rhy_FDR10       = sum(qvals < 0.10, na.rm = TRUE),
    rhy_FDR05       = sum(qvals < 0.05, na.rm = TRUE),
    rhy_FDR01       = sum(qvals < 0.01, na.rm = TRUE),
    rhy_p05         = sum(pvals < 0.05,  na.rm = TRUE),
    rhy_p01         = sum(pvals < 0.01,  na.rm = TRUE),
    rhy_p001        = sum(pvals < 0.001, na.rm = TRUE)
  ), stringsAsFactors = FALSE)
}

rows <- list()

# =====================================================================
# 1. GTEx — all tissues with n >= 80 (Human, passive)
# =====================================================================
cat("\n--- 1. GTEx ---\n")
gtex_path <- "/home/qtp1/Projects/Collaborative/GTEXdata/TOD Xiangning/data/CPM/CPM.all.norm.RData"
load(gtex_path)

extract_gtex <- function(tissue_name) {
  df   <- CPM.all.norm[[tissue_name]]
  if (is.null(df)) return(NULL)
  ids  <- as.character(colnames(df))
  hhmm <- sapply(strsplit(ids, "\\."), function(x) if (length(x) >= 3) x[3] else NA)
  hrs  <- as.numeric(substr(hhmm, 1, 2)) + as.numeric(substr(hhmm, 3, 4)) / 60
  ok   <- !is.na(hrs)
  list(expr = as.matrix(df[, ok]), times = hrs[ok], n = sum(ok))
}

MIN_N_GTEX   <- 80L
gtex_tissues <- if (SMOKE_TEST) names(CPM.all.norm)[1:3] else names(CPM.all.norm)

# Pre-extract all valid tissues before parallel fork
gtex_data <- lapply(gtex_tissues, function(tn) {
  t <- extract_gtex(tn)
  if (is.null(t) || t$n < MIN_N_GTEX) return(NULL)
  list(tn = tn, expr = t$expr, times = t$times)
})
gtex_data <- Filter(Negate(is.null), gtex_data)
rm(CPM.all.norm)
cat(sprintf("  %d GTEx tissues to fit (n >= %d)\n", length(gtex_data), MIN_N_GTEX))

gtex_rows <- parallel::mclapply(gtex_data, function(d) {
  tag <- sprintf("gtex_%s", gsub("[^A-Za-z0-9]", "_", d$tn))
  p   <- run_or_cache(tag, d$expr, d$times)
  make_row(tag, p, d$times, "Human", "GTEx", d$tn, "passive", "All")
}, mc.cores = N_CORES)
rows <- c(rows, Filter(Negate(is.null), gtex_rows))

# =====================================================================
# 2. Baboon GSE98965 — all 61 tissues (active, ZT design)
# =====================================================================
cat("\n--- 2. Baboon GSE98965 ---\n")
load("data/CAMO_PRC_hmb.RData")

bab_tissues <- if (SMOKE_TEST) names(baboon_withTOD$baboon)[1:3] else names(baboon_withTOD$baboon)

bab_data <- lapply(bab_tissues, function(tn) {
  list(tn = tn,
       expr  = baboon_withTOD$baboon[[tn]],
       times = baboon_withTOD$tod[[tn]])
})
rm(baboon_withTOD, gtex, mice)

bab_rows <- parallel::mclapply(bab_data, function(d) {
  tag <- sprintf("baboon_%s", d$tn)
  p   <- run_or_cache(tag, d$expr, d$times, input_type = "cpm")
  make_row(tag, p, d$times, "Baboon", "GSE98965", d$tn, "active", "Control")
}, mc.cores = N_CORES)
rows <- c(rows, Filter(Negate(is.null), bab_rows))

# =====================================================================
# 3. Mouse D1D2 — D1 and D2 (active) — small, sequential
# =====================================================================
cat("\n--- 3. Mouse D1D2 ---\n")
d1d2_expr <- read.csv("data/mouse_D1D2_logCPMfiltered_counts.csv",
                      row.names = 1, check.names = FALSE)
d1d2_meta <- read.csv("data/mouse_clinicalinfo_03082021_rmOutliers.csv")

for (ct in c("D1", "D2")) {
  idx   <- d1d2_meta$cell == ct
  cols  <- as.character(d1d2_meta$sample[idx])
  cols  <- cols[cols %in% colnames(d1d2_expr)]
  times <- d1d2_meta$time[d1d2_meta$cell == ct & d1d2_meta$sample %in% cols]
  expr  <- as.matrix(d1d2_expr[, cols])
  tag   <- sprintf("mouse_D1D2_%s", ct)
  p     <- run_or_cache(tag, expr, times)
  r     <- make_row(tag, p, times, "Mouse", "D1D2-Mouse", ct, "active", "Control")
  if (!is.null(r)) rows[[length(rows)+1]] <- r
}

# =====================================================================
# 4. Mouse GSE54651 — all 12 tissues (active)
# =====================================================================
cat("\n--- 4. Mouse GSE54651 ---\n")
dat_mouse   <- readRDS("data/mice_GSE54651_CPM.RData")
gse_tissues <- if (SMOKE_TEST) names(dat_mouse$count_clean)[1:3] else names(dat_mouse$count_clean)

gse_data <- lapply(gse_tissues, function(tn) {
  list(tn = tn,
       expr  = as.matrix(dat_mouse$count_clean[[tn]]),  # count_clean is a data.frame; coerce so per-gene rows are numeric
       times = dat_mouse$tod[[tn]])
})
rm(dat_mouse)

gse_rows <- parallel::mclapply(gse_data, function(d) {
  tag <- sprintf("mouse_GSE54651_%s", d$tn)
  p   <- run_or_cache(tag, d$expr, d$times)
  make_row(tag, p, d$times, "Mouse", "GSE54651", d$tn, "active", "All")
}, mc.cores = N_CORES)
rows <- c(rows, Filter(Negate(is.null), gse_rows))

# =====================================================================
# 5. Seney ACC (Human, passive) — controls only — small, sequential
# =====================================================================
cat("\n--- 5. Seney ACC ---\n")
library(readxl)
meta_s <- read_excel("data/MD5_MetaData_1-15-25.xlsx")
tod_s  <- read_excel("data/TOD.xlsx")
expr_s <- as.matrix(read.csv("data/ACC_RNA_filtered_normalized.csv",
                              row.names = 1, check.names = FALSE))

col_ids  <- gsub("[A-Za-z]+$", "", colnames(expr_s))
meta_idx <- match(col_ids, as.character(meta_s$HU_NUM))
tod_idx  <- match(col_ids, as.character(tod_s$HU_NUM))
tod_hr   <- as.numeric(format(tod_s$OFFC_TIME[tod_idx], "%H")) +
            as.numeric(format(tod_s$OFFC_TIME[tod_idx], "%M")) / 60
ctrl_ok  <- !is.na(meta_idx) & !is.na(tod_hr) &
            meta_s$DESCRIPT[meta_idx] == "unaffected_control"

tag <- "seney_ACC_ctrl"
p   <- run_or_cache(tag, expr_s[, ctrl_ok], tod_hr[ctrl_ok])
r   <- make_row(tag, p, tod_hr[ctrl_ok], "Human", "Seney-ACC", "ACC", "passive", "Control")
if (!is.null(r)) rows[[length(rows)+1]] <- r

# =====================================================================
# 6. BA11 / BA47 (Human, passive) — younger + older — small, sequential
# =====================================================================
cat("\n--- 6. BA11/BA47 ---\n")
ba_path  <- "/home/qtp1/Projects/Collaborative/Paper/Congruence/PNAS_aging/data/combined_data.rds"
comb     <- readRDS(ba_path)
snames   <- colnames(comb$expr)
ph_order <- match(snames, comb$pheno$sample_name)
valid    <- !is.na(ph_order)
pheno    <- comb$pheno[ph_order[valid], ]
expr_ba  <- comb$expr[, valid]
pheno$tod <- if ("TOD.x" %in% names(pheno)) pheno$TOD.x else pheno$TOD.y
pheno$age <- if ("AgeGroup" %in% names(pheno)) pheno$AgeGroup else pheno$age_group

for (reg in c("BA11", "BA47")) {
  for (ag in c("younger", "older")) {
    idx <- pheno$region == reg & pheno$age == ag & !is.na(pheno$tod)
    if (sum(idx) < 20) next
    tag <- sprintf("ba11ba47_%s_%s", reg, ag)
    p   <- run_or_cache(tag, expr_ba[, idx], pheno$tod[idx])
    r   <- make_row(tag, p, pheno$tod[idx], "Human", "BA11-BA47", reg, "passive", ag)
    if (!is.null(r)) rows[[length(rows)+1]] <- r
  }
}

# =====================================================================
# 7. GSE160521 — NAc, Caudate, Putamen, all 3 conditions (Human, passive)
# =====================================================================
cat("\n--- 7. GSE160521 ---\n")
kyle_dir <- "/home/qtp1/Projects/Circadian/Kyle/Kyle_multiBrainRegion"

region_files <- list(
  NAc     = list(expr = "NAc_CPMfiltered_logCPM_1215_rm97_rm231.csv",
                 clin = "NAc_clinical_1221_rm97_rm231_matchIndex34.csv"),
  Caudate = list(expr = "Caudate_CPMfiltered_logCPM_1215_rm97_rm231.csv",
                 clin = "DS_clinical_1221_rm97_rm231_matchIndex34.csv"),
  Putamen = list(expr = "Putamen_CPMfiltered_logCPM_1215_rm97_rm231.csv",
                 clin = "DS_clinical_1221_rm97_rm231_matchIndex34.csv")
)

# Build flat list of all (region, condition) jobs
kyle_jobs <- list()
for (reg in names(region_files)) {
  rf <- region_files[[reg]]
  cl <- read.csv(file.path(kyle_dir, rf$clin))
  expr <- read.csv(file.path(kyle_dir, rf$expr), row.names = 1, check.names = FALSE)
  for (diag in unique(cl$Diagnosis.3Grp)) {
    grp  <- cl[cl$Diagnosis.3Grp == diag, ]
    cols <- as.character(grp$pair)
    cols <- cols[cols %in% colnames(expr)]
    if (length(cols) < 10) next
    times <- grp$CorrectedTOD[as.character(grp$pair) %in% cols]
    kyle_jobs[[length(kyle_jobs)+1]] <- list(
      reg = reg, diag = diag,
      expr = as.matrix(expr[, cols]), times = times
    )
  }
}

kyle_rows <- parallel::mclapply(kyle_jobs, function(d) {
  tag <- sprintf("gse160521_%s_%s", d$reg, gsub("[^A-Za-z0-9]", "_", d$diag))
  p   <- run_or_cache(tag, d$expr, d$times)
  make_row(tag, p, d$times, "Human", "GSE160521", d$reg, "passive", d$diag)
}, mc.cores = N_CORES)
rows <- c(rows, Filter(Negate(is.null), kyle_rows))

# =====================================================================
# Assemble and save
# =====================================================================
result <- do.call(rbind, Filter(Negate(is.null), rows))
species_order <- c("Human", "Mouse", "Baboon")
result$species <- factor(result$species, levels = species_order)
result <- result[order(result$species, result$dataset, result$tissue, result$condition), ]
result$species <- as.character(result$species)

out_csv <- file.path(out_dir, sprintf("tissue_signal_summary_%s.csv", timestamp))
write.csv(result, out_csv, row.names = FALSE)

cat(sprintf("\n=== Done: %d tissue entries ===\n", nrow(result)))
cat(sprintf("Saved: %s\n", out_csv))
