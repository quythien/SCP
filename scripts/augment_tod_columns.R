#' Augment tissue_signal_summary.csv with TOD-pattern columns
#'
#' Adds four columns describing the sampling-time-of-day pattern of each pilot:
#'   tod_n_unique : number of distinct sampling phases (mod 24)
#'   tod_cycles   : total time span / 24, in cycles
#'   tod_step_hr  : common inter-phase step in hours if the grid is regular, else NA
#'   tod_phases   : comma-separated sorted unique phases (mod 24) when <= 12 distinct
#'                  phases; otherwise a range summary
#'
#' Strategy: walk each available source data file, extract the time-of-day
#' vector for the matching (species, dataset, tissue, condition) tuple, and
#' merge into the existing CSV. Rows without a locally-available source keep
#' NA; future ingestions populate these via the generator at
#' examples/publication/suppT1_pilot_dataset_summary.R (already patched).
#'
#' Usage:
#'   Rscript scripts/augment_tod_columns.R

`%||%` <- function(a, b) if (!is.null(a)) a else b
setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")

# ---------- helper -----------------------------------------------------------
tod_summary <- function(times) {
  if (length(times) == 0 || all(is.na(times))) {
    return(list(tod_n_unique = NA_integer_, tod_cycles = NA_real_,
                tod_step_hr  = NA_real_,    tod_phases  = NA_character_))
  }
  t_mod      <- times %% 24
  uniq_phases <- sort(unique(round(t_mod, 2)))
  n_unique    <- length(uniq_phases)
  span_full   <- max(times) - min(times)
  tod_cycles  <- if (span_full > 0) round(span_full / 24, 2) else 0
  step_diffs  <- if (n_unique >= 2) diff(uniq_phases) else numeric(0)
  is_regular  <- length(step_diffs) > 0 && (max(step_diffs) - min(step_diffs)) < 0.5
  tod_step    <- if (is_regular) round(median(step_diffs), 1) else NA_real_
  tod_phases  <- if (n_unique <= 12)
    paste(formatC(uniq_phases, format = "g", digits = 3), collapse = ",")
  else
    sprintf("%s..%s (%d phases)",
            formatC(uniq_phases[1],        format = "g", digits = 3),
            formatC(uniq_phases[n_unique], format = "g", digits = 3),
            n_unique)
  list(tod_n_unique = n_unique, tod_cycles = tod_cycles,
       tod_step_hr  = tod_step, tod_phases  = tod_phases)
}

# ---------- collect times by key ---------------------------------------------
times_by_key <- list()
key <- function(sp, ds, ts, cd) paste(sp, ds, ts, cd, sep = "|")

# Mouse GSE54651 (9 valid tissues, condition "All")
if (file.exists("data/mice_GSE54651_CPM.RData")) {
  d <- tryCatch(readRDS("data/mice_GSE54651_CPM.RData"), error = function(e) NULL)
  if (!is.null(d)) {
    for (tn in names(d$count_clean)) {
      tt <- d$tod[[tn]]
      if (length(tt) > 0)
        times_by_key[[key("Mouse", "GSE54651", tn, "All")]] <- tt
    }
  }
}

# GTEx single-tissue pilots locally available
gtex_local <- list(Adrenal = "data/gtex_AdrenalGland_single_pilot.rds",
                   Liver   = "data/gtex_Liver_single_pilot.rds")
for (nm in names(gtex_local)) {
  if (file.exists(gtex_local[[nm]])) {
    p <- tryCatch(readRDS(gtex_local[[nm]]), error = function(e) NULL)
    tt <- p$raw$times %||% p$times %||% NULL
    if (!is.null(tt))
      times_by_key[[key("Human", "GTEx", nm, "All")]] <- tt
  }
}

# GSE160521 striatum control pilots (3 regions x 1 cond = "CONTROL")
gse160521_paths <- list(
  list(rds="data/gse160521_nac_ctrl_pilot.rds",     tissue="NAc"),
  list(rds="data/gse160521_caudate_ctrl_pilot.rds", tissue="Caudate"),
  list(rds="data/gse160521_putamen_ctrl_pilot.rds", tissue="Putamen")
)
for (e in gse160521_paths) {
  if (file.exists(e$rds)) {
    p <- tryCatch(readRDS(e$rds), error = function(e) NULL)
    tt <- p$raw$times %||% p$times %||% NULL
    if (!is.null(tt))
      times_by_key[[key("Human", "GSE160521", e$tissue, "CONTROL")]] <- tt
  }
}

# Seney-ACC
if (file.exists("data/seney_ctrl_pilot.rds")) {
  p <- tryCatch(readRDS("data/seney_ctrl_pilot.rds"), error = function(e) NULL)
  tt <- p$raw$times %||% p$times %||% NULL
  if (!is.null(tt))
    times_by_key[[key("Human", "Seney-ACC", "ACC", "Control")]] <- tt
}

# BA11-BA47 (younger only locally; older needs CMC raw)
if (file.exists("data/ba11_ba47_younger.rds")) {
  p <- tryCatch(readRDS("data/ba11_ba47_younger.rds"), error = function(e) NULL)
  for (region in c("BA11", "BA47")) {
    tt <- p[[region]]$raw$times %||% p[[region]]$times %||% p$raw$times %||% NULL
    if (!is.null(tt))
      times_by_key[[key("Human", "BA11-BA47", region, "younger")]] <- tt
  }
}

# Baboon GSE98965 — Mure 2018, 61 tissues from CAMO_PRC_hmb.RData
if (file.exists("data/CAMO_PRC_hmb.RData")) {
  env <- new.env()
  ok <- tryCatch({load("data/CAMO_PRC_hmb.RData", envir = env); TRUE},
                 error = function(e) FALSE)
  if (ok) {
    objs <- ls(env)
    bb <- NULL
    for (o in objs) {
      v <- get(o, envir = env)
      if (is.list(v) && !is.null(v$tod) && is.list(v$tod)) { bb <- v; break }
    }
    if (!is.null(bb)) {
      for (tn in names(bb$tod))
        if (length(bb$tod[[tn]]) > 0)
          times_by_key[[key("Baboon", "GSE98965", tn, "Control")]] <- bb$tod[[tn]]
    }
  }
}

# ---------- generic walk of inst/extdata/pilots/<species>/*.rds --------------
# For every saved pilot rds, harvest p$raw$times (or p$times) and key by the
# manifest tuple (species, dataset, tissue, condition).  This makes the script
# generic to any future ingestion.
manifest_path <- "inst/extdata/pilots/manifest.csv"
if (file.exists(manifest_path)) {
  m <- read.csv(manifest_path, stringsAsFactors = FALSE)
  pilots_root <- "inst/extdata/pilots"
  for (i in seq_len(nrow(m))) {
    fpath <- file.path(pilots_root, m$file[i])
    if (!file.exists(fpath)) next
    p <- tryCatch(readRDS(fpath), error = function(e) NULL)
    if (is.null(p)) next
    tt <- p$raw$times %||% p$times %||% NULL
    if (is.null(tt) || length(tt) == 0) next
    sp <- m$species[i]; ds <- m$dataset[i]; ts <- m$tissue[i]; cd <- m$condition[i]
    # Title-case species for the CSV key compatibility
    sp_disp <- paste0(toupper(substr(sp,1,1)), substr(sp,2,nchar(sp)))
    times_by_key[[key(sp_disp, ds, ts, cd)]] <- as.numeric(tt)
  }
}

# ---------- apply ------------------------------------------------------------
csv_path <- "output/supp_tissue_summary/tissue_signal_summary.csv"
csv <- read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)
new_cols <- c("tod_n_unique", "tod_cycles", "tod_step_hr", "tod_phases")
for (c in new_cols) csv[[c]] <- if (c == "tod_phases") NA_character_ else NA

matched <- 0L
for (i in seq_len(nrow(csv))) {
  k  <- key(csv$species[i], csv$dataset[i], csv$tissue[i], csv$condition[i])
  tt <- times_by_key[[k]]
  if (!is.null(tt)) {
    s <- tod_summary(tt)
    for (c in new_cols) csv[[c]][i] <- s[[c]]
    matched <- matched + 1L
  }
}

# place new columns right after tod_sd
ord     <- names(csv)
idx_sd  <- which(ord == "tod_sd")
ord_new <- c(ord[1:idx_sd], new_cols,
             setdiff(ord[(idx_sd + 1):length(ord)], new_cols))
csv <- csv[, ord_new]

write.csv(csv, csv_path, row.names = FALSE)
file.copy(csv_path, "submission/tissue_signal_summary.csv", overwrite = TRUE)

cat(sprintf("Augmented %d / %d rows with TOD pattern columns.\n",
            matched, nrow(csv)))
cat("\nBreakdown by dataset (has_tod_pattern):\n")
print(table(csv$dataset, !is.na(csv$tod_n_unique),
            dnn = c("dataset", "has_tod_pattern")))
cat("\nSample rows with TOD pattern:\n")
have <- which(!is.na(csv$tod_n_unique))[c(1, 3, 5, 7)]
have <- have[!is.na(have)]
print(csv[have, c("species","dataset","tissue","condition","n",
                  "tod_n_unique","tod_cycles","tod_step_hr","tod_phases")])
