#' QA the SCP pilot database before merge to main
#'
#' Three-pass check on every pilot rds:
#'   1. Schema:    has the required CircadianBioOptions slots
#'   2. Values:    A > 0, sigma > 0, phase in [0, 24], r-tilde median finite
#'   3. Drive-test: a tiny runSimsSingleCohort sim returns values in [0, 1]
#'
#' Pilots that fail any step get marked status = "qa_failed_*" in the
#' manifest and excluded from summary/tissue_signal_summary.csv.
#'
#' Usage:
#'   Rscript scripts/qa_pilot_database.R

setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
suppressPackageStartupMessages(library(SCP))
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0L) a else b

manifest_path <- "inst/extdata/pilots/manifest.csv"
m <- read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
ing <- which(m$status == "ingested")
cat(sprintf("Auditing %d ingested pilots...\n\n", length(ing)))

# Required CircadianBioOptions slots for simulation
required <- c("ngenes", "amplitude", "sigma_rhythmic", "phase", "period")

# Tiny simulation settings (~1-2 sec per pilot)
sim_design   <- SCP::CircadianDesignOptions(
                   sample_sizes = c(40, 80), nsims = 5L,
                   cts = seq(0, 22, by = 4))
sim_analysis <- SCP::CircadianAnalysisOptions(alpha = 0.05)

results <- vector("list", length(ing))
names(results) <- m$file[ing]

for (k in seq_along(ing)) {
  i <- ing[k]
  fpath <- file.path("inst/extdata/pilots", m$file[i])
  rec <- list(file = m$file[i], status = "qa_ok", note = "")
  res <- tryCatch({
    if (!file.exists(fpath)) stop("file missing on disk")
    p <- readRDS(fpath)
    # Step 1: schema
    miss <- setdiff(required, names(p))
    if (length(miss) > 0) stop(sprintf("schema_invalid: missing %s",
                                       paste(miss, collapse = ",")))
    # Step 2: values
    amp <- p$amplitude; sig <- p$sigma_rhythmic; phi <- p$phase
    if (!any(is.finite(amp) & amp > 0))  stop("values_invalid: no positive A")
    if (!any(is.finite(sig) & sig > 0))  stop("values_invalid: no positive sigma")
    if (!any(is.finite(phi) & phi >= 0 & phi <= 24))
      stop("values_invalid: phase out of [0,24]")
    r <- amp[is.finite(amp) & amp > 0] / sig[is.finite(sig) & sig > 0]
    if (!is.finite(median(r, na.rm = TRUE)))
      stop("values_invalid: r-tilde median not finite")
    # Step 3: drive-test sim
    sim <- SCP::runSimsSingleCohort(
      bio.opts = p, design.opts = sim_design,
      analysis.opts = sim_analysis, mc.cores = 1L
    )
    mp <- sim$marginal_power
    if (is.null(mp) || all(is.na(mp))) stop("sim_failed: marginal_power all-NA")
    if (any(mp < 0 | mp > 1, na.rm = TRUE))
      stop("sim_failed: marginal_power out of [0,1]")
    "qa_ok"
  }, error = function(e) {
    msg <- conditionMessage(e)
    new_status <- if (grepl("schema_invalid", msg)) "qa_schema_invalid"
                  else if (grepl("values_invalid", msg)) "qa_values_invalid"
                  else if (grepl("sim_failed", msg)) "qa_sim_failed"
                  else "qa_failed"
    rec$status <<- new_status
    rec$note   <<- msg
    new_status
  })
  results[[k]] <- rec
  if ((k %% 20) == 0L)
    cat(sprintf("  %d / %d (status = %s)\n", k, length(ing), res))
}

# ---- write QA report ------------------------------------------------------
qa <- do.call(rbind, lapply(results, function(r) data.frame(r, stringsAsFactors = FALSE)))
write.csv(qa, "inst/extdata/pilots/qa_report.csv", row.names = FALSE)
cat("\n=== QA summary ===\n")
print(table(qa$status))

# ---- update manifest: demote failed pilots --------------------------------
demoted <- 0L
for (k in seq_along(ing)) {
  if (results[[k]]$status != "qa_ok") {
    m$status[ing[k]] <- results[[k]]$status
    demoted <- demoted + 1L
  }
}
write.csv(m, manifest_path, row.names = FALSE)
cat(sprintf("\nDemoted %d ingested pilots that failed QA.\n", demoted))

# ---- regenerate tissue_signal_summary.csv ---------------------------------
m_ok <- m[m$status == "ingested", , drop = FALSE]
cat(sprintf("Regenerating summary CSV for %d QA-passed pilots...\n",
            nrow(m_ok)))

safe_phases <- function(uniq_phases) {
  n <- length(uniq_phases)
  if (n == 0L || any(is.na(uniq_phases))) return(NA_character_)
  if (n <= 12) paste(formatC(uniq_phases, format = "g", digits = 3), collapse = ",")
  else sprintf("%s..%s (%d phases)",
               formatC(uniq_phases[1], format = "g", digits = 3),
               formatC(uniq_phases[n], format = "g", digits = 3), n)
}

rows <- vector("list", nrow(m_ok))
for (i in seq_len(nrow(m_ok))) {
  fpath <- file.path("inst/extdata/pilots", m_ok$file[i])
  p <- tryCatch(readRDS(fpath), error = function(e) NULL); if (is.null(p)) next
  amp <- p$amplitude %||% NA_real_
  sig <- p$sigma_rhythmic %||% NA_real_
  r   <- if (length(amp) && length(sig) && length(amp) == length(sig)) amp / sig else NA_real_
  r   <- r[is.finite(r) & r > 0]
  times <- p$times %||% p$raw$times %||% numeric(0)
  has_times <- length(times) > 0 && all(is.finite(times))
  tmod <- if (has_times) times %% 24 else numeric(0)
  uniq_phases <- if (length(tmod) > 0) sort(unique(round(tmod, 2))) else numeric(0)
  n_unique <- length(uniq_phases)
  span <- if (has_times) max(times) - min(times) else NA_real_
  step <- if (n_unique > 1) {
    dd <- diff(uniq_phases); if (max(dd) - min(dd) < 0.5) round(median(dd), 1) else NA_real_
  } else NA_real_
  rows[[i]] <- data.frame(
    species = m_ok$species[i], dataset = m_ok$dataset[i], tissue = m_ok$tissue[i],
    design = m_ok$design[i], condition = m_ok$condition[i],
    n = m_ok$n[i], ngenes = m_ok$ngenes[i],
    tod_min      = if (has_times) round(min(tmod), 1) else NA_real_,
    tod_max      = if (has_times) round(max(tmod), 1) else NA_real_,
    tod_sd       = if (has_times) round(sd(tmod), 2)  else NA_real_,
    tod_n_unique = if (has_times) n_unique else NA_integer_,
    tod_cycles   = if (has_times) round(span / 24, 2) else NA_real_,
    tod_step_hr  = step,
    tod_phases   = safe_phases(uniq_phases),
    r_median_top300 = if (length(r) > 0) round(median(r), 3) else NA_real_,
    r_q25_top300    = if (length(r) > 0) round(quantile(r, 0.25, names = FALSE), 3) else NA_real_,
    r_q75_top300    = if (length(r) > 0) round(quantile(r, 0.75, names = FALSE), 3) else NA_real_,
    prop_rhythmic   = if (!is.null(p$prop_rhythmic)) round(p$prop_rhythmic, 4) else NA_real_,
    stringsAsFactors = FALSE
  )
}
csv <- do.call(rbind, rows[!sapply(rows, is.null)])
write.csv(csv, "summary/tissue_signal_summary.csv", row.names = FALSE)
cat(sprintf("Wrote %d rows to summary/tissue_signal_summary.csv\n", nrow(csv)))
cat("\nBy species:\n"); print(table(csv$species))
cat("\nUnique tissues:\n"); print(length(unique(csv$tissue)))
