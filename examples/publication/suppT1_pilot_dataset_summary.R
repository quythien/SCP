# Current catalogue: single-harmonic circadian signal in every bundled pilot.
# Run from the SCP repository root:
#   Rscript examples/publication/suppT1_pilot_dataset_summary.R
# Optional SCP_SIGNAL_SUMMARY_OUT selects a separate review copy.
#
# r = amplitude / residual sigma uses the paired genes returned by
# scp_load_pilot(alpha_pilot = 0.01, adjust = "none", K = 1,
#                paired_sigma = TRUE): at most 300 genes, ordered by raw p.
# prop_rhythmic counts every gene passing that screen over all tested genes.
# It is distinct from the BH discovery counts in gtex_v10_signal_summary.csv.
# Quantiles use R's default type 7 and retain full numerical precision.

# 1. Use the current source loader and repository pilots -------------------
DIR <- normalizePath("inst/extdata/pilots", mustWork = TRUE)
OUT <- Sys.getenv("SCP_SIGNAL_SUMMARY_OUT", "summary/tissue_signal_summary.csv")

# A private environment keeps this script independent of an installed SCP
# version. Only the pilot-directory resolver is bound to this source tree.
pilot_api <- new.env(parent = globalenv())
pilot_api$`%||%` <- function(a, b) if (!is.null(a)) a else b
for (file in c("R/estimation.R", "R/simulation.R", "R/pilot_database.R"))
  sys.source(file, envir = pilot_api)
pilot_api$.scp_pilots_dir <- function() DIR

m <- pilot_api$scp_pilots(status = "ingested")
key_columns <- c("species", "dataset", "tissue", "condition", "time_basis")
stopifnot(nrow(m) > 0L, !anyDuplicated(m[key_columns]),
          all(file.exists(file.path(DIR, m$file))))

# 2. Summarize observed native sampling times, when available --------------
# GTEx contains an aggregate time distribution. Its grid does not represent
# observed samples, so the empirical time-summary columns remain NA.
tod_summary <- function(times) {
  if (!length(times))
    return(list(tod_min = NA_real_, tod_max = NA_real_, tod_sd = NA_real_,
                tod_n_unique = NA_integer_, tod_cycles = NA_real_,
                tod_step_hr = NA_character_, tod_phases = NA_character_))
  stopifnot(is.numeric(times), all(is.finite(times)))
  t_mod <- times %% 24
  uniq <- sort(unique(round(t_mod, 2)))
  nu <- length(uniq)
  span <- max(times) - min(times)
  d <- diff(uniq)
  regular <- length(d) > 0L && max(d) - min(d) < 0.5
  step <- if (regular) as.character(round(stats::median(d), 1)) else "irregular"
  phases <- if (nu <= 12L)
    paste(formatC(uniq, format = "g", digits = 3), collapse = ",")
  else
    paste0(formatC(uniq[1], format = "g", digits = 3), "..",
           formatC(uniq[nu], format = "g", digits = 3), " (", nu, " phases)")
  list(tod_min = round(min(t_mod), 1), tod_max = round(max(t_mod), 1),
       tod_sd = round(stats::sd(t_mod), 2), tod_n_unique = nu,
       tod_cycles = round(span / 24, 2), tod_step_hr = step,
       tod_phases = phases)
}

# 3. Load each pilot and verify the selected amplitude/sigma pairs ---------
rows <- vector("list", nrow(m))
for (i in seq_len(nrow(m))) {
  r <- m[i, , drop = FALSE]
  p <- pilot_api$scp_load_pilot(r$species, r$dataset, r$tissue, r$condition,
                               alpha_pilot = 0.01, adjust = "none", K = 1,
                               paired_sigma = TRUE)
  rf <- p$rhythm_fit
  stopifnot(is.data.frame(rf), nrow(rf) > 0L,
            all(c("pvalue", "A", "sigma") %in% names(rf)),
            all(is.finite(rf$pvalue)), all(diff(rf$pvalue) >= 0),
            isTRUE(p$pilot_cap >= 0.01), isTRUE(p$pilot_top_k == 300L),
            isTRUE(p$ngenes == r$ngenes),
            isTRUE(p$rhythm_denom == r$ngenes),
            is.finite(r$ngenes), r$ngenes > 0L)
  candidates <- rf[rf$pvalue < 0.01, , drop = FALSE]
  n_rhythmic <- nrow(candidates)
  n_selected <- min(300L, n_rhythmic)
  selected <- candidates[seq_len(n_selected), , drop = FALSE]
  stopifnot(length(p$amplitude) == n_selected,
            length(p$sigma_rhythmic) == n_selected,
            identical(p$amplitude, selected$A),
            identical(p$sigma_rhythmic, selected$sigma),
            all(is.finite(p$amplitude) & p$amplitude > 0),
            all(is.finite(p$sigma_rhythmic) & p$sigma_rhythmic > 0),
            isTRUE(all.equal(p$prop_rhythmic, n_rhythmic / r$ngenes)),
            n_rhythmic <= r$ngenes)
  rv <- p$amplitude / p$sigma_rhythmic
  stopifnot(all(is.finite(rv) & rv > 0))
  qs <- if (n_selected > 0L)
    stats::quantile(rv, c(0.25, 0.5, 0.75), names = FALSE, type = 7)
  else rep(NA_real_, 3L)

  # Check sample counts against both the manifest and pilot metadata.
  stopifnot(is.finite(r$n), r$n >= 2L, r$n == as.integer(r$n))
  if (identical(r$dataset, "GTEx")) {
    stopifnot(!length(p$times), !length(p$raw), !length(p$cts),
              isTRUE(p$n_samples == r$n),
              isTRUE(p$time_distribution$n_observed == r$n),
              identical(p$gtex_release, "v10"), identical(p$time_basis, "TOD"))
    times <- NULL
    time_summary_source <- "aggregate_only"
  } else {
    times <- p$times
    if (is.null(times)) times <- p$raw$times
    if (length(times)) {
      stopifnot(length(times) == r$n)
      time_summary_source <- "observed_times"
    } else {
      pilot_api$.scp_validate_time_distribution(p$time_distribution)
      stopifnot(isTRUE(p$n_samples == r$n),
                isTRUE(p$time_distribution$n_observed == r$n))
      time_summary_source <- "aggregate_only"
    }
  }
  td <- tod_summary(times)
  rows[[i]] <- data.frame(
    species = r$species, dataset = r$dataset, tissue = r$tissue,
    tissue_canonical = r$tissue_canonical, design = r$design,
    condition = r$condition, n = r$n, ngenes = r$ngenes,
    tod_min = td$tod_min, tod_max = td$tod_max, tod_sd = td$tod_sd,
    tod_n_unique = td$tod_n_unique, tod_cycles = td$tod_cycles,
    tod_step_hr = td$tod_step_hr, tod_phases = td$tod_phases,
    r_median_top300 = qs[2], r_q25_top300 = qs[1], r_q75_top300 = qs[3],
    prop_rhythmic = p$prop_rhythmic,
    n_selected = n_selected, n_rhythmic = n_rhythmic,
    alpha_pilot = 0.01, adjust_pilot = "none", K = 1L, paired_sigma = TRUE,
    gtex_release = r$gtex_release, time_basis = r$time_basis,
    time_summary_source = time_summary_source, stringsAsFactors = FALSE)
}

# 4. Save one row per current ingested pilot -------------------------------
result <- do.call(rbind, rows)
stopifnot(nrow(result) == nrow(m), !anyDuplicated(result[key_columns]))
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(result, OUT, row.names = FALSE)
cat(sprintf("Saved %s: %d pilots (%d GTEx v10 TOD, %d native).\n", OUT,
            nrow(result), sum(result$dataset == "GTEx"),
            sum(result$dataset != "GTEx")))
