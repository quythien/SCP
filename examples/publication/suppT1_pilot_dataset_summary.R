#' =======================================================================
#' Supplementary Table - Single-Tissue Circadian Signal Summary
#' =======================================================================
#'
#' Summarises the circadian signal strength of every pilot the package
#' actually bundles in inst/extdata/pilots/, so the output stays consistent
#' with what ships (paired top-K r-tilde = A/sigma, current manifest). One row
#' per ingested manifest pilot whose .rds exists on disk.
#'
#' This reads the BUNDLED CircadianBioOptions pilots (which carry the per-gene
#' rhythm_fit table and the paired amplitude/sigma/phase vectors); it does NOT
#' re-fit from raw expression. To rebuild the underlying pilots, see
#' pilot_builders/build_k1_pilots.R, build_2h_pilots.R, upgrade_geo_pilots.R.
#'
#' COLUMNS:
#'   species, dataset, tissue, design, condition, n, ngenes,
#'   tod_min, tod_max, tod_sd, tod_n_unique, tod_cycles, tod_step_hr,
#'   tod_phases, r_median_top300, r_q25_top300, r_q75_top300, prop_rhythmic
#'
#' USAGE:
#'   Rscript examples/publication/suppT1_pilot_dataset_summary.R

`%||%` <- function(a, b) if (!is.null(a)) a else b

DIR <- "inst/extdata/pilots"
OUT <- "submission/tissue_signal_summary.csv"

m <- utils::read.csv(file.path(DIR, "manifest.csv"), stringsAsFactors = FALSE)
# Only summarise pilots the package actually exposes: status == "ingested"
# (same filter as scp_pilots()/the Shiny app) and the .rds present on disk.
# This excludes ingest_incomplete / failed arms (e.g. tiny-n microarray arms
# whose residual sigma collapses and inflates r-tilde).
m <- m[m$status == "ingested" & file.exists(file.path(DIR, m$file)), , drop = FALSE]

# Time-of-day design summary from a pilot's empirical sampling times (mod 24).
tod_summary <- function(times) {
  if (is.null(times) || length(times) < 1L || all(!is.finite(times)))
    return(list(tod_min = NA, tod_max = NA, tod_sd = NA, tod_n_unique = NA,
                tod_cycles = NA, tod_step_hr = NA, tod_phases = NA_character_))
  times <- times[is.finite(times)]
  t_mod <- times %% 24
  uniq  <- sort(unique(round(t_mod, 2)))
  nu    <- length(uniq)
  span  <- max(times) - min(times)
  cyc   <- if (span > 0) round(span / 24, 2) else 0
  d     <- if (nu >= 2) diff(uniq) else numeric(0)
  reg   <- length(d) > 0 && (max(d) - min(d) < 0.5)
  step  <- if (reg) round(stats::median(d), 1) else NA_real_
  phases <- if (nu <= 12)
    paste(formatC(uniq, format = "g", digits = 3), collapse = ",")
  else
    paste0(formatC(uniq[1], format = "g", digits = 3), "..",
           formatC(uniq[nu], format = "g", digits = 3), " (", nu, " phases)")
  list(tod_min = round(min(t_mod), 1), tod_max = round(max(t_mod), 1),
       tod_sd = round(sd(t_mod), 2), tod_n_unique = nu, tod_cycles = cyc,
       tod_step_hr = step, tod_phases = phases)
}

rows <- lapply(seq_len(nrow(m)), function(i) {
  r <- m[i, ]
  p <- tryCatch(readRDS(file.path(DIR, r$file)), error = function(e) NULL)
  if (is.null(p)) return(NULL)
  # Paired top-K effect size r-tilde = A / sigma (gene-aligned as bundled).
  amp <- p$amplitude %||% numeric(0); sig <- p$sigma_rhythmic %||% numeric(0)
  k   <- min(length(amp), length(sig))
  rv  <- if (k > 0) { v <- amp[seq_len(k)] / sig[seq_len(k)]; v[is.finite(v) & v > 0] } else numeric(0)
  times <- p$times %||% p$raw$times %||% NULL
  td  <- tod_summary(times)
  data.frame(
    species = r$species, dataset = r$dataset, tissue = r$tissue,
    design = r$design, condition = r$condition,
    n = if (!is.null(times)) length(times) else r$n,
    ngenes = p$ngenes %||% r$ngenes,
    tod_min = td$tod_min, tod_max = td$tod_max, tod_sd = td$tod_sd,
    tod_n_unique = td$tod_n_unique, tod_cycles = td$tod_cycles,
    tod_step_hr = td$tod_step_hr, tod_phases = td$tod_phases,
    r_median_top300 = round(if (length(rv) >= 3) median(rv) else NA_real_, 3),
    r_q25_top300    = round(if (length(rv) >= 3) quantile(rv, 0.25, names = FALSE) else NA_real_, 3),
    r_q75_top300    = round(if (length(rv) >= 3) quantile(rv, 0.75, names = FALSE) else NA_real_, 3),
    prop_rhythmic   = round(p$prop_rhythmic %||% NA_real_, 4),
    stringsAsFactors = FALSE)
})

result <- do.call(rbind, Filter(Negate(is.null), rows))
utils::write.csv(result, OUT, row.names = FALSE)
cat(sprintf("Saved %s: %d rows (Seney: %d)\n", OUT, nrow(result),
            sum(grepl("Seney", result$dataset, ignore.case = TRUE))))
cat("r_median_top300 range:",
    paste(round(range(result$r_median_top300, na.rm = TRUE), 2), collapse = " - "), "\n")
