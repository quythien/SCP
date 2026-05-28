# Append newly-ingested pilots (from inst/extdata/pilots/manifest.csv) to
# output/supp_tissue_summary/tissue_signal_summary.csv using the same row
# schema as examples/publication/suppT1_pilot_dataset_summary.R::make_row.

setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
`%||%` <- function(a, b) if (!is.null(a)) a else b

make_row <- function(p, times, species, dataset, tissue, design, condition) {
  if (is.null(p)) return(NULL)
  pvals   <- p$raw$pvalue
  qvals   <- p.adjust(pvals, "BH")
  r_all   <- p$raw$r
  phi_all <- p$raw$phi

  ok  <- !is.na(r_all) & !is.na(pvals) & is.finite(r_all) & r_all > 0
  idx <- order(pvals[ok])[seq_len(min(300L, sum(ok)))]
  rv  <- r_all[ok][idx]; rv <- rv[is.finite(rv)]
  rhy_ok  <- !is.na(pvals) & pvals < 0.05 & !is.na(phi_all)
  phi_rhy <- phi_all[rhy_ok] %% 24

  t_mod <- times %% 24
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
    species         = species, dataset = dataset, tissue = tissue,
    design = design, condition = condition,
    n               = length(times), ngenes = length(pvals),
    tod_min         = round(min(t_mod), 1),
    tod_max         = round(max(t_mod), 1),
    tod_sd          = round(sd(t_mod),  2),
    tod_n_unique    = n_unique, tod_cycles = tod_cycles,
    tod_step_hr     = tod_step, tod_phases = tod_phases,
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

csv_path <- "output/supp_tissue_summary/tissue_signal_summary.csv"
csv <- read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)
existing_key <- with(csv, paste(species, dataset, tissue, condition, sep = "|"))

m <- read.csv("inst/extdata/pilots/manifest.csv", stringsAsFactors = FALSE)
ingested <- m[m$status == "ingested", ]

new_rows <- list()
n_skip <- 0L; n_norts <- 0L
for (i in seq_len(nrow(ingested))) {
  r <- ingested[i, ]
  sp_disp <- paste0(toupper(substr(r$species,1,1)), substr(r$species,2,nchar(r$species)))
  key_i <- paste(sp_disp, r$dataset, r$tissue, r$condition, sep = "|")
  if (key_i %in% existing_key) { n_skip <- n_skip + 1L; next }
  fpath <- file.path("inst/extdata/pilots", r$file)
  if (!file.exists(fpath)) next
  p <- tryCatch(readRDS(fpath), error = function(e) NULL)
  if (is.null(p)) next
  tt <- p$raw$times %||% p$times %||% NULL
  if (is.null(tt) || length(tt) == 0) { n_norts <- n_norts + 1L; next }
  row <- make_row(p, as.numeric(tt), sp_disp, r$dataset, r$tissue, r$design, r$condition)
  if (!is.null(row)) new_rows[[length(new_rows) + 1L]] <- row
}

cat(sprintf("Existing CSV rows: %d\n", nrow(csv)))
cat(sprintf("Already present (skipped): %d\n", n_skip))
cat(sprintf("New ingested w/o times: %d\n", n_norts))
cat(sprintf("New rows to append: %d\n", length(new_rows)))

if (length(new_rows) > 0) {
  add <- do.call(rbind, new_rows)
  # Align columns to existing CSV (pad missing as NA)
  for (c in setdiff(colnames(csv), colnames(add))) add[[c]] <- NA
  for (c in setdiff(colnames(add), colnames(csv))) csv[[c]] <- NA
  add <- add[, colnames(csv)]
  out <- rbind(csv, add)
  write.csv(out, csv_path, row.names = FALSE)
  dir.create("submission", showWarnings = FALSE)
  file.copy(csv_path, "submission/tissue_signal_summary.csv", overwrite = TRUE)
  cat(sprintf("Wrote %d rows total to %s\n", nrow(out), csv_path))
} else {
  cat("No new rows to append.\n")
}
