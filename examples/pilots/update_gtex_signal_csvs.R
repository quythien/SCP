# Update signal-summary rows for a supplied set of rebuilt GTEx pilots.
# Source this file from the rebuild script. It does not run on its own.

# 1. Identify the three existing summaries ------------------------------------

gtex_signal_csv_paths <- c(
  "summary/tissue_signal_summary.csv",
  "submission/tissue_signal_summary.csv",
  "submission/supplementary/tissue_signal_summary.csv"
)

gtex_signal_tissue_key <- function(tissue) {
  key <- tolower(gsub("[^[:alnum:]]", "", tissue))
  key[key == "adrenalgland"] <- "adrenal"
  key
}

gtex_signal_protected_keys <- c("adrenal", "liver", "muscleskeletal")


# 2. Summarize the sampling times using the existing CSV convention -----------

gtex_signal_tod_summary <- function(times) {
  if (!is.numeric(times) || !length(times) || any(!is.finite(times)))
    stop("Updated CSV rows require a finite numeric time vector.")
  if (any(times < 0 | times >= 24))
    stop("Updated sampling times must be hours in [0, 24).")

  t_mod <- times %% 24
  uniq <- sort(unique(round(t_mod, 2)))
  gaps <- diff(uniq)
  regular <- length(gaps) > 0L && diff(range(gaps)) < 0.5
  span <- diff(range(times))
  phases <- if (length(uniq) <= 12L) {
    paste(formatC(uniq, format = "g", digits = 3), collapse = ",")
  } else {
    paste0(formatC(uniq[1], format = "g", digits = 3), "..",
           formatC(tail(uniq, 1), format = "g", digits = 3),
           " (", length(uniq), " phases)")
  }

  list(
    tod_min = round(min(t_mod), 1),
    tod_max = round(max(t_mod), 1),
    tod_sd = round(stats::sd(t_mod), 2),
    tod_n_unique = length(uniq),
    tod_cycles = if (span > 0) round(span / 24, 2) else 0,
    tod_step_hr = if (regular) round(stats::median(gaps), 1) else NA_real_,
    tod_phases = phases
  )
}


# 3. Calculate signal and discovery statistics from the full gene fits -------

gtex_signal_update_statistics <- function(update) {
  required <- c("tissue", "pilot_file", "times", "ngenes", "full_fit",
                "prop_rhythmic")
  if (!all(required %in% names(update)))
    stop("An update is missing: ",
         paste(setdiff(required, names(update)), collapse = ", "))

  fit <- update$full_fit
  fit_columns <- c("gene", "pvalue", "A", "phi", "sigma", "mesor")
  if (!is.data.frame(fit) || !all(fit_columns %in% names(fit)))
    stop("full_fit must contain all per-gene cosinor columns.")
  if (length(update$ngenes) != 1L || !is.finite(update$ngenes) ||
      update$ngenes < 1 || nrow(fit) != update$ngenes)
    stop("full_fit must cover the complete tested gene universe.")
  if (!all(vapply(fit[c("pvalue", "A", "phi", "sigma", "mesor")],
                  is.numeric, logical(1))))
    stop("Per-gene parameter columns must be numeric.")
  finite_p <- is.finite(fit$pvalue)
  if (any(fit$pvalue[finite_p] < 0 | fit$pvalue[finite_p] > 1))
    stop("Per-gene p-values must be in [0, 1].")

  candidates <- which(finite_p & fit$pvalue < 0.01 &
                        is.finite(fit$A) & fit$A > 0 &
                        is.finite(fit$sigma) & fit$sigma > 0)
  selected <- head(candidates[order(fit$pvalue[candidates])], 300L)
  r <- fit$A[selected] / fit$sigma[selected]
  if (any(!is.finite(r))) stop("Selected signal/noise ratios must be finite.")
  phases <- fit$phi[selected]
  phases <- phases[is.finite(phases)]

  q <- stats::p.adjust(fit$pvalue, method = "BH", n = update$ngenes)
  q <- q[is.finite(q)]
  p <- fit$pvalue[finite_p]
  quantile_or_na <- function(x, probability) {
    if (!length(x)) return(NA_real_)
    unname(stats::quantile(x, probability, type = 7))
  }
  if (length(update$prop_rhythmic) != 1L ||
      !is.finite(update$prop_rhythmic) ||
      update$prop_rhythmic < 0 || update$prop_rhythmic > 1)
    stop("prop_rhythmic must be a finite proportion.")
  expected_prop <- length(candidates) / update$ngenes
  if (!isTRUE(all.equal(as.numeric(update$prop_rhythmic), expected_prop,
                       tolerance = 1e-10)))
    stop("prop_rhythmic disagrees with the full-fit rhythmic candidate set.")

  # Phase columns are arithmetic quantiles in hours.
  # The same top <= 300 genes define both the phase and signal summaries.
  # These are descriptive quantiles, not circular phase estimators.
  c(
    list(n = length(update$times), ngenes = update$ngenes),
    gtex_signal_tod_summary(update$times),
    list(
      r_median_top300 = round(quantile_or_na(r, 0.5), 3),
      r_q25_top300 = round(quantile_or_na(r, 0.25), 3),
      r_q75_top300 = round(quantile_or_na(r, 0.75), 3),
      prop_rhythmic = round(update$prop_rhythmic, 4),
      phi_median_rhy = round(quantile_or_na(phases, 0.5), 2),
      phi_q25_rhy = round(quantile_or_na(phases, 0.25), 2),
      phi_q75_rhy = round(quantile_or_na(phases, 0.75), 2),
      rhy_FDR20 = sum(q < 0.20),
      rhy_FDR10 = sum(q < 0.10),
      rhy_FDR05 = sum(q < 0.05),
      rhy_FDR01 = sum(q < 0.01),
      rhy_p05 = sum(p < 0.05),
      rhy_p01 = sum(p < 0.01),
      rhy_p001 = sum(p < 0.001)
    )
  )
}


# 4. Replace only the requested CSV rows -------------------------------------

gtex_signal_csv_row <- function(row) {
  output <- character()
  con <- textConnection("output", "w", local = TRUE)
  utils::write.table(row, con, sep = ",", row.names = FALSE,
                     col.names = FALSE, quote = TRUE, na = "NA",
                     qmethod = "double")
  close(con)
  if (length(output) != 1L)
    stop("Signal-summary rows must not contain embedded line breaks.")
  output
}

gtex_signal_patch_csv <- function(source, destination, updates, statistics) {
  original <- utils::read.csv(source, stringsAsFactors = FALSE,
                              check.names = FALSE)
  lines <- readLines(source, warn = FALSE)
  if (length(lines) != nrow(original) + 1L)
    stop("CSV must contain exactly one physical line per row: ", source)
  if (!all(c("species", "dataset", "tissue", "condition") %in%
           names(original)))
    stop("CSV is missing the pilot identity columns: ", source)

  updated <- original
  changed <- integer()
  for (j in seq_along(updates)) {
    target <- gtex_signal_tissue_key(updates[[j]]$tissue)
    hit <- which(tolower(original$species) == "human" &
                   tolower(original$dataset) == "gtex" &
                   gtex_signal_tissue_key(original$tissue) == target &
                   tolower(original$condition) == "all")
    if (length(hit) != 1L)
      stop("Expected one CSV row for ", updates[[j]]$tissue,
           " in ", source, ". Found ", length(hit), ".")
    columns <- intersect(names(statistics[[j]]), names(original))
    for (column in columns) updated[hit, column] <- statistics[[j]][[column]]
    lines[hit + 1L] <- gtex_signal_csv_row(updated[hit, , drop = FALSE])
    changed <- c(changed, hit)
  }
  if (anyDuplicated(changed)) stop("A CSV row was targeted more than once.")

  unchanged <- setdiff(seq_len(nrow(original)), changed)
  if (!identical(original[unchanged, , drop = FALSE],
                 updated[unchanged, , drop = FALSE]))
    stop("An untargeted CSV value changed: ", source)
  if (!identical(names(original), names(updated)))
    stop("CSV columns changed: ", source)

  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, destination, useBytes = TRUE)
  written <- utils::read.csv(destination, stringsAsFactors = FALSE,
                             check.names = FALSE)
  if (!isTRUE(all.equal(updated, written, check.attributes = FALSE)))
    stop("CSV round-trip validation failed: ", destination)
  written_lines <- readLines(destination, warn = FALSE)
  if (!identical(readLines(source, warn = FALSE)[c(1L, unchanged + 1L)],
                 written_lines[c(1L, unchanged + 1L)]))
    stop("An untargeted CSV line changed: ", destination)

  data.frame(file = source, rows = nrow(original), rows_updated = length(changed),
             rows_preserved = length(unchanged), stringsAsFactors = FALSE)
}


# 5. Write the three staged summaries ----------------------------------------

update_gtex_signal_csvs <- function(updates, output_root, package_root = ".",
                                    expected_tissues = 31L) {
  if (!is.list(updates) || length(updates) != expected_tissues)
    stop("Expected ", expected_tissues, " rebuilt tissue records.")
  if (normalizePath(output_root, mustWork = FALSE) ==
      normalizePath(package_root, mustWork = TRUE))
    stop("Write summaries to a separate staging directory first.")
  keys <- vapply(updates, function(x) gtex_signal_tissue_key(x$tissue),
                 character(1))
  if (anyNA(keys) || any(!nzchar(keys)) || anyDuplicated(keys))
    stop("Rebuilt tissue identities must be present and unique.")
  if (any(keys %in% gtex_signal_protected_keys))
    stop("Adrenal, Liver, and Muscle summary rows must remain unchanged.")
  statistics <- lapply(updates, gtex_signal_update_statistics)

  # The package summary is public. The two accompanying paper CSVs are local.
  # Update them when present without requiring them in a fresh GitHub checkout.
  available <- file.exists(file.path(package_root, gtex_signal_csv_paths))
  if (!available[1L]) stop("The package signal-summary CSV is missing.")
  results <- lapply(gtex_signal_csv_paths[available], function(relative_path) {
    gtex_signal_patch_csv(
      file.path(package_root, relative_path),
      file.path(output_root, relative_path),
      updates, statistics
    )
  })
  result <- do.call(rbind, results)
  rownames(result) <- NULL
  result
}
