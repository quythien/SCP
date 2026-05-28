## Build summary/tissue_signal_summary.csv
##
## Walks every "ingested" pilot RDS under inst/extdata/pilots/ and
## extracts per-tissue signal/timing summary fields. Writes a single
## CSV with no NA rows for the canonical columns.
##
## Run from the package root:
##   Rscript scripts/build_tissue_signal_summary.R

suppressWarnings(suppressMessages({
  library(stats)
  library(utils)
}))

manifest_path <- "inst/extdata/pilots/manifest.csv"
pilots_root   <- "inst/extdata/pilots"
out_csv       <- "summary/tissue_signal_summary.csv"

stopifnot(file.exists(manifest_path))

man <- read.csv(manifest_path, stringsAsFactors = FALSE)
man <- man[man$status == "ingested", , drop = FALSE]
cat("Ingested manifest rows:", nrow(man), "\n")

summarize_times <- function(t) {
  if (length(t) == 0L) {
    return(list(min = NA_real_, max = NA_real_, sd = NA_real_,
                n_unique = NA_integer_, cycles = NA_real_,
                step_hr = NA_real_, phases = NA_character_))
  }
  tt <- t %% 24
  u  <- sort(unique(round(tt, 2)))
  diffs <- diff(u)
  step <- if (length(diffs) > 0) round(median(diffs), 3) else NA_real_
  cycles <- (max(tt) - min(tt)) / 24
  list(
    min = round(min(tt), 3),
    max = round(max(tt), 3),
    sd  = round(sd(tt), 3),
    n_unique = length(u),
    cycles = round(cycles, 3),
    step_hr = step,
    phases = paste(u, collapse = ",")
  )
}

## parse manifest tod_phases string back to numeric vector
parse_phases <- function(s) {
  if (is.null(s) || is.na(s) || !nzchar(s)) return(numeric(0))
  v <- suppressWarnings(as.numeric(strsplit(s, "[,\\s]+", perl = TRUE)[[1]]))
  v[is.finite(v)]
}

rows <- vector("list", nrow(man))
skipped <- character()
for (i in seq_len(nrow(man))) {
  r <- man[i, ]
  rds <- file.path(pilots_root, r$file)
  if (!file.exists(rds)) {
    message("Skip missing rds: ", rds)
    next
  }
  p <- tryCatch(readRDS(rds), error = function(e) NULL)
  if (is.null(p)) {
    message("Skip unreadable rds: ", rds)
    next
  }

  ## Times vector: prefer pilot$times, fall back to manifest$tod_phases
  t <- if (!is.null(p$times)) p$times[is.finite(p$times)] else numeric(0)
  if (length(t) == 0L) {
    t <- parse_phases(r$tod_phases)
  }
  if (length(t) == 0L) {
    ## Differential-only pilot with no TOD; not a tissue signal sample.
    skipped <- c(skipped, r$file)
    next
  }

  ## n samples
  n_samp <- length(t)
  if (!is.null(p$times) && length(p$times) > 0L) n_samp <- length(p$times)
  if (is.na(n_samp) || n_samp == 0L) n_samp <- suppressWarnings(as.integer(r$n))

  ## ngenes
  ngenes <- if (!is.null(p$raw$A)) length(p$raw$A) else NA_integer_
  if (is.na(ngenes) || ngenes == 0L) ngenes <- suppressWarnings(as.integer(r$ngenes))

  ## time-of-day summary
  tod <- summarize_times(t)

  ## r = amplitude / sigma  (top-300 vectors after rebuild)
  a <- p$amplitude
  s <- p$sigma_rhythmic
  n_min <- min(length(a), length(s))
  rv <- if (!is.null(a) && !is.null(s) && n_min > 0)
          a[seq_len(n_min)] / s[seq_len(n_min)]
        else NA_real_
  rv <- rv[is.finite(rv)]
  if (length(rv) > 0) {
    r_med <- as.numeric(median(rv))
    r_q25 <- as.numeric(quantile(rv, 0.25))
    r_q75 <- as.numeric(quantile(rv, 0.75))
  } else {
    r_med <- r_q25 <- r_q75 <- NA_real_
  }

  ## prop rhythmic
  pr <- p$prop_rhythmic
  if (is.null(pr) || !is.finite(pr)) {
    if (!is.null(p$raw$qvalue))
      pr <- mean(p$raw$qvalue < 0.05, na.rm = TRUE)
    else
      pr <- NA_real_
  }

  rows[[i]] <- data.frame(
    species          = r$species,
    dataset          = r$dataset,
    tissue           = r$tissue,
    tissue_canonical = r$tissue_canonical,
    design           = r$design,
    condition        = r$condition,
    n                = n_samp,
    ngenes           = ngenes,
    tod_min          = tod$min,
    tod_max          = tod$max,
    tod_sd           = tod$sd,
    tod_n_unique     = tod$n_unique,
    tod_cycles       = tod$cycles,
    tod_step_hr      = tod$step_hr,
    tod_phases       = tod$phases,
    r_median_top300  = round(r_med, 4),
    r_q25_top300     = round(r_q25, 4),
    r_q75_top300     = round(r_q75, 4),
    prop_rhythmic    = round(pr, 5),
    stringsAsFactors = FALSE
  )
}

rows <- Filter(Negate(is.null), rows)
out  <- do.call(rbind, rows)

dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
write.csv(out, out_csv, row.names = FALSE)

cat("Wrote", nrow(out), "rows to", out_csv, "\n")
if (length(skipped) > 0) {
  cat("Skipped", length(skipped),
      "pilots with no time-of-day vector (differential-only):\n")
  cat(paste(" -", skipped), sep = "\n")
}

## QC: NA counts per column
na_counts <- vapply(out, function(col) sum(is.na(col) | (is.character(col) & col == "")), integer(1))
cat("NA counts per column:\n")
print(na_counts)

rows_with_any_na <- sum(apply(out, 1, function(x) any(is.na(x) | x == "")))
cat("Rows with any NA/empty:", rows_with_any_na, "\n")
