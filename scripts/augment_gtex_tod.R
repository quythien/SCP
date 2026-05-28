#' Populate manifest tod_phases for GTEx pilots from gtex_v10_with_tod
#'
#' The GTEx pilots were ingested as gene-level cosinor summaries (no raw
#' times). The source `/home/qtp1/Projects/Collaborative/GTEXdata/data/v10/
#' gtex_v10_with_tod_121325.rds` keeps the per-donor time-of-death vector
#' under `[[tissue]]$tod` for each of the 68 GTEx v10 tissues. This script
#' reads it, builds the comma-list / summary form of unique times mod 24,
#' and writes the result into the manifest's tod_phases column so the
#' Shiny TOD distribution panel renders.

setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")

src <- "/home/qtp1/Projects/Collaborative/GTEXdata/data/v10/gtex_v10_with_tod_121325.rds"
cat(sprintf("Loading %s (~1.5 GB)\n", basename(src)))
gtex <- readRDS(src)
cat(sprintf("  %d tissues\n", length(gtex)))

# safe phase formatter — same convention as augment_tod_columns.R
safe_phases <- function(uniq_phases) {
  n <- length(uniq_phases)
  if (n == 0L || any(is.na(uniq_phases))) return(NA_character_)
  if (n <= 100)
    paste(formatC(uniq_phases, format = "g", digits = 4), collapse = ",")
  else
    sprintf("%s..%s (%d phases)",
            formatC(uniq_phases[1], format = "g", digits = 4),
            formatC(uniq_phases[n], format = "g", digits = 4), n)
}

# Mapping from GTEx source tissue name -> manifest tissue name
# (we used 'Brain - Cortex' -> 'Brain_Cortex' style in ingestion)
to_manifest_tissue <- function(nm) {
  out <- nm
  out <- gsub(" - ", "_", out, fixed = TRUE)
  out <- gsub("\\(|\\)", "", out)
  out <- gsub(" ", "_", out)
  out <- gsub("__+", "_", out)
  out <- gsub("_+$", "", out)
  out
}

manifest_path <- "inst/extdata/pilots/manifest.csv"
m <- read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
new_cols <- c("tod_n_unique","tod_cycles","tod_step_hr","tod_phases")
for (c in new_cols) if (!c %in% names(m)) m[[c]] <- if (c == "tod_phases") NA_character_ else NA

updated <- 0L
for (src_name in names(gtex)) {
  ti <- to_manifest_tissue(src_name)
  rows <- which(m$species == "human" & m$dataset == "GTEx" & m$tissue == ti)
  if (length(rows) == 0L) next
  times <- tryCatch(as.numeric(gtex[[src_name]]$tod), error = function(e) NULL)
  if (is.null(times) || length(times) == 0L) next
  times <- times[is.finite(times)]
  if (length(times) == 0L) next
  t_mod <- times %% 24
  uniq <- sort(unique(round(t_mod, 2)))
  n_unique <- length(uniq)
  span <- max(times) - min(times)
  step <- if (n_unique > 1) {
    dd <- diff(uniq); if (max(dd) - min(dd) < 0.5) round(median(dd), 1) else NA_real_
  } else NA_real_
  for (j in rows) {
    m$tod_n_unique[j] <- n_unique
    m$tod_cycles[j]   <- round(span / 24, 2)
    m$tod_step_hr[j]  <- step
    m$tod_phases[j]   <- safe_phases(uniq)
  }
  updated <- updated + length(rows)
  cat(sprintf("  %s -> %s: n=%d, n_unique=%d\n", src_name, ti, length(times), n_unique))
}

write.csv(m, manifest_path, row.names = FALSE)
cat(sprintf("\nUpdated %d GTEx manifest rows with tod_phases.\n", updated))
cat(sprintf("Manifest tod_phases populated on %d / %d ingested rows total.\n",
            sum(!is.na(m$tod_phases) & m$status == "ingested"),
            sum(m$status == "ingested")))
