# Final summary after batch ingestion completes.
setwd("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")

m <- read.csv("inst/extdata/pilots/manifest.csv", stringsAsFactors = FALSE)
csv_path <- "output/supp_tissue_summary/tissue_signal_summary.csv"
csv <- if (file.exists(csv_path)) read.csv(csv_path, stringsAsFactors = FALSE) else NULL

cat("=== MANIFEST FINAL STATE ===\n")
print(table(m$status, useNA = "always"))

cat("\n=== INGESTED BY SPECIES (manifest) ===\n")
print(table(m$species, m$status))

# Restrict to the original 78 pending pool (anything that was pending at start)
# We approximate that by: not in original 84 ingested (GTEx + baboon + Seney + etc.)
# Just report counts on all non-controlled rows.
batch_pool <- m[!m$status %in% c("controlled", "skipped_low_priority"), ]
new_pool <- m[m$status %in% c("ingested", "download_failed", "no_tod_metadata", "ingest_failed"), ]

cat("\n=== PER-DATASET OUTCOME (newly processed only) ===\n")
# Identify newly processed by looking for non-empty note OR new ingestions w/ small n
new_acc <- unique(m$accession[m$status %in% c("download_failed","no_tod_metadata","ingest_failed")])
new_ok  <- m[m$accession %in% c(new_acc, unique(m$accession[!is.na(m$note) & m$note != ""])), ]
if (nrow(new_ok) > 0) {
  show <- new_ok[, c("accession","species","tissue","condition","n","ngenes","status","note")]
  show$note <- ifelse(is.na(show$note), "", substr(show$note, 1, 80))
  print(show, row.names = FALSE)
}

cat("\n=== CSV STATE ===\n")
if (!is.null(csv)) {
  cat("Rows:", nrow(csv), "\n")
  cat("By species:\n"); print(table(csv$species))
}

cat("\n=== TOP FAILURE MODES ===\n")
fails <- m[m$status %in% c("download_failed","no_tod_metadata","ingest_failed"), ]
if (nrow(fails) > 0) {
  fails$note_short <- substr(fails$note, 1, 60)
  print(table(fails$status))
  cat("\nNote distribution:\n")
  print(sort(table(fails$note_short), decreasing = TRUE))
}
