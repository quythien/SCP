# Batch driver: walks manifest, processes each pending accession,
# updates manifest row-by-row, writes ingestion.log.
# Run via: Rscript scripts/geo_ingest/run_batch.R

PILOTS_ROOT <- "/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim"
Sys.setenv(PILOTS_ROOT = PILOTS_ROOT)
setwd(PILOTS_ROOT)

source("scripts/geo_ingest/ingest_one.R")

MANIFEST <- "inst/extdata/pilots/manifest.csv"

# Add note column if missing
m <- read.csv(MANIFEST, stringsAsFactors = FALSE)
if (!"note" %in% colnames(m)) m$note <- ""

# Ordering: prefer Batch A (recent RNA-seq) accessions first per user spec
batch_A <- c("GSE201207","GSE190939","GSE198015","GSE197524","GSE197199","GSE131063",
"GSE104336","GSE145883","GSE165856","GSE159960","GSE164047","GSE149092","GSE182856",
"GSE193660","GSE201862","GSE70497","GSE83855","GSE38622","GSE39977","GSE36975",
"GSE104674","GSE76896","GSE45642","GSE114918","GSE113883","GSE182740","GSE125521",
"GSE48113","GSE39445","GSE99905","GSE158422","GSE193677","E-MTAB-9606","GSE143038",
"GSE127919","GSE161617","GSE188475","GSE140983","GSE164900","GSE89773","GSE169501","GSE46069")
batch_B <- c("GSE10045","GSE10644","GSE11516","GSE11922","GSE11923","GSE27366","GSE34018",
"GSE35026","GSE49638","GSE52333","GSE54650","GSE57830","GSE73222","GSE78215","GSE84580",
"GSE1654","GSE25612","GSE20635","GSE8988","GSE8989","GSE5417","GSE51220")
ordered_accs <- c(batch_A, batch_B)

pending <- m[m$status == "pending", ]
pending_accs <- unique(pending$accession)
# Order pending by batch
ord <- c(intersect(ordered_accs, pending_accs),
         setdiff(pending_accs, ordered_accs))

cat("[", format(Sys.time(), "%H:%M:%S"), "] starting batch:",
    length(ord), "accessions,", nrow(pending), "arms\n", sep="")
cat("[", format(Sys.time(), "%H:%M:%S"), "] batch start: ",
    length(ord), " accs / ", nrow(pending), " arms\n", sep="",
    file = file.path(PILOTS_ROOT, "inst/extdata/pilots/ingestion.log"),
    append = TRUE)

n_done <- 0
checkpoint_every <- 5

write_manifest <- function(m) {
  write.csv(m, MANIFEST, row.names = FALSE)
}

for (acc in ord) {
  # Skip controlled or anything not pending anymore (in case of resume)
  idx <- which(m$accession == acc & m$status == "pending")
  if (length(idx) == 0) next
  rows <- m[idx, ]

  res <- tryCatch(process_accession(acc, rows),
    error = function(e) {
      log_msg("FATAL ", acc, ": ", conditionMessage(e))
      lapply(seq_len(nrow(rows)), function(i) list(status="ingest_failed",
        note=substr(conditionMessage(e),1,180), n=NA, ngenes=NA))
    })

  for (k in seq_along(idx)) {
    r <- res[[k]]
    m$status[idx[k]] <- r$status
    m$note[idx[k]]   <- r$note
    if (!is.null(r$n) && !is.na(r$n)) m$n[idx[k]] <- r$n
    if (!is.null(r$ngenes) && !is.na(r$ngenes)) m$ngenes[idx[k]] <- r$ngenes
    if (r$status == "ingested") n_done <- n_done + 1
  }
  # Write manifest after every accession for resume safety
  write_manifest(m)

  if (n_done > 0 && (n_done %% checkpoint_every) == 0) {
    msg <- sprintf("checkpoint: %d arms ingested so far", n_done)
    cat("[", format(Sys.time(), "%H:%M:%S"), "] ", msg, "\n", sep="")
    cat("[", format(Sys.time(), "%H:%M:%S"), "] ", msg, "\n", sep="",
        file = file.path(PILOTS_ROOT, "inst/extdata/pilots/ingestion.log"),
        append = TRUE)
  }
}

# Final summary
cat("\n=== FINAL ===\n")
print(table(m$status, useNA="always"))
cat("\nBy species (pending pool only):\n")
pend_orig <- pending
final_pend <- m[m$file %in% pend_orig$file, ]
print(table(final_pend$species, final_pend$status))

cat("\nFAILURES:\n")
fail <- m[m$file %in% pend_orig$file & m$status != "ingested", c("accession","tissue","condition","status","note")]
if (nrow(fail) > 0) print(fail, row.names=FALSE)
write_manifest(m)
cat("\nManifest written:", MANIFEST, "\n")
