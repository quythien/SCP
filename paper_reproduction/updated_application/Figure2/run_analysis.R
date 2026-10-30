# 1. Use an isolated installation of the frozen production engine.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
bundle <- dirname(normalizePath(script))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 2L)
library_path <- normalizePath(args[1], mustWork = TRUE)
output <- normalizePath(args[2], mustWork = FALSE)
.libPaths(c(library_path, .libPaths()))
library(SCP)
.CPP_LOADED <- TRUE
stopifnot(normalizePath(find.package("SCP")) == normalizePath(file.path(library_path, "SCP")))
source(file.path(bundle, "analysis/02_production_functions.R"))
reference <- jsonlite::read_json(file.path(bundle, "analysis/v2_model_manifest.json"), simplifyVector = TRUE)
body_hash <- function(name) digest::digest(paste(deparse(body(getFromNamespace(name,"SCP"))),
  collapse="\n"), algo="sha256", serialize=FALSE)
installed <- vapply(names(reference$function_body_sha256), body_hash, character(1))
stopifnot(identical(unname(installed), unname(unlist(reference$function_body_sha256))))
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")

# 2. Select recorded cases. Omit case IDs to run all 32,400 experiments.
cases <- read.csv(file.path(bundle, "analysis/cases.csv"), stringsAsFactors=FALSE)
if (length(args) > 2L) cases <- cases[cases$case_id %in% args[-(1:2)], , drop=FALSE]
stopifnot(nrow(cases) > 0L)
pilot <- readRDS(file.path(bundle, "analysis/Caudate_Putamen_pilot.rds"))
job <- jsonlite::read_json(file.path(bundle, "analysis/job.json"), simplifyVector=TRUE)
job$runtime_DLL_sha256 <- digest::digest(getLoadedDLLs()[["SCP"]][["path"]], file=TRUE, algo="sha256")
dir.create(file.path(output,"cases"), recursive=TRUE, showWarnings=FALSE)

# 3. Run each experiment with its recorded truth and expression seeds.
status <- vector("list", nrow(cases))
for (i in seq_len(nrow(cases))) {
  status[[i]] <- run_recorded_case(as.list(cases[i, , drop=FALSE]), pilot, job, output)
  cat(cases$case_id[i], status[[i]]$status, "\n")
}
write.csv(do.call(rbind,status), file.path(output,"case_status.csv"), row.names=FALSE)
stopifnot(all(vapply(status, function(x) x$status == "complete", logical(1))))
