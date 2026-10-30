# 1. Read the recorded balanced-design jobs.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1])
bundle <- dirname(normalizePath(script))
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) >= 1L)
output <- normalizePath(args[1], mustWork=FALSE)
jobs <- read.csv(file.path(bundle,"analysis/job_manifest.csv"))
if (length(args) > 1L) jobs <- jobs[jobs$job_id %in% as.integer(args[-1]), , drop=FALSE]
stopifnot(nrow(jobs)>0L)

# 2. Recompute the selected Liver jobs using their frozen seeds.
for (id in jobs$job_id) {
  status <- system2(file.path(R.home("bin"),"Rscript"),
    c(shQuote(file.path(bundle,"analysis/02_run_job.R")), id, shQuote(output)))
  stopifnot(status == 0L)
}
cat("Caudate-Putamen C-E use Figure2/run_analysis.R and the active cases.\n")
