# 1. Locate this bundle when called from any working directory.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
root <- dirname(normalizePath(script, mustWork = TRUE))
setwd(root)

# 2. Redraw each figure in a separate process.
scripts <- c(paste0("Figure", 1:5, "/redraw.R"), "Supplement/redraw.R")
status <- vapply(scripts, function(script) system2(file.path(R.home("bin"), "Rscript"),
  shQuote(script)), integer(1))
stopifnot(all(status == 0L))
cat("Figures written to", file.path(root, "generated"), "\n")
