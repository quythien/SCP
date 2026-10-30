# Draw the current application figures from their compact plotting data.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1])
root <- normalizePath(file.path(dirname(script), "../.."))
status <- system2(file.path(R.home("bin"), "Rscript"),
  shQuote(file.path(root,"paper_reproduction/updated_application/redraw_all.R")))
stopifnot(status==0L)
