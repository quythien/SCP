# 1. Locate the authorized expression source and create a new analysis directory.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1])
bundle <- dirname(normalizePath(script))
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2L)
input <- normalizePath(args[1],mustWork=TRUE)
output <- normalizePath(args[2],mustWork=FALSE)
stopifnot(!dir.exists(output))
dir.create(file.path(output,"output/two_harmonic/results"),recursive=TRUE)
dir.create(file.path(output,"submission/figures"),recursive=TRUE)

# 2. Use the saved fit and enrichment results with the original observed-data script.
files <- c("tissue_Liver_fits.rds","tissue_Liver_venn.rds","kegg_table.rds")
targets <- c("tissue_Liver_fits.rds","tissue_Liver_venn.rds","tissue_Liver_kegg.rds")
stopifnot(all(file.copy(file.path(bundle,"data",files),
  file.path(output,"output/two_harmonic/results",targets))))
Sys.setenv(SCP_ROOT=output,GTEX_V10_CPM=input)
source(file.path(bundle,"code/observed_liver.R"))
