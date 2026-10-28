# Prepare the two local normalized expression atlases.
source("examples/pilots/native/00_setup.R")

# 1. Baboon GSE98965 -----------------------------------------------------
if (any(native_index$dataset == "GSE98965")) {
  source_file <- native_input("SCP_BABOON_INPUT")
  source_env <- new.env()
  load(source_file, envir = source_env)
  atlas <- source_env$baboon_withTOD
  stopifnot(is.list(atlas$baboon), is.list(atlas$tod))
  rows <- unique(native_index[native_index$dataset == "GSE98965",
                              c("tissue", "base_file")])
  for (i in seq_len(nrow(rows))) {
    tissue <- rows$tissue[i]
    expr <- atlas$baboon[[tissue]]
    times <- atlas$tod[[tissue]]
    native_save_input(expr, times, rows$base_file[i], source_file,
                      "Normalized atlas expression and supplied ZT. No second transformation.")
  }
}

# 2. Mouse GSE54651 ------------------------------------------------------
if (any(native_index$dataset == "GSE54651")) {
  source_file <- native_input("SCP_MOUSE_RNASEQ_INPUT")
  atlas <- readRDS(source_file)
  stopifnot(is.list(atlas$count_clean), is.list(atlas$tod))
  rows <- unique(native_index[native_index$dataset == "GSE54651",
                              c("tissue", "base_file")])
  for (i in seq_len(nrow(rows))) {
    tissue <- rows$tissue[i]
    native_save_input(atlas$count_clean[[tissue]], atlas$tod[[tissue]],
                      rows$base_file[i], source_file,
                      "Prepared log2 CPM matrix and supplied ZT. No second transformation.")
  }
}
