# GEO matrix and annotation helpers.
# Time and arm choices are supplied in the reviewed sample plan.

native_read_expression <- function(path) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    expr <- readRDS(path)
  } else {
    stopifnot(requireNamespace("data.table", quietly = TRUE))
    table <- data.table::fread(path, data.table = FALSE, check.names = FALSE)
    gene <- as.character(table[[1]])
    expr <- as.matrix(table[, -1, drop = FALSE])
    storage.mode(expr) <- "numeric"
    rownames(expr) <- gene
  }
  stopifnot(is.matrix(expr), is.numeric(expr), !anyDuplicated(rownames(expr)))
  expr
}

native_probe_symbols <- function(feature_data) {
  columns <- grep("^gene.?symbol$|^symbol$|gene_symbol|^Gene Symbol$",
                  names(feature_data), ignore.case = TRUE, value = TRUE)
  if (length(columns)) {
    symbols <- as.character(feature_data[[columns[1]]])
  } else {
    column <- grep("gene_assignment|gene assignment", names(feature_data),
                    ignore.case = TRUE, value = TRUE)
    if (!length(column)) return(NULL)
    raw <- as.character(feature_data[[column[1]]])
    symbols <- trimws(sub("^[^/]*//\\s*([^/]+?)\\s*//.*", "\\1", raw))
    symbols[symbols == raw] <- NA_character_
  }
  names(symbols) <- rownames(feature_data)
  symbols[!is.na(symbols) & nzchar(symbols)]
}

# This parser is for explicit time labels, never identifiers or arbitrary columns.
native_time_labels <- function(values) {
  text <- trimws(tolower(as.character(values)))
  numeric_only <- grepl("^-?[0-9]+(\\.[0-9]+)?$", text)
  labelled <- grepl("^(zt|ct|time|hour|hr|h)[ :=_-]*-?[0-9]+(\\.[0-9]+)?h?$", text)
  ampm <- grepl("^[0-9]{1,2}(:[0-9]{2})?[ ]*(am|pm)$", text)
  result <- rep(NA_real_, length(text))
  result[numeric_only] <- as.numeric(text[numeric_only])
  result[labelled] <- as.numeric(gsub("[^0-9.-]", "", text[labelled]))
  for (i in which(ampm)) {
    clock <- sub("[ ]*(am|pm)$", "", text[i])
    pieces <- strsplit(clock, ":", fixed = TRUE)[[1]]
    hour <- as.numeric(pieces[1])
    minute <- if (length(pieces) == 2L) as.numeric(pieces[2]) else 0
    if (hour < 1 || hour > 12 || minute > 59) next
    result[i] <- (hour %% 12) + minute / 60 + if (grepl("pm$", text[i])) 12 else 0
  }
  result
}
