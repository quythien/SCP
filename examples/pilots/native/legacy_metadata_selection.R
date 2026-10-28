# Earlier metadata-selection helpers retained for tracing candidate sample plans.
# Their results require review. The numbered builders do not call them automatically.
# Source identity and original hashes are recorded in builder_source_index.csv.

# Parse TOD from a vector of strings/numbers -> numeric hours, or NA
native_legacy_parse_time <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  x <- as.character(x)
  # Strip prefixes
  out <- suppressWarnings(as.numeric(gsub("[^0-9.\\-]", "",
           gsub("(?i)ZT|CT|time[^0-9]*|hour[^0-9]*|hr[^0-9]*|h$|am$|pm$|:[0-9]+", "", x, perl=TRUE))))
  # Handle AM/PM
  has_pm <- grepl("(?i)pm", x); has_am <- grepl("(?i)am", x)
  out[has_pm & !is.na(out) & out < 12] <- out[has_pm & !is.na(out) & out < 12] + 12
  out[has_am & !is.na(out) & out == 12] <- 0
  out
}

# Find TOD column in pData
native_legacy_find_time_column <- function(pd) {
  cols <- colnames(pd)
  patterns <- c("\\btime\\b", "\\btod\\b", "\\bZT\\b", "\\bCT\\b",
                "circadian", "\\bhour", "timepoint", "time.?point",
                "sampling.?time", "harvest.?time", "collection.?time",
                "\\bzt[0-9]", "time.of.day")
  for (p in patterns) {
    hits <- grep(p, cols, ignore.case = TRUE, value = TRUE)
    # Prefer characteristics_* cols which actually carry the value
    for (h in hits) {
      vals <- pd[[h]]
      parsed <- native_legacy_parse_time(vals)
      ok <- sum(!is.na(parsed) & parsed >= -12 & parsed <= 48)
      if (ok >= 4) return(h)
    }
  }
  # Fallback: scan every column
  for (h in cols) {
    vals <- pd[[h]]
    if (is.numeric(vals)) {
      ok <- sum(!is.na(vals) & vals >= -12 & vals <= 48)
      if (ok >= 6 && length(unique(vals)) >= 3) return(h)
    } else {
      sv <- as.character(vals)
      if (any(grepl("(?i)ZT[0-9]|CT[0-9]|hours?[: ]|h$", sv))) {
        parsed <- native_legacy_parse_time(sv)
        if (sum(!is.na(parsed)) >= 4) return(h)
      }
    }
  }
  NA_character_
}


# Parse TOD from sample column names like "ct18-young-muscle", "ZT06_rep1", "t12h_WT"
native_legacy_time_from_names <- function(cn) {
  # Try ZT/CT prefix
  hr1 <- suppressWarnings(as.numeric(sub("(?i).*?(zt|ct|t|h)[\\.\\-_ ]?([0-9]{1,2}(\\.[0-9]+)?).*", "\\2", cn, perl = TRUE)))
  # Also handle "hour6", "h6"
  hr2 <- suppressWarnings(as.numeric(sub("(?i).*?(hour|h)[\\.\\-_ ]?([0-9]{1,2}).*", "\\2", cn, perl = TRUE)))
  hr <- ifelse(!is.na(hr1) & hr1 >= 0 & hr1 < 48, hr1,
        ifelse(!is.na(hr2) & hr2 >= 0 & hr2 < 48, hr2, NA_real_))
  hr
}

# Map a manifest tissue+condition to a subset of samples.
# Heuristic: look at pData for matches; if not found, return all samples.
native_legacy_match_arm <- function(pd, tissue, condition) {
  n <- nrow(pd)
  keep <- rep(TRUE, n)
  # Tissue matching
  if (!is.na(tissue) && nzchar(tissue) && tolower(tissue) != "all") {
    tcols <- grep("source|tissue|organ|cell", colnames(pd), ignore.case=TRUE, value=TRUE)
    if (length(tcols) > 0) {
      tissue_terms <- unlist(strsplit(tissue, "[_/ ]"))
      tissue_terms <- tissue_terms[nchar(tissue_terms) >= 2]
      hits <- rep(FALSE, n)
      for (tc in tcols) {
        for (tt in tissue_terms) {
          hits <- hits | grepl(tt, pd[[tc]], ignore.case = TRUE)
        }
      }
      if (sum(hits) >= 6) keep <- keep & hits
    }
  }
  # Condition matching (genotype/treatment/condition columns)
  if (!is.na(condition) && nzchar(condition) && tolower(condition) != "all") {
    ccols <- grep("genotype|treatment|condition|group|strain|status|diagnosis",
                  colnames(pd), ignore.case=TRUE, value=TRUE)
    if (length(ccols) > 0) {
      cond_terms <- unlist(strsplit(condition, "[_/ ]"))
      cond_terms <- cond_terms[nchar(cond_terms) >= 2]
      hits <- rep(FALSE, n)
      for (cc in ccols) {
        for (ct in cond_terms) {
          hits <- hits | grepl(ct, pd[[cc]], ignore.case = TRUE)
        }
      }
      if (sum(hits) >= 6) keep <- keep & hits
    }
  }
  which(keep)
}
