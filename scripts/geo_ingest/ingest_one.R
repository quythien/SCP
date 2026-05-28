# Per-dataset ingestion helpers for SCP pilot DB.
# Sourced by run_batch.R. All paths absolute via PILOTS_ROOT env var.

suppressMessages({
  library(GEOquery)
  library(Biobase)
})

PILOTS_ROOT <- Sys.getenv("PILOTS_ROOT",
  "/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/PowerSim")
CACHE_DIR <- file.path(PILOTS_ROOT, "inst/extdata/pilots/cache_geo")
PILOT_DIR <- file.path(PILOTS_ROOT, "inst/extdata/pilots")
LOG_FILE  <- file.path(PILOT_DIR, "ingestion.log")
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)

# Source estimation code without loading full package machinery
source(file.path(PILOTS_ROOT, "R/estimation.R"))

log_msg <- function(...) {
  msg <- sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(..., collapse=""))
  cat(msg)
  cat(msg, file = LOG_FILE, append = TRUE)
}

# Parse TOD from a vector of strings/numbers -> numeric hours, or NA
parse_tod <- function(x) {
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
find_tod_col <- function(pd) {
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
      parsed <- parse_tod(vals)
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
        parsed <- parse_tod(sv)
        if (sum(!is.na(parsed)) >= 4) return(h)
      }
    }
  }
  NA_character_
}

# Build a pilot from expr + times, fitting circadian params
build_pilot <- function(expr, times, label = "") {
  expr <- as.matrix(expr)
  storage.mode(expr) <- "numeric"
  # Drop NA samples
  ok <- !is.na(times)
  expr <- expr[, ok, drop = FALSE]; times <- times[ok]
  if (ncol(expr) < 6) stop("too few samples: ", ncol(expr))
  # Drop genes with too many NAs / all-zero
  keep <- rowSums(is.na(expr)) < ncol(expr)/2 & matrixStats::rowVars(expr, na.rm=TRUE) > 0
  expr <- expr[keep, , drop = FALSE]
  # CPM normalize if looks like counts (integers, large range)
  if (max(expr, na.rm=TRUE) > 1e4 && all(expr[!is.na(expr)] >= 0)) {
    rs <- colSums(expr, na.rm=TRUE)
    if (all(rs > 0)) {
      expr <- t(t(expr) / rs) * 1e6
      expr <- log2(expr + 1)
    }
  }
  p <- estimate_circadian_params(expr, times, verbose = FALSE)
  if (is.null(p)) stop("estimate_circadian_params returned NULL")
  # Preserve times for downstream TOD augmentation
  if (is.null(p$raw)) p$raw <- list()
  p$raw$times <- times
  p$times <- times
  p
}

# Try to get an ExpressionSet
fetch_eset <- function(acc) {
  destdir <- file.path(CACHE_DIR, acc)
  dir.create(destdir, showWarnings = FALSE)
  res <- tryCatch(
    GEOquery::getGEO(acc, GSEMatrix = TRUE, destdir = destdir, AnnotGPL = FALSE),
    error = function(e) { log_msg("getGEO error ", acc, ": ", conditionMessage(e)); NULL }
  )
  if (is.null(res) || length(res) == 0) return(NULL)
  res
}

# Download supplementary files for an accession, return list of local paths
fetch_suppl <- function(acc) {
  destdir <- file.path(CACHE_DIR, acc, "suppl")
  dir.create(destdir, showWarnings = FALSE, recursive = TRUE)
  out <- tryCatch(
    GEOquery::getGEOSuppFiles(acc, baseDir = file.path(CACHE_DIR, acc), fetch_files = TRUE,
                              makeDirectory = FALSE),
    error = function(e) { log_msg("getGEOSuppFiles error ", acc, ": ", conditionMessage(e)); NULL }
  )
  if (is.null(out) || nrow(out) == 0) return(character(0))
  rownames(out)
}

# Parse a supplementary expression file into a numeric matrix (genes x samples)
parse_suppl_matrix <- function(path) {
  # Auto-detect delimiter from extension/content
  lower <- tolower(path)
  if (grepl("\\.rds$", lower)) {
    obj <- tryCatch(readRDS(path), error=function(e) NULL)
    if (is.matrix(obj) || is.data.frame(obj)) return(as.matrix(obj))
    return(NULL)
  }
  # Try a few delimiters/headers; use fread which handles dup row IDs gracefully
  if (!requireNamespace("data.table", quietly = TRUE)) return(NULL)
  for (sep in c("\t", ",", " ")) {
    df_full <- tryCatch(
      data.table::fread(path, header = TRUE, sep = sep,
                        check.names = FALSE, data.table = FALSE),
      error = function(e) NULL)
    if (is.null(df_full) || ncol(df_full) < 4) next
    # Take first column as row IDs (gene symbols / Ensembl IDs)
    row_ids <- as.character(df_full[[1]])
    df_full <- df_full[, -1, drop = FALSE]
    # Drop non-numeric columns (annotation cols)
    numcols <- sapply(df_full, function(x) {
      if (is.numeric(x)) return(TRUE)
      sv <- suppressWarnings(as.numeric(as.character(x[1:min(50, length(x))])))
      mean(!is.na(sv)) > 0.9
    })
    df_full <- df_full[, numcols, drop = FALSE]
    if (ncol(df_full) < 4) next
    mat <- suppressWarnings(as.matrix(sapply(df_full, function(x) as.numeric(as.character(x)))))
    # Handle duplicate row IDs by aggregating with max
    if (any(duplicated(row_ids))) {
      ord <- order(rowSums(mat, na.rm = TRUE), decreasing = TRUE)
      mat <- mat[ord, , drop = FALSE]; row_ids <- row_ids[ord]
      keep <- !duplicated(row_ids)
      mat <- mat[keep, , drop = FALSE]; row_ids <- row_ids[keep]
    }
    rownames(mat) <- row_ids
    if (ncol(mat) >= 4 && nrow(mat) >= 100) return(mat)
  }
  NULL
}

# Try to recover expression for a GEO series whose exprs() is empty by
# downloading suppl files and parsing the most likely one.
recover_expr_from_suppl <- function(acc, pd) {
  paths <- fetch_suppl(acc)
  if (length(paths) == 0) return(NULL)
  # Prefer counts/FPKM/TPM/CPM/normalized files
  pri <- grepl("(count|fpkm|tpm|cpm|rpkm|expression|matrix|normalized)", tolower(basename(paths)))
  paths <- c(paths[pri], paths[!pri])
  # Skip ones we know are non-matrix (BAM, fastq, BED, ...)
  paths <- paths[!grepl("\\.(bam|fastq|fq|bw|bed|wig|gtf|narrowpeak)(\\.gz)?$", tolower(paths))]
  for (p in paths) {
    if (file.size(p) > 5e8) next  # skip files > 500MB
    mat <- tryCatch(parse_suppl_matrix(p), error=function(e) NULL)
    if (!is.null(mat) && ncol(mat) >= 6) {
      # Try to align columns to pData rownames (GSM IDs)
      gsm_ids <- rownames(pd)
      cn <- colnames(mat)
      # Match by GSM prefix
      gsm_match <- sapply(cn, function(x) {
        hit <- gsm_ids[sapply(gsm_ids, function(g) grepl(g, x, fixed=TRUE) || grepl(x, g, fixed=TRUE))]
        if (length(hit) == 1) hit else NA_character_
      })
      if (sum(!is.na(gsm_match)) >= 6) {
        colnames(mat) <- gsm_match
        mat <- mat[, !is.na(colnames(mat)), drop = FALSE]
        return(list(expr = mat, file = basename(p)))
      }
      # Fall back: if column count matches sample count exactly, assume same order
      if (ncol(mat) == nrow(pd)) {
        colnames(mat) <- rownames(pd)
        return(list(expr = mat, file = basename(p)))
      }
      # Last resort: keep original column names; downstream code may parse TOD from them
      attr(mat, "tod_in_colnames") <- any(grepl("(?i)zt[0-9]|ct[0-9]|t[0-9]{1,2}[a-z_]", colnames(mat)))
      if (isTRUE(attr(mat, "tod_in_colnames"))) {
        return(list(expr = mat, file = basename(p), tod_in_colnames = TRUE))
      }
    }
  }
  NULL
}

# Parse TOD from sample column names like "ct18-young-muscle", "ZT06_rep1", "t12h_WT"
parse_tod_from_colnames <- function(cn) {
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
match_arm <- function(pd, tissue, condition) {
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

# Process all rows for one accession
process_accession <- function(acc, manifest_rows) {
  log_msg("=== ", acc, " (", nrow(manifest_rows), " arms) ===")
  results <- vector("list", nrow(manifest_rows))
  names(results) <- seq_len(nrow(manifest_rows))

  esets <- fetch_eset(acc)
  if (is.null(esets)) {
    for (i in seq_len(nrow(manifest_rows))) {
      results[[i]] <- list(status = "download_failed", note = "getGEO returned NULL", n=NA, ngenes=NA)
    }
    return(results)
  }

  # Use the first ExpressionSet (most series have one)
  eset <- esets[[1]]
  pd <- Biobase::pData(eset)
  expr <- Biobase::exprs(eset)

  # If exprs is empty (typical for RNA-seq series), try supplementary files
  if (nrow(expr) == 0 || all(is.na(expr))) {
    log_msg("  exprs() empty, trying suppl files for ", acc)
    rec <- tryCatch(recover_expr_from_suppl(acc, pd), error=function(e) {
      log_msg("  suppl recovery error: ", conditionMessage(e)); NULL })
    if (is.null(rec)) {
      for (i in seq_len(nrow(manifest_rows))) {
        results[[i]] <- list(status = "download_failed",
                             note = "exprs empty and no parseable suppl matrix", n=NA, ngenes=NA)
      }
      return(results)
    }
    expr <- rec$expr
    if (isTRUE(rec$tod_in_colnames)) {
      # Build synthetic pd from column names; TOD comes from the colname itself
      cn <- colnames(expr)
      tod_vec <- parse_tod_from_colnames(cn)
      keep <- !is.na(tod_vec)
      if (sum(keep) < 6) {
        for (i in seq_len(nrow(manifest_rows))) {
          results[[i]] <- list(status = "no_tod_metadata",
                               note = "suppl colnames non-GSM and no TOD parseable",
                               n=sum(keep), ngenes=NA)
        }
        return(results)
      }
      expr <- expr[, keep, drop = FALSE]
      tod_vec <- tod_vec[keep]
      pd <- data.frame(row.names = colnames(expr),
                       title = colnames(expr),
                       parsed_tod = tod_vec,
                       stringsAsFactors = FALSE)
      log_msg("  suppl ok (colname TOD): ", rec$file, " -> ", nrow(expr),
              " genes x ", ncol(expr), " samples; TOD range ",
              min(tod_vec), "-", max(tod_vec))
      # Force tod_col downstream
      tod_col <- "parsed_tod"
      all_times <- tod_vec
      # Proceed straight to arm matching below
      for (i in seq_len(nrow(manifest_rows))) {
        row <- manifest_rows[i, ]
        # Tissue/condition match on title col (use colnames as the text)
        arm_idx <- match_arm(pd, row$tissue, row$condition)
        if (length(arm_idx) < 6) {
          results[[i]] <- list(status = "ingest_failed",
                               note = sprintf("arm matched %d samples (<6) via colnames", length(arm_idx)),
                               n=length(arm_idx), ngenes=NA)
          next
        }
        expr_arm <- expr[, arm_idx, drop = FALSE]
        t_arm <- all_times[arm_idx]
        p <- tryCatch(build_pilot(expr_arm, t_arm, label = row$file),
                      error = function(e) { log_msg("build_pilot fail ", row$file, ": ", conditionMessage(e)); NULL })
        if (is.null(p)) {
          results[[i]] <- list(status = "ingest_failed", note = "build_pilot threw", n=length(t_arm), ngenes=NA)
          next
        }
        out_path <- file.path(PILOT_DIR, row$file)
        dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
        saveRDS(p, out_path)
        results[[i]] <- list(status = "ingested", note = "via colname TOD",
                             n = length(p$raw$times), ngenes = length(p$raw$pvalue))
        log_msg("  ok ", row$file, " n=", length(p$raw$times), " ngenes=", length(p$raw$pvalue))
      }
      return(results)
    }
    # Restrict pd to matched columns
    pd <- pd[colnames(expr)[colnames(expr) %in% rownames(pd)], , drop = FALSE]
    expr <- expr[, rownames(pd), drop = FALSE]
    log_msg("  suppl ok: ", rec$file, " -> ", nrow(expr), " genes x ", ncol(expr), " samples")
  }

  tod_col <- find_tod_col(pd)
  if (is.na(tod_col)) {
    for (i in seq_len(nrow(manifest_rows))) {
      results[[i]] <- list(status = "no_tod_metadata",
                           note = paste0("no TOD col among: ", paste(colnames(pd), collapse="|")),
                           n=NA, ngenes=NA)
    }
    return(results)
  }
  all_times <- parse_tod(pd[[tod_col]])

  for (i in seq_len(nrow(manifest_rows))) {
    row <- manifest_rows[i, ]
    arm_idx <- match_arm(pd, row$tissue, row$condition)
    if (length(arm_idx) < 6) {
      results[[i]] <- list(status = "ingest_failed",
                           note = sprintf("arm %s/%s matched %d samples (<6)", row$tissue, row$condition, length(arm_idx)),
                           n=length(arm_idx), ngenes=NA)
      next
    }
    expr_arm <- expr[, arm_idx, drop = FALSE]
    t_arm <- all_times[arm_idx]
    n_ok <- sum(!is.na(t_arm))
    if (n_ok < 6) {
      results[[i]] <- list(status = "no_tod_metadata",
                           note = sprintf("only %d samples with parseable TOD in arm", n_ok),
                           n=n_ok, ngenes=NA)
      next
    }
    p <- tryCatch(build_pilot(expr_arm, t_arm, label = row$file),
                  error = function(e) { log_msg("build_pilot fail ", row$file, ": ", conditionMessage(e)); NULL })
    if (is.null(p)) {
      results[[i]] <- list(status = "ingest_failed", note = "build_pilot threw", n=n_ok, ngenes=NA)
      next
    }
    out_path <- file.path(PILOT_DIR, row$file)
    dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
    saveRDS(p, out_path)
    results[[i]] <- list(status = "ingested", note = "",
                         n = length(p$raw$times), ngenes = length(p$raw$pvalue))
    log_msg("  ok ", row$file, " n=", length(p$raw$times), " ngenes=", length(p$raw$pvalue))
  }
  results
}
