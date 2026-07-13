#' Utility Functions for Circadian Power Analysis
#'



# QUICK ANALYTICAL POWER 
# ==============================================================================

#' Quick analytical power for single-group rhythmicity
#'
#' @description Fast closed-form power calculation based on F-test.
#' Useful for quick estimates before running full simulations.
#' Based on CircaPower's analytical approach.
#'
#' @param n Sample size (NULL to solve for n)
#' @param power Target power (NULL to solve for power)
#' @param r Effect size (A/sigma)
#' @param alpha Significance level
#' @param design_factor Sampling design factor (0.5 for even spacing)
#'
#' @return List with n, power, r, alpha
#' 
#' @examples
#' \dontrun{
#' # How many samples for 80% power?
#' CircaPower(n = NULL, power = 0.8, r = 1.5)
#'
#' # What power with n=24?
#' CircaPower(n = 24, power = NULL, r = 1.2)
#' }
#'
#' @keywords internal
CircaPower = function(n = NULL, power = NULL, r = NULL, 
                      alpha = 0.05, design_factor = 0.5) {
  
  # Exactly one parameter must be NULL
  if (sum(sapply(list(n, r, power, alpha), is.null)) != 1) {
    stop("Exactly one of n, r, power, and alpha must be NULL")
  }
  
  # Based on F-test with ncp = r^2 x n x d
  k = 2  # 2 parameters (cos, sin)
  df1 = k
  
  if (is.null(power)) {
    # Calculate power
    df2 = n - k - 1
    ncp = r^2 * n * design_factor
    F_crit = qf(1 - alpha, df1, df2, ncp = 0)
    power = 1 - pf(F_crit, df1, df2, ncp = ncp)
    
  } else if (is.null(n)) {
    # Solve for n
    n = uniroot(function(n) {
      df2 = n - k - 1
      ncp = r^2 * n * design_factor
      F_crit = qf(1 - alpha, df1, df2, ncp = 0)
      (1 - pf(F_crit, df1, df2, ncp = ncp)) - power
    }, c(4, 1e9))$root
    n = ceiling(n)
    
  } else if (is.null(r)) {
    # Solve for r
    r = uniroot(function(r) {
      df2 = n - k - 1
      ncp = r^2 * n * design_factor
      F_crit = qf(1 - alpha, df1, df2, ncp = 0)
      (1 - pf(F_crit, df1, df2, ncp = ncp)) - power
    }, c(1e-10, 1e5))$root
    
  } else if (is.null(alpha)) {
    # Solve for alpha
    alpha = uniroot(function(alpha) {
      df2 = n - k - 1
      ncp = r^2 * n * design_factor
      F_crit = qf(1 - alpha, df1, df2, ncp = 0)
      (1 - pf(F_crit, df1, df2, ncp = ncp)) - power
    }, c(1e-10, 1 - 1e-10))$root
  }
  
  return(list(n = n, power = power, r = r, alpha = alpha, 
              design_factor = design_factor))
}

# ==============================================================================
# DATA PREPARATION
# ==============================================================================

#' Prepare expression data for circadian power analysis
#'
#' @description
#' Converts raw expression data into the standard (matrix, times) format
#' expected by \code{estCircadianParam()}, \code{estCircadianParamTwoGroup()},
#' and \code{runBootstrapDesignGrid()}.
#'
#' The pipeline always expects:
#'   \itemize{
#'     \item \code{data}, numeric matrix, \strong{genes x samples}, log2-scale
#'     \item \code{times}, numeric vector, hours in [0, 24), same length as \code{ncol(data)}
#'   }
#'
#' This function handles the three most common raw formats:
#' \describe{
#'   \item{"counts"}{Raw integer counts (e.g. featureCounts, STAR).
#'     Converts to log2(CPM + 1): \code{log2(sweep(x, 2, colSums(x), "/") * 1e6 + 1)}.}
#'   \item{"cpm"}{Raw CPM values (e.g. edgeR \code{cpm()}, baboon data).
#'     Converts to log2(CPM + 1).}
#'   \item{"log2"}{Already log2-scale (e.g. GSE54651, Seney ACC).
#'     Wraps with \code{data.matrix()} only; no numeric transformation.}
#' }
#'
#' @param expr  A matrix, data.frame, or file path (CSV/TSV) of expression values.
#'              Rows = genes, columns = samples. Row names used as gene IDs.
#' @param times Numeric vector of time-of-day values (hours, 0-24) for each sample,
#'              OR a character column name in \code{pheno} to extract times from.
#' @param input_type One of \code{"counts"}, \code{"cpm"}, or \code{"log2"}.
#'              Default \code{"counts"}.
#' @param pheno Optional data.frame of sample metadata. Must have the same column
#'              order as \code{expr} if \code{times} is a column name string.
#' @param sample_col Optional: name of the column in \code{pheno} whose values match
#'              \code{colnames(expr)} for alignment. If NULL, assumes columns are
#'              already in the same order as rows of \code{pheno}.
#' @param time_col  Name of the TOD column in \code{pheno} (used when \code{times}
#'              is a column name string).
#' @param period  Circadian period in hours. Default 24. Used only for validation.
#' @param verbose Print a short summary of what was done. Default TRUE.
#'
#' @return A list with:
#'   \item{data}{Numeric matrix (genes x samples), log2-scale}
#'   \item{times}{Numeric vector of hours (0-24), length == ncol(data)}
#'   \item{n_genes}{Number of genes}
#'   \item{n_samples}{Number of samples}
#'   \item{input_type}{The \code{input_type} used}
#'
#' @examples
#' \dontrun{
#' # From raw counts matrix
#' prep <- prepCircadianData(counts_matrix, times = tod_vector, input_type = "counts")
#' bio  <- estCircadianParam(prep$data, prep$times)
#'
#' # From data.frame with metadata
#' prep <- prepCircadianData(
#'   expr       = expr_df,
#'   times      = "TOD",
#'   input_type = "log2",
#'   pheno      = pheno_df,
#'   sample_col = "SampleID"
#' )
#'
#' # Baboon raw CPM (load gives a data.frame)
#' prep <- prepCircadianData(baboon_df, times = tod_vector, input_type = "cpm")
#' }
#' @export
prepCircadianData <- function(expr,
                              times,
                              input_type = c("counts", "cpm", "log2"),
                              pheno      = NULL,
                              sample_col = NULL,
                              time_col   = NULL,
                              period     = 24,
                              verbose    = TRUE) {

  input_type <- match.arg(input_type)

  # --- 1. Load from file if a path was given ---
  if (is.character(expr) && length(expr) == 1 && file.exists(expr)) {
    if (verbose) cat(sprintf("Reading expression from: %s\n", expr))
    sep <- if (grepl("\\.tsv$|\\.txt$", expr, ignore.case=TRUE)) "\t" else ","
    expr <- as.matrix(read.csv(expr, row.names = 1, sep = sep, check.names = FALSE))
  }

  # --- 2. Coerce to numeric matrix ---
  mat <- data.matrix(expr)
  if (!is.numeric(mat))
    stop("Expression data must be numeric after coercion.")

  # --- 3. Align pheno and extract times ---
  if (is.character(times) && length(times) == 1) {
    # times is a column name, extract from pheno
    if (is.null(pheno))
      stop("'times' is a column name but 'pheno' was not provided.")
    time_col_use <- times

    if (!is.null(sample_col)) {
      # Align pheno rows to expression columns
      idx <- match(colnames(mat), as.character(pheno[[sample_col]]))
      if (any(is.na(idx)))
        warning(sprintf("%d expression columns could not be matched in pheno.", sum(is.na(idx))))
      pheno_aligned <- pheno[idx, , drop = FALSE]
    } else {
      pheno_aligned <- pheno
    }

    times_vec <- as.numeric(pheno_aligned[[time_col_use]])
  } else {
    times_vec <- as.numeric(times)
  }

  # --- 4. Validate alignment ---
  if (length(times_vec) != ncol(mat))
    stop(sprintf(
      "Length of times (%d) must equal ncol(expr) (%d).",
      length(times_vec), ncol(mat)))

  # Remove samples with missing times
  ok <- !is.na(times_vec)
  if (any(!ok)) {
    if (verbose)
      message(sprintf("Removing %d samples with NA times.", sum(!ok)))
    mat       <- mat[, ok, drop = FALSE]
    times_vec <- times_vec[ok]
  }

  # Warn if times are outside [0, 24)
  if (any(times_vec < 0 | times_vec >= 24, na.rm = TRUE))
    warning("Some times are outside [0, 24). Verify TOD extraction is correct.")

  # --- 5. Normalize ---
  if (input_type == "counts") {
    lib_size <- colSums(mat)
    if (any(lib_size == 0))
      stop("Some samples have zero library size. Check input.")
    mat <- log2(sweep(mat, 2, lib_size, "/") * 1e6 + 1)
    norm_desc <- "log2(CPM + 1) from raw counts"

  } else if (input_type == "cpm") {
    mat <- log2(mat + 1)
    norm_desc <- "log2(CPM + 1) from raw CPM"

  } else {
    # log2: already normalized, just ensure numeric matrix
    norm_desc <- "as-is (already log2-scale)"
  }

  # --- 6. Summary ---
  if (verbose) {
    message(sprintf(
      "prepCircadianData: %d genes x %d samples | input=%s | norm=%s",
      nrow(mat), ncol(mat), input_type, norm_desc))
    message(sprintf("  TOD range: [%.1f, %.1f] h  |  unique times: %d",
                    min(times_vec), max(times_vec), length(unique(times_vec))))
  }

  list(
    data       = mat,
    times      = times_vec,
    n_genes    = nrow(mat),
    n_samples  = ncol(mat),
    input_type = input_type
  )
}
