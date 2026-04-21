#' Find the sample size needed to achieve a target genome-wide power
#'
#' Works on output from either \code{runSimsSingleCohort()} or \code{runSimsDiff()}.
#' Recomputes marginal power at the requested FDR from raw simulation output,
#' then returns the smallest grid n that meets the target and a linearly
#' interpolated estimate between grid points.
#'
#' @param res          Output from \code{runSimsSingleCohort()} or \code{runSimsDiff()}.
#' @param target_power Numeric in (0, 1). Target marginal power (default 0.80).
#' @param fdr          FDR threshold (default 0.05).
#' @param endpoint     For \code{runSimsDiff()} output: one of \code{"DR"}, \code{"DP"},
#'                     \code{"DM"}. Ignored for single-cohort output.
#' @param p.adjust.method Multiple-testing correction (default \code{"BH"}).
#'                     Only applied for single-cohort (differential uses pre-computed FDR).
#' @param interpolate  Logical. Linearly interpolate between grid points (default TRUE).
#' @return A list:
#'   \item{n_grid}{Smallest sample size in the simulation grid achieving target power, or NA.}
#'   \item{n_interp}{Linearly interpolated n (ceiling), or NA if power never reached.}
#'   \item{power}{Named numeric vector: mean power at each sample size.}
#'   \item{target_power, fdr, endpoint}{Echo of inputs.}
#' @export
npower <- function(res,
                   target_power     = 0.80,
                   fdr              = 0.05,
                   endpoint         = NULL,
                   p.adjust.method  = "BH",
                   interpolate      = TRUE) {

  sample_sizes <- res$sample_sizes
  nsims        <- res$nsims
  n_sizes      <- length(sample_sizes)
  is_diff      <- !is.null(res$diff_type)

  # ------------------------------------------------------------------
  # Compute mean power curve at the requested FDR
  # ------------------------------------------------------------------
  if (!is_diff) {
    # Single-cohort: recompute BH-adjusted power from raw p-values
    power_sim <- matrix(NA_real_, nrow = n_sizes, ncol = nsims)
    for (j in seq_len(n_sizes)) {
      for (s in seq_len(nsims)) {
        pvals       <- res$pvalues[j, , s]
        r_values    <- res$r_values_list[[j]][[s]]
        is_rhythmic <- r_values > 0
        n_tgt       <- sum(is_rhythmic)
        if (n_tgt == 0) next
        pvals[is.na(pvals)] <- 1
        qvals <- p.adjust(pvals, method = p.adjust.method)
        power_sim[j, s] <- sum(qvals <= fdr & is_rhythmic) / n_tgt
      }
    }

  } else {
    # Differential: use pre-computed FDR arrays
    if (is.null(endpoint) || !endpoint %in% c("DR", "DP", "DM"))
      stop("For runSimsDiff() output, specify endpoint as 'DR', 'DP', or 'DM'.")

    target_types <- switch(endpoint, DR = c(2L, 3L), DP = 4L, DM = 5L)
    fdr_arr      <- res[[paste0("fdr_", endpoint)]]   # [ngenes, n_sizes, nsims]

    power_sim <- matrix(NA_real_, nrow = n_sizes, ncol = nsims)
    for (s in seq_len(nsims)) {
      is_target <- res$diff_type[[s]] %in% target_types
      n_tgt     <- sum(is_target)
      if (n_tgt == 0) next
      for (j in seq_len(n_sizes)) {
        fdr_g <- fdr_arr[, j, s]
        power_sim[j, s] <- sum(!is.na(fdr_g) & fdr_g <= fdr & is_target) / n_tgt
      }
    }
  }

  mean_power          <- rowMeans(power_sim, na.rm = TRUE)
  names(mean_power)   <- sample_sizes

  # ------------------------------------------------------------------
  # Smallest grid n meeting target
  # ------------------------------------------------------------------
  above  <- which(mean_power >= target_power)
  n_grid <- if (length(above) > 0) sample_sizes[min(above)] else NA_integer_

  # ------------------------------------------------------------------
  # Linear interpolation between the two bracketing grid points
  # ------------------------------------------------------------------
  n_interp <- NA_real_
  if (interpolate && !is.na(n_grid)) {
    j2 <- min(above)
    if (j2 == 1L) {
      n_interp <- sample_sizes[1L]
    } else {
      j1 <- j2 - 1L
      p1 <- mean_power[j1]; p2 <- mean_power[j2]
      n1 <- sample_sizes[j1]; n2 <- sample_sizes[j2]
      n_interp <- ceiling(n1 + (target_power - p1) / (p2 - p1) * (n2 - n1))
    }
  }

  structure(
    list(
      n            = n_interp,
      n_grid       = n_grid,
      power        = mean_power,
      target_power = target_power,
      fdr          = fdr,
      endpoint     = endpoint
    ),
    class = "npower"
  )
}

#' @export
print.npower <- function(x, ...) {
  ep  <- if (!is.null(x$endpoint)) sprintf(" [%s]", x$endpoint) else ""
  cat(sprintf("npower%s — target: %.0f%% power at FDR %.0f%%\n",
              ep, x$target_power * 100, x$fdr * 100))
  if (is.na(x$n_grid)) {
    cat(sprintf("  Recommended n : not reached within simulation grid\n"))
    cat(sprintf("  Max power     : %.1f%% at n = %s\n",
                max(x$power, na.rm = TRUE) * 100,
                names(which.max(x$power))))
  } else {
    cat(sprintf("  Recommended n : %s\n",
                if (is.na(x$n)) as.character(x$n_grid) else as.character(x$n)))
  }
  cat("\n  Power curve:\n")
  for (nm in names(x$power))
    cat(sprintf("    n = %-4s  %.1f%%\n", nm, x$power[nm] * 100))
  invisible(x)
}
