#' ====================================================================
#' Two-harmonic pilot fit helper (GTEx Adrenal Gland)
#'
#' Loads raw GTEx CPM for Adrenal Gland and runs estCircadianParam2H().
#' Caches the resulting CircadianBioOptions (carrying paired (A1, phi1,
#' A2, phi2, sigma) tuples) so the three K=2 figure scripts can share
#' a single empirical pilot.
#'
#' Sourced by fig4 / fig5 / fig6 two_harmonic scripts; not intended to
#' be Rscript'd directly.
#'
#' Exports:
#'   psi2_2h  CircadianBioOptions with paired_2h = TRUE
#'   pilot_expr  raw expression matrix (genes x samples) for Fig 4 Panel A
#'   pilot_times sample times (hours mod 24)
#' ====================================================================

.PILOT_RDS_2H <- "output/two_harmonic/results/pilot_2h_AdrenalGland.rds"
.PILOT_EXPR   <- "output/two_harmonic/results/pilot_2h_AdrenalGland_expr.rds"
.GTEX_CPM     <- "/home/qtp1/Projects/Collaborative/GTEXdata/TOD Xiangning/data/CPM/CPM.all.norm.RData"
.TOP_K_2H     <- 1500L   # broader pool dilutes signal so power curves don't saturate
.NGENES_TGT   <- 2000L
.MIN_RHYTHM_P <- 0.05    # less stringent pre-screen, matches a broader r distribution

if (!exists("psi2_2h", inherits = FALSE) || !exists("pilot_expr", inherits = FALSE)) {

  if (file.exists(.PILOT_RDS_2H) && file.exists(.PILOT_EXPR)) {
    cat(sprintf("Loading cached two-harmonic pilot: %s\n", .PILOT_RDS_2H))
    psi2_2h    <- readRDS(.PILOT_RDS_2H)
    .exprlist  <- readRDS(.PILOT_EXPR)
    pilot_expr  <- .exprlist$expr
    pilot_times <- .exprlist$times
  } else {
    cat(sprintf("Loading raw GTEx CPM: %s\n", .GTEX_CPM))
    load(.GTEX_CPM)

    df   <- CPM.all.norm[["Adrenal Gland"]]
    ids  <- as.character(colnames(df))
    hhmm <- sapply(strsplit(ids, "\\."), function(x) if (length(x) >= 3) x[3] else NA)
    hrs  <- as.numeric(substr(hhmm, 1, 2)) + as.numeric(substr(hhmm, 3, 4)) / 60
    ok   <- !is.na(hrs)
    pilot_expr  <- as.matrix(df[, ok])
    pilot_times <- hrs[ok]
    rm(CPM.all.norm, df)

    cat(sprintf("GTEx Adrenal: %d genes x %d samples\n",
                nrow(pilot_expr), ncol(pilot_expr)))

    cat("Fitting estCircadianParam2H() ...\n")
    t0 <- Sys.time()
    psi2_2h <- estCircadianParam2H(
      data            = pilot_expr,
      times           = pilot_times,
      period          = 24,
      min_rhythm_pval = .MIN_RHYTHM_P,
      top_k           = .TOP_K_2H,
      prop_DR         = 0,
      prop_DP         = 0,
      prop_DM         = 0,
      verbose         = TRUE
    )
    cat(sprintf("estCircadianParam2H runtime: %.1f s\n",
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))

    psi2_2h$ngenes <- .NGENES_TGT
    saveRDS(psi2_2h, .PILOT_RDS_2H)
    saveRDS(list(expr = pilot_expr, times = pilot_times), .PILOT_EXPR)
    cat(sprintf("Cached pilot -> %s\n", .PILOT_RDS_2H))
  }

  # Always reset ngenes to the figure-budget value (cached pilot may
  # have been fit before the override; this is cheap and idempotent).
  psi2_2h$ngenes <- .NGENES_TGT

  # ----------------------------------------------------------------
  # IMPORTANT: restore 2H pairing.
  # CircadianBioOptions() resamples `amplitude` to length n_rhythmic via
  # setAmplitude() (when paired_sigma=FALSE), independent of the 2H
  # vectors. The constructor leaves `amplitude2`/`phase2` at their
  # original top-K length (300). This silently breaks the joint
  # (A1, phi1, A2, phi2, sigma) pairing required by the 2H simulator,
  # so the runner draws `amplitude2[ji]` with ji > 300 -> NAs.
  # Workaround: overwrite the constructor-expanded fields with the
  # paired top-K vectors stored on $diagnostics. The runner then sees
  # length(amplitude)==length(amplitude2)==length(sigma_rhythmic) and
  # uses a single shared `sample.int(L_pilot, n_rhythmic, replace=TRUE)`
  # to draw paired tuples correctly.
  d <- psi2_2h$diagnostics
  stopifnot(!is.null(d), length(d$A1_emp) == length(d$A2_emp))
  psi2_2h$amplitude       <- d$A1_emp
  psi2_2h$amplitude2      <- d$A2_emp
  psi2_2h$phase           <- d$phi1_emp
  psi2_2h$phase2          <- d$phi2_emp
  psi2_2h$sigma_rhythmic  <- d$sigma_emp
  psi2_2h$paired_2h       <- TRUE

  # Cap simulated rhythmic proportion so FDR retains a non-trivial null
  # set. A pre-screen at p<0.05 inflates the empirical estimate.
  if (psi2_2h$prop_rhythmic > 0.30) psi2_2h$prop_rhythmic <- 0.30

  cat(sprintf("Pilot: ngenes(sim)=%d  L_pilot=%d  prop_rhythmic=%.3f  median r=%.2f\n",
              psi2_2h$ngenes, length(psi2_2h$amplitude),
              psi2_2h$prop_rhythmic,
              stats::median(psi2_2h$amplitude / psi2_2h$sigma_rhythmic)))
}
