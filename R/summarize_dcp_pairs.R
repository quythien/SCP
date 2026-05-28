#' Summarize DCP pair results into a wide gene-level table
#'
#' Loads saved DCP_Analyze RDS files and computes a wide-format table with:
#'   - Sample and gene counts per tissue
#'   - TOD coverage per tissue
#'   - Single-tissue rhythmicity at multiple p/FDR thresholds
#'   - Joint rhythmicity (TOJR) breakdown
#'   - DR/DP/DA/DM counts at multiple thresholds
#'   - r = A/sigma for top-300 rhythmic genes (cross-tissue comparable)
#'   - Effect sizes: median |delta_phi| (DP), median |delta_M| (DM)
#'   - Recommendation flags (>= min_n_recommend genes per endpoint at FDR 5%)
#'
#' @param rds_paths   Named character vector: names = pair tags, values = RDS paths.
#' @param pair_meta   Data frame with columns: tag, tissue1, tissue2.
#' @param out_csv     Path for output CSV (NULL = return only).
#' @param min_n_recommend Minimum gene count to flag endpoint as sufficient (default 50).
#' @param p_thresholds  Raw p-value thresholds (default c(0.05, 0.01, 0.001)).
#' @param fdr_thresholds FDR thresholds (default c(0.20, 0.10, 0.05, 0.01)).
#' @param top_k       Number of top rhythmic genes (by p-value) for r statistics (default 300).
#' @param species     Character string for species column.
#' @param dataset     Character string for dataset column.
#' @return Data frame, one row per pair, sorted by sufficient_all3 then score_all3.
#' @export
summarizeDCPPairs <- function(rds_paths,
                              pair_meta        = NULL,
                              out_csv          = NULL,
                              min_n_recommend  = 50L,
                              p_thresholds     = c(0.05, 0.01, 0.001),
                              fdr_thresholds   = c(0.20, 0.10, 0.05, 0.01),
                              top_k            = 300L,
                              species          = "Human",
                              dataset          = "") {

  sidak2 <- 1 - (1 - 0.05)^(1 / 2)   # Sidak correction for 2 simultaneous post-hoc tests

  p_tags   <- sub("0\\.", "", as.character(p_thresholds))
  fdr_tags <- sub("0\\.", "", as.character(fdr_thresholds))

  rows <- lapply(seq_along(rds_paths), function(i) {
    tag      <- names(rds_paths)[i]
    rds_path <- rds_paths[i]

    if (!file.exists(rds_path)) {
      message(sprintf("[%s] RDS not found, skipping.", tag)); return(NULL)
    }
    message(sprintf("[%s] Processing...", tag))
    res <- readRDS(rds_path)

    rhy   <- res$rhythm
    dr    <- res$DR
    dp    <- res$DP
    clf   <- res$classification
    ngenes <- nrow(clf)

    rhy1  <- rhy$x1$rhythm
    rhy2  <- rhy$x2$rhythm
    joint <- rhy$rhythm.joint

    # ------------------------------------------------------------------
    # 1. Sample info & TOD coverage
    # ------------------------------------------------------------------
    n1 <- ncol(rhy$x1$data)
    n2 <- ncol(rhy$x2$data)
    t1 <- rhy$x1$time %% 24
    t2 <- rhy$x2$time %% 24

    # ------------------------------------------------------------------
    # 2. Per-tissue rhythmicity at multiple thresholds
    # ------------------------------------------------------------------
    rhy_counts <- function(pvals, label) {
      qvals <- p.adjust(pvals, "BH")
      out   <- c()
      for (k in seq_along(p_thresholds))
        out[sprintf("%s_p%s", label, p_tags[k])] <- sum(pvals < p_thresholds[k], na.rm = TRUE)
      for (k in seq_along(fdr_thresholds))
        out[sprintf("%s_FDR%s", label, fdr_tags[k])] <- sum(qvals < fdr_thresholds[k], na.rm = TRUE)
      out
    }
    rc1 <- rhy_counts(rhy1$pvalue, "rhy_g1")
    rc2 <- rhy_counts(rhy2$pvalue, "rhy_g2")

    # r = A/sigma for top-K rhythmic genes (comparable across tissues/species)
    top_r <- function(pvals, amps, sigmas) {
      ok  <- !is.na(pvals) & !is.na(amps) & !is.na(sigmas) & sigmas > 0 & amps > 0
      if (sum(ok) < 3) return(c(NA, NA, NA))
      idx <- order(pvals[ok])[seq_len(min(top_k, sum(ok)))]
      r   <- (amps[ok] / sigmas[ok])[idx]
      r   <- r[is.finite(r)]
      if (length(r) < 3) return(c(NA, NA, NA))
      c(median(r), quantile(r, 0.25, names=FALSE), quantile(r, 0.75, names=FALSE))
    }
    sig1 <- if ("sigma" %in% names(rhy1)) rhy1$sigma else rep(NA, nrow(rhy1))
    sig2 <- if ("sigma" %in% names(rhy2)) rhy2$sigma else rep(NA, nrow(rhy2))
    amp1 <- if ("A" %in% names(rhy1)) rhy1$A else if ("amp" %in% names(rhy1)) rhy1$amp else rep(NA, nrow(rhy1))
    amp2 <- if ("A" %in% names(rhy2)) rhy2$A else if ("amp" %in% names(rhy2)) rhy2$amp else rep(NA, nrow(rhy2))
    r1 <- top_r(rhy1$pvalue, amp1, sig1)
    r2 <- top_r(rhy2$pvalue, amp2, sig2)

    # ------------------------------------------------------------------
    # 3. TOJR joint rhythmicity breakdown
    # ------------------------------------------------------------------
    tojr_raw <- joint$TOJR
    tojr_fdr <- joint$TOJR.FDR
    n_both_raw  <- sum(tojr_raw == "both",  na.rm = TRUE)
    n_g1only    <- sum(tojr_raw == "rhyI",  na.rm = TRUE)
    n_g2only    <- sum(tojr_raw == "rhyII", na.rm = TRUE)
    n_arrhy     <- sum(tojr_raw == "arrhy", na.rm = TRUE)
    n_both_fdr  <- sum(tojr_fdr == "both",  na.rm = TRUE)

    # ------------------------------------------------------------------
    # 4. DR counts at multiple thresholds + directionality
    # ------------------------------------------------------------------
    dr_pvals <- if (!is.null(dr$p.R2)) dr$p.R2 else rep(NA_real_, ngenes)
    dr_qvals <- if (!is.null(dr$q.R2)) dr$q.R2 else p.adjust(dr_pvals, "BH")
    dr_counts <- c()
    for (k in seq_along(p_thresholds))
      dr_counts[sprintf("DR_p%s", p_tags[k])] <- sum(dr_pvals < p_thresholds[k], na.rm = TRUE)
    for (k in seq_along(fdr_thresholds))
      dr_counts[sprintf("DR_FDR%s", fdr_tags[k])] <- sum(dr_qvals < fdr_thresholds[k], na.rm = TRUE)

    dr_sig_idx <- which(dr_qvals < 0.05)
    tojr_at_dr <- tojr_raw[dr_sig_idx]
    dr_g1only_fdr05 <- sum(tojr_at_dr == "rhyI",  na.rm = TRUE)
    dr_g2only_fdr05 <- sum(tojr_at_dr == "rhyII", na.rm = TRUE)

    # ------------------------------------------------------------------
    # 5. DP counts at multiple thresholds + median |delta_phi|
    # ------------------------------------------------------------------
    dp_counts <- setNames(rep(NA_integer_, length(p_thresholds) + length(fdr_thresholds)),
                          c(sprintf("DP_p%s", p_tags), sprintf("DP_FDR%s", fdr_tags)))
    delta_phase_med <- NA_real_

    # DP uses "peak" terminology in DCP (peak time = phase)
    dp_phase_q  <- if ("q.delta.peak"  %in% names(dp)) dp$q.delta.peak  else dp$q.delta.phase
    dp_phase_p  <- if ("p.delta.peak"  %in% names(dp)) dp$p.delta.peak  else dp$p.delta.phase
    dp_delta_ph <- if ("delta.peak"    %in% names(dp)) dp$delta.peak    else dp$delta.phase

    if (!is.null(dp) && !is.null(dp$q.overall) && !is.null(dp_phase_q)) {
      q_ov <- dp$q.overall; p_ov <- dp$p.overall
      for (k in seq_along(fdr_thresholds))
        dp_counts[sprintf("DP_FDR%s", fdr_tags[k])] <-
          sum(!is.na(q_ov) & q_ov < fdr_thresholds[k] & dp_phase_q < sidak2, na.rm = TRUE)
      for (k in seq_along(p_thresholds))
        dp_counts[sprintf("DP_p%s", p_tags[k])] <-
          sum(!is.na(p_ov) & p_ov < p_thresholds[k] & dp_phase_p < sidak2, na.rm = TRUE)
      dp_sig <- !is.na(q_ov) & q_ov < 0.05 & dp_phase_q < sidak2
      if (sum(dp_sig, na.rm = TRUE) > 0 && !is.null(dp_delta_ph))
        delta_phase_med <- median(abs(dp_delta_ph[dp_sig]), na.rm = TRUE)
    }

    # ------------------------------------------------------------------
    # 6. DA (differential amplitude) counts
    # ------------------------------------------------------------------
    da_counts <- setNames(rep(NA_integer_, length(p_thresholds) + length(fdr_thresholds)),
                          c(sprintf("DA_p%s", p_tags), sprintf("DA_FDR%s", fdr_tags)))

    if (!is.null(dp) && !is.null(dp$q.overall) && !is.null(dp$q.delta.A)) {
      q_ov <- dp$q.overall; p_ov <- dp$p.overall
      q_A  <- dp$q.delta.A; p_A  <- dp$p.delta.A
      for (k in seq_along(fdr_thresholds))
        da_counts[sprintf("DA_FDR%s", fdr_tags[k])] <-
          sum(!is.na(q_ov) & q_ov < fdr_thresholds[k] & q_A < sidak2, na.rm = TRUE)
      for (k in seq_along(p_thresholds))
        da_counts[sprintf("DA_p%s", p_tags[k])] <-
          sum(!is.na(p_ov) & p_ov < p_thresholds[k] & p_A < sidak2, na.rm = TRUE)
    }

    # ------------------------------------------------------------------
    # 7. DM (differential mesor) — run DCP_DiffPar Par="M" on joint genes
    # ------------------------------------------------------------------
    dm_counts <- setNames(rep(NA_integer_, length(p_thresholds) + length(fdr_thresholds)),
                          c(sprintf("DM_p%s", p_tags), sprintf("DM_FDR%s", fdr_tags)))
    delta_M_med <- NA_real_

    if (n_both_raw >= 5) {
      dm_res <- tryCatch(
        DCP_DiffPar(rhy, Par = "M", alpha = 0.05, parallel.ncores = 1),
        error = function(e) NULL
      )
      if (!is.null(dm_res) && "pvalue" %in% names(dm_res)) {
        p_dm <- dm_res$pvalue
        q_dm <- p.adjust(p_dm, "BH")
        for (k in seq_along(p_thresholds))
          dm_counts[sprintf("DM_p%s", p_tags[k])] <- sum(p_dm < p_thresholds[k], na.rm = TRUE)
        for (k in seq_along(fdr_thresholds))
          dm_counts[sprintf("DM_FDR%s", fdr_tags[k])] <- sum(q_dm < fdr_thresholds[k], na.rm = TRUE)
        dm_sig <- !is.na(q_dm) & q_dm < 0.05
        if (sum(dm_sig) > 0 && "delta.Par" %in% names(dm_res))
          delta_M_med <- median(abs(dm_res$delta.Par[dm_sig]), na.rm = TRUE)
      }
    }

    # ------------------------------------------------------------------
    # 8. Recommendation flags (>= min_n_recommend at FDR 5%)
    # ------------------------------------------------------------------
    n_DR_rec <- dr_counts["DR_FDR05"];  n_DR_rec[is.na(n_DR_rec)] <- 0L
    n_DP_rec <- dp_counts["DP_FDR05"];  n_DP_rec[is.na(n_DP_rec)] <- 0L
    n_DM_rec <- dm_counts["DM_FDR05"];  n_DM_rec[is.na(n_DM_rec)] <- 0L

    flag_DR  <- n_DR_rec >= min_n_recommend
    flag_DP  <- n_DP_rec >= min_n_recommend
    flag_DM  <- n_DM_rec >= min_n_recommend

    # ------------------------------------------------------------------
    # 9. Assemble row
    # ------------------------------------------------------------------
    meta <- if (!is.null(pair_meta) && tag %in% pair_meta$tag)
              pair_meta[pair_meta$tag == tag, ]
            else
              data.frame(tag=tag, tissue1=NA_character_, tissue2=NA_character_,
                         stringsAsFactors=FALSE)

    as.data.frame(c(
      list(
        species  = species,  dataset = dataset,
        pair     = tag,
        tissue1  = meta$tissue1[[1]],
        tissue2  = meta$tissue2[[1]],
        n1 = n1, n2 = n2, ngenes = ngenes,
        # TOD coverage
        tod_min_g1  = round(min(t1), 1), tod_max_g1 = round(max(t1), 1),
        tod_sd_g1   = round(sd(t1), 2),
        tod_min_g2  = round(min(t2), 1), tod_max_g2 = round(max(t2), 1),
        tod_sd_g2   = round(sd(t2), 2),
        # r = A/sigma — top-300 rhythmic genes (cross-tissue comparable)
        r_median_top300_g1 = round(r1[1], 3),
        r_q25_top300_g1    = round(r1[2], 3),
        r_q75_top300_g1    = round(r1[3], 3),
        r_median_top300_g2 = round(r2[1], 3),
        r_q25_top300_g2    = round(r2[2], 3),
        r_q75_top300_g2    = round(r2[3], 3),
        # TOJR joint breakdown
        n_arrhy     = n_arrhy,
        n_g1only    = n_g1only,
        n_g2only    = n_g2only,
        n_both      = n_both_raw,
        n_both_fdr  = n_both_fdr,
        pct_jointly_rhythmic = round(100 * n_both_raw / ngenes, 2),
        # DR directionality at FDR 5%
        DR_FDR05_g1only = dr_g1only_fdr05,
        DR_FDR05_g2only = dr_g2only_fdr05,
        # Effect sizes
        delta_phase_median_h = round(delta_phase_med, 2),
        delta_M_median       = round(delta_M_med, 3),
        # Recommendation flags
        sufficient_DR   = flag_DR,
        sufficient_DP   = flag_DP,
        sufficient_DM   = flag_DM,
        sufficient_all3 = flag_DR & flag_DP & flag_DM
      ),
      as.list(rc1), as.list(rc2),
      as.list(dr_counts), as.list(dp_counts),
      as.list(da_counts), as.list(dm_counts)
    ), stringsAsFactors = FALSE)
  })

  rows   <- Filter(Negate(is.null), rows)
  result <- do.call(rbind, lapply(rows, function(r) as.data.frame(r, stringsAsFactors=FALSE)))

  if (nrow(result) > 0) {
    dr_col <- if ("DR_FDR05" %in% names(result)) result$DR_FDR05 else 0
    dp_col <- if ("DP_FDR05" %in% names(result)) result$DP_FDR05 else 0
    dm_col <- if ("DM_FDR05" %in% names(result)) result$DM_FDR05 else 0
    result$score_all3 <- ifelse(is.na(dr_col)|is.na(dp_col)|is.na(dm_col), 0L,
                                dr_col * dp_col * dm_col)
    result <- result[order(-result$sufficient_all3, -result$score_all3), ]
  }

  if (!is.null(out_csv)) {
    write.csv(result, out_csv, row.names = FALSE)
    message(sprintf("Saved: %s", out_csv))
  }
  invisible(result)
}
