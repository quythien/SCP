make_panels <- function(index, display) {
  marginal <- power[power$index == index, ]
  s <- strata[strata$index == index, ]
  sizes <- sort(unique(marginal$N))
  thresholds <- c(.01, .05, .10, .20)
  strata_ids <- sort(unique(s$stratum))
  marginal_mean <- marginal_se <- matrix(NA_real_, length(sizes), 4)
  for (j in seq_along(sizes)) for (k in seq_along(thresholds)) {
    x <- marginal[marginal$N == sizes[j] & marginal$FDR_target == thresholds[k], ]
    marginal_mean[j, k] <- x$power
    marginal_se[j, k] <- x$power_se
  }
  selected <- which(sizes %in% display)
  p <- p_se <- td <- td_se <- matrix(NA_real_, length(selected), length(strata_ids))
  for (j in seq_along(selected)) for (k in seq_along(strata_ids)) {
    x <- s[s$N == sizes[selected[j]] & s$stratum == strata_ids[k], ]
    p[j, k] <- x$power; p_se[j, k] <- x$power_se
    td[j, k] <- x$TD; td_se[j, k] <- x$TD_se
  }
  list(ss_a = sizes, marg_mean = marginal_mean, marg_se = marginal_se,
    vline_n = n80$first_tested_N_at80[n80$index == index & n80$FDR_target == .05],
    thresh_cols = c("darkgreen", "steelblue", "orange", "red"),
    thresh_labels = c("FDR 1%", "FDR 5%", "FDR 10%", "FDR 20%"),
    fdr_label = "FDR 5%", strata = s$label[match(strata_ids, s$stratum)],
    ss_bc = sizes[selected], size_cols = rainbow(length(sizes), s = .6, v = .8)[selected],
    pow = p, pow_se = p_se, td = td, td_se = td_se,
    gene_counts = vapply(strata_ids, function(k) mean(s$targets[s$stratum == k]), numeric(1)))
}
