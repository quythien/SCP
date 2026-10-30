# 1. Read the complete recorded NAc and Caudate experiments.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1])
bundle <- dirname(dirname(normalizePath(script)))
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2L)
input <- normalizePath(args[1])
out <- args[2]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
power <- read.csv(file.path(input,"power_summary.csv"))
n80 <- read.csv(file.path(input,"N80_summary.csv"))
old_strata <- read.csv(file.path(input,"stratified_summary.csv"))
# Use the first internal empty interval to define a shared upper-tail bin.
pair_indices <- unique(power$index)
upper_limit <- 2
for (index in pair_indices) {
  a <- old_strata[old_strata$index == index, ]
  support <- tapply(a$targets, a$stratum, mean)
  occupied <- which(support > 0)
  gaps <- which(support == 0 & seq_along(support) > min(occupied) &
    seq_along(support) < max(occupied))
  if (length(gaps)) upper_limit <- min(upper_limit, (min(gaps) - 1) * .25)
}
tail_bin <- as.integer(upper_limit / .25) + 1L
labels <- c(sprintf("(%g,%g]", seq(0, upper_limit - .25, .25),
  seq(.25, upper_limit, .25)), paste0(">", upper_limit))
paths <- list.files(file.path(input,"replicates"),pattern="^[12]_N[0-9]+[.]rds$",full.names=TRUE)
stopifnot(length(paths) == 30L)
summaries <- checks <- list()

# 2. Combine the upper-tail counts within each saved experiment -----------
for (path in paths) {
  saved <- readRDS(path)
  d <- saved$strata
  stopifnot(length(unique(d$replicate)) == 600L)
  d$bin <- pmin(d$stratum, tail_bin)
  pooled <- aggregate(cbind(targets, TD) ~ replicate + bin, d, sum)
  pooled$power <- ifelse(pooled$targets > 0, pooled$TD / pooled$targets, NA_real_)
  original_totals <- aggregate(cbind(targets, TD) ~ replicate, d, sum)
  pooled_totals <- aggregate(cbind(targets, TD) ~ replicate, pooled, sum)
  stopifnot(identical(original_totals, pooled_totals))
  checks[[length(checks) + 1L]] <- data.frame(
    index = saved$index, N = saved$N, experiments = 600L,
    target_totals_identical = TRUE, discovery_totals_identical = TRUE,
    source_sha256 = digest::digest(path, file = TRUE, algo = "sha256"))
  for (bin in seq_along(labels)) {
    a <- pooled[pooled$bin == bin, ]
    valid <- is.finite(a$power)
    summaries[[length(summaries) + 1L]] <- data.frame(
      index = saved$index, key = saved$key, N = saved$N,
      stratum = bin, label = labels[bin],
      power = if (any(valid)) mean(a$power[valid]) else NA_real_,
      power_se = if (sum(valid) > 1L) sd(a$power[valid]) / sqrt(sum(valid)) else NA_real_,
      TD = mean(a$TD), TD_se = sd(a$TD) / sqrt(nrow(a)),
      targets = mean(a$targets), contributing_experiments = sum(valid))
  }
}
strata <- do.call(rbind, summaries)
unchanged <- merge(strata[strata$stratum < tail_bin, ], old_strata,
  by = c("index", "N", "stratum"), suffixes = c(".new", ".old"))
stopifnot(nrow(unchanged) == length(paths) * (tail_bin - 1L))
for (field in c("power", "power_se", "TD", "TD_se", "targets")) {
  a <- unchanged[[paste0(field, ".new")]]
  b <- unchanged[[paste0(field, ".old")]]
  stopifnot(identical(is.na(a), is.na(b)),
    all(abs(a[is.finite(a)] - b[is.finite(b)]) < 1e-10))
}
write.csv(strata, file.path(out, "stratified_summary.csv"), row.names = FALSE)
write.csv(power, file.path(out, "power_summary.csv"), row.names = FALSE)
write.csv(n80, file.path(out, "N80_summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, checks), file.path(out, "pooling_checks.csv"), row.names = FALSE)


# 3. Build the six panel matrices and apply the paper's N80 rule.
source(file.path(bundle,"Figure1/code/assemble_panels.R"))
left <- make_panels(2L,c(20,40,80,120))
right <- make_panels(1L,c(20,40,80,120))
full <- list(left,right)
for(i in seq_along(full)) {
  d <- full[[i]]
  p <- d$marg_mean[,2]
  j <- which(p>=.8)[1]
  d$vline_n <- if(is.na(j)) NA_real_ else if(j==1L) d$ss_a[1] else
    ceiling(d$ss_a[j-1]+(.8-p[j-1])/(p[j]-p[j-1])*(d$ss_a[j]-d$ss_a[j-1]))
  d$size_cols <- c("#CC5252","#C1CC52","#52ABCC","#7E52CC")
  full[[i]] <- d
}

# 4. Keep the selected leading limits for power and discovery panels.
power_panels <- discovery <- full
for(i in seq_along(full)) {
  for(kind in c("power","discovery")) {
    d <- full[[i]]
    first <- if(kind=="power") which(d$gene_counts>0)[1] else c(3L,1L)[i]
    keep <- seq.int(first,length(d$strata))
    d$strata <- d$strata[keep]
    d$gene_counts <- d$gene_counts[keep]
    for(field in c("pow","pow_se","td","td_se")) d[[field]] <- d[[field]][,keep,drop=FALSE]
    if(kind=="power") power_panels[[i]] <- d else discovery[[i]] <- d
  }
}
names(power_panels) <- c("NAc","Caudate")
saveRDS(power_panels,file.path(out,"panel_data.rds"))
saveRDS(discovery,file.path(out,"discovery_panel_data.rds"))
