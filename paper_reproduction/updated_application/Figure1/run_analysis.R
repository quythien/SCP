# 1. Read the saved NAc and Caudate calibration and time distributions.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1])
setwd(dirname(dirname(normalizePath(script))))
source("code/exact_F.R")
pilots <- readRDS("Figure1/data/calibrations.rds")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
out <- args[1]

G <- 5000L
nsims <- as.integer(Sys.getenv("SCP_REPLICATES", "600"))
thresholds <- c(.01, .05, .10, .20)
breaks <- c(seq(0, 3, .25), Inf)
labels <- c(sprintf("(%g,%g]", head(breaks, -2), breaks[2:(length(breaks) - 1)]), ">3")
active_N <- c(6L, 12L, 18L, 24L, 30L, 36L, 48L, 60L, 72L, 96L)
passive_N <- c(12L, 18L, 20L, 24L, 30L, 36L, 40L, 48L, 60L, 66L, 72L, 80L, 96L, 120L, 160L)
if (nzchar(Sys.getenv("SCP_N"))) passive_N <- as.integer(strsplit(Sys.getenv("SCP_N"), ",")[[1]])
dir.create(out, recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(list(G = G, nsims = nsims, FDR = thresholds,
  active_B = 6L, active_N = active_N, passive_N = passive_N,
  pilot_p = .01, top_k = 300L, paired_A_sigma = TRUE,
  input_sha256 = digest::digest(file = "Figure1/data/calibrations.rds", algo = "sha256"),
  interval_scope = "Fixed-pilot Monte Carlo error, separate from subject bootstrap uncertainty"),
  file.path(out, "settings.json"), pretty = TRUE, auto_unbox = TRUE)

# 2. Keep the same parameter law under a balanced active schedule ---------
draw_active <- function(bio, N, seed) {
  set.seed(seed, kind = "L'Ecuyer-CMRG")
  statistic_stream <- parallel::nextRNGStream(parallel::nextRNGStream(.Random.seed))
  nr <- round(bio$ngenes * bio$prop_rhythmic)
  truth <- logical(bio$ngenes); truth[sample(bio$ngenes, nr)] <- TRUE
  index <- sample(length(bio$amplitude), nr, replace = TRUE)
  A <- pmax(bio$amplitude[index], .05)
  sigma <- pmax(bio$sigma_rhythmic[index], 1e-6)
  phase <- sample(bio$phase, nr, replace = TRUE)
  times <- rep(seq(0, 20, 4), each = N / 6)
  lambda <- scp_f_noncentrality(A, phase, sigma, times, 24)
  stopifnot(max(abs(lambda - N * (A / sigma)^2 / 2)) < 1e-8)
  assign(".Random.seed", statistic_stream, envir = .GlobalEnv)
  p <- numeric(bio$ngenes)
  p[!truth] <- stats::runif(bio$ngenes - nr)
  p[truth] <- stats::pf(stats::rf(nr, 2, N - 3, ncp = lambda), 2, N - 3, lower.tail = FALSE)
  list(pvalue = p, truth = truth, amplitude = A, sigma = sigma, n_targets = nr)
}

# 3. Save the four-FDR and signal-stratum summaries for each experiment ---
dir.create(file.path(out, "replicates"), showWarnings = FALSE)
for (pilot in pilots) {
  bio <- list(ngenes = G, period = 24, prop_rhythmic = pilot$prop_rhythmic,
    amplitude = pilot$amplitude, sigma_rhythmic = pilot$sigma, phase = pilot$phase)
  sizes <- if (pilot$design == "active") active_N else passive_N
  for (N in sizes) {
    marg <- strat <- vector("list", nsims)
    for (i in seq_len(nsims)) {
      seed <- as.integer(20261029 + N * 100000 + i)
      value <- if (pilot$design == "active") draw_active(bio, N, seed) else
        simulate_scp_f_replicate(bio, NULL, N, seed, retain_parameters = TRUE,
          time_sampler = list(x = pilot$sampler$grid, probability = pilot$sampler$weights))
      q <- stats::p.adjust(value$pvalue, "BH")
      bin <- cut(value$amplitude / value$sigma, breaks, include.lowest = TRUE, labels = FALSE)
      count <- tabulate(bin, nbins = length(labels))
      marg[[i]] <- do.call(rbind, lapply(thresholds, function(alpha) {
        discovered <- q <= alpha
        td <- sum(discovered & value$truth); fd <- sum(discovered & !value$truth)
        data.frame(FDR_target = alpha, replicate = i, seed = seed,
          power = td / value$n_targets, TD = td, FD = fd, FDP = fd / max(td + fd, 1))
      }))
      td <- tabulate(bin[q[value$truth] <= .05], nbins = length(labels))
      strat[[i]] <- data.frame(replicate = i, stratum = seq_along(labels), label = labels,
        targets = count, TD = td, power = ifelse(count > 0, td / count, NA_real_))
    }
    file <- file.path(out, "replicates", paste0(pilot$index, "_N", N, ".rds"))
    stopifnot(!file.exists(file))
    saveRDS(list(key = pilot$key, index = pilot$index, N = N,
      marginal = do.call(rbind, marg), strata = do.call(rbind, strat)), file, compress = "xz")
    cat(pilot$key, "N", N, "complete\n"); flush.console()
  }
}

# 4. Summarize the completed experiments without broadening their effects --
summary <- list(); strat_summary <- list(); n80 <- list()
for (pilot in pilots) {
  paths <- list.files(file.path(out, "replicates"), pattern = paste0("^", pilot$index, "_N.*[.]rds$"), full.names = TRUE)
  for (path in paths) {
    x <- readRDS(path)
    for (alpha in thresholds) {
      a <- x$marginal[x$marginal$FDR_target == alpha, ]
      summary[[paste(pilot$index, x$N, alpha)]] <- data.frame(index = pilot$index, key = pilot$key,
        N = x$N, FDR_target = alpha, power = mean(a$power), power_se = sd(a$power) / sqrt(nrow(a)),
        mean_TD = mean(a$TD), mean_FDP = mean(a$FDP), nsims = nrow(a))
    }
    for (s in seq_along(labels)) {
      a <- x$strata[x$strata$stratum == s, ]
      valid <- is.finite(a$power)
      strat_summary[[paste(pilot$index, x$N, s)]] <- data.frame(index = pilot$index, key = pilot$key,
        N = x$N, stratum = s, label = labels[s],
        power = if (any(valid)) mean(a$power[valid]) else NA_real_,
        power_se = if (sum(valid) > 1) sd(a$power[valid]) / sqrt(sum(valid)) else NA_real_,
        TD = mean(a$TD), TD_se = sd(a$TD) / sqrt(nrow(a)), targets = mean(a$targets))
    }
  }
}
summary <- do.call(rbind, summary)
strat_summary <- do.call(rbind, strat_summary)
for (pilot in pilots) for (alpha in thresholds) {
  x <- summary[summary$index == pilot$index & summary$FDR_target == alpha, ]
  x <- x[order(x$N), ]; high <- which(x$power >= .8)[1]
  n80[[paste(pilot$index, alpha)]] <- data.frame(index = pilot$index, key = pilot$key,
    design = pilot$design, FDR_target = alpha,
    first_tested_N_at80 = if (is.finite(high)) x$N[high] else NA_integer_)
}
utils::write.csv(summary, file.path(out, "power_summary.csv"), row.names = FALSE)
utils::write.csv(strat_summary, file.path(out, "stratified_summary.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, n80), file.path(out, "N80_summary.csv"), row.names = FALSE)
