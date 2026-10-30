# Exact paired K=1 and K=2 Gaussian statistics for the same simulated dataset.

# 1. Draw the source parameters using the installed generator's order ------
draw_matched_parameters <- function(bio, G) {
  n_targets <- round(G * bio$prop_rhythmic)
  target_id <- sample(G, n_targets)
  truth <- logical(G)
  truth[target_id] <- TRUE
  mesor <- bio$lBaselineExpr
  sigma <- exp(bio$lOD)
  if (length(mesor) != G) {
    index <- sample(length(mesor), G, replace = TRUE)
    mesor <- mesor[index]
    sigma <- sigma[index]
  }
  A1 <- A2 <- phi1 <- phi2 <- numeric(G)
  joint_index <- sample.int(length(bio$amplitude), n_targets, replace = TRUE)
  A1[target_id] <- pmax(bio$amplitude[joint_index], .05)
  A2[target_id] <- pmax(bio$amplitude2[joint_index], 0)
  phi1[target_id] <- bio$phase[joint_index]
  phi2[target_id] <- bio$phase2[joint_index]
  sigma[target_id] <- pmax(bio$sigma_rhythmic[joint_index], 1e-6)
  stopifnot(length(mesor) == G, length(sigma) == G,
            all(is.finite(sigma)), all(sigma > 0))
  list(truth = truth, n_targets = n_targets, mesor = mesor, sigma = sigma,
       A1 = A1, A2 = A2, phi1 = phi1, phi2 = phi2)
}

# 2. Preserve passive sampling and use six exactly balanced phases --------
matched_schedule <- function(N, design, sampler) {
  stopifnot(length(N) == 1L, N >= 12L, N %% 6L == 0L)
  if (design == "passive") {
    return((sample(sampler$x, N, replace = TRUE, prob = sampler$probability) +
      rnorm(N, 0, .1)) %% 24)
  }
  stopifnot(design == "balanced6")
  rep(seq(0, 20, by = 4), each = N / 6L)
}

# 3. Retain the ordered nested harmonic subspaces -------------------------
matched_projection <- function(times) {
  omega <- 2 * pi / 24
  X <- cbind(1, cos(omega * times), sin(omega * times),
             cos(2 * omega * times), sin(2 * omega * times))
  decomposition <- qr(X, LAPACK = FALSE)
  stopifnot(length(times) > 5L, decomposition$rank == 5L,
            identical(decomposition$pivot, 1:5))
  list(X = X, Q = qr.Q(decomposition), R = qr.R(decomposition),
       decomposition = decomposition)
}

matched_mean_coordinates <- function(parameters, projection) {
  omega <- 2 * pi / 24
  beta <- cbind(parameters$A1 * cos(omega * parameters$phi1),
                parameters$A1 * sin(omega * parameters$phi1),
                parameters$A2 * cos(2 * omega * parameters$phi2),
                parameters$A2 * sin(2 * omega * parameters$phi2)) / parameters$sigma
  beta %*% t(projection$R[2:5, 2:5, drop = FALSE])
}

# 4. Test both detectors using the same four coordinates and residual -----
matched_pvalues <- function(coordinates, residual_ss, N) {
  S1 <- rowSums(coordinates[, 1:2, drop = FALSE]^2)
  S2 <- rowSums(coordinates[, 3:4, drop = FALSE]^2)
  F1 <- (S1 / 2) / ((S2 + residual_ss) / (N - 3))
  F2 <- ((S1 + S2) / 4) / (residual_ss / (N - 5))
  p1 <- pf(F1, 2, N - 3, lower.tail = FALSE)
  p2 <- pf(F2, 4, N - 5, lower.tail = FALSE)
  stopifnot(all(is.finite(p1)), all(is.finite(p2)))
  list(p1 = p1, p2 = p2)
}

# 5. Give parameters, passive times and each design's noise separate streams -
matched_streams <- function(seed) {
  set.seed(seed, kind = "L'Ecuyer-CMRG")
  parameters <- .Random.seed
  passive_times <- parallel::nextRNGStream(parameters)
  passive_noise <- parallel::nextRNGStream(passive_times)
  balanced6_noise <- parallel::nextRNGStream(passive_noise)
  list(parameters = parameters, passive_times = passive_times,
       passive_noise = passive_noise, balanced6_noise = balanced6_noise)
}

matched_metrics <- function(p, parameters, design, N) {
  q1 <- p.adjust(p$p1, "BH")
  q2 <- p.adjust(p$p2, "BH")
  rows <- lapply(c(.01, .05, .10, .20), function(alpha) {
    D1 <- q1 <= alpha
    D2 <- q2 <= alpha
    TD1 <- sum(D1 & parameters$truth)
    TD2 <- sum(D2 & parameters$truth)
    FD1 <- sum(D1 & !parameters$truth)
    FD2 <- sum(D2 & !parameters$truth)
    data.frame(design = design, N = N, FDR_target = alpha,
      n_targets = parameters$n_targets,
      power_K1 = TD1 / parameters$n_targets, power_K2 = TD2 / parameters$n_targets,
      paired_power_difference = (TD2 - TD1) / parameters$n_targets,
      TD_K1 = TD1, TD_K2 = TD2, FD_K1 = FD1, FD_K2 = FD2,
      FDP_K1 = FD1 / max(TD1 + FD1, 1L), FDP_K2 = FD2 / max(TD2 + FD2, 1L))
  })
  do.call(rbind, rows)
}

# 6. Share truth and parameter draws across designs, and data across detectors -
simulate_matched_pair <- function(bio, sampler, N, seed, G = 5000L,
                                  keep_fixture = FALSE) {
  streams <- matched_streams(seed)
  assign(".Random.seed", streams$parameters, envir = .GlobalEnv)
  parameters <- draw_matched_parameters(bio, G)
  stopifnot(parameters$n_targets > 0L, parameters$n_targets < G)
  rows <- fixtures <- list()
  for (design in c("passive", "balanced6")) {
    assign(".Random.seed", streams$passive_times, envir = .GlobalEnv)
    times <- matched_schedule(N, design, sampler)
    projection <- matched_projection(times)
    means <- matched_mean_coordinates(parameters, projection)
    assign(".Random.seed", streams[[paste0(design, "_noise")]], envir = .GlobalEnv)
    coordinates <- means + matrix(rnorm(G * 4L), nrow = G, ncol = 4L)
    residual_ss <- rchisq(G, df = N - 5L)
    p <- matched_pvalues(coordinates, residual_ss, N)
    rows[[design]] <- matched_metrics(p, parameters, design, N)
    if (keep_fixture) fixtures[[design]] <- list(times = times, truth = parameters$truth,
      parameters = parameters, coordinates = coordinates, residual_ss = residual_ss,
      pvalue_K1 = p$p1, pvalue_K2 = p$p2)
  }
  result <- list(metrics = do.call(rbind, rows))
  if (keep_fixture) result$fixture <- list(seed = seed, N = N,
    RNG_streams = streams, designs = fixtures)
  result
}

# 7. Verify portable files with SHA256 on macOS or Linux --------------------
matched_sha256 <- function(path) {
  stopifnot(file.exists(path))
  executable <- Sys.which("sha256sum")
  args <- shQuote(path)
  if (!nzchar(executable)) {
    executable <- Sys.which("shasum")
    args <- c("-a", "256", shQuote(path))
  }
  stopifnot(nzchar(executable))
  output <- system2(executable, args, stdout = TRUE)
  stopifnot(is.null(attr(output, "status")), length(output) == 1L)
  hash <- strsplit(output, "[[:space:]]+")[[1]][1]
  stopifnot(grepl("^[a-f0-9]{64}$", hash))
  hash
}

matched_bundle_files <- function() c("00_paired_engine.R", "01_validate_engine.R",
  "02_run_jobs.R", "03_summarize.R", "calibration_and_sampler.rds",
  "settings.rds", "job_manifest.csv", "source_provenance.csv", "validation_receipt.json")

matched_verify_gate <- function(bundle) {
  gate_path <- file.path(bundle, "independent_gate.csv")
  if (!file.exists(gate_path)) stop("Independent root review is required before production.")
  gate <- read.csv(gate_path, stringsAsFactors = FALSE)
  required <- matched_bundle_files()
  stopifnot(all(c("file", "sha256", "status", "reviewer") %in% names(gate)),
            !anyDuplicated(gate$file), setequal(gate$file, required),
            all(gate$status == "PASS"), all(nzchar(gate$reviewer)))
  current <- vapply(file.path(bundle, gate$file), matched_sha256, character(1))
  stopifnot(identical(unname(current), gate$sha256))
  list(gate_sha256 = matched_sha256(gate_path), hashes = gate)
}
