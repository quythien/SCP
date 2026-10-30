# Exact Gaussian statistics for Figure 2.5 A/B balanced designs.

# 1. Preserve each panel's own calibrated parameter law --------------------
active_B <- function(panel) {
  stopifnot(panel %in% c("A", "B"))
  if (panel == "A") c(3L, 4L, 6L, 8L, 12L, 24L) else c(6L, 8L, 12L, 24L)
}

active_parameters_K1 <- function(bio, G) {
  n_targets <- round(G * bio$prop_rhythmic)
  target_id <- sample(G, n_targets)
  truth <- logical(G)
  truth[target_id] <- TRUE
  index <- sample(length(bio$amplitude), n_targets, replace = TRUE)
  amplitude <- pmax(bio$amplitude[index], .05)
  sigma_target <- pmax(bio$sigma_rhythmic[index], 1e-6)
  phase <- sample(bio$phase, n_targets, replace = TRUE)
  A1 <- phi1 <- numeric(G)
  sigma <- rep(1, G)
  A1[target_id] <- amplitude
  phi1[target_id] <- phase
  sigma[target_id] <- sigma_target
  list(truth = truth, n_targets = n_targets, mesor = numeric(G), sigma = sigma,
       A1 = A1, phi1 = phi1, A2 = numeric(G), phi2 = numeric(G),
       target_id = target_id, source_index = index)
}

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

# 2. Verify actual balanced design matrices against their common information -
active_geometry <- function(N, B, K) {
  stopifnot(N >= 24L, N %% 24L == 0L, N %% B == 0L, B > 2L * K,
            K %in% 1:2)
  times <- rep(seq(0, 24 - 24 / B, length.out = B), each = N / B)
  omega <- 2 * pi / 24
  X <- cbind(1, cos(omega * times), sin(omega * times))
  if (K == 2L) X <- cbind(X, cos(2 * omega * times), sin(2 * omega * times))
  expected <- diag(c(N, rep(N / 2, 2L * K)))
  error <- max(abs(crossprod(X) - expected))
  decomposition <- qr(X, LAPACK = FALSE)
  stopifnot(decomposition$rank == 2L * K + 1L, error < N * 1e-12,
            all(table(times) == N / B))
  list(times = times, X = X, decomposition = decomposition,
       rank = decomposition$rank, information_max_error = error)
}

active_coefficients <- function(parameters, K) {
  omega <- 2 * pi / 24
  coefficients <- cbind(parameters$A1 * cos(omega * parameters$phi1),
                         parameters$A1 * sin(omega * parameters$phi1))
  if (K == 2L) coefficients <- cbind(coefficients,
    parameters$A2 * cos(2 * omega * parameters$phi2),
    parameters$A2 * sin(2 * omega * parameters$phi2))
  coefficients
}

active_pvalues <- function(coordinates, residual_ss, N, K) {
  statistic <- (rowSums(coordinates^2) / (2L * K)) /
    (residual_ss / (N - 2L * K - 1L))
  p <- stats::pf(statistic, 2L * K, N - 2L * K - 1L, lower.tail = FALSE)
  stopifnot(all(is.finite(p)), all(p >= 0 & p <= 1))
  p
}

# 3. Couple the canonical Gaussian coordinates across all B within a panel -
simulate_active_panel <- function(bio, panel, N, seed, G = 5000L, keep_fixture = FALSE) {
  stopifnot(N >= 24L, N %% 24L == 0L)
  K <- if (panel == "A") 1L else 2L
  B_values <- active_B(panel)
  set.seed(seed, kind = "L'Ecuyer-CMRG")
  parameter_stream <- .Random.seed
  reserved_design_stream <- parallel::nextRNGStream(parameter_stream)
  noise_stream <- parallel::nextRNGStream(reserved_design_stream)
  parameters <- if (panel == "A") active_parameters_K1(bio, G) else draw_matched_parameters(bio, G)
  stopifnot(parameters$n_targets > 0L, parameters$n_targets < G)
  means <- sqrt(N / 2) * active_coefficients(parameters, K) / parameters$sigma
  assign(".Random.seed", noise_stream, envir = .GlobalEnv)
  coordinates <- means + matrix(rnorm(G * 2L * K), nrow = G)
  residual_ss <- rchisq(G, N - 2L * K - 1L)
  pvalue <- active_pvalues(coordinates, residual_ss, N, K)
  qvalue <- p.adjust(pvalue, "BH")
  rows <- list()
  for (B in B_values) for (alpha in c(.01, .05, .10, .20)) {
    discovery <- qvalue <= alpha
    TD <- sum(discovery & parameters$truth)
    FD <- sum(discovery & !parameters$truth)
    rows[[paste(B, alpha)]] <- data.frame(panel = panel, K = K, N = N, B = B,
      m = N / B, FDR_target = alpha, n_targets = parameters$n_targets,
      TD = TD, FD = FD, power = TD / parameters$n_targets,
      FDP = FD / max(TD + FD, 1L))
  }
  result <- list(metrics = do.call(rbind, rows))
  if (keep_fixture) result$fixture <- list(panel = panel, K = K, N = N, seed = seed,
    parameters = parameters, parameter_stream = parameter_stream, noise_stream = noise_stream,
    means = means, coordinates = coordinates, residual_ss = residual_ss, pvalue = pvalue)
  result
}

# 4. Require reviewed code, source references and frozen calibrations -------
active_sha256 <- function(path) {
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

active_gate_files <- function() c("00_exact_active_engine.R", "01_validate.R", "02_run_job.R",
  "03_summarize.R", "04_draw_panels.R", "calibrations.rds", "settings.rds",
  "job_manifest.csv", "source_provenance.csv", "validation_receipt.json")

active_verify_gate <- function(bundle) {
  path <- file.path(bundle, "independent_gate.csv")
  if (!file.exists(path)) stop("Independent root review is required before production.")
  gate <- read.csv(path, stringsAsFactors = FALSE)
  stopifnot(all(c("file", "sha256", "status", "reviewer") %in% names(gate)),
            !anyDuplicated(gate$file), setequal(gate$file, active_gate_files()),
            all(gate$status == "PASS"), all(nzchar(gate$reviewer)),
            identical(unname(vapply(file.path(bundle, gate$file), active_sha256, character(1))), gate$sha256))
  list(gate_sha256 = active_sha256(path), files = gate)
}
