// cosinor_fast.cpp
//
// Rcpp/RcppArmadillo implementations of the two computational hot paths
// in the SCP PowerSim framework:
//
//   1. cosinor_pvals_cpp()
//      Vectorised cosinor OLS over G genes: fits y = M + A*cos(ωt) + B*sin(ωt)
//      for each gene via the 3x3 normal equations (identical maths to
//      one_cosinor_OLS() in estimation.R / fitCosinorAll() in bootstrap_sim.R).
//      Returns a numeric vector of F-test p-values, length G.
//
//   2. sim_cosinor_expr_cpp()
//      Generates a G x N expression matrix under the cosinor model with optional
//      2nd/3rd harmonics (identical maths to the for-g loop in simulation.R
//      simCircadianDiff() and runner.R runSimsSingleCohort()).
//      Returns a G x N numeric matrix.
//
// Usage from R (after Rcpp::sourceCpp("code/src/cosinor_fast.cpp")):
//
//   pvals <- cosinor_pvals_cpp(expr, times, period = 24)
//
//   expr  <- sim_cosinor_expr_cpp(mesor, amplitude, phase, sigma,
//                                  times, period = 24,
//                                  alpha2 = 0.0, alpha3 = 0.0)
//
// Design notes:
//   - Single-threaded (no OpenMP): parallelism is handled by mclapply at the
//     R level, so each worker calls these functions sequentially.
//   - The 3x3 solve is done analytically (closed form) to avoid LAPACK overhead
//     on a problem that is always exactly 3x3.
//   - NA/NaN inputs: if any element of y is non-finite, the gene returns pval=1.
//   - Degenerate case (RSS <= 0, n <= 3): returns pval = NA (as R code does).

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// ---------------------------------------------------------------------------
// Analytical 3x3 symmetric matrix inverse via cofactor expansion.
// Returns false if the matrix is singular (|det| < tol).
// ---------------------------------------------------------------------------
static inline bool solve3x3(const double S[9], const double d[3],
                              double out[3], double tol = 1e-14) {
  // S is stored row-major: S[3*row + col]
  // Cofactors of S (symmetric so we only need 6 unique entries)
  const double s00 = S[0], s01 = S[1], s02 = S[2];
  const double s11 = S[4], s12 = S[5];
  const double s22 = S[8];

  // Cofactors (upper-triangle only, using symmetry)
  const double C00 = s11 * s22 - s12 * s12;
  const double C01 = -(s01 * s22 - s12 * s02);
  const double C02 = s01 * s12 - s11 * s02;
  const double C11 = s00 * s22 - s02 * s02;
  const double C12 = -(s00 * s12 - s01 * s02);
  const double C22 = s00 * s11 - s01 * s01;

  const double det = s00 * C00 + s01 * C01 + s02 * C02;
  if (std::abs(det) < tol) return false;

  const double inv_det = 1.0 / det;
  // inv(S) * d  (S is symmetric so inv is also symmetric)
  out[0] = (C00 * d[0] + C01 * d[1] + C02 * d[2]) * inv_det;
  out[1] = (C01 * d[0] + C11 * d[1] + C12 * d[2]) * inv_det;
  out[2] = (C02 * d[0] + C12 * d[1] + C22 * d[2]) * inv_det;
  return true;
}

// ---------------------------------------------------------------------------
// cosinor_pvals_cpp
//
// Arguments:
//   expr   — G x N numeric matrix (genes x samples), row-major in R
//   times  — numeric vector of length N (sample collection times)
//   period — circadian period in hours (default 24)
//
// Returns:
//   Numeric vector of length G with cosinor F-test p-values.
//   Returns 1.0 for genes with non-finite data or degenerate fits.
//   Returns NA_real_ for genes where n <= 3 (insufficient df).
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector cosinor_pvals_cpp(const NumericMatrix& expr,
                                 const NumericVector& times,
                                 double period = 24.0) {
  const int G = expr.nrow();
  const int N = expr.ncol();

  if (N != times.size())
    stop("ncol(expr) must equal length(times)");
  if (N <= 3)
    return NumericVector(G, NA_REAL);

  // Pre-compute cos/sin basis and their cross-products (constant across genes)
  const double omega = 2.0 * M_PI / period;

  std::vector<double> c1(N), s1(N);
  double sum_c1 = 0, sum_s1 = 0;
  double sum_c1c1 = 0, sum_c1s1 = 0, sum_s1s1 = 0;

  for (int j = 0; j < N; ++j) {
    c1[j] = std::cos(omega * times[j]);
    s1[j] = std::sin(omega * times[j]);
    sum_c1   += c1[j];
    sum_s1   += s1[j];
    sum_c1c1 += c1[j] * c1[j];
    sum_c1s1 += c1[j] * s1[j];
    sum_s1s1 += s1[j] * s1[j];
  }

  // Normal-equations matrix (3x3, row-major, symmetric)
  // [ N       sum_c1   sum_s1   ]
  // [ sum_c1  sum_c1c1 sum_c1s1 ]
  // [ sum_s1  sum_c1s1 sum_s1s1 ]
  const double Smat[9] = {
    (double)N,  sum_c1,   sum_s1,
    sum_c1,     sum_c1c1, sum_c1s1,
    sum_s1,     sum_c1s1, sum_s1s1
  };

  NumericVector pvals(G);

  for (int g = 0; g < G; ++g) {
    // Accumulate y-related sums
    double sum_y = 0, sum_yc1 = 0, sum_ys1 = 0;
    double sum_y2 = 0;
    bool has_nan = false;

    for (int j = 0; j < N; ++j) {
      const double y = expr(g, j);
      if (!R_finite(y)) { has_nan = true; break; }
      sum_y   += y;
      sum_yc1 += y * c1[j];
      sum_ys1 += y * s1[j];
      sum_y2  += y * y;
    }

    if (has_nan) { pvals[g] = 1.0; continue; }

    const double dvec[3] = { sum_y, sum_yc1, sum_ys1 };
    double est[3];
    if (!solve3x3(Smat, dvec, est)) { pvals[g] = 1.0; continue; }

    // est[0] = m_hat, est[1] = beta1, est[2] = beta2 (referenced by index below)
    const double mean_y   = sum_y / N;

    // TSS = sum(y - mean_y)^2 = sum_y2 - N * mean_y^2
    const double TSS = sum_y2 - (double)N * mean_y * mean_y;

    // RSS = sum(y - yhat)^2, computed via
    //   RSS = sum_y2 - 2*(m_hat*sum_y + beta1*sum_yc1 + beta2*sum_ys1)
    //         + (normal equations evaluated at est)
    // Faster: RSS = TSS - MSS, where MSS = d' * est - N * mean_y^2 - m_hat^2 * ...
    // Simplest numerically stable form: RSS = TSS - beta1*sum_yc1 - beta2*sum_ys1
    //   (after centring; derivation via projection)
    // We use the direct form for clarity and numerical accuracy:
    //   RSS = sum_y2 - est' * d
    //   (since X'X * est = d  =>  est' * X'X * est = est' * d
    //    and RSS = y'y - est'X'Xy = sum_y2 - est'*d)
    const double RSS = sum_y2 - (est[0]*dvec[0] + est[1]*dvec[1] + est[2]*dvec[2]);
    const double MSS = TSS - RSS;

    if (RSS <= 0.0 || TSS <= 0.0) { pvals[g] = NA_REAL; continue; }

    const double Fstat = (MSS / 2.0) / (RSS / (double)(N - 3));
    if (!R_finite(Fstat) || Fstat < 0.0) { pvals[g] = 1.0; continue; }

    pvals[g] = R::pf(Fstat, 2.0, (double)(N - 3), 0, 0);  // lower.tail=FALSE
  }

  return pvals;
}


// ---------------------------------------------------------------------------
// sim_cosinor_expr_cpp
//
// Generates a G x N expression matrix under:
//   y_g(t) = M_g + A_g * [cos(ωt - ωφ_g)
//                         + alpha2 * cos(2ωt - 2ωφ_g)
//                         + alpha3 * cos(3ωt - 3ωφ_g)]
//             + epsilon_g,    epsilon_g ~ N(0, sigma_g^2)
//
// Arguments:
//   mesor     — length-G numeric vector of mesor values M_g
//   amplitude — length-G numeric vector of amplitudes A_g (0 for non-rhythmic)
//   phase     — length-G numeric vector of phases phi_g (hours)
//   sigma     — length-G numeric vector of residual SDs sigma_g
//   times     — length-N numeric vector of sample times
//   period    — circadian period (default 24)
//   alpha2    — 2nd-harmonic coefficient (default 0, i.e., no 2nd harmonic)
//   alpha3    — 3rd-harmonic coefficient (default 0, i.e., no 3rd harmonic)
//
// Returns:
//   G x N numeric matrix of simulated expression values.
//   Random draws use the R RNG (so set.seed() in R controls reproducibility).
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
NumericMatrix sim_cosinor_expr_cpp(const NumericVector& mesor,
                                    const NumericVector& amplitude,
                                    const NumericVector& phase,
                                    const NumericVector& sigma,
                                    const NumericVector& times,
                                    double period = 24.0,
                                    double alpha2 = 0.0,
                                    double alpha3 = 0.0) {
  const int G = mesor.size();
  const int N = times.size();

  if (amplitude.size() != G || phase.size() != G || sigma.size() != G)
    stop("mesor, amplitude, phase, sigma must all have length G");

  const double omega = 2.0 * M_PI / period;
  const bool use_h2 = (alpha2 != 0.0);
  const bool use_h3 = (alpha3 != 0.0);

  // Pre-compute trig terms per time point (shared across all genes)
  std::vector<double> cos1(N), sin1(N), cos2(N), cos3(N);
  for (int j = 0; j < N; ++j) {
    const double ot = omega * times[j];
    cos1[j] = std::cos(ot);
    sin1[j] = std::sin(ot);
    if (use_h2) cos2[j] = std::cos(2.0 * ot);
    if (use_h3) cos3[j] = std::cos(3.0 * ot);
  }

  NumericMatrix expr(G, N);

  for (int g = 0; g < G; ++g) {
    const double M  = mesor[g];
    const double A  = amplitude[g];
    const double ph = phase[g];
    const double sg = sigma[g];

    // Phase rotation: cos(ωt - ωφ) = cos(ωt)cos(ωφ) + sin(ωt)sin(ωφ)
    const double cos_phi1 = std::cos(omega * ph);
    const double sin_phi1 = std::sin(omega * ph);

    // 2nd harmonic phase rotation
    double cos_phi2 = 0.0, sin_phi2 = 0.0;
    if (use_h2) {
      cos_phi2 = std::cos(2.0 * omega * ph);
      sin_phi2 = std::sin(2.0 * omega * ph);
    }

    // 3rd harmonic phase rotation
    double cos_phi3 = 0.0;
    if (use_h3) {
      cos_phi3 = std::cos(3.0 * omega * ph);
      // sin(3ωφ) is not needed because cos(3ωt - 3ωφ) = cos3*cos_phi3 + sin3*sin_phi3,
      // but we store sin3 inside the loop only if needed. Recompute inline to avoid
      // the extra array allocation at cost of one trig call per gene per time point,
      // which is negligible compared to the N RNG calls.
    }

    for (int j = 0; j < N; ++j) {
      // Fundamental: cos(ωt - ωφ) = cos1*cos_phi1 + sin1*sin_phi1
      double mu = M + A * (cos1[j] * cos_phi1 + sin1[j] * sin_phi1);

      if (use_h2)
        mu += A * alpha2 * (cos2[j] * cos_phi2 + std::sin(2.0 * omega * times[j]) * sin_phi2);
      if (use_h3)
        mu += A * alpha3 * (cos3[j] * cos_phi3 + std::sin(3.0 * omega * times[j]) * (std::sin(3.0 * omega * ph)));

      expr(g, j) = R::rnorm(mu, sg);
    }
  }

  return expr;
}
