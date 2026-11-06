// C++ kernels for the multidimensional hierarchical model (mhgrm).
//
// Latent vector theta_i (d <= 3) with prior N(mu_i, Sigma_i), where Sigma_i is
// parameterized by the modified Cholesky decomposition T_i Sigma_i T_i' = D_i:
//   eps_1 = e_1,  eps_2 = phi_21 eps_1 + e_2,  eps_3 = phi_31 eps_1 + phi_32 eps_2 + e_3,
//   e_k ~ N(0, sig2_ik), eps = theta - mu.
// Items load on one dimension each (simple structure), so the likelihood on the
// product grid is the sum of per-dimension log-likelihoods on the 1-D grids.
//
// Quadrature is adaptive: on dimension a, person i's nodes are
// G_a(k, i) = c_ia + h_ia * g_k, where g_k are standard nodes with log weights
// lqw_k, and c_ia and h_ia are the current posterior mean and SD. With c = 0 and
// h = 1 the grid is the same for everyone (non-adaptive quadrature).

#include <Rcpp.h>
using namespace Rcpp;

static inline double plogis_md(double x) {
  if (x == R_PosInf) return 1.0;
  if (x == R_NegInf) return 0.0;
  return 1.0 / (1.0 + std::exp(-x));
}

static inline double dlogis_md(double x) {
  if (!R_FINITE(x)) return 0.0;
  double p = plogis_md(x);
  return p * (1.0 - p);
}

// Log-likelihood of a set of items on person-specific nodes G (K x N):
// LL(i, k) = sum_j log P(y_ij | theta = G(k, i)). Missing responses are skipped.
// [[Rcpp::export]]
NumericMatrix loglik_nodes_cpp(const IntegerMatrix& Y, const NumericMatrix& alpha,
                               const NumericVector& beta, const NumericMatrix& G) {
  const int N = Y.nrow(), J = Y.ncol(), K = G.nrow();
  NumericMatrix LL(N, K);
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < J; ++j) {
      const int y = Y(i, j);
      if (y == NA_INTEGER) continue;
      const double lo = alpha(j, y - 1), hi = alpha(j, y);
      for (int k = 0; k < K; ++k) {
        const double u = beta[j] * G(k, i);
        LL(i, k) += std::log(plogis_md(lo + u) - plogis_md(hi + u));
      }
    }
  return LL;
}

// E-step on the (person-specific) product grid.
//   LL   : list of d matrices (N x K), log-likelihood of each dimension's items
//   mu   : N x d prior means;  sig2: N x d innovation variances
//   phi  : N x 3 Cholesky coefficients in the order (phi_21, phi_31, phi_32)
//   G    : list of d node matrices (K x N);  lsc: N x d log node scales h_ia
//   lqw  : log weights of the standard nodes
// Returns posterior means (N x d), posterior second moments E[theta_a theta_b]
// (N x d(d+1)/2, pairs (1,1),(1,2),(2,2),(1,3),(2,3),(3,3) in column-major
// upper-triangular order), marginal weights (list of d K x N matrices), and the
// log marginal likelihood.
// [[Rcpp::export]]
List estep_md_cpp(const List& LL, const NumericMatrix& mu, const NumericMatrix& sig2,
                  const NumericMatrix& phi, const List& G, const NumericMatrix& lsc,
                  const NumericVector& lqw) {
  const int d = LL.size(), N = mu.nrow(), K = lqw.size();
  if (d < 1 || d > 3) stop("d must be 1, 2, or 3");
  const int K2 = d >= 2 ? K : 1, K3 = d >= 3 ? K : 1;
  NumericMatrix L1 = LL[0];
  NumericMatrix L2 = d >= 2 ? as<NumericMatrix>(LL[1]) : NumericMatrix(N, 1);
  NumericMatrix L3 = d >= 3 ? as<NumericMatrix>(LL[2]) : NumericMatrix(N, 1);
  NumericMatrix G1 = G[0];
  NumericMatrix G2 = d >= 2 ? as<NumericMatrix>(G[1]) : NumericMatrix(1, N);
  NumericMatrix G3 = d >= 3 ? as<NumericMatrix>(G[2]) : NumericMatrix(1, N);
  const int npair = d * (d + 1) / 2;
  NumericMatrix m(N, d), S2(N, npair);
  NumericMatrix W1(K, N), W2(K2, N), W3(K3, N);
  std::vector<double> lw(K * K2 * K3);
  const double log2pi = std::log(2.0 * M_PI);
  double loglik = 0.0;

  for (int i = 0; i < N; ++i) {
    const double s1 = sig2(i, 0), m1 = mu(i, 0);
    const double s2 = d >= 2 ? sig2(i, 1) : 1.0, m2 = d >= 2 ? mu(i, 1) : 0.0;
    const double s3 = d >= 3 ? sig2(i, 2) : 1.0, m3 = d >= 3 ? mu(i, 2) : 0.0;
    const double p21 = phi(i, 0), p31 = phi(i, 1), p32 = phi(i, 2);
    double cst = -0.5 * (log2pi + std::log(s1)) + lsc(i, 0);
    if (d >= 2) cst += -0.5 * (log2pi + std::log(s2)) + lsc(i, 1);
    if (d >= 3) cst += -0.5 * (log2pi + std::log(s3)) + lsc(i, 2);
    double mx = R_NegInf;
    int idx = 0;
    for (int k1 = 0; k1 < K; ++k1) {
      const double e1 = G1(k1, i) - m1;
      const double a1 = cst - 0.5 * e1 * e1 / s1 + L1(i, k1) + lqw[k1];
      for (int k2 = 0; k2 < K2; ++k2) {
        double e2 = 0.0, a2 = a1;
        if (d >= 2) {
          e2 = G2(k2, i) - m2;
          const double r2 = e2 - p21 * e1;
          a2 += -0.5 * r2 * r2 / s2 + L2(i, k2) + lqw[k2];
        }
        for (int k3 = 0; k3 < K3; ++k3) {
          double a3 = a2;
          if (d >= 3) {
            const double r3 = (G3(k3, i) - m3) - p31 * e1 - p32 * e2;
            a3 += -0.5 * r3 * r3 / s3 + L3(i, k3) + lqw[k3];
          }
          lw[idx++] = a3;
          if (a3 > mx) mx = a3;
        }
      }
    }
    // nodes with relative weight below exp(-40) contribute nothing measurable;
    // setting them to zero also avoids slow subnormal arithmetic
    double tot = 0.0;
    for (int t = 0; t < idx; ++t) {
      const double v = lw[t] - mx;
      lw[t] = v < -40.0 ? 0.0 : std::exp(v);
      tot += lw[t];
    }
    loglik += mx + std::log(tot);
    // moments and marginals
    double e[3] = {0, 0, 0}, ee[6] = {0, 0, 0, 0, 0, 0};
    idx = 0;
    for (int k1 = 0; k1 < K; ++k1)
      for (int k2 = 0; k2 < K2; ++k2)
        for (int k3 = 0; k3 < K3; ++k3) {
          const double w = lw[idx++] / tot;
          if (w == 0.0) continue;
          const double t1 = G1(k1, i), t2 = d >= 2 ? G2(k2, i) : 0.0, t3 = d >= 3 ? G3(k3, i) : 0.0;
          W1(k1, i) += w;
          e[0] += w * t1; ee[0] += w * t1 * t1;
          if (d >= 2) {
            W2(k2, i) += w;
            e[1] += w * t2; ee[1] += w * t1 * t2; ee[2] += w * t2 * t2;
          }
          if (d >= 3) {
            W3(k3, i) += w;
            e[2] += w * t3; ee[3] += w * t1 * t3; ee[4] += w * t2 * t3; ee[5] += w * t3 * t3;
          }
        }
    for (int a = 0; a < d; ++a) m(i, a) = e[a];
    for (int r = 0; r < npair; ++r) S2(i, r) = ee[r];
  }
  List marg(d);
  marg[0] = W1;
  if (d >= 2) marg[1] = W2;
  if (d >= 3) marg[2] = W3;
  return List::create(_["m"] = m, _["S2"] = S2, _["marg"] = marg, _["loglik"] = loglik);
}

// Item scores for one dimension: by the Fisher identity, the score of the
// marginal likelihood is sum_k w_ik d log P(y_ij | G(k, i)) / d psi, where w is
// the marginal posterior on that dimension's nodes G (K x N). Rows per item:
// a_2..a_H, b.
// [[Rcpp::export]]
NumericMatrix item_scores_md_cpp(const IntegerMatrix& Y, const NumericMatrix& alpha,
                                 const NumericVector& beta, const NumericMatrix& G,
                                 const NumericMatrix& w, const IntegerVector& H) {
  const int N = Y.nrow(), J = Y.ncol(), K = G.nrow();
  int n_rows = 0;
  std::vector<int> start(J);
  for (int j = 0; j < J; ++j) { start[j] = n_rows; n_rows += H[j]; }
  NumericMatrix s(n_rows, N);
  for (int j = 0; j < J; ++j) {
    const int nA = H[j] - 1, rb = start[j] + nA;
    for (int i = 0; i < N; ++i) {
      const int y = Y(i, j);
      if (y == NA_INTEGER) continue;
      double sl = 0.0, su = 0.0, sb = 0.0;
      for (int k = 0; k < K; ++k) {
        if (w(k, i) == 0.0) continue;
        const double t = G(k, i);
        const double el = alpha(j, y - 1) + beta[j] * t, eu = alpha(j, y) + beta[j] * t;
        const double P = plogis_md(el) - plogis_md(eu);
        const double dl = (y == 1) ? 0.0 : dlogis_md(el);
        const double du = (y == H[j]) ? 0.0 : -dlogis_md(eu);
        const double c = w(k, i) / P;
        sl += c * dl; su += c * du; sb += c * t * (dl + du);
      }
      if (y >= 2) s(start[j] + y - 2, i) += sl;
      if (y <= nA) s(start[j] + y - 1, i) += su;
      s(rb, i) = sb;
    }
  }
  return s;
}
