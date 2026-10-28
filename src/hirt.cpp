// C++ kernels for hIRT.
//
// Binary items are handled as graded items with two categories, so the same
// kernels serve hltm/hltm2 and hgrm/hgrm2/hgrmDIF. Throughout,
//   Y      : N x J integer matrix of responses coded 1..H_j (NA allowed)
//   alpha  : J x (Hmax + 1) matrix; row j is c(Inf, a_j2, ..., a_jH, -Inf, <pad>)
//   beta   : J-vector of discrimination parameters
//   off0, off1 : N x J offsets for DIF items (0 x 0 when there is no DIF); the
//            linear predictor is alpha + beta * theta + off0 + theta * off1
// and P(Y_ij = h | theta) = F(alpha[h] + util) - F(alpha[h + 1] + util), with
// F the logistic CDF and alpha indexed from 1 as in the R code.

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

static inline double plogis_(double x) {
  if (x == R_PosInf) return 1.0;
  if (x == R_NegInf) return 0.0;
  return 1.0 / (1.0 + std::exp(-x));
}

// density of the logistic distribution, F(1 - F); zero at +/- Inf
static inline double dlogis_(double x) {
  if (!R_FINITE(x)) return 0.0;
  double p = plogis_(x);
  return p * (1.0 - p);
}

// N x K matrix of sum_j log P(y_ij | theta_k) over observed items
// [[Rcpp::export]]
NumericMatrix loglik_nk_cpp(const IntegerMatrix& Y, const NumericMatrix& alpha,
                            const NumericVector& beta, const NumericVector& theta,
                            const NumericMatrix& off0, const NumericMatrix& off1,
                            const IntegerVector& dif) {
  const int N = Y.nrow(), J = Y.ncol(), K = theta.size();
  const int W = alpha.ncol();
  const bool has_off = off0.nrow() > 0;
  NumericMatrix LL(N, K);
  std::vector<double> logtab(K * W);

  for (int j = 0; j < J; ++j) {
    const bool dif_j = has_off && dif[j];
    if (!dif_j) {
      // no person-specific offsets: tabulate log P(h | theta_k) once per item
      for (int k = 0; k < K; ++k) {
        double u = beta[j] * theta[k];
        for (int h = 1; h < W; ++h) {
          double lo = alpha(j, h - 1), hi = alpha(j, h);
          if (ISNAN(lo) || ISNAN(hi)) { logtab[k * W + h] = NA_REAL; continue; }
          logtab[k * W + h] = std::log(plogis_(lo + u) - plogis_(hi + u));
        }
      }
      for (int k = 0; k < K; ++k) {
        double* LLk = &LL(0, k);
        const double* lt = &logtab[k * W];
        for (int i = 0; i < N; ++i) {
          int y = Y(i, j);
          if (y == NA_INTEGER) continue;
          LLk[i] += lt[y];
        }
      }
    } else {
      for (int k = 0; k < K; ++k) {
        double* LLk = &LL(0, k);
        for (int i = 0; i < N; ++i) {
          int y = Y(i, j);
          if (y == NA_INTEGER) continue;
          double u = beta[j] * theta[k] + off0(i, j) + theta[k] * off1(i, j);
          LLk[i] += std::log(plogis_(alpha(j, y - 1) + u) - plogis_(alpha(j, y) + u));
        }
      }
    }
  }
  return LL;
}

// Posterior weights on the quadrature grid. Returns w (K x N, columns sum to
// one) and the log marginal likelihood. Uses log-sum-exp, so it does not
// underflow when the number of items is large.
// [[Rcpp::export]]
List estep_cpp(const NumericMatrix& LL, const NumericVector& fitted_mean,
               const NumericVector& fitted_var, const NumericVector& theta,
               const NumericVector& qw) {
  const int N = LL.nrow(), K = LL.ncol();
  NumericMatrix w(K, N);
  std::vector<double> lw(K), logqw(K);
  for (int k = 0; k < K; ++k) logqw[k] = std::log(qw[k]);
  double loglik = 0.0;
  for (int i = 0; i < N; ++i) {
    double sd = std::sqrt(fitted_var[i]);
    double mx = R_NegInf;
    for (int k = 0; k < K; ++k) {
      lw[k] = LL(i, k) + R::dnorm(theta[k] - fitted_mean[i], 0.0, sd, 1) + logqw[k];
      if (lw[k] > mx) mx = lw[k];
    }
    double s = 0.0;
    for (int k = 0; k < K; ++k) { lw[k] = std::exp(lw[k] - mx); s += lw[k]; }
    for (int k = 0; k < K; ++k) w(k, i) = lw[k] / s;
    loglik += mx + std::log(s);
  }
  return List::create(_["w"] = w, _["loglik"] = loglik);
}

// Pseudo counts for item j: K x H_j matrix, tab(k, h) = sum_i w(k, i) 1(y_ij = h)
// [[Rcpp::export]]
List pseudo_tab_cpp(const IntegerMatrix& Y, const NumericMatrix& w, const IntegerVector& H) {
  const int N = Y.nrow(), J = Y.ncol(), K = w.nrow();
  List out(J);
  for (int j = 0; j < J; ++j) {
    NumericMatrix tab(K, H[j]);
    for (int i = 0; i < N; ++i) {
      int y = Y(i, j);
      if (y == NA_INTEGER) continue;
      const double* wi = &w(0, i);
      double* col = &tab(0, y - 1);
      for (int k = 0; k < K; ++k) col[k] += wi[k];
    }
    out[j] = tab;
  }
  return out;
}

// Weighted cumulative logit (proportional odds) model by Newton-Raphson:
//   logit P(y >= h) = a_h + x'b,  h = 2..H,  a_2 > ... > a_H.
// y coded 1..H; rows with weight <= tol or missing y are dropped (as lrm_fit).
// Returns c(a_2, ..., a_H, b). The log likelihood is concave, so Newton with
// step halving converges to the same MLE as rms::lrm.fit.
static double cumlogit_ll(const arma::vec& par, const IntegerVector& y, const arma::mat& X,
                          const std::vector<int>& rows, const NumericVector& wt, int H,
                          arma::vec* grad, arma::mat* hess) {
  const int m = X.n_cols, P = H - 1 + m;
  double ll = 0.0;
  if (grad) grad->zeros(P);
  if (hess) hess->zeros(P, P);
  for (int r : rows) {
    int yr = y[r];
    double xb = 0.0;
    for (int c = 0; c < m; ++c) xb += X(r, c) * par[H - 1 + c];
    // indices of a_{y} and a_{y+1} in par (a_h sits at h - 2); -1 if infinite
    int il = (yr >= 2) ? yr - 2 : -1;
    int iu = (yr <= H - 1) ? yr - 1 : -1;
    double el = (il >= 0) ? par[il] + xb : R_PosInf;
    double eu = (iu >= 0) ? par[iu] + xb : R_NegInf;
    double Fl = plogis_(el), Fu = plogis_(eu);
    double Pr = Fl - Fu;
    if (!(Pr > 0)) return R_NegInf;
    double w = wt[r];
    ll += w * std::log(Pr);
    if (!grad) continue;
    double gl = Fl * (1 - Fl), gu = Fu * (1 - Fu);            // f(eta)
    double hl = gl * (1 - 2 * Fl), hu = gu * (1 - 2 * Fu);    // f'(eta)
    double gd = (gl - gu) / Pr;
    if (il >= 0) (*grad)[il] += w * gl / Pr;
    if (iu >= 0) (*grad)[iu] -= w * gu / Pr;
    for (int c = 0; c < m; ++c) (*grad)[H - 1 + c] += w * X(r, c) * gd;
    // Hessian
    double Pr2 = Pr * Pr;
    if (il >= 0) (*hess)(il, il) += w * (hl / Pr - gl * gl / Pr2);
    if (iu >= 0) (*hess)(iu, iu) += w * (-hu / Pr - gu * gu / Pr2);
    if (il >= 0 && iu >= 0) {
      double v = w * gl * gu / Pr2;
      (*hess)(il, iu) += v; (*hess)(iu, il) += v;
    }
    double hbb = (hl - hu) / Pr - gd * gd;
    double hlb = hl / Pr - gl * gd / Pr;
    double hub = -hu / Pr + gu * gd / Pr;
    for (int c = 0; c < m; ++c) {
      double xc = X(r, c);
      if (il >= 0) { (*hess)(il, H - 1 + c) += w * xc * hlb; (*hess)(H - 1 + c, il) += w * xc * hlb; }
      if (iu >= 0) { (*hess)(iu, H - 1 + c) += w * xc * hub; (*hess)(H - 1 + c, iu) += w * xc * hub; }
      for (int d = 0; d <= c; ++d) {
        double v = w * xc * X(r, d) * hbb;
        (*hess)(H - 1 + c, H - 1 + d) += v;
        if (d != c) (*hess)(H - 1 + d, H - 1 + c) += v;
      }
    }
  }
  return ll;
}

// [[Rcpp::export]]
NumericVector cumlogit_fit_cpp(const IntegerVector& y, const NumericMatrix& Xr,
                               const NumericVector& wt, int H, NumericVector init,
                               double tol = 1e-16, int maxit = 100, double eps = 1e-10) {
  arma::mat X(const_cast<double*>(Xr.begin()), Xr.nrow(), Xr.ncol(), false, true);
  std::vector<int> rows;
  rows.reserve(y.size());
  for (int r = 0; r < y.size(); ++r)
    if (y[r] != NA_INTEGER && wt[r] > tol) rows.push_back(r);
  const int P = H - 1 + Xr.ncol();
  arma::vec par(init.begin(), P);
  arma::vec g; arma::mat Hs;
  double ll = cumlogit_ll(par, y, X, rows, wt, H, &g, &Hs);
  if (!R_FINITE(ll)) {
    // fall back to evenly spaced intercepts and zero slopes
    for (int h = 0; h < H - 1; ++h) par[h] = 1.0 - 2.0 * h / std::max(H - 2, 1);
    for (int c = H - 1; c < P; ++c) par[c] = 0.0;
    ll = cumlogit_ll(par, y, X, rows, wt, H, &g, &Hs);
  }
  for (int it = 0; it < maxit; ++it) {
    arma::vec step;
    bool ok = arma::solve(step, -Hs, g, arma::solve_opts::no_approx);
    if (!ok) step = g / (arma::abs(Hs.diag()).max() + 1.0);
    double t = 1.0, ll_new = R_NegInf;
    arma::vec cand;
    for (int half = 0; half < 30; ++half) {
      cand = par + t * step;
      ll_new = cumlogit_ll(cand, y, X, rows, wt, H, nullptr, nullptr);
      if (R_FINITE(ll_new) && ll_new >= ll - 1e-12 * std::abs(ll)) break;
      t *= 0.5;
    }
    if (!R_FINITE(ll_new)) break;
    double change = arma::abs(cand - par).max();
    par = cand;
    ll = cumlogit_ll(par, y, X, rows, wt, H, &g, &Hs);
    if (change < eps) break;
  }
  return NumericVector(par.begin(), par.end());
}

// Scores for the outer-product-of-gradients covariance matrix.
// Returns
//   s_item : (sum_j n_par_j) x N, for each item rows (a_j2..a_jH, eta_j, b_j)
//   s_gamma: p x N, s_lambda: q x N, loglik, and Li (scaled marginal likelihood)
// X0: N x p0 DIF covariates (0 columns if none); dif: J-vector of 0/1 flags;
// nonuniform: whether DIF also enters the slope (adds theta * x0 columns).
// [[Rcpp::export]]
List scores_cpp(const IntegerMatrix& Y, const NumericMatrix& alpha, const NumericVector& beta,
                const NumericVector& theta, const NumericVector& qw,
                const NumericVector& fitted_mean, const NumericVector& fitted_var,
                const NumericMatrix& x, const NumericMatrix& z,
                const NumericMatrix& off0, const NumericMatrix& off1,
                const NumericMatrix& X0, const IntegerVector& dif, bool nonuniform,
                const IntegerVector& H, bool items = true) {
  const int N = Y.nrow(), J = Y.ncol(), K = theta.size();
  const int p = x.ncol(), q = z.ncol(), p0 = X0.ncol();
  NumericMatrix LL = loglik_nk_cpp(Y, alpha, beta, theta, off0, off1, dif);

  // pik * Lik, scaled per person by exp(-max_k LL) for numerical stability
  NumericMatrix PL(N, K);
  NumericVector Li(N);
  double loglik = 0.0;
  for (int i = 0; i < N; ++i) {
    double sd = std::sqrt(fitted_var[i]);
    double mx = R_NegInf;
    for (int k = 0; k < K; ++k) if (LL(i, k) > mx) mx = LL(i, k);
    double s = 0.0;
    for (int k = 0; k < K; ++k) {
      PL(i, k) = R::dnorm(theta[k], fitted_mean[i], sd, 0) * qw[k] * std::exp(LL(i, k) - mx);
      s += PL(i, k);
    }
    Li[i] = s;
    loglik += mx + std::log(s);
  }

  // gamma and lambda scores
  NumericMatrix s_gamma(p, N), s_lambda(q, N);
  for (int i = 0; i < N; ++i) {
    double a = 0.0, b = 0.0, v = fitted_var[i];
    for (int k = 0; k < K; ++k) {
      double d = theta[k] - fitted_mean[i];
      a += PL(i, k) * d;
      b += 0.5 * PL(i, k) * (d * d / v - 1.0);
    }
    a /= v * Li[i]; b /= Li[i];
    for (int c = 0; c < p; ++c) s_gamma(c, i) = a * x(i, c);
    for (int c = 0; c < q; ++c) s_lambda(c, i) = b * z(i, c);
  }
  if (!items) return List::create(_["s_gamma"] = s_gamma, _["s_lambda"] = s_lambda,
                                  _["loglik"] = loglik, _["Li"] = Li);

  // item scores
  int n_rows = 0;
  std::vector<int> start(J);
  const int n_eta = nonuniform ? 2 * p0 : p0;
  for (int j = 0; j < J; ++j) {
    start[j] = n_rows;
    n_rows += H[j] - 1 + (dif[j] ? n_eta : 0) + 1;
  }
  NumericMatrix s_item(n_rows, N);
  const bool has_off = off0.nrow() > 0;
  std::vector<double> Ftab, Dtab;
  for (int j = 0; j < J; ++j) {
    const int nA = H[j] - 1;
    const int rb = start[j] + nA + (dif[j] ? n_eta : 0);   // row of b_j
    const bool dif_j = has_off && dif[j];
    // without person-specific offsets, F and f at each (theta_k, alpha_jm) are
    // shared by all respondents: tabulate them once per item
    const int Wj = H[j] + 1;
    if (!dif_j) {
      Ftab.assign(K * Wj, 0.0); Dtab.assign(K * Wj, 0.0);
      for (int k = 0; k < K; ++k)
        for (int m = 0; m < Wj; ++m) {
          double e = alpha(j, m) + beta[j] * theta[k];
          Ftab[k * Wj + m] = plogis_(e);
          Dtab[k * Wj + m] = dlogis_(e);
        }
    }
    for (int i = 0; i < N; ++i) {
      int y = Y(i, j);
      if (y == NA_INTEGER) continue;
      double sa_l = 0.0, sa_u = 0.0, sb = 0.0, se0 = 0.0, se1 = 0.0;
      double o0 = dif_j ? off0(i, j) : 0.0, o1 = dif_j ? off1(i, j) : 0.0;
      for (int k = 0; k < K; ++k) {
        double Pr, dl, du;
        if (!dif_j) {
          const double* Fk = &Ftab[k * Wj];
          const double* Dk = &Dtab[k * Wj];
          Pr = Fk[y - 1] - Fk[y];
          dl = (y == 1) ? 0.0 : Dk[y - 1];
          du = (y == H[j]) ? 0.0 : -Dk[y];
        } else {
          double u = beta[j] * theta[k] + o0 + theta[k] * o1;
          double el = alpha(j, y - 1) + u, eu = alpha(j, y) + u;
          Pr = plogis_(el) - plogis_(eu);
          dl = (y == 1) ? 0.0 : dlogis_(el);
          du = (y == H[j]) ? 0.0 : -dlogis_(eu);
        }
        double comp = PL(i, k) / Pr;
        sa_l += comp * dl;
        sa_u += comp * du;
        sb += comp * theta[k] * (dl + du);
        se0 += comp * (dl + du);
        se1 += comp * theta[k] * (dl + du);
      }
      // a_{y} sits in row start + y - 2, a_{y+1} in row start + y - 1
      if (y >= 2) s_item(start[j] + y - 2, i) += sa_l / Li[i];
      if (y <= nA) s_item(start[j] + y - 1, i) += sa_u / Li[i];
      s_item(rb, i) = sb / Li[i];
      if (dif[j]) {
        for (int c = 0; c < p0; ++c) {
          s_item(start[j] + nA + c, i) = X0(i, c) * se0 / Li[i];
          if (nonuniform) s_item(start[j] + nA + p0 + c, i) = X0(i, c) * se1 / Li[i];
        }
      }
    }
  }
  return List::create(_["s_item"] = s_item, _["s_gamma"] = s_gamma, _["s_lambda"] = s_lambda,
                      _["loglik"] = loglik, _["Li"] = Li);
}

// Gamma GLM with log link by IRLS, mirroring glm.fit(z, y, family = Gamma("log"))
// (starts from mu = y; working weights are all one under the log link; stops
// when the relative change in deviance is below epsilon, as glm.control).
// [[Rcpp::export]]
List gamma_glm_cpp(const arma::mat& Z, const arma::vec& y, double epsilon = 1e-8,
                   int maxit = 25) {
  arma::vec eta = arma::log(y), mu = y, coef;
  auto dev_fun = [&](const arma::vec& m) {
    return -2.0 * arma::accu(arma::log(y / m) - (y - m) / m);
  };
  double devold = dev_fun(mu), dev = devold;
  for (int it = 0; it < maxit; ++it) {
    arma::vec zw = eta + (y - mu) / mu;
    coef = arma::solve(Z, zw);
    eta = Z * coef;
    mu = arma::exp(eta);
    dev = dev_fun(mu);
    if (std::abs(dev - devold) / (std::abs(dev) + 0.1) < epsilon) break;
    devold = dev;
  }
  return List::create(_["coefficients"] = coef, _["fitted.values"] = mu);
}

// weighted least squares coefficients, as lm.wfit(x, y, w)$coefficients
// [[Rcpp::export]]
arma::vec wls_cpp(const arma::mat& X, const arma::vec& y, const arma::vec& w) {
  arma::vec sw = arma::sqrt(w);
  return arma::solve(X.each_col() % sw, y % sw);
}
