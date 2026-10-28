#' @useDynLib hIRT, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

# Wrappers around the C++ kernels in src/hirt.cpp. Binary items are passed as
# graded items with two categories (responses coded 1, 2).

# alpha list (each element c(Inf, a_2, ..., a_H, -Inf)) -> J x (Hmax + 1) matrix
alpha_mat <- function(alpha) {
  W <- max(lengths(alpha))
  t(vapply(alpha, function(a) c(unname(a), rep(NA_real_, W - length(a))), double(W)))
}

# integer response matrix coded 1..H_j
y_mat <- function(y) {
  m <- as.matrix(y)
  storage.mode(m) <- "integer"
  unname(m)
}

dbl_mat <- function(m) {
  storage.mode(m) <- "double"
  m
}

no_off <- matrix(0, 0, 0)

# E-step: posterior weights (K x N) and log marginal likelihood
estep_grm <- function(Y, alpha, beta, theta_ls, qw_ls, fitted_mean, fitted_var,
                      off0 = no_off, off1 = no_off, dif = integer(ncol(Y))) {
  LL <- loglik_nk_cpp(Y, alpha_mat(alpha), as.double(beta), theta_ls, off0, off1,
                      as.integer(dif))
  estep_cpp(LL, as.double(fitted_mean), as.double(fitted_var), theta_ls, qw_ls)
}

# weighted cumulative logit fit with named coefficients, as from rms::lrm.fit
cumlogit_fit <- function(y, X, weights, H, init) {
  X <- dbl_mat(as.matrix(X))
  coefs <- cumlogit_fit_cpp(as.integer(y), X, as.double(weights), as.integer(H),
                            as.double(init))
  # rms::lrm.fit calls the single intercept of a binary response "Intercept"
  int_names <- if (H == 2) "Intercept" else paste0("y>=", seq(2, H))
  setNames(coefs, c(int_names, colnames(X)))
}

# item M-step: one grouped-data cumulative logit per item, warm-started
mstep_grm <- function(Y, w, H, theta_ls, alpha, beta) {
  K <- length(theta_ls)
  tabs <- pseudo_tab_cpp(Y, w, as.integer(H))
  fits <- lapply(seq_along(H), function(j) {
    Hj <- H[[j]]
    cumlogit_fit(rep(seq_len(Hj), each = K), cbind(x = rep(theta_ls, Hj)),
                 as.double(tabs[[j]]), Hj, c(alpha[[j]][2:Hj], beta[[j]]))
  })
  beta <- setNames(vapply(fits, function(f) f[[length(f)]], double(1L)), names(H))
  alpha <- setNames(lapply(fits, function(f) c(Inf, f[-length(f)], -Inf)), names(H))
  list(alpha = alpha, beta = beta)
}

# scores for the outer product of gradients
scores_grm <- function(Y, alpha, beta, theta_ls, qw_ls, fitted_mean, fitted_var, x, z, H,
                       items = TRUE, off0 = no_off, off1 = no_off, X0 = matrix(0, nrow(Y), 0),
                       dif = integer(ncol(Y)), nonuniform = FALSE) {
  scores_cpp(Y, alpha_mat(alpha), as.double(beta), theta_ls, qw_ls,
             as.double(fitted_mean), as.double(fitted_var), dbl_mat(x), dbl_mat(z),
             off0, off1, dbl_mat(X0), as.integer(dif), nonuniform, as.integer(H), items)
}

# hltm stores intercepts as a numeric vector; convert to the graded-item layout
ltm_alpha <- function(alpha) lapply(alpha, function(a) c(Inf, a, -Inf))

# Level-II update: conditional maximization of the mean and variance
# regressions given the posterior means and variances (Supplementary Material B)
var_regression <- function(x, z, theta_eap, theta_vap, lm_opr, con) {
  gamma <- lm_opr %*% theta_eap
  r2 <- as.double((theta_eap - x %*% gamma)^2 + theta_vap)
  if (ncol(z) == 1) return(list(gamma = gamma, lambda = log(mean(r2))))
  x <- dbl_mat(x); z <- dbl_mat(z); theta_eap <- as.double(theta_eap)
  s2 <- as.double(gamma_glm_cpp(z, r2)[["fitted.values"]])
  LL0 <- sum(-0.5 * (log(s2) + r2/s2))
  for (m in seq(1, con[["max_iter2"]])) {
    gamma <- setNames(as.double(wls_cpp(x, theta_eap, 1/s2)), colnames(x))
    r2 <- as.double((theta_eap - x %*% gamma)^2 + theta_vap)
    var_reg <- gamma_glm_cpp(z, r2)
    s2 <- as.double(var_reg[["fitted.values"]])
    LL_temp <- sum(-0.5 * (log(s2) + r2/s2))
    if (LL_temp - LL0 < con[["eps2"]]) break
    LL0 <- LL_temp
  }
  list(gamma = gamma, lambda = setNames(as.double(var_reg[["coefficients"]]), colnames(z)))
}
