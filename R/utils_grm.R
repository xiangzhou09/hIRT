# check if a vector has at least two valid responses
invalid_grm <- function(x) max(x, na.rm = TRUE) < 2

lrm_fit <- function(x, y, weights, tol = 1e-16, ...){
  valid <- weights>tol & !is.na(y)
  lrm.fit(x[valid, , drop = FALSE], y[valid], weights = weights[valid], ...)
}

# log likelihood function (return N * J matrix) y: N*J data frame alpha:
# length J list beta: length J numeric vector theta: length N numeric
# vector
loglik_grm <- function(alpha, beta, theta) {
    util <- outer(theta, beta)
    alpha_l <- simplify2array(unname(Map(function(x, y) x[y], alpha, y)))
    alpha_h <- simplify2array(unname(Map(function(x, y) x[y + 1L], alpha, y)))
    log(plogis(util + alpha_l) - plogis(util + alpha_h))
}

# posterior of theta (unnormalized) (returns N-vector)
# y: N*J data frame
# x: N*p model matrix
# z: N*q model matrix
# alpha: length J list
# beta: length J numeric vector
# gamma: p-vector
# lambda: q-vector
# theta_k: numeric scalar
# qw_k numeric scalar
theta_post_grm <- function(theta_k, qw_k) {
    wt_k <- dnorm(theta_k - fitted_mean, sd = sqrt(fitted_var)) * qw_k  # prior density * quadrature weight
    loglik <- rowSums(loglik_grm(alpha, beta, rep(theta_k, N)), na.rm = TRUE)
    logPop <- log(wt_k)
    exp(loglik + logPop)
}

theta_prior_grm <- function(theta_k, qw_k) {
  wt_k <- dnorm(theta_k - fitted_mean, sd = sqrt(fitted_var)) * qw_k  # prior density * quadrature weight
  # loglik <- rowSums(loglik_grm(alpha, beta, rep(theta_k, N)), na.rm = TRUE)
  logPop <- log(wt_k)
  exp(logPop)
}

# pseudo tabulated data for item J (returns K*H_j matrix)
# y_j: N-vector
# H_j: number of response categories for item j
# w: K*N matrix
dummy_fun_grm <- function(y_j, H_j) {
    dummy_mat <- outer(y_j, 1:H_j, "==")  # N*H_j matrix
    dummy_mat[is.na(dummy_mat)] <- 0
    w %*% dummy_mat
}

# pseudo tabulated data to pseudo data frame
# tab: K*H_j matrix
# theta_ls: K-vector
tab2df_grm <- function(tab, theta_ls) {
    H_j <- ncol(tab)
    theta <- rep(theta_ls, H_j)
    y <- rep(1:H_j, each = K)
    data.frame(y = factor(y), x = theta, wt = as.double(tab))
}

# score function of alpha and beta (returns an H_j*N matrix) Lik: N*K
# matrix pik: N*K matrix alpha: J-list beta: J-vector theta_ls: K-vector
sj_ab_grm <- function(j) {
    temp2 <- array(0, c(N, K, H[[j]] + 1))
    h <- .subset2(y, j)
    drv_h <- vapply(theta_ls, function(theta_k) exp(alpha[[j]][h] + beta[[j]] *
        theta_k)/(1 + exp(alpha[[j]][h] + beta[[j]] * theta_k))^2, double(N))
    drv_h_plus_one <- -vapply(theta_ls, function(theta_k) exp(alpha[[j]][h +
        1L] + beta[[j]] * theta_k)/(1 + exp(alpha[[j]][h + 1L] + beta[[j]] *
        theta_k))^2, double(N))
    drv_h[h == 1, ] <- 0
    drv_h_plus_one[h == H[[j]], ] <- 0
    for (i in seq_len(N)) {
        if (is.na(h[[i]])) next
        temp2[i, , h[[i]]] <- drv_h[i, ]
        temp2[i, , h[[i]] + 1L] <- drv_h_plus_one[i, ]
    }
    comp_a <- pik * Lik/vapply(Lijk, `[`, 1:N, j, FUN.VALUE = double(N))  # N*K matrix
    comp_a[is.na(comp_a)] <- 0
    s_alpha <- vapply(1:N, function(i) comp_a[i, ] %*% temp2[i, , 2:H[[j]]],
        double(H[[j]] - 1L))  # (H[j]-1)*N matrix
    temp2_beta <- drv_h + drv_h_plus_one
    s_beta <- rowSums(comp_a * matrix(theta_ls, N, K, byrow = TRUE) * temp2_beta)  # N-vector
    s <- sweep(rbind(s_alpha, s_beta), 2, rowSums(Lik * pik), FUN = "/")
}


x0DIF <- function(theta_k){
  if(form_dif == "uniform") x0 else cbind(x0, theta_k * x0)
}

# log likelihood function (return N * J matrix) y for DIF, returns N*J data frame
# alpha: length J list
# beta: length J numeric vector
# theta: length N numeric vector
loglik_grmDIF <- function(alpha, beta, theta, eta) {
  util <- outer(theta, beta)
  tmp_x <- x0DIF(theta)
  util2 <- vapply(1:length(eta), function(j) if(j %in% items_dif) tmp_x %*% eta[[j]] else double(N), double(N))
  alpha_l <- simplify2array(unname(Map(function(x, y) x[y], alpha, y)))
  alpha_h <- simplify2array(unname(Map(function(x, y) x[y + 1L], alpha, y)))
  log(plogis(util + util2 + alpha_l) - plogis(util + util2 + alpha_h))
}

# posterior of theta (unnormalized) (returns N-vector)
# y: N*J data frame
# x: N*p model matrix
# z: N*q model matrix
# alpha: length J list
# beta: length J numeric vector
# gamma: p-vector
# lambda: q-vector
# theta_k: numeric scalar
# qw_k numeric scalar
theta_post_grmDIF <- function(theta_k, qw_k) {
  wt_k <- dnorm(theta_k - fitted_mean, sd = sqrt(fitted_var)) * qw_k  # prior density * quadrature weight
  loglik <- rowSums(loglik_grmDIF(alpha, beta, rep(theta_k, N), eta), na.rm = TRUE)
  logPop <- log(wt_k)
  exp(loglik + logPop)
}

# score function of alpha and beta and eta (returns an ?*N matrix)
# Lik: N*K matrix pik: N*K matrix alpha: J-list beta: J-vector theta_ls: K-vector
sj_ab_grmDIF <- function(j) {

  temp2 <- array(0, c(N, K, H[[j]] + 1))
  h <- .subset2(y, j)

  if(j %in% items_dif){

    drv_h <- vapply(theta_ls, function(theta_k) exp(alpha[[j]][h] + beta[[j]] * theta_k +
                                                      x0DIF(theta_k) %*% eta[[j]])/
                        (1 + exp(alpha[[j]][h] + beta[[j]] * theta_k + x0DIF(theta_k) %*% eta[[j]]))^2,
                    double(N))

    drv_h_plus_one <- -vapply(theta_ls, function(theta_k) exp(alpha[[j]][h + 1L] + beta[[j]] * theta_k +
                                                                x0DIF(theta_k) %*% eta[[j]])/
                                (1 + exp(alpha[[j]][h + 1L] + beta[[j]] * theta_k + x0DIF(theta_k) %*% eta[[j]]))^2,
                              double(N))
    } else{

      drv_h <- vapply(theta_ls, function(theta_k) exp(alpha[[j]][h] + beta[[j]] * theta_k)/
                        (1 + exp(alpha[[j]][h] + beta[[j]] * theta_k))^2,
                      double(N))

      drv_h_plus_one <- -vapply(theta_ls, function(theta_k) exp(alpha[[j]][h + 1L] + beta[[j]] * theta_k)/
                                  (1 + exp(alpha[[j]][h + 1L] + beta[[j]] * theta_k))^2,
                                double(N))
    }
  drv_h[h == 1, ] <- 0
  drv_h_plus_one[h == H[[j]], ] <- 0

  for (i in seq_len(N)) {
    if (is.na(h[[i]])) next
    temp2[i, , h[[i]]] <- drv_h[i, ]
    temp2[i, , h[[i]] + 1L] <- drv_h_plus_one[i, ]
  }

  comp_a <- pik * Lik/vapply(Lijk, `[`, 1:N, j, FUN.VALUE = double(N))  # N*K matrix
  comp_a[is.na(comp_a)] <- 0

  s_alpha <- vapply(1:N, function(i) comp_a[i, ] %*% temp2[i, , 2:H[[j]]],
                    double(H[[j]] - 1L))  # (H[j]-1)*N matrix

  temp2_beta <- drv_h + drv_h_plus_one
  s_beta <- rowSums(comp_a * matrix(theta_ls, N, K, byrow = TRUE) * temp2_beta)  # N-vector

  if(j %in% items_dif){
    s_eta <- vapply(1:p0, function(i) rowSums(comp_a * matrix(x0[, i, drop = FALSE], N, K, byrow = FALSE) * temp2_beta), double(N))
    if(form_dif == "non-uniform"){
      s_eta <- cbind(s_eta, vapply(1:p0, function(i) rowSums(comp_a * outer(x0[, i], theta_ls) * temp2_beta), double(N)))
    }
    s <- sweep(rbind(s_alpha, t(s_eta), s_beta), 2, rowSums(Lik * pik), FUN = "/")
  } else{
    s <- sweep(rbind(s_alpha, s_beta), 2, rowSums(Lik * pik), FUN = "/")
  }
}

