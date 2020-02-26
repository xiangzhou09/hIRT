# check if a vector is dichotomous
invalid_ltm <- function(x) max(x, na.rm = TRUE) != 1

glm_fit <- function(x, y, weights, tol = 1e-16, ...){
    glm.fit(x[weights>tol, , drop = FALSE], y[weights>tol], weights = weights[weights>tol], ...)
}

# log likelihood function (return N * J matrix) y: N*J data frame alpha:
# length J numeric vector beta: length J numeric vector theta: length N
# numeric vector
loglik_ltm <- function(alpha, beta, theta) {
    util <- matrix(alpha, N, J, byrow = TRUE) + outer(theta, beta)
    log(exp(as.matrix(y) * util)/(1 + exp(util)))
}

# posterior of theta (unnormalized) (returns N-vector) y: N*J data frame
# x: N*p model matrix z: N*q model matrix alpha: length J list beta:
# length J numeric vector gamma: p-vector lambda: q-vector theta_k:
# numeric scalar qw_k numeric scalar
theta_post_ltm <- function(theta_k, qw_k) {
    N <- nrow(y)
    wt_k <- dnorm(theta_k - fitted_mean, sd = sqrt(fitted_var)) * qw_k  # prior density * quadrature weight
    loglik <- rowSums(loglik_ltm(alpha, beta, rep(theta_k, N)), na.rm = TRUE)
    logPop <- log(wt_k)
    exp(loglik + logPop)
}

# pseudo tabulated data for item J (returns K*2 matrix) y_j: N-vector w:
# K*N matrix
dummy_fun_ltm <- function(y_j) {
    dummy_mat <- outer(y_j, c(0, 1), "==")  # N*H_j matrix
    dummy_mat[is.na(dummy_mat)] <- 0
    w %*% dummy_mat
}

# pseudo tabulated data to pseudo data frame tab: K*2 matrix theta_ls:
# K-vector
tab2df_ltm <- function(tab, theta_ls) {
    theta <- rep(theta_ls, 2)
    y <- rep(c(0, 1), each = K)
    data.frame(y = factor(y), x = theta, wt = as.double(tab))
}

# derivative of likelihood wrt alpha, given theta_k
dalpha_ltm <- function(alpha, beta) {
    putil <- plogis(matrix(alpha, K, J, byrow = TRUE) + outer(theta_ls, beta))
    putil * (1 - putil)
}

# score function of alpha and beta (return a N*2 matrix) Lik: N*K matrix
# pik: N*K matrix alpha: J-vector beta: J-vector theta_ls: K-vector
sj_ab_ltm <- function(j) {
    tmp_mat <- (pik * Lik/vapply(Lijk, `[`, 1:N, j, FUN.VALUE = double(N)))  # N*K matrix
    dalpha_j <- dalpha[, j, drop = FALSE]  # K*1 matrix
    dbeta_j <- dalpha_j * theta_ls  # K*1 matrix
    sgn <- .subset2(y, j) * 2 - 1
    sgn[is.na(sgn)] <- 0  # N-vector (1, -1, or 0)
    drv_alpha <- sgn * (tmp_mat %*% dalpha_j)/Li
    drv_beta <- sgn * (tmp_mat %*% dbeta_j)/Li
    cbind(drv_alpha, drv_beta)
}
