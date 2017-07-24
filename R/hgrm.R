#' Fitting Hierarchical Graded Response Models
#'
#' \code{hgrm} fits a hierarchical graded response model in which both
#' the mean and the variance of the latent preference (ability parameter)
#' may depend on person-specific covariates (\code{x} and \code{z}).
#' Specifically, the mean is specified as a linear combination of \code{x}
#' and the log of the variance is specified as a linear combination of
#' \code{z}. Nonresponses are treated as missing at random.
#'
#' @param y A data frame or matrix of item responses.
#' @param x A matrix of covariates, including intercept, that predict the
#'   mean of the latent preference.
#' @param z A matrix of covariates, including intercept, that predict the
#'   variance of the latent preference.
#' @param beta_set The index of the item whose discrimination parameter is
#'   restricted to be positive (or negative). It may take an integer value from
#'   1 to \code{ncol(y)}.
#' @param sign_set Logical. Should the discrimination parameter of
#'   the corresponding item (indexed by \code{beta_set}) be positive
#'   (if \code{TRUE}) or negative (if \code{FALSE})?
#' @param control A list of control values
#' \describe{
#'  \item{max_iter}{The maximum number of iterations of the EM algorithm.
#'    Default 150.}
#'  \item{eps}{Tolerance parameter used to determine convergence of the
#'   EM algorithm. Specifically, iterations continue until the Euclidean
#'   distance between \eqn{\beta_{n}} and \eqn{\beta_{n-1}} falls under \code{eps},
#'   where \eqn{\beta} is the vector of item discrimination parameters.
#'   \code{eps}=1e-4 by default.}
#'  \item{max_iter2}{The maximum number of iterations of the conditional
#'    maximization procedures for updating \eqn{\gamma} and \eqn{\lambda}.
#'    Default 15.}
#'  \item{eps2}{Tolerance parameter used to determine convergence of the
#'    conditional maximization procedures for updating \eqn{\gamma} and
#'    \eqn{\lambda}. Specifically, iterations continue until the Euclidean
#'   distance between two consecutive log likelihoods falls under \code{eps2}.
#'   \code{eps2}=1e-3 by default.}
#'  \item{K}{Number of Gauss-Legendre quadrature points for the E-step. Default 21.}
#'  \item{C}{[-C, C] sets the range of integral in the E-step. \code{C}=5 by default.}
#' }
#'
#' @return An object of class \code{hgrm}.
#'  \item{coefficients}{A data frame of parameter estimates, standard errors,
#'   z values and p values.}
#'  \item{scores}{A data frame of EAP estimates of latent preferences and
#'   their approximate standard errors.}
#'  \item{vcov}{Variance-covariance matrix of parameter estimates.}
#'  \item{log_Lik}{The log-likelihood value at convergence.}
#'  \item{H}{A vector denoting the number of response categories for each item.}
#'  \item{p}{The number of predictors for the mean equation.}
#'  \item{q}{The number of predictors for the variance equation.}
#'  \item{item_names}{Names of items.}
#'  \item{call}{The matched call.}
#' @references Zhou, Xiang. 2017. "Hierarchical Item Response Models for Analyzing Public Opinion."
#'   Working Paper.
#' @importFrom rms lrm.fit
#' @importFrom pryr compose
#' @importFrom pryr partial
#' @import stats
#' @export
#' @examples
#' y <- nes_econ2008[, -(1:3)]
#' x <- model.matrix( ~ party * educ, nes_econ2008)
#' z <- model.matrix( ~ party, nes_econ2008)
#' nes_m1 <- hgrm(y, x, z)
#' print(nes_m1)

hgrm <- function(y, x = matrix(1, nrow(y), 1), z = matrix(1, nrow(y), 1),
    beta_set = 1, sign_set = TRUE, control = list()) {

    # match call
    cl <- match.call()

    # check y
    if ((!is.data.frame(y) & !is.matrix(y)) || ncol(y) == 1L)
        stop("'y' must be either a data.frame or a matrix with at least two columns.")

    # check missing columns
    y <- as.data.frame(y)
    N <- nrow(y)
    J <- ncol(y)
    for (j in seq(1, J)) y[[j]] <- fac2int(y[[j]])
    tmp <- match(TRUE, vapply(y, invalid_grm, logical(1L)))
    if (!is.na(tmp))
      stop(paste(names(y)[tmp], "does not have at least two valid responses"))
    if(max(y, na.rm = TRUE)==2)
      stop("All items are dichotomous. Use 'hltm' ")

    # check x and z (x and z should contain an intercept column)
    if (is.null(nrow(x)))
        x <- as.matrix(x)
    if (is.null(nrow(x)))
        z <- as.matrix(z)
    if (nrow(x) != N || nrow(z) != N)
        stop("both 'x' and 'z' must have the same number of rows as 'y'")
    x <- `colnames<-`(model.matrix(~0 + x), colnames(x))
    z <- `colnames<-`(model.matrix(~0 + z), colnames(z))

    # check beta_set and sign_set
    stopifnot(beta_set %in% 1:J, is.logical(sign_set))

    # control parameters
    con <- list(max_iter = 150, max_iter2 = 15, eps = 1e-04, eps2 = 0.001,
        K = 21, C = 5)  # control parameters
    con[names(control)] <- control

    # dimensions, response categories, etc.
    H <- vapply(y, max, numeric(1), na.rm = TRUE)
    p <- ncol(x)
    q <- ncol(z)
    y_names <- names(y)
    x_names <- colnames(x)
    z_names <- colnames(z)
    environment(loglik_grm) <- environment(theta_post_grm) <- environment(dummy_fun_grm) <- environment(tab2df_grm) <- environment()

    # GL points
    K <- con[["K"]]
    theta_ls <- con[["C"]] * GLpoints[[K]][["x"]]
    qw_ls <- con[["C"]] * GLpoints[[K]][["w"]]

    # initialization
    theta_eap <- {
        tmp <- rowMeans(y, na.rm = TRUE)
        (tmp - mean(tmp, na.rm = TRUE))/sd(tmp, na.rm = TRUE)
    }
    theta_eap[is.na(theta_eap)] <- 0
    alpha <- lapply(H, function(x) c(Inf, seq(1, -1, length.out = x - 1),
        -Inf))
    beta <- vapply(y, function(y) cov(y, theta_eap, use = "complete.obs")/var(theta_eap),
        numeric(1L))
    gamma <- solve(t(x) %*% x) %*% t(x) %*% theta_eap
    lambda <- rep(0, q)

    # EM algorithm
    for (iter in seq(1, con[["max_iter"]])) {

        # store previous parameters
        alpha_prev <- alpha
        beta_prev <- beta
        gamma_prev <- gamma
        lambda_prev <- lambda

        # construct w_ik
        posterior <- Map(theta_post_grm, theta_ls, qw_ls)
        w <- {
            tmp <- matrix(unlist(posterior), N, K)
            t(sweep(tmp, 1, rowSums(tmp), FUN = "/"))
        }

        # maximization
        pseudo_tab <- Map(dummy_fun_grm, y, H)
        pseudo_y <- lapply(pseudo_tab, tab2df_grm, theta_ls = theta_ls)
        pseudo_lrm <- lapply(pseudo_y, function(df) lrm.fit(df[["x"]],
            df[["y"]], weights = df[["wt"]])[["coefficients"]])
        beta <- vapply(pseudo_lrm, function(x) x[length(x)], numeric(1L))
        alpha <- lapply(pseudo_lrm, function(x) c(Inf, x[-length(x)], -Inf))
        theta_eap <- t(theta_ls %*% w)
        theta_vap <- t(theta_ls^2 %*% w) - theta_eap^2

        # variance regression
        gamma <- solve(t(x) %*% x) %*% t(x) %*% theta_eap
        r2 <- (theta_eap - x %*% gamma)^2 + theta_vap
        s2 <- glm.fit(x = z, y = r2, intercept = FALSE, family = Gamma(link = "log"))[["fitted.values"]]
        loglik <- -0.5 * (log(s2) + r2/s2)
        LL0 <- sum(loglik)
        dLL <- 1
        for (m in seq(1, con[["max_iter2"]])) {
            gamma <- lm(theta_eap ~ 0 + x, weights = 1/s2)[["coefficients"]]
            r2 <- (theta_eap - x %*% gamma)^2 + theta_vap
            var_reg <- glm.fit(x = z, y = r2, intercept = FALSE, family = Gamma(link = "log"))
            s2 <- var_reg[["fitted.values"]]
            loglik <- -0.5 * (log(s2) + r2/s2)
            LL_temp <- sum(loglik)
            dLL <- LL_temp - LL0
            if (dLL < con[["eps2"]])
                break
            LL0 <- LL_temp
        }
        lambda <- var_reg[["coefficients"]]

        # location constraint
        tmp <- mean(x %*% gamma)
        alpha <- Map(function(x, y) x + tmp * y, alpha, beta)
        gamma[1L] <- gamma[1L] - tmp

        # scale constraint
        tmp <- mean(z %*% lambda)
        gamma <- gamma/exp(tmp/2)
        beta <- beta * exp(tmp/2)
        lambda[1L] <- lambda[1L] - tmp

        # direction contraint
        if (sign_set == (beta[beta_set] < 0)) {
            gamma <- -gamma
            beta <- -beta
        }
        cat(".")

        # check convergence
        if (sqrt(sum((beta - beta_prev)^2)) < con[["eps"]]) {
            cat("\n converged at iteration", iter, "\n")
            gamma <- setNames(gamma, paste("x", colnames(x), sep = "_"))
            lambda <- setNames(lambda, paste("z", colnames(z), sep = "_"))
            break
        } else if (iter == con[["max_iter"]]) {
            stop("algorithm did not converge; try increasing max_iter.")
            break
        } else next
    }

    # inference
    pik <- matrix(unlist(Map(partial(dnorm, x = theta_ls), mean = as.vector(x %*%
        gamma), sd = as.vector(exp(z %*% lambda)))), N, K, byrow = TRUE) *
        matrix(qw_ls, N, K, byrow = TRUE)
    Lijk <- lapply(theta_ls, function(theta_k) exp(loglik_grm(alpha = alpha,
        beta = beta, rep(theta_k, N))))  # K-list
    Lik <- vapply(Lijk, compose(exp, partial(rowSums, na.rm = TRUE),
        log), numeric(N))
    Li <- rowSums(Lik * pik)

    # log likelihood
    log_Lik <- sum(log(Li))

    # outer product of gradients
    environment(sj_ab_grm) <- environment(si_gamma) <- environment(si_lambda) <- environment()
    s_ab <- unname(Reduce(rbind, lapply(1:J, sj_ab_grm)))
    s_lambda <- s_gamma <- NULL
    if (p > 1)
        s_gamma <- vapply(1:N, si_gamma, numeric(p - 1))
    if (q > 1)
        s_lambda <- vapply(1:N, si_lambda, numeric(q - 1))
    s_all <- rbind(s_ab, s_gamma, s_lambda)
    s_all[is.na(s_all)] <- 0
    covmat <- solve(s_all %*% t(s_all))
    se_all <- sqrt(diag(covmat))

    # reorganize se_all
    sH <- sum(H)
    lambda_indices <- gamma_indices <- NULL
    if (p > 1)
        gamma_indices <- (sH + 1):(sH + p - 1)
    if (q > 1)
        lambda_indices <- (sH + p):(sH + p + q - 2)
    se_all <- c(se_all[1:sH], NA, se_all[gamma_indices], NA, se_all[lambda_indices])

    # name se_all and covmat
    names_ab <- unlist(lapply(names(alpha), function(x) {
        tmp <- alpha[[x]]
        paste(x, c(names(tmp)[-c(1, length(tmp))], "Dscrmn"))
    }))
    names(se_all) <- c(names_ab, names(gamma), names(lambda))
    rownames(covmat) <- colnames(covmat) <- c(names_ab, names(gamma)[-1L],
        names(lambda)[-1L])

    # item coefficients
    coef_item <- Map(function(a, b) c(a[-c(1L, length(a))], Dscrmn = b),
        alpha, beta)

    # all coefficients
    coef_all <- c(unlist(coef_item), gamma, lambda)
    coefs <- data.frame(Estimate = coef_all, Std_Error = se_all, z_value = coef_all/se_all,
        p_value = 2 * (1 - pnorm(abs(coef_all/se_all))))
    rownames(coefs) <- names(se_all)

    # ability parameter estimates
    theta <- data.frame(est = theta_eap, se = sqrt(theta_vap))

    # output
    out <- list(coefficients = coefs, scores = theta, vcov = covmat, log_Lik = log_Lik,
        H = H, p = p, q = q, item_names = names(y), call = cl)
    class(out) <- c("hgrm", "hIRT")
    out
}

