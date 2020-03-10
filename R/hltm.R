#' Fitting Hierarchical Latent Trait Models (for Binary Responses).
#'
#' \code{hltm} fits a hierarchical latent trait model in which both
#' the mean and the variance of the latent preference (ability parameter)
#' may depend on person-specific covariates (\code{x} and \code{z}).
#' Specifically, the mean is specified as a linear combination of \code{x}
#' and the log of the variance is specified as a linear combination of
#' \code{z}.
#'
#' @inheritParams hgrm
#'
#' @return An object of class \code{hltm}.
#'  \item{coefficients}{A data frame of parameter estimates, standard errors,
#'   z values and p values.}
#'  \item{scores}{A data frame of EAP estimates of latent preferences and
#'   their approximate standard errors.}
#'  \item{vcov}{Variance-covariance matrix of parameter estimates.}
#'  \item{log_Lik}{The log-likelihood value at convergence.}
#'  \item{N}{Number of units.}
#'  \item{J}{Number of items.}
#'  \item{H}{A vector denoting the number of response categories for each item.}
#'  \item{ylevels}{A list showing the levels of the factorized response categories.}
#'  \item{p}{The number of predictors for the mean equation.}
#'  \item{q}{The number of predictors for the variance equation.}
#'  \item{control}{List of control values.}
#'  \item{call}{The matched call.}
#' @references Zhou, Xiang. 2019. "\href{https://doi.org/10.1017/pan.2018.63}{Hierarchical Item Response Models for Analyzing Public Opinion.}" Political Analysis.
#' @importFrom rms lrm.fit
#' @importFrom pryr compose
#' @importFrom pryr partial
#' @import stats
#' @export
#' @examples
#' y <- nes_econ2008[, -(1:3)]
#' x <- model.matrix( ~ party * educ, nes_econ2008)
#' z <- model.matrix( ~ party, nes_econ2008)
#'
#' dichotomize <- function(x) findInterval(x, c(mean(x, na.rm = TRUE)))
#' y[] <- lapply(y, dichotomize)
#' nes_m1 <- hltm(y, x, z)
#' nes_m1

hltm <- function(y, x = NULL, z = NULL, constr = c("latent_scale", "items"),
                 beta_set = 1L, sign_set = TRUE, init = c("naive", "glm", "irt"),
                 control = list()) {

  # match call
  cl <- match.call()

  # check y and convert y into data.frame if needed
  if(missing(y)) stop("`y` must be provided.")
  if ((!is.data.frame(y) && !is.matrix(y)) || ncol(y) == 1L)
    stop("'y' must be either a data.frame or a matrix with at least two columns.")
  if(is.matrix(y)) y <- as.data.frame(y)

  # number of units and items
  N <- nrow(y)
  J <- ncol(y)

  # convert each y_j into an integer vector
  y[] <- lapply(y, factor, exclude = c(NA, NaN))
  ylevels <- lapply(y, levels)
  y[] <- lapply(y, function(x) as.integer(x) - 1)
  if (!is.na(invalid <- match(TRUE, vapply(y, invalid_ltm, logical(1L)))))
    stop(paste(names(y)[invalid], "is not a dichotomous variable"))
  H <- vapply(y, max, double(1L), na.rm = TRUE) + 1

  # check x and z (x and z should contain an intercept column)
  x <- x %||% as.matrix(rep(1, N))
  z <- z %||% as.matrix(rep(1, N))
  if (!is.matrix(x)) stop("`x` must be a matrix.")
  if (!is.matrix(z)) stop("`z` must be a matrix.")
  if (nrow(x) != N || nrow(z) != N) stop("both 'x' and 'z' must have the same number of rows as 'y'")
  p <- ncol(x)
  q <- ncol(z)
  colnames(x) <- colnames(x) %||% paste0("x", 1:p)
  colnames(z) <- colnames(z) %||% paste0("x", 1:q)

  # check beta_set and sign_set
  stopifnot(beta_set %in% 1:J, is.logical(sign_set))

  # check constraint
  constr <- match.arg(constr)
  init <- match.arg(init)

  # control parameters
  con <- list(max_iter = 150, max_iter2 = 15, eps = 1e-03, eps2 = 1e-03, K = 25, C = 4)
  con[names(control)] <- control

  # set environments for utility functions
  environment(loglik_ltm) <- environment(theta_post_ltm) <- environment(dummy_fun_ltm) <- environment(tab2df_ltm) <- environment()

  # GL points
  K <- con[["K"]]
  theta_ls <- con[["C"]] * GLpoints[[K]][["x"]]
  qw_ls <- con[["C"]] * GLpoints[[K]][["w"]]

  # imputation
  y_imp <- y
  if(anyNA(y)) y_imp[] <- lapply(y, impute)

  # pca for initial values of theta_eap
  theta_eap <- {
    tmp <- princomp(y_imp, cor = TRUE)$scores[, 1]
    (tmp - mean(tmp, na.rm = TRUE))/sd(tmp, na.rm = TRUE)
  }

  # initialization of alpha and beta parameters
  if (init == "naive"){
    alpha <- rep(0, J)
    beta <- vapply(y, function(y) cov(y, theta_eap, use = "complete.obs")/var(theta_eap), double(1L))
  } else if (init == "glm"){
    pseudo_logit <- lapply(y_imp, function(y) glm.fit(cbind(1, theta_eap), y, family = binomial("logit"))[["coefficients"]])
    beta <- vapply(pseudo_logit, function(x) x[2L], double(1L))
    alpha <- vapply(pseudo_logit, function(x) x[1L], double(1L))
  } else {
    ltm_coefs <- ltm(y ~ z1)[["coefficients"]]
    beta <- ltm_coefs[, 2, drop = TRUE]
    alpha <- ltm_coefs[, 1, drop = TRUE]
  }

  # initial values of gamma and lambda
  lm_opr <- tcrossprod(solve(crossprod(x)), x)
  gamma <- lm_opr %*% theta_eap
  lambda <- rep(0, q)
  fitted_mean <- as.double(x %*% gamma)
  fitted_var <- rep(1, N)

  # EM algorithm
  for (iter in seq(1, con[["max_iter"]])) {

      # store previous parameters
      alpha_prev <- alpha
      beta_prev <- beta
      gamma_prev <- gamma
      lambda_prev <- lambda

      # construct w_ik
      posterior <- Map(theta_post_ltm, theta_ls, qw_ls)
      w <- {
          tmp <- matrix(unlist(posterior), N, K)
          t(sweep(tmp, 1, rowSums(tmp), FUN = "/"))
      }

      # maximization
      pseudo_tab <- lapply(y, dummy_fun_ltm)
      pseudo_y <- lapply(pseudo_tab, tab2df_ltm, theta_ls = theta_ls)
      pseudo_logit <- lapply(pseudo_y, function(df) glm.fit(cbind(1, df[["x"]]),
          df[["y"]], weights = df[["wt"]], family = quasibinomial("logit"))[["coefficients"]])
      beta <- vapply(pseudo_logit, function(x) x[2L], double(1L))
      alpha <- vapply(pseudo_logit, function(x) x[1L], double(1L))

      # EAP and VAP estimates of latent preferences
      theta_eap <- t(theta_ls %*% w)
      theta_vap <- t(theta_ls^2 %*% w) - theta_eap^2

      # variance regression
      gamma <- lm_opr %*% theta_eap
      r2 <- (theta_eap - x %*% gamma)^2 + theta_vap

      if (ncol(z)==1) lambda <- log(mean(r2)) else{
        s2 <- glm.fit(x = z, y = r2, intercept = FALSE, family = Gamma(link = "log"))[["fitted.values"]]
        loglik <- -0.5 * (log(s2) + r2/s2)
        LL0 <- sum(loglik)
        dLL <- 1
        for (m in seq(1, con[["max_iter2"]])) {
            gamma <- lm.wfit(x, theta_eap, w = 1/s2)[["coefficients"]]
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
      }

      # location constraint
      tmp <- mean(x %*% gamma)
      alpha <- unlist(Map(function(x, y) x + tmp * y, alpha, beta))
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
      fitted_mean <- as.double(x %*% gamma)
      fitted_var <- exp(as.double(z %*% lambda))

      cat(".")

      # check convergence
      if (sqrt(mean((beta - beta_prev)^2)) < con[["eps"]]) {
          cat("\n converged at iteration", iter, "\n")
          break
      } else if (iter == con[["max_iter"]]) {
          stop("algorithm did not converge; try increasing max_iter.")
          break
      } else next
  }

  gamma <- setNames(as.double(gamma), paste("x", colnames(x), sep = ""))
  lambda <- setNames(as.double(lambda), paste("z", colnames(z), sep = ""))

  # inference
  pik <- matrix(unlist(Map(partial(dnorm, x = theta_ls), mean = fitted_mean, sd = sqrt(fitted_var))),
                N, K, byrow = TRUE) * matrix(qw_ls, N, K, byrow = TRUE)
  Lijk <- lapply(theta_ls, function(theta_k) exp(loglik_ltm(alpha = alpha, beta = beta, rep(theta_k, N))))  # K-list
  Lik <- vapply(Lijk, compose(exp, partial(rowSums, na.rm = TRUE), log), double(N))
  Li <- rowSums(Lik * pik)

  # log likelihood
  log_Lik <- sum(log(Li))

  # outer product of gradients
  environment(dalpha_ltm) <- environment(sj_ab_ltm) <- environment(si_gamma) <- environment(si_lambda) <- environment()
  dalpha <- dalpha_ltm(alpha, beta)  # K*J matrix
  s_ab <- unname(Reduce(cbind, lapply(1:J, sj_ab_ltm)))
  s_gamma <- vapply(1:N, si_gamma, double(p))
  s_lambda <- vapply(1:N, si_lambda, double(q))

  s_all <- rbind(t(s_ab)[-c(1L, ncol(s_ab)), , drop = FALSE], s_gamma, s_lambda)
  s_all[is.na(s_all)] <- 0
  covmat <- tryCatch(solve(tcrossprod(s_all)),
                     error = function(e) {warning("The information matrix is singular; SE calculation failed.");
                       matrix(NA, nrow(s_all), nrow(s_all))})
  se_all <- sqrt(diag(covmat))

  # reorganize se_all
  sH <- 2 * J
  gamma_indices <- (sH - 1):(sH + p - 2)
  lambda_indices <- (sH + p - 1):(sH + p + q - 2)
  se_all <- c(NA, se_all[1:(sH-2)], NA, se_all[gamma_indices], se_all[lambda_indices])

  # name se_all and covmat
  names_ab <- paste(rep(names(alpha), each = 2), c("Diff", "Dscrmn"))
  names(se_all) <- c(names_ab, names(gamma), names(lambda))
  rownames(covmat) <- colnames(covmat) <- c(names_ab[-c(1L, length(names_ab))], names(gamma), names(lambda))

  # item coefficients
  coefs_item <- Map(function(a, b) c(Diff = a, Dscrmn = b), alpha, beta)

  # all coefficients
  coef_all <- c(unlist(coefs_item), gamma, lambda)
  coefs <- data.frame(Estimate = coef_all, Std_Error = se_all, z_value = coef_all/se_all,
                      p_value = 2 * (1 - pnorm(abs(coef_all/se_all))))
  rownames(coefs) <- names(se_all)

  # item constraints
  if (constr == "items"){

    gamma0_prev <- gamma[1L]

    # location constraint
    alpha_sum <- sum(alpha)
    beta_sum <- sum(beta)
    c1 <- alpha_sum/beta_sum
    gamma[1L] <- gamma[1L] + c1  # adjust gamma0
    alpha <- unlist(Map(function(x, y) x - c1 * y, alpha, beta))

    # scale constraint
    c2 <- 2 * mean(log(abs(beta)))
    gamma <- gamma * exp(c2/2)
    lambda[1L] <- lambda[1L] + c2
    beta <- beta / exp(c2/2)

    # fitted means and variances
    fitted_mean <- as.double(x %*% gamma)
    fitted_var <- exp(as.double(z %*% lambda))

    # theta_eap and theta_vap
    theta_eap <- (theta_eap - gamma0_prev) * exp(c2/2) + gamma[1L]
    theta_vap <- theta_vap * (exp(c2/2))^2

    # covmat for new parameterization
    tmp_fun <- function(d) {
      mat <- diag(d)
      mat[d, d] <- exp(-c2/2)
      mat[1:(d-1), d] <- rep(-c1, d-1)
      mat
    }
    A <- Reduce(Matrix::bdiag, lapply(H, tmp_fun))
    A2 <- A[seq(2, nrow(A)-1), seq(2, ncol(A)-1)]
    B <- Matrix::bdiag(exp(c2/2) * diag(p), diag(q))
    C <- Matrix::bdiag(A2, B)
    covmat <- C %*% Matrix::tcrossprod(covmat, C)

    se_all <- sqrt(Matrix::diag(covmat))

    # reorganize se_all
    sH <- 2 * J
    lambda_indices <- gamma_indices <- NULL
    gamma_indices <- (sH - 1):(sH + p - 2)
    lambda_indices <- (sH + p - 1):(sH + p + q - 2)
    se_all <- c(NA, se_all[1:(sH-2)], NA, se_all[gamma_indices], se_all[lambda_indices])

    # name se_all and covmat
    names_ab <- paste(rep(names(alpha), each = 2), c("Diff", "Dscrmn"))
    names(se_all) <- c(names_ab, names(gamma), names(lambda))
    rownames(covmat) <- colnames(covmat) <- c(names_ab[-c(1L, length(names_ab))], names(gamma), names(lambda))

    # item coefficients
    coefs_item <- Map(function(a, b) c(Diff = a, Dscrmn = b), alpha, beta)

    # all coefficients
    coef_all <- c(unlist(coefs_item), gamma, lambda)
    coefs <- data.frame(Estimate = coef_all, Std_Error = se_all, z_value = coef_all/se_all,
                        p_value = 2 * (1 - pnorm(abs(coef_all/se_all))))
    rownames(coefs) <- names(se_all)
  }

  # ability parameter estimates
  theta <- data.frame(post_mean = theta_eap, post_sd = sqrt(theta_vap),
                      prior_mean = fitted_mean, prior_sd = sqrt(fitted_var))

  # output
  out <- list(coefficients = coefs, scores = theta, vcov = covmat, log_Lik = log_Lik, constr = constr,
              N = N, J = J, H = H, ylevels = ylevels, p = p, q = q, control = con, call = cl)
  class(out) <- c("hltm", "hIRT")
  out
}

