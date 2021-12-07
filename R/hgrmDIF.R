#' Test Differential Item Functioning in Hierarchical Graded Response Models
#'
#' \code{hgrmDIF} fits a hierarchical graded response model similar to hgrm(), but person-specific
#' covariates \code{x} are allowed to affect item responses directly (not via the latent preference).
#' This model can be used to test for the presence of differential item functioning.
#'
#' @param y A data frame or matrix of item responses.
#' @param x An optional model matrix, including the intercept term, that predicts the
#'   mean of the latent preference. If not supplied, only the intercept term is included.
#' @param z An optional model matrix, including the intercept term, that predicts the
#'   variance of the latent preference. If not supplied, only the intercept term is included.
#' @param items_dif The indices of the items for which differential item functioning is tested.
#' @param form_dif Form of differential item functioning being tested. Either "uniform" or "non-uniform."
#' @param constr The type of constraints used to identify the model: "latent_scale",
#'   or "items". The default, "latent_scale" constrains the mean of latent preferences
#'   to zero and the geometric mean of prior variance to one; "items" places constraints
#'   on item parameters instead and sets the mean of item difficulty parameters to zero
#'   and the geometric mean of the discrimination parameters to one. Currently, only "latent_scale"
#'   is supported in hgrmDIF().
#' @param beta_set The index of the item for which the discrimination parameter is
#'   restricted to be positive (or negative). It may take any integer value from
#'   1 to \code{ncol(y)}.
#' @param sign_set Logical. Should the discrimination parameter of
#'   the corresponding item (indexed by \code{beta_set}) be positive
#'   (if \code{TRUE}) or negative (if \code{FALSE})?
#' @param init A character string indicating how item parameters are initialized. It can be
#'   "naive", "glm", or "irt".
#' @param control A list of control values
#' \describe{
#'  \item{max_iter}{The maximum number of iterations of the EM algorithm.
#'    The default is 150.}
#'  \item{eps}{Tolerance parameter used to determine convergence of the
#'   EM algorithm. Specifically, iterations continue until the Euclidean
#'   distance between \eqn{\beta_{n}} and \eqn{\beta_{n-1}} falls under \code{eps},
#'   where \eqn{\beta} is the vector of item discrimination parameters.
#'   \code{eps}=1e-4 by default.}
#'  \item{max_iter2}{The maximum number of iterations of the conditional
#'    maximization procedures for updating \eqn{\gamma} and \eqn{\lambda}.
#'    The default is 15.}
#'  \item{eps2}{Tolerance parameter used to determine convergence of the
#'    conditional maximization procedures for updating \eqn{\gamma} and
#'    \eqn{\lambda}. Specifically, iterations continue until the Euclidean
#'   distance between two consecutive log likelihoods falls under \code{eps2}.
#'   \code{eps2}=1e-3 by default.}
#'  \item{K}{Number of Gauss-Legendre quadrature points for the E-step. The default is 21.}
#'  \item{C}{[-C, C] sets the range of integral in the E-step. \code{C}=3 by default.}
#' }
#'
#' @return An object of class \code{hgrm}.
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
#'  \item{coef_item}{Item coefficient estimates.}
#'  \item{control}{List of control values.}
#'  \item{call}{The matched call.}
#' @importFrom rms lrm.fit
#' @importFrom pryr compose
#' @importFrom pryr partial
#' @importFrom ltm grm
#' @importFrom ltm ltm
#' @import stats
#' @export
#' @examples
#' y <- nes_econ2008[, -(1:3)]
#' x <- model.matrix( ~ party * educ, nes_econ2008)
#' nes_m2 <- hgrmDIF(y, x, items_DIF = 1:2)
#' coef_item(nes_m2)

hgrmDIF <- function(y, x = NULL, z = NULL, items_dif = 1L, form_dif = c("uniform", "non-uniform"),
                    constr = c("latent_scale"), beta_set = 1L, sign_set = TRUE,
                    init = c("naive", "glm", "irt"), control = list()) {

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
  y[] <- lapply(y, as.integer)
  if (!is.na(invalid <- match(TRUE, vapply(y, invalid_grm, logical(1L)))))
    stop(paste(names(y)[invalid], "does not have at least two valid responses"))
  H <- vapply(y, max, integer(1L), na.rm = TRUE)

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

  # check item, beta_set and sign_set
  stopifnot(beta_set %in% 1:J, is.logical(sign_set))

  # check constraint
  constr <- match.arg(constr)
  init <- match.arg(init)
  form_dif <- match.arg(form_dif)

  # control parameters
  con <- list(max_iter = 150, max_iter2 = 15, eps = 1e-03, eps2 = 1e-03, K = 25, C = 4)
  con[names(control)] <- control

  # set environments for utility functions
  environment(xDIF) <- environment(loglik_grmDIF) <- environment(theta_post_grmDIF) <- environment(dummy_fun_grm) <- environment(tab2df_grm) <- environment()

  # GL points
  K <- con[["K"]]
  theta_ls <- con[["C"]] * GLpoints[[K]][["x"]]
  qw_ls <- con[["C"]] * GLpoints[[K]][["w"]]

  # imputation
  y_imp <- y
  if(anyNA(y)) y_imp[] <- lapply(y, impute)

  # pca
  theta_eap <- {
    tmp <- princomp(y_imp, cor = TRUE)$scores[, 1]
    (tmp - mean(tmp, na.rm = TRUE))/sd(tmp, na.rm = TRUE)
  }

  # initialization of alpha and beta parameters
  if (init == "naive"){

    alpha <- lapply(H, function(x) c(Inf, seq(1, -1, length.out = x - 1), -Inf))
    beta <- vapply(y, function(y) cov(y, theta_eap, use = "complete.obs")/var(theta_eap), double(1L))

  } else if (init == "glm"){

    pseudo_lrm <- lapply(y_imp, function(y) lrm.fit(theta_eap, y)[["coefficients"]])
    beta <- vapply(pseudo_lrm, function(x) x[[length(x)]], double(1L))
    alpha <- lapply(pseudo_lrm, function(x) c(Inf, x[-length(x)], -Inf))

  } else {

    grm_coefs <- grm(y)[["coefficients"]]
    beta <- vapply(grm_coefs, function(x) x[[length(x)]], double(1L))
    alpha <- lapply(grm_coefs, function(x) c(Inf, rev(x[-length(x)]), -Inf))

  }

  # initialization of xDIFnames, x_lrm, and eta
  tmp_theta <- rep(theta_ls, N)
  tmp_x <- apply(x[, -1, drop = FALSE], 2, rep, each = K)
  if(form_dif == "uniform") {
    xDIFnames <- colnames(x)[-1]
    x_lrm <- cbind(tmp_theta, tmp_x)
    eta <- lapply(1:J, function(j) if(j %in% items_dif) rep(0, p - 1) else double(0L))
  } else{
    xDIFnames <- c(colnames(x)[-1], paste0("Dscrmn * ", colnames(x)[-1]))
    x_lrm <- cbind(tmp_theta, tmp_x, tmp_theta * tmp_x)
    eta <- lapply(1:J, function(j) if(j %in% items_dif) rep(0, 2 * (p - 1)) else double(0L))
  }
  colnames(x_lrm) <- c("theta", xDIFnames)
  names(eta) <- names(alpha) <- names(H)

  # initial values of gamma and lambda
  lm_opr <- tcrossprod(solve(crossprod(x)), x)
  gamma <- lm_opr %*% theta_eap
  lambda <- rep(0, q)
  fitted_mean <- as.double(x %*% gamma)
  fitted_var <- rep(1, N)
  pseudo_lrm <- vector(mode = "list", length = J)

  # EM algorithm
  for (iter in seq(1, con[["max_iter"]])) {

    # store previous parameters
    alpha_prev <- alpha
    beta_prev <- beta
    gamma_prev <- gamma
    lambda_prev <- lambda
    eta_prev <- eta

    # construct w_ik
    # list of length K, each element an N-vector
    posterior <- Map(theta_post_grmDIF, theta_ls, qw_ls)
    # K-by-N matrix
    w <- {
      tmp <- matrix(unlist(posterior), N, K)
      t(sweep(tmp, 1, rowSums(tmp), FUN = "/"))
    }

    # maximization with DIF
    w_lrm <- as.vector(w)
    for (j in seq(1, J)){
      y_lrm <- rep(y[[j]], each = K)
      pseudo_lrm[[j]] <- if(j %in% items_dif) lrm_fit(x_lrm, y_lrm, weights = w_lrm)[["coefficients"]] else
        lrm_fit(cbind(theta = tmp_theta), y_lrm, weights = w_lrm)[["coefficients"]]
    }
    beta <- vapply(pseudo_lrm, function(xx) xx[["theta"]], double(1L))
    eta <- lapply(1:J, function(j) if(j %in% items_dif) pseudo_lrm[[j]][-(1:H[[j]])] else double(0L))
    alpha <- Map(function(xx, H_j) c(Inf, xx[1:(H_j-1)], -Inf), pseudo_lrm, H)
    names(alpha) <- names(H)

    # # maximization in hgrm()
    # # list of length J, each element a K-by-H_j matrix containing
    # # the pseudo number of people with theta^k choosing h
    # pseudo_tab <- Map(dummy_fun_grm, y, H)
    # # list of length J, each element a data frame with K*H_j units and 3 columns: y, x(theta), and wt
    # pseudo_y <- lapply(pseudo_tab, tab2df_grm, theta_ls = theta_ls)
    # # weighted lrm on each element of pseudo_y
    # pseudo_lrm <- lapply(pseudo_y, function(df) lrm_fit(df[["x"]], df[["y"]], weights = df[["wt"]])[["coefficients"]])
    # beta <- vapply(pseudo_lrm, function(x) x[[length(x)]], double(1L))
    # alpha <- lapply(pseudo_lrm, function(x) c(Inf, x[-length(x)], -Inf))

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
    alpha <- Map(function(x, y) x + tmp * y, alpha, beta)
    gamma[[1L]] <- gamma[[1L]] - tmp

    # scale constraint
    tmp <- mean(z %*% lambda)
    gamma <- gamma/exp(tmp/2)
    beta <- beta * exp(tmp/2)
    lambda[[1L]] <- lambda[[1L]] - tmp

    # direction contraint
    if (sign_set == (beta[[beta_set]] < 0)) {
      gamma <- -gamma
      beta <- -beta
    }

    fitted_mean <- as.double(x %*% gamma)
    fitted_var <- exp(as.double(z %*% lambda))
    # cat(beta, "\n")
    # cat(abs(beta - beta_prev), "\n")

    cat(".")

    # check convergence
    if (sqrt(mean((beta - beta_prev)^2)) < con[["eps"]]) {
      cat("\n converged at iteration", iter, "\n")
      break
    } else if (iter == con[["max_iter"]]) {
      stop("algorithm did not converge; try increasing `max_iter` or decreasing `eps`")
      break
    } else next
  }

  gamma <- setNames(as.double(gamma), paste("x", colnames(x), sep = ""))
  lambda <- setNames(as.double(lambda), paste("z", colnames(z), sep = ""))

  # inference
  # pik equals p_ik * w_k in Zhou2019supp
  pik <- matrix(unlist(Map(partial(dnorm, x = theta_ls), mean = fitted_mean, sd = sqrt(fitted_var))),
                N, K, byrow = TRUE) * matrix(qw_ls, N, K, byrow = TRUE)
  Lijk <- lapply(theta_ls, function(theta_k) exp(loglik_grmDIF(alpha = alpha, beta = beta, eta = eta, rep(theta_k, N))))  # K-list
  Lik <- vapply(Lijk, compose(exp, partial(rowSums, na.rm = TRUE), log), double(N))
  Li <- rowSums(Lik * pik)

  # log likelihood
  log_Lik <- sum(log(Li))

  # outer product of gradients
  environment(sj_ab_grmDIF) <- environment(sj_ab_grm) <- environment(si_gamma) <- environment(si_lambda) <- environment()
  s_abe <- unname(Reduce(rbind, lapply(1:J, sj_ab_grmDIF)))
  s_lambda <- s_gamma <- NULL
  s_gamma <- vapply(1:N, si_gamma, double(p))
  s_lambda <- vapply(1:N, si_lambda, double(q))

  # covariance matrix and standard errors
  s_all <- rbind(s_abe[-c(1L, nrow(s_abe)), , drop = FALSE], s_gamma, s_lambda)
  s_all[is.na(s_all)] <- 0
  covmat <- tryCatch(solve(tcrossprod(s_all)),
                     error = function(e) {warning("The information matrix is singular; SE calculation failed.");
                       matrix(NA, nrow(s_all), nrow(s_all))})
  se_all <- sqrt(diag(covmat))

  # reorganize se_all
  sH <- sum(H) + (p - 1) * length(items_dif) + (p - 1) * as.double(form_dif == "non-uniform") * length(items_dif)
  gamma_indices <- (sH - 1):(sH + p - 2)
  lambda_indices <- (sH + p - 1):(sH + p + q - 2)
  se_all <- c(NA, se_all[1:(sH-2)], NA, se_all[gamma_indices], se_all[lambda_indices])

  # name se_all and covmat
  names(alpha) <- names(H)
  names_abe <- unlist(lapply(1:J, function(j) {
    itemname <- names(alpha)[[j]]
    tmp <- alpha[[j]]
    if(j %in% items_dif) paste(itemname, c(names(tmp)[-c(1L, length(tmp))], xDIFnames, "Dscrmn")) else{
      paste(itemname, c(names(tmp)[-c(1L, length(tmp))], "Dscrmn"))
    }
  }))
  names(se_all) <- c(names_abe, names(gamma), names(lambda))
  rownames(covmat) <- colnames(covmat) <- c(names_abe[-c(1L, length(names_abe))], names(gamma), names(lambda))

  # item coefficients
  coef_item <- Map(function(a, b, c) c(a[-c(1L, length(a))], c, Dscrmn = b), alpha, beta, eta)

  # all coefficients
  coef_all <- c(unlist(coef_item), gamma, lambda)
  coefs <- data.frame(Estimate = coef_all, Std_Error = se_all, z_value = coef_all/se_all,
                      p_value = 2 * (1 - pnorm(abs(coef_all/se_all))))

  # ability parameter estimates
  theta <- data.frame(post_mean = theta_eap, post_sd = sqrt(theta_vap),
                      prior_mean = fitted_mean, prior_sd = sqrt(fitted_var))

  # output
  out <- list(coefficients = coefs, scores = theta, vcov = covmat, log_Lik = log_Lik,
              items_dif = items_dif,  form_dif = form_dif, constr = constr,
              N = N, J = J, H = H, ylevels = ylevels, p = p, q = q, coef_item = coef_item,
              control = con, call = cl)
  class(out) <- c("hgrmDIF")
  out
}

