#' Predicting from Hierarchical IRT Models
#'
#' \code{predict.hIRT} predicts latent preferences for new data (\eqn{y}, \eqn{x}, \eqn{z})
#'   using a fitted \code{hIRT} object. Both posterior (given \eqn{y}, \eqn{x}, \eqn{z}) and
#'   prior (given \eqn{x}, \eqn{z}) estimates are returned.
#'
#' @param object An object of class \code{hIRT}.
#' @param y A data frame or matrix of item responses.
#' @param x Model matrix, including the intercept term, representing the covariates
#'   used to predict the means.
#' @param z Model matrix, including the intercept term, representing the covariates
#'   used to predict the variances.
#'
#' @return A data frame.
#'  \item{post_mean}{Posterior mean estimates of latent preferences, i.e., \eqn{E[\theta|y,x,z]}}
#'  \item{post_sd}{Posterior standard deviation estimates of latent preferences, i.e., \eqn{\sqrt{V[\theta|y,x,z]}}}
#'  \item{prior_mean}{Prior mean estimates of latent preferences, i.e., \eqn{E[\theta|x]}}
#'  \item{prior_sd}{Prior standard deviation estimates of latent preferences, i.e., \eqn{\sqrt{V[\theta|z]}}}
#'
#' @export
#' @examples
#' n <- nrow(nes_econ2008)
#' n_train <- sample.int(n, 3*n/4)
#' n_test <- setdiff(1:n, n_train)
#'
#' nes_econ2008_train <- nes_econ2008[n_train, , drop = FALSE]
#' nes_econ2008_test <- nes_econ2008[n_test, , drop = FALSE]
#'
#' y <- nes_econ2008_train[, -(1:3)]
#' x <- model.matrix( ~ party * educ, nes_econ2008_train)
#' z <- model.matrix( ~ party, nes_econ2008_train)
#' mod_train <- hgrm(y, x, z)
#'
#' y2 <- nes_econ2008_test[, -(1:3)]
#' x2 <- model.matrix( ~ party * educ, nes_econ2008_test)
#' z2 <- model.matrix( ~ party, nes_econ2008_test)
#' preds <- predict_hIRT(mod_train, y2, x2, z2)
#' summary(preds)
#'
predict_hIRT <- function(object, y, x, z) {

  # check object
  if (!inherits(object, "hIRT"))
    stop("Use only with 'hIRT' objects.\n")

  # check y
  if ((!is.data.frame(y) & !is.matrix(y)))
    stop("'y' must be either a data.frame or a matrix")
  y <- as.data.frame(y)
  N <- nrow(y)
  J <- ncol(y)
  if (J != object$J)
    stop("'y' must have the same number of columns as 'object$J'")

  # check x and z (x and z should contain an intercept column)
  if (missing(x)) x <- as.matrix(rep(1, N))
  if (missing(z)) z <- as.matrix(rep(1, N))
  if (is.null(nrow(x))) x <- as.matrix(x)
  if (is.null(nrow(x))) z <- as.matrix(z)
  if (nrow(x) != N || nrow(z) != N)
    stop("both 'x' and 'z' must have the same number of rows as 'y'")
  if (ncol(x) != object$p || ncol(z) != object$q)
    stop("'x' and 'z' must have the same number of columns as 'object$p' and 'object$q', respectively")
  x <- `colnames<-`(model.matrix(~ 0 + x), colnames(x) %||% paste("x", 1:object$p, sep = ""))
  z <- `colnames<-`(model.matrix(~ 0 + z), colnames(z) %||% paste("z", 1:object$q, sep = ""))

  # obtain lists of theta_k and w_k parameters
  con <- object$control
  K <- con$K
  theta_ls <- con[["C"]] * GLpoints[[K]][["x"]]
  qw_ls <- con[["C"]] * GLpoints[[K]][["w"]]

  # get parameters gamma, lambda
  gamma <- coef_mean(object, digits = 7)[["Estimate"]]
  lambda <- coef_var(object, digits = 7)[["Estimate"]]

  # fitted mean and var
  fitted_mean <- as.vector(x %*% gamma)
  fitted_var <- exp(as.vector(z %*% lambda))

  if(inherits(object, "hgrm")){

    # preprocess y
    for (j in seq(1, J)) {
      y[[j]] <- as.integer(factor(y[[j]], levels = object$ylevels[[j]], exclude = c(NA, NaN)))
    }

    # get parameters alpha, beta
    tmp <- coef_item(object, digits = 7)
    alpha <- lapply(tmp, function(x) c(Inf, x[-nrow(x), "Estimate"], -Inf))
    beta <- vapply(tmp, function(x) x[nrow(x), "Estimate"], numeric(1L))

    # set environments
    environment(loglik_grm) <- environment(theta_post_grm) <- environment()

    # posterior densities
    posterior <- Map(theta_post_grm, theta_ls, qw_ls)

  } else if (inherits(object, "hltm")){

    # preprocess y
    for (j in seq(1, J)) {
      y[[j]] <- as.integer(factor(y[[j]], levels = object$ylevels[[j]], exclude = c(NA, NaN))) - 1
    }

    # get parameters alpha, beta
    tmp <- coef_item(object, digits = 7)
    alpha <- vapply(tmp, function(x) x[1L, "Estimate"], numeric(1L))
    beta <- vapply(tmp, function(x) x[2L, "Estimate"], numeric(1L))

    # set environments
    environment(loglik_ltm) <- environment(theta_post_ltm) <- environment()

    # posterior densities
    posterior <- Map(theta_post_ltm, theta_ls, qw_ls)

  } else stop("Use only with 'hgrm' or 'hltm' objects.\n")

  # EAP and VAP estimates
  tmp <- matrix(unlist(posterior), N, K)
  w <- t(sweep(tmp, 1, rowSums(tmp), FUN = "/"))
  theta_eap <- t(theta_ls %*% w)
  theta_vap <- t(theta_ls^2 %*% w) - theta_eap^2

  # output
  out <- data.frame(post_mean = theta_eap, post_sd = sqrt(theta_vap),
                    prior_mean = fitted_mean, prior_sd = sqrt(fitted_var))
}
