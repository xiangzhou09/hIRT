#' Hierarchical Graded Response Models with Known Item Parameters
#'
#' \code{hgrm2} fits a hierarchical graded response model where the item parameters
#'   are known and supplied by the user.
#'
#' @param y A data frame or matrix of item responses.
#' @param x An optional model matrix, including the intercept term, that predicts the
#'   mean of the latent preference. If not supplied, only the intercept term is included.
#' @param z An optional model matrix, including the intercept term, that predicts the
#'   variance of the latent preference. If not supplied, only the intercept term is included.
#' @param item_coefs A list of known item parameters. The parameters of item \eqn{j} are given
#'   by the \eqn{j}th element, which should be a vector of length \eqn{H_j}, containing
#'   \eqn{H_j - 1} item difficulty parameters (in descending order) and one item discrimination
#'   parameter.
#' @param control A list of control values
#' \describe{
#'  \item{max_iter}{The maximum number of iterations of the EM algorithm.
#'    The default is 150.}
#'  \item{eps}{Tolerance parameter used to determine convergence of the
#'   EM algorithm. Specifically, iterations continue until the Euclidean
#'   distance between \eqn{\beta_{n}} and \eqn{\beta_{n-1}} falls under \code{eps},
#'   where \eqn{\beta} is the vector of item discrimination parameters.
#'   \code{eps}=1e-3 by default.}
#'  \item{max_iter2}{The maximum number of iterations of the conditional
#'    maximization procedures for updating \eqn{\gamma} and \eqn{\lambda}.
#'    The default is 15.}
#'  \item{eps2}{Tolerance parameter used to determine convergence of the
#'    conditional maximization procedures for updating \eqn{\gamma} and
#'    \eqn{\lambda}. Specifically, iterations continue until the Euclidean
#'   distance between two consecutive log likelihoods falls under \code{eps2}.
#'   \code{eps2}=1e-3 by default.}
#'  \item{K}{Number of Gauss-Legendre quadrature points for the E-step. The default is 25.}
#'  \item{C}{[-C, C] sets the range of integral in the E-step. \code{C}=4 by default.}
#'  \item{verbose}{Logical. Should the progress of the EM algorithm be printed
#'   to the console? \code{FALSE} by default.}
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
#'  \item{control}{List of control values.}
#'  \item{call}{The matched call.}
#' @importFrom rms lrm.fit
#' @import stats
#' @export
#' @examples
#'
#' y <- nes_econ2008[, -(1:3)]
#' x <- model.matrix( ~ party * educ, nes_econ2008)
#' z <- model.matrix( ~ party, nes_econ2008)
#'
#' n <- nrow(nes_econ2008)
#' id_train <- sample.int(n, n/4)
#' id_test <- setdiff(1:n, id_train)
#'
#' y_train <- y[id_train, ]
#' x_train <- x[id_train, ]
#' z_train <- z[id_train, ]
#'
#' mod_train <- hgrm(y_train, x_train, z_train)
#'
#' y_test <- y[id_test, ]
#' x_test <- x[id_test, ]
#' z_test <- z[id_test, ]
#'
#' item_coefs <- lapply(coef_item(mod_train), `[[`, "Estimate")
#'
#' model_test <- hgrm2(y_test, x_test, z_test, item_coefs = item_coefs)

hgrm2 <- function(y, x = NULL, z = NULL, item_coefs, control = list()) {

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
  Y <- y_mat(y)

  # extract item parameters
  if(missing(item_coefs))
    stop("`item_coefs` must be supplied.")
  if(!is.list(item_coefs) || length(item_coefs) != J)
    stop("`item_coefs` must be a list of `ncol(y)` elements")
  item_coefs_H <- vapply(item_coefs, length, integer(1L))
  if(!all.equal(item_coefs_H, H))
    stop("`item_coefs` do not match the number of response categories in `y`")
  alpha <- lapply(item_coefs, function(x) c(Inf, x[-length(x)], -Inf))
  beta <- vapply(item_coefs, function(x) x[[length(x)]], double(1L))

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

  # control parameters
  con <- list(max_iter = 150, max_iter2 = 15, eps = 1e-03, eps2 = 1e-03, K = 25, C = 4,
              verbose = FALSE)
  con[names(control)] <- control

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

  # initial values of gamma and lambda
  lm_opr <- tcrossprod(solve(crossprod(x)), x)
  gamma <- lm_opr %*% theta_eap
  lambda <- rep(0, q)
  fitted_mean <- as.double(x %*% gamma)
  fitted_var <- rep(1, N)

  # EM algorithm
  for (iter in seq(1, con[["max_iter"]])) {

    # store previous parameters
    # alpha_prev <- alpha
    # beta_prev <- beta
    gamma_prev <- gamma
    lambda_prev <- lambda

    # construct w_ik
    w <- estep_grm(Y, alpha, beta, theta_ls, qw_ls, fitted_mean, fitted_var)[["w"]]

    # # maximization
    # pseudo_tab <- Map(dummy_fun_grm, y, H)
    # pseudo_y <- lapply(pseudo_tab, tab2df_grm, theta_ls = theta_ls)
    # pseudo_lrm <- lapply(pseudo_y, function(df) lrm.fit(df[["x"]], df[["y"]], weights = df[["wt"]])[["coefficients"]])
    # beta <- vapply(pseudo_lrm, function(x) x[[length(x)]], double(1L))
    # alpha <- lapply(pseudo_lrm, function(x) c(Inf, x[-length(x)], -Inf))

    # EAP and VAP estimates of latent preferences
    theta_eap <- t(theta_ls %*% w)
    theta_vap <- t(theta_ls^2 %*% w) - theta_eap^2

    # variance regression
    vr <- var_regression(x, z, theta_eap, theta_vap, lm_opr, con)
    gamma <- vr[["gamma"]]
    lambda <- vr[["lambda"]]

    fitted_mean <- as.double(x %*% gamma)
    fitted_var <- exp(as.double(z %*% lambda))
    if (con[["verbose"]]) cat(".")

    if (sqrt(mean((gamma/gamma_prev - 1)^2)) < con[["eps"]]) {
      if (con[["verbose"]]) cat("\n converged at iteration", iter, "\n")
      break
    } else if (iter == con[["max_iter"]]) {
      stop("algorithm did not converge; try increasing `max_iter` or decreasing `eps`")
      break
    } else next
  }

  gamma <- setNames(as.double(gamma), paste("x", colnames(x), sep = ""))
  lambda <- setNames(as.double(lambda), paste("z", colnames(z), sep = ""))

  # inference
  sc <- scores_grm(Y, alpha, beta, theta_ls, qw_ls, fitted_mean, fitted_var, x, z, H,
                   items = FALSE)

  # log likelihood
  log_Lik <- sc[["loglik"]]

  # outer product of gradients
  s_gamma <- sc[["s_gamma"]]
  s_lambda <- sc[["s_lambda"]]

  # covariance matrix
  s_all <- rbind(s_gamma, s_lambda)
  s_all[is.na(s_all)] <- 0
  covmat <- tryCatch(solve(tcrossprod(s_all)),
                     error = function(e) {warning("The information matrix is singular; SE calculation failed.");
                       matrix(NA, nrow(s_all), nrow(s_all))})

  # reorganize se_all
  sH <- sum(H)
  gamma_indices <- (sH - 1):(sH + p - 2)
  lambda_indices <- (sH + p - 1):(sH + p + q - 2)
  se_all <- c(rep(0, sH), sqrt(diag(covmat)))

  # name se_all and covmat
  names_ab <- unlist(lapply(names(alpha), function(x) {
    tmp <- alpha[[x]]
    paste(x, c(paste0("y>=", seq(2, length(tmp)-1)), "Dscrmn"))
  }))
  names(se_all) <- c(names_ab, names(gamma), names(lambda))
  rownames(covmat) <- colnames(covmat) <- c(names(gamma), names(lambda))

  # item coefficients
  coef_item <- Map(function(a, b) c(a[-c(1L, length(a))], Dscrmn = b), alpha, beta)

  # all coefficients
  coef_all <- c(unlist(coef_item), gamma, lambda)
  coefs <- data.frame(Estimate = coef_all, Std_Error = se_all, z_value = coef_all/se_all,
                      p_value = 2 * (1 - pnorm(abs(coef_all/se_all))))
  rownames(coefs) <- names(se_all)

  # ability parameter estimates
  theta <- data.frame(post_mean = theta_eap, post_sd = sqrt(theta_vap),
                      prior_mean = fitted_mean, prior_sd = sqrt(fitted_var))

  # output
  out <- list(coefficients = coefs, scores = theta, vcov = covmat, log_Lik = log_Lik,
              N = N, J = J, H = H, ylevels = ylevels, p = p, q = q, control = con, call = cl)
  class(out) <- c("hgrm", "hIRT")
  out
}
