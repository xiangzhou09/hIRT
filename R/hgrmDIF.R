#' Hierarchical Graded Response Models with Differential Item Functioning
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
#' @param x0 A matrix specifying the covariates by which differential item functioning operates. If not supplied,
#'   \code{x0} is taken to be a matrix containing all predictors in \code{x} except the intercept.
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
#'  \item{p0}{The number of predictors for items with DIF.}
#'  \item{coef_item}{Item coefficient estimates.}
#'  \item{control}{List of control values.}
#'  \item{call}{The matched call.}
#' @importFrom rms lrm.fit
#' @importFrom ltm grm
#' @importFrom ltm ltm
#' @import stats
#' @export
#' @examples
#' \donttest{
#' y <- nes_econ2008[, -(1:3)]
#' x <- model.matrix( ~ party * educ, nes_econ2008)
#' nes_m2 <- hgrmDIF(y, x, items_dif = 1:2)
#' coef_item(nes_m2)
#' }

hgrmDIF <- function(y, x = NULL, z = NULL, x0 = x[, -1, drop = FALSE],
                    items_dif = 1L, form_dif = c("uniform", "non-uniform"),
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
  Y <- y_mat(y)

  # check x0
  if(!is.matrix(x0)) x0 <- as.matrix(x0)

  # check x and z (x and z should contain an intercept column)
  x <- x %||% as.matrix(rep(1, N))
  z <- z %||% as.matrix(rep(1, N))
  if (!is.matrix(x)) stop("`x` must be a matrix.")
  if (!is.matrix(z)) stop("`z` must be a matrix.")
  if (nrow(x) != N || nrow(z) != N) stop("both 'x' and 'z' must have the same number of rows as 'y'")
  p <- ncol(x)
  q <- ncol(z)
  p0 <- ncol(x0)
  colnames(x) <- colnames(x) %||% paste0("x", 1:p)
  colnames(z) <- colnames(z) %||% paste0("z", 1:q)
  colnames(x0) <- colnames(x0) %||% paste0("x0", 1:p0)

  # check item, beta_set and sign_set
  stopifnot(beta_set %in% 1:J, is.logical(sign_set))

  # check constraint
  constr <- match.arg(constr)
  init <- match.arg(init)
  form_dif <- match.arg(form_dif)

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

  # initialization of x0DIFnames and eta
  if(form_dif == "uniform") {
    x0DIFnames <- colnames(x0)
    eta <- lapply(1:J, function(j) if(j %in% items_dif) rep(0, p0) else double(0L))
  } else{
    x0DIFnames <- c(colnames(x0), paste0("Dscrmn * ", colnames(x0)))
    eta <- lapply(1:J, function(j) if(j %in% items_dif) rep(0, 2 * p0) else double(0L))
  }
  lrm_names <- c("theta", x0DIFnames)
  dif_flag <- as.integer(seq_len(J) %in% items_dif)
  x0_key <- {
    tmp <- do.call(paste, c(unname(as.data.frame(x0)), sep = "\r"))
    match(tmp, unique(tmp))
  }
  x0_u <- x0[!duplicated(x0_key), , drop = FALSE]
  # N x J offsets entering the linear predictor: off0 + theta * off1
  dif_offsets <- function(eta) {
    off0 <- off1 <- matrix(0, N, J)
    for (j in items_dif) {
      off0[, j] <- x0 %*% eta[[j]][seq_len(p0)]
      if (form_dif == "non-uniform") off1[, j] <- x0 %*% eta[[j]][p0 + seq_len(p0)]
    }
    list(off0 = off0, off1 = off1)
  }
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
    off <- dif_offsets(eta)
    w <- estep_grm(Y, alpha, beta, theta_ls, qw_ls, fitted_mean, fitted_var,
                   off[["off0"]], off[["off1"]], dif_flag)[["w"]]

    # maximization with DIF
    # Respondents with the same x0 pattern and response contribute identical
    # rows at each node, so their weights are pooled before fitting; this
    # leaves the weighted likelihood unchanged.
    tw <- t(w)
    for (j in seq(1, J)){
      init_j <- c(alpha[[j]][2:H[[j]]], beta[[j]], eta[[j]])
      is_dif <- j %in% items_dif
      obs <- which(!is.na(y[[j]]))
      cell <- (if (is_dif) x0_key[obs] - 1L else 0L) * H[[j]] + y[[j]][obs]
      agg <- rowsum(tw[obs, , drop = FALSE], cell)          # cells x K
      cg <- as.integer(rownames(agg))
      y_lrm <- rep((cg - 1L) %% H[[j]] + 1L, each = K)
      theta_lrm <- rep(theta_ls, length(cg))
      X_lrm <- if (is_dif) {
        x0_cell <- x0_u[rep((cg - 1L) %/% H[[j]] + 1L, each = K), , drop = FALSE]
        if (form_dif == "uniform") cbind(theta_lrm, x0_cell) else
          cbind(theta_lrm, x0_cell, theta_lrm * x0_cell)
      } else cbind(theta_lrm)
      colnames(X_lrm) <- if (is_dif) lrm_names else "theta"
      pseudo_lrm[[j]] <- cumlogit_fit(y_lrm, X_lrm, as.double(t(agg)), H[[j]], init_j)
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
    vr <- var_regression(x, z, theta_eap, theta_vap, lm_opr, con)
    gamma <- vr[["gamma"]]
    lambda <- vr[["lambda"]]

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

    if (con[["verbose"]]) cat(".")

    # check convergence
    if (sqrt(mean((beta - beta_prev)^2)) < con[["eps"]]) {
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
  # pik equals p_ik * w_k in Zhou2019supp
  off <- dif_offsets(eta)
  sc <- scores_grm(Y, alpha, beta, theta_ls, qw_ls, fitted_mean, fitted_var, x, z, H,
                   off0 = off[["off0"]], off1 = off[["off1"]], X0 = x0, dif = dif_flag,
                   nonuniform = form_dif == "non-uniform")

  # log likelihood
  log_Lik <- sc[["loglik"]]

  # outer product of gradients
  s_abe <- sc[["s_item"]]
  s_gamma <- sc[["s_gamma"]]
  s_lambda <- sc[["s_lambda"]]

  # covariance matrix and standard errors
  s_all <- rbind(s_abe[-c(1L, nrow(s_abe)), , drop = FALSE], s_gamma, s_lambda)
  s_all[is.na(s_all)] <- 0
  covmat <- tryCatch(solve(tcrossprod(s_all)),
                     error = function(e) {warning("The information matrix is singular; SE calculation failed.");
                       matrix(NA, nrow(s_all), nrow(s_all))})
  se_all <- sqrt(diag(covmat))

  # reorganize se_all
  sH <- sum(H) + p0 * length(items_dif) + p0 * as.double(form_dif == "non-uniform") * length(items_dif)
  gamma_indices <- (sH - 1):(sH + p - 2)
  lambda_indices <- (sH + p - 1):(sH + p + q - 2)
  se_all <- c(NA, se_all[1:(sH-2)], NA, se_all[gamma_indices], se_all[lambda_indices])

  # name se_all and covmat
  names(alpha) <- names(H)
  names_abe <- unlist(lapply(1:J, function(j) {
    itemname <- names(alpha)[[j]]
    tmp <- alpha[[j]]
    if(j %in% items_dif) paste(itemname, c(names(tmp)[-c(1L, length(tmp))], x0DIFnames, "Dscrmn")) else{
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
              N = N, J = J, H = H, ylevels = ylevels, p = p, q = q, p0 = p0, coef_item = coef_item,
              control = con, call = cl)
  class(out) <- c("hgrmDIF")
  out
}

