#' Hierarchical Latent Trait Models with Known Item Parameters.
#'
#' \code{hltm2} fits a hierarchical latent trait model where the item parameters
#'   are known and supplied by the user.
#'
#' @inheritParams hgrm2
#' @param item_coefs A list of known item parameters. The parameters of item \eqn{j} are given
#'   by the \eqn{j}th element, which should be a vector of length 2, containing
#'   the item difficulty parameter and item discrimination parameter.
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
#' @importFrom rms lrm.fit
#' @importFrom pryr compose
#' @importFrom pryr partial
#' @import stats
#' @export
#' @examples
#' y <- nes_econ2008[, -(1:3)]
#' x <- model.matrix( ~ party * educ, nes_econ2008)
#' z <- model.matrix( ~ party, nes_econ2008)
#' dichotomize <- function(x) findInterval(x, c(mean(x, na.rm = TRUE)))
#' y_bin <- y
#' y_bin[] <- lapply(y, dichotomize)
#'
#' n <- nrow(nes_econ2008)
#' id_train <- sample.int(n, n/4)
#' id_test <- setdiff(1:n, id_train)
#'
#' y_bin_train <- y_bin[id_train, ]
#' x_train <- x[id_train, ]
#' z_train <- z[id_train, ]
#'
#' mod_train <- hltm(y_bin_train, x_train, z_train)
#'
#' y_bin_test <- y_bin[id_test, ]
#' x_test <- x[id_test, ]
#' z_test <- z[id_test, ]
#'
#' item_coefs <- lapply(coef_item(mod_train), `[[`, "Estimate")
#'
#' model_test <- hltm2(y_bin_test, x_test, z_test, item_coefs = item_coefs)

hltm2 <- function(y, x = NULL, z = NULL, item_coefs, control = list()) {

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

  # extract item parameters
  if(missing(item_coefs))
    stop("`item_coefs` must be supplied.")
  if(!is.list(item_coefs) || length(item_coefs) != J)
    stop("`item_coefs` must be a list of `ncol(y)` elements")
  item_coefs_H <- vapply(item_coefs, length, integer(1L))
  if(!all.equal(item_coefs_H, H))
    stop("`item_coefs` do not match the number of response categories in `y`")
  alpha <- vapply(item_coefs, function(x) x[[1L]], double(1L))
  beta <- vapply(item_coefs, function(x) x[[2L]], double(1L))

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
  con <- list(max_iter = 150, max_iter2 = 15, eps = 1e-04, eps2 = 1e-03, K = 21, C = 3)
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
    posterior <- Map(theta_post_ltm, theta_ls, qw_ls)
    w <- {
      tmp <- matrix(unlist(posterior), N, K)
      t(sweep(tmp, 1, rowSums(tmp), FUN = "/"))
    }

    # # maximization
    # pseudo_tab <- lapply(y, dummy_fun_ltm)
    # pseudo_y <- lapply(pseudo_tab, tab2df_ltm, theta_ls = theta_ls)
    # pseudo_logit <- lapply(pseudo_y, function(df) glm.fit(cbind(1, df[["x"]]),
    #                                                       df[["y"]], weights = df[["wt"]], family = quasibinomial("logit"))[["coefficients"]])
    # beta <- vapply(pseudo_logit, function(x) x[2L], double(1L))
    # alpha <- vapply(pseudo_logit, function(x) x[1L], double(1L))

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

    fitted_mean <- as.double(x %*% gamma)
    fitted_var <- exp(as.double(z %*% lambda))
    cat(".")

    # check convergence
    if (sqrt(mean((gamma/gamma_prev - 1)^2)) < con[["eps"]]) {
      cat("\n converged at iteration", iter, "\n")
      break
    } else if (iter == con[["max_iter"]]) {
      stop("algorithm did not converge; try increasing max_iter.")
      break
    } else next
  }

  gamma <- setNames(as.double(gamma), paste("x", colnames(x), sep = "_"))
  lambda <- setNames(as.double(lambda), paste("z", colnames(z), sep = "_"))

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
  # s_ab <- unname(Reduce(cbind, lapply(1:J, sj_ab_ltm)))
  s_gamma <- vapply(1:N, si_gamma, double(p))
  s_lambda <- vapply(1:N, si_lambda, double(q))

  s_all <- rbind(s_gamma, s_lambda)
  s_all[is.na(s_all)] <- 0
  covmat <- tryCatch(solve(tcrossprod(s_all)),
                     error = function(e) {warning("The information matrix is singular; SE calculation failed.");
                       matrix(NA, nrow(s_all), nrow(s_all))})
  se_all <- sqrt(diag(covmat))

  # reorganize se_all
  sH <- 2 * J
  gamma_indices <- (sH - 1):(sH + p - 2)
  lambda_indices <- (sH + p - 1):(sH + p + q - 2)
  se_all <- c(rep(0, sH), sqrt(diag(covmat)))

  # name se_all and covmat
  names_ab <- paste(rep(names(alpha), each = 2), c("Diff", "Dscrmn"))
  names(se_all) <- c(names_ab, names(gamma), names(lambda))
  rownames(covmat) <- colnames(covmat) <- c(names(gamma), names(lambda))

  # item coefficients
  coefs_item <- Map(function(a, b) c(Diff = a, Dscrmn = b), alpha, beta)

  # all coefficients
  coef_all <- c(unlist(coefs_item), gamma, lambda)
  coefs <- data.frame(Estimate = coef_all, Std_Error = se_all, z_value = coef_all/se_all,
                      p_value = 2 * (1 - pnorm(abs(coef_all/se_all))))
  rownames(coefs) <- names(se_all)

  # ability parameter estimates
  theta <- data.frame(post_mean = theta_eap, post_sd = sqrt(theta_vap),
                      prior_mean = fitted_mean, prior_sd = sqrt(fitted_var))

  # output
  out <- list(coefficients = coefs, scores = theta, vcov = covmat, log_Lik = log_Lik,
              N = N, J = J, H = H, ylevels = ylevels, p = p, q = q, control = con, call = cl)
  class(out) <- c("hltm", "hIRT")
  out
}
