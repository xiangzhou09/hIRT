#' Fitting Multidimensional Hierarchical Graded Response Models
#'
#' \code{mhgrm} fits a hierarchical graded response model with up to three
#' latent dimensions (for example, economic, civil rights, and moral
#' attitudes). Each item measures one dimension (simple structure). The latent
#' vector follows a multivariate normal prior whose mean depends on \code{x}
#' and whose covariance matrix, including the covariances between dimensions,
#' depends on \code{z}.
#'
#' The covariance matrix is parameterized by its modified Cholesky
#' decomposition: with \eqn{\epsilon_i = \theta_i - \mu_i},
#' \deqn{\epsilon_{i1} = e_{i1}, \quad
#'   \epsilon_{ik} = \sum_{l<k} \phi_{kl}(z_i)\epsilon_{il} + e_{ik}, \quad
#'   e_{ik} \sim N(0, \exp(z_i'\gamma_k)),}
#' where \eqn{\phi_{kl}(z_i) = z_i'\lambda_{kl}}. Any values of \eqn{\gamma}
#' and \eqn{\lambda} give a positive definite covariance matrix. The
#' coefficients \eqn{\lambda_{kl}} are slopes from regressing dimension
#' \eqn{k} on earlier dimensions, which depend on the order of the
#' dimensions; use \code{\link{mhgrm_cov}} to obtain the implied variances,
#' covariances, and correlations at given covariate values.
#'
#' The model is fitted by an EM algorithm. By default the E-step uses a fixed
#' product grid of Gauss-Legendre nodes shared by all respondents, as in
#' \code{\link{hgrm}}. A fixed grid can be too coarse when dimensions have many
#' items (so that posteriors are narrow relative to the grid spacing) or when a
#' latent correlation is close to one. With \code{adapt = TRUE}, the E-step
#' uses adaptive Gauss-Hermite quadrature instead: each respondent's product
#' grid is centered on the current posterior mean and scaled by the posterior
#' standard deviation on each dimension. This is more accurate per node but
#' slower per iteration, because the nodes differ across respondents. The grid
#' adapts along each axis but is not rotated, so it cannot resolve a posterior
#' concentrated along a diagonal; if a latent correlation approaches one, check
#' the log-likelihood at a larger \code{K}.
#'
#' On a fixed grid the update of the level-II parameters is the M-step for a
#' continuous normal prior rather than for its discretized version, so when
#' prior distributions drift off the grid (typically when a dimension is weakly
#' identified, with few or weakly related items) the log-likelihood can fall.
#' Brief dips are tolerated, but if the log-likelihood falls far below its
#' maximum or stays below it for 50 iterations, or if the level-II update
#' fails, the algorithm stops with a warning and returns the estimates with the
#' highest log-likelihood; it does the same, with a warning, if the final
#' estimates are worse than earlier ones. A warning
#' is also given when more than 1\% of some respondent's prior distribution
#' lies outside the grid \eqn{[-C, C]}.
#'
#' The M-step fits one cumulative logit model per item and
#' updates the level-II parameters by conditional maximization. Standard errors come from the outer product of gradients,
#' computed under the identification constraints.
#'
#' @param y A data frame or matrix of item responses.
#' @param dims A vector of length \code{ncol(y)} assigning each item to a
#'   dimension. It can be numeric (1, 2, 3) or character/factor (dimension
#'   names, used in the output, in order of first appearance unless a factor).
#'   Each dimension needs at least two items.
#' @param x An optional model matrix, including the intercept term, that
#'   predicts the means of the latent dimensions. If not supplied, only the
#'   intercept term is included.
#' @param z An optional model matrix, including the intercept term as its
#'   first column, that predicts the covariance matrix of the latent
#'   dimensions. If not supplied, only the intercept term is included.
#' @param constr The type of constraints used to identify the scale of each
#'   dimension: "latent_scale" (default) sets the geometric mean of the prior
#'   variance of each dimension to one; "items" sets the geometric mean of the
#'   discrimination parameters within each dimension to one. In both cases the
#'   mean of the prior means of each dimension is set to zero.
#' @param beta_set The indices (in \code{y}) of one item per dimension whose
#'   discrimination parameter is restricted to be positive (or negative). By
#'   default, the first item of each dimension.
#' @param sign_set Logical. Should the discrimination parameters of the items
#'   in \code{beta_set} be positive (\code{TRUE}) or negative (\code{FALSE})?
#' @param control A list of control values
#' \describe{
#'  \item{max_iter}{The maximum number of EM iterations. The default is 500.}
#'  \item{eps}{Convergence tolerance: iterations stop when the largest absolute
#'    change in any parameter falls under \code{eps}. The default is 1e-4.}
#'  \item{adapt}{Logical. Use adaptive Gauss-Hermite quadrature (\code{TRUE})
#'    or a fixed Gauss-Legendre grid (\code{FALSE}, the default)?}
#'  \item{K}{Number of quadrature points per dimension. With a fixed grid the
#'    default is 25 for one dimension and 21 for two or three dimensions; with
#'    adaptive quadrature it is 11.}
#'  \item{C}{Fixed grid only: [-C, C] sets the range of integration on each
#'    dimension. The default is 4 for one dimension and 5 otherwise.}
#'  \item{verbose}{Logical. Should the progress of the EM algorithm be
#'    printed? \code{FALSE} by default.}
#' }
#'
#' @return An object of class \code{mhgrm}, a list with
#'  \item{coefficients}{A data frame of parameter estimates, standard errors,
#'   z values, and p values.}
#'  \item{vcov}{Variance-covariance matrix of parameter estimates.}
#'  \item{scores}{A list with the posterior means (\code{post_mean}) and
#'   standard deviations (\code{post_sd}) of the latent dimensions, their
#'   posterior correlations (\code{post_cor}, pairs in the order 1-2, 1-3,
#'   2-3), and the prior means (\code{prior_mean}).}
#'  \item{log_Lik}{The log-likelihood at convergence.}
#'  \item{alpha, beta, B, gamma, lambda}{Parameter estimates in matrix form.}
#'  \item{dims, dim_names, d, N, J, H, p, q, iter, control, call}{Model details.}
#' @seealso \code{\link{mhgrm_cov}}, \code{\link{hgrm}}
#' @export
#' @examples
#' \donttest{
#' # economic, civil rights, and moral attitudes in ANES 2008
#' y <- nes2008[, -(1:3)]
#' dims <- rep(c("econ", "civil", "moral"), c(10, 13, 9))
#' x <- model.matrix( ~ party * educ, nes2008)
#' z <- model.matrix( ~ party + educ, nes2008)
#' # anchor each dimension on an item whose higher codes are conservative
#' anchors <- match(c("health_ins7", "blacks_aid7", "women_role7"), names(y))
#' fit <- mhgrm(y, dims, x, z, beta_set = anchors)
#' fit
#' # variances, covariances, and correlations for Democrats and Republicans
#' # without a college degree
#' mhgrm_cov(fit, z = rbind(c(1, 0, 0, 0), c(1, 0, 1, 0)))
#' }
mhgrm <- function(y, dims, x = NULL, z = NULL, constr = c("latent_scale", "items"),
                  beta_set = NULL, sign_set = TRUE, control = list()) {

  cl <- match.call()
  constr <- match.arg(constr)

  # responses
  if ((!is.data.frame(y) && !is.matrix(y)) || ncol(y) < 2L)
    stop("'y' must be either a data.frame or a matrix with at least two columns.")
  if (is.matrix(y)) y <- as.data.frame(y)
  N <- nrow(y)
  J <- ncol(y)
  y[] <- lapply(y, factor, exclude = c(NA, NaN))
  ylevels <- lapply(y, levels)
  y[] <- lapply(y, as.integer)
  if (!is.na(invalid <- match(TRUE, vapply(y, invalid_grm, logical(1L)))))
    stop(paste(names(y)[invalid], "does not have at least two valid responses"))
  H <- vapply(y, max, integer(1L), na.rm = TRUE)
  Y <- y_mat(y)

  # dimensions
  if (length(dims) != J) stop("`dims` must have one entry per column of `y`.")
  dim_f <- if (is.factor(dims)) droplevels(dims) else factor(dims, levels = unique(dims))
  dim_names <- levels(dim_f)
  dim_idx <- as.integer(dim_f)
  d <- length(dim_names)
  if (d > 3L) stop("At most three dimensions are supported.")
  items <- lapply(seq_len(d), function(k) which(dim_idx == k))
  if (any(lengths(items) < 2L)) stop("Each dimension needs at least two items.")

  # covariates
  x <- x %||% as.matrix(rep(1, N))
  z <- z %||% as.matrix(rep(1, N))
  if (!is.matrix(x) || !is.matrix(z)) stop("`x` and `z` must be matrices.")
  if (nrow(x) != N || nrow(z) != N)
    stop("both 'x' and 'z' must have the same number of rows as 'y'")
  if (any(x[, 1] != 1) || any(z[, 1] != 1))
    stop("the first column of `x` and `z` must be the intercept (all ones).")
  x <- dbl_mat(x)
  z <- dbl_mat(z)
  p <- ncol(x)
  q <- ncol(z)
  colnames(x) <- fill_names(colnames(x), p, "x")
  colnames(z) <- fill_names(colnames(z), q, "z")

  # sign restrictions
  beta_set <- beta_set %||% vapply(items, `[`, integer(1L), 1L)
  if (length(beta_set) != d || any(dim_idx[beta_set] != seq_len(d)))
    stop("`beta_set` must give one item from each dimension, in dimension order.")
  stopifnot(is.logical(sign_set))

  # control
  adapt <- if (is.null(control[["adapt"]])) FALSE else control[["adapt"]]
  con <- list(max_iter = 500, eps = 1e-4, adapt = adapt,
              K = if (adapt) 11 else if (d == 1) 25 else 21,
              C = if (d == 1) 4 else 5, verbose = FALSE)
  con[names(control)] <- control
  K <- con[["K"]]
  if (adapt) {
    gh <- gh_rule(K)
    g <- gh$x
    lqw <- gh$lw
  } else {
    g <- con[["C"]] * GLpoints[[K]][["x"]]
    lqw <- log(con[["C"]] * GLpoints[[K]][["w"]])
  }

  # ---- initial values: separate unidimensional fits ----
  alpha <- vector("list", J)
  beta <- double(J)
  B <- matrix(0, p, d)
  eap <- psd <- matrix(0, N, d)
  for (k in seq_len(d)) {
    yk <- y[, items[[k]], drop = FALSE]
    bk <- match(beta_set[k], items[[k]])
    # these fits only supply starting values: if one does not converge (as with
    # weakly related items), retry with a loose tolerance before giving up
    fit_k <- tryCatch(
      hgrm(yk, x, beta_set = bk, sign_set = sign_set),
      error = function(e) tryCatch(
        hgrm(yk, x, beta_set = bk, sign_set = sign_set,
             control = list(eps = 1e-2, max_iter = 1000)),
        error = function(e2)
          stop("starting values for dimension '", dim_names[k], "' failed: ",
               conditionMessage(e2), call. = FALSE)))
    est <- fit_k$coefficients$Estimate
    pos <- 0
    for (j in items[[k]]) {
      a <- est[pos + seq_len(H[[j]])]
      alpha[[j]] <- c(Inf, a[-H[[j]]], -Inf)
      beta[j] <- a[H[[j]]]
      pos <- pos + H[[j]]
    }
    B[, k] <- est[pos + seq_len(p)]
    eap[, k] <- fit_k$scores$post_mean
    psd[, k] <- fit_k$scores$post_sd
  }
  names(alpha) <- names(beta) <- names(y)
  Sigma0 <- if (d == 1) matrix(1) else cor(eap - x %*% B)
  ch <- chol_params(Sigma0)
  gam <- matrix(0, q, d)
  gam[1, ] <- log(ch$sig2)
  lam <- matrix(0, q, 3)          # columns: lambda_21, lambda_31, lambda_32
  lam[1, ] <- ch$phi
  # node centers and scales: posterior means and SDs (adaptive), or a common grid
  ctr <- if (adapt) eap else matrix(0, N, d)
  scl <- if (adapt) psd else matrix(1, N, d)

  # ---- EM algorithm ----
  par_vec <- function() c(unlist(alpha[]) [is.finite(unlist(alpha))], beta, B, gam, lam)
  # Guards. On a fixed grid the level-II update is the M-step for a continuous
  # normal prior, not for its discretized version, so when prior distributions
  # drift off the grid (typically with weakly identified dimensions) the
  # log-likelihood can fall and the fit can run away. Brief dips are tolerated;
  # the fit stops, returning the estimates with the highest log-likelihood, if
  # the log-likelihood falls far below its maximum or stays below it for
  # `patience` iterations, or if the level-II update fails.
  ll_best <- -Inf
  best <- NULL
  stop_msg <- NULL
  patience <- 50L
  for (iter in seq_len(con[["max_iter"]])) {
    par_prev <- par_vec()

    # E-step (es$loglik is the log-likelihood of the current estimates)
    es <- estep_md(Y, items, alpha, beta, g, lqw, ctr, scl, x %*% B, exp(z %*% gam), z %*% lam)
    if (es$loglik > ll_best) {
      ll_best <- es$loglik
      best <- list(alpha = alpha, beta = beta, B = B, gam = gam, lam = lam,
                   ctr = ctr, scl = scl, iter = iter)
    } else if (!adapt) {
      fall <- ll_best - es$loglik
      if (fall > max(10, 1e-3 * abs(ll_best))) {
        stop_msg <- sprintf("the log-likelihood fell from %.2f (iteration %d) to %.2f (iteration %d)",
                            ll_best, best$iter, es$loglik, iter)
        break
      }
      if (iter - best$iter >= patience && fall > 1e-6 * abs(ll_best)) {
        stop_msg <- sprintf(paste("the log-likelihood has stayed below its maximum of %.2f",
                                  "(iteration %d) for %d iterations"),
                            ll_best, best$iter, patience)
        break
      }
    }

    # M-step: items
    for (k in seq_len(d)) {
      jk <- items[[k]]
      ab <- if (adapt) {
        mstep_grm_nodes(Y[, jk, drop = FALSE], es$marg[[k]], H[jk], es$G[[k]], alpha[jk], beta[jk])
      } else mstep_grm(Y[, jk, drop = FALSE], es$marg[[k]], H[jk], g, alpha[jk], beta[jk])
      alpha[jk] <- ab$alpha
      beta[jk] <- ab$beta
    }

    # M-step: level II (one ECM cycle)
    lv2 <- tryCatch(level2_update(es$m, post_cov(es), x, z, B, gam, lam, d),
                    error = function(e) e)
    if (inherits(lv2, "error")) {
      stop_msg <- sprintf("the update of the level-II parameters failed at iteration %d (%s)",
                          iter, conditionMessage(lv2))
      break
    }
    B <- lv2$B
    gam <- lv2$gam
    lam <- lv2$lam

    # identification: location, scale, and sign of each dimension
    nm <- normalize_md(alpha, beta, B, gam, lam, x, z, items, d, constr, beta_set, sign_set)
    alpha <- nm$alpha
    beta <- nm$beta
    B <- nm$B
    gam <- nm$gam
    lam <- nm$lam
    if (adapt) {
      # re-center the grids on the new posteriors, expressed on the normalized scale
      ctr <- sweep(sweep(es$m, 2, nm$shift), 2, nm$cvec, `*`)
      scl <- sweep(post_sd_md(es), 2, abs(nm$cvec), `*`)
    }

    if (con[["verbose"]]) cat(sprintf("iter %3d  log-lik %.4f\n", iter, es$loglik))
    if (max(abs(par_vec() - par_prev)) < con[["eps"]]) break
    if (iter == con[["max_iter"]])
      warning("algorithm did not converge; try increasing `max_iter`.")
  }
  grid_hint <- if (adapt) "" else " or a larger `C` (or `adapt = TRUE`)"
  if (!is.null(stop_msg)) {
    alpha <- best$alpha; beta <- best$beta; B <- best$B; gam <- best$gam; lam <- best$lam
    ctr <- best$ctr; scl <- best$scl
    warning(stop_msg, "; returning the estimates with the highest log-likelihood (iteration ",
            best$iter, "). The model may be weakly identified (e.g., dimensions with few or ",
            "weakly related items); try more items per dimension, a simpler `z`", grid_hint, ".",
            call. = FALSE)
  }

  # ---- final E-step and inference ----
  mu <- x %*% B
  sig2 <- exp(z %*% gam)
  phi <- z %*% lam
  es <- estep_md(Y, items, alpha, beta, g, lqw, ctr, scl, mu, sig2, phi)
  if (adapt) es <- estep_md(Y, items, alpha, beta, g, lqw, es$m, post_sd_md(es), mu, sig2, phi)
  if (!adapt && is.null(stop_msg) && es$loglik < ll_best - 1e-6 * abs(ll_best)) {
    # the final estimates are worse than earlier ones: return the best
    warning(sprintf(paste("the log-likelihood was not monotone on the fixed grid; returning the",
                          "estimates with the highest log-likelihood (%.2f at iteration %d,",
                          "against %.2f at the last iteration)."), ll_best, best$iter, es$loglik),
            call. = FALSE)
    alpha <- best$alpha; beta <- best$beta; B <- best$B; gam <- best$gam; lam <- best$lam
    mu <- x %*% B
    sig2 <- exp(z %*% gam)
    phi <- z %*% lam
    es <- estep_md(Y, items, alpha, beta, g, lqw, ctr, scl, mu, sig2, phi)
  }
  if (!adapt) {
    # does the fixed grid cover the prior distributions?
    S <- sigma_entries(sig2, phi, d)
    out_mass <- max(vapply(seq_len(d), function(k) {
      s <- sqrt(S[[k]][[k]])
      max(pnorm(-con[["C"]], mu[, k], s) + pnorm(con[["C"]], mu[, k], s, lower.tail = FALSE))
    }, double(1L)))
    if (out_mass > 0.01)
      warning(sprintf(paste0("the fixed grid [-%g, %g] leaves up to %.1f%% of some respondents' ",
                             "prior distribution outside its range, so the estimates may be ",
                             "inaccurate; refit with a larger `C` or `adapt = TRUE`."),
                      con[["C"]], con[["C"]], 100 * out_mass), call. = FALSE)
  }
  Spost <- post_cov(es)

  # item scores (rows: a_j2..a_jH, b_j for each item, in the order of y)
  s_item <- matrix(0, sum(H), N)
  row_start <- cumsum(c(0, H))[seq_len(J)]
  for (k in seq_len(d)) {
    jk <- items[[k]]
    sk <- item_scores_md_cpp(Y[, jk, drop = FALSE], alpha_mat(alpha[jk]), as.double(beta[jk]),
                             es$G[[k]], es$marg[[k]], as.integer(H[jk]))
    rows <- unlist(lapply(jk, function(j) row_start[j] + seq_len(H[[j]])))
    s_item[rows, ] <- sk
  }
  s_lv2 <- level2_scores(es$m, Spost, x, z, mu, sig2, phi, d)

  # Outer-product information for all parameters. It is singular along the
  # 2d directions of location and scale invariance, so the covariance matrix
  # of the constrained estimates is the top-left block of the inverse of the
  # information bordered by the Jacobian of the identification constraints.
  s_all <- rbind(s_item, s_lv2$s)
  s_all[is.na(s_all)] <- 0
  Cmat <- constraint_jacobian(beta, B, gam, lam, x, z, items, H, d, constr)
  n_par <- nrow(s_all)
  bordered <- rbind(cbind(tcrossprod(s_all), t(Cmat)),
                    cbind(Cmat, matrix(0, nrow(Cmat), nrow(Cmat))))
  covmat <- tryCatch(solve(bordered)[seq_len(n_par), seq_len(n_par)],
                     error = function(e) {
                       warning("The information matrix is singular; SE calculation failed.")
                       matrix(NA_real_, n_par, n_par)
                     })
  se_all <- sqrt(pmax(diag(covmat), 0))
  se_item <- se_all[seq_len(sum(H))]
  se_lv2 <- se_all[-seq_len(sum(H))]

  # ---- assemble output ----
  item_names <- unlist(lapply(seq_len(J), function(j)
    paste(names(y)[j], c(if (H[[j]] == 2) "Intercept" else paste0("y>=", seq(2, H[[j]])),
                         "Dscrmn"))))
  est_item <- unlist(lapply(seq_len(J), function(j) c(alpha[[j]][2:H[[j]]], beta[[j]])))
  pairs <- lam_pairs(d)
  names_lv2 <- c(
    paste0("mean[", rep(dim_names, each = p), "]: ", rep(colnames(x), d)),
    paste0("logvar[", rep(dim_names, each = q), "]: ", rep(colnames(z), d)),
    if (d > 1) paste0("chol[", rep(dim_names[pairs[, 1]], each = q), " ~ ",
                      rep(dim_names[pairs[, 2]], each = q), "]: ",
                      rep(colnames(z), nrow(pairs))))
  est_lv2 <- c(B, gam, lam[, pairs[, "col"], drop = FALSE])
  est_all <- c(est_item, est_lv2)
  se_all <- c(se_item, se_lv2)
  coefs <- data.frame(Estimate = est_all, Std_Error = se_all, z_value = est_all/se_all,
                      p_value = 2 * (1 - pnorm(abs(est_all/se_all))),
                      row.names = c(item_names, names_lv2))
  dimnames(covmat) <- list(c(item_names, names_lv2), c(item_names, names_lv2))

  post_sd <- sqrt(vapply(seq_len(d), function(k) Spost[[k]][[k]], double(N)))
  post_sd <- matrix(post_sd, N, d)
  post_cor <- if (d > 1) {
    vapply(seq_len(nrow(pairs)), function(r) {
      a <- pairs[r, 2]; b <- pairs[r, 1]
      Spost[[a]][[b]] / (post_sd[, a] * post_sd[, b])
    }, double(N))
  } else NULL
  colnames(es$m) <- colnames(post_sd) <- colnames(mu) <- dim_names
  if (!is.null(post_cor))
    colnames(post_cor) <- paste(dim_names[pairs[, 2]], dim_names[pairs[, 1]], sep = "-")
  rownames(B) <- colnames(x)
  rownames(gam) <- rownames(lam) <- colnames(z)
  colnames(B) <- colnames(gam) <- dim_names
  colnames(lam) <- c("21", "31", "32")

  out <- list(coefficients = coefs, vcov = covmat,
              scores = list(post_mean = es$m, post_sd = post_sd, post_cor = post_cor,
                            prior_mean = mu),
              log_Lik = es$loglik, alpha = alpha, beta = beta, B = B, gamma = gam,
              lambda = lam[, pairs[, "col"], drop = FALSE],
              dims = dim_idx, dim_names = dim_names, d = d, constr = constr,
              N = N, J = J, H = H, ylevels = ylevels, p = p, q = q, iter = iter,
              z_names = colnames(z), control = con, call = cl)
  class(out) <- "mhgrm"
  out
}

#' Implied Covariance Structure of a Multidimensional Hierarchical Model
#'
#' Computes the prior variances, covariances, and correlations of the latent
#' dimensions implied by a fitted \code{\link{mhgrm}} model at given values of
#' the covariance predictors, with delta-method standard errors.
#'
#' @param object A fitted \code{mhgrm} object.
#' @param z A matrix of covariance predictors (with the intercept as its first
#'   column and the same columns as in the fitted model). Each row gives one
#'   covariate profile.
#' @return A data frame with one row per covariate profile and quantity:
#'   \code{profile} (row of \code{z}), \code{quantity} (e.g. \code{"var[econ]"},
#'   \code{"cov[econ, social]"}, \code{"cor[econ, social]"}), \code{estimate},
#'   and \code{std_error}.
#' @export
mhgrm_cov <- function(object, z) {
  if (!inherits(object, "mhgrm")) stop("Use only with 'mhgrm' objects.")
  z <- dbl_mat(rbind(z))
  if (ncol(z) != object$q) stop("`z` must have ", object$q, " columns.")
  d <- object$d
  q <- object$q
  dn <- object$dim_names
  pairs <- lam_pairs(d)
  # parameters of the covariance model and their covariance matrix
  theta0 <- c(object$gamma, object$lambda)
  idx <- match(c(paste0("logvar[", rep(dn, each = q), "]: ", rep(object$z_names, d)),
                 if (d > 1) paste0("chol[", rep(dn[pairs[, 1]], each = q), " ~ ",
                                   rep(dn[pairs[, 2]], each = q), "]: ",
                                   rep(object$z_names, nrow(pairs)))),
               rownames(object$vcov))
  V <- as.matrix(object$vcov)[idx, idx, drop = FALSE]
  labels <- c(paste0("var[", dn, "]"),
              if (d > 1) c(paste0("cov[", dn[pairs[, 2]], ", ", dn[pairs[, 1]], "]"),
                           paste0("cor[", dn[pairs[, 2]], ", ", dn[pairs[, 1]], "]")))
  f <- function(th, zr) {
    gam <- matrix(th[seq_len(q * d)], q, d)
    lam <- matrix(0, q, 3)
    if (d > 1) lam[, pairs[, "col"]] <- th[-seq_len(q * d)]
    S <- sigma_entries(exp(zr %*% gam), zr %*% lam, d)
    vars <- vapply(seq_len(d), function(k) S[[k]][[k]], double(1L))
    if (d == 1) return(vars)
    covs <- vapply(seq_len(nrow(pairs)), function(r) S[[pairs[r, 2]]][[pairs[r, 1]]], double(1L))
    cors <- covs / sqrt(vars[pairs[, 2]] * vars[pairs[, 1]])
    c(vars, covs, cors)
  }
  out <- lapply(seq_len(nrow(z)), function(r) {
    zr <- z[r, , drop = FALSE]
    est <- f(theta0, zr)
    # numerical Jacobian (central differences)
    Jm <- vapply(seq_along(theta0), function(t) {
      h <- 1e-6 * max(1, abs(theta0[t]))
      e <- replace(numeric(length(theta0)), t, h)
      (f(theta0 + e, zr) - f(theta0 - e, zr)) / (2 * h)
    }, double(length(est)))
    Jm <- matrix(Jm, length(est))
    se <- sqrt(pmax(diag(Jm %*% V %*% t(Jm)), 0))
    data.frame(profile = r, quantity = labels, estimate = est, std_error = se)
  })
  do.call(rbind, out)
}

#' @export
print.mhgrm <- function(x, digits = 3, ...) {
  cat("\nCall:\n", paste(deparse(x[["call"]]), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Dimensions:", paste0(x$dim_names, " (", tabulate(x$dims, x$d), " items)",
                            collapse = ", "), "\n\n")
  cf <- x$coefficients
  lv2 <- !grepl("( y>=| Intercept$| Dscrmn$)", rownames(cf))
  cat("Level-II coefficients:\n")
  print(round(cf[lv2, , drop = FALSE], digits), ...)
  cat("\nLog Likelihood: ", round(x$log_Lik, digits), " (", x$iter, " EM iterations)\n\n", sep = "")
  invisible(x)
}

# ---- internal helpers ----

# column names for a design matrix whose first column is the intercept;
# missing or empty names are filled in
fill_names <- function(nms, n, prefix) {
  default <- c("(Intercept)", paste0(prefix, seq_len(n - 1)))[seq_len(n)]
  if (is.null(nms)) return(default)
  ifelse(is.na(nms) | nms == "", default, nms)
}

# rows of the lambda matrix: (2,1) -> col 1, (3,1) -> col 2, (3,2) -> col 3
lam_pairs <- function(d) {
  all <- cbind(k = c(2L, 3L, 3L), l = c(1L, 1L, 2L), col = 1:3)
  all[all[, "k"] <= d, , drop = FALSE]
}

# modified Cholesky parameters of a fixed covariance matrix
chol_params <- function(S) {
  d <- nrow(S)
  Lc <- t(chol(S))                       # S = Lc Lc'
  Dh <- diag(Lc)
  Lu <- Lc %*% diag(1 / Dh, d)           # unit lower triangular, S = Lu D Lu'
  Tm <- solve(Lu)                        # T S T' = D
  phi <- c(if (d >= 2) -Tm[2, 1] else 0, if (d >= 3) -Tm[3, 1] else 0,
           if (d >= 3) -Tm[3, 2] else 0)
  list(sig2 = Dh^2, phi = phi)
}

# entries of Sigma_i = L D L' (L = T^{-1}) as a nested list of N-vectors
sigma_entries <- function(sig2, phi, d) {
  N <- nrow(sig2)
  one <- rep(1, N); zero <- rep(0, N)
  Lm <- list(list(one, zero, zero))
  if (d >= 2) Lm[[2]] <- list(phi[, 1], one, zero)
  if (d >= 3) Lm[[3]] <- list(phi[, 2] + phi[, 3] * phi[, 1], phi[, 3], one)
  S <- vector("list", d)
  for (a in seq_len(d)) {
    S[[a]] <- vector("list", d)
    for (b in seq_len(d)) {
      v <- 0
      for (k in seq_len(d)) v <- v + Lm[[a]][[k]] * Lm[[b]][[k]] * sig2[, k]
      S[[a]][[b]] <- v
    }
  }
  S
}

# entries of the precision matrix T' D^{-1} T
prec_entries <- function(sig2, phi, d) {
  N <- nrow(sig2)
  one <- rep(1, N); zero <- rep(0, N)
  Tm <- list(list(one, zero, zero))
  if (d >= 2) Tm[[2]] <- list(-phi[, 1], one, zero)
  if (d >= 3) Tm[[3]] <- list(-phi[, 2], -phi[, 3], one)
  P <- vector("list", d)
  for (a in seq_len(d)) {
    P[[a]] <- vector("list", d)
    for (b in seq_len(d)) {
      v <- 0
      for (k in seq_len(d)) v <- v + Tm[[k]][[a]] * Tm[[k]][[b]] / sig2[, k]
      P[[a]][[b]] <- v
    }
  }
  P
}

# Gauss-Hermite rule for the standard normal (Golub-Welsch), returned as nodes x
# and log weights lw for integrating against Lebesgue measure:
# int f(t) dt ~ sum_k exp(lw_k) f(x_k)
gh_rule <- function(K) {
  i <- seq_len(K - 1)
  Jm <- matrix(0, K, K)
  Jm[cbind(i, i + 1)] <- Jm[cbind(i + 1, i)] <- sqrt(i / 2)
  e <- eigen(Jm, symmetric = TRUE)
  x <- rev(e$values) * sqrt(2)
  w <- rev(e$vectors[1, ]^2)
  list(x = x, lw = log(w) - dnorm(x, log = TRUE))
}

# E-step wrapper. On dimension a, person i's nodes are ctr[i, a] + scl[i, a] * g,
# stored as K x N matrices (returned as es$G).
estep_md <- function(Y, items, alpha, beta, g, lqw, ctr, scl, mu, sig2, phi) {
  d <- length(items)
  K <- length(g)
  G <- lapply(seq_len(d), function(a) outer(g, scl[, a]) + rep(ctr[, a], each = K))
  LL <- lapply(seq_len(d), function(a) {
    jk <- items[[a]]
    loglik_nodes_cpp(Y[, jk, drop = FALSE], alpha_mat(alpha[jk]), as.double(beta[jk]), G[[a]])
  })
  phi3 <- matrix(0, nrow(mu), 3)
  phi3[, seq_len(ncol(phi))] <- phi
  es <- estep_md_cpp(LL, dbl_mat(mu), dbl_mat(sig2), phi3, G, dbl_mat(log(scl)), lqw)
  es$G <- G
  es
}

# posterior SDs (N x d) from an E-step result
post_sd_md <- function(es) {
  d <- ncol(es$m)
  diag_col <- c(1, 3, 6)[seq_len(d)]
  sqrt(pmax(es$S2[, diag_col, drop = FALSE] - es$m^2, 1e-12))
}

# item M-step on person-specific nodes: one cumulative logit per item, fitted
# to the N x K expanded data (response y_ij at node G[k, i], weight w[k, i])
mstep_grm_nodes <- function(Y, w, H, G, alpha, beta) {
  K <- nrow(G)
  xg <- cbind(x = as.vector(G))
  wv <- as.vector(w)
  fits <- lapply(seq_along(H), function(j) {
    Hj <- H[[j]]
    cumlogit_fit(rep(Y[, j], each = K), xg, wv, Hj, c(alpha[[j]][2:Hj], beta[[j]]))
  })
  beta <- setNames(vapply(fits, function(f) f[[length(f)]], double(1L)), names(H))
  alpha <- setNames(lapply(fits, function(f) c(Inf, f[-length(f)], -Inf)), names(H))
  list(alpha = alpha, beta = beta)
}

# posterior covariance entries from an E-step result (nested list of N-vectors)
post_cov <- function(es) {
  d <- ncol(es$m)
  pair_col <- matrix(c(1, 2, 4, 2, 3, 5, 4, 5, 6), 3)[seq_len(d), seq_len(d), drop = FALSE]
  lapply(seq_len(d), function(a) lapply(seq_len(d), function(b)
    es$S2[, pair_col[a, b]] - es$m[, a] * es$m[, b]))
}

# one conditional-maximization cycle for the level-II parameters
level2_update <- function(m, Spost, x, z, B, gam, lam, d) {
  p <- ncol(x); q <- ncol(z)
  # (a) means by GLS given the current covariance matrices
  P <- prec_entries(exp(z %*% gam), z %*% lam, d)
  A <- matrix(0, p * d, p * d)
  rhs <- double(p * d)
  for (a in seq_len(d)) {
    ra <- (a - 1) * p + seq_len(p)
    for (b in seq_len(d)) {
      A[ra, (b - 1) * p + seq_len(p)] <- crossprod(x, P[[a]][[b]] * x)
      rhs[ra] <- rhs[ra] + crossprod(x, P[[a]][[b]] * m[, b])
    }
  }
  B <- matrix(solve(A, rhs), p, d)
  # (b) expected outer products of eps = theta - mu
  R <- m - x %*% B
  M <- lapply(seq_len(d), function(a) lapply(seq_len(d), function(b)
    Spost[[a]][[b]] + R[, a] * R[, b]))
  # (c) Cholesky rows: lambda by WLS, then gamma by a gamma GLM
  pairs <- lam_pairs(d)
  for (k in seq_len(d)) {
    s2 <- exp(drop(z %*% gam[, k]))
    if (k > 1) {
      L <- seq_len(k - 1)
      cols <- pairs[pairs[, "k"] == k, "col"]
      Ak <- matrix(0, q * (k - 1), q * (k - 1))
      bk <- double(q * (k - 1))
      for (l in L) {
        rl <- (l - 1) * q + seq_len(q)
        for (l2 in L) Ak[rl, (l2 - 1) * q + seq_len(q)] <- crossprod(z, M[[l]][[l2]] / s2 * z)
        bk[rl] <- crossprod(z, M[[l]][[k]] / s2)
      }
      lam[, cols] <- solve(Ak, bk)
      phik <- z %*% lam[, cols, drop = FALSE]
      qk <- M[[k]][[k]]
      for (l in L) {
        qk <- qk - 2 * phik[, l] * M[[l]][[k]]
        for (l2 in L) qk <- qk + phik[, l] * phik[, l2] * M[[l]][[l2]]
      }
    } else qk <- M[[1]][[1]]
    gam[, k] <- if (q == 1) log(mean(qk)) else
      as.double(gamma_glm_cpp(z, qk)[["coefficients"]])
  }
  list(B = B, gam = gam, lam = lam)
}

# scores of the level-II parameters (Fisher identity): rows B, gamma, lambda
level2_scores <- function(m, Spost, x, z, mu, sig2, phi, d) {
  N <- nrow(m); p <- ncol(x); q <- ncol(z)
  R <- m - mu
  P <- prec_entries(sig2, phi, d)
  sB <- do.call(rbind, lapply(seq_len(d), function(a) {
    PR <- Reduce(`+`, lapply(seq_len(d), function(b) P[[a]][[b]] * R[, b]))
    t(PR * x)
  }))
  M <- lapply(seq_len(d), function(a) lapply(seq_len(d), function(b)
    Spost[[a]][[b]] + R[, a] * R[, b]))
  pairs <- lam_pairs(d)
  sG <- NULL; sL <- vector("list", 3)
  for (k in seq_len(d)) {
    L <- seq_len(k - 1)
    cols <- pairs[pairs[, "k"] == k, "col"]
    # E[e_k eps_l] for l < k and E[e_k^2]
    Eke <- lapply(seq_len(d), function(l) {
      v <- M[[k]][[l]]
      for (j in L) v <- v - phi[, cols[j]] * M[[j]][[l]]
      v
    })
    qk <- Eke[[k]]
    for (j in L) qk <- qk - phi[, cols[j]] * Eke[[j]]
    sG <- rbind(sG, t(0.5 * (qk / sig2[, k] - 1) * z))
    for (j in L) sL[[cols[j]]] <- t(Eke[[j]] / sig2[, k] * z)
  }
  list(s = rbind(sB, sG, do.call(rbind, sL[pairs[, "col"]])))
}

# Jacobian (2d x n_par) of the identification constraints, in the parameter
# order (item thresholds and discriminations, B, gamma, lambda):
#   location_k: mean_i mu_ik = 0
#   scale_k   : mean_i log Var(theta_ik) = 0 ("latent_scale"), or
#               mean_{j in k} log|beta_j| = 0 ("items")
constraint_jacobian <- function(beta, B, gam, lam, x, z, items, H, d, constr) {
  p <- ncol(x); q <- ncol(z)
  pairs <- lam_pairs(d)
  n_item <- sum(H)
  n_par <- n_item + p * d + q * d + q * nrow(pairs)
  row_start <- cumsum(c(0, H))[seq_along(H)]
  Cm <- matrix(0, 2 * d, n_par)
  xbar <- colMeans(x)
  for (k in seq_len(d)) Cm[k, n_item + (k - 1) * p + seq_len(p)] <- xbar
  if (constr == "items") {
    for (k in seq_len(d)) for (j in items[[k]])
      Cm[d + k, row_start[j] + H[[j]]] <- 1 / (length(items[[k]]) * beta[[j]])
  } else {
    th0 <- c(gam, lam[, pairs[, "col"]])
    f <- function(th) {
      g2 <- matrix(th[seq_len(q * d)], q, d)
      l2 <- matrix(0, q, 3)
      if (d > 1) l2[, pairs[, "col"]] <- th[-seq_len(q * d)]
      S <- sigma_entries(exp(z %*% g2), z %*% l2, d)
      vapply(seq_len(d), function(k) mean(log(S[[k]][[k]])), double(1L))
    }
    Jm <- vapply(seq_along(th0), function(t) {
      h <- 1e-6 * max(1, abs(th0[t]))
      e <- replace(numeric(length(th0)), t, h)
      (f(th0 + e) - f(th0 - e)) / (2 * h)
    }, double(d))
    Cm[d + seq_len(d), n_item + p * d + seq_along(th0)] <- matrix(Jm, d)
  }
  Cm
}

# move along the invariance orbit so that each dimension satisfies the
# location, scale, and sign constraints
normalize_md <- function(alpha, beta, B, gam, lam, x, z, items, d, constr, beta_set, sign_set) {
  mu <- x %*% B
  cvec <- shift <- double(d)
  for (k in seq_len(d)) {
    # location: mean of prior means is zero
    shift[k] <- s <- mean(mu[, k])
    for (j in items[[k]]) alpha[[j]] <- alpha[[j]] + beta[[j]] * s
    B[1, k] <- B[1, k] - s
  }
  if (constr == "latent_scale") {
    S <- sigma_entries(exp(z %*% gam), z %*% lam, d)
    for (k in seq_len(d)) cvec[k] <- exp(-mean(log(S[[k]][[k]])) / 2)
  } else {
    for (k in seq_len(d)) cvec[k] <- exp(mean(log(abs(beta[items[[k]]]))))
  }
  for (k in seq_len(d))
    if (sign_set == (beta[[beta_set[k]]] < 0)) cvec[k] <- -cvec[k]
  # theta_k -> c_k theta_k
  for (k in seq_len(d)) {
    beta[items[[k]]] <- beta[items[[k]]] / cvec[k]
    B[, k] <- B[, k] * cvec[k]
    gam[1, k] <- gam[1, k] + 2 * log(abs(cvec[k]))
  }
  pairs <- lam_pairs(d)
  for (r in seq_len(nrow(pairs)))
    lam[, pairs[r, "col"]] <- lam[, pairs[r, "col"]] * cvec[pairs[r, "k"]] / cvec[pairs[r, "l"]]
  # theta_k -> c_k (theta_k - shift_k)
  list(alpha = alpha, beta = beta, B = B, gam = gam, lam = lam, shift = shift, cvec = cvec)
}
