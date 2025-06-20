#' Parameter Estimates from Hierarchical IRT Models.
#'
#' Parameter estimates from either \code{hltm} or \code{hgrm} or \code{hgrmDIF} models. \code{code_item}
#' reports estimates of item parameters. \code{coef_mean} reports results for the mean equation.
#' \code{coef_var} reports results for the variance equation.
#'
#' @inheritParams print.hIRT
#' @param by_item Logical. Should item parameters be stored item by item
#'   (if \code{TRUE}) or put together in a data frame (if \code{FALSE})?
#'
#' @return Parameter estimates, standard errors, z values, and p values
#'   organized as a data frame (if \code{by_item = TRUE}) or a list (if \code{
#'   by_item = FALSE}).
#'
#' @export
#' @examples
#' y <- nes_econ2008[, -(1:3)]
#' x <- model.matrix( ~ party * educ, nes_econ2008)
#' z <- model.matrix( ~ party, nes_econ2008)
#' nes_m1 <- hgrm(y, x, z)
#' coef_item(nes_m1)
coef_item <- function(x, by_item = TRUE, digits = 3) {

  if(inherits(x, "hgrm")){
    H <- unname(x[["H"]])
    xitem <- x[["coefficients"]][1:sum(H), , drop = FALSE]
    if (by_item == FALSE) return(round(xitem, digits))
    index <- findInterval(1:sum(H), c(1, cumsum(H)[-length(H)] + 1))
    out <- split(xitem, index)
    for (i in seq_along(out)) {
      tmp <- strsplit(rownames(out[[i]]), "\\.")
      rownames(out[[i]]) <- vapply(tmp, function(x) x[length(x)], FUN.VALUE = character(1L))
      out[[i]] <- round(out[[i]], digits)
    }
  } else if (inherits(x, "hltm")){
    J <- x[["J"]]
    xitem <- x[["coefficients"]][1:(2 * J), , drop = FALSE]
    if (by_item == FALSE)
      return(round(xitem, digits))
    out <- split(xitem, rep(1:J, each = 2))
    for (i in 1:J) rownames(out[[i]]) <- c("Diff", "Dscrmn")
  } else if (inherits(x, "hgrmDIF")){
    H <- unname(x[["H"]])
    p <- x[["p"]]
    ncoefs <- vapply(x[["coef_item"]], length, integer(1L))
    sH <- sum(ncoefs)
    xitem <- x[["coefficients"]][1:sH, , drop = FALSE]
    if (by_item == FALSE) return(round(xitem, digits))
    index <- findInterval(1:sH, c(1, cumsum(ncoefs[-length(ncoefs)]) + 1))
    out <- split(xitem, index)
    for (i in seq_along(out)) {
      tmp <- strsplit(rownames(out[[i]]), "\\.")
      rownames(out[[i]]) <- vapply(tmp, function(x) x[length(x)], FUN.VALUE = character(1L))
      out[[i]] <- round(out[[i]], digits)
    }
  } else stop("Use only with 'hgrm' or 'hltm' or `hgrmDIF` objects.\n")

  stats::setNames(out, names(x[["H"]]))
}

#' @inheritParams print.hIRT
#'
#' @export
#' @rdname coef_item
#' @examples
#' coef_mean(nes_m1)
coef_mean <- function(x, digits = 3) {

  # if (x[["p"]] < 2) return(NULL)
  sH <- if (inherits(x, "hltm"))
    2 * x[["J"]] else if (inherits(x, "hgrm"))
      sum(x[["H"]]) else if (inherits(x, "hgrmDIF"))
        sum(vapply(x[["coef_item"]], length, integer(1L))) else{
          stop("Use only with 'hgrm' or 'hltm' or `hgrmDIF` objects.\n")
        }
  gamma_indices <- (sH + 1):(sH + x[["p"]])
  round(x[["coefficients"]][gamma_indices, , drop = FALSE], digits = digits)
}

#' @inheritParams print.hIRT
#'
#' @export
#' @rdname coef_item
#' @examples
#' coef_var(nes_m1)
coef_var <- function(x, digits = 3) {

  # if (x[["q"]] < 2) return(NULL)
  sH <- if (inherits(x, "hltm"))
    2 * x[["J"]] else if (inherits(x, "hgrm"))
      sum(x[["H"]]) else if (inherits(x, "hgrmDIF"))
        sum(vapply(x[["coef_item"]], length, integer(1L))) else{
          stop("Use only with 'hgrm' or 'hltm' or `hgrmDIF` objects.\n")
        }
  lambda_indices <- (sH + x[["p"]] + 1):(sH + x[["p"]] + x[["q"]])
  round(x[["coefficients"]][lambda_indices, , drop = FALSE], digits = digits)
}
