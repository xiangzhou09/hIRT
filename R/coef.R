#' Extracting Parameter Estimates from Hierarchical IRT Models.
#'
#' Parameter estimates from either \code{hltm} or \code{hgrm} models. \code{coef_mean}
#' reports results only for the mean equation. \code{coef_var} reports results only for
#' the variance equation.
#'
#' @param object An object of class \code{hIRT}
#' @param digits The number of significant digits to use when printing
#'
#' @return A data frame of parameter estimates, standard errors, z values and p values.
#'
#' @export
#'
#' @examples
#' y <- nes_econ2008[, -(1:3)]
#' x <- model.matrix( ~ party * educ, nes_econ2008)
#' z <- model.matrix( ~ party, nes_econ2008)
#' nes_m1 <- hgrm(y, x, z)
#' coef(nes_m1)
coef.hIRT <- function(object, digits = 3, ...) {
  round(object[["coefficients"]], digits)
}

#' @inheritParams print.hIRT
#'
#' @export
#' @rdname coef.hIRT
#' @examples
#' coef_mean(nes_m1)
coef_mean <- function(x, digits = 3) {
  if (!inherits(x, "hIRT"))
    stop("Use only with 'hIRT' objects.\n")
  if (x[["p"]] < 2) return(NULL)
  sH <- if (inherits(x, "hltm"))
    2 * x[["J"]] else sum(x[["H"]])
  gamma_indices <- (sH + 1):(sH + x[["p"]])
  round(x[["coefficients"]][gamma_indices, , drop = FALSE], digits = digits)
}

#' @inheritParams print.hIRT
#'
#' @export
#' @rdname coef.hIRT
#' @examples
#' coef_var(nes_m1)
coef_var <- function(x, digits = 3) {
  if (!inherits(x, "hIRT"))
    stop("Use only with 'hIRT' objects.\n")
  if (x[["q"]] < 2)return(NULL)
  sH <- if (inherits(x, "hltm"))
    2 * x[["J"]] else sum(x[["H"]])
  lambda_indices <- (sH + x[["p"]] + 1):(sH + x[["p"]] + x[["q"]])
  round(x[["coefficients"]][lambda_indices, , drop = FALSE], digits = digits)
}

#' Extracting Estimates of Item Parameters from Hierarchical IRT Models.
#'
#' Item parameter estimates from either \code{hltm} or \code{hgrm} models.
#'
#' @inheritParams print.hIRT
#' @param by_item Logical. Should item parameters be stored item by item
#'   (if \code{TRUE}) or put together in a data frame (if \code{FALSE})?
#'
#' @return Item parameter estimates, standard errors, z values, and p values
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

coef_item <- function(x, ...) UseMethod("coef_item")

#' @export
#' @rdname coef_item
coef_item.hgrm <- function(x, by_item = TRUE, digits = 3) {
  H <- unname(x[["H"]])
  index <- findInterval(1:sum(H), c(1, cumsum(H)[-length(H)] + 1))
  xitem <- x[["coefficients"]][1:sum(H), , drop = FALSE]
  if (by_item == FALSE) return(round(xitem, digits))
  out <- split(xitem, index)
  for (i in seq_along(out)) {
    tmp <- strsplit(rownames(out[[i]]), " ")
    rownames(out[[i]]) <- vapply(tmp, function(x) x[length(x)], FUN.VALUE = character(1L))
    out[[i]] <- round(out[[i]], digits)
  }
  stats::setNames(out, x[["item_names"]])
}

#' @export
#' @rdname coef_item
coef_item.hltm <- function(x, by_item = TRUE, digits = 3) {
  J <- x[["J"]]
  xitem <- x[["coefficients"]][1:(2 * J), , drop = FALSE]
  if (by_item == FALSE)
    return(round(xitem, digits))
  out <- split(xitem, rep(1:J, each = 2))
  for (i in 1:J) rownames(out[[i]]) <- c("Diff", "Dscrmn")
  stats::setNames(out, x[["item_names"]])
}
