#' Summarizing Hierarchical Item Response Theory Models
#'
#' Summarizing the fit of either \code{hltm} or \code{hgrm}.
#'
#' @param object An object of class \code{hIRT}.
#' @param by_item Logical. Should item parameters be stored item by item
#'   (if \code{TRUE}) or put together in a data frame (if \code{FALSE})?
#' @param digits the number of significant digits to use when printing.
#'
#' @return An object of class \code{summary_hIRT}.
#'  \item{call}{The matched call.}
#'  \item{model}{Model fit statistics: Log likelihood, AIC, and BIC.}
#'  \item{item_coefs}{Item parameter estimates, standard errors,
#'   z values, and p values.}
#'  \item{mean_coefs}{Parameter estimates for the mean equation.}
#'  \item{var_coefs}{Parameter estimates for the variance equation.}
#'
#' @export
#' @examples
#' y <- nes_econ2008[, -(1:3)]
#' x <- model.matrix( ~ party * educ, nes_econ2008)
#' z <- model.matrix( ~ party, nes_econ2008)
#' nes_m1 <- hgrm(y, x, z)
#' summary(nes_m1, by_item = TRUE)

summary.hIRT <- function(object, by_item = FALSE, digits = 3, ...) {
  item_coefs <- coef_item(object, by_item = by_item, digits)
  mean_coefs <- coef_mean(object, digits)
  var_coefs <- coef_var(object, digits)
  log_Lik = object[["log_Lik"]]
  df <- sum(object[["H"]]) + sum(object[["p"]]) + sum(object[["q"]]) - 2
  N <- length(object[["scores"]])
  AIC <- -2 * log_Lik + 2 * df
  BIC <- -2 * log_Lik + df * log(N)
  model <- list(log_Lik = log_Lik, AIC = AIC, BIC = BIC)
  out <- list(call = object[["call"]], model = model, item_coefs = item_coefs,
              mean_coefs = mean_coefs, var_coefs = var_coefs)
  class(out) <- c("summary_hIRT")
  out
}

#' @inheritParams print.hIRT
#' @export
#' @rdname summary.hIRT
print.summary_hIRT <- function(x, digits = 3, ...) {
  cat("\nCall:\n", paste(deparse(x[["call"]]), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Model Summary:\n")
  model_sum <- as.data.frame(x[["model"]], row.names = "")
  print(model_sum)
  cat("\n Item Coefficients: \n")
  print(x[["item_coefs"]], ...)
  cat("\n Mean Regression Coefficients:\n")
  print(x[["mean_coefs"]], ...)
  cat("\n Variance Regression Coefficients:\n")
  print(x[["var_coefs"]], ...)
  cat("\n\n")
  invisible(x)
}
