#' Printing an object of class \code{hIRT}
#' @param x An object of class \code{hIRT}
#' @param digits The number of significant digits to use when printing
#' @param ... further arguments passed to \code{\link{print}}.
#' @return The input object \code{x}, returned invisibly. The function is called
#'   for its side effect of printing the call, the mean and variance regression
#'   coefficients, and the log likelihood of the fitted model.
#' @export
print.hIRT <- function(x, digits = 3, ...) {
  cat("\nCall:\n", paste(deparse(x[["call"]]), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Mean Regression:\n")
  print(coef_mean(x, digits), ...)
  cat("\n")
  cat("Variance Regression:\n")
  print(coef_var(x, digits), ...)
  cat("\nLog Likelihood:", round(x[["log_Lik"]], digits))
  cat("\n\n")
  invisible(x)
}
