#' Printing an object of class \code{hIRT}
#' @param x An object of class \code{hIRT}
#' @param digits the number of significant digits to use when printing
#' @param ... further arguments passed to or from other methods
#' @export
print.hIRT <- function(x, digits = 3, ...) {
  cat("\nCall:\n", paste(deparse(x[["call"]]), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Mean Regression:\n")
  print(coef_mean(x, digits), print.gap = 2, quote = FALSE)
  cat("\n")
  cat("Variance Regression:\n")
  print(coef_var(x, digits), print.gap = 2, quote = FALSE)
  cat("\nLog Likelihood:", round(x[["log_Lik"]], digits))
  cat("\n\n")
  invisible(x)
}
