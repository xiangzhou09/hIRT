utils::globalVariables(c("H", "J", "K", "Li", "Lijk", "Lik", "N",
                  "alpha", "dalpha", "lambda", "p", "q", "pik",
                  "theta_ls", "w", "x", "y", "z", "fitted_mean", "fitted_var"))

# median impute
impute <- function(vec){
  vec[is.na(vec)] <- median(vec, na.rm = TRUE)
  vec
}

# logical or infix function
`%||%` <- function(a, b) if (!is.null(a)) a else b

# calculate gamma gradient for case i
si_gamma <- function(i) {
    sum(pik[i, ] * Lik[i, ] * (theta_ls - fitted_mean[[i]]))/
    fitted_var[[i]]/Li[[i]] * x[i, 1:p]
}

# calculate lambda gradient for case i
si_lambda <- function(i) {
    sum(0.5 * pik[i, ] * Lik[i, ] * ((theta_ls - fitted_mean[[i]])^2/fitted_var[[i]] - 1))/Li[[i]] * z[i, 1:q]
}
