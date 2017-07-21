utils::globalVariables(c("H", "J", "K", "Li", "Lijk", "Lik", "N",
                  "alpha", "dalpha", "lambda", "p", "q", "pik",
                  "theta_ls", "w", "x", "y", "z"))

# convert original response data to integers
fac2int <- function(x) as.integer(factor(x, exclude = c(NA, NaN)))

# calculate gamma gradient for case i
si_gamma <- function(i) {
    sum(pik[i, ] * Lik[i, ] * (theta_ls - as.vector(x[i, ] %*% gamma)))/exp(as.vector(z[i, ] %*%
        lambda))/Li[i] * x[i, 2:p]
}

# calculate lambda gradient for case i
si_lambda <- function(i) {
    sum(0.5 * pik[i, ] * Lik[i, ] * ((theta_ls - as.vector(x[i, ] %*% gamma))^2/exp(as.vector(z[i,
        ] %*% lambda)) - 1))/Li[i] * z[i, 2:q]
}
