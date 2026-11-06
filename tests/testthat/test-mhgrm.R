# Multidimensional hierarchical graded response model

test_that("mhgrm with one dimension reproduces hgrm estimates", {
  y <- nes_econ2008[, -(1:3)]
  x <- model.matrix(~ party * educ, nes_econ2008)
  z <- model.matrix(~ party, nes_econ2008)
  ref <- hgrm(y, x, z, control = list(eps = 1e-6))
  # same fixed grid as hgrm
  fit <- mhgrm(y, dims = rep(1, ncol(y)), x, z, control = list(eps = 1e-6, adapt = FALSE))
  expect_equal(fit$coefficients$Estimate, ref$coefficients$Estimate, tolerance = 1e-4)
  expect_equal(fit$log_Lik, ref$log_Lik, tolerance = 1e-6)
  expect_equal(fit$scores$post_mean[, 1], ref$scores$post_mean, tolerance = 1e-4)
  expect_true(all(is.finite(fit$coefficients$Std_Error)))
  # adaptive quadrature matches a wide, fine fixed grid; hgrm's default range
  # [-4, 4] is too narrow here, so it is not the reference
  fit_ad <- mhgrm(y, dims = rep(1, ncol(y)), x, z, control = list(eps = 1e-6, adapt = TRUE))
  fine <- mhgrm(y, dims = rep(1, ncol(y)), x, z,
                control = list(eps = 1e-6, adapt = FALSE, K = 41, C = 7))
  expect_equal(fit_ad$coefficients$Estimate, fine$coefficients$Estimate, tolerance = 1e-4)
  expect_equal(fit_ad$log_Lik, fine$log_Lik, tolerance = 1e-6)
})

sim_md <- function(N, rho0, rho1, seed = 1, n_items = 4) {
  set.seed(seed)
  z1 <- rbinom(N, 1, 0.5)
  theta <- matrix(rnorm(2 * N), N)
  r <- ifelse(z1 == 1, rho1, rho0)
  theta[, 2] <- r * theta[, 1] + sqrt(1 - r^2) * theta[, 2]
  y <- sapply(rep(1:2, each = n_items), function(k)
    findInterval(theta[, k] + rlogis(N), c(-1.5, -0.5, 0.5, 1.5)))
  list(y = y, x = cbind(1, z1), z = cbind(1, z1))
}

test_that("mhgrm recovers covariate-dependent correlations in two dimensions", {
  s <- sim_md(1500, rho0 = 0.1, rho1 = 0.7)
  fit <- mhgrm(s$y, dims = rep(c("a", "b"), each = 4), x = s$x, z = s$z,
               control = list(K = 15, adapt = TRUE))
  expect_s3_class(fit, "mhgrm")
  expect_equal(dim(fit$scores$post_mean), c(1500L, 2L))
  cc <- mhgrm_cov(fit, z = rbind(c(1, 0), c(1, 1)))
  cors <- cc$estimate[cc$quantity == "cor[a, b]"]
  expect_length(cors, 2)
  expect_lt(abs(cors[1] - 0.1), 0.15)
  expect_lt(abs(cors[2] - 0.7), 0.15)
  expect_true(all(cc$std_error > 0))
  # location constraint: prior means average to zero on each dimension
  expect_equal(unname(colMeans(fit$scores$prior_mean)), c(0, 0), tolerance = 1e-8)
})

test_that("adaptive quadrature matches a fine fixed grid", {
  s <- sim_md(1000, rho0 = 0.3, rho1 = 0.3)
  dims <- rep(c("a", "b"), each = 4)
  ad <- mhgrm(s$y, dims = dims, control = list(eps = 1e-6, adapt = TRUE))
  fine <- mhgrm(s$y, dims = dims, control = list(eps = 1e-6, adapt = FALSE, K = 41, C = 6))
  expect_equal(ad$log_Lik, fine$log_Lik, tolerance = 1e-6)
  expect_equal(ad$coefficients$Estimate, fine$coefficients$Estimate, tolerance = 1e-3)
})

test_that("long, highly correlated dimensions are estimated without degenerating", {
  # with many items per dimension the posterior is much narrower than a fixed
  # grid's spacing; adaptive quadrature keeps the correlation away from 1
  s <- sim_md(1000, rho0 = 0.9, rho1 = 0.9, n_items = 25)
  fit <- mhgrm(s$y, dims = rep(c("a", "b"), each = 25), control = list(adapt = TRUE))
  cc <- mhgrm_cov(fit, z = matrix(1, 1, 1))
  rho <- cc$estimate[cc$quantity == "cor[a, b]"]
  expect_lt(abs(rho - 0.9), 0.03)
  expect_true(all(is.finite(fit$coefficients$Std_Error)))
})

test_that("the fixed grid is the default", {
  s <- sim_md(300, rho0 = 0.3, rho1 = 0.3)
  dims <- rep(c("a", "b"), each = 4)
  def <- mhgrm(s$y, dims = dims)
  fixed <- mhgrm(s$y, dims = dims, control = list(adapt = FALSE))
  expect_equal(def$log_Lik, fixed$log_Lik)
  expect_equal(def$coefficients$Estimate, fixed$coefficients$Estimate)
})

test_that("ordinary fits raise no warnings", {
  s <- sim_md(800, rho0 = 0.3, rho1 = 0.6)
  expect_no_warning(mhgrm(s$y, dims = rep(c("a", "b"), each = 4), x = s$x, z = s$z))
  y <- nes_econ2008[, -(1:3)]
  expect_no_warning(mhgrm(y, dims = rep(1, ncol(y)), model.matrix(~ party * educ, nes_econ2008),
                          model.matrix(~ party, nes_econ2008)))
})

test_that("a fixed grid too narrow for the prior distributions is flagged", {
  set.seed(3)
  N <- 1500
  g <- rbinom(N, 1, 0.5)
  theta <- 3.5 * g + rnorm(N)                       # group means far apart
  y <- sapply(1:4, function(j) findInterval(theta + rlogis(N), c(-1, 1, 3, 5)))
  x <- cbind(1, g)
  expect_warning(mhgrm(y, dims = rep(1, 4), x = x), "outside its range")
  expect_no_warning(mhgrm(y, dims = rep(1, 4), x = x, control = list(C = 8, K = 41)))
})

test_that("mhgrm rejects unsupported inputs", {
  s <- sim_md(200, 0.2, 0.2)
  expect_error(mhgrm(s$y, dims = rep(1:4, each = 2)), "At most three")
  expect_error(mhgrm(s$y, dims = c(1, rep(2, 7))), "at least two items")
  expect_error(mhgrm(s$y, dims = rep(1:2, 3)), "one entry per column")
})
