# C++ kernels against their R counterparts

test_that("cumlogit_fit matches rms::lrm.fit (ordinal and binary, with weights)", {
  set.seed(1)
  n <- 400
  X <- cbind(theta = rnorm(n), x0 = rbinom(n, 1, 0.5))
  for (H in c(2L, 5L)) {
    eta <- 1.2 * X[, 1] - 0.5 * X[, 2]
    cuts <- seq(1, -1, length.out = H - 1)
    y <- 1L + rowSums(outer(eta, cuts, "+") > qlogis(runif(n)))
    wt <- runif(n, 0.2, 2)
    ref <- rms::lrm.fit(X, y, weights = wt)$coefficients
    got <- hIRT:::cumlogit_fit(y, X, wt, H, c(seq(1, -1, length.out = H - 1), 0, 0))
    expect_equal(unname(got), unname(ref), tolerance = 1e-6)
    expect_identical(names(got), names(ref))
  }
})

test_that("gamma_glm_cpp and wls_cpp match glm.fit and lm.wfit", {
  set.seed(2)
  n <- 500
  Z <- cbind(1, rnorm(n), rbinom(n, 1, 0.4))
  y <- rgamma(n, shape = 2, rate = 2 / exp(Z %*% c(0.2, 0.3, -0.4)))
  ref <- glm.fit(Z, y, family = Gamma(link = "log"), intercept = FALSE)
  got <- hIRT:::gamma_glm_cpp(Z, y)
  expect_equal(as.double(got$coefficients), unname(ref$coefficients), tolerance = 1e-8)
  w <- runif(n)
  expect_equal(as.double(hIRT:::wls_cpp(Z, y, w)), unname(lm.wfit(Z, y, w)$coefficients),
               tolerance = 1e-10)
})

test_that("E-step weights match a direct R computation", {
  set.seed(3)
  N <- 60; J <- 4; K <- 9
  H <- c(2L, 3L, 4L, 5L)
  Y <- sapply(H, function(h) sample(seq_len(h), N, replace = TRUE))
  Y[sample(N * J, 20)] <- NA
  storage.mode(Y) <- "integer"
  alpha <- lapply(H, function(h) c(Inf, seq(1, -1, length.out = h - 1), -Inf))
  beta <- c(0.8, 1.2, -0.5, 1)
  theta <- seq(-4, 4, length.out = K); qw <- rep(8 / K, K)
  m <- rnorm(N, sd = 0.3); v <- exp(rnorm(N, sd = 0.2))
  logP <- function(j, th) {
    a <- alpha[[j]]; y <- Y[, j]
    out <- log(plogis(a[y] + beta[j] * th) - plogis(a[y + 1] + beta[j] * th))
    out[is.na(out)] <- 0
    out
  }
  post <- sapply(seq_len(K), function(k)
    exp(Reduce(`+`, lapply(seq_len(J), logP, th = theta[k]))) *
      dnorm(theta[k], m, sqrt(v)) * qw[k])
  ref_w <- t(post / rowSums(post))
  got <- hIRT:::estep_grm(Y, alpha, beta, theta, qw, m, v)
  expect_equal(got$w, ref_w, tolerance = 1e-12)
  expect_equal(got$loglik, sum(log(rowSums(post))), tolerance = 1e-10)
})
