# End-to-end agreement with hIRT 0.4.0 (pure R implementation). The reference
# values were produced by the CRAN/GitHub version on the same calls.

ref <- readRDS(test_path("ref_hIRT_0.4.0.rds"))
y <- nes_econ2008[, -(1:3)]
x <- model.matrix(~ party * educ, nes_econ2008)
z <- model.matrix(~ party, nes_econ2008)
ybin <- y
ybin[] <- lapply(y, function(v) as.integer(v > median(v, na.rm = TRUE)))

check_fit <- function(fit, r, tol) {
  expect_identical(rownames(fit$coefficients), rownames(r$coefficients))
  expect_equal(fit$coefficients$Estimate, r$coefficients$Estimate, tolerance = tol)
  expect_equal(fit$coefficients$Std_Error, r$coefficients$Std_Error, tolerance = tol)
  expect_equal(fit$log_Lik, r$log_Lik, tolerance = tol)
  expect_equal(fit$scores$post_mean, r$scores$post_mean, tolerance = tol)
}

test_that("hgrm reproduces hIRT 0.4.0", {
  check_fit(hgrm(y, x, z), ref$hgrm_nes, 1e-6)
  check_fit(hgrm(y, x, z, constr = "items"), ref$hgrm_nes_items, 1e-6)
})

test_that("hltm reproduces hIRT 0.4.0", {
  check_fit(hltm(ybin, x, z), ref$hltm_nes, 1e-5)
})

test_that("hgrm2 reproduces hIRT 0.4.0", {
  item_coefs <- lapply(coef_item(hgrm(y, x, z)), `[[`, "Estimate")
  check_fit(hgrm2(y, x, z, item_coefs = item_coefs), ref$hgrm2_nes, 1e-6)
})

test_that("hgrmDIF reproduces hIRT 0.4.0", {
  check_fit(hgrmDIF(y, x, items_dif = 1:2), ref$hgrmDIF_unif, 1e-5)
  check_fit(hgrmDIF(y, x, items_dif = 1, form_dif = "non-uniform"), ref$hgrmDIF_nonunif, 1e-5)
})
