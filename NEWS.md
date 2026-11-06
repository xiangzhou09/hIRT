# hIRT (development version)

* New `mhgrm()` fits a multidimensional hierarchical graded response model with
  up to three latent dimensions. Each item measures one dimension. The means
  of the latent dimensions depend on `x`, and their variances and covariances
  depend on `z` through a modified Cholesky parameterization. Standard errors
  are computed under the identification constraints. By default the E-step
  uses a fixed Gauss-Legendre grid as in `hgrm()`; `control = list(adapt =
  TRUE)` uses adaptive Gauss-Hermite quadrature (grids centered and scaled on
  each respondent's posterior), which is slower per iteration but stays
  accurate when dimensions have many items or are highly correlated. If the
  log-likelihood falls well below its maximum on a fixed grid (as when weakly
  identified dimensions push the prior distributions off the grid) or the
  level-II update fails, `mhgrm()` stops with a warning and returns the estimates with the highest
  log-likelihood; it also warns when the fixed grid does not cover the prior
  distributions. If a unidimensional fit used for starting values does not
  converge, it is retried with a looser tolerance.
* New `mhgrm_cov()` reports the implied variances, covariances, and
  correlations of the latent dimensions at given covariate values, with
  delta-method standard errors.
* New dataset `nes2008` adds 13 civil rights and 9 moral items to the
  respondents and economic items of `nes_econ2008`.

# hIRT 0.5.0

* The computational core is now written in C++ (Rcpp/RcppArmadillo): the E-step,
  the item-level M-step (a Newton solver for the weighted cumulative logit model,
  replacing `rms::lrm.fit()` and `glm.fit()` inside the EM loop), the level-II
  variance regression, and the scores used for standard errors. Estimates,
  standard errors, and log-likelihoods agree with version 0.4.0 to numerical
  precision, and fits are roughly 10 to 100 times faster.
* The E-step is computed on the log scale, so it no longer underflows when a
  respondent answers many items.
* `hgrmDIF()` pools respondents who share the same DIF covariate pattern in the
  M-step, which makes it much faster when `x0` is categorical.

# hIRT 0.4.0

* Removed dependency on the archived package pryr; function composition is now done in base R.
* Added `hgrmDIF()` for differential item functioning in the graded response model.

# hIRT 0.3.1

* fixed bug in README.md




