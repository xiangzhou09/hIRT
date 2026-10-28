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




