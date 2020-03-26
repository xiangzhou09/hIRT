
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hIRT: hierarchical item response theory (IRT) models

hIRT is an R package that implements a class of hierarchical item
response theory (IRT) models where both the mean and the variance of the
latent “ability parameters” may depend on observed covariates. The
current implementation includes both the two-parameter latent trait
model for binary data (`hltm()` and `hltm2()`) and the graded response
model for ordinal data (`hgrm()` and `hgrm2()`). Both are fitted via the
Expectation-Maximization (EM) algorithm. Asymptotic standard errors are
derived from the observed information matrix.

**Main Reference**: Zhou, Xiang. 2019. “Hierarchical Item Response
Models for Analyzing Public Opinion.” Political Analysis, 27(4):
481-502. Available at: <https://doi.org/10.1017/pan.2018.63>

Full paper with technical appendix is available at:
<https://scholar.harvard.edu/files/xzhou/files/Zhou2019_hIRT.pdf>

## Installation

You can install the released version of hIRT from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("hIRT")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("xiangzhou09/hIRT")
```

## Example

The following example illustrates how the `hgrm()` function can be used
to examine the effects of education and party affiliation on economic
ideology, a latent variable gauged by a number of survey items in the
American National Election Studies (ANES), 2008. Documentation of the
dataset `nes_econ2008` can be accessed by running `?nes_econ2008` in R
after loading the `hIRT` package.

``` r
library(hIRT)
#> Registered S3 method overwritten by 'pryr':
#>   method      from
#>   print.bytes Rcpp

# survey items used to measure economic ideology
y <- nes_econ2008[, -(1:3)]

# predictors for the mean of economic ideology
x <- model.matrix( ~ party * educ, nes_econ2008)

# predictors for the variance of economic ideology
z <- model.matrix( ~ party, nes_econ2008)

# fitting a hierarhical graded response model
nes_m1 <- hgrm(y, x, z)
#> ............
#>  converged at iteration 12

nes_m1
#> 
#> Call:
#> hgrm(y = y, x = x, z = z)
#> 
#> Mean Regression:
#>                          Estimate Std_Error z_value p_value
#> x_(Intercept)              -0.480     0.105  -4.570   0.000
#> x_partyindependent          0.386     0.086   4.473   0.000
#> x_partyRepublican           1.133     0.135   8.408   0.000
#> x_educ2                     0.037     0.079   0.467   0.641
#> x_partyindependent:educ2    0.235     0.117   2.007   0.045
#> x_partyRepublican:educ2     0.428     0.148   2.886   0.004
#> 
#> Variance Regression:
#>                    Estimate Std_Error z_value p_value
#> z_(Intercept)        -0.097     0.139  -0.697   0.486
#> z_partyindependent    0.166     0.100   1.661   0.097
#> z_partyRepublican     0.172     0.126   1.373   0.170
#> 
#> Log Likelihood: -16259.16
```

The output from `hgrm` is an object of class `hIRT`. The `print()`
method for `hIRT` outputs the regression tables for the mean regression
and the variance regression.

## Extracting coefficients

The `coef_item()`, `coef_mean()`, and `coef_var()` functions can be used
to extract coefficient tables for item parameters, the mean regression,
and the variance regression respectively.

``` r
coef_item(nes_m1)
#> $health_ins7
#>        Estimate Std_Error z_value p_value
#> y>=2      1.279        NA      NA      NA
#> y>=3      0.541     0.063   8.542   0.000
#> y>=4     -0.075     0.083  -0.898   0.369
#> y>=5     -1.047     0.107  -9.826   0.000
#> y>=6     -1.852     0.124 -14.901   0.000
#> y>=7     -2.684     0.149 -17.990   0.000
#> Dscrmn    1.016     0.096  10.569   0.000
#> 
#> $jobs_guar7
#>        Estimate Std_Error z_value p_value
#> y>=2      2.136     0.173  12.377       0
#> y>=3      1.352     0.153   8.860       0
#> y>=4      0.607     0.141   4.299       0
#> y>=5     -0.520     0.137  -3.797       0
#> y>=6     -1.611     0.141 -11.429       0
#> y>=7     -2.785     0.163 -17.043       0
#> Dscrmn    1.305     0.114  11.448       0
#> 
#> $gov_services7
#>        Estimate Std_Error z_value p_value
#> y>=2      3.950     0.222  17.760   0.000
#> y>=3      2.859     0.182  15.707   0.000
#> y>=4      1.831     0.158  11.592   0.000
#> y>=5      0.247     0.147   1.679   0.093
#> y>=6     -1.001     0.154  -6.490   0.000
#> y>=7     -2.020     0.169 -11.947   0.000
#> Dscrmn   -1.363     0.116 -11.715   0.000
#> 
#> $FS_poor3
#>        Estimate Std_Error z_value p_value
#> y>=2     -1.180     0.179  -6.601       0
#> y>=3     -4.459     0.243 -18.357       0
#> Dscrmn    1.918     0.164  11.679       0
#> 
#> $FS_childcare3
#>        Estimate Std_Error z_value p_value
#> y>=2     -0.808     0.148  -5.474       0
#> y>=3     -4.051     0.192 -21.132       0
#> Dscrmn    1.608     0.128  12.535       0
#> 
#> $FS_crime3
#>        Estimate Std_Error z_value p_value
#> y>=2     -0.845     0.066 -12.866       0
#> y>=3     -3.150     0.108 -29.048       0
#> Dscrmn    0.516     0.059   8.823       0
#> 
#> $FS_publicschools3
#>        Estimate Std_Error z_value p_value
#> y>=2     -1.790     0.136 -13.197       0
#> y>=3     -4.144     0.188 -22.022       0
#> Dscrmn    1.302     0.111  11.751       0
#> 
#> $FS_welfare3
#>        Estimate Std_Error z_value p_value
#> y>=2      1.054     0.117   8.970       0
#> y>=3     -1.355     0.116 -11.650       0
#> Dscrmn    1.178     0.099  11.937       0
#> 
#> $FS_envir3
#>        Estimate Std_Error z_value p_value
#> y>=2     -0.855     0.106  -8.071       0
#> y>=3     -3.499     0.159 -22.023       0
#> Dscrmn    1.101     0.092  11.953       0
#> 
#> $FS_socsec3
#>        Estimate Std_Error z_value p_value
#> y>=2     -1.091     0.104 -10.535       0
#> y>=3     -4.278     0.178 -24.033       0
#> Dscrmn    1.028        NA      NA      NA

coef_mean(nes_m1)
#>                          Estimate Std_Error z_value p_value
#> x_(Intercept)              -0.480     0.105  -4.570   0.000
#> x_partyindependent          0.386     0.086   4.473   0.000
#> x_partyRepublican           1.133     0.135   8.408   0.000
#> x_educ2                     0.037     0.079   0.467   0.641
#> x_partyindependent:educ2    0.235     0.117   2.007   0.045
#> x_partyRepublican:educ2     0.428     0.148   2.886   0.004

coef_var(nes_m1)
#>                    Estimate Std_Error z_value p_value
#> z_(Intercept)        -0.097     0.139  -0.697   0.486
#> z_partyindependent    0.166     0.100   1.661   0.097
#> z_partyRepublican     0.172     0.126   1.373   0.170
```

## Latent scores

The `latent_scores()` function can be used to extract the Expected A
Posteriori (EAP) estimates of the latent ability parameters, along with
their “prior” estimates (without the random effects). In this example,
the latent ability estimates can be interpreted as the estimated
ideological positions of ANES respondents on economic issues.

``` r

pref <- latent_scores(nes_m1)

summary(pref)
#>    post_mean            post_sd         prior_mean            prior_sd    
#>  Min.   :-2.082000   Min.   :0.3940   Min.   :-0.4800000   Min.   :0.953  
#>  1st Qu.:-0.751000   1st Qu.:0.4788   1st Qu.:-0.4440000   1st Qu.:0.953  
#>  Median :-0.104000   Median :0.5280   Median :-0.0950000   Median :1.035  
#>  Mean   :-0.000147   Mean   :0.5469   Mean   :-0.0001561   Mean   :1.001  
#>  3rd Qu.: 0.629500   3rd Qu.:0.6090   3rd Qu.: 0.1770000   3rd Qu.:1.035  
#>  Max.   : 3.359000   Max.   :0.9780   Max.   : 1.1170000   Max.   :1.039
```

## Identification constraints.

The `constr` parameter in the `hgrm()` and `hltm()` function can be used
to specify the type of constraints used to identify the model. The
default option, `"latent_scale"`, constrains the mean of the latent
ability parameters to zero and the geometric mean of their prior
variance to one; Alternatively, `"items"` sets the mean of the item
difficulty parameters to zero and the geometric mean of the
discrimination parameters to one.

In practice, one may want to interpret the effects of the mean
predictors (in the above example, education and party affiliation) on
the standard deviation scale of the latent trait. This can be easily
achieved through rescaling their point estimates and standard errors.

``` r

library(dplyr)

total_sd <- sqrt(var(pref$post_mean) + mean(pref$post_sd^2))

coef_mean_sd_scale <- coef_mean(nes_m1) %>%
  mutate(`Estimate` = `Estimate`/total_sd,
         `Std_Error` = `Std_Error`/total_sd)

coef_mean_sd_scale
#>      Estimate  Std_Error z_value p_value
#> 1 -0.42437486 0.09283200  -4.570   0.000
#> 2  0.34126812 0.07603383   4.473   0.000
#> 3  1.00170150 0.11935543   8.408   0.000
#> 4  0.03271223 0.06984503   0.467   0.641
#> 5  0.20776686 0.10344137   2.007   0.045
#> 6  0.37840092 0.13084892   2.886   0.004
```

## hIRT with fixed item parameters

Sometimes, the researcher might want to fit the hIRT models using a set
of fixed item parameters, for example, to make results comparable across
different studies. The `hgrm2()` and `hltm2()` functions can be used for
this purpose. They are illustrated in more detail in the package
documentation.
