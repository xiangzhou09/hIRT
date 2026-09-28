## Resubmission

This is a resubmission of hIRT, which was archived on 2026-01-30 because it
required the archived package 'pryr'. The pryr dependency has been removed and
replaced with base R code; results are numerically identical to version 0.3.0's
implementation.

In response to the reviewer's comments on the previous submission:

* Added a reference describing the methods to the Description field of the
  DESCRIPTION file: Zhou (2019) <doi:10.1017/pan.2018.63>.
* Added a \value section to print.hIRT.Rd, documenting that the function
  invisibly returns its input and is called for the side effect of printing.
* Removed unsuppressible console output from hgrm(), hgrm2(), hgrmDIF(),
  hltm(), and hltm2(). EM progress messages are now printed only when the
  user sets `control = list(verbose = TRUE)`; the default is `verbose = FALSE`.

We also corrected several documented default values of control parameters
(eps, K, C) so that they match the code.

## Test environments
* local Windows 11, R 4.5.1
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

* New submission (package was archived on CRAN).
* Possibly misspelled words in DESCRIPTION: Zhou. This is the author's surname
  in the added method reference, Zhou (2019), and is spelled correctly.

## Downstream dependencies

There are currently no downstream dependencies for this package.
