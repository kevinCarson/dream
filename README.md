
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dream: An R Package for Dynamic Relational Event Analysis and Modeling

<!-- badges: start -->

[![](https://www.r-pkg.org/badges/version/dream?color=blue)](https://cran.r-project.org/package=dream)
[![](http://cranlogs.r-pkg.org/badges/grand-total/dream?color=red)](https://cran.r-project.org/package=dream)
[![R build
status](https://github.com/kevinCarson/dream/workflows/R-CMD-check/badge.svg)](https://github.com/kevinCarson/dream/actions)
[![CRAN
checks](https://badges.cranchecks.info/summary/dream.svg)](https://cran.r-project.org/web/checks/check_results_dream.html)

[![R-CMD-check](https://github.com/kevinCarson/dream/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/kevinCarson/dream/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The dream package provides users with helpful functions for relational
event modeling/analysis. In particular, dream provides users with helper
functions for large relational event analysis, such as recently proposed
sampling procedures for creating relational risk sets. Alongside the set
of functions for relational event analysis, this package includes
functions for the structural analysis of one- and two-mode networks,
such as network constraint and effective size measures.

This package was developed with support from the National Science
Foundation’s (NSF) Human Networks and Data Science Program (HNDS) under
award number 2241536 (PI: Diego F. Leal). Any opinions, findings, and
conclusions, or recommendations expressed in this material are those of
the authors and do not necessarily reflect the views of the NSF.

## Authors

## Installation

You can install the stable verison of `dream` from
[CRAN](https://cran.r-project.org/package=dream) via:

``` r
install.packages("dream")
```

You can install the development version of dream from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("kevinCarson/dream")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(dream)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
