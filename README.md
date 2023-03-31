
<!-- README.md is generated from README.Rmd. Please edit that file -->


# fda.usc: Functional Data Analysis and Utilities for Statistical Computing 

<!-- ![](inst/figures/fda.usc.png)
pkgdown::build_site()
-->

[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/fda.usc)](https://cran.r-project.org/package=fda.usc)
[![Licence](https://img.shields.io/badge/licence-GPL--2-blue.svg)](https://www.gnu.org/licenses/gpl-2.0.en.html)
[![](https://cranlogs.r-pkg.org/badges/fda.usc)](https://cran.r-project.org/package=fda.usc)

## Package overview

**fda.usc** package carries out exploratory and descriptive analysis of
functional data exploring its most important features such as depth
measurements or functional outliers detection, among others. It also
helps to explain and model the relationship between a dependent variable
and independent (regression models) and make predictions. Methods for
supervised or unsupervised classification of a set of functional data
regarding a feature of the data are also included. It can perform
functional ANOVA, hypothesis testing, functional response models and
many others.

## Installation

You can install the current fda.usc version from CRAN with:

``` r
install.packages("fda.usc")
```

or the latest patched version from Github with:

``` r
library(devtools)
devtools::install_github("moviedo5/fda.usc")
```

## Issues & Feature Requests

For issues, bugs, feature requests etc. please use the [Github
Issues](https://github.com/moviedo5/fda.usc/issues). Input is always
welcome.

## Documentation

A hands on introduction to  can be found in the reference
[vignette](https://www.jstatsoft.org/article/view/v051i04/).

Details on specific functions are in the [reference
manual](https://cran.r-project.org/package=fda.usc/fda.usc.pdf).

Cheatsheet [fda.usc reference
card](https://zenodo.org/record/3386752/files/RefCard_fda.usc_v1.pdf?download=1).

## References

Febrero-Bande, M. and Oviedo de la Fuente, M. (2012). Statistical
Computing in Functional Data Analysis: The R Package fda.usc. *Journal
of Statistical Software*, 51(4):1-28.
<https://dx.doi.org/10.18637/jss.v051.i04>

<!-- 
<https://www.jstatsoft.org/v51/i04/>
library(roxygen2)
setwd("C:/Users/moviedo/github/fda.usc/")
getwd()
pkgbuild::compile_dll()
roxygenize()
devtools::document()

library(devtools)

# devtools::install_github("moviedo5/fda.usc",auth_user="moviedo5")

R CMD build fda.usc
R CMD check fda.usc_2.1.0.tar.gz --as-cran
R CMD INSTALL fda.usc_2.1.0.tar.gz --build

Manuel Oviedo PhD thesis [Advances in functional regression and classification models](https://hdl.handle.net/10347/18236)

-->


