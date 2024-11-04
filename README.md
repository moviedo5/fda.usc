
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
manual](docs/fda.usc-manual.pdf).

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
# setwd("D:/Users/moviedo/github/fda.usc/")
getwd()
pkgbuild::compile_dll()
roxygenize()
devtools::document()

tools::checkRd("man/Outliers.fdata.Rd")
tools::checkRd("man/LMDC.select.Rd")
tools::checkRd("man/accuracy.Rd")
tools::checkRd("man/classif.DD.Rd")
tools::checkRd("man/classif.depth.Rd")
tools::checkRd("man/cond.mode.Rd")
tools::checkRd("man/depth.mdata.Rd")
tools::checkRd("man/dfv.test.Rd")
tools::checkRd("man/fEqDistrib.test.Rd")
tools::checkRd("man/fanova.RPm.Rd")
tools::checkRd("man/fanova.hetero.Rd")
tools::checkRd("man/fanova.onefactor.Rd")
tools::checkRd("man/fdata.bootstrap.Rd")
tools::checkRd("man/fdata2basis.Rd")
tools::checkRd("man/fdata2pc.Rd")
tools::checkRd("man/flm.Ftest.Rd")
tools::checkRd("man/flm.test.Rd")
tools::checkRd("man/fregre.basis.Rd")
tools::checkRd("man/fregre.basis.cv.Rd")
tools::checkRd("man/fregre.basis.fr.Rd")
tools::checkRd("man/fregre.bootstrap.Rd")
tools::checkRd("man/fregre.gkam.Rd")
tools::checkRd("man/fregre.glm.Rd")
tools::checkRd("man/fregre.gls.Rd")
tools::checkRd("man/fregre.gsam.Rd")
tools::checkRd("man/fregre.igls.Rd")
tools::checkRd("man/fregre.lm.Rd")
tools::checkRd("man/fregre.np.Rd")
tools::checkRd("man/fregre.np.cv.Rd")
tools::checkRd("man/fregre.pc.Rd")
tools::checkRd("man/fregre.pc.cv.Rd")
tools::checkRd("man/fregre.pls.Rd")
tools::checkRd("man/fregre.pls.cv.Rd")
tools::checkRd("man/fregre.plm.Rd")
tools::checkRd("man/influence.fregre.fd.Rd")
tools::checkRd("man/influence_quan.Rd")
tools::checkRd("man/optim.basis.Rd")
tools::checkRd("man/optim.np.Rd")
tools::checkRd("man/predict.classif.DD.Rd")
tools::checkRd("man/predict.fregre.lm.Rd")
tools::checkRd("man/summary.fregre.fd.Rd")
tools::checkRd("man/rp.flm.test.Rd")
tools::checkRd("man/GCCV.S.Rd")

roxygenize()
devtools::document()




library(devtools)
devtools::build()
devtools::check()
devtools::build_win()

# devtools::install_github("moviedo5/fda.usc",auth_user="moviedo5")
R CMD check --as-cran and R-wind-builder 
 
R CMD build fda.usc
R CMD check fda.usc_2.2.0.tar.gz --as-cran  R-wind-builder 
R CMD INSTALL fda.usc_2.2.0.tar.gz --build

Manuel Oviedo PhD thesis [Advances in functional regression and classification models](https://hdl.handle.net/10347/18236)

-->


