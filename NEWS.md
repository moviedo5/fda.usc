# fda.usc 2.2.0 

* Modification in several Rd fiels: Lost braces in itemize.
* classif.DD.aux.r, `solve.ab`  has been renamed to `solve_ab`.
* classif.kfold.R, `all.vars1`  has been renamed to `all_vars1`.
* kmeans.fd.dist.R, `predict.kmeans.fd`  has been renamed to `predict_kmeans.fd`.
* kmeans.fd.R, `kmeans.assig.groups`  has been renamed to `kmeans_assig_groups`.
* kmeans.fd.R, `kmeans.centers.update`  has been renamed to `kmeans_centers_update`.
* kmeans.assig.groups.R, `kmeans.assig.groups`  has been renamed to `kmeans_assig_groups`.
* kmeans.assig.groups.R, `kmeans.assig.groups`  has been renamed to `kmeans_assig_groups`.
* plot.fdata.R, `image.scale`  has been renamed to `image_scale`.
* predict.gls.nlme.r, `predict.gls`  has been renamed to `predict_gls`.
* predict.gls.nlme.r, `predict.mregre`  has been renamed to `predict_._mregre`.
* predict.mfregre.r, `predict.gls`  has been renamed to `predict_gls`.
* quantile.outliers.pond.r, `quantile_outliers.pond`  has been renamed to `quantile_outliers_pond`.
* quantile.outliers.trim.r, `quantile_outliers.trim`  has been renamed to `quantile_outliers_trim`.
* the `plot.lfdata` function is exported, no longer internal


# fda.usc 2.1.0 

* fda.usc 2.1.0 is a major release with several new feature and fixed bugs.

* `fdata2basis()` always return centred fdataobj and mean. The mean is computed 
  using the basis.

* New funtions: `fEqMoments.test()`, `fmean.test.fdata()`, `cov.test.fdata()` for checking 
  the equality of means and/or covariance between two populations under gaussianity. 

* New funtions: `fEqDistrib.test()`, `XYRP.test()`, `MMD.test()`, `MMDA.test()`, 
  `fEqDistrib.test()` for checking the equality of distributions between two functional populations. 
  
* The wavelength units in the Tecator dataset are labeled as nm in dataset description, but they are nm in the original description. (bug detected by vnmabus)

  

# fda.usc 2.0.3

* Several changes on `summary.fdata.comp()`
* Now, we uses `data.matrix` (instead of as.matrix) to convert data.table in a matrix class object, option recommended when data.frame contains characters.
* `classif.cv.glmnet()` and `classif.gbm()`, functional basis classsification using `cv.glmnet()`, require glmnet package, and `gbm()`, require gbm package.
* `h.default()`,  new argument 'Ker'
* `mfdata()`, new class object for multivariate functional data
* `fregre.basis.cv()`, `fregre.basis.cv()` and `fregre.pc()` return df.residual object
* Bug corrected in `S.LPR()` and `S.LLR()`
* `classif.gsam.vs`, new function for variable selection in additive classifier
* `fdata2basis()` is used in `fregre.lm()` and `predict.fregre.lm()`
* summary for `fdata2basis()`
 
# fda.usc 2.0.2

* `fdata.bootstrap()` and `fregre.bootstrap()` functions addapted to parallel backend.
* Corrected bug in `kmeans.fd()` function.
* Bug corrected in `predict.fregre.glm()`, `predict.fregre.lm()` and `predict.gsam()`, now works with type="effects".
* The "main", "xlab" and "ylab" arguments can be used in the `plot.bifd()` function.
* `plot.ldata()` draws each curve according to the factor indicated in the argument "var.name".
* Minor changes in `ldata()`, "na.rm = T" is removed in the sweep function.
* The documentation for the argument '...' that is not included in the "usage" was deleted.
* Modification in internal function `wmestadis()` used in `fanova.onefactor()`, now it replicates the statistic defined by Cuevas, 2004.
* Deleted arguments "y" and "corplot" in `summary.fdata.comp()` function.
* Deleted `dev.new()` in code of `summary.fdata.comp()` function.

# fda.usc 2.0.1

* Modification in `fdata()` function to avoid class()== and class()!= instead use `is()`,
* `kmeans.fd()` function:

    a. "par.ini" argument is depreciated, the user can use "method" argument. 
    b. "cluster.size" argument are added.
    c. New internal function `predict.kmeans.fd`.

* Bug corrected in internal function `pred2glm2boost()`, it is used for predictions of `classiff.DD()` outputs
* New functions: `Ops.ldata()`, `Math.ldata()`, `Summary.ldata()`, `mean.ldata()` and `mean.fdata()` (deprecated `ldata.mean()`, `mfdata.mean()`)

* Modifications in `ldata.cen()`

* A bug in `S.LPR()` has been fixed.

* A bug in internal function `wmestadis()` used in `fanova.onefactor()` has been fixed.


# fda.usc 2.0.0

Version 2.0.0 is a major release with several new features, including:

* `inprod.fdata()` and `metric.lp` funcitons addapted to parallel backend.

* zzz.R file includes `.onAttach()` function (welcome package message)

* ops.fda.usc.R file includes `ops.fda.usc()` function that control general parameters of packages such as ncores argument.

* par.fda.usc.R file is deleted: par.fda.usc is now an internal object created and modified by `ops.fda.usc()`

* New function `metric.DTW()` computes distances between functional data using dynamic time warping (DTW)
`metric.WDTW()` and `metric.TWED()` are extended version (not parallelized yet, pending to completed the Rd document)

* New functions `S.LPR()` and `S.LCR()` for computing smoothing matrix S by nonparametric method.

* The functions anova.hetero(), anova.onefactor(), anova.RPm(), influence.fdata(), influence.quan(), min.basis(), min.np() and unlist.fdata() are renamed fanova.hetero(), fanova.onefactor(), fanova.RPm(), `influence.fregre.fd()`, influence_quan(), `optim.basis()`, `optim.np()` and `unlist_fdata()`

* optim.np() (deprecated min.np())  allows Local polynomial regression with correlated errors using the new parameter (correl=TRUE) 

* `Kernel.correlated()` new functions

* New class: "ldata":

1. `ldata()` class definition.

2.  Redefined `metric.ldata()`, it computes distance for ldata object: list with m functional data `mfdata()` and univariate data included in a data frame called "df"

3.  New function `metric.mfdata()`: compute distance for mfdata class object:  list with m functional data

4. `plot.ldata()`: plots for ldata object, it allows drawing using a color bar.
5. `plot.mfdata()`: plot formfdata object (internal function, pending to completed the Rd document)
6. `depth.modep()`, depth.mode()  call `metric.lp()` and `metric.ldata()` propperly
7. New functions: `subset.ldata()`, `is.lfdata()`, `[.lfdata()`, `[.ldata`, `is.ldata()`,
`names.ldata()` and `c.ldata()`

* `classic.tree()` is replaced by the `classic.rpart()` (which requires the rpart library to be installed). 
The internal function `classif.tree2boost()` and the dependency of the rpart package are also removed

* New functions and utilities in accuracy.r file.  

* New functions related with Machine Learning procedures (rdepend on packages not included in "fda.usc"):
  + `classif.svm()` and `classif.naiveBayes()` (e1071 pkg), `classif.ksvm()` (personalized pkg), 
  + `classif.rpart()` (rpart pkg), `classif.nnet()` (nnet pkg), `classif.multinom()` (nnet pkg),
  + `classif.randomForest()` (randomForest)
  + `clasiff.univariante()` is used in classif.DD and allow multiclass labels
  + `classif.kfold()` selects the parameters using k-fold cross-validation 

* Minor changes in `classif.gkam()` and `fregre.gkam()`

* Settings in `fregre.np()`, `fregre.np.cv()` with type.S = S.KNN

* New script file: FDA_REviewClasif_V2 classification example

* Bug corrected in `h.default()` (specially using k-nearest neighbors smoothing, 
type.S="S.KNN")

* "type.CV" and "par.CV()" arguments are removed in `classif.np()`, 
`classif.kernel()` and `classif.knn()`

* dcor.xy.r/.Rd includes Rdnames

* fdata2model.R shortcut to use in classif and fregre method 


# fda.usc 1.5.0

* This version was released in Jan. 2019 to accompany Manuel Oviedo de la Fuente PhD Thesis, 
see Minerva (University of Santiago de Compostela) repository.
 
* New function implemented: fregre.gsam.vs() accompany paper: 
Febrero-Bande, M., Gonz\'{a}lez-Manteiga, W. and Oviedo de la Fuente, M. 
Variable selection in functional additive regression models, (2018).  
Computational Statistics, 1-19. DOI: 10.1007/s00180-018-0844-5

* The current function `fregre.basis.cv()` returns an object called fregre.basis 
(same output as if the `fregre.basis()` function had been used) that uses the selected parameters
according to the indicated criteria (see example below).  The previous function version
(up to version 1.5.0) has been renamed in the function "fregre.basis.cv.old". 
It is marked as deprecated in the current version and will be deleted in the next version of the package, thanks to Beatriz Bueno.

* New functions: plot.fregre.lm and summary.fregre.lm() solve errors in the summary of the 
in `fregre.lm()` function, thanks to Prof. Andros Kourtellos.
* A bug in `fregre.pc()` has been fixed (thanks to Prof. Eduardo Garcia-Portugues).


# fda.usc 1.4.0

* This was published in December 2017 to accompany the document: 

 + Ordonez, C., Oviedo de la Fuente, M., Roca-Pardinas, J.,  Rodriguez-Perez, J. R. (2017).
 Determining optimum wavelengths for leaf water content estimation from reflectance: A 
 distance correlation approach. \emph{Chemometrics and Intelligent Laboratory Systems}, (2018)  173,41-50 DOI: 10.1016/j.chemolab.2017.12.001.
 + New functions implemented: `LMDC.select()` and `LMDC.regression()`.

* Oviedo de la Fuente M, Febrero-Bande M, Muñoz MP, Domínguez À (2018) Predicting seasonal influenza transmission using functional regression models with temporal dependence. PLoS ONE 13(4): e0194250. DOI: 10.1371/journal.pone.0194250
  + The functions fregre.gls and fregre.igls implement the functional linear model
with dependent errors (functional GLs and functional iterative GLS respectively).
  + The "fda.usc/inst/script/" folder contains the code that reproduces the results
of the paper (see the scripts whose prefix is "PLOS").

* This package version also companion for the paper:

+ "Goodness-of-fit tests for the functional linear model
based on randomly projected empirical processes" 
Cuesta-Albertos et al., 2017). The package implements 
goodness-of-fit tests for the functional linear model with scalar response.

* A bug in functional derivative by raw derivation 
(function `fdata.deriv()` with method="diff") has been fixed, thanks to Marcos Matabuena.
* A bug in `classif.knn()` and  predict.classif() has been fixed, thanks to Ricardo Recarey.
* A bug in `CV.S()` function has been fixed, thanks to Miquel Carbajo.
* A bug in `anova.hetero()` has been fixed, thanks to Beatriz Bueno.


# fda.usc 1.3.0

* Beta version functions to Fit Functional Linear Model Using Generalized Least Squares:
`fregre.gls()`, `fregre.igls()`, `GCCV.S()`, `predict.fregre.gls()` and `predict.fregre.igls()`. Internal function "auxiliar", "corSigma()", "corStruct()".

* The functionality of the functions "+.fdata()", "-.fdata()", "*.fdata()"
an "/.fdata" has been improved.

* S3 functions for fdata class calculations: `is.na.fdata()` and `anyNA.fdata`.
Function "count.na.fdata()"  returns a vector with the number of "NA" of each curve.

* Internal function "count.na" is deprecated.

* fdata function converts "xtab" and "ftable" class object into "fdata" class object.


# fda.usc 1.2.3

* Modification of `classif.DD()` function for DDk classifier.
* New functions `length.fdata()`, `NROW.fdata()`, `NCOL.fdata()`, `gridfdata()` and `rcombfdata()`.
`depth.KFSD()`  function implements a depth measure based on Kernelized Functional Spatial Depth.   
`depth.FSD()`  function implements a depth measure based on Functional Spatial Depth.   
* A bug in `fregre.pc()` function has been fixed.
* A bug in `depth.RPD()` function has been fixed.


# fda.usc 1.2.2

* Warning message in `fregre.basis.cv()`, `fregre.pc.cv()` and `fregre.pls.cv()`, 
`fregre.basis()`, `fregre.pc()` and `fregre.pls()` functions when system is computationally singular.
* A bug in predict.classif() function using a fitted object by classif.knn() has been fixed.
* The argument "trim" is modified from 0.1 to 0.25 in the function quantile.outliers.trim
* The internal distance functions: euclidean, manhattan, minkowski and maximum allow a vector of weights.
* In `classif.DD()` function, the  polynomial classifier ("DD1", "DD2" and "DD3") uses the original procedure
 proposed by Li et al. (2012), rotating the DD-plot (to exchange abscise and ordinate). The procedure 
 extend to multi-class problems by incorporating the method of majority voting in the case of polynomial
 classifier and the method One vs the Rest in the logistic case ("glm" and "gam").
* The `fregre.gkam()` function only considers  functional covariates (not implemented for non-functional covariates). 
* In the `depth.FM()` function it has been renamed the argument "dfunc" by "dfunc2".
`subset.fdata()` is a wrapper function of subset function.


# fda.usc 1.2.1

* The functions `dcor.xy`, `dcor.test()`, `bcdcor.dist()` and `dcor.dist()` 
(wrapper function of energy package) are added.
`fregre.gsam()` function can be used without smoothed vairables.
* A bug in "selec" argument of `summary.fregre.gkam()` has been fixed.
* A bug in "h" argument of fregre.plm() has been fixed.
* A bug in sigma="vexponential" of `rproc2fdata()` has been fixed. 
The default values depth.RPp, depth.RP and `rproc2fdata()` have been modified.
* A `classif.DD()` function uses the same bandwidth "h" for k groups in modal depth
and same projections "proj" for k groups in RP depth.
* A bug in `predict.fregre.gsam()` when PLS are previously estimated using norm=TRUE has been fixed.


# fda.usc 1.2.0

New functions:

 + `classif.DD()`, fits Nonparametric Classification Procedure Based on DD-plot (depth-versus-depth plot) for G groups.
 + `depth.FMp()`, `depth.modep` and `depth.RPp()` functions provide the depth measure for a list of p--functional data objects.
 + `metric.ldata()`, computes distance for a list of p--functional data objects.
 + `metric.hausdorff()`, computes hausdorff distance.
 + Multivariate depth functions have been renamed: "depth." to "mdepth.".


# fda.usc 1.1.0

New functions:

 + `fregre.basis.fr()` fits functional response model.
 + `metric.kl()` computes Kullback--Leibler distance.
 + `anova.onefactor()`: tests one--way anova model for functional data.
 + `split.fdata()`, `unlist.fdata()`: A wrapper functions of the split and unlist function for functional data.
 + `func.mean.formula()` computes the mean curve for the each level of grouping variable.

New dataset: Mithochondiral calcium overload (MCO) data set.

New utilities:

 + `fdata()` converts arrays of 3 dimension in a functional data of 2 dimension `plot.fdata()` allows functional data of 2 dimension.
 + The functions fdata2ppc(), fdata2ppls(), fregre.ppc(), fregre.ppls(), 
 fregre.ppc.cv(), fregre.ppls.cv() are deprecated in favor of `fdata2pc()`, `fdata2pls()`, `fregre.pc()`, `fregre.pc.cv()`, `fregre.pls()`, `fregre.pls.cv()`. These latter functions include penalty arguments. 


# fda.usc 1.0.5

* "pls" package dependency has been removed.
* A bug in `outlier.ltr()` function has been fixed in the case of the rownames (of fdata) can not be converted to numeric values.
* A bug in `fregre.lm()` function has been fixed in the case of one of the covariates is a factor and penalization argument is required (rn or lambda greater than zero). 
* The new argument "lambda" in `fregre.lm()` penalizes the derivative of second order of the functional data.
* New arguments in `predict.fregre.fd()` and `predict.fregre.lm()` produce confidence or prediction intervals at the specified level mimicking `predict.lm()`.
* New argument "verbose" in min.basis(), min.np() and `rproc2fdata()`.
* The argument "mu" in `rproc2fdata()` allows vector and also fdata class object.


# fda.usc 1.0.4

* A bug in `fregre.pls()` function has been fixed in the case of the FPLS basis are created with the argument norm is TRUE (the curves are centred and scaled). 
* A bug in `outliers.depth.trim()` function has been fixed in the case of the procedure requires more than one iteration step.


# fda.usc 1.0.3 

* This version introduces new function `classif.depth()` that fits a nonparametric classification procedure based on maximum depth measure.
* Penalized FPC an FPLS basis are computed in `create.pc.basis()` and `create.pls.basis()` by the new arguments "lambda" and "P".
* A bug in `predict.fregre.fd()` function has been fixed (in the case of the "object" is fitted using funtional partial least square basis). 


# fda.usc 1.0.2 

* Release 1.0.2 introduces new functions:
* New argument "se.fit" in predict.fregre.fd() and predict.fregre.lm() function.
* A bug in `CV.S()` function when "y" argument is a fdata object has been fixed. 
* This "bug" has involved the min.np() function.
* In `metric.lp()` function NA values are returned, if the fdata has NA's values. 
* In `metric.lp()` function supremum distance is computed, if "lp" argument is 0.


# fda.usc 1.0.1 

New functions:

* New depth functions and its corresponding shortcut functions (see `help(Descriptive)` form more details):
 + `depth.SD()` provides the simplicial depth measure for bivariate data.
 + `depth.PD()` provides the depth measure using random projections for multivariate data.
 + `depth.MhD()` provides the Mahalanobis depth measure for multivariate data.
 + `depth.HD()` provides the half-space depth measure for multivariate data.

* It introduces a new functions for functional PC and PLS regression:
+  fregre.ppc, fregre.ppls, fregre.ppc.cv, fregre.ppls.cv,
and the auxiliary functions:  fdata2ppc, fdata2ppls, P.penalty.
* The function rber.gold() has been renamed by `rwild()` function. 
* ow, r`wild()` contructs the Wild bootstrap residuals. 
* `order.fdata()` is a wrapper function of order function.

New arguments and options:

*New arguments "wild" and "type.wild" in `fregre.bootstrap()`.
* In `fregre.glm()`, `fregre.gsam()`, classif.glm2boost(), classif.gsam2boost() the
"fdataobj" argument allows a multivariate data or functional data.
* `fregre.lm()` allows penalization by "rn" parameter (ridge regression).
* `fregre.pc()` and fregre.basis() allow weighted least squares by "weights" argument.

* Correction of bugs:
 + A bug in "draw" argument of `fdata.bootstrap()` has been fixed.
 + A bug in `predict.classif()` function ussing a fitted object by `classif.knn()` has been fixed.
 + A bug in `create.pc.basis()` function when "l" argument has length 1 has been fixed. 
 + This "bug" has involved the following functions: `fregre.lm()`, `fregre.glm()`
 and `fregre.gsam()`.


# fda.usc 1.0.0

* Release 1.0.0 was released in Oct. 2012 as the working version to
accompany 'Febrero-Bande, M. and  Oviedo de la Fuente, M. (2012).
'Statistical Computing in Functional Data Analysis: The R Package fda.usc.
Journal of Statistical Software, 51(4), 1-28., URL https://www.jstatsoft.org/article/view/v051i04

* New functions added: 
  +  `depth.RT()` implements a random Tukey depth (RT) and its corresponding  some shortcut functions.
    `metric.dist()` methods as wrappers for dist(), 
  + New arguments in `depth.FM()`, `depth.RP()`, `depth.mode()` and `fdata()`.

* The functions: `fregre.glm()`, `fregre.gsam()`, `fregre.gkam()`, `classif.np()`,
`classif.glm()`, `classi.glm()` allow functional and multivariate analysis. 


# fda.usc 0.9.8.1

Release 0.9.8.1 introduces new functions flm.Ftest() and 
dfv.test(). The first performs a functional F-test and the 
second implements the test of Delsol, Ferraty and Vieu (2010). 

Function `flm.test()` now has a better computational performance 
and function Aijr() has been replaced by `Adot()`.

New argument "lambda" in `fdata2fd()` function.

New argument "rn" in `create.pc.basis()` function.

fregre.kgam() has been renamed to `fregre.gkam()`.


# fda.usc 0.9.8

Release 0.9.8 introduces a new function `flm.test()` that 
allows to test for the Functional Linear Model with scalar 
response for a given dataset. Is based on the new functions 
`PCvM.statistic()`, Aijr() and rber.gold().

A bug in fregre.kgam() has been fixed.


# fda.usc 0.9.7

* New functions: 
 + fregre.kgam(), classif.kgam(), dev.S(),
 + predict.fregre.kgam(), print.fregre.kgam(), 
 + summary.fregre.kgam(), `fregre.gsam()`, 
 + `classif.np()`, classif.kgam(), `classif.gsam()`.

* New argument "par.S" in: `fregre.np()`, `fregre.np.cv()`, `fregre.plm()`,
`S.NW()`, `S.KNN()`, `S.LLR()`.

* New attributes for: `metric.lp()`, `semimetric.basis()` and `semimetric.NPFDA()`

# fda.usc 0.9.6

* Release 0.9.6 renames the functions: 
 + pc.fdata()-->`fdata2pc()`
 +   pls.fdata()-->`fdata2pls()`
 + pc.cor()-->`summary.fdata.comp()`
 + pc.fdata()-->`summary.fdata.comp()`

* It added `create.pls.basis()`, `Math.fdata()`,  `Ops.fdata()`, `Summary.fdata()` and
`dis.cos.cor()` function.

* New argument par.S in: `fregre.np()`, `fregre.np.cv()`, `fregre.plm()`,
New argument cv in: `S.NW()`, `S.KNN()`, `S.LLR()`

* In `metric.lp()` the argument p now is called lp.


# fda.usc 0.9.5

Release 0.9.5 improves `fdata.bootstrap()` function
(better computational efficiency). It introduces a new 
functions: for Partial Linear Square (pls.fdata(), `fregre.pls() `
and `fregre.pls.cv()`) and Simpson integration (int.simpson()
and int.simpson2()). It modifies the functions `metric.lp()`, 
`inprod.fdata()`, `summary.fregre.fd()` and `predict.fregre.fd()`.

# fda.usc 0.9.4

Release 0.9.4 added 3 script files: Outliers_fdata.R, 
flm_beta_estimation_brownian_data.R and Classif_phoneme.R.
It has introduced the functions `fregre.glm()` and 
`predict.fregre.glm()` which allow fit and predict
respectively Functional Generalized Linear Models.
It has introduced the functions create.pc.basis and 
`create.fdata.basis()`  which allow to create basis objects
for functional data of class "fdata".

# fda.usc 0.9

Release 0.9 introduces a new function `h.default()` that 
simplifies the calculation of the bandwidth parameter
"h" in the functions: `fregre.np()`, `fregre.np.cv()` and `fregre.plm()`.  
In most of the functions has added a stop control when the
dataset has missing data (NA's). It adds the attribute "call"
to the distance  matrix calculated in `metric.lp()`, 
`semimetric.basis()`  and `semimetric.NPFDA()` functions.
