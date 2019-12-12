#' @title Data-driven sampling of random directions guided by sample of functional
#' data
#' 
#' @description Generation of random directions based on the principal components \eqn{\hat
#' e_1,\ldots,\hat e_k}{\hat e_1,...,\hat e_k} of a sample of functional data
#' \eqn{X_1,\ldots,X_n}{X_1,...,X_n}. The random directions are sampled as
#' \deqn{h=\sum_{j=1}^kh_j\hat e_j,}{h=\sum_{j=1}^kh_j\hat e_j,} with
#' \eqn{h_j\sim\mathcal{N}(0, \sigma_j^2)}{h_j~N(0, \sigma_j^2)},
#' \eqn{j=1,\ldots,k}{j=1,...,k}. Useful for sampling non-orthogonal random
#' directions \eqn{h}{h} such that they are non-orthogonal for the random
#' sample. 
#' 
#' @param n number of curves to be generated.
#' @param X.fdata an \code{\link{fdata}} object used to compute the
#' functional principal components.
#' @param ncomp if an integer vector is provided, the index for the principal
#' components to be considered. If a threshold between \code{0} and \code{1} is
#' given, the number of components \eqn{k}{k} is determined automatically as
#' the minimum number that explains at least the \code{ncomp} proportion of the
#' total variance of \code{X.fdata}.
#' @param fdata2pc.obj output of \code{\link{fdata2pc}} containing as
#' many components as the ones to be selected by \code{ncomp}. Otherwise, it is
#' computed internally.
#' @param sd if \code{0}, the standard deviations \eqn{\sigma_j} are estimated
#' by the standard deviations of the scores for \eqn{e_j}. If not, the
#' \eqn{\sigma_j}'s are set to \code{sd}.
#' @param zero.mean whether the projections should have zero mean. If not, the
#' mean is set to the mean of \code{X.fdata}.
#' @param norm whether the samples should be L2-normalized or not.
#' @return A \code{\link{fdata}} object with the sampled directions.
#' @author Eduardo Garcia-Portugues (\email{edgarcia@@est-econ.uc3m.es}) and
#' Manuel Febrero-Bande (\email{manuel.febrero@@usc.es}).
#' @examples
#' \dontrun{
#' # Simulate some data
#' set.seed(345673)
#' X.fdata <- r.ou(n = 200, mu = 0, alpha = 1, sigma = 2, t = seq(0, 1, l = 201), 
#'                 x0 = rep(0, 200))
#' pc <- fdata2pc(X.fdata, ncomp = 20)
#' 
#' # Samples
#' set.seed(34567)
#' rdir.pc(n = 5, X.fdata = X.fdata, zero.mean = FALSE)$data[, 1:5]
#' set.seed(34567)
#' rdir.pc(n = 5, X.fdata = X.fdata, fdata2pc.obj = pc)$data[, 1:5]
#' 
#' # Comparison for the variance type
#' set.seed(456732)
#' n.proj <- 100
#' set.seed(456732)
#' samp1 <- rdir.pc(n = n.proj, X.fdata = X.fdata, sd = 1, norm = FALSE, ncomp = 0.99)
#' set.seed(456732)
#' samp2 <- rdir.pc(n = n.proj, X.fdata = X.fdata, sd = 0, norm = FALSE, ncomp = 0.99)
#' set.seed(456732)
#' samp3 <- rdir.pc(n = n.proj, X.fdata = X.fdata, sd = 1, norm = TRUE, ncomp = 0.99)
#' set.seed(456732)
#' samp4 <- rdir.pc(n = n.proj, X.fdata = X.fdata, sd = 0, norm = TRUE, ncomp = 0.99)
#' par(mfrow = c(1, 2))
#' plot(X.fdata, col = gray(0.85), lty = 1)
#' lines(samp1[1:10], col = 2, lty = 1)
#' lines(samp2[1:10], col = 4, lty = 1)
#' legend("topleft", legend = c("Data", "Different variances", "Equal variances"), 
#'        col = c(gray(0.85), 2, 4), lwd = 2)
#' plot(X.fdata, col = gray(0.85), lty = 1)
#' lines(samp3[1:10], col = 5, lty = 1)
#' lines(samp4[1:10], col = 6, lty = 1)
#' legend("topleft", legend = c("Data", "Different variances, normalized", 
#'        "Equal variances, normalized"), col = c(gray(0.85), 5:6), lwd = 2)
#' 
#' # Correlations (stronger with different variances and unnormalized; 
#' # stronger with lower ncomp)
#' ind <- lower.tri(matrix(nrow = n.proj, ncol = n.proj))
#' median(abs(cor(sapply(1:n.proj, function(i) inprod.fdata(X.fdata, samp1[i]))))[ind])
#' median(abs(cor(sapply(1:n.proj, function(i) inprod.fdata(X.fdata, samp2[i]))))[ind])
#' median(abs(cor(sapply(1:n.proj, function(i) inprod.fdata(X.fdata, samp3[i]))))[ind])
#' median(abs(cor(sapply(1:n.proj, function(i) inprod.fdata(X.fdata, samp4[i]))))[ind])
#' 
#' # Comparison for the threshold
#' samp1 <- rdir.pc(n = 100, X.fdata = X.fdata, ncomp = 0.25, fdata2pc.obj = pc)
#' samp2 <- rdir.pc(n = 100, X.fdata = X.fdata, ncomp = 0.50, fdata2pc.obj = pc)
#' samp3 <- rdir.pc(n = 100, X.fdata = X.fdata, ncomp = 0.90, fdata2pc.obj = pc)
#' samp4 <- rdir.pc(n = 100, X.fdata = X.fdata, ncomp = 0.95, fdata2pc.obj = pc)
#' samp5 <- rdir.pc(n = 100, X.fdata = X.fdata, ncomp = 0.99, fdata2pc.obj = pc)
#' cols <- rainbow(5, alpha = 0.25)
#' par(mfrow = c(3, 2))
#' plot(X.fdata, col = gray(0.75), lty = 1, main = "Data")
#' plot(samp1, col = cols[1], lty = 1, main = "Threshold = 0.25")
#' plot(samp2, col = cols[2], lty = 1, main = "Threshold = 0.50")
#' plot(samp3, col = cols[3], lty = 1, main = "Threshold = 0.90")
#' plot(samp4, col = cols[4], lty = 1, main = "Threshold = 0.95")
#' plot(samp5, col = cols[5], lty = 1, main = "Threshold = 0.99")
#' 
#' # Normalizing
#' samp1 <- rdir.pc(n = 100, X.fdata = X.fdata, ncomp = 0.50, fdata2pc.obj = pc,
#'                  norm = TRUE)
#' samp2 <- rdir.pc(n = 100, X.fdata = X.fdata, ncomp = 0.90, fdata2pc.obj = pc,
#'                  norm = TRUE)
#' samp3 <- rdir.pc(n = 100, X.fdata = X.fdata, ncomp = 0.95, fdata2pc.obj = pc,
#'                  norm = TRUE)
#' samp4 <- rdir.pc(n = 100, X.fdata = X.fdata, ncomp = 0.99, fdata2pc.obj = pc,
#'                  norm = TRUE)
#' samp5 <- rdir.pc(n = 100, X.fdata = X.fdata, ncomp = 0.999, fdata2pc.obj = pc,
#'                  norm = TRUE)
#' cols <- rainbow(5, alpha = 0.25)
#' par(mfrow = c(3, 2))
#' plot(X.fdata, col = gray(0.75), lty = 1, main = "Data")
#' plot(samp1, col = cols[1], lty = 1, main = "Threshold = 0.50")
#' plot(samp2, col = cols[2], lty = 1, main = "Threshold = 0.90")
#' plot(samp3, col = cols[3], lty = 1, main = "Threshold = 0.95")
#' plot(samp4, col = cols[4], lty = 1, main = "Threshold = 0.99")
#' plot(samp5, col = cols[5], lty = 1, main = "Threshold = 0.999")
#' }
#' @rdname rdir.pc
#' @export 
rdir.pc <- function(n, X.fdata, ncomp = 0.95, fdata2pc.obj = 
                      fdata2pc(X.fdata, ncomp = min(length(X.fdata$argvals), 
                                                             nrow(X.fdata))), 
                    sd = 0, zero.mean = TRUE, norm = FALSE) {
  
  # Check fdata
  if (class(X.fdata) != "fdata") {
    
    stop("X.fdata must be of class fdata")
    
  } 

  # Consider PCs up to a threshold of the explained variance or up to max(ncomp)
  ej <- fdata2pc.obj
  m <- switch((ncomp < 1) + 1,
              max(ncomp),
              max(2, min(which(cumsum(ej$d^2) / sum(ej$d^2) > ncomp))))
  if (ncomp < 1) {
    
    ncomp <- 1:m
    
  }
  
  # Compute PCs with fdata2pc if ej contains less eigenvectors than m
  # The problem is that fdata2pc computes all the PCs and then returns 
  # the eigenvalues (d) for all the components but only the eigenvectors (rotation)
  # for the ncomp components.
  if (nrow(ej$rotation) < m) {
    
    ej <- fdata2pc(X.fdata, ncomp = m)
    
  }
    
  # Standard deviations of the normal coefficients
  if (sd == 0) {
    
    # Standard deviations of scores of X.fdata on the eigenvectors
    sdarg <- apply(ej$x[, ncomp], 2, sd)
    
  } else {
    
    # Constant standard deviation
    sdarg <- rep(sd, length(ncomp))
    
  }

  # Eigenvectors
  eigv <- ej$rotation[ncomp]
  
  # Compute linear combinations of the eigenvectors with coefficients sampled 
  # from a centred normal with standard deviations sdarg
  x <- matrix(rnorm(n * m), ncol = m)
  x <- t(t(x) * sdarg)
  rprojs <- fdata(mdata = x %*% eigv$data, 
                           argvals = argvals(X.fdata))
  
  # Normalize
  if (norm) {

    rprojs$data <- rprojs$data / drop(fda.usc::norm.fdata(rprojs))

  }
  
  # Add mean
  if (!zero.mean) {
    
    rprojs$data <- t(t(rprojs$data) + drop(ej$mean$data))

  }
  
  return(rprojs)
  
}




#' @rdname rp.flm.statistic 
#' @title Statistics for testing the functional linear model using random projections
#' 
#' @description Computes the Cramer-von Mises (CvM) and Kolmogorv-Smirnov (kS) statistics on
#' the projected process \deqn{T_{n,h}(u)=\frac{1}{n}\sum_{i=1}^n (Y_i-\langle
#' X_i,\hat \beta\rangle)1_{\{\langle X_i, h\rangle\leq u\}},}{T_{n,
#' h}(u)=1/n\sum_{i = 1}^n (Y_i - <X_i, \hat \beta>)1_{<X_i, h> \le u},}
#' designed to test the goodness-of-fit of a functional linear model with
#' scalar response.
#' \code{NA}'s are not allowed neither in the functional covariate nor in the
#' scalar response.
#' 
#' @param proj.X matrix of size \code{c(n, n.proj)} containing, for each
#' column, the projections of the functional data \eqn{X_1,\ldots,X_n} into a
#' random direction \eqn{h}. Not required if \code{proj.X.ord} is provided.
#' @param residuals the residuals of the fitted funtional linear model,
#' \eqn{Y_i-\langle X_i,\hat \beta\rangle}{Y_i - <X_i, \hat \beta, Y_i>}.
#' Either a vector of length \code{n} (same residuals for all projections) or a
#' matrix of size \code{c(n.proj, n)} (each projection has an associated set
#' residuals).
#' @param proj.X.ord matrix containing the row permutations of \code{proj.X}
#' which rearranges them increasingly, for each column. So, for example
#' \code{proj.X[proj.X.ord[, 1], 1]} equals \code{sort(proj.X[, 1])}. If not
#' provided, it is computed internally.
#' @param F.code whether to use faster \code{FORTRAN} code or \code{R} code.
#' @return A list containing: 
#' \itemize{ 
#' \item {list("statistic")}{ a matrix of size \code{c(n.proj, 2)} with the the CvM (first column) and KS (second)
#' statistics, for the \code{n.proj} different projections.}
#' \item {list("proj.X.ord")}{the computed row permutations of \code{proj.X},
#' useful for recycling in subsequent calls to \code{rp.flm.statistic} with the
#' same projections but different residuals.} 
#' }
#' @author Eduardo Garcia-Portugues (\email{edgarcia@@est-econ.uc3m.es}) and
#' Manuel Febrero-Bande (\email{manuel.febrero@@usc.es}).
#' @references Cuesta-Albertos, J.A., Garcia-Portugues, E., Febrero-Bande, M.
#' and Gonzalez-Manteiga, W. (2017). Goodness-of-fit tests for the functional
#' linear model based on randomly projected empirical processes.
#' arXiv:1701.08363. \url{https://arxiv.org/abs/1701.08363}
#' @examples
#' \dontrun{
#' # Simulated example
#' set.seed(345678)
#' t <- seq(0, 1, l = 101)
#' n <- 100
#' X <- r.ou(n = n, t = t)
#' beta0 <- fdata(mdata = cos(2 * pi * t) - (t - 0.5)^2, argvals = t,
#'                rangeval = c(0,1))
#' Y <- inprod.fdata(X, beta0) + rnorm(n, sd = 0.1)
#' 
#' # Linear model
#' mod <- fregre.pc(fdataobj = X, y = Y, l = 1:3)
#' 
#' # Projections
#' proj.X1 <- inprod.fdata(X, r.ou(n = 1, t = t))
#' proj.X2 <- inprod.fdata(X, r.ou(n = 1, t = t))
#' proj.X12 <- cbind(proj.X1, proj.X2)
#' 
#' # Statistics
#' t1 <- rp.flm.statistic(proj.X = proj.X1, residuals = mod$residuals)
#' t2 <- rp.flm.statistic(proj.X = proj.X2, residuals = mod$residuals)
#' t12 <- rp.flm.statistic(proj.X = proj.X12, residuals = mod$residuals)
#' t1$statistic
#' t2$statistic
#' t12$statistic
#' 
#' # Recycling proj.X.ord
#' rp.flm.statistic(proj.X.ord = t1$proj.X.ord, residuals = mod$residuals)$statistic
#' t1$statistic
#' 
#' # Sort in the columns
#' cbind(proj.X12[t12$proj.X.ord[, 1], 1], proj.X12[t12$proj.X.ord[, 2], 2]) -
#' apply(proj.X12, 2, sort)
#' 
#' # FORTRAN and R code
#' rp.flm.statistic(proj.X = proj.X1, residuals = mod$residuals)$statistic -
#' rp.flm.statistic(proj.X = proj.X1, residuals = mod$residuals, 
#'                  F.code = FALSE)$statistic
#' 
#' # Matrix and vector residuals
#' rp.flm.statistic(proj.X = proj.X12, residuals = mod$residuals)$statistic
#' rp.flm.statistic(proj.X = proj.X12, 
#'                  residuals = rbind(mod$residuals, mod$residuals * 2))$statistic
#' }
#' 
#' @export 
rp.flm.statistic <- function(proj.X, residuals, proj.X.ord = NULL, F.code = TRUE) {
  
  # Number of projections
  n.proj <- ifelse(is.null(proj.X.ord), ncol(proj.X), ncol(proj.X.ord))
  
  # Residuals as a matrix
  if (!is.matrix(residuals)) {
    
    residuals <- matrix(residuals, nrow = n.proj, ncol = length(residuals),
                        byrow = TRUE)
    
  }
  n <- ncol(residuals)
  if (nrow(residuals) != n.proj) {
    
    stop("The number of rows in residuals must be the number of projections")
    
  }
  
  # Matrix of statistics (columns) projected in n.proj projections (rows)
  rp.stat <- matrix(0, nrow = n.proj, ncol = 2)
  
  # Order projections if not provided
  if (is.null(proj.X.ord)) {
    
    proj.X.ord <- apply(proj.X, 2, order)
    
  }
  
  # Compute statistics
  if (F.code) {
    
    # Statistic
    rp.stat <- .Fortran("rp_stat", proj_X_ord = proj.X.ord, residuals = residuals,
                        n_proj = n.proj, n = n, rp_stat_proj = rp.stat,
                        PACKAGE = "fda.usc")$rp_stat_proj
    
  } else {
    
    # R implementation
    for (i in 1:n.proj) {
      
      # Empirical process
      y <- cumsum(residuals[i, proj.X.ord[, i]])
      
      # Statistics (CvM and KS, rows)
      CvM <- sum(y^2)
      KS <- max(abs(y))
      rp.stat[i, ] <- c(CvM, KS)
      
    }
    
    # Standardize
    rp.stat[, 1] <- rp.stat[, 1] / (n^2)
    rp.stat[, 2] <- rp.stat[, 2] / sqrt(n)
    
  }
  
  # Return both statistics
  colnames(rp.stat) <- c("CvM", "KS")
  return(list(statistic = rp.stat, proj.X.ord = proj.X.ord))
  
}



#' @rdname rp.flm.test
#' @title Goodness-of fit test for the functional linear model using random
#' projections
#' 
#' @description Tests the composite null hypothesis of a Functional Linear Model with scalar
#' response (FLM), \deqn{H_0:\,Y=\langle
#' X,\beta\rangle+\epsilon\quad\mathrm{vs}\quad H_1:\,Y\neq\langle
#' X,\beta\rangle+\epsilon.}{H_0: Y = <X, \beta> + \epsilon vs H_1: Y != <X,
#' \beta> + \epsilon.} If \eqn{\beta=\beta_0}{\beta=\beta_0} is provided, then
#' the simple hypothesis \eqn{H_0:\,Y=\langle X,\beta_0\rangle+\epsilon}{H_0: Y
#' = <X, \beta_0> + \epsilon} is tested. The way of testing the null hypothesis
#' is via a norm (Cramer-von Mises or Kolmogorov-Smirnov) in the empirical
#' process indexed by the projections.
#' 
#' No NA's are allowed neither in the functional covariate nor in the scalar
#' response.
#' 
#' @param X.fdata functional observations in the class
#' \code{\link{fdata}}.
#' @param Y scalar responses for the FLM. Must be a vector with the same number
#' of elements as functions are in \code{X.fdata}.
#' @param beta0.fdata functional parameter for the simple null hypothesis, in
#' the \code{\link{fdata}} class. The \code{argvals} and
#' \code{rangeval} arguments of \code{beta0.fdata} must be the same of
#' \code{X.fdata}. If \code{beta0.fdata=NULL} (default), the function will test
#' for the composite null hypothesis.
#' @param B number of bootstrap replicates to calibrate the distribution of the
#' test statistic.
#' @param n.proj vector with the number of projections to consider.
#' @param est.method estimation method for \eqn{\beta}{\beta}, only used in the
#' composite case. There are three methods: 
#' \itemize{ 
#' \item {list("\"pc\"")}{ if \code{p} is given, then \eqn{\beta}{\beta} is estimated by
#' \code{\link{fregre.pc}}. Otherwise, \code{p} is chosen using \code{\link{fregre.pc.cv}} and the \code{p.criterion} criterion.}
#' \item {list("\"pls\"")}{ if \code{p} is given, \eqn{\beta}{\beta} is estimated by \code{\link{fregre.pls}}. 
#' Otherwise, \code{p} is chosen using \code{\link{fregre.pls.cv}} and the \code{p.criterion} criterion.}
#' \item {list("\"basis\"")}{ if \code{p} is given, \eqn{\beta}{\beta} is estimated by \code{\link{fregre.basis}}. 
#' Otherwise, \code{p} is' chosen using \code{\link{fregre.basis.cv}} and the \code{p.criterion} criterion. 
#' Both in \code{\link{fregre.basis}} and \code{\link{fregre.basis.cv}}, the same basis for
#' \code{basis.x} and \code{basis.b} is considered.} 
#' }
#' @param p number of elements for the basis representation of
#' \code{beta0.fdata} and \code{X.fdata} with the \code{est.method} (only
#' composite hypothesis). If not supplied, it is estimated from the data.
#' @param p.criterion for \code{est.method} equal to \code{"pc"} or
#' \code{"pls"}, either \code{"SIC"}, \code{"SICc"} or one of the criterions
#' described in \code{\link{fregre.pc.cv}}. For \code{"basis"} a value
#' for \code{type.CV} in \code{\link{fregre.basis.cv}} such as
#' \code{GCV.S}.
#' @param pmax maximum size of the basis expansion to consider in when using
#' \code{p.criterion}.
#' @param type.basis type of basis if \code{est.method = "basis"}.
#' @param projs a \code{\link{fdata}} object containing the random
#' directions employed to project \code{X.fdata}. If numeric, the convenient
#' value for \code{ncomp} in \code{\link{rdir.pc}}.
#' @param verbose whether to show or not information about the testing
#' progress.
#' @param same.rwild wether to employ the same wild bootstrap residuals for
#' different projections or not.
#' @param ... further arguments passed to \link[fda]{create.basis} (not
#' \code{rangeval} that is taken as the \code{rangeval} of \code{X.fdata}).
#' @return An object with class \code{"htest"} whose underlying structure is a
#' list containing the following components: 
#' \itemize{
#' \item {list("p.values.fdr")}{ a matrix of size \code{c(n.proj, 2)}, containing
#' in each row the FDR p-values of the CvM and KS tests up to that projection.}
#' \item {list("proj.statistics")}{ a matrix of size \code{c(max(n.proj), 2)}
#' with the value of the test statistic on each projection.}
#' \item {list("boot.proj.statistics")}{ an array of size \code{c(max(n.proj), 2,
#' B)} with the values of the bootstrap test statistics for each projection.}
#' \item {list("proj.p.values")}{ a matrix of size \code{c(max(n.proj), 2)}}
#' \item {list("method")}{ information about the test performed and the kind of
#' estimation performed.} 
#' \item {list("B")}{ number of bootstrap replicates used.} 
#' \item {list("n.proj")}{ number of projections specified}
#' \item {list("projs")}{ random directions employed to project \code{X.fdata}.}
#' \item {list("type.basis")}{ type of basis for \code{est.method = "basis"}.}
#' \item {list("beta.est")}{ estimated functional parameter \eqn{\hat \beta}{\hat
#' \beta} in the composite hypothesis. For the simple hypothesis, \code{beta0.fdata}.} 
#' \item {list("p")}{ number of basis elements considered for estimation of \eqn{\beta}{\beta}.} 
#' \item {list("p.criterion")}{ criterion employed for selecting \code{p}.} 
#' \item {list("data.name")}{ the character string "Y = <X, b> + e"} 
#' }
#' @author Eduardo Garcia-Portugues (\email{edgarcia@@est-econ.uc3m.es}) and
#' Manuel Febrero-Bande (\email{manuel.febrero@@usc.es}).
#' 
#' @references Cuesta-Albertos, J.A., Garcia-Portugues, E., Febrero-Bande, M.
#' and Gonzalez-Manteiga, W. (2017). Goodness-of-fit tests for the functional
#' linear model based on randomly projected empirical processes.
#' arXiv:1701.08363. \url{https://arxiv.org/abs/1701.08363}
#' 
#' Garcia-Portugues, E., Gonzalez-Manteiga, W. and Febrero-Bande, M. (2014). A
#' goodness-of-fit test for the functional linear model with scalar response.
#' Journal of Computational and Graphical Statistics, 23(3), 761--778.
#' \url{http://dx.doi.org/10.1080/10618600.2013.812519}
#' 
#' @examples
#' \dontrun{
#' # Simulated example
#' 
#' set.seed(345678)
#' t <- seq(0, 1, l = 101)
#' n <- 100
#' X <- r.ou(n = n, t = t, alpha = 2, sigma = 0.5)
#' beta0 <- fdata(mdata = cos(2 * pi * t) - (t - 0.5)^2, argvals = t,
#'                rangeval = c(0,1))
#' Y <- inprod.fdata(X, beta0) + rnorm(n, sd = 0.1)
#' 
#' # Test all cases
#' rp.flm.test(X.fdata = X, Y = Y, est.method = "pc")
#' rp.flm.test(X.fdata = X, Y = Y, est.method = "pls")
#' rp.flm.test(X.fdata = X, Y = Y, est.method = "basis", 
#'             p.criterion = fda.usc::GCV.S)
#' rp.flm.test(X.fdata = X, Y = Y, est.method = "pc", p = 5)
#' rp.flm.test(X.fdata = X, Y = Y, est.method = "pls", p = 5)
#' rp.flm.test(X.fdata = X, Y = Y, est.method = "basis", p = 5)
#' rp.flm.test(X.fdata = X, Y = Y, beta0.fdata = beta0)
#' 
#' # Composite hypothesis: do not reject FLM
#' rp.test <- rp.flm.test(X.fdata = X, Y = Y, est.method = "pc")
#' rp.test$p.values.fdr
#' pcvm.test <- flm.test(X.fdata = X, Y = Y, est.method = "pc", B = 1e3,
#'                       plot.it = FALSE)
#' pcvm.test
#' 
#' # Estimation of beta
#' par(mfrow = c(1, 3))
#' plot(X, main = "X")
#' plot(beta0, main = "beta")
#' lines(rp.test$beta.est, col = 2)
#' lines(pcvm.test$beta.est, col = 3)
#' plot(density(Y), main = "Density of Y", xlab = "Y", ylab = "Density")
#' rug(Y)
#' 
#' # Simple hypothesis: do not reject beta = beta0
#' rp.flm.test(X.fdata = X, Y = Y, beta0.fdata = beta0)$p.values.fdr
#' flm.test(X.fdata = X, Y = Y, beta0.fdata = beta0, B = 1e3, plot.it = FALSE)
#' 
#' # Simple hypothesis: reject beta = beta0^2
#' rp.flm.test(X.fdata = X, Y = Y, beta0.fdata = beta0^2)$p.values.fdr
#' flm.test(X.fdata = X, Y = Y, beta0.fdata = beta0^2, B = 1e3, plot.it = FALSE)
#' 
#' # Tecator dataset
#' 
#' # Load data
#' data(tecator)
#' absorp <- tecator$absorp.fdata
#' ind <- 1:129 # or ind <- 1:215
#' x <- absorp[ind, ]
#' y <- tecator$y$Fat[ind]
#' 
#' # Composite hypothesis
#' rp.tecat <- rp.flm.test(X.fdata = x, Y = y, est.method = "pc")
#' pcvm.tecat <- flm.test(X.fdata = x, Y = y, est.method = "pc", B = 1e3,
#'                        plot.it = FALSE)
#' rp.tecat$p.values.fdr[c(5, 10), ]
#' pcvm.tecat
#' 
#' # Simple hypothesis
#' zero <- fdata(mdata = rep(0, length(x$argvals)), argvals = x$argvals,
#'               rangeval = x$rangeval)
#' rp.flm.test(X.fdata = x, Y = y, beta0.fdata = zero)
#' flm.test(X.fdata = x, Y = y, beta0.fdata = zero, B = 1e3)
#' 
#' # With derivatives
#' rp.tecat <- rp.flm.test(X.fdata = fdata.deriv(x, 1), Y = y, est.method = "pc")
#' rp.tecat$p.values.fdr
#' rp.tecat <- rp.flm.test(X.fdata = fdata.deriv(x, 2), Y = y, est.method = "pc")
#' rp.tecat$p.values.fdr
#' 
#' # AEMET dataset
#' 
#' # Load data
#' data(aemet)
#' wind.speed <- apply(aemet$wind.speed$data, 1, mean)
#' temp <- aemet$temp
#' 
#' # Remove the 5% of the curves with less depth (i.e. 4 curves)
#' par(mfrow = c(1, 1))
#' res.FM <- depth.FM(temp, draw = TRUE)
#' qu <- quantile(res.FM$dep, prob = 0.05)
#' l <- which(res.FM$dep <= qu)
#' lines(aemet$temp[l], col = 3)
#' 
#' # Data without outliers
#' wind.speed <- wind.speed[-l]
#' temp <- temp[-l]
#' 
#' # Composite hypothesis
#' rp.aemet <- rp.flm.test(X.fdata = temp, Y = wind.speed, est.method = "pc")
#' pcvm.aemet <- flm.test(X.fdata = temp, Y = wind.speed, B = 1e3,
#'                        est.method = "pc", plot.it = FALSE)
#' rp.aemet$p.values.fdr
#' apply(rp.aemet$p.values.fdr, 2, range)
#' pcvm.aemet
#' 
#' # Simple hypothesis
#' zero <- fdata(mdata = rep(0, length(temp$argvals)), argvals = temp$argvals,
#'               rangeval = temp$rangeval)
#' flm.test(X.fdata = temp, Y = wind.speed, beta0.fdata = zero, B = 1e3,
#'          plot.it = FALSE)
#' rp.flm.test(X.fdata = temp, Y = wind.speed, beta0.fdata = zero)
#' }
#' @useDynLib fda.usc, .registration = TRUE
#' @export 
rp.flm.test <- function(X.fdata, Y, beta0.fdata = NULL, B = 1000, n.proj = 10, 
                        est.method = "pc", p = NULL, p.criterion = "SICc", 
                        pmax = 20, type.basis = "bspline", projs = 0.95, 
                        verbose = TRUE, same.rwild = FALSE, ...) {
  
  # Sample size
  n <- dim(X.fdata)[1]
  
  # p data driven flag
  p.data.driven <- is.null(p)
  
  # Display progress
  if (verbose) {
    
    cat("Computing estimation of beta... ")
    
  }
  
  # Truncate maximum basis expansion
  pmax <- min(pmax, n)
  
  ## Estimation of beta
  
  # Composite hypothesis: optimal estimation of beta and the basis expansion
  if (is.null(beta0.fdata)) {
    
    # Center the data first
    X.fdata <- fdata.cen(X.fdata)$Xcen
    Y <- Y - mean(Y)
    
    # Method
    meth <- "Random projection based test for the functional linear model using"
    
    # PC
    if (est.method == "pc") {
      
      # Optimal p by p.criterion
      if (p.data.driven) {
        
        # Method
        meth <- paste(meth, "optimal PC basis representation")
        
        # Choose the number of basis elements
        mod <- fregre.pc.cv(fdataobj = X.fdata, y = Y, kmax = 1:pmax,
                                     criteria = p.criterion)
        p.opt <- length(mod$pc.opt)
        ord.opt <- mod$pc.opt
        
        # Return the best model
        mod <- mod$fregre.pc
        pc.comp <- mod$fdata.comp
        
      # Fixed p
      } else {
        
        # Method
        meth <- paste(meth, " a representation in a PC basis of ", p, "elements")
        
        # Estimation of beta on the given fixed basis
        mod <- fregre.pc(fdataobj = X.fdata, y = Y, l = 1:p)
        pc.comp <- mod$fdata.comp
        p.opt <- p
        ord.opt <- mod$l
        
      }
      
    # PLS
    } else if (est.method == "pls") {
      
      # Optimal p by p.criterion
      if (p.data.driven) {
        
        # Method
        meth <- paste(meth, "optimal PLS basis representation")
        
        # Choose the number of the basis: SIC is probably the best criteria
        mod <- fregre.pls.cv(fdataobj = X.fdata, y = Y, kmax = pmax,
                                      criteria = p.criterion)
        p.opt <- length(mod$pls.opt)
        ord.opt <- mod$pls.opt
        
        # Return the best model
        mod <- mod$fregre.pls
        
      # Fixed p
      } else {
        
        # Method
        meth <- paste(meth, "a representation in a PLS basis of ", p, "elements")
        
        # Estimation of beta on the given fixed basis
        mod <- fregre.pls(fdataobj = X.fdata, y = Y, l = 1:p)
        p.opt <- p
        ord.opt <- mod$l
        
      }
      
    # Deterministic basis
    } else if (est.method == "basis") {
      
      # Optimal p by p.criterion
      if (p.data.driven) {
        
        # Method
        meth <- paste(meth, "optimal", type.basis, "basis representation")
        
        # Choose the number of the bspline basis with GCV.S
        if (type.basis == "fourier") {
          
          basis.x <- seq(1, pmax, by = 2)
          
        } else {
          
          basis.x <- 5:max(pmax, 5)
          
        }
        mod <- fregre.basis.cv(fdataobj = X.fdata, y = Y, 
                                        basis.x = basis.x, basis.b = NULL,
                                        type.basis = type.basis, 
                                        type.CV = p.criterion, verbose = FALSE, 
                                        ...)
        p.opt <- mod$basis.x.opt$nbasis
        ord.opt <- 1:p.opt

      # Fixed p
      } else {
        
        # Method
        meth <- paste(meth, "a representation in a", type.basis, "basis of ",
                      p, "elements")
        
        # Estimation of beta on the given fixed basis
        basis.opt <- do.call(what = paste("create.", type.basis,
                                          ".basis", sep = ""),
                             args = list(rangeval = X.fdata$rangeval,
                                         nbasis = p, ...))
        mod <- fregre.basis(fdataobj = X.fdata, y = Y, 
                                     basis.x = basis.opt, basis.b = basis.opt)
        p.opt <- p
        ord.opt <- 1:p.opt
        
      }
      
    } else {
      
      stop(paste("Estimation method", est.method, "not implemented."))
      
    }
    
    # Estimated beta
    beta.est <- mod$beta.est
    
    # Hat matrix
    H <- mod$H
    
    # Residuals
    e <- mod$residuals
    
  # Simple hypothesis
  } else {
    
    # Method
    meth <- "Random projection based test for the simple hypothesis in a functional linear model"
    
    # Do not need to estimate beta
    beta.est <- beta0.fdata
    p.opt <- NA
    
    # Compute the residuals
    e <- drop(Y - inprod.fdata(X.fdata, beta.est))
    
  }
  
  ## Computation of the statistic
  
  # Fix seed for projections and bootstrap samples
  old <- .Random.seed
  set.seed(987654321)
  on.exit({.Random.seed <<- old})
  
  # Sample random directions
  if (verbose) {
    
    cat("Done.\nComputing projections... ")
    
  }
  if (is.numeric(projs)) {
    
    # Compute PCs if not done yet
    if (est.method != "pc" | !is.null(beta0.fdata)) {
      
      pc.comp <- fdata2pc(X.fdata, ncomp = min(length(X.fdata$argvals), 
                                                        nrow(X.fdata)))
      
      
    }
    
    # Random directions
    if (length(n.proj) > 1) {
      
      vec.nproj <- sort(n.proj)
      n.proj <- max(n.proj)
      
    } else {
      
      vec.nproj <- 1:n.proj
      
    }
    projs <- rdir.pc(n = n.proj, X.fdata = X.fdata, ncomp = projs, 
                     fdata2pc.obj = pc.comp, sd = 1)
    
  } else {
    
    n.proj <- length(projs)
    vec.nproj <- 1:n.proj
    
  }
  
  # Compute projections for the statistic and the bootstrap replicates
  proj.X <- inprod.fdata(X.fdata, projs) # A matrix n x n.proj
  
  # Statistic
  rp.stat <- rp.flm.statistic(proj.X = proj.X, residuals = e, F.code = TRUE)
  
  ## Bootstrap calibration
  
  # Define required objects
  rp.stat.star <- array(NA, dim = c(n.proj, 2, B))
  if (verbose) {
    
    cat("Done.\nBootstrap calibration...\n ")
    pb <- txtProgressBar(style = 3)
    
  }
  
  # Composite hypothesis
  if (is.null(beta0.fdata)) {
    
    # Calculate the matrix that gives the residuals of the linear model 
    # from the observed response. This allows to resample efficiently the 
    # residuals without re-estimating again the beta
    Y.to.residuals.matrix <- diag(rep(1, n)) - H
    
    # Bootstrap resampling
    for (i in 1:B) {
      
      # Generate bootstrap errors
      if (same.rwild) {
        
        e.star <- matrix(rwild(e, "golden"), nrow = n.proj, ncol = n,
                        byrow = TRUE)
        
      } else {
        
        e.star <- matrix(rwild(rep(e, n.proj), "golden"), nrow = n.proj,
                        ncol = n, byrow = TRUE)
        
      }

      # Residuals from the bootstrap estimated model (implicit column recycling)
      Y.star <- t(t(e.star) + drop(Y - e))
      e.hat.star <- Y.star %*% Y.to.residuals.matrix
      
      # Calculate the bootstrap statistics
      rp.stat.star[, , i] <- rp.flm.statistic(residuals = e.hat.star,
                                              proj.X.ord = rp.stat$proj.X.ord,
                                              F.code = TRUE)$statistic
      
      # Display progress
      if (verbose) {
        
        setTxtProgressBar(pb, i / B)
        
      }
      
    }
    
  # Simple hypothesis
  } else {
    
    # Bootstrap resampling
    for (i in 1:B) {
      
      # Generate bootstrap errors
      if (same.rwild) {
        
        e.hat.star <- matrix(rwild(e, "golden"), nrow = n.proj,
                             ncol = n, byrow = TRUE)
        
      } else {
        
        e.hat.star <- matrix(rwild(rep(e, n.proj), "golden"),
                             nrow = n.proj, ncol = n, byrow = TRUE)
        
      }
      
      # Calculate the bootstrap statistics
      rp.stat.star[, , i] <- rp.flm.statistic(residuals = e.hat.star,
                                              proj.X.ord = rp.stat$proj.X.ord,
                                              F.code = TRUE)$statistic
      
      # Display progress
      if (verbose) {
        
        setTxtProgressBar(pb, i/B)
        
      }
      
    }
    
  }
  
  # Compute the p-values of the projected tests
  positiveCorrection <- FALSE
  pval <- t(sapply(1:n.proj, function(i) {
    
    c(sum(rp.stat$statistic[i, 1] <= rp.stat.star[i, 1, ]), 
      sum(rp.stat$statistic[i, 2] <= rp.stat.star[i, 2, ]))
    
    }) + positiveCorrection) / (B + positiveCorrection)
  
  # Compute p-values depending for the vector of projections
  rp.pvalue <- t(sapply(seq_along(vec.nproj), function(k) {
    
    apply(pval[1:vec.nproj[k], , drop = FALSE], 2, function(x) {
      
      l <- length(x)
      return(min(l / (1:l) * sort(x)))
      
    })
    
  }))
  colnames(rp.pvalue) <- colnames(pval) <- c("CvM", "KS")
  rownames(rp.pvalue) <- vec.nproj
  
  # Return result
  if (verbose) {
    
    cat("\nDone.\n")
    
  }
  options(warn = -1)
  mean.stats <- colMeans(rp.stat$statistic)
  names(mean.stats) <- paste("Mean", names(mean.stats))
  
  result <- structure(list(statistics.mean = mean.stats,
                           p.values.fdr = rp.pvalue,
                           proj.statistics = rp.stat$statistic,
                           boot.proj.statistics = rp.stat.star,
                           proj.p.values = pval, method = meth, B = B,
                           n.proj = vec.nproj, projs = projs,
                           type.basis = type.basis, beta.est = beta.est, 
                           p = p.opt, p.criterion = p.criterion,
                           data.name = "Y = <X, b> + e"))
  class(result) <- "htest"
  return(result)
  
}

