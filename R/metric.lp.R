#' Approximates Lp-metric distances for functional data.
#' 
#' @description  Measures the proximity between the functional data and curves approximating
#' Lp-metric. If \code{w = 1} approximates the Lp-metric by Simpson's rule. By
#' default it uses \code{lp = 2} and weights \code{w = 1}.
#' 
#' @details By default it uses the L2-norm with \code{lp = 2}.  \deqn{Let \ \ f(x)=
#' fdata1(x)-fdata2(x)}{f(x)= fdata1(x)-fdata2(x)}
#' \deqn{\left\|f\right\|_p=\left ( \frac{1}{\int_{a}^{b}w(x)dx} \int_{a}^{b}
#' \left|f(x)\right|^{p}w(x)dx \right)^{1/p}}{} \cr 
#' The observed points on each curve are equally spaced (by default) or not.
#' 
#' The L\eqn{\infty}-norm is computed with \code{lp = 0}.
#' \deqn{d(fdata1(x),fdata2(x))_{\infty}=sup
#' \left|fdata1(x)-fdata2(x)\right|}
#' 
#' @param fdata1 Functional data 1 or curve 1. If \code{fdata} class, the
#' dimension of \code{fdata1$data} object is (\code{n1} x \code{m}), where
#' \code{n1} is the number of curves and \code{m} are the points observed in
#' each curve. 
#' @param fdata2 Functional data 2 or curve 2. If \code{fdata} class, the
#' dimension of \code{fdata2$data} object is (\code{n2} x \code{m}), where
#' \code{n2} is the number of curves and \code{m} are the points observed in
#' each curve.
#' @param lp Lp norm, by default it uses \code{lp = 2}
#' @param w Vector of weights with length \code{m}, If \code{w = 1}
#' approximates the metric Lp by Simpson's rule. By default it uses \code{w =
#' 1}
#' @param dscale If scale is a numeric, the distance matrix is divided by the
#' scale value. If scale is a function (as the mean for example) the distance
#' matrix is divided by the corresponding value from the output of the
#' function.
#' @param \dots Further arguments passed to or from other methods.
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@usc.es}
#' @seealso See also \code{\link{semimetric.basis}} and
#' \code{\link{semimetric.NPFDA}}
#' @references Febrero-Bande, M., Oviedo de la Fuente, M. (2012).
#' \emph{Statistical Computing in Functional Data Analysis: The R Package
#' fda.usc.} Journal of Statistical Software, 51(4), 1-28.
#' \url{http://www.jstatsoft.org/v51/i04/}
#' @keywords cluster
#' @examples
#' \dontrun{
#' #	INFERENCE PHONDAT
#' data(phoneme)
#' mlearn<-phoneme$learn[1:100]
#' mtest<-phoneme$test[1:100]
#' glearn<-phoneme$classlearn[1:100]
#' gtest<-phoneme$classtest[1:100]
#' # Matrix of distances of curves of DATA1
#' mdist1<-metric.lp(mlearn)
#' 
#' # Matrix of distances between curves of DATA1 and curves of DATA2
#' mdist2<-metric.lp(mlearn,mtest,lp=2)
#' # mdist with L1 norm and weigth=v
#' v=dnorm(seq(-3,3,len=dim(mlearn)[2]))
#' mdist3<-metric.lp(mlearn,mtest,lp=1,w=v)
#' plot(1:100,mdist2[1,],type="l",ylim=c(1,max(mdist3[1,])))
#' lines(mdist3[1,],type="l",col="2")
#' 
#' # mdist with mlearn with different discretization points.
#' # mlearn2=mlearn
#' # mlearn2[["argvals"]]=seq(0,1,len=150)
#' # mdist5<-metric.lp(mlearn,mlearn2)
#' # mdist6<-metric.lp(mlearn2,mlearn) 
#' # sum(mdist5-mdist6)
#' # sum(mdist1-mdist6)
#' 
#' x<-seq(0,2*pi,length=1001)
#' fx<-fdata(sin(x)/sqrt(pi),x)
#' fx0<-fdata(rep(0,length(x)),x)
#' metric.lp(fx,fx0)
#' # The same
#' integrate(function(x){(abs(sin(x)/sqrt(pi))^2)},0,2*pi)
#' }
#' 
#' @export

metric.lp=function (fdata1, fdata2 = NULL, lp = 2, w = 1, dscale=1,...) {
  p <- lp
  C1 <- match.call()
  same <- FALSE
  if (is.fdata(fdata1)) {
    fdat <- TRUE
    tt <- tt1 <- fdata1[["argvals"]]
    rtt <- fdata1[["rangeval"]]
    nas1 <- is.na.fdata(fdata1)
    if (any(nas1)) 
      warning("fdata1 contain ", sum(nas1), " curves with some NA value \n")
    if (is.null(fdata2)) {
      fdata2 <- fdata1
      same <- TRUE
    }
    else if (!is.fdata(fdata2)) {
      fdata2 <- fdata(fdata2, tt1, rtt, fdata1$names)
    }
    nas2 <- is.na.fdata(fdata2)
    if (any(nas2)) {
      warning("fdata2 contain ", sum(nas2), " curves with some NA value \n")
    }
    DATA1 <- fdata1[["data"]]
    DATA2 <- fdata2[["data"]]
    tt2 <- fdata2[["argvals"]]
    if (!same) {
      if (sum(tt1 != tt2) != 0) {
        stop("Error: different discretization points in the input data.\n")
      }
    }
  }  else {
    fdat <- FALSE
    DATA1 <- fdata1
    if (is.null(fdata2)) {
      fdata2 <- fdata1
      same <- TRUE
    }
    DATA2 <- fdata2  }
  numgr = nrow(DATA1)
  numgr2 = nrow(DATA2)
  np <- ncol(DATA1)
  inumgr <- 1:numgr
  inumgr2 <- 1:numgr2
  if ((length(w) != np) & (length(w) != 1)) {
    stop("DATA ERROR: The weight vector hasn't the length of the functions\n")
  }
  testfordim <- sum(dim(DATA1) == dim(DATA2)) == 2
  twodatasets <- TRUE
  etiq1=rownames(DATA1)
  etiq2=rownames(DATA2)    
  if (testfordim) 
    twodatasets <- sum(DATA1 - DATA2, na.rm = TRUE) == 0
 
  # requireNamespace("foreach", quietly = TRUE)
  # requireNamespace("parallel", quietly = TRUE)
  # requireNamespace("doParallel", quietly = TRUE)
  
  par.fda.usc <- eval(parse(text="fda.usc:::par.fda.usc"), envir=.GlobalEnv)
  ncores <- par.fda.usc$ncores
  int.method <- par.fda.usc$int.method
  #if (ncores==1) {
   # foreach::registerDoSEQ()
  #}  
  # if (!foreach::getDoParRegistered()){
  #   if (ncores==1) {
  #     foreach::registerDoSEQ()
  #   }  else{
  #     if (foreach::getDoParWorkers()!=ncores){
  #       # cat("getDoParWorkers != ncores")
  #       cl <-  suppressWarnings(parallel::makePSOCKcluster(ncores))
  #       doParallel::registerDoParallel(cl)
  #     }
  #   }
  #   #ops.fda.usc()
  # }
  inumgr <- icount(numgr)
  if (fdat) {
    #supremum distance    metric.dist(iris[1:20,1:4],method="maximum",
    if (lp==0) mdist<-metric.dist(DATA1,DATA2,method="maximum",...)
    else {  
      equi = argvals.equi(tt)
      mdist = array(0, dim = c(numgr, numgr2))
      predi <- TRUE
      if (testfordim) {
        if (twodatasets) {
          predi <- FALSE
         # print("2 datasets")
            mdist <-foreach(icnt = inumgr, .combine = 'rbind') %dopar% {
            ii = icnt + 1
            aux<-numeric(numgr2)
            for (ii in icnt:numgr2) {
              f = w * abs(DATA1[icnt, ] - DATA2[ii, ])^p
              aux[ii] =int.simpson2(tt, f, equi,int.method)^(1/p)
            }
            aux           }
          mdist = t(mdist) + mdist
        } 
      }
      if (predi) {
        #print("entra parallel1")
        aux<-numeric(numgr2)
        mdist <- foreach(icnt =  inumgr, .combine = 'rbind') %dopar% {
          for (ii in inumgr2) {
            f<-w * abs(DATA1[icnt, ] - DATA2[ii, ])^p
            aux[ii]<-int.simpson2(tt, f, equi,int.method)^(1/p)  }
          aux           
        }
      }
    }
  }
  else {
    if (lp==0) mdist<-metric.dist(DATA1,DATA2,method="maximum",...)
    else {
      mdist = array(0, dim = c(numgr, numgr2))
      aux<-numeric(numgr2)
      #print("entra parallel3")
      #dummy =  function() {
       # globalVariables('icnt')
#dummyfun<-function(iter){iter}
       # registerDoMC(4)
       #  results = foreach(i = 1:10) %dopar% dummyfun(iter = i)
      #}
     #globalVariables('icnt')
     mdist <- foreach(icnt = 1:numgr, .combine = 'rbind') %dopar% {
       for (ii in  inumgr2) {
          #f <- w * abs(DATA1[dummyfun(icnt), ] - DATA2[ii, ])^p
         f <- w * abs(DATA1[icnt, ] - DATA2[ii, ])^p
          aux[ii]<- (sum(f))^(1/p) 
       }
       aux
      }     }
  }
  #  if (stp) stopCluster(cl)
  dim(mdist)=c(numgr,numgr2)
  mdist2<-mdist
  #diag(mdist2)<-NA
  if (is.function(dscale)) 
    dscale<-dscale(mdist2,na.rm=T)
  mdist<-mdist/dscale                           
  rownames(mdist)<-etiq1
  colnames(mdist)<-etiq2    
  attr(mdist, "call") <- "metric.lp"
  attr(mdist, "par.metric") <- list(lp = p, w = w,dscale=dscale)
#if (stp) suppressWarnings(stopCluster(cl))
  return(mdist)
}
###########################################
