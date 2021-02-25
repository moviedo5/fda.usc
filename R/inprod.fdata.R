#' Inner products of Functional Data Objects o class (fdata)
#' 
#' @description Computes a inner products of functional data objects of class fdata.
#' 
#' @details By default it uses weights \code{w=1}.  \deqn{ \left\langle fdata1,fdata2
#' \right\rangle=\frac{1}{\int_{a}^{b}w(x)dx} \int_{a}^{b} fdata1(x) *
#' fdata2(x)w(x) dx }{} The observed points on each curve are equally spaced
#' (by default) or not.
#' 
#' @param fdata1 Functional data 1 or curve 1. \code{fdata1$data} with
#' dimension (\code{n1} x \code{m}), where \code{n1} is the number of curves
#' and \code{m} are the points observed in each curve.
#' @param fdata2 Functional data 2 or curve 2. \code{fdata2$data} with
#' dimension (\code{n2} x \code{m}), where \code{n2} is the number of curves
#' and \code{m} are the points observed in each curve.
#' @param w Vector of weights with length \code{m}, If \code{w} = 1
#' approximates the metric Lp by Simpson's rule. By default it uses \code{w} =
#' 1
#' @param \dots Further arguments passed to or from other methods.
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' @seealso See also \link[fda]{inprod} and \code{\link{norm.fdata}}
#' @keywords cluster
#' @examples
#' \dontrun{
#' x<-seq(0,2*pi,length=1001)
#' fx1<-sin(x)/sqrt(pi)
#' fx2<-cos(x)/sqrt(pi)
#' argv<-seq(0,2*pi,len=1001)
#' fdat0<-fdata(rep(0,len=1001),argv,range(argv))
#' fdat1<-fdata(fx1,x,range(x))
#' inprod.fdata(fdat1,fdat1)
#' inprod.fdata(fdat1,fdat1)
#' metric.lp(fdat1)
#' metric.lp(fdat1,fdat0)
#' norm.fdata(fdat1)
#' # The same
#' integrate(function(x){(abs(sin(x)/sqrt(pi))^2)},0,2*pi)
#' integrate(function(x){(abs(cos(x)/sqrt(pi))^2)},0,2*pi)
#' }
#' @export
inprod.fdata=function (fdata1,fdata2=NULL, w = 1, ...)   {
  if (!inherits(fdata1,"fdata")) stop("No fdata class")
  tt1<-fdata1[["argvals"]]
  DATA1<-fdata1[["data"]]
  nas1<-is.na.fdata(fdata1)
  if (any(nas1)) {
     stop("fdata1 contain ",sum(nas1)," curves with some NA value \n")
     }
  if (is.null(fdata2)) {fdata2<-fdata1
  }   else  if (!inherits(fdata2,"fdata")) stop("No fdata class")
  nas2<-is.na.fdata(fdata2)
  if (any(nas2)) {
     stop("fdata2 contain ",sum(nas2)," curves with some NA value \n")
  }
  DATA2<-fdata2[["data"]]
  tt2<-fdata2[["argvals"]]
  if  (sum(tt1!=tt2)!=0) 
    stop("Error: different discretization points in the input data.\n")
  numgr = NROW(DATA1) 
  numgr2 = NROW(DATA2)
  inumgr2 <- 1:numgr2
  equi <- argvals.equi(tt1) # Check if data is equispaced
  np <- length(tt1)
  if ((length(w)!=np) & (length(w) != 1)) {
      stop("DATA ERROR: The weight vector hasn't the length of the functions\n")
  }
  mdist = array(0, dim = c(numgr, numgr2))
  #if (is.null(getDefaultCluster())){
  #if (!("doParallel" %in% rownames(installed.packages()))) {
  #  stop("Please install package 'doParallel'") }
  #requireNamespace("foreach", quietly = TRUE)
  # if (!getDoParRegistered())  ops.fda.usc()
  
  aux <- numeric()          
  par.fda.usc <-      eval(parse(text="fda.usc:::par.fda.usc"), envir=.GlobalEnv)
  ncores <- par.fda.usc$ncores
  int.method <- par.fda.usc$int.method
  inumgr <- icount(numgr)
 # if (ncores==1) {    foreach::registerDoSEQ()  }  
  #globalVariables('i') 
  #globalVariables(names=c("i"), package="fda.usc", add=F)
  mdist <- foreach(icnt =  inumgr, .combine = 'rbind') %dopar% {
                  d0 <-DATA1[icnt,]
                  for (ii in  inumgr2) {
                      f=w*(d0 * DATA2[ii,])
                      aux[ii]=int.simpson2(tt1,f,equi,int.method)
                      }         
                     aux
  }
 #mdist<-data.matrix(mdist,numgr,numgr2)
 dim(mdist)<-c(numgr,numgr2)
 rownames(mdist)=rownames(DATA1)
 colnames(mdist)=rownames(DATA2)
 return(mdist)
} 
# ops.fda.usc()
# inprod.fdata(a1,a2)
# aa<-inprod.fdata(a1,a2)
# aa<-metric.lp(a1,a2)
# #aa<-inprod.fdata2(a1,a2)
# dim(aa)
# Y=inprod.fdata(X,beta0)+rnorm(100,sd=0.1)
# 
# 
# cores <- detectCores()
# cl <- makeCluster(cores)
# registerDoParallel(cl)
