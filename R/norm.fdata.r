#' @title Approximates Lp-norm for functional data.
#' 
#' @description Approximates Lp-norm for functional data (fdata) object using metric or
#' semimetric functions. Norm for functional data using by default Lp-metric.
#' 
#' @details By default it computes the L2-norm with \code{p = 2} and weights \code{w}
#' with length=\code{(m-1)}.  \deqn{Let \ \ f(x)= fdataobj(x)\ }{}
#' \deqn{\left\|f\right\|_p=\left ( \frac{1}{\int_{a}^{b}w(x)dx} \int_{a}^{b}
#' \left|f(x)\right|^{p}w(x)dx \right)^{1/p}}{}
#' 
#' The observed points on each curve are equally spaced (by default) or not.
#' 
#' @aliases norm.fdata norm.fd
#' @param fdataobj \code{\link{fdata}} class object.
#' @param fdobj Functional data or curves of \link[fda]{fd} class.
#' @param metric Metric function, by default \code{\link{metric.lp}}.
#' @param \dots Further arguments passed to or from other methods.
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' @seealso See also \code{\link{metric.lp}} and \code{\link{norm}}\cr
#' Alternative method: \code{\link{inprod}} of fda-package
#' @keywords math
#' @examples
#' \dontrun{
#' x<-seq(0,2*pi,length=1001)
#' fx1<-sin(x)/sqrt(pi)
#' fx2<-cos(x)/sqrt(pi)
#' argv<-seq(0,2*pi,len=1001)
#' fdat0<-fdata(rep(0,len=1001),argv,range(argv))
#' fdat1<-fdata(fx1,x,range(x))
#' metric.lp(fdat1)
#' metric.lp(fdat1,fdat0)
#' norm.fdata(fdat1)
#' # The same
#' integrate(function(x){(abs(sin(x)/sqrt(pi))^2)},0,2*pi)
#' integrate(function(x){(abs(cos(x)/sqrt(pi))^2)},0,2*pi)
#' 
#' bspl1<- create.bspline.basis(c(0,2*pi),21)
#' fd.bspl1 <- fd(basisobj=bspl1)
#' fd.bspl2<-fdata2fd(fdat1,nbasis=21)
#' norm.fd(fd.bspl1)
#' norm.fd(fd.bspl2)
#' }
#' 
#' @rdname  norm.fdata
#' @export 
norm.fdata <- function(fdataobj,metric=metric.lp,...){
if (!inherits(fdataobj,"fdata")) stop("No fdata class")
if (is.vector(fdataobj$data))    fdataobj$data=matrix(fdataobj$data,nrow=1)
z0<-matrix(0,ncol=ncol(fdataobj),nrow=1)
z0<-fdata(z0,fdataobj[["argvals"]],fdataobj[["rangeval"]],fdataobj[["names"]])
n.lp<-metric(fdataobj,z0,...)
res<-n.lp[,1]
names(res)<-rownames(fdataobj$data)
res
}

#' @rdname  norm.fdata
#' @export 
norm.fd <- function(fdobj){
if (is.fd(fdobj)) rng<- fdobj[[2]]$rangeval
else if (is.basis(fdobj)) rng<- fdobj$rangeval
else stop("Non fd or basis class")
sqrt(inprod(fdobj,fdobj))#/diff(rng))

}


