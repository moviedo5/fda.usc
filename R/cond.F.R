#' @title Conditional Distribution Function
#' 
#' @description Calculate the conditional distribution function of a scalar response with
#' functional data.
#' 
#' @details If \code{x.dist=NULL} the distance matrix between \code{fdata} objects is
#' calculated by function passed in \code{metric} argument.
#' 
#' @param fdata0 Conditional explanatory functional data of \code{\link{fdata}}
#' class.
#' @param y0 Vector of conditional response with length \code{n}.
#' @param fdataobj \code{\link{fdata}} class object.
#' @param y Vector of scalar response with length \code{nn}.
#' @param h Smoothing parameter or bandwidth of response \code{y}.
#' @param g Smoothing parameter or bandwidth of explanatory functional data
#' \code{fdataobj}.
#' @param metric Metric function, by default \code{\link{metric.lp}}.
#' @param Ker List of 2 arguments. The fist argument is a character string that
#' determines the type of asymetric kernel (see
#' \code{\link{Kernel.asymmetric}}). Asymmetric Epanechnikov kernel is selected
#' by default. The second argumentis a string that determines the type of
#' integrated kernel(see \code{\link{Kernel.integrate}}). Integrate
#' Epanechnikov kernel is selected by default.\cr.
#' @param \dots Further arguments passed to or from other methods.
#' @aliases cond.F
#' @return
#' \itemize{
#' \item {Fc}{ Conditional distribution function.} 
#' \item {y0}{ Vector of conditional response.} 
#' \item {g}{ Smoothing parameter or bandwidth of explanatory functional data (\code{fdataobj}).} 
#' \item {h}{ Smoothing parameter or bandwidth of respone, \code{y}.} 
#' \item {x.dist}{ Distance matrix between curves of \code{fdataobj} object.} 
#' \item {xy.dist}{ Distance matrix between cuves of \code{fdataobj} and \code{fdata0} objects.}
#' }
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@usc.es}
#' @seealso See Also as: \code{\link{cond.mode}} and
#' \code{\link{cond.quantile}}.
#' @references Ferraty, F. and Vieu, P. (2006). \emph{Nonparametric functional
#' data analysis.} Springer Series in Statistics, New York.
#' @keywords distribution
#' @examples
#' \dontrun{
#' # Read data
#' n= 500
#' t= seq(0,1,len=101)
#' beta = t*sin(2*pi*t)^2
#' x = matrix(NA, ncol=101, nrow=n)
#' y=numeric(n)
#' x0<-rproc2fdata(n,seq(0,1,len=101),sigma="wiener")
#' x1<-rproc2fdata(n,seq(0,1,len=101),sigma=0.1)
#' x<-x0*3+x1
#' fbeta = fdata(beta,t)
#' y<-inprod.fdata(x,fbeta)+rnorm(n,sd=0.1) 
#' prx=x[1:100];pry=y[1:100]
#' ind=101;ind2=102:110
#' pr0=x[ind];pr10=x[ind2,]
#' ndist=61
#' gridy=seq(-1.598069,1.598069, len=ndist)
#' 
#' # Conditional Function
#' res1 = cond.F(pr10, gridy, prx, pry,p=1)
#' res2 = cond.F(pr10, gridy, prx, pry,h=0.3)
#' res3 = cond.F(pr10, gridy, prx, pry,g=0.25,h=0.3)
#' 
#' plot(res1$Fc[,1],type="l",ylim=c(0,1))
#' lines(res2$Fc[,1],type="l",col=2)
#' lines(res3$Fc[,1],type="l",col=3)
#' }
#' 
#' @export
cond.F=function(fdata0,y0,fdataobj,y,h=0.15,g=0.15,metric=metric.lp,
Ker=list(AKer=AKer.epa,IKer=IKer.epa),...){
 if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
  nas1<-is.na.fdata(fdataobj)
tt<-fdataobj$argvals
 if (any(nas1))  stop("fdataobj contain ",sum(nas1)," curves with some NA value \n")
if (!is.fdata(fdata0))  fdata0=fdata(fdata0,tt)
  nas2<-is.na.fdata(fdata0)
 if (any(nas2))  stop("fdata0 contain ",sum(nas2)," curves with some NA value \n")
 if (any(is.na(y0)))  stop("y0 contain ",sum(is.na(y0)),"  NA values \n")
 if (any(is.na(y)))   stop("y contain ",sum(is.na(y)),"  NA values \n")
data<-fdataobj[["data"]]
n = nrow(data)
m = ncol(data)
nn = nrow(fdata0)
ndist=length(y0)
  xy.dist=metric(fdataobj,fdata0,...)
  x.dist=metric(fdataobj,...)
  h3=quantile(x.dist, probs = h, na.rm = TRUE,type=4)
  h=h3
  W =apply(xy.dist/h,2,Ker$AKer)
  g=g*(max(y)-min(y))
  A=outer(y,y0,"-")
  Wy = apply(A/g,2,Ker$IKer) ###
  WW= W
  if (nn==1)  {
      s=sum(WW)
      Fc =1-drop(t(Wy)%*%WW/s)
 }
 else  {
#    s=apply(WW,2,sum,na.rm=TRUE)
    s=colSums(WW,na.rm=TRUE)
    s=diag(1/s,nrow=length(s),ncol=length(s))
    Fc =1-drop((t(Wy)%*%WW)%*%s)
    }
  out=list("Fc"=Fc,"y0"=y0,"g"=g,"h"=h,"x.dist"=x.dist,"xy.dist"=xy.dist)
return(out)
}
