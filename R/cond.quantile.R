#' @name cond.quantile
#' @title Conditional quantile
#' 
#' @description Computes the quantile for conditional distribution function.
#' 
#' @param qua Quantile value, by default the median (\code{qua}=0.5).
#' @param fdata0 Conditional functional explanatory data of \code{\link{fdata}}
#' class object.
#' @param fdataobj Functional explanatory data of \code{\link{fdata}} class
#' object.
#' @param y Scalar Response.
#' @param fn Conditional distribution function.
#' @param a Lower limit.
#' @param b Upper limit.
#' @param tol Tolerance.
#' @param iter.max Maximum iterations allowed, by default \code{100}.
#' @param \dots Further arguments passed to or from other methods.
#' 
#' @return Return the quantile for conditional distribution function.
#' 
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' 
#' @seealso See Also as: \code{\link{cond.F}} and \code{\link{cond.mode}}.
#' 
#' @aliases cond.quantile
#' 
#' @references Ferraty, F. and Vieu, P. (2006). \emph{Nonparametric functional
#' data analysis.} Springer Series in Statistics, New York.
#' @keywords distribution
#' @examples
#' \dontrun{
#' n= 100
#' t= seq(0,1,len=101)
#' beta = t*sin(2*pi*t)^2
#' x = matrix(NA, ncol=101, nrow=n)
#' y=numeric(n)
#' x0<-rproc2fdata(n,seq(0,1,len=101),sigma="wiener")
#' x1<-rproc2fdata(n,seq(0,1,len=101),sigma=0.1)
#' x<-x0*3+x1
#' fbeta = fdata(beta,t)
#' y<-inprod.fdata(x,fbeta)+rnorm(n,sd=0.1)
#' 
#' prx=x[1:50];pry=y[1:50]
#' ind=50+1;ind2=51:60
#' pr0=x[ind];pr10=x[ind2]
#' ndist=161
#' gridy=seq(-1.598069,1.598069, len=ndist)
#' ind4=5
#' y0 = gridy[ind4]
#' 
#' # Conditional median
#' med=cond.quantile(qua=0.5,fdata0=pr0,fdataobj=prx,y=pry,fn=cond.F,h=1)
#' 
#' # Conditional CI 95% conditional
#' lo=cond.quantile(qua=0.025,fdata0=pr0,fdataobj=prx,y=pry,fn=cond.F,h=1)
#' up=cond.quantile(qua=0.975,fdata0=pr0,fdataobj=prx,y=pry,fn=cond.F,h=1)
#' print(c(lo,med,up))
#' }
#' 
#' @rdname cond.quantile
#' @export
cond.quantile<-function(qua=0.5,fdata0,fdataobj,y,fn, a=min(y),
                        b=max(y), tol=10^floor(log10(max(y)-min(y))-3),
                        iter.max=100,...){

 if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
  nas1<-is.na.fdata(fdataobj)
tt<-fdataobj$argvals
 if (any(nas1))  stop("fdataobj contain ",sum(nas1)," curves with some NA value \n")
if (!is.fdata(fdata0))  fdata0=fdata(fdata0,tt)
  nas2<-is.na.fdata(fdata0)
 if (any(nas2))  stop("fdata0 contain ",sum(nas2)," curves with some NA value \n")
 if (any(is.na(y)))   stop("y contain ",sum(is.na(y)),"  NA values \n")
  i<-0
  tol.up=qua+tol
  tol.lo=qua-tol
  medio<-(a+b)/2
  rmed=fn(fdata0,medio, fdataobj, y,...)$Fc
  rlo=fn(fdata0,a, fdataobj, y, ...)$Fc
  rup=fn(fdata0,b, fdataobj, y, ...)$Fc
  res<-c(i,medio,rmed)
  if ((rup < qua) | (rlo >qua) | (a>b)) {
     print('Error in input data')
     return(NA)
     }
  else {
    while (abs(rmed-qua)>tol & i<iter.max) {
      i<-i+1
      if (rmed<qua) {
        a<-medio
        medio<-(a+b)/2
        rlo=rmed
      }
      else {
        b<-medio
        medio<-(a+b)/2
        rup=rmed
      }
      res<-rbind(res,c(i,medio,rmed))
      rmed=fn(fdata0,medio, fdataobj, y,...)$Fc
  }
  }
print(paste("Fc=",medio,"with Conditional quantile=",qua,", tol=",tol," and iter=",i) )
return(medio)
}

