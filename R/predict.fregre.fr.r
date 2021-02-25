#' Predict method for functional response model
#' 
#' Computes predictions for regression between functional explanatory variables
#' and functional response.
#' 
#' @aliases predict.fregre.fr
#' @param object \code{fregre.fr} object.
#' @param new.fdataobj New functional explanatory data of \code{fdata} class.
#' @param \dots Further arguments passed to or from other methods.
#' @return Return the predicted functional data.
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' @seealso See Also as: \code{\link{fregre.basis.fr}}
#' @keywords regression
#' @examples 
#' \dontrun{ 
#' # CV prediction for CandianWeather data
#' rtt<-c(0, 365)
#' basiss  <- create.bspline.basis(rtt,7)
#' basist  <- create.bspline.basis(rtt,9)
#' nam<-dimnames(CanadianWeather$dailyAv)[[2]]
#' 
#' # fdata class (raw data)
#' tt<-1:365
#' tempfdata<-fdata(t(CanadianWeather$dailyAv[,,1]),tt,rtt)
#' log10precfdata<-fdata(t(CanadianWeather$dailyAv[,,3]),tt,rtt)
#' rng<-range(log10precfdata) 
#' for (ind in 1:35){
#'  res1<-  fregre.basis.fr(tempfdata[-ind], log10precfdata[-ind],
#'  basis.s=basiss,basis.t=basist)
#'  pred1<-predict(res1,tempfdata[ind])
#'  plot( log10precfdata[ind],col=1,ylim=rng,main=nam[ind])
#'  lines(pred1,lty=2,col=2)
#'  Sys.sleep(1)
#' }
#' 
#' # fd class  (smooth data)
#' basis.alpha  <- create.constant.basis(rtt)
#' basisx  <- create.bspline.basis(rtt,65)
#' 
#' dayfd<-Data2fd(day.5,CanadianWeather$dailyAv,basisx)
#' tempfd<-dayfd[,1]
#' log10precfd<-dayfd[,3]
#' for (ind in 1:35){
#'  res2 <-  fregre.basis.fr(tempfd[-ind], log10precfd[-ind],
#'  basis.s=basiss,basis.t=basist)
#'  pred2<-predict(res2,tempfd[ind])
#'  plot(log10precfd[ind],col=1,ylim=range(log10precfd$coef),main=nam[ind]) 
#'  lines(pred2,lty=2,col=2)
#'  Sys.sleep(.5)
#' }
#' }
#' 
#' @export
predict.fregre.fr<-function(object,new.fdataobj=NULL,...){
if (is.null(object)) stop("No fregre.fd object entered")
if (is.null(new.fdataobj)) return(object$fitted.values)
if (object$call[[1]]=="fregre.basis.fr" || object$call[[1]]=="fregre.basis.fr.cv"){
beta.est<-object$coefficients
isfdx<-is.fd(new.fdataobj)
if (isfdx) {
  xcoef<-new.fdataobj$coef
  ncurves<-ncol(xcoef)
  }
else {
 xfdobj<-Data2fd(argvals =new.fdataobj$argvals, y = t(new.fdataobj$data), basisobj = object$basis.s)
 xcoef<-xfdobj$coef
 ncurves<-ncol(xcoef)
 if (any(new.fdataobj$argvals!=object$x$argvals)) stop("Incorrect argvals")
}
H = t(xcoef) %*% object$H 
beta.xest = beta.est %*% t(H)
beta.xfd   = fd(beta.xest, object$basis.t)
if (isfdx) {
 #fitted.values  <- fd(coef=yhat, basisobj=object$y$basis, fdnames=object$y$fdnames)
  yhat = eval.fd(object$argvals.y,object$alpha.est) %*% matrix(1,1,ncurves) + eval.fd(object$argvals.y, beta.xfd)  
  fitted.values  <-smooth.basis(object$argvals.y, yhat, object$y$basis)$fd 
}
else {
  yhat = eval.fd(object$y$argvals,object$alpha.est) %*% matrix(1,1,ncurves) + eval.fd(object$y$argvals, beta.xfd)
  fitted.values<-fdata(t(yhat),new.fdataobj$argvals,new.fdataobj$rangeval,new.fdataobj$names)
}
return(fitted.values)
}
else stop("predict is only implemented for fregre.basis.fr output object")
}
