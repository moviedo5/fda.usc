#' Sampling Variance estimates
#' 
#' Sampling variance or error variance estimates for regression estimates.
#' 
#' @aliases Var.y Var.e
#' @param y \code{\link{fdata}} class object.
#' @param S Smoothing matrix calculated by \code{\link{S.basis}} or
#' \code{\link{S.NW}} functions.
#' @param Var.e Error Variance Estimates. If \code{Var.e}=NULL, Var.e is
#' calculated.
#' @return Var.y: returns the sampling variance of the functional data. Var.e:
#' returns the sampling error variance of the functional data.
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' @seealso See Also as \code{\link{Var.e}}
#' @references Ramsay, James O., and Silverman, Bernard W. (2006), \emph{
#' Functional Data Analysis}, 2nd ed., Springer, New York.
#' @examples
#' a1<-seq(0,1,by=.01)
#' a2=rnorm(length(a1),sd=0.2)
#' f1<-(sin(2*pi*a1))+rnorm(length(a1),sd=0.2)
#' nc<-50
#' np<-length(f1)
#' tt=1:101
#' mdata<-matrix(NA,ncol=np,nrow=nc)
#' for (i in 1:nc) mdata[i,]<- (sin(2*pi*a1))+rnorm(length(a1),sd=0.2)
#' mdata<-fdata(mdata,tt)
#' S=S.NW(tt,h=0.15)
#' var.e<-Var.e(mdata,S)
#' var.y<-Var.y(mdata,S)
#' var.y2<-Var.y(mdata,S,var.e) #the same
#' @export
Var.y=function(y,S,Var.e=NULL){
  df=fdata.trace(S)
  if (is.null(Var.e)) {
     if (is.vector(y)) {
     	  n=length(y)
     	  y.est=S%*%y
        se=sum((y-y.est)^2,na.rm=TRUE)/(n-df)
	      Var.e=se*diag(ncol(S))
       }
     else {
     if (!is.fdata(y)) y<-fdata(y)
        y<-y[["data"]]
        n=ncol(S)
        y.est<-t(S%*%t(y))
	      Var.e<-t(y-y.est)%*%(y-y.est)/(n-df)
  	}
  }
  var.y=S%*%Var.e%*%t(S)
  return(var.y)
}
