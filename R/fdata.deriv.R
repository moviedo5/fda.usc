#' Computes the derivative of functional data object.
#' 
#' Computes the derivative of functional data.
#' \itemize{ 
#' \item If method =\emph{"bspline"}, \emph{"exponential"}, \emph{"fourier"}, 
#' \emph{"monomial"} or \emph{"polynomial"}.
#'  \code{fdata.deriv} function creates a basis to represent the functional
#' data. 
#' The functional data are converted to class \code{fd} using the
#' \code{\link{Data2fd}} function and the basis indicated in the \code{method}.
#' Finally, the function calculates the derivative of order \code{nderiv} of
#' curves using \code{\link{deriv.fd}} function.\cr 
#' \item If \code{method}=\emph{"fmm"}, \emph{"periodic"}, \emph{"natural"} or
#' \emph{"monoH.FC"} is used \code{\link{splinefun}} function.
#' \item If \code{method}=\emph{"diff"}, raw derivation is applied.  Not recommended to
#' use this method when the values are not equally spaced.\cr 
#' }
#' 
#' @param fdataobj \code{\link{fdata}} class object.
#' @param nderiv Order of derivation, by defalult \code{nderiv}=1.
#' @param method Type of derivative method, for more information see
#' \bold{details}.
#' @param class.out Class of functional data returned: \code{fdata} or
#' \code{fd} class.
#' @param nbasis Number of Basis for \code{fdatataobj\$DATA}. It is only used
#' if method =\emph{"bspline"}, \emph{"exponential"}, \emph{"fourier"},
#' \emph{"monomial"} or \emph{"polynomial"}
#' @param \dots Further arguments passed to or from other methods.
#' 
#' @return Returns the derivative of functional data of \code{fd} class if
#' \code{class.out}="\emph{fd}" or  \code{fdata} class if \code{class.out}="\emph{fdata}".
#' @seealso See also \code{\link{deriv.fd}}, \code{\link{splinefun}} and
#' \code{\link{fdata}}
#' @keywords manip
#' @examples
#' data(tecator)
#' absorp=tecator$absorp.fdata
#' tecator.fd1=fdata2fd(absorp)
#' tecator.fd2=fdata2fd(absorp,"fourier",9)
#' tecator.fd3=fdata2fd(absorp,"fourier",nbasis=9,nderiv=1)
#' #tecator.fd1;tecator.fd2;tecator.fd3
#' tecator.fdata1=fdata(tecator.fd1)
#' tecator.fdata2=fdata(tecator.fd2)
#' tecator.fdata3=fdata(tecator.fd3)
#' tecator.fdata4=fdata.deriv(absorp,nderiv=1,method="bspline",
#' class.out='fdata',nbasis=9)
#' tecator.fd4=fdata.deriv(tecator.fd3,nderiv=0,class.out='fd',nbasis=9)
#' plot(tecator.fdata4)
#' plot(fdata.deriv(absorp,nderiv=1,method="bspline",class.out='fd',nbasis=11))
#' 
#' @export
fdata.deriv<-function(fdataobj,nderiv=1,method="bspline",
                      class.out='fdata',nbasis=NULL,...) {
 if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
 nas1<-is.na(fdataobj)
 if (any(nas1))  stop("fdataobj contain ",sum(nas1)," curves with some NA value \n")
 DATA<-fdataobj[["data"]]
 tt=fdataobj[["argvals"]]
 rtt=fdataobj[["rangeval"]]
 labels=fdataobj[["names"]]
 nr <- nrow(DATA)
 nc <- ncol(DATA)
 if (method=="diff") {
  res=matrix(NA, nrow=nr, ncol = nc)
  for (i in 1:nr) {
   a=diff(DATA[i,],differences=nderiv)/(tt[2:nc]-tt[1:(nc-1)])^nderiv
   ab=matrix(NA,ncol=nc,nrow=2)
   ab[1,2:nc]=a
   ab[2,1:(nc-1)]=a
   res[i,]=colMeans(ab,na.rm=TRUE)
  }
  labels$main<-paste("d(",labels$main,",",nderiv,")",sep="")
  res<-fdata(res,tt,rtt,names=labels)
 }
  else {
      if (any(method==c("fmm", "periodic", "natural", "monoH.FC"))) {
       res=matrix(NA,nrow=nrow(DATA),ncol=ncol(DATA))
       for (i in 1:nrow(DATA)) {
         f1<-splinefun(x=tt,y=DATA[i,],method=method)
         res[i,]=f1(tt,deriv=nderiv)
        }
          labels$main<-paste("d(",labels$main,",",nderiv,")",sep="")
         res<-fdata(res,tt,rtt,names=labels)
       }
  else{
    if (any(method==c("bspline","exponential", "fourier",
      "monomial","polynomial"))) {
      #no run  "constant","polygonal","power"
      res=fdata2fd(fdataobj=fdataobj,type.basis=method,nbasis=nbasis,nderiv=nderiv,...)
      if (class.out=='fdata') {
         ff<-eval.fd(tt,res)
         labels$ylab<-paste("d(",labels$ylab,",",nderiv,")",sep="")
         res=fdata(t(ff),tt,rtt,names=labels)
         }
      }
  }
}
res
}
