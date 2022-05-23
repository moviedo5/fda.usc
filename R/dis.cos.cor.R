#' Proximities between functional data
#' 
#' Computes the cosine correlation distance between two functional dataset of
#' class \code{fdata}.
#' 
#' @param fdata1 Functional data 1 or curve 1.
#' @param fdata2 Functional data 2 or curve 2.
#' @param as.dis Returns the distance matrix como 1-cor(fdata1,fdata2).
#' @return Returns a proximity/distance matrix (depending on \code{as.dis})  between functional data.
#' @seealso See also \code{\link{metric.lp}} and \code{\link{semimetric.NPFDA}}
#' @references Kemmeren P, van Berkum NL, Vilo J, et al. (2002). \emph{Protein
#' Interaction Verification and Functional Annotation by Integrated Analysis of
#' Genome-Scale Data }. Mol Cell. 2002 9(5):1133-43.
#' @keywords cluster
#' @examples
#' \dontrun{
#'  r1<-rnorm(1001,sd=.01)
#'  r2<-rnorm(1001,sd=.01)
#'  x<-seq(0,2*pi,length=1001)
#'  fx<-fdata(sin(x)/sqrt(pi)+r1,x)
#'  dis.cos.cor(fx,fx)
#'  dis.cos.cor(c(fx,fx),as.dis=TRUE)
#'  fx0<-fdata(rep(0,length(x))+r2,x)
#'  plot(c(fx,fx0))
#'  dis.cos.cor(c(fx,fx0),as.dis=TRUE)
#'  }
#'  
#' @export
dis.cos.cor<-function(fdata1,fdata2=NULL,as.dis=FALSE){
if (is.null(fdata2)) {
     a<-inprod.fdata(fdata1)# cambiar por as.vector(norm.fdata(fdata1))
     b1=b2=sqrt(diag(a))
#   b<-norm.fdata(fdata1)
#   print(a/(outer(b[,1],b[,1])))

   }
else {
   a<-inprod.fdata(fdata1,fdata2)
   b1<-norm.fdata(fdata1)
   b2<-norm.fdata(fdata2)
#   print(a/(outer(b1[,1],b2[,1])))
}
   if (as.dis) {
   rho<-1-abs(a/outer(b1,b2))} else {
   rho<-a/outer(b1,b2)}
return(rho)
}
