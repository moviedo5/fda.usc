#' @title Functional regression with scalar response using non-parametric kernel
#' estimation
#' 
#' @description Computes functional regression between functional explanatory variables and
#' scalar response using kernel estimation. 
#' 
#' @details The non-parametric functional regression model can be written as follows \deqn{y_i =r(X_i)+\epsilon_i}{ y
#' = r(X) + \epsilon} where the unknown smooth real function \eqn{r} is
#' estimated using kernel estimation by means of
#' \deqn{\hat{r}(X)=\frac{\sum_{i=1}^{n}{K(h^{-1}d(X,X_{i}))y_{i}}}{\sum_{i=1}^{n}{K(h^{-1}d(X,X_{i}))}}}{\hat{r}(X)=(\sum_i
#' K(d(X,X_i))y_i/h) / (\sum_i K(d(X,X_i)/h)) i=1,...,n} where \eqn{K} is an
#' kernel function (see \code{Ker} argument), \code{h} is the smoothing
#' parameter and \eqn{d} is a metric or a semi-metric (see \code{metric}
#' argument).
#' 
#' The distance between curves is calculated using the \code{\link{metric.lp}}
#' although any other semimetric could be used (see
#' \code{\link{semimetric.basis}} or \code{\link{semimetric.NPFDA}} functions).
#' The kernel is applied to a metric or semi-metrics that provides non-negative
#' values, so it is common to use asymmetric kernels. Different asymmetric
#' kernels can be used, see \code{\link{Kernel.asymmetric}}.\cr
#' 
#' @param fdataobj \code{\link{fdata}} class object.
#' @param y Scalar response with length \code{n}.
#' @param h Bandwidth, \code{h>0}. Default argument values are provided as the
#' 5\%--quantile of the distance between \code{fdataobj} curves, see
#' \code{\link{h.default}}.
#' @param Ker Type of asymmetric kernel used, by default asymmetric normal
#' kernel.
#' @param metric Metric function, by default \code{\link{metric.lp}}.
#' @param type.S Type of smothing matrix \code{S}. By default \code{S} is
#' calculated by Nadaraya-Watson kernel estimator (\code{S.NW}).
#' @param par.S List of parameters for \code{type.S}: \code{w}, the weights.
#' @param \dots Arguments to be passed for \code{\link{metric.lp}} o other
#' metric function.
#' @return Return:
#' \itemize{
#' \item {call}{ The matched call.} 
#' \item {fitted.values}{ Estimated scalar response.} 
#' \item {H}{ Hat matrix.} 
#' \item {residuals}{ \code{y} minus \code{fitted values}.} 
#' \item {df.residual}{ The residual degrees of freedom.} 
#' \item {r2}{ Coefficient of determination.} 
#' \item {sr2}{ Residual variance.} 
#' \item {y}{ Response.} 
#' \item {fdataobj}{ Functional explanatory data.} 
#' \item {mdist}{ Distance matrix between \code{x} and \code{newx}.}
#' \item {Ker}{ Asymmetric kernel used.} 
#' \item {h.opt}{ smoothing parameter or' bandwidth.}
#' }
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' @seealso See Also as: \code{\link{fregre.np.cv}},
#' \code{\link{summary.fregre.fd}} and \code{\link{predict.fregre.fd}} .\cr
#' Alternative method: \code{\link{fregre.basis}},cand \code{\link{fregre.pc}}.
#' @references Ferraty, F. and Vieu, P. (2006). \emph{Nonparametric functional
#' data analysis.} Springer Series in Statistics, New York. \cr
#' 
#' Febrero-Bande, M., Oviedo de la Fuente, M. (2012).  \emph{Statistical
#' Computing in Functional Data Analysis: The R Package fda.usc.} Journal of
#' Statistical Software, 51(4), 1-28. \url{https://www.jstatsoft.org/v51/i04/}
#' 
#' Hardle, W. \emph{Applied Nonparametric Regression}. Cambridge University
#' Press, 1994.
#' @keywords regression
#' @examples
#' \dontrun{
#' data(tecator)
#' absorp=tecator$absorp.fdata
#' ind=1:129
#' x=absorp[ind,]
#' y=tecator$y$Fat[ind]
#' 
#' res.np=fregre.np(x,y,Ker=AKer.epa)
#' summary(res.np)
#' res.np2=fregre.np(x,y,Ker=AKer.tri)
#' summary(res.np2)
#' 
#' # with other semimetrics.
#' res.pca1=fregre.np(x,y,Ker=AKer.tri,metri=semimetric.pca,q=1)
#' summary(res.pca1)
#' res.deriv=fregre.np(x,y,metri=semimetric.deriv)
#' summary(res.deriv)
#' x.d2=fdata.deriv(x,nderiv=1,method="fmm",class.out='fdata')
#' res.deriv2=fregre.np(x.d2,y)
#' summary(res.deriv2)
#' x.d3=fdata.deriv(x,nderiv=1,method="bspline",class.out='fdata')
#' res.deriv3=fregre.np(x.d3,y)
#' summary(res.deriv3)
#' }
#' 
#' @export
fregre.np<-function(fdataobj,y,h=NULL,Ker=AKer.norm,
                    metric=metric.lp,type.S=S.NW,par.S=list(w=1),...){

  if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
isfdata<-is.fdata(y)
nas<-is.na.fdata(fdataobj)
nas.g<-is.na(y)
if (is.null(names(y))) names(y)<-1:length(y)
if (any(nas) & !any(nas.g)) {
   bb<-!nas
   cat("Warning: ",sum(nas)," curves with NA are omited\n")
   fdataobj$data<-fdataobj$data[bb,]
  y<-y[bb]
   }
else {
if (!any(nas) & any(nas.g)) {
   cat("Warning: ",sum(nas.g)," values of group with NA are omited \n")
   bb<-!nas.g
   fdataobj$data<-fdataobj$data[bb,]
     y<-y[bb]
   }
else {
if (any(nas) & any(nas.g))  {
   bb<-!nas & !nas.g
   cat("Warning: ",sum(!bb)," curves  and values of group with NA are omited \n")
   fdataobj$data<-fdataobj$data[bb,]
   y<-y[bb]
   }
}}                              
x<-fdataobj[["data"]]
tt<-fdataobj[["argvals"]]
rtt<-fdataobj[["rangeval"]]
C<-match.call()
mf <- match.call(expand.dots = FALSE)
m<-match(c("fdataobj", "y","h","Ker","metric","type.S","par.S"),names(mf),0L)
#    if (is.vector(x))         x <- t(data.matrix(x))
n = nrow(x)
np <- ncol(x)   
   if (!isfdata) {
   if (n != (length(y)))         stop("ERROR IN THE DATA DIMENSIONS")
   if (is.null(rownames(x)))         rownames(x) <- 1:n
   if (is.null(colnames(x)))         colnames(x) <- 1:np
   if (is.vector(y)) y.mat<-matrix(y,ncol=1)
   else y.mat<-data.matrix(y)
   ny = nrow(y.mat)
   npy <- ncol(y.mat)
   }
   else {
     tty<-y$argvals
     rtty<-y$rangeval
     y.mat<-y$data
     ny = nrow(y.mat)
     npy <- ncol(y.mat)
     if (n != ny | npy!=np)         stop("ERROR IN THE DATA DIMENSIONS")
      }

if (is.matrix(metric)) mdist<-metric
else mdist=metric(fdataobj,fdataobj,...)
ke<-deparse(substitute(Ker))
ty<-deparse(substitute(type.S))
attr(par.S, "call") <- ty
#print(h)
if (is.null(h)) h=h.default(fdataobj,prob=0.05,len=1,metric = mdist, type.S = ty,...)
#     H =type.S(mdist,h,Ker,cv=FALSE)
#     par.S$w<-y
#S.NW2<-function (tt, h, Ker = Ker.norm,cv=FALSE,weights=rep(1,len=length(tt)))
#    print(H[1:2,1:3]);    print(ty)
    par.S$tt<-mdist
    if (is.null(par.S$Ker))  par.S$Ker<-Ker
    if (is.null(par.S$h))  par.S$h<-h
    #if  (type.S=="S.KNN")  par.S$cv<-TRUE
    H=do.call(type.S,par.S)
    par.S$cv<-TRUE
    H.cv=do.call(type.S,par.S)
#    for (j in 1:npy) {
#        y.est[,j]=H%*%y.mat[,j]
#        y.est.cv[,j]=H.cv%*%y.mat[,j]
#        }
   df=trace.matrix(H)
   yp=H%*%y.mat
   yp2<-H.cv%*%y.mat^2-(H.cv%*%y.mat)^2
   if (isfdata) {
      yp<-fdata(yp,tty,rtty,names=y$names)
#      yp.cv<-fdata(y.est.cv,tty,rtty,names=y$names)
      rownames(yp$data)=rownames(y$data)
#      rownames(yp.cv$data)=rownames(y$data)
      ydif<-y-yp
#      ydif.cv<-y-yp.cv
      e<-y-yp
#      ecv<-y-yp.cv
#      sr2=sum(e^2)/(n-df)
      ycen=fdata.cen(y)$Xcen
#  	  r2=1-sum(e^2)/sum(ycen^2)
    norm.e<-norm.fdata(e,metric=metric,...)[,1]^2
    sr2=sum(norm.e)/(n-df)
    ycen=fdata.cen(y)$Xcen
 	  r2=1-sum(norm.e)/sum(ycen^2)
      out<-list("call"=C,"fitted.values"=yp,"H"=H,"residuals"=e,"df.residual"=df,"r2"=r2,
"sr2"=sr2,"y"=y,"fdataobj"=fdataobj,"mdist"=mdist,"Ker"=Ker,
"metric"=metric,"type.S"=type.S,"par.S"=par.S,"h.opt"=h,"m"=m)
      }
else {
    yp<-drop(yp)
    y<-drop(y)
    e<-y-yp
    names(e)<-rownames(x)
    sr2=sum(e^2)/(n-df)
    ycen=y-mean(y)
	  r2=1-sum(e^2)/sum(ycen^2)
	  out<-list("call"=C,"fitted.values"=yp,"H"=H,"residuals"=e,"df.residual"=df,"r2"=r2,
"sr2"=sr2,"y"=drop(y),"fdataobj"=fdataobj,"mdist"=mdist,"Ker"=Ker,
"metric"=metric,"type.S"=type.S,"par.S"=par.S,"h.opt"=h,"m"=m,var.y=yp2)
  }
class(out) <- "fregre.fd"
return(out)
}
