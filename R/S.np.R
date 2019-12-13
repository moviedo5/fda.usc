#' @name S.np 
#' @title Smoothing matrix by nonparametric methods
#' 
#' @description Provides the smoothing matrix \code{S} for the discretization points \code{tt}
#' 
#' @details Options: 
#' \itemize{
#'  \item Nadaraya-Watson kernel estimator (S.NW) with bandwidth parameter \code{h}. 
#'  \item Local Linear Smoothing (S.LLR) with bandwidth parameter \code{h}.
#'  \item K nearest neighbors estimator (S.KNN) with parameter \code{knn}.
#'  \item Polynomial Local Regression Estimator (S.LCR) with parameter of polynomial \code{p} and of kernel \code{Ker}.
#'  \item Local Cubic Regression Estimator (S.LPR) with kernel \code{Ker}.
#'  }
#' @aliases S.np S.LLR S.KNN S.LPR S.LCR S.NW

#' @param tt Vector of discretization points or distance matrix \code{mdist}
#' @param h Smoothing parameter or bandwidth. In S.KNN, number of k-nearest neighbors.
#' @param Ker Type of kernel used, by default normal kernel.
#' @param w Optional case weights.
#' @param cv If \code{TRUE}, cross-validation is done.
#' @param p Polynomial degree.
#' @param \dots Further arguments passed to or from other methods. Arguments to
#' be passed by default to \link[fda]{create.basis}
#' @return Return the smoothing matrix \code{S}.
#' \itemize{
#'  \item {\code{S.LLR}}{ return the smoothing matrix by  Local Linear Smoothing.}
#'  \item { \code{S.NW}}{ return the smoothing matrix by Nadaraya-Watson kernel estimator.}
#'  \item {\code{S.KNN}}{ return the smoothing matrix by k nearest neighbors estimator.}
#'  \item {\code{S.LPR}}{ return the smoothing matrix by Local Polynomial Regression Estimator.}
#'  \item {\code{S.LCR}}{ return the smoothing matrix by Cubic Polynomial Regression.}
#'  }
#'   
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente \email{manuel.oviedo@@usc.es}
#' @seealso See Also as \code{\link{S.basis}}
#' @references
#' Ferraty, F. and Vieu, P. (2006). \emph{Nonparametric functional data analysis.} Springer Series in Statistics, New York. 
#'   
#' Wasserman, L. \emph{All of Nonparametric Statistics}. Springer Texts in Statistics, 2006.
#'   
#' Opsomer, J. D., and Ruppert, D. (1997). Fitting a bivariate additive model by local polynomial regression. \emph{The Annals of Statistics}, 25(1), 186-211.
#' @keywords smooth
#' @examples
#' \dontrun{
#'   tt=1:101
#'   S=S.LLR(tt,h=5)
#'   S2=S.LLR(tt,h=10,Ker=Ker.tri)
#'   S3=S.NW(tt,h=10,Ker=Ker.tri)
#'   S4=S.KNN(tt,h=5,Ker=Ker.tri)
#'   par(mfrow=c(2,3))
#'   image(S)
#'   image(S2)
#'   image(S3)
#'   image(S4)
#'   S5=S.LPR(tt,h=10,p=1, Ker=Ker.tri)
#'   S6=S.LCR(tt,h=10,Ker=Ker.tri)
#'   image(S5)
#'   image(S6)
#' }
#' @rdname S.np
#' @export 
S.LLR<-function (tt, h, Ker = Ker.norm,w=NULL,cv=FALSE)
{
 if (is.matrix(tt)) {
    if (ncol(tt)!=nrow(tt)) {
      if (ncol(tt)==1) {
         tt=as.vector(tt)
         tt=abs(outer(tt,tt, "-"))}
#      else stop("Error: incorrect arguments passed")
    }}
 else if (is.vector(tt))    tt=abs(outer(tt,tt, "-"))
 else stop("Error: incorrect arguments passed")
 if (cv)  diag(tt)=Inf
 k=Ker(tt/h)
 if (cv)   diag(k)=0
# S1=apply(k*tt,1,sum,na.rm=TRUE)
# S2=apply(k*(tt^2),1,sum,na.rm=TRUE)
 S1=rowSums(k*tt,na.rm=TRUE)
 S2=rowSums(k*(tt^2),na.rm=TRUE)
 b=k*(S2-tt*S1)
 if (cv)   diag(b)=0
if (is.null(w)) w<-rep(1,nrow(b))
 b<-sweep(b,1,w,FUN="*")
# res =b/apply(b,1,sum,na.rm=TRUE)
 res =b/rowSums(b,na.rm=TRUE)
return(res)}

#' @rdname S.np
#' @export 
S.LPR<-function (tt, h, p=1, Ker = Ker.norm, w = NULL, cv = FALSE) 
{
  if (is.matrix(tt)) {
    if (ncol(tt) != nrow(tt)) {
      if (ncol(tt) == 1) {
        tt = as.vector(tt)
        tt = abs(outer(tt, tt, "-"))
      }
    }
  }
  else if (is.vector(tt)) 
    tt = abs(outer(tt, tt, "-"))
  else stop("Error: incorrect arguments passed")
  if (is.null(w)) 
    w <- rep(1, nrow(tt))
  k = Ker(tt/h)/h	
  if (cv) diag(k) = 0
  xx=outer(tt/h,0:p,function(a,b) a^b)
  e=matrix(c(1,rep(0,dim(xx)[3]-1)),ncol=1)
  S=matrix(NA,nrow=nrow(tt),ncol=ncol(tt))
  for (i in 1:nrow(tt)){
    Xx=xx[i,,]
    W=diag(k[i,])
    Saux=t(Xx)%*%W%*%Xx
    S[i,]=t(e)%*%solve(Saux)%*%t(Xx)%*%W
  }
  if (cv) 
    diag(S) = 0
  S <- sweep(S, 1, w, FUN = "*")
  res = S/rowSums(S, na.rm = TRUE)
  return(res)
}


#' @rdname S.np
#' @export 
S.LCR<-function(tt, h, Ker=Ker.norm, w=NULL, cv=FALSE){
  res=S.LPR(tt=tt,h=h,p=3,Ker=Ker.norm,w=w,cv=cv)
  return(res)
}

#' @rdname S.np
#' @export 
S.NW<-function (tt, h=NULL, Ker = Ker.norm,w=NULL,cv=FALSE) {
if (is.matrix(tt)) {
  if (ncol(tt)!=nrow(tt)) {
    if (ncol(tt)==1) {
      tt=as.vector(tt)
      tt=abs(outer(tt,tt, "-"))}
    #else stop("Error: incorrect arguments passed")
  }}
else if (is.vector(tt))    tt=abs(outer(tt,tt, "-"))
else stop("Error: incorrect arguments passed")
if (is.null(h)) {
  h=quantile(tt,probs=0.15,na.rm=TRUE)
  while(h==0) {
    h=quantile(tt,probs=pp,na.rm=TRUE)
    pp<-pp+.05
  }
}
if (cv)  diag(tt)=Inf
tt2<-as.matrix(sweep(tt,1,h,FUN="/"))
k<-Ker(tt2)
#print(any(is.na(tt2)))  
if (is.null(w)) w<-rep(1,len=ncol(tt))
k1<-sweep(k,2,w,FUN="*")
#  S =k1/apply(k1,1,sum)
rw<-rowSums(k1,na.rm = TRUE)
rw[rw==0]<-1e-28
S =k1/rw
return(S)
}

#' @rdname S.np
#' @export 
S.KNN<-function(tt,h=NULL,Ker=Ker.unif,w=NULL,cv=FALSE)      {
  if (is.matrix(tt)) {
    if (ncol(tt)!=nrow(tt)) {
      if (ncol(tt)==1) {
        tt=as.vector(tt)
        tt=abs(outer(tt,tt, "-"))}
      #      else stop("Error: incorrect arguments passed")
    }}
  else if (is.vector(tt))    tt=abs(outer(tt,tt, "-"))
  else stop("Error: incorrect arguments passed")
  numgr=ncol(tt)
  if (is.null(h)) h=floor(quantile(1:numgr,probs=0.05,na.rm=TRUE,type=4))
  else if (h<=0 ) stop("Error: incorrect knn value")
  tol=1e-19
  tol=diff(range(tt)*tol)
  tol=1e-19
  if (cv) diag(tt)=Inf
  vec=apply(tt,1,quantile,probs=((h)/numgr),type=4)+tol
  rr=sweep(tt,1,vec,"/")
  rr=Ker(rr)
  #if (cv) diag(rr)=0
  if (!is.null(w)){ #w<-rep(1,ncol(rr))
    rr<-sweep(rr,2,w,FUN="*") }  ## antes un 2
  #print(colSums(rr,na.rm=TRUE))
  rr=rr/rowSums(rr,na.rm=TRUE)
  return(rr)
}