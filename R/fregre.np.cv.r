#' Cross-validation functional regression with scalar response using kernel
#' estimation.
#' 
#' @description Computes functional regression between functional explanatory variables and
#' scalar response using asymmetric kernel estimation by cross-validation
#' method.
#' 
#' @details The non-parametric functional regression model can be written as follows
#' \deqn{ y_i =r(X_i) + \epsilon_i } where the unknown smooth real function
#' \eqn{r} is estimated using kernel estimation by means of
#' \deqn{\hat{r}(X)=\frac{\sum_{i=1}^{n}{K(h^{-1}d(X,X_{i}))y_{i}}}{\sum_{i=1}^{n}{K(h^{-1}d(X,X_{i}))}}}
#' where \eqn{K} is an kernel function (see \code{Ker} argument), \code{h} is
#' the smoothing parameter and \eqn{d} is a metric or a semi-metric (see
#' \code{metric} argument).
#' 
#' The function estimates the value of smoothing parameter (also called
#' bandwidth) \code{h} through Generalized Cross-validation \code{GCV}
#' criteria, see \code{\link{GCV.S}} or \code{\link{CV.S}}.
#' 
#' The function estimates the value of smoothing parameter or the bandwidth
#' through the cross validation methods: \code{\link{GCV.S}} or
#' \code{\link{CV.S}}. It computes the distance between curves using the
#' \code{\link{metric.lp}}, although any other semimetric could be used (see
#' \code{\link{semimetric.basis}} or \code{\link{semimetric.NPFDA}} functions).
#' Different asymmetric kernels can be used, see
#' \code{\link{Kernel.asymmetric}}.\cr
#' 
#' @param fdataobj \code{\link{fdata}} class object.
#' @param y Scalar response with length \code{n}.
#' @param h Bandwidth, \code{h>0}. Default argument values are provided as the
#' sequence of length 25 from 2.5\%--quantile to 25\%--quantile of the distance
#' between \code{fdataobj} curves, see \code{\link{h.default}}.
#' @param Ker Type of asymmetric kernel used, by default asymmetric normal
#' kernel.
#' @param metric Metric function, by default \code{\link{metric.lp}}.
#' @param type.CV Type of cross-validation. By default generalized
#' cross-validation \code{\link{GCV.S}} method.
#' @param type.S Type of smothing matrix \code{S}. By default \code{S} is
#' calculated by Nadaraya-Watson kernel estimator (\code{S.NW}).
#' @param par.CV List of parameters for \code{type.CV}: \code{trim}, the alpha
#' of the trimming\cr and \code{draw=TRUE}.
#' @param par.S List of parameters for \code{type.S}: \code{w}, the weights.
#' @param \dots Arguments to be passed for \code{\link{metric.lp}} o other
#' metric function.
#' @return Return:
#' \itemize{
#' \item \code{call}{ The matched call.} 
#' \item \code{residuals}{ \code{y} minus \code{fitted values}.} 
#' \item \code{fitted.values}{ Estimated scalar response.} 
#' \item \code{df.residual}{ The residual degrees of freedom.} 
#' \item \code{r2}{ Coefficient of determination.} 
#' \item \code{sr2}{ Residual variance.} 
#' \item \code{H}{ Hat matrix.} 
#' \item \code{y}{ Response.} 
#' \item \code{fdataobj}{ Functional explanatory data.}
#' \item \code{mdist}{ Distance matrix between \code{x} and \code{newx}.} 
#' \item \code{Ker}{ Asymmetric kernel used.} 
#' \item \code{gcv}{ CV or GCV values.} 
#' \item \code{h.opt}{ smoothing parameter or bandwidth that minimizes CV or GCV method.} 
#' \item \code{h}{ Vector of smoothing parameter or bandwidth.} 
#' \item \code{cv}{ List with the fitted values and residuals estimated by CV, without the same curve.}
#' }
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' @seealso See Also as: \code{\link{fregre.np}},
#' \code{\link{summary.fregre.fd}} and \code{\link{predict.fregre.fd}} .\cr
#' Alternative method: \code{\link{fregre.basis.cv}} and
#' \code{\link{fregre.np.cv}}.
#' @references Ferraty, F. and Vieu, P. (2006). \emph{Nonparametric functional
#' data analysis.} Springer Series in Statistics, New York.
#' 
#' Hardle, W. \emph{Applied Nonparametric Regression}. Cambridge University
#' Press, 1994.
#' 
#' Febrero-Bande, M., Oviedo de la Fuente, M. (2012).  \emph{Statistical
#' Computing in Functional Data Analysis: The R Package fda.usc.} Journal of
#' Statistical Software, 51(4), 1-28. \url{https://www.jstatsoft.org/v51/i04/}
#' @keywords regression
#' @examples 
#' \dontrun{
#' data(tecator)
#' absorp=tecator$absorp.fdata
#' ind=1:129
#' x=absorp[ind,]
#' y=tecator$y$Fat[ind]
#' Ker=AKer.tri
#' res.np=fregre.np.cv(x,y,Ker=Ker)
#' summary(res.np)
#' res.np2=fregre.np.cv(x,y,type.CV=GCV.S,criteria="Shibata")
#' summary(res.np2)
#' 
#' ## Example with other semimetrics (not run)
#' res.pca1=fregre.np.cv(x,y,Ker=Ker,metric=semimetric.pca,q=1)
#' summary(res.pca1)
#' res.deriv=fregre.np.cv(x,y,Ker=Ker,metric=semimetric.deriv)
#' summary(res.deriv)
#' 
#' x.d2=fdata.deriv(x,nderiv=1,method="fmm",class.out='fdata')
#' res.deriv2=fregre.np.cv(x.d2,y,Ker=Ker)
#' summary(res.deriv2)
#' x.d3=fdata.deriv(x,nderiv=1,method="bspline",class.out='fdata')
#' res.deriv3=fregre.np.cv(x.d3,y,Ker=Ker)
#' summary(res.deriv3)
#' }
#' 
#' @export
fregre.np.cv=function(fdataobj,y,h=NULL,Ker=AKer.norm,metric=metric.lp,
type.CV = GCV.S,type.S=S.NW,par.CV=list(trim=0),par.S=list(w=1),...){
#print("np.CV")
if (is.function(type.CV)) tcv<-deparse(substitute(type.CV))
else tcv<-type.CV
if (is.function(type.S)) ty<-deparse(substitute(type.S))
else ty<-type.S
#print(tcv)
#print(ty)
if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
isfdata<-is.fdata(y)
nas<-is.na.fdata(fdataobj)
nas.g<-is.na(y)
if (is.null(names(y))) names(y) <- seq_len(length(y))
if (any(nas) & !any(nas.g)) {
   bb<-!nas
   if (ops.fda.usc()$warning)
    warning(sum(nas)," curves with NA are omited\n")
   fdataobj$data<-fdataobj$data[bb,]
  y<-y[bb]
   }
else {
if (!any(nas) & any(nas.g)) {
   if (ops.fda.usc()$warning)   
     warning(sum(nas.g)," values of group with NA are omited \n")
   bb<-!nas.g
   fdataobj$data<-fdataobj$data[bb,]
   y<-y[bb]
   }
else {
if (any(nas) & any(nas.g))  {
   bb<-!nas & !nas.g
   if (ops.fda.usc()$warning) 
     warning(sum(!bb)," curves  and values of group with NA are omited \n")
   fdataobj$data<-fdataobj$data[bb,]
   y<-y[bb]
   }
}}
x<-fdataobj[["data"]]
tt<-fdataobj[["argvals"]]
rtt<-fdataobj[["rangeval"]]
   C<-match.call()
   mf <- match.call(expand.dots = FALSE)
   m<-match(c("x", "y","h","Ker","metric","type.CV","type.S","par.CV","par.S"),names(mf),0L)
#   if (is.vector(x))         x <- t(data.matrix(x))
   n = nrow(x)
   np <- ncol(x)
   if (!isfdata) {
   if (n != (length(y)))         stop("ERROR IN THE DATA DIMENSIONS")
   if (is.null(rownames(x)))         rownames(x) <- 1:n
   if (is.null(colnames(x)))         colnames(x) <- 1:np
   if (is.vector(y)) y.mat<-matrix(y,ncol=1)
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
types=FALSE

if (is.matrix(metric)) mdist<-metric
else mdist=metric(fdataobj,fdataobj,...)

attr(par.S, "call") <- ty
#if (!is.function(Ker)) Ker<-get(Ker)


# if (is.character(Ker)){  nker <- function(u,mik=Ker){get(mik)(u)}
# } else {  nker <- function(u,mik=Ker){mik(u)} }

if (is.null(h)) {
  #aa <- strsplit(deparse(substitute(Ker)),"[.]")
#  print(aa)
  #print(deparse(substitute(Ker)))
#nker=get(paste0("Ker.",unlist(aa)[2]))
  h = do.call(h.default,c(list(fdataobj=fdataobj,metric=mdist,
                             prob=c(0.025,0.25),type.S=ty,Ker=Ker),...))
}
else {if   (any(h<=0)) stop("Error: Invalid range for h")}


lenh <- length(h)
cv=gcv1=gcv=cv.error <- array(NA, dim = c(lenh))
par.S2<-par.S
if (is.null(par.S2$h))  par.S$h<-h
if (is.null(par.S$Ker))  par.S$Ker<-Ker
y.est.cv<-y.est<-matrix(NA,nrow=nrow(y.mat),ncol=ncol(y.mat))
par.S$tt<-mdist
par.CV$metric<-metric
  for (i in 1:lenh) {
#print(i)
#print("h")  
#     H2=type.S(mdist,h[i],Ker,cv=FALSE)
    par.S$h<-h[i]
    par.S$cv=TRUE
#    H.cv=do.call(type.S,par.S)
    H.cv=do.call(ty,par.S)
    par.S$cv=FALSE
#    H=do.call(type.S,par.S)
    H=do.call(ty,par.S)

#     gcv[i] <- type.CV(y, H,trim=par.CV$trim,draw=par.CV$draw,...)
    par.CV$S<-switch(tcv,CV.S=H.cv,GCV.S=H,dev.S=H,GCCV.S=H)
#    if (tcv=="CV.S")  par.CV$S<-H.cv
#    if (tcv=="GCV.S") par.CV$S<-H
#    if (tcv=="dev.S") par.CV$S<-H
    for (j in 1:npy) {
        par.CV$y<-y.mat[,j]
#        if (!isfdata) gcv[i]<- do.call(type.CV,par.CV)    #si es fdata no hace falta!!!
#print(j)
#print(tcv)
#print(names(par.CV))
        if (!isfdata) gcv[i]<- do.call(tcv,par.CV)    #si es fdata no hace falta!!!
        y.est[,j]=H%*%y.mat[,j]
        y.est.cv[,j]=H.cv%*%y.mat[,j]
        }
   if (isfdata) {
      par.CV$y<-y
########################
#      gcv[i]<- do.call(type.CV,par.CV)
      gcv[i]<- do.call(tcv,par.CV)
########################
#      calculo directo del CV y GCV (respuesta funcional)
#      yp<-fdata(y.est,tty,rtty,names=y$names)
#      yp.cv<-fdata(y.est.cv,tty,rtty,names=y$names)
#      ydif<-y-yp
#      ydif.cv<-y-yp.cv
#      nmdist1<-norm.fdata(ydif,metric=metric,...)^2
#      gcv1[i]<-sum(nmdist1)/(n*(1-fdata.trace(H)/(n))^2)
#      nmdist2<-norm.fdata(ydif.cv,metric=metric,...)^2
#      cv[i]<- sum(nmdist2)
      }
#     e2=y-H%*%y
#     cv.error[i]=sum(e2^2)/(n-fdata.trace(H))
   }
if (all(is.infinite(gcv)) & ops.fda.usc()$warning
    ) warning(" Warning: Invalid range for h")
   l = which.min(gcv)
   h.opt <- h[l] ######################################################### arreglar
   if (h.opt==min(h) & ops.fda.usc()$warning)
    warning(" Warning: h.opt is the minimum value of bandwidths
   provided, range(h)=",range(h),"\n")
   else if (h.opt==max(h) &    ops.fda.usc()$warning)
    warning(" Warning: h.opt is the maximum value of bandwidths
   provided, range(h)=",range(h),"\n")
   #H =type.S(mdist,h.opt,Ker,cv=FALSE)
  par.S$tt<-mdist
  par.S$h=h.opt
  par.S$cv=FALSE
  H<-do.call(ty,par.S)
 	yp<-H%*%y.mat
  par.S$cv<-TRUE
  Hcv<-do.call(ty,par.S)
 	ypcv<-Hcv%*%y.mat
  df=fdata.trace(H)
	names(gcv)<-h
  if (isfdata) {
  	names(cv)<-h
  	yp<-fdata(yp,tty,rtty)
  	rownames(yp$data)<-rownames(y$data)
    e<-y-yp
    ypcv<-fdata(ypcv,tty,rtty)
   	rownames(ypcv$data)<-rownames(y$data)
    ecv<-y-ypcv
                                                                  
    norm.e<-norm.fdata(e,metric=metric,...)^2
    sr2=sum(norm.e)/(n-df)
    ycen=fdata.cen(y)$Xcen
#	  r2=1-sum(e^2)/sum(ycen^2)
	  r2=1-sum(norm.e)/sum(ycen^2)
    yp2<-Hcv%*%y.mat^2-(Hcv%*%y.mat)^2
    out<-list("call"=C,"fitted.values"=yp,"H"=H,"residuals"=e,"df.residual"=df,"r2"=r2,
"sr2"=sr2,"var.y"=yp2,"y"=y,"fdataobj"=fdataobj,"mdist"=mdist,"Ker"=Ker,
"metric"=metric,"type.S"=type.S,"par.S"=par.S,"gcv"=gcv,"h.opt"=h.opt,"h"=h,"m"=m,
"fit.CV"=list("fitted.values"=ypcv,"residuals"=ecv))
 }
else {
    e<-y-drop(yp)
    names(e)<-rownames(x)
    ecv<-y-drop(ypcv)

    sr2=sum(e^2)/(n-df)
    ycen=y-mean(y)
	  r2=1-sum(e^2)/sum(ycen^2)      	  
    yp2<-Hcv%*%y.mat[,1]^2-(Hcv%*%y.mat[,1])^2
	  out<-list("call"=C,"fitted.values"=yp,"H"=H,"residuals"=e,"df.residual"=df,"r2"=r2,
"sr2"=sr2,"var.y"=yp2,"y"=y,"fdataobj"=fdataobj,"mdist"=mdist,"Ker"=Ker,
"metric"=metric,"type.S"=type.S,"par.S"=par.S,"gcv"=gcv,"h.opt"=h.opt,"h"=h,"m"=m,
"fit.CV"=list("fitted.values"=ypcv,"residuals"=ecv))
  }
class(out)="fregre.fd"
return(out)
}

