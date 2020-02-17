#' @title Semi-functional partially linear model with scalar response.
#' 
#' @description Computes functional regression between functional (and non functional)
#' explanatory variables and scalar response using asymmetric kernel
#' estimation.
#' 
#' @details An extension of the non-parametric functional regression models is the
#' semi-functional partial linear model proposed in Aneiros-Perez and Vieu
#' (2005). This model uses a non-parametric kernel procedure as that described
#' in \code{\link{fregre.np}}. The output \eqn{y} is scalar. A functional
#' covariate \eqn{X} and a multivariate non functional covariate \eqn{Z} are
#' considered.
#' 
#' \deqn{y =\emph{r(X)}+\sum_{j=1}^{p}{Z_j\beta_j}+\epsilon}{y =
#' r(X)+\sum_(j=1:p) Z_j \beta_j+\epsilon}
#' 
#' The unknown smooth real function \eqn{r} is estimated by means of
#' \deqn{\hat{r}_{h}(X)=\sum_{i=1}^{n}{w_{n,h}(X,X_{i})(Y_{i}-Z_{i}^{T}\hat{\beta}_{h})}}{\hat{r}_{h}(X)=\sum_(i=1:n)
#' w_{n,h}(X,X_i) (Y_i - Z_i^{T} \beta.est_h)} where \eqn{W_h} is the weight
#' function:
#' 
#' \eqn{w_{n,h}(X,X_{i})=\frac{K(d(X,X_i)/h)}{\sum_{j=1}^{n}K(d(X,X_j)/h)}}{w_{n,h}(X,X_i)=
#' K( d( X , X_i )/h ) / \sum_(j=1:n) K( d( X, X_j )/h )} with smoothing
#' parameter \eqn{h}, an asymmetric kernel \eqn{K} and a metric or semi-metric
#' \eqn{d}.  In \code{fregre.plm()} by default \eqn{W_h} is a functional
#' version of the Nadaraya-Watson-type weights (\code{type.S=S.NW}) with
#' asymmetric normal kernel (\code{Ker=AKer.norm}) in \eqn{L_2}\cr
#' (\code{metric=metric.lp} with \code{p=2}). The unknown parameters
#' \eqn{\beta_j} for the multivariate non functional covariates are estimated
#' by means of
#' \eqn{\hat{\beta}_j=(\tilde{Z}_{h}^{T}\tilde{Z}_{h})^{-1}\tilde{Z}_{h}^{T}\tilde{Z}_{h}}{\beta.est_j=(Z_h'Z_h)^{-1}
#' Z_h^{T}Z_h} where \eqn{\tilde{Z}_{h}=(I-W_{h})Z}{Z_h=(I-W_h)Z} with the
#' smoothing parameter \eqn{h}. The errors \eqn{\epsilon} are independent, with
#' zero mean, finite variance \eqn{\sigma^2} and
#' \eqn{E[\epsilon|Z_1,\ldots,Z_p,X(t)]=0}{E[\epsilon|Z_1,...,Z_p,X(t)]=0}.\cr
#' 
#' The first item in the \code{data} list is called \emph{"df"} and is a data
#' frame with the response and non functional explanatory variables, as
#' \code{link{lm}}. If non functional data into the formula then
#' \code{\link{lm}} regression is performed.\cr Functional variable
#' (\code{fdata} or \code{fd} class) is introduced in the second item in the
#' \code{data} list.  If only functional variable into the formula then
#' \code{\link{fregre.np.cv}} is performed.\cr
#' 
#' The function estimates the value of smoothing parameter or the bandwidth
#' \code{h} through Generalized Cross-validation \code{GCV} criteria. It
#' computes the distance between curves using the \code{\link{metric.lp}},
#' although you can also use other metric function. \cr Different asymmetric
#' kernels can be used, see \code{\link{Kernel.asymmetric}}.\cr
#' 
#' @param formula an object of class \code{formula} (or one that can be coerced
#' to that class): a symbolic description of the model to be fitted. The
#' details of model specification are given under \code{Details}.
#' @param data List that containing the variables in the model.
#' @param Ker Type of asymmetric kernel used, by default asymmetric normal
#' kernel.
#' @param h Bandwidth, \code{h>0}. Default argument values are provided as the
#' sequence of length 51 from 2.5\%--quantile to 25\%--quantile of the distance
#' between the functional data, see \code{\link{h.default}}.
#' @param metric Metric function, by default \code{\link{metric.lp}}.
#' @param type.CV Type of cross-validation. By default generalized
#' cross-validation \code{\link{GCV.S}} method.
#' @param type.S Type of smothing matrix \code{S}. By default \code{S} is
#' calculated by Nadaraya-Watson kernel estimator (\code{S.NW}).
#' @param par.CV List of parameters for \code{type.CV}: \code{trim}, the alpha
#' of the trimming\cr and \code{draw=TRUE}.
#' @param par.S List of parameters for \code{type.S}: \code{w}, the weights.
#' @param \dots Further arguments passed to or from other methods.
#' @return 
#' \itemize{
#' \item \code{call}{ The matched call.} 
#' \item \code{fitted.values}{ Estimated scalar response.} 
#' \item \code{residuals}{ \code{y} minus \code{fitted values}.}
#' \item \code{df}{ The residual degrees of freedom.} 
#' \item \code{H}{ Hat matrix.}
#' \item \code{r2}{ Coefficient of determination.} 
#' \item \code{sr2}{ Residual variance.}
#' \item \code{y}{ Scalar response.} 
#' \item \code{fdataobj}{ Functional explanatory data.}
#' \item \code{XX}{ Non functional explanatory data.} 
#' \item \code{mdist}{ Distance matrix between curves.} 
#' \item \code{betah}{ beta coefficient estimated} 
#' \item \code{data}{ List that containing the variables in the model.} 
#' \item \code{Ker}{ Asymmetric kernel used.} 
#' \item \code{h.opt}{ Value that minimizes CV or GCV method.} 
#' \item \code{h}{ Smoothing parameter or bandwidth.}
#' \item \code{data}{ List that containing the
#' variables in the model.} 
#' \item \code{gcv}{ GCV values.} 
#' \item \code{formula}{ formula.}
#' }
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@usc.es}
#' @seealso See Also as: \code{\link{predict.fregre.plm}} and
#' \code{\link{summary.fregre.fd}}\cr Alternative methods:
#' \code{\link{fregre.lm}}, \code{\link{fregre.np}} and
#' \code{\link{fregre.np.cv}}
#' @references Aneiros-Perez G. and Vieu P. (2005). \emph{Semi-functional
#' partial linear regression}.  Statistics & Probability Letters, 76:1102-1110.
#' 
#' Ferraty, F. and Vieu, P. (2006). \emph{Nonparametric functional data
#' analysis.} Springer Series in Statistics, New York.
#' 
#' Hardle, W. \emph{Applied Nonparametric Regression}. Cambridge University
#' Press, 1994.
#' 
#' Febrero-Bande, M., Oviedo de la Fuente, M. (2012).  \emph{Statistical
#' Computing in Functional Data Analysis: The R Package fda.usc.} Journal of
#' Statistical Software, 51(4), 1-28. \url{http://www.jstatsoft.org/v51/i04/}
#' @keywords regression
#' @examples
#' \dontrun{
#' data(tecator)
#' x=tecator$absorp.fdata[1:129]
#' dataf=tecator$y[1:129,]
#' 
#' f=Fat~Water+x
#' ldata=list("df"=dataf,"x"=x)
#' res.plm=fregre.plm(f,ldata)
#' summary(res.plm)
#' 
#' # with 2nd derivative of functional data
#' x.fd=fdata.deriv(x,nderiv=2)
#' f2=Fat~Water+x.fd
#' ldata2=list("df"=dataf,"x.fd"=x.fd)
#' res.plm2=fregre.plm(f2,ldata2)
#' summary(res.plm2)
#' }
#' @export
fregre.plm=function(formula,data,h=NULL,Ker=AKer.norm,metric=metric.lp,
type.CV = GCV.S,type.S=S.NW,par.CV=list(trim=0,draw=FALSE),par.S=list(w=1),...){
 C<-match.call()
 mf <- match.call(expand.dots = FALSE)
 m<-match(c("formula","data","h","Ker","metric","type.CV","type.S","par.CV"),names(mf),0L)
 tf <- terms.formula(formula)
 terms <- attr(tf, "term.labels")
 nt <- length(terms)
 if (attr(tf, "response") > 0) {
        response <- as.character(attr(tf, "variables")[2])
        pf <- rf <- paste(response, "~", sep = "")
    } else pf <- rf <- "~"
 vtab<-rownames(attr(tf,"factors"))
 vnf=intersect(terms,names(data$df))
# vnf2=intersect(vtab[-1],names(data$df)[-1])
 vfunc2=setdiff(terms,vnf)
 vint=setdiff(terms,vtab)
 vfunc=setdiff(vfunc2,vint)
# vnf=c(vnf2,vint)
 off<-attr(tf,"offset")
 kterms=1
 z=list()
 lenvnf=length(vnf)
 ty<-deparse(substitute(type.S))
 ke<-deparse(substitute(Ker))
 if (lenvnf>0) {
# cat(" Non functional variables: ",vnf,"\n")
 for ( i in 1:length(vnf)){
    if (kterms > 1)   pf <- paste(pf, "+", vnf[i], sep = "")
     else pf <- paste(pf, vnf[i], sep = "")
     kterms <- kterms + 1
 }
 if   (attr(tf,"intercept")==0) {pf<- paste(pf,-1,sep="")}
 y=as.matrix(data[[1]][,response],ncol=1)
 n=nrow(y)
 if (length(vfunc)>0) {
  if (!is.fdata(data[[vfunc[1]]])) fdataobj=fdata(data[[vfunc[1]]])
  else fdataobj=data[[vfunc[1]]]
  x.fd<-fdataobj[["data"]]
  tt<-fdataobj[["argvals"]]
  rtt<-fdataobj[["rangeval"]]
  mdist=metric(fdataobj,fdataobj,...)
  if (is.null(h))  h<-h.default(data[[vfunc[1]]],type.S=ty,metric=mdist,Ker=ke)
  lenh <- length(h)
  df=gcv<- array(NA, dim = c(lenh))
  yph <- array(NA, dim = c(nrow(y),lenh))
  H <- array(NA, dim = c(nrow(yph),nrow(y),lenh))
  I=diag(1,ncol=nrow(x.fd),nrow=nrow(x.fd))
  pb=txtProgressBar(min=0,max=lenh,style=3)
  XX=as.matrix(data[[1]][,vnf])
  colnames(XX)=vnf
  for (i in 1:lenh) {
        setTxtProgressBar(pb,i-0.5)
####  antes:
##    kmdist=Ker(mdist/h[i])
##    ww=kmdist/apply(kmdist, 1, sum)
####
#    ww=type.S(mdist,h[i],Ker,cv=FALSE)
       par.S$tt<-mdist
    if (is.null(par.S$Ker))  par.S$Ker<-Ker
    if (is.null(par.S$h))  par.S$h<-h[i]
    ww=do.call(type.S,par.S)
    wh=(I-ww)
    yh=wh%*%y
    xh=wh%*%XX
    betah=solve(t(xh)%*%xh)%*%t(xh)%*%yh
    mh=ww%*%(y-XX%*%betah)
    yph[,i]=XX%*%betah+mh
    e=yph[,i]-drop(y)
    c1=solve(t(xh)%*%xh)%*%t(xh)
    H[,,i]=wh%*%XX%*%c1%*%wh+ww
    df[i]=fdata.trace(H[,,i])
    y.pred3=H[,,i]%*%y
#    gcv[i] <- type.CV(y,H[,,i],trim=par.CV$trim,draw=par.CV$draw,...)
     par.CV$S<-H[,,i]
     par.CV$y<-y
     gcv[i]<- do.call(type.CV,par.CV)
     setTxtProgressBar(pb,i)
    }
close(pb)
  if (all(is.infinite(gcv))) print(" Warning: Invalid range for h")
  l = which.min(gcv)
  df=df[l]+ lenvnf
  h.opt <- h[l]
 	names(gcv)<-h
  yph=yph[,l]
  HH=H[,,l]
  e=drop(y)-drop(yph)
#  names(e)<-rownames(fdataobj)
  sr2 = sum(e^2)/(n - df)
  ycen = y - mean(y)
  r2 = 1 - sum(e^2)/sum(ycen^2)
vcov2=sr2*solve(t(xh)%*%xh)
std.error=sqrt(diag(vcov2))
t.value=betah/std.error
p.value= 2 * pt(abs(t.value),df, lower.tail = FALSE)
result<-cbind(betah,std.error,t.value,p.value)
colnames(result) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
if (lenh>1) {
  if (h.opt==min(h))  cat(" Warning: h.opt is the minimum value of bandwidths
   provided, range(h)=",range(h),"\n")
  else if (h.opt==max(h)) cat(" Warning: h.opt is the maximum value of bandwidths
   provided, range(h)=",range(h),"\n")}
z=list(coefficients=result,vcov=vcov2,r2=r2,residuals=e,sr2=sr2,
formula=formula,h.opt=h.opt,h=h,fdataobj=fdataobj,XX=XX,xh=xh,yh=yh,wh=wh,mdist=mdist,y=y,betah=betah,H=HH,data=data,call=C,fitted.values=yph,gcv=gcv,df=df,m=m,metric=metric,Ker=Ker,type.S=type.S)
class(z)="fregre.plm"
}
else {
  XX=data[[1]][,c(response,vnf)]
  print("Warning: lm regession done, non functional data in the formula")
  if (!is.data.frame(XX)) XX=data.frame(XX)
      z=lm(formula=formula,data=XX,x=TRUE,y=TRUE,...)
      z$formula=formula
      z$data=data    }
}
else {
 print("Warning: fregre.np.cv done, only functional data in the formula")
 if (m[5]==0) {
  if (is.null(h)) h<-h.default(data[[vfunc[1]]],type.S=ty,Ker=ke)
  
  z=fregre.np.cv(data[[vfunc[1]]],data[[1]][,response],h=h, 
  Ker=Ker,metric=metric,type.CV=deparse(substitute(type.CV)),type.S=deparse(substitute(type.S)),par.CV=par.CV,...)
 }
 else {
  a1<-match.fun(C[[m[5]]])
  for (i in 1:length(z$call)) {if (z$call[[i]]=="metric") {z$call[[i]]=metric}}
  z$metric=metric
  }
 z$formula=formula
 z$data=data
}
class(z)<-c("fregre.plm","fregre.fd")
z
}

