#' @title Functional Regression with scalar response using basis representation.
#' 
#' @description Computes functional regression between functional explanatory variable
#' \eqn{X(t)} and scalar response \eqn{Y} using basis representation.
#' 
#' @details
#' \deqn{Y=\big<X,\beta\big>+\epsilon=\int_{T}{X(t)\beta(t)dt+\epsilon}}{Y=<X,\beta>+\epsilon} 
#' where \eqn{ \big< \cdot , \cdot \big>}{<.,.>} denotes the inner product on
#' \eqn{L_2} and \eqn{\epsilon} are random errors with mean zero, finite
#' variance \eqn{\sigma^2} and \eqn{E[X(t)\epsilon]=0}{E[X(t)\epsilon]=0}.
#' 
#' The function uses the basis representation proposed by Ramsay and Silverman (2005) to model the
#' relationship between the scalar response and the functional covariate by
#' basis representation of the observed functional data
#' \eqn{X(t)\approx\sum_{k=1}^{k_{n1}} c_k \xi_k(t)}{X(t)} and the unknown
#' functional parameter \eqn{\beta(t)\approx\sum_{k=1}^{k_{n2}} b_k
#' \phi_k(t)}{\beta(t)}. \cr
#' 
#' The functional linear models estimated by the expression: \deqn{\hat{y}=
#' \big< X,\hat{\beta} \big> =
#' C^{T}\psi(t)\phi^{T}(t)\hat{b}=\tilde{X}\hat{b}}{y.est= < X,\beta.est > =
#' C'\psi\phi' \beta.est=Z\beta.est} where
#' \eqn{\tilde{X}(t)=C^{T}\psi(t)\phi^{T}(t)}{Z=C'\psi\phi'}, and
#' \eqn{\hat{b}=(\tilde{X}^{T}\tilde{X})^{-1}\tilde{X}^{T}y}{\beta.est=(Z'Z)^{-1}Z'y}
#' and so,
#' \eqn{\hat{y}=\tilde{X}\hat{b}=\tilde{X}(\tilde{X}^{T}\tilde{X})^{-1}\tilde{X}^{T}y=Hy}{y.est=Z(Z'Z)^{-1}Z'y}
#' where \eqn{H} is the hat matrix with degrees of freedom: \eqn{df=tr(H)}.\cr
#' 
#' If \eqn{\lambda>0}{\lambda>0} then \code{fregre.basis} incorporates a
#' roughness penalty: \cr
#' \eqn{\hat{y}=\tilde{X}\hat{b}=\tilde{X}(\tilde{X}^{T}\tilde{X}+\lambda
#' R_0)^{-1}\tilde{X}^{T}y= H_{\lambda}y}{y.est=Z(Z'Z+\lambda R_0)^{-1} Z'y=
#' H_{\lambda}y} where \eqn{R_0} is the penalty matrix.\cr
#' 
#' This function allows covariates of class \code{fdata}, \code{matrix},
#' \code{data.frame} or directly covariates of class \code{fd}.  The function
#' also gives default values to arguments \code{basis.x} and \code{basis.b} for
#' representation on the basis of functional data \eqn{X(t)} and the functional
#' parameter \eqn{\beta(t)}, respectively.
#' 
#' If \code{basis=}\code{NULL} creates the \code{bspline} basis by
#' \code{\link{create.bspline.basis}}. \cr If the functional covariate
#' \code{fdataobj} is a matrix or data.frame, it creates an object of class
#' "fdata" with default attributes, see \code{\link{fdata}}.\cr If
#' \code{basis.x$type=``fourier''} and \code{basis.b$type=``fourier''}, the
#' basis are orthonormal and the function decreases the number of fourier basis
#' elements on the \eqn{min(k_{n1},k_{n2})}{min(k.x,k_b)}, where
#' \eqn{k_{n1}}{k.x} and \eqn{k_{n2}}{k.b} are the number of basis element of
#' \code{basis.x} and \code{basis.b} respectively.
#' 
#' @param fdataobj \code{\link{fdata}} class object.
#' @param y Scalar response with length \code{n}.
#' @param basis.x Basis for functional explanatory data \code{fdataobj}.
#' @param basis.b Basis for functional beta parameter.
#' @param lambda A roughness penalty. By default, no penalty \code{lambda=0}.
#' @param Lfdobj See \link[fda]{eval.penalty}.
#' @param weights weights
#' @param \dots Further arguments passed to or from other methods.
#' @return Return:
#' \itemize{
#' \item {call:}{ The matched call.} 
#' \item {coefficients:}{ A named vector of coefficients}
#' \item {residuals:}{ \code{y} minus \code{fitted values}.} 
#' \item {fitted.values:}{ Estimated scalar response.} 
#' \item {beta.est:}{ beta parameter estimated of class \code{fd}} 
#' \item {weights:}{(only for' weighted fits) the specified weights.} 
#' \item {df.residual:}{ The residual degrees of' freedom.} 
#' \item {r2:}{ Coefficient of determination.} 
#' \item {sr2:}{ Residual' variance.} 
#' \item {Vp:}{ Estimated covariance matrix for the parameters.}
#' \item {H:}{ Hat matrix.} 
#' \item {y:}{ Response.}
#' \item {fdataobj:}{ Functional explanatory data of class \code{fdata}.}
#' \item {a.est:}{ Intercept parameter estimated} 
# \item {x.fd:}{ Centered functional explanatory data of class \code{fd}.} 
#' \item {basis.b:}{ Basis used' for beta parameter estimation.} 
#' \item {lambda.opt:}{ A roughness penalty.}
#' \item {Lfdobj:}{ Order of a derivative or a linear differential operator.}
#' \item {P:}{ Penalty matrix.} 
#' \item {lm:}{ Return \code{lm} object }
#' }
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' @seealso See Also as: \code{\link{fregre.basis.cv}},
#' \code{\link{summary.fregre.fd}} and \code{\link{predict.fregre.fd}}.\cr
#' Alternative method: \code{\link{fregre.pc}} and \code{\link{fregre.np}}.
#' @references Ramsay, James O., and Silverman, Bernard W. (2006), \emph{
#' Functional Data Analysis}, 2nd ed., Springer, New York.
#' 
#' Febrero-Bande, M., Oviedo de la Fuente, M. (2012).  \emph{Statistical
#' Computing in Functional Data Analysis: The R Package fda.usc.} Journal of
#' Statistical Software, 51(4), 1-28. \url{https://www.jstatsoft.org/v51/i04/}
#' @keywords regression
#' @examples 
#' \dontrun{
#' # fregre.basis
#' data(tecator)
#' names(tecator)
#' absorp=tecator$absorp.fdata
#' ind=1:129
#' x=absorp[ind,]
#' y=tecator$y$Fat[ind]
#' tt=absorp[["argvals"]]
#' res1=fregre.basis(x,y)
#' summary(res1)
#' basis1=create.bspline.basis(rangeval=range(tt),nbasis=19)
#' basis2=create.bspline.basis(rangeval=range(tt),nbasis=9)
#' res5=fregre.basis(x,y,basis1,basis2)
#' summary(res5)
#' x.d2=fdata.deriv(x,nbasis=19,nderiv=1,method="bspline",class.out="fdata")
#' res7=fregre.basis(x.d2,y,basis1,basis2)
#' summary(res7)
#' }
#' @export
fregre.basis=function(fdataobj,y,basis.x=NULL,basis.b=NULL,lambda=0,
Lfdobj=vec2Lfd(c(0,0),rtt),weights= rep(1,n),...){
R<-NULL
if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
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
  call<-match.call()
  n = nrow(x)
  np <- ncol(x)
  if (n != (length(y))) stop("ERROR IN THE DATA DIMENSIONS")
  if (is.null(rownames(x)))        rownames(x) <- 1:n
  if (is.null(colnames(x)))        colnames(x) <- 1:np
	if (is.null(basis.x))  {
	          nbasis.x=max(floor(np/6),7)
            basis.x=create.bspline.basis(rangeval=rtt,nbasis=nbasis.x,...)
            }
	if (is.null(basis.b))  {
	   if (basis.x$type=="fourier") basis.b<-basis.x
	   else {
     nbasis.b=min(basis.x$nbasis,max(floor(np/10),5))
     basis.b=create.bspline.basis(rangeval=rtt,nbasis=nbasis.b)
     }
  }
  xx<-fdata.cen(fdataobj)
	xmean=xx[[2]]
  xcen=xx[[1]]
	ymean=mean(y)
  ycen=y-ymean
  if (basis.x$type=="fourier" &  basis.x$type=="fourier")
  fou=TRUE
  else fou=FALSE
  if (fou)   {
   if (basis.x$nbasis<basis.b$nbasis) {
 #    cat("Warning: The number of fourier basis elements in basis.b \n
 #    decrease to number of fourier basis elements in the basis.x \n")
     cat("Warning: The nbasis=",basis.b$nbasis," of basis.b must be the same
      of nbasis=",basis.x$nbasis," of basis.x; will be decreased by ",basis.b$nbasis-basis.x$nbasis,sep="","\n")
      basis.b<-basis.x
     }
   if (basis.x$nbasis>basis.b$nbasis) {
    cat("Warning: The nbasis=",basis.x$nbasis," of basis.x must be the same
      of nbasis=",basis.b$nbasis," of basis.b; will be decreased by ",basis.x$nbasis-basis.b$nbasis,sep="","\n")
     basis.x<-basis.b
    }
  }
  J=inprod(basis.x,basis.b)
  ###################  codigo viejo
  # x.fd=Data2fd(argvals=tt,y=t(xcen$data),basisobj=basis.x)
  # C=t(x.fd$coefs)
  # Z<-C%*%J
  ################### codigo nuevo 
   xaux <- fdata2basis(fdataobj,basis.x)
   name.coef <- colnames(xaux$coefs) <- paste(gsub( " ", "", fdataobj$names$main),".",colnames(xaux$coefs),sep="")
   Z <- xaux$coefs
   if (!is.null(basis.b)){
     #mean.list[[vfunc[i]]] <- mean.fd(x.fd);          x.fd <- center.fd(x.fd)
     colnam <- colnames(Z)
     Z <- Z %*% J
     name.coef <- colnames(Z) <- colnam[1:NCOL(Z)]
   }  
   ###################
  vfunc=call[[2]]
  #  XX<-Z
  #Z=cbind(rep(1,len=n),Z)
  #colnames(Z)<-1:ncol(Z)
  #colnames(Z)[2:ncol(Z)]= paste(vfunc,".",basis.b$names, sep = "")
  W<-diag(weights)
if (lambda==0) {
#       S=Z%*%solve(t(Z)%*%Z)%*%t(Z)
    S<-t(Z)%*%W%*%Z
    Lmat    <- chol(S)          
    Lmatinv <- solve(Lmat)
    SS<- Lmatinv %*% t(Lmatinv) 
    S<-Z%*%SS%*%t(Z)%*%W
       yp2=S%*%y
       response="y"
       pf <- paste(response, "~", sep = "")
       for ( i in 1:length(colnames(Z))) pf <- paste(pf, "+", colnames(Z)[i], sep = "")
       # print(pf)
       object.lm=lm(formula=pf,data=data.frame(y,Z,weights),x=TRUE,y=TRUE,weights=weights,...)
       yp=object.lm$fitted.values
       e<-object.lm$residuals
       b.est<-object.lm$coefficients[-1]
       beta.est=fd(b.est,basis.b)
       a.est<-object.lm$coefficients[1]
       df=basis.b$nbasis+1
       rdf<-n-df
       sr2 <- sum(e^2)/ rdf
       Vp<-sr2*SS 
       r2 <- 1 - sum(e^2)/sum(ycen^2)
       r2.adj<- 1 - (1 - r2) * ((n -    1)/ rdf)
       GCV <- sum(e^2)/(n - df)^2
	     coefficients<-object.lm$coefficients
}
else {
    R=diag(0,ncol= basis.b$nbasis+1,nrow=basis.b$nbasis+1)
    #       R[-1,-1]<-getbasispenalty(basis.b,Lfdobj) ############
    R[-1,-1]<-eval.penalty(basis.b,Lfdobj)
    Z=cbind(rep(1,len=n),Z)
    Sb= t(Z)%*%W%*%Z + lambda*R
#    fda:::eigchk(Sb)    
    Lmat <- chol(Sb)          
    Lmatinv <- solve(Lmat)
    Cinv <- Lmatinv %*% t(Lmatinv)   
#       Cinv<-solve(Sb)
       Sb2 <- Cinv%*%t(Z)
       DD <- t(Z)%*%W%*%y
       S <- Z%*%Sb2%*%W
       yp <- S%*%y
       e <- y - yp                   
       b.est=Sb2%*%y
       beta.est=fd(b.est[-1,1],basis.b)
       df=basis.b$nbasis+1     
       rdf <- n-df
       sr2 <- sum(e^2)/ rdf
       r2 <- 1 - sum(e^2)/sum(ycen^2)
       Vp <- sr2*Cinv       
#       bet<-Cinv%*%DD
       a.est=b.est[1,1]
       #beta.est2=fd(b.est2[-1,1]*diff(rtt),basis.b)

       r2.adj<- 1 - (1 - r2) * ((n -    1)/ rdf)
       GCV <- sum(e^2)/(n - df)^2
       object.lm=list()
       object.lm$coefficients<-drop(b.est)
       object.lm$residuals<-drop(e)
       object.lm$fitted.values<-yp
       object.lm$y<-y
       object.lm$rank<-df
       object.lm$df.residual<-n-df
       colnames(Z)[1]="(Intercept)"
       vcov2=sr2*Cinv
       std.error=sqrt(diag(vcov2))
       t.value=b.est/std.error
       p.value= 2 * pt(abs(t.value),n-df, lower.tail = FALSE)
       coefficients<-cbind(b.est,std.error,t.value,p.value)
       colnames(coefficients) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
       rownames(coefficients)[1]="(Intercept)"
       class(object.lm) <- "lm"
	   b.est <- b.est[-1]
	   names(b.est) <- rownames(coefficients)[-1]
     }
    #  GCV <- sum(e^2)/(n - df)^2      #GCV"=GCV,

#hat<-diag(hat(Z, intercept = TRUE),ncol=n)
out<-list("call"=call, coefficients=coefficients, "residuals"=e, "fitted.values"=yp
,"beta.est"=beta.est,weights= weights,"df.residual"=n-df,"r2"=r2,"sr2"=sr2,
"Vp"=Vp,"H"=S,"y"=y,"fdataobj"=fdataobj, "basis.x.opt"=basis.x, #x.fd=x.fd,
"basis.b.opt"=basis.b,"J"=J,"lambda.opt"=lambda,P=R, Lfdobj=Lfdobj,
  lm=object.lm,"mean"=xmean, "b.est"=b.est,"a.est"=a.est,XX=Z)
class(out) <- "fregre.fd"
return(out)
}
