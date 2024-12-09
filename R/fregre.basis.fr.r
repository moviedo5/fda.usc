#' Functional Regression with functional response using basis representation.
#' 
#' @description Computes functional regression between functional explanatory variable
#' \eqn{X(s)} and functional response \eqn{Y(t)} using basis representation.
#' 
#' @details 
#' \deqn{Y(t)=\alpha(t)+\int_{T}{X(s)\beta(s,t)ds+\epsilon(t)}}{Y=\alpha(t)+\int
#' X(s)\beta(s,t)ds+\epsilon(t)}
#' 
#' where \eqn{\alpha(t)} is the intercept function, \eqn{\beta(s,t)} is the
#' bivariate resgression function and \eqn{\epsilon(t)} are the error term with
#' mean zero.
#'  
#' The function is a wrapped of \link[fda]{linmod} function proposed by
#' Ramsay and Silverman (2005) to model the relationship between the functional
#' response \eqn{Y(t)} and the functional covariate \eqn{X(t)} by basis
#' representation of both.
#' 
#' The unknown bivariate functional parameter \eqn{\beta(s,t)}{\beta(s,t)} can
#' be expressed as a double expansion in terms of \eqn{K} basis function
#' \eqn{\nu_k} and \eqn{L} basis functions \eqn{\theta_l},
#' \deqn{\beta(s,t)=\sum_{k=1}^{K}\sum_{l=1}^{L} b_{kl}
#' \nu_{k}(s)\theta_{l}(t)=\nu(s)^{\top}\bold{B}\theta(t)}{\beta(s,t)=\sum_k
#' \sum_l b_{kl} \nu_k(s)\theta_l(t)=\nu(s)'B\theta(t)} Then, the model can be
#' re--written in a matrix version as,
#' \deqn{Y(t)=\alpha(t)+\int_{T}{X(s)\nu(s)^{\top}\bold{B}\theta(t)ds+\epsilon(t)}=\alpha(t)+\bold{XB}\theta(t)+\epsilon(t)}{Y(t)=\alpha(t)+\int
#' X(s)\nu(s)'B\theta(t)ds+\epsilon(t)=\alpha(t)+XB\theta(t)+\epsilon(t)} where
#' \eqn{\bold{X}=\int X(s)\nu^{\top}(t)ds}{X=\int X(s)\nu'(t)ds} \cr
#' 
#' This function allows objects of class \code{fdata} or directly covariates of
#' class \code{fd}.  If \code{x} is a \code{fdata} class, \code{basis.s} is
#' also the basis used to represent \code{x} as \code{fd} class object. If
#' \code{y} is a \code{fdata} class, \code{basis.t} is also the basis used to
#' represent \code{y} as \code{fd} class object. The function also gives
#' default values to arguments \code{basis.s} and \code{basis.t} for construct
#' the bifd class object used in the estimation of \eqn{\beta(s,t)}.  If
#' \code{basis.s=}\code{NULL} or \code{basis.t=}\code{NULL} the function
#' creates a \code{bspline} basis by \link[fda]{create.bspline.basis}.
#' 
#' \code{fregre.basis.fr} incorporates a roughness penalty using an appropiate
#' linear differential operator; \code{lambda.s}, \code{Lfdobj.s} for
#' penalization of \eqn{\beta}{\beta}'s variations with respect to \eqn{s} and
#' \cr \code{lambda.t}, \code{Lfdobj.t} for penalization of
#' \eqn{\beta}{\beta}'s variations with respect to \eqn{t}.\cr
#' 
#' @param x Functional explanatory variable.
#' @param y Functional response variable.
#' @param basis.s Basis related with \code{s} and it is used in the estimation
#' of \eqn{\beta(s,t)}.
#' @param basis.t Basis related with \code{t} and it is used in the estimation
#' of \eqn{\beta(s,t)}.
#' @param lambda.s A roughness penalty with respect to \code{s} to be applied
#' in the estimation of \eqn{\beta(s,t)}. By default, no penalty
#' \code{lambda.s=0}.
#' @param lambda.t A roughness penalty with respect to \code{t} to be applied
#' in the estimation of \eqn{\beta(s,t)}.  By default, no penalty
#' \code{lambda.t=0}.
#' @param Lfdobj.s A linear differential operator object with respect to
#' \code{s} . See \link[fda]{eval.penalty}.
#' @param Lfdobj.t A linear differential operator object with respect to
#' \code{t}. See \link[fda]{eval.penalty}.
#' @param weights Weights.
#' @param \dots Further arguments passed to or from other methods.
#' @return Return:
#' \itemize{
#' \item \code{call}: The matched call. 
#' \item \code{a.est}: Intercept parameter estimated. 
#' \item \code{coefficients}: The matrix of the coefficients. 
#' \item \code{beta.est}: A bivariate functional data object of class \code{bifd} with the estimated parameters of \eqn{\beta(s,t)}{\beta(s,t)}. 
#' \item \code{fitted.values}: Estimated response. 
#' \item \code{residuals}: \code{y} minus \code{fitted values}. 
#' \item \code{y}: Functional response. 
#' \item \code{x}: Functional explanatory data. 
#' \item \code{lambda.s}: A roughness penalty with respect to \code{s}. 
#' \item \code{lambda.t}: A roughness penalty with respect to \code{t}. 
#' \item \code{Lfdobj.s}: A linear differential operator with respect to \code{s}. 
#' \item \code{Lfdobj.t}: A linear differential operator with respect to \code{t}. 
#' \item \code{weights}: Weights. 
#' }
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' @seealso See Also as: \code{\link{predict.fregre.fr}}.
#' Alternative method: \link[fda]{linmod}.
#' @references Ramsay, James O., and Silverman, Bernard W. (2006), \emph{
#' Functional Data Analysis}, 2nd ed., Springer, New York.
#' @keywords regression
#' @examples
#' \dontrun{
#' rtt<-c(0, 365)
#' basis.alpha  <- create.constant.basis(rtt)
#' basisx  <- create.bspline.basis(rtt,11)
#' basisy  <- create.bspline.basis(rtt,11)
#' basiss  <- create.bspline.basis(rtt,7)
#' basist  <- create.bspline.basis(rtt,9)
#' 
#' # fd class
#' dayfd<-Data2fd(day.5,CanadianWeather$dailyAv,basisx)
#' tempfd<-dayfd[,1]
#' log10precfd<-dayfd[,3]
#' res1 <-  fregre.basis.fr(tempfd, log10precfd,
#' basis.s=basiss,basis.t=basist)
#' 
#' # fdata class
#' tt<-1:365
#' tempfdata<-fdata(t(CanadianWeather$dailyAv[,,1]),tt,rtt)
#' log10precfdata<-fdata(t(CanadianWeather$dailyAv[,,3]),tt,rtt)
#' res2<-fregre.basis.fr(tempfdata,log10precfdata,
#' basis.s=basiss,basis.t=basist)
#' 
#' # penalization
#' Lfdobjt <- Lfdobjs <- vec2Lfd(c(0,0), rtt)
#' Lfdobjt <- vec2Lfd(c(0,0), rtt)
#' lambdat<-lambdas <- 100
#' res1.pen <- fregre.basis.fr(tempfdata,log10precfdata,basis.s=basiss,
#' basis.t=basist,lambda.s=lambdas,lambda.t=lambdat,
#' Lfdobj.s=Lfdobjs,Lfdobj.t=Lfdobjt)
#' 
#' res2.pen <- fregre.basis.fr(tempfd, log10precfd,
#' basis.s=basiss,basis.t=basist,lambda.s=lambdas,
#' lambda.t=lambdat,Lfdobj.s=Lfdobjs,Lfdobj.t=Lfdobjt)
#' 
#' plot(log10precfd,col=1)
#' lines(res1$fitted.values,col=2)
#' plot(res1$residuals)
#' plot(res1$beta.est,tt,tt)
#' plot(res1$beta.est,tt,tt,type="persp",theta=45,phi=30)
#' }
#' @aliases fregre.basis.fr
#' @export
fregre.basis.fr <- function(x,y,basis.s=NULL,basis.t=NULL,
lambda.s=0,lambda.t=0,Lfdobj.s=vec2Lfd(c(0,0),range.s),
Lfdobj.t=vec2Lfd(c(0,0),range.t),weights=NULL,...){
call<-match.call()
isfdx<-is.fd(x)
x.orig<-x
isfdy<-is.fd(y)
y.orig<-y
if (isfdx) {
  xfdobj<-x
  basis.x<- x$basis
  nbasis.x<- basis.x$nbasis
  range.x<- basis.x$rangeval
  if (is.null(basis.s))  basis.s<-basis.x
}
else {
  if (!is.fdata(x))  stop("x is not a functional data object of class fd or fdata.")
    range.x<-x$rangeval
    if (is.null(basis.s))  {
        np<-ncol(x)
        #nbasis.b = min(51,floor(np/10))
        nbasis.b =min(7,floor(np/2))
        basis.s<-basis.x<-create.bspline.basis(rangeval=range.x,nbasis=nbasis.b)
   }
   else   basis.x<-basis.s
   xfdobj<-Data2fd(argvals =x$arg, y = t(x$data), basisobj = basis.s,...)
}
if (isfdy) {
  yfdobj<-y
  basis.y<-y$basis
  nbasis.y<-basis.y$nbasis
  range.y<-basis.y$rangeval
  np.y<-nbasis.y
  nfine = min(c(201, 10 * np.y+ 1))
  tty <- seq(range.y[1],range.y[2],len=nfine)
#tfine = seq(ranget[1], ranget[2], len = nfine)

  if (is.null(basis.t))  basis.t<-basis.y
}
else {
  if (!is.fdata(y))  stop("y is not a functional data object of class fd or fdata.")
  range.y<-y$rangeval
  np.y<-ncol(y)
  nbasis.y =min(7,floor(np.y/2))
#  tty<-y$argvals
  if (is.null(basis.t))  {
        #nbasis.y = min(51,floor(np/10))
        basis.y<-basis.t<-create.bspline.basis(rangeval=range.y,nbasis=nbasis.y)
   }
   else   basis.y<-basis.t
  nfine = min(c(201, 10 * np.y+ 1))
  tty <- seq(range.y[1],range.y[2],len=nfine)
   yfdobj<-Data2fd(argvals =y$argvals, y = t(y$data), basisobj = basis.t,...)
}
#xLfdobj<-vec2Lfd(c(0,0,0), c(range.x[1],range.x[2]))
coefy   = yfdobj$coef
coefx   = xfdobj$coef
coefdx  = dim(coefx)
coefdy  = dim(coefy)
n=ncurves = coefdx[2]
################################################################################
#alphabasis<-basis.t#create.constant.basis(range.y)
range.t = basis.t$rangeval
range.s = basis.s$rangeval
###fdPar(smallbasis3, xLfdobj, xlambda)$Lfd
alphafdPar = fdPar(basis.s,Lfdobj.s, lambda.s)
alphalambda = alphafdPar$lambda
alphabasis<-alphafdPar$fd$basis
alphanbasis = alphabasis$nbasis
Finprod = inprod(basis.y, alphabasis)
alphattmat = diff(range.y) #cambiar
alphattmat = inprod(alphabasis, alphabasis)
if (range.s[1] != range.x[1] || range.s[2] != range.x[2]) {
    stop("Range of beta.s and x not compatible.")
}
nbasis.s = basis.s$nbasis
Hinprod = inprod(basis.x, basis.s)
xcoef = xfdobj$coef
basis.ss  = inprod(basis.s, basis.s)
Q<-S<-R<- NULL
range.t = basis.t$rangeval
if (range.t[1] != range.y[1] || range.t[2] != range.y[2]) {
    stop("Range of BETATFD coefficient and YFD not compatible.")
}
nbasis.t = basis.t$nbasis
Ginprod = inprod(basis.y, basis.t)
ycoef = yfdobj$coef
basis.tt  = inprod(basis.t, basis.t)
basis.talpha= inprod(basis.t, alphabasis)
Fmat = t(ycoef) %*% Finprod
Gmat = t(ycoef) %*% Ginprod
Hmat = t(xcoef) %*% Hinprod
if (is.null(weights)) {
    HHCP = t(Hmat) %*% Hmat
    HGCP = t(Hmat) %*% Gmat
    H1CP = as.matrix(colSums(Hmat))
    F1CP = as.matrix(colSums(Fmat))
} else {
    HHCP = t(Hmat) %*% (outer(weights,rep(nbasis.s))*Hmat)
    HGCP = t(Hmat) %*% (outer(weights,rep(nbasis.t))*Gmat)
    H1CP = t(Hmat) %*% weights
    F1CP = t(Fmat) %*% weights
}
#############################
alphattmat = inprod(alphabasis, alphabasis)
betalttmat = inprod(basis.t, alphabasis)
betassmat = inprod(basis.s, basis.s)
betattmat = inprod(basis.t, basis.t)

Beta1fd  = bifd(matrix(0,nbasis.s,nbasis.t), basis.s, basis.t)#
Beta1Par = bifdPar(Beta1fd, Lfdobj.s, Lfdobj.t, lambda.s, lambda.t)
betaslambda<-Beta1Par$lambdas
betatlambda<-Beta1Par$lambdat

if (alphalambda > 0) {   Q = eval.penalty(alphabasis, alphafdPar$Lfd) }
if (lambda.s > 0) {   R = eval.penalty(basis.s,Beta1Par$Lfds)}
if (lambda.t > 0) {  S = eval.penalty(basis.t,Beta1Par$Lfdt)}

betan = nbasis.s*nbasis.t
ncoef = alphanbasis+betan
Cmat = matrix(0, ncoef, ncoef)
Dmat = matrix(0, ncoef, 1)
ind1 = 1:alphanbasis
ind2 = ind1
Cmat[ind1, ind2] = ncurves * alphattmat
if (alphalambda > 0) {    Cmat[ind1, ind2] = Cmat[ind1, ind2] + alphalambda * Q  }
ind2 = alphanbasis + (1:betan)
Cmat[ind1,ind2] = t(kronecker(H1CP,basis.talpha))
Dmat[ind1] = F1CP
ind1 = alphanbasis + (1:betan)
ind2 = 1:alphanbasis
Cmat[ind1, ind2] = t(Cmat[ind2, ind1])
ind2 = ind1
Cmat[ind1, ind2] = kronecker(HHCP, betattmat)

if (betaslambda> 0) {
  Cmat[ind1,ind2] = Cmat[ind1,ind2] + lambda.s*kronecker(R,basis.tt)
}
if (betatlambda > 0) {
    Cmat[ind1,ind2] = Cmat[ind1,ind2] + lambda.t*kronecker(basis.ss,S)
}
Dmat[ind1] = matrix(t(HGCP), betan, 1)
coefvec = symsolve(Cmat, Dmat)
ind1 = 1:alphanbasis
alpha.est = coefvec[ind1]
alphafdnames = yfdobj$fdnames
alphafdnames[[3]] = "Intercept"
alphafd = fd(alpha.est, alphabasis, alphafdnames)
ind1 = alphanbasis + (1:betan)
beta.est = matrix(coefvec[ind1],nbasis.t,nbasis.s)
betafdnames = xfdobj$fdnames
betafdnames[[3]] = "Reg. Coefficient"
betafd = bifd(t(beta.est), basis.s, basis.t, betafdnames)
beta.xest = beta.est %*% t(Hmat)
xbetafd = fd(beta.xest, basis.t)
if (isfdy) {
   yhatmat = eval.fd(tty, alphafd) %*% matrix(1, 1, ncurves)+   eval.fd(tty,  xbetafd)
   yhatfd = smooth.basis(tty, yhatmat, basis.y)$fd
   fitted.values  <-yhatfd # fd(coef=yhatmat, basisobj=basis.y, fdnames=yfdobj$fdnames)
 }
else {
   yhatmat = eval.fd(y$argvals, alphafd) %*% matrix(1, 1, ncurves)+   eval.fd(y$argvals,  xbetafd)
   fitted.values<-fdata(t(yhatmat),y$argvals,y$rangeval,y$names)
   #betafd = fdata(eval.bifd(t(beta.est), basis.s, basis.t, betafdnames)
   #      betafd = fdata(xbetafd,fdata2d)
   #plot(fdata(eval.bifd(argvals.s,argvals.t,x),tt,rtt,fdata2d=TRUE),... )
}

residuals<-y-fitted.values
out = list(call = call, alpha.est = alphafd, coefficients = beta.est, 
             beta.estbifd = betafd, fitted.values = fitted.values, 
             residuals = residuals, lambda.s = lambda.s, lambda.t = lambda.t, 
             Lfdobj.s = Lfdobj.s, Lfdobj.t = Lfdobj.t, weights = weights, 
             x = x, y = y, H = Hinprod, basis.s = basis.s, basis.t = basis.t, 
             argvals.y = tty
            , betaspenmat= R, betatpenmat= S,alphapenmat=Q)
class(out)<-"fregre.fr"
return(out)
}
