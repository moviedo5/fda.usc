#' @title Partial least squares components for functional data.
#' 
#' @description Compute penalized partial least squares (PLS) components for functional
#' data. 
#' 
#' @details If \code{norm=TRUE}, computes the PLS by
#' \code{NIPALS} algorithm and the Degrees of Freedom using the Krylov
#' representation of PLS, see Kraemer and Sugiyama (2011).\cr 
#' If \code{norm=FALSE}, computes the PLS by Orthogonal Scores Algorithm and
#' the Degrees of Freedom are the number of components \code{ncomp}, see
#' Martens and Naes (1989).
#' 
#' @aliases fdata2pls 
#' @param fdataobj \code{\link{fdata}} class object.
#' @param y Scalar response with length \code{n}.
#' @param ncomp The number of components to include in the model.
#' @param lambda Amount of penalization. Default value is 0, i.e. no
#' penalization is used.
#' @param P If P is a vector: coefficients to define the penalty matrix object.
#' By default \eqn{P=c(0,0,1)} penalizes the second derivative (curvature) or
#' acceleration.  If P is a matrix: the penalty matrix object.
#' @param norm If \code{TRUE} the \code{fdataobj} are centered and scaled.
#' @param \dots Further arguments passed to or from other methods.
#' @return \code{fdata2pls} function return: 
#' \itemize{
#' \item {df}{ degree of freedom}
#' \item {rotation}{ \code{\link{fdata}} class object.} 
#' \item {x}{ Is true the value of the rotated data (the centred data multiplied by the rotation matrix) is returned.}
#' \item {fdataobj.cen}{ The centered \code{fdataobj} object.} 
#' \item {mean}{ mean of \code{fdataobj}.} 
#' \item {l}{Vector of index of principal components.} 
#' \item {C}{ The matched call.}
#' \item {lambda}{ Amount of penalization.} 
#' \item {P}{ Penalty matrix.}
#' }
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente  \email{manuel.oviedo@@usc.es}
#' @seealso  Used in:
#' \code{\link{fregre.pls}}, \code{\link{fregre.pls.cv}}.
#' Alternative method: \code{\link{fdata2pc}}.
#' @references
#'  Kraemer, N., Sugiyama M. (2011). \emph{The Degrees of Freedom of
#' Partial Least Squares Regression}. Journal of the American Statistical
#' Association. Volume 106, 697-705.
#' 
#' Febrero-Bande, M., Oviedo de la Fuente, M. (2012).  \emph{Statistical
#' Computing in Functional Data Analysis: The R Package fda.usc.} Journal of
#' Statistical Software, 51(4), 1-28. \url{http://www.jstatsoft.org/v51/i04/}
#' 
#' Martens, H., Naes, T. (1989) \emph{Multivariate calibration.} Chichester:
#' Wiley.
#' @keywords multivariate
#' @examples 
#' \dontrun{
#' n= 500;tt= seq(0,1,len=101)
#' x0<-rproc2fdata(n,tt,sigma="wiener")
#' x1<-rproc2fdata(n,tt,sigma=0.1)
#' x<-x0*3+x1
#' beta = tt*sin(2*pi*tt)^2
#' fbeta = fdata(beta,tt)
#' y<-inprod.fdata(x,fbeta)+rnorm(n,sd=0.1)
#' pls1=fdata2pls(x,y)
#' norm.fdata(pls1$rotation)
#' }
#' @export
fdata2pls<-function(fdataobj,y,ncomp = 2,lambda=0,P=c(0,0,1),norm=TRUE,...) {
  if (!is.fdata(fdataobj)) fdataobj<-fdataobj(fdataobj)
  C <- match.call()
  X<-fdataobj$data
  tt<-fdataobj[["argvals"]]
  rtt<-fdataobj[["rangeval"]]
  nam<-fdataobj[["names"]]
  J <- ncol(X);n <- nrow(X)
  Jmin<-min(c(J,n))
  Jmax = min(J+1,n-1)
  Beta <- matrix(, J,ncomp)
  W <- V <- Beta
  dW <- dBeta <- dV <- array(dim = c(ncomp, J, n))
  X0 <- X;  y0 <- y
  mean.y <- mean(y)
  y <- scale(y, scale = FALSE)
  center<-fdata.cen(fdataobj)
  mean.X<-center$meanX
  repnX<-rep(1, n)
  X<-center$Xcen$data
  
  if (norm)    {
    sd.X <- sqrt(apply(X, 2, var))
    X <- X/(rep(1, nrow(X)) %*% t(sd.X))
    X2<-fdata(X,tt,rtt,nam)
    dcoefficients = NULL
    tX<-t(X)    
    A <- crossprod(X)
    b <- crossprod(X,y)
    if (lambda>0) {
      if (is.vector(P))  {     P<-P.penalty(tt,P)          }
      #     else {
      dimp<-dim(P)
      if (!(dimp[1]==dimp[2] & dimp[1]==J))
        stop("Incorrect matrix dimension P")
      #         }
      M <- solve( diag(J) + lambda*P)
      W[, 1]<- M %*%b
    }
    else    {
      M<-NULL
      W[, 1] <- b
    }
    dV[1, , ] <- dW[1, , ] <- dA(W[, 1], A, tX)
    W[, 1] <- W[, 1]/sqrt((sum((W[, 1]) * (A %*% W[,1]))))
    V[, 1] <- W[, 1]
    Beta[, 1] <- sum(V[, 1] * b) * V[, 1]
    dBeta[1, , ] <-dvvtz(V[, 1], b, dV[1, , ],tX)
    if (ncomp>1) {
      for (i in 2:ncomp) {
        vsi <-b - A %*% Beta[, i - 1]
        if (!is.null(M))  vsi<- M %*% vsi
        W[, i]<-vsi
        dW[i, , ] <- t(X) - A %*% dBeta[i - 1, , ]
        V[, i] <- W[, i] - vvtz(V[, 1:(i - 1), drop = FALSE],A %*% W[, i])
        dV[i, , ] = dW[i, , ] - dvvtz(V[, 1:(i - 1),drop = FALSE],
                                      A %*% W[, i], dV[1:(i - 1), , , drop = FALSE], A %*% dW[i, , ])
        dV[i, , ] <- dA(V[, i], A, dV[i, , ])
        V[, i] <- V[, i]/sqrt((sum(t(V[, i]) %*% A %*% V[,i])))
        Beta[, i] = Beta[, i - 1] + sum(V[, i] * b) * V[,i]
        dBeta[i,,]<-dBeta[i-1,,]+dvvtz(V[,i],b,dV[i,,],tX)
      }
    }
    dcoefficients <- NULL
    dcoefficients <- array(0, dim = c(ncomp + 1, J, n))
    dcoefficients[2:(ncomp + 1), , ] = dBeta
    sigmahat <- RSS <- yhat <- vector(length = ncomp + 1)
    DoF <- 1:(ncomp + 1)
    Yhat <- matrix(, n, ncomp + 1)
    dYhat <- array(dim = c(ncomp + 1, n, n))
    coefficients <- matrix(0, J, ncomp + 1)
    coefficients[, 2:(ncomp + 1)] = Beta/(sd.X %*% t(rep(1, ncomp)))
    intercept <- rep(mean.y, ncomp + 1) - t(coefficients) %*% t(mean.X$data)
    covariance <- array(0, dim = c(ncomp + 1, J, J))
    DD <- diag(1/sd.X)
    for (i in 1:(ncomp + 1)) {
      Yhat[, i] <- X0 %*% coefficients[, i] + intercept[i]
      res <- y0 - Yhat[, i]
      yhat[i] <- sum((Yhat[, i])^2)
      RSS[i] <- sum(res^2)
      dYhat[i, , ]<-X %*%dcoefficients[i, , ] + matrix(1,n, n)/n
      DoF[i] <- sum(diag(dYhat[i, , ]))
      #        dummy <- (diag(n) - dYhat[i, , ]) %*% (diag(n) - t(dYhat[i, , ]))
      #        sigmahat[i] <- sqrt(RSS[i]/sum(diag(dummy)))
      #        if (i > 1) {
      #                covariance[i, , ] <- sigmahat[i]^2 * DD %*% dcoefficients[i,
      #                  , ] %*% t(dcoefficients[i, , ]) %*% DD
      #            }
    }
    V2<- fdata(t(V)*(rep(1, nrow(t(V))) %*% t(sd.X)),tt,rtt,nam)
    V2$data<-sweep(V2$data,1,norm.fdata(V2),"/")
    #   V2<-fdata(t(V),tt,rtt,nam)
    #   V2$data<-sweep(V2$data,1,norm.fdata(V2),"/")
    #    W2<-fdata(t(W),tt,rtt,nam)
    #    X3<-fdata(X,tt,rtt,nam)
    #    beta.est<-fdata(t(Beta),tt,rtt,nam)
    DoF[DoF > Jmax] = Jmax
    #    intercept <- as.vector(intercept)    
  }
  else {
    plsr<-mplsr(X,y,ncomp=ncomp,lambda=lambda,P=P,...)   
    V2<-fdata(t(plsr$loading.weights),tt,rtt,nam)
    X2<-fdata(X,tt,rtt,nam)
    DoF<-1:ncomp
    Yhat<-plsr$fitted.values    
    yhat<-sum(Yhat)^2          
  }    
  scores<-inprod.fdata(X2,V2,...)
  #   TT = X %*% V  ## son los scores
  l<-1:ncomp
  colnames(scores) <- paste("PLS", l, sep = "")
  outlist = list(call=C,df = DoF, rotation=V2, x=scores, lambda=lambda,P=P,norm=norm,type="pls",
                 fdataobj=fdataobj,y=y0,l=l,fdataobj.cen=center$Xcen,mean=mean.X,Yhat = Yhat,yhat = yhat)
  #    Yhat = Yhat, coyhat = yhat,efficients = coefficients,intercept = intercept,
  #     RSS = RSS,TT=TT, sigmahat = sigmahat,covariance = covariance,
  #     W2=W2,beta.est=beta.est,
  class(outlist) <- "fdata.comp"
  return(outlist)
}



#' Principal components for functional data
#' 
#' @description Compute (penalized) principal components for functional data. 
#' 
#' @details Smoothing is achieved by penalizing the integral of the square of the
#' derivative of order m over rangeval: \itemize{ \item m = 0 penalizes the
#' squared difference from 0 of the function \item m = 1 penalize the square of
#' the slope or velocity \item m = 2 penalize the squared acceleration \item m
#' = 3 penalize the squared rate of change of acceleration }
#' 
#' @aliases fdata2pc 
#' @param fdataobj \code{\link{fdata}} class object.
#' @param ncomp Number of principal comoponents.
#' @param norm =TRUE the norm of eigenvectors \code{(rotation)} is 1.
#' @param lambda Amount of penalization. Default value is 0, i.e. no
#' penalization is used.
#' @param P If P is a vector: coefficients to define the penalty matrix object.
#' By default P=c(0,0,1) penalize the second derivative (curvature) or
#' acceleration.  If P is a matrix: the penalty matrix object.
#' @param \dots Further arguments passed to or from other methods.
#' @return
#' \itemize{
#' \item {d}{ The standard deviations of the functional principal components.} 
#' \item {rotation}{ are also known as loadings.  A \code{fdata} class object whose rows contain the eigenvectors.} 
#' \item {x}{ are also known as scores. The value of the rotated functional data is returned.}
#' \item {fdataobj.cen}{ The centered \code{fdataobj} object.} 
#' \item {mean}{ The functional mean of \code{fdataobj} object.} 
#' \item {l}{ Vector of index of principal components.} 
#' \item {C}{ The matched call.} 
#' \item {lambda}{ Amount of penalization.} 
#' \item {P}{ Penalty matrix.}
#' }
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@usc.es}
#' @seealso See Also as \link[base]{svd} and \link[stats]{varimax}.
#' @references Venables, W. N. and B. D. Ripley (2002). \emph{Modern Applied
#' Statistics with S}. Springer-Verlag.
#' 
#' N. Kraemer, A.-L. Boulsteix, and G. Tutz (2008). Penalized Partial Least
#' Squares with Applications to B-Spline Transformations and Functional Data.
#' Chemometrics and Intelligent Laboratory Systems, 94, 60 - 69.
#' \url{http://dx.doi.org/10.1016/j.chemolab.2008.06.009}
#' 
#' Febrero-Bande, M., Oviedo de la Fuente, M. (2012).  \emph{Statistical
#' Computing in Functional Data Analysis: The R Package fda.usc.} Journal of
#' Statistical Software, 51(4), 1-28. \url{http://www.jstatsoft.org/v51/i04/}
#' @keywords multivariate
#' @examples
#'  \dontrun{
#'  n= 100;tt= seq(0,1,len=51)
#'  x0<-rproc2fdata(n,tt,sigma="wiener")
#'  x1<-rproc2fdata(n,tt,sigma=0.1)
#'  x<-x0*3+x1
#'  pc=fdata2pc(x,lambda=1)
#'  summary(pc)
#'  }
#' @export
fdata2pc<-function (fdataobj,  ncomp = 2,norm = TRUE,lambda=0,P=c(0,0,1),...)
{
  C <- match.call()
  if (!is.fdata(fdataobj))
    stop("No fdata class")
  nas1 <- is.na.fdata(fdataobj)
  if (any(nas1))
    stop("fdataobj contain ", sum(nas1), " curves with some NA value \n")
  X <- fdataobj[["data"]]
  tt <- fdataobj[["argvals"]]
  rtt <- fdataobj[["rangeval"]]
  nam <- fdataobj[["names"]]
  mm <- fdata.cen(fdataobj)
  xmean <- mm$meanX
  Xcen.fdata <- mm$Xcen
  dimx<-dim(X)
  n <- dimx[1]
  J <- dimx[2]
  Jmin <- min(c(J, n))
  if (lambda>0) {
    if (is.vector(P))  {     P<-P.penalty(tt,P)          }
    #     else {
    dimp<-dim(P)
    if (!(dimp[1]==dimp[2] & dimp[1]==J))
      stop("Incorrect matrix dimension P")
    #          }
    M <- solve( diag(J) + lambda*P)
    Xcen.fdata$data<- Xcen.fdata$data %*%t(M)
  }
  eigenres <- svd(Xcen.fdata$data)
  v <- eigenres$v
  u <- eigenres$u
  d <- eigenres$d
  D <- diag(d)
  vs <- fdata(t(v), tt, rtt, list(main = "fdata2pc", xlab = "t",
                                  ylab = "rotation"))
  scores <- matrix(0, ncol = J, nrow = n)
  if (norm) {
    dtt <- diff(tt)
    drtt <- diff(rtt)
    eps <- as.double(.Machine[[1]] * 10)
    inf <- dtt - eps
    sup <- dtt + eps
    if (all(dtt > inf) & all(dtt < sup))
      delta <- sqrt(drtt/(J - 1))
    else delta <- 1/sqrt(mean(1/dtt))
    no <- norm.fdata(vs)
    vs <- vs/delta
    newd <- d * delta
    scores[, 1:Jmin] <- inprod.fdata(Xcen.fdata, vs, ...)
  }
  else {
    scores[, 1:Jmin] <- inprod.fdata(Xcen.fdata, vs, ...)
    newd <- d
  }
  colnames(scores) <- paste("PC", 1:J, sep = "")
  l <- 1:ncomp
  out <- list(call = C,d = newd, rotation = vs[1:ncomp],x = scores,
              lambda = lambda,P=P, fdataobj.cen = Xcen.fdata,norm=norm,type="pc",
              mean = xmean, fdataobj = fdataobj,l=l,u=u[,1:ncomp,drop=FALSE])
  class(out) = "fdata.comp"
  return(out)
}


D.penalty=function(tt) {
  rtt<-diff(tt)
  hh<--(1/mean(1/rtt))/rtt
  Dmat1<-diag(hh)
  p<-length(tt)
  mat<-matrix(rep(0,p-1),ncol=1)
  Dmat12<-cbind(Dmat1,mat)
  Dmat22<-cbind(mat,-Dmat1)
  Dmat<-Dmat12+Dmat22
  return(Dmat)
}

################################################################
################################################################
# P.penalty <-function(tt,order=2) {
#   D.pen2<-D.pen <- D.penalty(tt)
#  p<-ncol(D.pen)
#   if (order > 1) {
#      for (k in 2:order) D.pen <- D.pen2[1:(p-k),1:(p-k+1)] %*% D.pen
#   }
#   D.pen <- t(D.pen) %*% D.pen
# }
################################################################
################################################################


#' Penalty matrix for higher order differences
#' 
#' This function computes the matrix that penalizes the higher order
#' differences.
#' 
#' For example, if \code{P}=c(0,1,2), the function return the penalty matrix
#' the second order difference of a vector \eqn{tt}. That is \deqn{v^T P_j tt=
#' \sum_{i=3} ^{n} (\Delta tt_i) ^2} where \deqn{\Delta tt_i= tt_i -2 tt_{i-1}
#' + tt_{i-2}} is the second order difference. More details can be found in
#' Kraemer, Boulesteix, and Tutz (2008).
#' 
#' @param tt vector of the \code{n} discretization points or argvals.
#' @param P vector of coefficients with the order of the differences. Default
#' value \code{P}=c(0,0,1) penalizes the second order difference.
#' @return penalty matrix of size \code{sum(n)} x \code{sum(n)}
#' @note The discretization points can be equidistant or not.
#' @author This version is created by Manuel Oviedo de la Fuente modified the
#' original version created by Nicole Kramer in \code{ppls} package.
#' @seealso \code{\link{fdata2pls}}
#' @references N. Kraemer, A.-L. Boulsteix, and G. Tutz (2008). \emph{Penalized
#' Partial Least Squares with Applications to B-Spline Transformations and
#' Functional Data}. Chemometrics and Intelligent Laboratory Systems, 94, 60 -
#' 69. \url{http://dx.doi.org/10.1016/j.chemolab.2008.06.009}
#' @keywords math
#' @examples
#' 
#' P.penalty((1:10)/10,P=c(0,0,1))
#' # a more detailed example can be found under script file 
#' @export P.penalty
P.penalty <- function(tt,P=c(0,0,1)) {
  lenp<-length(P)
  D.pen2<-D.pen <- D.penalty(tt)
  p<-ncol(D.pen)
  D.pen3<-matrix(0,p,p)
  if  (sum(P)!=0){
    D.pen3<-D.pen3+P[1]*diag(p)
    if (lenp > 1) {
      P<-P[-1]
      for (i in 2:lenp-1){
        D.pen2<-D.pen <- D.penalty(tt)
        if (i > 1) {
          for (k in 2:i) D.pen <- D.pen2[1:(p-k),1:(p-k+1)] %*% D.pen
        }
        D.pen3 <- D.pen3+P[i]*(t(D.pen) %*% D.pen)
      }
    }}
  D.pen3
}
################
################
dA<-function (w, A, dw){
  wa <- sqrt(sum((w * (A %*% w))))
  dummy <- (1/wa) * (diag(length(w)) - w %*% t(w) %*% A/(wa^2)) %*%   dw
  return(dummy)
}
################
################
vvtz<-function (v, z){as.vector(v %*% (t(v) %*% z))}
################
################
dvvtz<-function (v, z, dv, dz) {
  if (is.matrix(v) == FALSE) {
    v <- matrix(v, ncol = 1)
    dv <- array(dv, dim = c(1, nrow(dv), ncol(dv)))
  }
  k = ncol(v)
  p <- nrow(v)
  n <- dim(dv)[3]
  dummy <- matrix(0, dim(dv)[2], dim(dv)[3])
  for (i in 1:k) {
    D <- (v[, i] %*% t(z) + sum(v[, i] * z) * diag(p)) %*%
      dv[i, , ] + v[, i] %*% t(v[, i]) %*% dz
    dummy <- dummy + D
  }
  return(dummy)
}