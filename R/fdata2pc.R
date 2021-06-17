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
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente  \email{manuel.oviedo@@udc.es}
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
#' Statistical Software, 51(4), 1-28. \url{https://www.jstatsoft.org/v51/i04/}
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
#' pls1$call
#' summary(pls1)
#' pls1$l
#' norm.fdata(pls1$rotation)
#' }
#' @export
fdata2pls<-function(fdataobj, y, ncomp = 2,lambda = 0,
                    P=c(0,0,1), norm = TRUE,...) {
  if (!is.fdata(fdataobj)) fdataobj<-fdataobj(fdataobj)
  C <- match.call()
  # X <- fdataobj$data
  # tt <- fdataobj[["argvals"]]
  # rtt<- fdataobj[["rangeval"]]
  # nam <- fdataobj[["names"]]
  # J <- ncol(X);n <- nrow(X)
  # Jmin <- min(c(J,n))
  # Jmax <- min(J+1,n-1)
  # Beta <- matrix(, J,ncomp)
  # W <- V <- Beta
  # dW <- dBeta <- dV <- array(dim = c(ncomp, J, n))
  # X0 <- X;  y0 <- y
  # mean.y <- mean(y)
  # y <- scale(y, scale = FALSE)
  # center <- fdata.cen(fdataobj)
  # mean.X <- center$meanX
  # repnX <- rep(1, n)
  # X <- center$Xcen$data
  
  plsr <- plsfit(fdataobj, y, ncomp=ncomp, lambda=lambda, P=P,...)   
  
  
 
  # scores[, 1:Jmin] <- inprod.fdata(Xcen.fdata, vs, ...)
  # V2 <- plsr$loading.weights
  # X2 <- fdata(X,tt,rtt,nam)
  # scores <- inprod.fdata(X2,V2,...) # se hace dentro de plsfit
  
  l <- DoF <- 1:ncomp
  # Yhat <- plsr$fitted.values    
  # yhat <- sum(Yhat)^2          
  # }    
  # colnames(scores) <- paste("PLS", l, sep = "")
  # outlist = list(call = C, df = DoF, rotation=V2, x=scores, lambda=lambda,P=P,
  #                norm=norm, type="pls", fdataobj=fdataobj,y=y0, l=l,
  #                fdataobj.cen=center$Xcen, mean=mean.X, Yhat = Yhat, yhat = yhat
  #                ,basis.y=basis.y, basis=basis, basis.x=basis.x
  # )
  colnames(plsr$coefs) <- paste("PLS", l, sep = "")
  plsr$rotation <- plsr$basis
  plsr$x <- plsr$coefs
  plsr$call <- C
  plsr$df <- DoF
  plsr$norm <- norm 
  plsr$y <- y
  plsr$fdataobj <- fdataobj
  plsr$basis$names$main <-  plsr$type <- "pls"
  class(plsr) <- "fdata.comp"
  return(plsr)
}



fdata2pls.old<-function(fdataobj, y, ncomp = 2,lambda = 0,
                    P=c(0,0,1), norm = TRUE,...) {
  if (!is.fdata(fdataobj)) fdataobj<-fdataobj(fdataobj)
  C <- match.call()
  X <- fdataobj$data
  tt <- fdataobj[["argvals"]]
  rtt<- fdataobj[["rangeval"]]
  nam <- fdataobj[["names"]]
  J <- ncol(X);n <- nrow(X)
  Jmin <- min(c(J,n))
  Jmax <- min(J+1,n-1)
  Beta <- matrix(, J,ncomp)
  W <- V <- Beta
  dW <- dBeta <- dV <- array(dim = c(ncomp, J, n))
  X0 <- X;  y0 <- y
  mean.y <- mean(y)
  y <- scale(y, scale = FALSE)
  center <- fdata.cen(fdataobj)
  mean.X <- center$meanX
  repnX <- rep(1, n)
  X <- center$Xcen$data
  
  #   if (norm)    {
  #     sd.X <- sqrt(apply(X, 2, var))
  #     X <- X/(rep(1, nrow(X)) %*% t(sd.X))
  #     X2 <- fdata(X,tt,rtt,nam)
  #     dcoefficients = NULL
  #     tX <- t(X)    
  #     A <- crossprod(X)
  #     b <- crossprod(X,y)
  # print(dim(A))    
  # print(dim(b))    
  #     if (lambda>0) {
  #       if (is.vector(P))  {     P<-P.penalty(tt,P)          }
  #       #     else {
  #       dimp<-dim(P)
  #       if (!(dimp[1]==dimp[2] & dimp[1]==J))
  #         stop("Incorrect matrix dimension P")
  #       #         }
  #       M <- solve( diag(J) + lambda*P)
  #       W[, 1]<- M %*%b
  #     }
  #     else    {
  #       M <- NULL
  #       W[, 1] <- b
  #     }
  #     dV[1, , ] <- dW[1, , ] <- dA(W[, 1], A, tX)
  #     W[, 1] <- W[, 1]/sqrt((sum((W[, 1]) * (A %*% W[,1]))))
  #     V[, 1] <- W[, 1]
  #     Beta[, 1] <- sum(V[, 1] * b) * V[, 1]
  #     dBeta[1, , ] <-dvvtz(V[, 1], b, dV[1, , ],tX)
  #     if (ncomp>1) {
  #       for (i in 2:ncomp) {
  #         vsi <-b - A %*% Beta[, i - 1]
  #         if (!is.null(M))  vsi<- M %*% vsi
  #         W[, i] <- vsi
  #         dW[i, , ] <- t(X) - A %*% dBeta[i - 1, , ]
  #         V[, i] <- W[, i] - vvtz(V[, 1:(i - 1), drop = FALSE],A %*% W[, i])
  #         dV[i, , ] = dW[i, , ] - dvvtz(V[, 1:(i - 1),drop = FALSE],
  #                                       A %*% W[, i], dV[1:(i - 1), , , drop = FALSE], A %*% dW[i, , ])
  #         dV[i, , ] <- dA(V[, i], A, dV[i, , ])
  #         V[, i] <- V[, i]/sqrt((sum(t(V[, i]) %*% A %*% V[,i])))
  #         Beta[, i] = Beta[, i - 1] + sum(V[, i] * b) * V[,i]
  #         dBeta[i,,]<-dBeta[i-1,,]+dvvtz(V[,i],b,dV[i,,],tX)
  #       }
  #     }
  #     dcoefficients <- NULL
  #     dcoefficients <- array(0, dim = c(ncomp + 1, J, n))
  #     dcoefficients[2:(ncomp + 1), , ] = dBeta
  #     sigmahat <- RSS <- yhat <- vector(length = ncomp + 1)
  #     DoF <- 1:(ncomp + 1)
  #     Yhat <- matrix(, n, ncomp + 1)
  #     dYhat <- array(dim = c(ncomp + 1, n, n))
  #     coefficients <- matrix(0, J, ncomp + 1)
  #     coefficients[, 2:(ncomp + 1)] = Beta/(sd.X %*% t(rep(1, ncomp)))
  #     intercept <- rep(mean.y, ncomp + 1) - t(coefficients) %*% t(mean.X$data)
  #     covariance <- array(0, dim = c(ncomp + 1, J, J))
  #     DD <- diag(1/sd.X)
  #     for (i in 1:(ncomp + 1)) {
  #       Yhat[, i] <- X0 %*% coefficients[, i] + intercept[i]
  #       res <- y0 - Yhat[, i]
  #       yhat[i] <- sum((Yhat[, i])^2)
  #       RSS[i] <- sum(res^2)
  #       dYhat[i, , ]<-X %*%dcoefficients[i, , ] + matrix(1,n, n)/n
  #       DoF[i] <- sum(diag(dYhat[i, , ]))
  #       #        dummy <- (diag(n) - dYhat[i, , ]) %*% (diag(n) - t(dYhat[i, , ]))
  #       #        sigmahat[i] <- sqrt(RSS[i]/sum(diag(dummy)))
  #       #        if (i > 1) {
  #       #                covariance[i, , ] <- sigmahat[i]^2 * DD %*% dcoefficients[i,
  #       #                  , ] %*% t(dcoefficients[i, , ]) %*% DD
  #       #            }
  #     }
  #     V2 <- fdata(t(V)*(rep(1, nrow(t(V))) %*% t(sd.X)),tt,rtt,nam)
  #     V2$data <- sweep(V2$data,1,norm.fdata(V2),"/")
  #     #   V2<-fdata(t(V),tt,rtt,nam)
  #     #   V2$data<-sweep(V2$data,1,norm.fdata(V2),"/")
  #     #    W2<-fdata(t(W),tt,rtt,nam)
  #     #    X3<-fdata(X,tt,rtt,nam)
  #     #    beta.est<-fdata(t(Beta),tt,rtt,nam)
  #     DoF[DoF > Jmax] = Jmax
  #     #    intercept <- as.vector(intercept)    
  #   }
  #   else {
  plsr <- mplsr(X,y,ncomp=ncomp,lambda=lambda,P=P,...)   
  V2 <- fdata(t(plsr$loading.weights),tt,rtt,nam)
  X2 <- fdata(X,tt,rtt,nam)
  DoF <- 1:ncomp
  Yhat <- plsr$fitted.values    
  yhat <- sum(Yhat)^2          
  # }    
  scores <- inprod.fdata(X2,V2,...)
  #   TT = X %*% V  ## son los scores
  l <- 1:ncomp
  colnames(scores) <- paste("PLS", l, sep = "")
  outlist = list(call = C, df = DoF, rotation=V2, x=scores, lambda=lambda,P=P,
                 norm=norm,type="pls", fdataobj=fdataobj,y=y0, l=l,
                 fdataobj.cen=center$Xcen,mean=mean.X,Yhat = Yhat,yhat = yhat
                 #,basis.y=basis.y,basis=basis,basis.x=basis.x
  )
  #    Yhat = Yhat, coyhat = yhat,efficients = coefficients,intercept = intercept,
  #     RSS = RSS,TT=TT, sigmahat = sigmahat,covariance = covariance,
  #     W2=W2,beta.est=beta.est,
  class(outlist) <- "fdata.comp"
  return(outlist)
}

#######################
mplsr <- function(X, Y, ncomp = 2, lambda=0, P=c(0,0,1),...)
{
  #
  # Orthogonal Scores Algorithm for PLS (Martens and Naes, pp. 121--123)
  #
  # X: predictors (matrix)
  #
  # Y: multivariate response (matrix)
  #
  # K: The number of PLS factors in the model which must be less than or
  #    equal to the  rank of X.
  #
  # Returned Value is the vector of PLS regression coefficients
  #
  dx <- dim(X)
  J<-dx[2]
  if (is.fdata(X)) {
    X<-X$data
    arg<-X$argvals
  }
  else arg<-1:J
  K<-ncomp
  tol <- 1e-10
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  nbclass <- ncol(Y)
  #	xbar <- apply(X, 2, sum)/dx[1]
  #	ybar <- apply(Y, 2, sum)/dx[1]
  xbar <- colSums(X)/dx[1]
  ybar <- colSums(Y)/dx[1]
  X0 <- X - outer(rep(1, dx[1]), xbar)
  Y0 <- Y - outer(rep(1, dx[1]), ybar)
  Xtotvar <- sum(X0 * X0)
  PP<-W <- matrix(0, dx[2], K)
  Q <- matrix(0, nbclass, K)
  #	sumofsquaresY <- apply(Y0^2, 2, sum)
  sumofsquaresY <-colSums(Y0^2)
  u <- Y0[, order(sumofsquaresY)[nbclass]]
  TT <- matrix(NA, dx[1], K)
  tee <- 0
  cee<-numeric(K)
  M<-NULL
  if (lambda>0) { 
    if (is.vector(P))  {     P <- P.penalty(arg,P)          }
    M <- solve( diag(J) + lambda*P)
  }                              	
  for(i in 1:K) {
    test <- 1 + tol
    while(test > tol) {
      w <- crossprod(X0, u)
      if (!is.null(M)) { w <- M %*% w}		
      w <- w/sqrt(crossprod(w)[1])
      W[, i] <- w
      teenew <- X0 %*% w
      test <- sum((tee- teenew)^2)       #norm.fdata(tee,teenew)
      tee<- teenew			
      TT[,i] <-tee
      cee[i] <- crossprod(tee)[1]
      p <- crossprod(X0, (tee/cee[i]))
      PP[, i] <- p
      q <- crossprod(Y0, tee)[, 1]/cee[i]
      u <- Y0 %*% q
      u <- u/crossprod(q)[1]
    }
    Q[, i] <- q
    X0 <- X0 - tee %*% t(p)
    Y0 <- Y0 - tee %*% t(q)
  }
  tQ=solve(crossprod(PP, W)) %*% t(Q)	
  COEF <- W %*% tQ
  b0 <- ybar - t(COEF) %*% xbar   
  
  #	fitted <- drop(tee) *drop(tQ) + rep(ybar, each = dx[1])
  residuals <- Y0     
  fitted <- drop(Y-Y0)
  #Yscores = u,Yloadings=t(q),  
  #list(b0 = b0, COEF = COEF,scores=tee,loadings=p,loading.weights = W,
  list(b0 = b0,  "Yloadings" = COEF, scores=TT, loadings=PP, 
       loading.weights = W, projection=W, 
       Xmeans=xbar, Ymeans=ybar, fitted.values = fitted,
       residuals = residuals,cee=cee,Xvar=colSums(PP*PP)*cee,Xtotvar=Xtotvar)
}

#' @title Principal components for functional data
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
#' @param ncomp Number of principal components.
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
#' \email{manuel.oviedo@@udc.es}
#' @seealso See Also as \link[base]{svd} and \link[stats]{varimax}.
#' @references Venables, W. N. and B. D. Ripley (2002). \emph{Modern Applied
#' Statistics with S}. Springer-Verlag.
#' 
#' N. Kraemer, A.-L. Boulsteix, and G. Tutz (2008). Penalized Partial Least
#' Squares with Applications to B-Spline Transformations and Functional Data.
#' Chemometrics and Intelligent Laboratory Systems, 94, 60 - 69.
#' \doi{10.1016/j.chemolab.2008.06.009}
#' 
#' Febrero-Bande, M., Oviedo de la Fuente, M. (2012).  \emph{Statistical
#' Computing in Functional Data Analysis: The R Package fda.usc.} Journal of
#' Statistical Software, 51(4), 1-28. \url{https://www.jstatsoft.org/v51/i04/}
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
    d <- d * delta
  }  #  else {    newd <- d  }
  scores[, 1:Jmin] <- inprod.fdata(Xcen.fdata, vs, ...)
  colnames(scores) <- paste("PC", 1:J, sep = "")
  l <- 1:ncomp
  vs$names$main <-   "pls"
  out <- list(call = C,d = d, basis = vs[1:ncomp], rotation = vs[1:ncomp],
              coefs = scores, x = scores,
              lambda = lambda,P=P, fdataobj.cen = Xcen.fdata,norm=norm, type="pc",
              mean = xmean, fdataobj = fdataobj,l=l, u=u[,1:ncomp,drop=FALSE])
  class(out) = "fdata.comp"
  return(out)
}


D.penalty=function(tt) {
  #tt <- 1:2
  rtt <- diff(tt)
  p <- length(tt)
  hh <- -(1/mean(1/rtt))/rtt
  if (p==2)  Dmat1 <- -diag(p-1)    else Dmat1 <- diag(hh)
  mat <- matrix(rep(0,p-1),ncol=1)
  Dmat12 <- cbind(Dmat1,mat)
  Dmat22 <- cbind(mat,-Dmat1)
  Dmat <- Dmat12 + Dmat22
  Dmat
  return(Dmat)
}


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
#' 69. \doi{10.1016/j.chemolab.2008.06.009}
#' @keywords math
#' @examples
#' 
#' P.penalty((1:10)/10,P=c(0,0,1))
#' # a more detailed example can be found under script file 
#' @export 
P.penalty <- function(tt,P=c(0,0,1)) {
  lenp <- length(P)
  if (length(tt)==1) return(matrix(1,1,1))
  D.pen2 <- D.pen <- D.penalty(tt)
  p <- ncol(D.pen)
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




#########################################
# Wrapper of Orthogonal scores PLSR
# Fits a PLSR model with the orthogonal scores algorithm (aka the NIPALS algorithm).
#
# X: predictors (matrix)
#
# Y: multivariate response (matrix)
#
# K: The number of PLS factors in the model which must be less than or
#    equal to the  rank of X.
#
# Returned Value is the vector of PLS regression coefficients
# 
# Pendiente centrar Y y X usando fdata.cen


plsfit <- function (X, Y, ncomp=3, lambda= 0, P=c(0,0,1),norm=TRUE, maxit = 100) 
{
  tol = .Machine$double.eps^0.5
  dx <- dim(X)
  nobj <- dx[1]
  npred <- dx[2]
  isxfdata <- is.fdata(X)
  if (isxfdata) {
  aux=fdata.cen(X)
  Xmean=aux$meanX
  X0=aux$Xcen$data }
  else {
  X=data.matrix(X)
  Xmean=colMeans(X)
  X0=sweep(X,2,Xmean,"-")
  }
 
  
#  Y <- data.matrix(Y)
  dy <- dim(Y)
  nresp <- dy[2]
  
  isyfdata <- is.fdata(Y)
  if (isyfdata) {
      aux <- fdata.cen(Y)
      Ymean <- aux$meanX
      Y0 <- aux$Xcen$data
	  Ymat=Y$data
    }
	else {
	Ymat=data.matrix(Y)
	Ymean=colMeans(Ymat)
	Y0=sweep(Ymat,2,Ymean,"-")
	}

  dy <- dim(Ymat)
  nresp <- dy[2]
    
  M <- NULL
  #penaltyMatrix <- matrix(0,npred,npred)
  if (lambda>0) {
    if (is.vector(P))  {     P <- P.penalty(X$argvals,P)          }
    dimp <- dim(P)
    if (!(dimp[1]==dimp[2] & dimp[1]==npred))
      stop("Incorrect matrix dimension P")
    # M <- solve( diag(npred) + lambda * P)
    M <- Minverse( diag(npred) + lambda * P)
    penaltyMatrix <- P 
  } else  penaltyMatrix <- matrix(0,npred,npred) # o M
  
  
  dnX <- dimnames(X0)
  dnY <- dimnames(Y0)

  
  tQ <- matrix(0, nrow = ncomp, ncol = nresp)
  B <- array(0, dim = c(npred, nresp, ncomp))
  TT <- U <- matrix(0, nrow = nobj, ncol = ncomp)
  tsqs <- numeric(ncomp)
  fitted <- residuals <- array(0, dim = c(nobj, nresp, 
                                            ncomp))

  W <- P <- matrix(0, nrow = npred, ncol = ncomp)
  Xtotvar <- sum(X0 * X0)
  for (a in 1:ncomp) {
    if (nresp == 1) {
      u.a <- Y0
    }
    else {
      u.a <- Y0[, which.max(colSums(Y0 * Y0))]
      t.a.old <- 0
    }
    nit <- 0
    repeat {
      nit <- nit + 1
      w.a <- crossprod(X0, u.a)
      # Penalization
      if (!is.null(M)) { w.a <- M %*% w.a}		
      
      w.a <- w.a/sqrt(c(crossprod(w.a)))
      t.a <- X0 %*% w.a
      tsq <- c(crossprod(t.a))
      t.tt <- t.a/tsq
      q.a <- crossprod(Y0, t.tt)
      if (nresp == 1) 
        break
      if (sum(abs((t.a - t.a.old)/t.a), na.rm = TRUE) < 
          tol) 
        break
      else {
        u.a <- Ymat %*% q.a/c(crossprod(q.a))
        t.a.old <- t.a
      }
      if (nit >= maxit) {
        warning("No convergence in ", maxit, " iterations\n")
        break
      }
    }
    p.a <- crossprod(X0, t.tt)
    X0 <- X0 - t.a %*% t(p.a)
    Y0 <- Y0 - t.a %*% t(q.a)
    W[, a] <- w.a
    P[, a] <- p.a
    tQ[a, ] <- q.a
    TT[, a] <- t.a
    U[, a] <- u.a
    tsqs[a] <- tsq
    fitted[, , a] <- TT[, 1:a] %*% tQ[1:a, , drop = FALSE]
    residuals[, , a] <- Y0
   }
  if (ncomp == 1) {
    R <- W
  }
  else {
    PW <- crossprod(P, W)
    if (nresp == 1) {
      PWinv <- diag(ncomp)
      bidiag <- -PW[row(PW) == col(PW) - 1]
      for (a in 1:(ncomp - 1)) PWinv[a, (a + 1):ncomp] <- cumprod(bidiag[a:(ncomp - 1)])
    }
    else {
      PWinv <- backsolve(PW, diag(ncomp))
    }
    R <- W %*% PWinv
  }
  for (a in 1:ncomp) {
    B[, , a] <- R[, 1:a, drop = FALSE] %*% tQ[1:a, , drop = FALSE]
  }
#   fitted <- fitted + rep(Ymean, each = nobj)
    if (isyfdata) {
	fitted <- sweep(fitted,2,Ymean$data,"+")
	} else {
   fitted <- sweep(fitted,2,Ymean,"+")
	}
    objnames <- dnX[[1]]
    if (is.null(objnames)) objnames <- dnY[[1]]
    prednames <- dnX[[2]]
    respnames <- dnY[[2]]
    compnames <- paste("Comp", 1:ncomp)
    nCompnames <- paste(1:ncomp, "comps")
    dimnames(TT) <- dimnames(U) <- list(objnames, compnames)
    dimnames(R) <- dimnames(W) <- dimnames(P) <- list(prednames, 
                                                      compnames)
    dimnames(tQ) <- list(compnames, respnames)
    dimnames(B) <- list(prednames, respnames, nCompnames)
    dimnames(fitted) <- dimnames(residuals) <- list(objnames,respnames, nCompnames)

    Pvar <- colSums(P * P)
    P <- t(P)
    W <- t(W)
    R <- t(R)
#    Yloadings <- t(tQ)

    
    if (isxfdata){
       P <- fdata(P,X$argvals,X$rangeval,X$names)
       W <- fdata(W,X$argvals,X$rangeval,X$names)
	   R <- fdata(R,X$argvals,X$rangeval,X$names)
	   TT <-t(Minverse(inprod.fdata(P))%*%inprod.fdata(P,X-Xmean))
    }
    
    if (isyfdata){
      Yloadings <- fdata(tQ,Y$argvals,Y$rangeval,Y$names)
	  U <- t(Minverse(inprod.fdata(Yloadings))%*%inprod.fdata(Yloadings,Y-Ymean))
    } else {
	Yloadings=tQ
	}
	if (norm) {
	if (isxfdata){
		P$data=sweep(P$data,1,norm.fdata(P),"/")
		TT=t(Minverse(inprod.fdata(P))%*%inprod.fdata(P,X-Xmean))
		# Xtilde=gridfdata(TT,P,Xmean)
	} else {
		P=P/sqrt(apply(P,1,function(v){sum(v^2)}))
		TT=t(Minverse(tcrossprod(P,P))%*%P%*%t(sweep(X,2,Xmean,"-")))
		#Xtilde=sweep(TT%*%P,2,Xmean,"+")
	}
	if (isyfdata){
		Yloadings$data=sweep(Yloadings$data,1,norm.fdata(Yloadings),"/")
		U=t(Minverse(inprod.fdata(Yloadings))%*%inprod.fdata(Yloadings,Y-Ymean))
		#Ytilde=gridfdata(U,Yloadings,Ymean)
		dimnames(U) <- list(objnames, compnames)
	} else {
		Yloadings=tQ/sqrt(apply(tQ,1,function(v){sum(v^2)}))
		U=t(Minverse(tcrossprod(Yloadings,Yloadings))%*%Yloadings%*%t(sweep(Ymat,2,Ymean,"-")))
		dimnames(U) <- list(objnames, compnames)
		#Ytilde=sweep(U%*%Yloadings,2,Ymean,"+")
	}	
	}

    list(Regcoefs = B, coefs = TT, basis = P,
		 loading.weigths = W, 
         coefs.y = U, basis.y = Yloadings, projection = R, 
         mean = Xmean, mean.Y = Ymean, l = 1:ncomp,
         # Xmeans = Xmeans, 
         # fitted.values = fitted, residuals = residuals, 
         # Xvar = Pvar * tsqs,         Xtotvar = Xtotvar,
         penaltyMatrix=penaltyMatrix)
}

