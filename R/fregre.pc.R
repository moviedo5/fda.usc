#' Functional Regression with scalar response using Principal Components
#' Analysis
#' 
#' Computes functional (ridge or penalized) regression between functional
#' explanatory variable \eqn{X(t)} and scalar response \eqn{Y} using Principal
#' Components Analysis.\cr \deqn{Y=\big<X,\beta\big>+\epsilon=\int_{T}{X(t)\beta(t)dt+\epsilon}}{Y=<X,\beta>+\epsilon} 
#' where \eqn{ \big< \cdot , \cdot \big>}{<.,.>} denotes the inner product on
#' \eqn{L_2} and \eqn{\epsilon} are random errors with mean zero , finite
#' variance \eqn{\sigma^2} and \eqn{E[X(t)\epsilon]=0}{E[X(t)\epsilon]=0}.\cr
#'  
#' @details The function computes the \eqn{\left\{\nu_k\right\}_{k=1}^{\infty}}{\nu_1,...,\nu_\infty} orthonormal
#' basis of functional principal components to represent the functional data as
#' \eqn{X_i(t)=\sum_{k=1}^{\infty}\gamma_{ik}\nu_k}{X(t)=\sum_(k=1:\infty)
#' \gamma_k \nu_k} and the functional parameter as
#' \eqn{\beta(t)=\sum_{k=1}^{\infty}\beta_k\nu_k}{\beta(t)=\sum_(k=1:\infty)
#' \beta_k \nu_k}, where \eqn{\gamma_{ik}=\Big< X_i(t),\nu_k\Big>}{\gamma_k= < X(t), \nu_k >} and
#' \eqn{\beta_{k}=\Big<\beta,\nu_k\Big>}{\beta_k=<\beta,\nu_k>}. \cr 
#' The response can be fitted by: 
#' \itemize{ 
#' \item \eqn{\lambda=0}{\lambda=0}, no penalization,
#' \deqn{\hat{y}=\nu_k^{\top}(\nu_k^{\top}\nu_k)^{-1}\nu_k^{\top}y}{y.est= \nu'(\nu'\nu)^{-1}\nu'y} 
#' \item Ridge regression, \eqn{\lambda>0}{\lambda>0} and \eqn{P=1}{P=1},
#'  \deqn{\hat{y}=\nu_k^{\top}(\nu_k\top \nu_k+\lambda I)^{-1}\nu_k^{\top}y}{y.est=\nu'(\nu'\nu+\lambda I)^{-1}\nu'y}
#' \item Penalized regression, \eqn{\lambda>0}{\lambda>0} and \eqn{P\neq0}{P!=0}. For example, \eqn{P=c(0,0,1)}{P=c(0,0,1)} penalizes the second derivative (curvature) by \code{P=P.penalty(fdataobj["argvals"],P)},
#' \deqn{\hat{y}=\nu_k^{\top}(\nu_k\top \nu_k+\lambda \nu_k^{\top}
#' \textbf{P}\nu_k)^{-1}\nu_k^{\top}y}{y.est=\nu'(\nu'\nu+\lambda v'Pv)^{-1}\nu'y}
#' } 
#' @aliases fregre.pc
#' @param fdataobj \code{\link{fdata}} class object or \code{fdata.comp} class
#' object created\cr by \code{\link{create.pc.basis}} function.
#' @param y Scalar response with length \code{n}.
#' @param l Index of components to include in the model.If is null \code{l} (by
#' default), \code{l=1:3}.
#' @param lambda Amount of penalization. Default value is 0, i.e. no
#' penalization is used.
#' @param P If \code{P} is a vector: \code{P} are coefficients to define the
#' penalty matrix object, see \code{\link{P.penalty}}. If \code{P} is a matrix:
#' P is the penalty matrix object.
#' @param weights weights
#' @param \dots Further arguments passed to or from other methods.
#' @return Return:
#' \itemize{
#' \item \code{call}{ The matched call of \code{\link{fregre.pc}} function.} 
#' \item \code{coefficients}{ A named vector of coefficients.}
#' \item \code{residuals}{ \code{y}-\code{fitted values}.} 
#' \item \code{fitted.values}{ Estimated scalar response.} 
#' \item \code{beta.est}{ beta coefficient estimated of class \code{fdata}}
#' \item \code{df.residual}{ The residual degrees of freedom. In ridge regression, \code{df(rn)} is the effective degrees of freedom.} 
#' \item \code{r2}{ Coefficient of determination.}
#' \item \code{sr2}{ Residual variance.} 
#' \item \code{Vp}{ Estimated covariance matrix for the parameters.} 
#' \item \code{H}{ Hat matrix.}
#' \item \code{l}{ Index of principal components selected.} 
#' \item \code{lambda}{ Amount of shrinkage.} 
#' \item \code{P}{ Penalty matrix.} 
#' \item \code{fdata.comp}{ Fitted object in \code{\link{fdata2pc}} function.} 
#' \item \code{lm}{ \code{lm} object.}
#' \item \code{fdataobj}{ Functional explanatory data.} 
#' \item \code{y}{ Scalar response.}
#' }
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' @seealso See Also as: \code{\link{fregre.pc.cv}},
#' \code{\link{summary.fregre.fd}} and \code{\link{predict.fregre.fd}}.
#' 
#' Alternative method: \code{\link{fregre.basis}} and \code{\link{fregre.np}}.
#' @references Cai TT, Hall P. 2006. \emph{Prediction in functional linear
#' regression}. Annals of Statistics 34: 2159-2179.
#' 
#' Cardot H, Ferraty F, Sarda P. 1999. \emph{Functional linear model}.
#' Statistics and Probability Letters 45: 11-22.
#' 
#' Hall P, Hosseini-Nasab M. 2006. \emph{On properties of functional principal
#' components analysis}. Journal of the Royal Statistical Society B 68:
#' 109-126.
#' 
#' Febrero-Bande, M., Oviedo de la Fuente, M. (2012).  \emph{Statistical
#' Computing in Functional Data Analysis: The R Package fda.usc.} Journal of
#' Statistical Software, 51(4), 1-28. \url{https://www.jstatsoft.org/v51/i04/}
#' 
#' N. Kraemer, A.-L. Boulsteix, and G. Tutz (2008). \emph{Penalized Partial
#' Least Squares with Applications to B-Spline Transformations and Functional
#' Data}. Chemometrics and Intelligent Laboratory Systems, 94, 60 - 69.
#' \doi{10.1016/j.chemolab.2008.06.009}
#' @keywords regression
#' @examples
#' \dontrun{
#' data(tecator)
#' absorp <- tecator$absorp.fdata
#' ind <- 1:129
#' x <- absorp[ind,]
#' y <- tecator$y$Fat[ind]
#' res <- fregre.pc(x,y)
#' summary(res)
#' res2 <- fregre.pc(x,y,l=c(1,3,4))
#' summary(res2)
#' # Functional Ridge Regression
#' res3 <- fregre.pc(x,y,l=c(1,3,4),lambda=1,P=1)
#' summary(res3)
#' # Functional Regression with 2nd derivative penalization
#' res4 <- fregre.pc(x,y,l=c(1,3,4),lambda=1,P=c(0,0,1))
#' summary(res4)
#' betas <- c(res$beta.est,res2$beta.est,
#'            res3$beta.est,res4$beta.est)
#' plot(betas)
#' } 
#' 
#' @export fregre.pc
fregre.pc=function (fdataobj, y, l =NULL,lambda=0,P=c(0,0,1),
                    weights=rep(1,len=n),...){
  if (is(fdataobj,"fdata.comp")) {
    pc <- fdataobj
    fdataobj <- pc$fdataobj
    if (is.null(l))    {
      l <- pc$l
    }
    else if (length(l)>nrow(pc$basis)) stop("Incorrect value for  argument l")
    x <- fdataobj[["data"]]
    tt <- fdataobj[["argvals"]]                                                    
  }
  else {
    if (is.null(l)) l <- 1:3
    if (!is.fdata(fdataobj))    fdataobj <- fdata(fdataobj)
    #  omit<-omit.fdata(fdataobj,y)
    #  fdataobj<-omit[[1]]
    #  y<-omit[[2]]
    tt <- fdataobj[["argvals"]]
    x <- fdataobj[["data"]]
    pc <- fdata2pc(fdataobj,ncomp=max(l),lambda=lambda,P=P)
  }  
  rtt <- fdataobj[["rangeval"]]
  names <- fdataobj[["names"]]
  n <- nrow(x)
  np <- ncol(x)
  lenl = length(l)
  if (is.null(rownames(x)))        rownames(x) <- 1:n
  X <- xcen <- pc$fdataobj.cen
  if (n != (length(y)))   stop("ERROR IN THE DATA DIMENSIONS")
  C <- match.call()
  ycen <- y - mean(y)
  vs <-t (pc$basis$data[l,,drop=F])
  scores <- Z <- (pc$coefs[,l,drop=F])
  cnames <- colnames(pc$coefs)[l]
  df <- lenl+1
  J <- min(np,lenl)
  ymean <- mean(y)
  ycen <-  y - ymean
  W <- diag(weights) 
  if (is.logical(lambda)) {
    #   val <- log(.25*(pc$d[1]^2),base=2)
    lambda <- .25*(pc$d[1]^2)#lambda <- c(0,2^seq(0,val,len=10))
  }
  if (lambda>0) {
    xmean <- pc$mean
    d <- pc$newd[l]
    D <- diag(d)
    diagJ <- diag(J)
    #    lenrn <- length(rn)
    scores <- cbind(rep(1,n),pc$coefs[,l])
    order.deriv <- 0     
    if (!is.matrix(P)){
      if (is.vector(P)) {
        for (i in 1:length(P))   {      if (P[i]!=0)         order.deriv <- i}
        P <- P.penalty(tt,P)
        
        P <- t(vs)%*%P%*%vs 
        P <- P*(diff(rtt)/(np -1))^(order.deriv*2-1)
      }         }
    mat <- diag(J+1)
    mat[-1,-1] <- lambda*P
    mat[1,1] <- 0
    Sb <- t(scores)%*%W%*%scores + mat
    # S <- solve(Sb)       
    #    S=solve(t(Z)%*%W%*%Z)    
    S <- Minverse(Sb) 
    Cinv <- S%*%t(scores)%*%W         
    coefs <- Cinv%*%y
    yp <- drop(scores%*%coefs)
    H <- scores%*%Cinv
    df <- fdata.trace(H)
    coefs <- drop(coefs)
    names(coefs) <- c("Intercept",cnames)
    beta.est <- coefs[-1]*pc$basis[l]
    beta.est$data <- colSums(beta.est$data)
    beta.est$names$main <- "beta.est"
    beta.est$data  <-  matrix(as.numeric(beta.est$data),nrow=1)
    e <- y-yp
    rdf <- n-df
    sr2 <- sum(e^2)/ rdf
    r2 <- 1 - sum(e^2)/sum(ycen^2)
    r2.adj<- 1 - (1 - r2) * ((n -    1)/ rdf)
    #    GCV <- sum(e^2)/(n - df)^2
    object.lm <- list()
    object.lm$coefficients <- coefs
    object.lm$residuals <- drop(e)
    object.lm$fitted.values <- yp
    object.lm$x<-scores 
    object.lm$y <- y
    object.lm$rank <- df
    object.lm$df.residual <-  rdf
    Z <- cbind(rep(1,len=n),Z)
    colnames(Z)[1] <- "(Intercept)"
    std.error <- sqrt(diag(S) *sr2)
    Vp <- sr2*S 
    t.value <- coefs/std.error
    p.value <- 2 * pt(abs(t.value), rdf, lower.tail = FALSE)
    coefficients <- cbind(coefs, std.error, t.value, p.value)
    colnames(coefficients) <- c("Estimate", "Std. Error",
                                "t value", "Pr(>|t|)")
    class(object.lm) <- "lm"
    out <- list(call = C, beta.est = beta.est,coefficients=coefs,
                fitted.values =yp,residuals = e,H=H,df.residual = rdf,r2=r2,#GCV=GCV,
                sr2 = sr2,Vp=Vp,l = l,lambda=lambda,fdata.comp=pc,lm=object.lm,
                scoefs=coefficients,fdataobj = fdataobj,y = y)
    ##################################
  }
  else {
    #print("no rn")
    response <- "y"
    dataf <- data.frame(y,Z,weights)
    colnames(dataf) <- c("y",cnames,"weights")
    pf <- paste(response, "~", sep = "")
    for (i in 1:length(cnames)) 
      pf <- paste(pf,"+",cnames[i],sep="")
    object.lm <- lm(formula = pf,data=data.frame(dataf),weights=weights,x=TRUE, y=TRUE)
    beta.est <- object.lm$coefficients[2:(lenl+1)]*pc$basis[l]
    beta.est$data <- colSums(beta.est$data)
    beta.est$names$main<-"beta.est"
    beta.est$data <- matrix(as.numeric(beta.est$data),nrow=1)
    Z=cbind(rep(1,len=n),Z)
    #    S=solve(t(Z)%*%W%*%Z)    
    S <- t(Z)%*%W%*%Z
    S <- Minverse(S) 
    H <- Z%*%S%*%t(Z)
    e <- object.lm$residuals
    df <- fdata.trace(df)#n- object.lm$df
    sr2 <- sum(e^2)/(n - df)
    Vp <- sr2*S 
    r2 <- 1 - sum(e^2)/sum(ycen^2)
    #     r2.adj<- 1 - (1 - r2) * ((n -    1)/(n-df))
    #     GCV <- sum(e^2)/(n - df)^2
    out <- list(call = C, coefficients=object.lm$coefficients,residuals = e,
                fitted.values =object.lm$fitted.values,weights=weights,beta.est = beta.est,
                df.residual = n-df,r2=r2,sr2 = sr2,Vp=Vp,H=H, l = l,lambda=lambda,P=P,fdata.comp=pc,
                lm=object.lm,XX=Z, fdataobj = fdataobj,y = y)
  }
  class(out) = "fregre.fd"
  return(out)
}

#' Functional Penalized PLS regression with scalar response
#' 
#' @description
#' Computes functional linear regression between functional explanatory variable \eqn{X(t)} and scalar response \eqn{Y} using penalized Partial
#' Least Squares (PLS) \deqn{Y=\big<\tilde{X},\beta\big>+\epsilon=\int_{T}{\tilde{X}(t)\beta(t)dt+\epsilon}}{Y=<\tilde{X},\beta>+\epsilon} where \eqn{ \big< \cdot , \cdot \big>}{<.,.>} denotes the inner product on
#' \eqn{L_2} and \eqn{\epsilon} are random errors with mean zero , finite variance \eqn{\sigma^2} and \eqn{E[\tilde{X}(t)\epsilon]=0}{E[X(t)\epsilon]=0}.\cr
#' \eqn{\left\{\nu_k\right\}_{k=1}^{\infty}}{\nu_1,...,\nu_\infty} orthonormal basis of PLS to represent the functional data as \eqn{X_i(t)=\sum_{k=1}^{\infty}\gamma_{ik}\nu_k}{X(t)=\sum_(k=1:\infty)
#' \gamma_k \nu_k}.%, where \eqn{\tilde{X}=MX} with \eqn{M=(I+\lambda P)^{-1}},\eqn{\gamma_{ik}=\Big< \tilde{X}_i(t),\nu_k\Big>}{\gamma_k= < \tilde{X}_i(t), \nu_k >}. 
#' 
#' 
#' @details Functional (FPLS) algorithm maximizes the covariance between \eqn{X(t)} and the scalar response \eqn{Y} via the partial least squares (PLS) components.
#' The functional penalized PLS are calculated in \code{\link{fdata2pls}} by alternative formulation of the NIPALS algorithm proposed by Kraemer and
#' Sugiyama (2011).\cr 
#' Let \eqn{\left\{\tilde{\nu}_k\right\}_{k=1}^{\infty}}{{\nu_k}_k=1:\infty} the functional PLS components and \eqn{\tilde{X}_i(t)=\sum_{k=1}^{\infty}\tilde{\gamma}_{ik}\tilde{\nu}_k}{X_i(t)=\sum{k=1:\infty}
#' \gamma_{ik} \nu_k} and \eqn{\beta(t)=\sum_{k=1}^{\infty}\tilde{\beta}_k\tilde{\nu}_k}{\beta(t)=\sum{k=1:\infty}
#' \beta_k \nu_k}. The functional linear model is estimated by: \deqn{\hat{y}=\big< X,\hat{\beta} \big> \approx \sum_{k=1}^{k_n}\tilde{\gamma}_{k}\tilde{\beta}_k }{ y.est=< X,\beta.est >
#' \approx \sum{k=1:k_n} \gamma_k \beta_k }\cr 
#' The response can be fitted by: 
#' \itemize{ \item \eqn{\lambda=0}{\lambda=0}, no penalization,
#' \deqn{\hat{y}=\nu_k^{\top}(\nu_k^{\top}\nu_k)^{-1}\nu_k^{\top}y}{y.est= \nu'(\nu'\nu)^{-1}\nu'y}
#' \itemize{
#' \item Penalized regression, \eqn{\lambda>0}{\lambda>0} and \eqn{P\neq0}{P!=0}. For example, \eqn{P=c(0,0,1)}{P=c(0,0,1)} penalizes the
#' second derivative (curvature) by \code{P=P.penalty(fdataobj["argvals"],P)},
#' \deqn{\hat{y}=\nu_k^{\top}(\nu_k\top \nu_k+\lambda \nu_k^{\top} \textbf{P}\nu_k)^{-1}\nu_k^{\top}y}{y.est=\nu'(\nu'\nu+\lambda v'Pv)^{-1}\nu'y} }
#' }
#' 
#' @param fdataobj \code{\link{fdata}} class object.
#' @param y Scalar response with length \code{n}.
#' @param l Index of components to include in the model.
#' @param lambda Amount of penalization. Default value is 0, i.e. no
#' penalization is used.
#' @param P If \code{P} is a vector: \code{P} are coefficients to define the
#' penalty matrix object. By default \code{P=c(0,0,1)} penalize the second
#' derivative (curvature) or acceleration.  If \code{P} is a matrix: P is the
#' penalty matrix object.
#' @param \dots Further arguments passed to or from other methods.
#' 
#' @return Return: 
#' \itemize{
#' \item \code{call}{ The matched call of \code{\link{fregre.pls}} function.} 
#' \item \code{beta.est}{ Beta coefficient estimated of class \code{fdata}.} 
#' \item \code{coefficients}{ A named vector of coefficients.}
#' \item \code{fitted.values}{ Estimated scalar response.} 
#' \item \code{residuals}{\code{y}-\code{fitted values}.} 
#' \item \code{H}{ Hat matrix.} 
#' \item \code{df.residual}{ The residual degrees of freedom.} 
#' \item \code{r2}{ Coefficient of determination.}
#' \item \code{GCV}{ GCV criterion.} 
#' \item \code{sr2}{ Residual variance.} 
#' \item \code{l}{ Index of components to include in the model.} 
#' \item \code{lambda}{ Amount of shrinkage.}
#' \item \code{fdata.comp}{ Fitted object in \code{\link{fdata2pls}} function.}
#' \item \code{lm}{ Fitted object in \code{\link{lm}} function} 
#' \item \code{fdataobj}{ Functional explanatory data.} 
#' \item \code{y}{ Scalar response.}
#' }
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' @seealso See Also as: \code{\link{P.penalty}} and
#' \code{\link{fregre.pls.cv}}.\cr Alternative method: \code{\link{fregre.pc}}.
#' @references Preda C. and Saporta G. \emph{PLS regression on a stochastic
#' process}. Comput. Statist. Data Anal. 48 (2005): 149-158.
#' 
#' N. Kraemer, A.-L. Boulsteix, and G. Tutz (2008). \emph{Penalized Partial
#' Least Squares with Applications to B-Spline Transformations and Functional
#' Data}. Chemometrics and Intelligent Laboratory Systems, 94, 60 - 69.
#' \doi{10.1016/j.chemolab.2008.06.009}
#' 
#' Martens, H., Naes, T. (1989) \emph{Multivariate calibration.} Chichester:
#' Wiley.
#' 
#' Kraemer, N., Sugiyama M. (2011). \emph{The Degrees of Freedom of Partial
#' Least Squares Regression}. Journal of the American Statistical Association.
#' Volume 106, 697-705.
#' 
#' Febrero-Bande, M., Oviedo de la Fuente, M. (2012).  \emph{Statistical
#' Computing in Functional Data Analysis: The R Package fda.usc.} Journal of
#' Statistical Software, 51(4), 1-28. \url{https://www.jstatsoft.org/v51/i04/}
#' @keywords regression
#' @examples
#' \dontrun{
#' data(tecator)
#' x <- tecator$absorp.fdata
#' y <- tecator$y$Fat
#' res <- fregre.pls(x,y,c(1:4))
#' summary(res)
#' res1 <- fregre.pls(x,y,l=1:4,lambda=100,P=c(1))
#' res4 <- fregre.pls(x,y,l=1:4,lambda=1,P=c(0,0,1))
#' summary(res4)#' plot(res$beta.est)
#' lines(res1$beta.est,col=4)
#' lines(res4$beta.est,col=2)

#' }
#' @export fregre.pls
fregre.pls=function(fdataobj, y=NULL, l = NULL,
                    lambda=0,P=c(0,0,1),...){
  if (is(fdataobj,"fdata.comp")) {
    pc <- fdataobj
    fdataobj <- pc$fdataobj
    if  (is.null(l)) l<-1:nrow(pc$basis)
    if  (is.null(y)) y <- pc$y
    else if (all(y!=pc$y)) warning("y is different from that calculated on the pls basis")
  }
  else {
    
    if (is.null(l)) l<- 1:3
    # omit<-omit.fdata(fdataobj,y)
    # fdataobj<-omit[[1]]
    # y<-omit[[2]]
    pc<-fdata2pls(fdataobj,y,ncomp=max(l),lambda=lambda,P=P,...)
  }
  if (length(l)==1) l<-1:l       
  x <- fdataobj[["data"]]
  tt <- fdataobj[["argvals"]]
  rtt <- fdataobj[["rangeval"]]
  names <- fdataobj[["names"]]
  n <- nrow(x)
  np <- ncol(x)
  lenl <- length(l)
  if (n != (length(y)))   
    stop("ERROR IN THE DATA DIMENSIONS")
  C <- match.call()
  if (is.null(rownames(x)))        rownames(x) <- 1:n
  ycen <- y - mean(y)
  vs <- pc$basis$data[,,drop=FALSE]
  Z <- pc$coefs[,l,drop=F]
  xcen <- pc$fdataobj.cen
  cnames <- colnames(pc$coefs)[l]
  response <- "y"
  df <- data.frame(y,Z)
  colnames(df) <- c("y",cnames)
  pf <- paste(response, "~", sep = "")
  for (i in 1:length(cnames)) 
    pf <- paste(pf,"+",cnames[i],sep="")
  object.lm <- lm(formula = pf, data =df , x = TRUE,y = TRUE)
  beta.est <- object.lm$coefficients[2:(lenl+1)]*pc$basis[l]
  beta.est$data <- apply(beta.est$data,2,sum)
  beta.est$names$main <- "beta.est"
  beta.est$data <- matrix(as.numeric(beta.est$data),nrow=1)
  #            if  (pc$type=="pls") {
  if (pc$norm)  {
    sd.X <- sqrt(apply(fdataobj$data, 2, var))
    beta.est$data<-  beta.est$data/sd.X
  }      
  #            }     
  #    H<-diag(hat(Z, intercept = TRUE),ncol=n)
  # H2<-lm.influence(object.lm, do.coef = T)$hat# o bien
  #    I <- diag(1/(n*pc$lambdas[l]), ncol = lenl) #1/n
  Z <- cbind(rep(1,len=n),Z)
  order.deriv <- 0 
  if (lambda==0) mat<-0
  else {      
    if (!is.matrix(P)){
      if (is.vector(P)) {
        for (i in 1:length(P))   {
          if (P[i]!=0)         order.deriv<-i
          }
        P <- P.penalty(tt,P)                
        P <- vs%*%P%*%t(vs) 
        P <- P*(diff(rtt)/(np -1))^(order.deriv*2-1)
      }         }
    mat<-diag(lenl+1)
    mat[-1,-1] <- lambda*P
    mat[1,1] <- 0   
  }
  print(mat)
  S <- t(Z)%*%Z+mat
  S <- Minverse(S)
#  H <- Z%*%S%*%t(Z)
  
#Sb <- t(scores)%*%W%*%scores + mat
# S <- Minverse(Sb) 
 Cinv <- S%*%t(Z)         
 coefs <- Cinv%*%y
 yp <- drop(Z%*%coefs)
 H <- Z%*%Cinv

 
 coefs <- drop(coefs)
 names(coefs) <- c("Intercept",cnames)
 beta.est <- coefs[-1]*pc$basis[l]
 beta.est$data <- colSums(beta.est$data)
 beta.est$names$main <- "beta.est"
 beta.est$data  <-  matrix(as.numeric(beta.est$data),nrow=1)
 
  e <- y - yp
  df <- max(fdata.trace(H),pc$df[lenl]+1)
  rdf <- n-df
  sr2 <- sum(e^2)/rdf
  Vp <- sr2*S 
  r2 <- 1 - sum(e^2)/sum(ycen^2)
  #    r2.adj<- 1 - (1 - r2) * ((n -    1)/ rdf)
  #    GCV <- sum(e^2)/rdf^2             #GCV=GCV,
  std.error <- sqrt(diag(S) *sr2)
  t.value  <- coefs/std.error
  p.value <- 2 * pt(abs(t.value), rdf, lower.tail = FALSE)
  coefficients <- cbind(coefs, std.error, t.value, p.value)
  colnames(coefficients) <- c("Estimate", "Std. Error","t value", "Pr(>|t|)")
  
  out <- list(call = C,coefficients=coefs, residuals = e,
              fitted.values =yp, beta.est = beta.est, scoefs=coefficients,
              H=H,df.residual = rdf,r2=r2, sr2 = sr2, Vp=Vp,l = l,lambda=lambda,P=P, fdata.comp=pc,
              lm=object.lm,fdataobj = fdataobj,y = y)
  class(out) = "fregre.fd"
  return(out)
}

#' Functional penalized PLS regression with scalar response using selection of
#' number of PLS components
#' 
#' @description Functional Regression with scalar response using selection of number of
#' penalized principal componentes PPLS through cross-validation. The algorithm
#' selects the PPLS components with best estimates the response. The selection
#' is performed by cross-validation (CV) or Model Selection Criteria (MSC).
#' After is computing functional regression using the best selection of PPLS
#' components.
#' 
#' @details The algorithm selects the best principal components
#' \code{pls.opt} from the first \code{kmax} PLS and (optionally) the best
#' penalized parameter \code{lambda.opt} from a sequence of non-negative
#' numbers \code{lambda}. 
#' \itemize{ 
#' \item The method selects the best principal components with
#' minimum MSC criteria by stepwise regression using \code{\link{fregre.pls}}
#' in each step.  
#' \item The process (point 1) is repeated for each \code{lambda} value.  
#' \item The method selects the principal components (\code{pls.opt}=\code{pls.order[1:k.min]}) and (optionally) the lambda parameter with minimum MSC criteria.
#' } 
#' Finally, is computing functional PLS regression between functional explanatory variable \eqn{X(t)} and scalar response \eqn{Y} using the best selection of PLS \code{pls.opt} and ridge parameter \code{rn.opt}.  
#' The criteria selection is done by cross-validation (CV) or Model Selection Criteria (MSC).   
#' \itemize{ 
#' \item Predictive Cross-Validation: \eqn{PCV(k_n)=\frac{1}{n}\sum_{i=1}^{n}{\Big(y_i -\hat{y}_{(-i,k_n)}\Big)^2}}{PCV(k_n)=1/n \sum_(i=1:n) (y_i - \hat{y}_{-i})^2}, \code{criteria}=``CV'' 
#' \item Model Selection Criteria: 
#' \eqn{MSC(k_n)=log \left[ \frac{1}{n}\sum_{i=1}^{n}{\Big(y_i-\hat{y}_i\Big)^2} \right] +p_n\frac{k_n}{n} }{MSC(k_n)=log [ 1/n \sum_(i=1:n){ (y_i- \hat{y}_i )^2} ] +p_n k_n/n } \cr
#'  \eqn{p_n=\frac{log(n)}{n}}{p_n=log(n)/n}, \code{criteria}=``SIC'' (by default)\cr
#' \eqn{p_n=\frac{log(n)}{n-k_n-2}}{p_n=log(n)/(n-k_n-2)},
#' \code{criteria}=``SICc''\cr \eqn{p_n=2}, \code{criteria}=``AIC''\cr
#' \eqn{p_n=\frac{2n}{n-k_n-2}}{p_n=2n/(n-k_n-2)}, \code{criteria}=``AICc''\cr
#' \eqn{p_n=\frac{2log(log(n))}{n}}{p_n=2log(log(n))/(n)},
#' \code{criteria}=``HQIC''\cr
#' where \code{criteria} is an argument that controls the
#' type of validation used in the selection of the smoothing parameter
#' \code{kmax}\eqn{=k_n} and penalized parameter \code{lambda}\eqn{=\lambda}.
#' }
#' 
#' @param fdataobj \code{\link{fdata}} class object.
#' @param y Scalar response with length \code{n}.
#' @param kmax The number of components to include in the model.
#' @param lambda Vector with the amounts of penalization. Default value is 0,
#' i.e. no penalization is used.  If \code{lambda=TRUE} the algorithm computes
#' a sequence of lambda values.
#' @param P The vector of coefficients to define the penalty matrix object. For
#' example, if \code{P=c(0,0,1)}, penalized regression is computed penalizing
#' the second derivative (curvature).
#' @param criteria Type of cross-validation (CV) or Model Selection Criteria
#' (MSC) applied. Possible values are \emph{"CV"}, \emph{"AIC"}, \emph{"AICc"},
#' \emph{"SIC"}, \emph{"SICc"}, \emph{"HQIC"}.
#' @param \dots Further arguments passed to \code{\link{fregre.pls}}.
#' @return Return:
#' \itemize{
#' \item \code{fregre.pls}{ Fitted regression object by the best (\code{pls.opt}) components.} 
#' \item \code{pls.opt}{ Index of PLS components' selected.} 
#' \item \code{MSC.min}{ Minimum Model Selection Criteria (MSC) value for
#' the (\code{pls.opt} components.} 
#' \item \code{MSC}{ Minimum Model Selection Criteria (MSC) value for \code{kmax} components.}
#' }
#' @note \code{criteria=``CV''} is not recommended: time-consuming.
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' @seealso See also as:\code{\link{fregre.pc}} .
#' @references Preda C. and Saporta G. \emph{PLS regression on a stochastic
#' process}. Comput. Statist. Data Anal. 48 (2005): 149-158.
#' 
#' Febrero-Bande, M., Oviedo de la Fuente, M. (2012).  \emph{Statistical
#' Computing in Functional Data Analysis: The R Package fda.usc.} Journal of
#' Statistical Software, 51(4), 1-28. \url{https://www.jstatsoft.org/v51/i04/}
#' @keywords regression
#' @examples
#' \dontrun{
#' data(tecator)
#' x<-tecator$absorp.fdata[1:129]
#' y<-tecator$y$Fat[1:129]
#' # no penalization
#' pls1<- fregre.pls.cv(x,y,8)
#' # 2nd derivative penalization
#' pls2<-fregre.pls.cv(x,y,8,lambda=0:5,P=c(0,0,1))
#' }
#' 
#' @export fregre.pls.cv
fregre.pls.cv=function (fdataobj, y, kmax=8,lambda=0,P=c(0,0,1),
                        criteria = "SIC",...) {
  if (is(fdataobj,"fdata.comp")) {
    pc<-fdataobj
    fdataobj<-pc$fdataobj
    kmax<-nrow(pc$basis)
  }
  else {
    if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
    omit<-omit.fdata(fdataobj,y)
    fdataobj<-omit[[1]]
    y<-omit[[2]]
  }
  x<-fdataobj[["data"]]
  tt<-fdataobj[["argvals"]]
  rtt<-fdataobj[["rangeval"]]
  n <- nrow(x);    nc <- ncol(x)
  if (min(lambda)>0) {
    if (is.vector(P)) {    P2<-P.penalty(tt,P=P)   }
    if (is.logical(lambda[1]))   {
      normx<-sqrt(sum(abs((t(x)%*%x)^2)))
      normp<-sqrt(sum(abs(P2^2)))
      normip<-sqrt(sum(abs(diag(nc)+P2)^2))
      pii <-fdata.trace(P2)
      lambda0<-(-2*pii+sqrt(4*pii^2-4*(nc-normx^2)*normp^2))/(2*normp^2)
      #  lambda1<-normx/normp
      lambda<-seq(0,sqrt(lambda0),len=10)
    }
  }
  tol<-sqrt(.Machine$double.eps)
  lenrn<-length(lambda)
  ind =1:kmax
  l = l2 = list()
  ck = 1
  tab = list("AIC", "AICc","SIC", "SICc","HQIC","CV")
  type.i = pmatch(criteria, tab)
  MSC.min<-Inf
  cv.AIC <- matrix(NA,nrow=lenrn,ncol=kmax)
  if (is.na(type.i))     stop("Error: incorrect criteria")
  else {
    if (type.i < 6) {
      #        cv.AIC <- rep(NA, kmax)
      for (r in 1:lenrn) {    
        pls<-fdata2pls(fdataobj,y,ncomp=kmax,lambda=lambda[r],P=P,...)
        for (j in 1:kmax) {
          pls2<-pls
          pls2$basis<-pls$basis[1:j]
          out = fregre.pls(pls2,y,lambda=lambda[r],P=P,...)
          ck<-n-out$df.residual
          s2 <- sum(out$residuals^2)/n  #(n-ck)
          cv.AIC[r,j]<-switch(criteria,
                              "AIC"=log(s2) + 2 * (ck)/n,
                              "AICc"=log(s2) + 2 * (ck)/(n - ck - 2),
                              "SIC"=log(s2) + log(n) * ck/n,
                              "SICc"=log(s2) + log(n) * ck/(n-ck-2),
                              "HQIC"=log(s2) + 2*log(log(n)) * ck/n
          )
          if ( MSC.min>(cv.AIC[r,j]+tol)) {
            rn.opt<-r
            pc.opt<-j
            MSC.min= cv.AIC[r,j]
          }
        }
        #    min.AIC = min(cv.AIC)
        #    pc.opt <- which.min(cv.AIC)
      }
    }
    # CV criteria
    else {
      #        pc2<-pc
      for (j in 1:kmax) {
        residuals2<-rep(NA,n)
        for (r in 1:lenrn) {
          for (i in 1:n){
            out = fregre.pls(fdataobj[-i], y[-i],lambda=lambda[r],P=P,...)
            ck<-n-out$df.residual
            a1<-out$coefficients[1]
            out$beta.est$data<-matrix(out$beta.est$data,nrow=1)
            b1<-inprod.fdata(fdata.cen(fdataobj[i],out$fdata.comp$mean)[[1]],out$beta.est)
            yp<- a1+b1
            #            residuals[i] <- y[i] - yp
            residuals2[i] <- ((y[i] - yp)/(n-ck))^2
          }
          #         cv.AIC[j] <- mean(residuals^2)/(n-j)^2###
          cv.AIC[r,j] <-sum(residuals2)/n
          
          if ( MSC.min>cv.AIC[r,j]) {
            rn.opt<-r
            pc.opt<-j
            MSC.min= cv.AIC[r,j]
          }
        }   }
    }   }
  colnames(cv.AIC) = paste("PLS",1:kmax , sep = "")
  rownames(cv.AIC) = paste("lambda=",signif(lambda,4) , sep = "")
  #    pc2$basis<-pc$basis[1:pc.opt]
  fregre=fregre.pls(fdataobj,y,l=1:pc.opt,lambda=lambda[rn.opt],P=P,...) #B.B bug detected
  MSC.min = cv.AIC[rn.opt,pc.opt]
  return(list("fregre.pls"=fregre,pls.opt = 1:pc.opt,lambda.opt=lambda[rn.opt],
              MSC.min = MSC.min,MSC = cv.AIC))
}

#' Functional penalized PC regression with scalar response using selection of
#' number of PC components
#' 
#' @description Functional Regression with scalar response using selection of number of
#' (penalized) principal components PC through cross-validation. The algorithm
#' selects the PC with best estimates the response. The selection is performed
#' by cross-validation (CV) or Model Selection Criteria (MSC). After is
#' computing functional regression using the best selection of principal
#' components.
#' 
#' @details The algorithm selects the best principal components \code{pc.opt} from the first \code{kmax} PC and (optionally) the best penalized parameter \code{lambda.opt} from a sequence of non-negative
#' numbers \code{lambda}. \cr  
#' If \code{kmax} is a integer (by default and recomended) the procedure is as follows (see example 1):
#' \itemize{
#'  \item Calculate the best principal component (\emph{pc.order[1]}) between \code{kmax} by
#' \code{\link{fregre.pc}}.
#'  \item Calculate the second-best principal component (\code{pc.order [2]}) between the \code{(kmax-1)} by
#' \code{\link{fregre.pc}} and calculate the criteria value of the two
#' principal components.  
#' \item The process (point 1 and 2) is repeated until \code{kmax} principal component (\emph{pc.order[kmax]}). 
#' \item The proces (point 1, 2 and 3) is repeated for each \code{lambda} value.
#' \item The method selects the principal components (\code{pc.opt}=\code{pc.order[1:k.min]}) and (optionally) the lambda
#' parameter with minimum MSC criteria.
#' }
#' If \code{kmax} is a sequence of integer the procedure is as follows (see example 2): 
#' \itemize{
#'  \item The method selects the best principal components with minimum MSC criteria by
#' stepwise regression using \code{\link{fregre.pc}} in each step.  
#' \item The process (point 1) is repeated for each \code{lambda} value.
#' \item The method selects the principal components (\code{pc.opt}=\code{pc.order[1:k.min]}) and (optionally) the lambda
#' parameter with minimum MSC criteria. 
#' }
#' Finally, is computing functional PC regression between functional explanatory variable \eqn{X(t)} and scalar
#' response \eqn{Y} using the best selection of PC \code{pc.opt} and ridge
#' parameter \code{rn.opt}.  \cr  
#' The criteria selection is done by cross-validation (CV) or Model Selection
#' Criteria (MSC).  
#' \itemize{
#' \item Predictive Cross-Validation:
#' \eqn{PCV(k_n)=\frac{1}{n}\sum_{i=1}^{n}{\Big(y_i -\hat{y}_{(-i,k_n)}
#' \Big)^2}}{PCV(k_n)=1/n \sum_(i=1:n) (y_i - \hat{y}_{-i})^2},\cr
#' \code{criteria}=``CV''
#' \item Model Selection Criteria: \eqn{MSC(k_n)=log \left[
#' \frac{1}{n}\sum_{i=1}^{n}{\Big(y_i-\hat{y}_i\Big)^2} \right]
#' +p_n\frac{k_n}{n} }{MSC(k_n)=log [ 1/n \sum_(i=1:n){ (y_i- \hat{y}_i )^2} ]
#' +p_n k_n/n } 
#' }
#' \eqn{p_n=\frac{log(n)}{n}}{p_n=log(n)/n}, \code{criteria}=``SIC'' (by default)\cr
#' \eqn{p_n=\frac{log(n)}{n-k_n-2}}{p_n=log(n)/(n-k_n-2)}, \code{criteria}=``SICc''\cr 
#' \eqn{p_n=2}, \code{criteria}=``AIC''\cr
#' \eqn{p_n=\frac{2n}{n-k_n-2}}{p_n=2n/(n-k_n-2)}, \code{criteria}=``AICc''\cr
#' \eqn{p_n=\frac{2log(log(n))}{n}}{p_n=2log(log(n))/(n)}, \code{criteria}=``HQIC''\cr
#'  where \code{criteria} is an argument that controls the
#' type of validation used in the selection of the smoothing parameter
#' \code{kmax}\eqn{=k_n} and penalized parameter
#' \code{lambda}\eqn{=\lambda}.
#' 
#' @param fdataobj \code{\link{fdata}} class object.
#' @param y Scalar response with length \code{n}.
#' @param kmax The number of components to include in the model.
#' @param lambda Vector with the amounts of penalization. Default value is 0,
#' i.e. no penalization is used.  If \code{lambda=TRUE} the algorithm computes
#' a sequence of lambda values.
#' @param P The vector of coefficients to define the penalty matrix object. For
#' example, if \code{P=c(1,0,0)}, ridge regresion is computed and if
#' \code{P=c(0,0,1)}, penalized regression is computed penalizing the second
#' derivative (curvature).
#' @param criteria Type of cross-validation (CV) or Model Selection Criteria
#' (MSC) applied. Possible values are \emph{"CV"}, \emph{"AIC"}, \emph{"AICc"},
#' \emph{"SIC"}, \emph{"SICc"}, \emph{"HQIC"}.
#' @param weights weights
#' @param \dots Further arguments passed to \code{\link{fregre.pc}} or
#' \code{\link{fregre.pls}}
#' @return Return:
#' \itemize{
#' \item \code{fregre.pc}{ Fitted regression object by the best (\code{pc.opt}) components.} 
#' \item \code{pc.opt}{ Index of PC components selected.} 
#' \item \code{MSC.min}{ Minimum Model Selection Criteria (MSC) value for the (\code{pc.opt} components.} 
#' \item \code{MSC}{ Minimum Model Selection Criteria (MSC) value for \code{kmax} components.}
#' }
#' @note \code{criteria=``CV''} is not recommended: time-consuming.
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' @seealso See also as:\code{\link{fregre.pc}} .
#' @references Febrero-Bande, M., Oviedo de la Fuente, M. (2012).
#' \emph{Statistical Computing in Functional Data Analysis: The R Package
#' fda.usc.} Journal of Statistical Software, 51(4), 1-28.
#' \url{https://www.jstatsoft.org/v51/i04/}
#' @keywords regression
#' @examples
#' \dontrun{
#' data(tecator)
#' x<-tecator$absorp.fdata[1:129]
#' y<-tecator$y$Fat[1:129]
#' # no penalization
#'  res.pc1=fregre.pc.cv(x,y,8)
#' # 2nd derivative penalization
#'  res.pc2=fregre.pc.cv(x,y,8,lambda=TRUE,P=c(0,0,1))
#' # Ridge regression
#' res.pc3=fregre.pc.cv(x,y,1:8,lambda=TRUE,P=1) 
#' }
#' 
#' @export fregre.pc.cv
fregre.pc.cv = function (fdataobj, y, kmax=8,lambda=0,P=c(0,0,1),criteria = "SIC",weights=rep(1,len=n),...) {
  sequen=FALSE
  if (length(kmax)>1) {
    sequen=TRUE
    l<-kmax
    kmax<-max(kmax)
  }
  if (is(fdataobj,"fdata.comp")) {
    fdataobj<-fdataobj$fdataobj
    if (min(lambda)!=0 | !is.null(P)) warning("The arguments lambda and P are not used")
  }
  else {
    if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
    tt<-fdataobj[["argvals"]]
    x<-fdataobj[["data"]]
    X<-fdata.cen(x)[[1]]$data
    np<-ncol(x)
    if (min(lambda,na.rm=TRUE)>0) {
      if (is.logical(lambda[1]))   {
        if (is.vector(P)) {    P2<-P.penalty(tt,P=P)   }
        normx<-sqrt(sum(abs((t(x)%*%x)^2)))
        normp<-sqrt(sum(abs(P2^2)))
        normip<-sqrt(sum(abs(diag(np)+P2)^2))
        pii <-fdata.trace(P2)
        lambda0<-(-2*pii+sqrt(4*pii^2-4*(np-normx^2)*normp^2))/(2*normp^2)
        #  lambda1<-normx/normp
        #print(lambda0)
        #print(lambda1)
        lambda<-seq(0,sqrt(lambda0),len=10)
      }
    }
  }
  n=nrow(fdataobj)
  #  pc<-fdata2ppc(fdataobj,ncomp=kmax,lambda=lambda,P=P,...)
  if (is.null(names(y))) names(y)<-1:length(y)
  rtt<-fdataobj[["rangeval"]]
  #  n <- nrow(x);    #nc <- ncol(x)
  cv.opt1 = Inf;    pc.opt1 = NA
  c1 = matrix(1:kmax, nrow = 1)
  num.pc = nrow(c1)
  l = l2 = list()
  max.c = length(c1)
  c0 = 1:kmax
  use = rep(FALSE, kmax)
  tab = list("AIC", "AICc","SIC", "SICc","HQIC","CV")
  type.i = pmatch(criteria, tab)
  #pc2<-pc
  lenrn<-length(lambda)
  MSC3<-list()
  pc.opt2 <- matrix(NA,nrow=lenrn,ncol=kmax)
  rownames(pc.opt2)<-paste("lambda=",signif(lambda,4),sep="")
  colnames(pc.opt2)<-paste("PC(",1:kmax,")",sep="")
  MSC2<-pc.opt2
  MSC.min<-Inf
  min.rn<-lambda[1]
  if (is.na(type.i))     stop("Error: incorrect criteria")
  else {
    if (type.i < 6) {
      for (r in 1:lenrn) {
        pc<-fdata2pc(fdataobj,ncomp=kmax,lambda=lambda[r],P=P)
        #       pc<-fdata2pc(fdataobj,ncomp=kmax,...)
        #       if (!is.matrix(P)) if (is.vector(P)) P<-P.penalty(tt,P)
        cv.opt1 = Inf
        pc.opt1 = NA
        l = l2 = list()
        c1 = matrix(1:kmax, nrow = 1)
        num.pc = nrow(c1)
        max.c = length(c1)
        c0 = 1:kmax
        use = rep(FALSE, kmax)
        pc2<-pc
        for (k in 1:kmax) {
          cv.AIC <- rep(NA, max.c)
          if (sequen) {
            max.c=1
            c1<-matrix(pc$l[1:k],ncol=1)
          }
          for (j in 1:max.c) {
            pc2$basis <- pc$basis#[c1[, j]]
            pc2$l <- pc$l[c1[, j]]
            out = fregre.pc(pc2, y,l=c1[, j],lambda=lambda[r],P=P,weights=weights,...)
            ck<-n-out$df.residual
            s2 <- sum(out$residuals^2)/n
            cv.AIC[j]<-switch(criteria,
                              "AIC"=log(s2) + 2 * (ck)/n,
                              "AICc"=log(s2) + 2 * (ck)/(n - ck - 2),
                              "SIC"=log(s2) + log(n) * ck/n,
                              "SICc"=log(s2) + log(n) * ck/(n-ck-2),
                              "HQIC"=log(s2) + 2*log(log(n)) * ck/n
            )
          }
          if (!sequen){
            min.AIC = min(cv.AIC)
            pc.opt1 <- c1[, which.min(cv.AIC)[1]]
            l[[k]] = pc.opt1[k]
            l2[[k]] = min.AIC
            use[pc.opt1[k]] = TRUE
            l[[k + 1]] = c0[use == FALSE]
            c1 = t(expand.grid(l))
            ck = nrow(c1) + 1
            max.c = ncol(c1)
          }
          else {
            pc.opt1 <- 1:k
            l[[k]] = k
            l2[[k]] =drop(cv.AIC[1])
          }
        }
        mn = which.min(l2)[1]
        MSC = as.numeric(l2)
        if ( MSC.min>MSC[mn]) {
          min.rn<-r
          MSC.min = MSC[mn]
          pc.opt3<-pc.opt1[1:mn]
        }
        pc.opt = pc.opt1[1:mn]
        MSC2[r,]<-MSC
        pc.opt2[r,]<-pc.opt1
      }
    }
    #### CV criteria
    else {
      pb=txtProgressBar(min=0,max=lenrn,width=50,style=3)
      pcl<-list()
      for (r in 1:lenrn) {
        setTxtProgressBar(pb,r-0.5)
        for (i in 1:n) {
          pcl[[i]]<-fdata2pc(fdataobj[-i,],ncomp=kmax,lambda=lambda[r],P=P)
        }
        cv.opt1 = Inf
        pc.opt1 = NA
        l = l2 = list()
        c1 = matrix(1:kmax, nrow = 1)
        num.pc = nrow(c1)
        max.c = length(c1)
        c0 = 1:kmax
        use = rep(FALSE, kmax)
        for (k in 1:kmax) {
          if (sequen) {   max.c=1;   c1<-matrix(pc$l[1:k],ncol=1)   }
          cv.AIC <- rep(NA, max.c)
          cv.AIC2 <- matrix(NA,nrow=max.c,ncol=lenrn)
          
          rownames(cv.AIC2)<-1:max.c
          for (j in 1:max.c) {
            residuals2<- rep(NA, n)
            maxk<-max(c1[, j])
            for (i in 1:n){
              pc2<-pcl[[i]]
              pc2$basis<-pcl[[i]]$basis#[c1[,j]]
              pc2$l<-pcl[[i]]$l[c1[,j]]
              out = fregre.pc(pc2,y[-i],l=c1[,j],weights=weights[-i],...) #####
              ck<-n-out$df.residual
              residuals2[i] <- ((y[i] - predict(out,fdataobj[i,]))/(n-ck))^2
            }
            cv.AIC[j] <-sum(residuals2)/n
          }
          if (!sequen){
            min.AIC = min(cv.AIC)
            pc.opt1 <- c1[, which.min(cv.AIC)[1]]
            l[[k]] = pc.opt1[k]
            l2[[k]] = min.AIC
            use[pc.opt1[k]] = TRUE
            l[[k + 1]] = c0[use == FALSE]
            c1 = t(expand.grid(l))
            ck = nrow(c1) + 1
            max.c = ncol(c1)
          }
          else { pc.opt1 <- 1:k; l[[k]] = k;l2[[k]] =drop(cv.AIC[1]) }
        }
        mn = which.min(l2)[1]
        MSC = as.numeric(l2)
        if ( MSC.min>MSC[mn]) {
          min.rn<-r
          MSC.min = MSC[mn]
          pc.opt3<-pc.opt1[1:mn]
        }
        pc.order<-names(MSC)
        pc.opt = pc.opt1[1:mn]
        MSC2[r,]<-MSC
        pc.opt2[r,]<-pc.opt1
        setTxtProgressBar(pb,r)
      }
      close(pb)
    }
    if (all(is.na(MSC))) stop("System is computationally singular: try with other number of basis elements")
    mn = which.min(l2)[1]
    MSC = as.numeric(l2)
    names(pc.opt3)<-paste("PC", pc.opt3, sep = "")
    rn.opt<-lambda[min.rn]
    fregre=fregre.pc(fdataobj,y,l=drop(pc.opt3),lambda=rn.opt,P=P,weights=weights,...)
    return(list("fregre.pc"=fregre,pc.opt = pc.opt3,lambda.opt=rn.opt,
                PC.order=pc.opt2,MSC.order=MSC2))
  }
}

