#' @title Fitting Functional Linear Models
#' 
#' @description Computes functional regression between functional (and non functional)
#' explanatory variables and scalar response using basis representation.
#' 
#' @details This section is presented as an extension of the linear regression models:
#' \code{\link{fregre.pc}}, \code{\link{fregre.pls}} and
#' \code{\link{fregre.basis}}. Now, the scalar response \eqn{Y} is estimated by
#' more than one functional covariate \eqn{X^j(t)} and also more than one non
#' functional covariate \eqn{Z^j}. The regression model is given by:
#' \deqn{E[Y|X,Z]=\alpha+\sum_{j=1}^{p}\beta_{j}Z^{j}+\sum_{k=1}^{q}\frac{1}{\sqrt{T_k}}\int_{T_k}{X^{k}(t)\beta_{k}(t)dt}
#' }{E[Y|X,Z]=\alpha+\sum_j \beta_j Z^j + \sum_k <X^k,\beta_k>}
#' 
#' where \eqn{Z=\left[ Z^1,\cdots,Z^p \right]}{Z=[Z^1,...,Z^p]} are the non
#' functional covariates, \eqn{X(t)=\left[ X^{1}(t_1),\cdots,X^{q}(t_q)
#' \right]}{X(t)=[X^1(t),...,X^q(t)]} are the functional ones and
#' \eqn{\epsilon} are random errors with mean zero , finite variance
#' \eqn{\sigma^2} and \eqn{E[X(t)\epsilon]=0}{E[X(t)\epsilon]=0}.  
#' 
#' The first item in the \code{data} list is called \emph{"df"} and is a data
#' frame with the response and non functional explanatory variables, as
#' \code{\link{lm}}. Functional covariates of class \code{fdata} or \code{fd}
#' are introduced in the following items in the \code{data} list.\cr
#' 
#' \code{basis.x} is a list of basis for represent each functional covariate.
#' The basis object can be created by the function:
#' \code{\link{create.pc.basis}}, \code{\link{pca.fd}}
#' \code{\link{create.pc.basis}}, \code{\link{create.fdata.basis}} or
#' \code{\link{create.basis}}.\cr \code{basis.b} is a list of basis for
#' represent each functional \eqn{\beta_k} parameter. If \code{basis.x} is a
#' list of functional principal components basis (see
#' \code{\link{create.pc.basis}} or \code{\link{pca.fd}}) the argument
#' \code{basis.b} \emph{(is unnecessary and)} is ignored.\cr
#' 
#' Penalty options are under development, not guaranteed to work properly.
#' The user can penalty the basis elements by: (i) \code{lambda} is a list of
#' rough penalty values of each functional covariate, see
#'  \code{\link{P.penalty}} for more details. 
#' 
#' @param formula an object of class \code{formula} (or one that can be coerced
#' to that class): a symbolic description of the model to be fitted. The
#' details of model specification are given under \code{Details}.
#' @param data List that containing the variables in the model. 
#' Functional covariates are recommended to be of class fdata. 
#' Objects of class "fd" can be used at the user's own risk.
#' @param basis.x List of basis for functional explanatory data estimation.
#' @param basis.b List of basis for functional beta parameter estimation.
#' @param lambda List, indexed by the names of the functional covariates, 
#' which contains the Roughness penalty parameter. 
#' @param P List, indexed by the names of the functional covariates, which contains the parameters for the creation of the penalty matrix.
#' @param weights weights
#' @param \dots Further arguments passed to or from other methods.
#' @return Return \code{lm} object plus:
#' \itemize{
#' \item \code{sr2}{ Residual variance.}
#' \item \code{Vp}{ Estimated covariance matrix for the parameters.} 
#' \item \code{lambda}{ A roughness penalty.} 
#' \item \code{basis.x}{ Basis used for \code{fdata} or \code{fd} covariates.} 
#' \item \code{basis.b}{ Basis used for beta parameter estimation.}
#' \item \code{beta.l}{ List of estimated beta parameter of functional covariates.}
#' \item \code{data}{ List that containing the variables in the model.}
#' \item \code{formula}{ formula.}
#' }
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@usc.es}
#' @seealso See Also as: \code{\link{predict.fregre.lm}} and
#' \code{\link{summary.lm}}.\cr Alternative method: \code{\link{fregre.glm}}.
#' @references Ramsay, James O., and Silverman, Bernard W. (2006), \emph{
#' Functional Data Analysis}, 2nd ed., Springer, New York.
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
#' tt <- x[["argvals"]]
#' dataf <- as.data.frame(tecator$y)
#' 
#' nbasis.x <- 11
#' nbasis.b <- 5
#' basis1 <- create.bspline.basis(rangeval=range(tt),nbasis=nbasis.x)
#' basis2 <- create.bspline.basis(rangeval=range(tt),nbasis=nbasis.b)
#' basis.x <- list("x"=basis1)
#' basis.b <- list("x"=basis2)
#' f <- Fat ~ Protein + x
#' ldat <- ldata("df"=dataf,"x"=x)
#' res <- fregre.lm(f,ldat,  basis.b=basis.b)
#' summary(res)
#' f2 <- Fat ~ Protein + xd +xd2
#' xd <- fdata.deriv(x,nderiv=1,class.out='fdata', nbasis=nbasis.x)
#' xd2 <- fdata.deriv(x,nderiv=2,class.out='fdata', nbasis=nbasis.x)
#' ldat2 <- list("df"=dataf,"xd"=xd,"x"=x,"xd2"=xd2)
#' basis.x2 <- NULL#list("xd"=basis1)
#' basis.b2 <- NULL#list("xd"=basis2)
#' basis.b2 <- list("xd"=basis2,"xd2"=basis2,"x"=basis2)

#' res2 <- fregre.lm(f2, ldat2,basis.b=basis.b2)
#' summary(res2)
#' par(mfrow=c(2,1))
#' plot(res$beta.l$x,main="functional beta estimation")
#' plot(res2$beta.l$xd,col=2)
#' }
#' @rdname fregre.lm
#' @export
fregre.lm <- function(formula, data, basis.x = NULL, basis.b = NULL
                       #,rn, lambda
                       , lambda=NULL, P=NULL
                       , weights=rep(1,n),...){
  
  # print("fregre.lm")
  tf <- terms.formula(formula)
  terms <- attr(tf, "term.labels")
  nt <- length(terms)
  if (attr(tf, "response") > 0) {
    response <- as.character(attr(tf, "variables")[2])
    pf <- rf <- paste(response, "~", sep = "")
  } else pf <- rf <- "~"
  vtab <- rownames(attr(tf,"factors"))                                                    
  vnf <- intersect(terms,names(data$df))
  vfunc2 <- setdiff(terms,vnf)
  vint <- setdiff(terms,vtab)
  vfunc <- setdiff(vfunc2,vint)
  off <- attr(tf,"offset")
  name.coef <-  nam  <-  par.fregre <- beta.l <- list()
  kterms <- 1
  n <- length(data[["df"]][,response])
  XX <- data.frame(data[["df"]][,c(response)],weights)
  
  namxx <- names(XX) <- c(response,"weights")
  lenvnf <- length(vnf)
  df <- 0
  intercept <- attr(tf,"intercept")==1
  # if (missing(P))    {    P  <-rn0=FALSE
  #       rn=list()} else rn0<-TRUE
  # if (missing(lambda))    {    lambda0=FALSE
  #       lambda=list()} else lambda0<-TRUE
  
  # if (missing(P))  P <- list()
  # if (missing(lambda))  lambda <- list()
  
  lenvfunc <- length(vfunc)
  hay.pls<-FALSE
  
  # out <- fdata2model(vfunc, vnf, response, data, basis.x,basis.b,pf,tf)
  # P <- list("x"=1,"x2"=1); lambda <- list("x"=1,"x2"=1)
  #P <- lambda <- NULL
  #out <- fdata2model.penalty(vfunc, vnf, response, data, basis.x,basis.b,pf,tf,lambda,P)
  #print("entra fdata2model")
  out <- fdata2model.penalty(vfunc = vfunc, vnf=vnf, response=response,
							 data=data, basis.x = basis.x, basis.b=basis.b, 
                             pf=pf, tf=tf, lambda=lambda, P=P)
  #print("sale fdata2model")
  mean.list <- out$mean.list
  name.coef <- out$name.coef
#  bsp1<-out$bsp1 # No tiene sentido que sea común a todas las variables.
  pf <- out$pf
  XX <- out$XX
  basis.x <- out$basis.x
  basis.b <- out$basis.b
  penalty <- out$penalty
  # print("penalty")
  # print(penalty)
  # print(head(XX))
  #scores<-as.matrix(XX[,-(1:2)])  
  #if (!is.data.frame(XX)) XX=data.frame(XX)
  # print("pasa 1")
  
  par.fregre$formula=pf
  par.fregre$data=XX
  y <- XX[,1] 
  ycen = y - mean(y)
  
  #Z <- as.matrix(XX[,-1])     
  W <- diag(weights)  
  if (!penalty) {
   if (lenvfunc==0 & length(vnf)==0)      {
      z=lm(formula=pf,data=XX,x=TRUE,...)   
      class(z)<-c("lm","fregre.lm")
      return(z)
    }       else       z <- lm(formula = pf,data=XX,x=TRUE,...)
#  print(summary(res))
    e <- z$residuals
    z$coefs<- summary(z)$coefficients
    z$r2 <- 1 - sum(z$residuals^2)/sum(ycen^2)  
    z$lambda <- FALSE
    # if (intercept)   {
    #   Z <- cbind(rep(1,len=n),Z)      
    #   colnames(Z)[1]<-"(Intercept)"
    # }
    
    
    #cbind(z$coefficients,coefs[,1],coefs3,coef.qr)
    class(z) <- c(class(z),"fregre.lm")
  } 
  else {
    # print("lm.penaltyyyy")
    # z <- lm.penalty(Z, W, scores, XX, basis.x, lambda0, rn0)
    # print(out$lpenalty)
    # print("out$lpenalty")
    # print(out$ipenalty)
    z<- lm.penalty( XX,W, vfunc, basis.x, out$lpenalty,out$ipenalty)
    z$lambda <- TRUE
    #  z$fitted.values
    # out$ipenalty
    #  z$mat
    # z$fitted.values-XX[,1]
  }       
  #    z$call<-z$call[1:2]
  if (length(vfunc)>0){
  for (i in 1:length(vfunc)) {
    z$coefficients[is.na(z$coefficients)]<-0
     if (inherits(basis.b[[vfunc[i]]],"basisfd")) 
      beta.l[[vfunc[i]]]=fd(z$coefficients[name.coef[[vfunc[i]]]],basis.b[[vfunc[i]]])
    else{
      if(!is.null(basis.b[[vfunc[i]]]$basis)) {
        #     beta.est<-z$coefficients[name.coef[[vfunc[i]]]]*vs.list[[vfunc[i]]]
#        beta.est <- z$coefficients[ name.coef[[vfunc[i]]]] * basis.list[[vfunc[i]]]
        beta.est <- gridfdata(matrix(z$coefficients[name.coef[[vfunc[i]]]],nrow=1),basis.b[[vfunc[i]]]$basis)
#        beta.est$data<-colSums(beta.est$data)
     
         beta.est$names$main<-bquote(paste(hat(beta),"(",.(vfunc[i]),")",sep=""))
#        beta.est$data <- matrix(as.numeric(beta.est$data),nrow=1)
#        beta.est$names$main<-"beta.est"
#        beta.est$data <- matrix(as.numeric(beta.est$data),nrow=1)

#        if  (basis.b[[vfunc[i]]]$type=="pls") {
#          if (basis.b[[vfunc[i]]]$norm)  {
#            sd.X <- sqrt(apply(data[[vfunc[i]]]$data, 2, var))
#            beta.est$data<-  beta.est$data/sd.X
#          }      
#        }  

        beta.l[[vfunc[i]]]<-beta.est     
      }
      else {
#        beta.est<-z$coefficients[name.coef[[vfunc[i]]]]*t(basis.list[[vfunc[i]]])
        #     beta.est<-apply(beta.est,2,sum)
        beta.est<-fd(z$coefficients[name.coef[[vfunc[i]]]],basis.b[[vfunc[i]]]$harmonics$basis)
#        beta.est<-colSums(beta.est)
#        beta.l[[vfunc[i]]]<-fd(beta.est,basis.x[[vfunc[i]]]$harmonics$basis)
        beta.l[[vfunc[i]]]<-beta.est
      }
    }
  }
  }
  # z$H <- design2hat(Z,W)  # usarla en fdata2model
  # print("design2hat")
  # print(names(z))
  # z$yp <- z$H %*% y
  #    z$coefficients <- design2coefs(Z,W,y)  # usarla en fdata2model
  # print("design2coefs"); print(z$coefficients)
  # print(name.coef)
  e <- z$residuals
  z$sr2 <- sum(e^2)/z$df.residual
  ###################### z$Vp <- z$sr2*S
  z$Vp <- z$sr2 * z$H  # 20210321

  z$formula <- pf
  z$formula.ini <- formula
  z$basis.x <- basis.x
  z$basis.b <- basis.b
  z$mean <- mean.list
  if (length(vfunc)>0){
      z$beta.l <- beta.l
  } else {
  z$beta.l=NULL
  }
  z$vs.list <- out$vs.list 
  #
  #z$data <- z$data
  z$XX <- XX
  # print("pasa 3") 
  z$data <- data
  #z$fdataobj <- data[[vfunc[1]]] #Por qué es necesario esto?
  #z$rn <- rn0
  if (is.list(z$lambda.opt)) z$lambda <- TRUE
  #z$JJ <- vs.list   
  
  class(z) <- c("fregre.lm","lm")
  # class(z$beta.l) <- c("mfdata","list")
  z
}     
################# auxiliary functions: 
# design2coefs(); design2hat() lm.penalty(); fdata2model.penalty (in fdata2model.R)
design2coefs <- function(Z,W,y){
  tZ <- t(Z)      
  # option 1
  # A  <- tZ%*%sqrt(W)
  #         qr.solve(t(A),y) },
  
  # option 2
  # A0 <- tZ %*% W %*% Z 
  # S0 <- solve(A0)
  # Cinv <- S0%*%tZ %*% W      
  # coefs <- Cinv %*% y
  
  # option 3
  # S<-diag(coef(summary(z))[,2])
  # qr0<-qr(Z, LAPACK = F)
  # coef.qr <- qr.coef(qr0, y)
  
  qr0<-qr(Z, LAPACK = F)
  coef <- qr.coef(qr0, y)         
  return(coef)
}
############################     
design2hat <- function(Z,W){
  tZ <- t(Z)      
  Sb = tZ %*% W %*% Z# + lambda * R
  Lmat <- chol(Sb)
  Lmatinv <- solve(Lmat)
  Cinv <- Lmatinv %*% t(Lmatinv)
  Sb2 = Cinv %*% tZ
  H = Z %*% Sb2 %*% W
  return(H)
}


############################      
lm.penalty <- function(  XX,W,vfunc
                         #Z, W, scores, XX
                         , basis.x, lpenalty,ipenalty){
  
  # print("lm.penalty")
  #lpenalty <-out$lpenalty
  
  #ipenalty <-out$ipenalty
  # print("si rn0 o si lambda0")
  #      qr0<-solve(qr(t(scores)%*%W%*%scores+mat2, LAPACK = TRUE),ycen)
  #      mat2<-0
  y <- XX[,1]
  n <- length(y)
  x <-  data.matrix(cbind(rep(1,NROW(XX)),XX[,-1]))
  #x <- data.matrix(XX[,-1])
  colnames(x)[1]<-"Intercept"
  head(x)
  rownames(x)<-1:nrow(x)
  #rownames(XX)[1]<-"Intercept"
  lenvfunc <- length(vfunc)
  pp <- ncol(x)
  #        mat <- rep(0,len=pp-2)
  mat <- diag(0,nrow=pp)
  # print(pp);  print(dim(x));  print(ipenalty)
  for (ifun in 1:lenvfunc){
    imat0 <- ipenalty[[vfunc[ifun]]]
    mat[imat0,imat0] <- lpenalty[[vfunc[ifun]]]
  }
  ###################### pc
  
 # mat <- t(x)%*%P%*%P%*%x #R
 #  X<- res1$XX[,3:9]
#  dim(X);dim(P)
 # P <- P.penalty(1:7,c(0,1))
  #P
  #mat <- X%*%P%*%t(X)
  Sb<- t(x)%*%W%*%x + mat
  
  # S<-solve(Sb)       
  #    S=solve(t(Z)%*%W%*%Z)
  #S=solve(t(Z)%*%W%*%Z + lambda*mat)    
  S<-Minverse(Sb) 
  Cinv<-S%*%t(x)%*%W         
  coefs<-Cinv%*%y
  yp<-drop(x%*%coefs)
  H<-x%*%Cinv
   df<-fdata.trace(H)
  #if  (!hay.pls) df <- fdata.trace(H)  
  coefs<-drop(coefs)
  ######################
  ycen <- y-mean(y)
  z <- list()  
  z$coefficients <- coefs
  z$fitted.values <- yp   
  e <- z$residuals <- y - z$fitted.values      
  #z$mean.list<-mean.list
  z$df.residual<-n-df
  z$H <- H
  # z$r2 <- 1 - sum(z$residuals^2)/sum(ycen^2)   # no se pasa ycen!!!
  if  (inherits(basis.x[[vfunc[1]]],"basisfd")) {
    z$call[[1]] = "fregre.basis"
  }  else {
    if  (basis.x[[vfunc[1]]]$type=="pc")   z$call[[1]] = "fregre.pc"
    if  (basis.x[[vfunc[1]]]$type=="pls")  z$call[[1]] = "fregre.pls"        
  }             
  class(z) <- c("fregre.fd","fregre.lm")
  rdf <- n-df
  sr2 <- sum(e^2)/ rdf
  r2 <- 1 - sum(e^2)/sum(ycen^2)
  r2.adj <- 1 - (1 - r2) * ((n -    1)/ rdf)
  GCV <- sum(e^2)/(n - df)^2
  ####################
  z$residuals <- drop(e)
  z$fitted.values <- yp
  z$y <- y
  z$rank <- df
  z$df.residual <-  rdf
  # Z=cbind(rep(1,len=n),Z)
  # colnames(Z)[1] = "(Intercept)"
  std.error = sqrt(diag(S) *sr2)
  t.value = coefs/std.error
  p.value = 2 * pt(abs(t.value), n - df, lower.tail = FALSE)
  coefficients <- cbind(coefs, std.error, t.value, p.value)
  colnames(coefficients) <- c("Estimate", "Std. Error",
                              "t value", "Pr(>|t|)")
  z$coefs<-coefficients    
  z$terms <- terms
  # z$lambda.opt <- lambda
  # z$lambda <- lambda
  z$mat <- mat
  # print(mat)
  class(z) <- "lm"
  return(z)
}
