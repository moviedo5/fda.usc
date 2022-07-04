#' Fitting Functional Generalized Linear Models
#' 
#' Computes functional generalized linear model between functional covariate
#' \eqn{X^j(t)}{X_j(t)} (and non functional covariate \eqn{Z^j}{Z_j}) and
#' scalar response \eqn{Y} using basis representation.
#' 
#' This function is an extension of the linear regression models:
#' \code{\link{fregre.lm}} where the \eqn{E[Y|X,Z]} is related to the linear
#' prediction \eqn{\eta} via a link function \eqn{g(.)}.
#' 
#' \deqn{E[Y|X,Z]=\eta=g^{-1}(\alpha+\sum_{j=1}^{p}\beta_{j}Z^{j}+\sum_{k=1}^{q}\frac{1}{\sqrt{T_k}}\int_{T_k}{X^{k}(t)\beta_{k}(t)dt})}{E[Y|X,Z]=
#' \eta = g^{-1}(\alpha + \sum \beta_j Z_j+\sum < X_k(t) , \beta_k(t) >)}
#' 
#' where \eqn{Z=\left[ Z^1,\cdots,Z^p \right]}{ Z = [Z_1 ,..., Z_p]} are the
#' non functional covariates and \eqn{X(t)=\left[ X^{1}(t_1),\cdots,X^{q}(t_q)
#' \right]}{X(t) = [ X_1(t_1) ,..., X_q(t_q)]} are the functional ones.
#' 
#' The first item in the \code{data} list is called \emph{"df"} and is a data
#' frame with the response and non functional explanatory variables, as
#' \code{\link{glm}}.\cr
#' 
#' Functional covariates of class \code{fdata} or \code{fd} are introduced in
#' the following items in the \code{data} list.\cr \code{basis.x} is a list of
#' basis for represent each functional covariate. The basis object can be
#' created by the function: \code{\link{create.pc.basis}}, \code{\link{pca.fd}}
#' \code{\link{create.pc.basis}}, \code{\link{create.fdata.basis}} o
#' \code{\link{create.basis}}.\cr \code{basis.b} is a list of basis for
#' represent each \eqn{\beta(t)} parameter. If \code{basis.x} is a list of
#' functional principal components basis (see \code{\link{create.pc.basis}} or
#' \code{\link{pca.fd}}) the argument \code{basis.b} is ignored.
#' 
#' %When using functional data derived recommend using a number of basis to
#' represent beta lower than the number of basis used to represent the
#' functional data.
#' 
#' @param formula an object of class \code{formula} (or one that can be coerced
#' to that class): a symbolic description of the model to be fitted. The
#' details of model specification are given under \code{Details}.
#' @param family a description of the error distribution and link function to
#' be used in the model. This can be a character string naming a family
#' function, a family function or the result of a call to a family function.
#' (See \code{\link{family}} for details of family functions.)
#' @param data List that containing the variables in the model.
#' @param basis.x List of basis for functional explanatory data estimation.
#' @param basis.b List of basis for \eqn{\beta(t)} parameter estimation.
# @param CV =TRUE, Cross-validation (CV) is done . Este parámetro no está incluido.
#' @param subset an optional vector specifying a subset of observations to be
#' used in the fitting process.
#' @param weights weights
#' @param \dots Further arguments passed to or from other methods.
#' @return Return \code{glm} object plus:
#' \itemize{ 
#' \item \code{basis.x}{ Basis used for \code{fdata} or \code{fd} covariates.} 
#' \item \code{basis.b}{ Basis used for beta parameter estimation.} 
#' \item \code{beta.l}{ List of estimated beta parameter of functional covariates.} 
#' \item \code{data}{ List that containing the variables in the model.} 
#' \item \code{formula}{ formula.} 
# \item \code{CV}{ predicted response by cross-validation.}
#' }
#' @note If the formula only contains a non functional explanatory variables
#' (multivariate covariates), the function compute a standard \code{\link{glm}}
#' procedure.
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' @seealso See Also as: \code{\link{predict.fregre.glm}} and
#' \code{\link{summary.glm}}.\cr Alternative method if
#' \code{family}=\emph{gaussian}: \code{\link{fregre.lm}}.
#' @references Ramsay, James O., and Silverman, Bernard W. (2006), \emph{
#' Functional Data Analysis}, 2nd ed., Springer, New York.
#' 
#' McCullagh and Nelder (1989), \emph{Generalized Linear Models} 2nd ed.
#' Chapman and Hall.
#' 
#' Venables, W. N. and Ripley, B. D. (2002) \emph{Modern Applied Statistics
#' with S}, New York: Springer.
#' @keywords regression
#' @examples
#' \dontrun{ 
#' data(tecator)
#' x=tecator$absorp.fdata
#' y=tecator$y$Fat
#' tt=x[["argvals"]]
#' dataf=as.data.frame(tecator$y)
#' nbasis.x=11
#' nbasis.b=7
#' basis1=create.bspline.basis(rangeval=range(tt),nbasis=nbasis.x)
#' basis2=create.bspline.basis(rangeval=range(tt),nbasis=nbasis.b)
#' f=Fat~Protein+x
#' basis.x=list("x"=basis1)
#' basis.b=list("x"=basis2)
#' ldata=list("df"=dataf,"x"=x)
#' res=fregre.glm(f,family=gaussian(),data=ldata,basis.x=basis.x,
#' basis.b=basis.b)
#' summary(res)
#' }

#' @export
fregre.glm=function(formula,family = gaussian(), data, 
                    basis.x=NULL, basis.b=NULL,# CV=FALSE, 
                    subset = NULL,
                    weights= NULL,
                    ...)
  # 
  {
    # print("fregre.glm")
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
    # XX <- data.frame(data[["df"]][,c(response)],weights)
    # namxx <- names(XX) <- c(response,"weights")
    XX <- data.frame(data[["df"]][,c(response),drop=FALSE])
    lenvnf <- length(vnf)
    df <- 0
    intercept <- attr(tf,"intercept")==1
     
    lenvfunc <- length(vfunc)
    hay.pls<-FALSE
    
    
    out <- fdata2model.penalty(vfunc = vfunc, vnf=vnf, response=response
                               , data=data, basis.x = basis.x, basis.b=basis.b, 
                               pf=pf, tf=tf)
    
    # print("sale fdata2model")
    
#    vs.list <- out$vs.list
    mean.list <- out$mean.list
    name.coef <- out$name.coef
#    bsp1<-out$bsp1
    pf <- out$pf
    XX <- out$XX
    basis.x <- out$basis.x
    basis.b <- out$basis.b
    par.fregre$formula=pf
    par.fregre$data=XX
    y <- XX[,1] 
    ycen = y - mean(y)
    if (missing(weights)) weights=rep(1,len=n) 
    
    
    #Z <- as.matrix(XX[,-1])     
    if (lenvfunc==0 & length(vnf)==0)      {
      XX$weights<-weights
        z=glm(formula=formula(pf),data=XX,family=family,
              x=TRUE,y=TRUE,subset=subset
              ,weights = weights)# ,...)
        class(z)<-c("fregre.glm",class(z))
        return(z)
      }       else      {
        XX$weights<-weights
        z=glm(formula=pf,data=XX,
                               family=family,
                               weights = weights,
                               #subset=subset,
                               x=TRUE,y=TRUE)
                                #,...)
                               
                               
       
      }
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
    #} 
    # if (missing(weights)) weights <- rep(1,len=n)
    # W <- diag(weights)  
      
  for (i in 1:length(vfunc)) {
     z$coefficients[is.na(z$coefficients)]<-0
     if (inherits(basis.b[[vfunc[i]]],"basisfd")) 
      beta.l[[vfunc[i]]]=fd(z$coefficients[name.coef[[vfunc[i]]]],basis.b[[vfunc[i]]])
    else{
      if(!is.null(basis.b[[vfunc[i]]]$basis)) {
        #     beta.est<-z$coefficients[name.coef[[vfunc[i]]]]*vs.list[[vfunc[i]]]
#        beta.est <- z$coefficients[ name.coef[[vfunc[i]]]] * basis.list[[vfunc[i]]]
        z$coefficients[name.coef[[vfunc[i]]]]
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

    #  z$H <- design2hat(Z,W)  # usarla en fdata2model
    # z$yp <- z$H %*% y
    #    z$coefficients <- design2coefs(Z,W,y)  # usarla en fdata2model
    # print("design2coefs"); print(z$coefficients)
    # print(name.coef)
    e <- z$residuals
    z$sr2 <- sum(e^2)/z$df.residual
    ###################### z$Vp <- z$sr2*S
    z$Vp <- z$sr2 * z$H  # 20210321
    z$beta.l <- beta.l
    z$formula <- pf
    z$mean <- mean.list
    z$formula.ini <- formula
    z$basis.x <- basis.x
    z$basis.b <- basis.b
    z$vs.list <- out$vs.list
    z$data <- z$data
    z$XX <- XX
    # print("pasa 3") 
    z$data <- data
    z$fdataobj <- data[[vfunc[1]]]
    #z$rn <- rn0
    
    #z$JJ <- vs.list   
#    z$basis.list <- basis.list   
    class(z) <- c("fregre.glm",class(z))
    #class(z)<-c(class(z),"fregre.glm")
    # class(z$beta.l) <- c("mfdata","list")
    z
  } 

