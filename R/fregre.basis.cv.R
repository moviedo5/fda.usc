#' Cross-validation Functional Regression with scalar response using basis
#' representation.
#' 
#' @description  Computes functional regression between functional explanatory variables and
#' scalar response using basis representation.  
#' 
#' @details
#' The function \code{fregre.basis.cv()}  uses validation criterion defined by argument
#' \code{type.CV} to estimate the number of basis elements and/or the penalized
#' parameter (\code{lambda}) that best predicts the response.
#' 
#' If \code{basis = NULL} creates bspline basis.\cr
#' 
#' If the functional covariate \code{fdataobj} is in a format \code{raw data},
#' such as matrix or data.frame, creates an object of class \code{fdata} with
#' default attributes, see \code{\link{fdata}}.\cr
#' 
#' If \code{basis.x} is a vector of number of basis elements and
#' \code{basis.b=NULL}, the function force the same number of elements in the
#' basis of \code{x} and \code{beta}.\cr
#' 
#' If \code{basis.x$type=``fourier''} and \code{basis.b$type=``fourier''}, the
#' function decreases the number of fourier basis elements on the
#' \eqn{min(k_{n1},k_{n2})}{min(k.x,k_b)}, where \eqn{k_{n1}}{k.x} and
#' \eqn{k_{n2}}{k.b} are the number of basis element of \code{basis.x} and
#' \code{basis.b} respectively.
#' 
#' @param fdataobj \code{\link{fdata}} class object.
#' @param y Scalar response with length \code{n}.
#' @param basis.x Basis for functional explanatory data \code{fdataobj}.
#' @param basis.b Basis for functional beta parameter.
#' @param type.basis A vector of character string which determines type of
#' basis. By default \emph{"bspline"}. It is only used when \code{basis.x} or
#' \code{basis.b} are a vector of number of basis considered.
#' @param lambda A roughness penalty. By default, no penalty \code{lambda=0}.
#' @param Lfdobj See \link[fda]{eval.penalty}.
#' @param type.CV Type of cross-validation. By default generalized
#' cross-validation \code{\link{GCV.S}} method.
#' @param par.CV List of parameters for \code{type.CV}: \code{trim}, the alpha
#' of the trimming and \code{draw}.
#' @param weights weights
#' @param verbose If \code{TRUE} information about the procedure is printed.
#' Default is \code{FALSE}.
#' @param \dots Further arguments passed to or from other methods.
#' @return Return:
#' \itemize{
#' \item \code{fregre.basis}: Fitted regression object by the best parameters (basis elements for data and beta and lambda penalty). 
#' \item \code{basis.x.opt}: Basis used for functional explanatory data estimation \code{fdata}. 
#' \item \code{basis.b.opt}: Basis used for functional \code{beta} parameter estimation. 
#' \item \code{lambda.opt}: \code{lambda} value that minimizes CV or GCV method. 
#' \item \code{gcv.opt}: Minimum value of CV or GCV method. 
#' }
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' @seealso See Also as: \code{\link{fregre.basis}},
#' \code{\link{summary.fregre.fd}} and \code{\link{predict.fregre.fd}} .\cr
#' Alternative method: \code{\link{fregre.pc.cv}} and
#' \code{\link{fregre.np.cv}}.
#' @references Ramsay, James O. and Silverman, Bernard W. (2006),
#' \emph{Functional Data Analysis}, 2nd ed., Springer, New York.
#' 
#' Febrero-Bande, M., Oviedo de la Fuente, M. (2012).  \emph{Statistical
#' Computing in Functional Data Analysis: The R Package fda.usc.} Journal of
#' Statistical Software, 51(4), 1-28. \url{https://www.jstatsoft.org/v51/i04/}
#' @keywords regression
#' @examples 
 #' \dontrun{
#' data(tecator)
#' x<-tecator$absorp.fdata[1:129]
#' y=tecator$y$Fat[1:129]
#' b1<-c(15,21,31)
#' b2<-c(7,9)
#' res1=fregre.basis.cv(x,y,basis.x=b1)
#' res2=fregre.basis.cv(x,y,basis.x=b1,basis.b=b2)
#' res1$gcv
#' res2$gcv
#' l=2^(-4:10)
#' res3=fregre.basis.cv(x,y,basis.b=b1,type.basis="fourier",
#' lambda=l,type.CV=GCV.S,par.CV=list(trim=0.15))
#' res3$gcv
#' }
#' 
#' @export
fregre.basis.cv <- function(fdataobj,y,basis.x=NULL,basis.b=NULL,
                           type.basis=NULL,lambda=0,Lfdobj=vec2Lfd(c(0,0),rtt),
                           type.CV=GCV.S,par.CV=list(trim=0),weights= rep(1,n),
                           verbose=FALSE,...){
  call<-match.call()
  if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
  #        omit <- omit.fdata(fdataobj, y)
  #        fdataobj <- omit[[1]]
  #        y <- omit[[2]]
  tol<-sqrt(.Machine$double.eps)
  x<-fdataobj[["data"]]
  tt<-fdataobj[["argvals"]]
  rtt<-fdataobj[["rangeval"]]
  n = nrow(x)
  np <- ncol(x)
  W<-diag(weights)
  if (n != (length(y))) stop("ERROR IN THE DATA DIMENSIONS")
  if (is.null(rownames(x)))        rownames(x) <- 1:n
  if (is.null(colnames(x)))        colnames(x) <- 1:np
  if (is.matrix(y)) y=as.vector(y)
  fou<-FALSE
  if (is.null(basis.x))  {
    #nbasis1=seq(max(floor(np/10),11),max(floor(np/5),11),by=2)
    nbasis1=seq(5,max(floor(np/5),11),by=2)
    lenbasis.x=length(nbasis1)
    basis.x=list()
    for (nb.x in 1:lenbasis.x) {
    if (!is.null(type.basis))  {
      aa1 <- paste("create.", type.basis[1], ".basis", sep = "")
      as <- list()
      as[[1]] <- range(tt)
      names(as)[[1]] <- "rangeval"
      as[[2]] <- nbasis1[nb.x]
      names(as)[[2]] <- "nbasis"
      if (verbose) basis.x[[nb.x]]=do.call(aa1, as)
      else  basis.x[[nb.x]]=suppressWarnings(do.call(aa1,as))
      }
    else      basis.x[[nb.x]]=create.bspline.basis(rangeval=rtt,nbasis=nbasis1[nb.x],...)
     }
  }
  else nbasis1<-basis.x
  if (is.null(basis.b))  {
  basis.b<-basis.x
  nbasis2<-nbasis1
  fou<-TRUE
   }
  else nbasis2<-basis.b
  lenlambda=length(lambda)
  a1=list()  ;a2=list()
  if (!is.null(type.basis)){
     if (type.basis=="fourier") fou<-TRUE}
  if (!is.list(basis.x))  {
   lenbasis.x=length(basis.x)
   for (nb.x in 1:lenbasis.x) {
     if (!is.null(type.basis))  {
       aa1 <- paste("create.", type.basis[1], ".basis", sep = "")
       as <- list()
       as[[1]] <- rtt
       names(as)[[1]] <- "rangeval"
       as[[2]] <- basis.x[nb.x]
       names(as)[[2]] <- "nbasis"
       if (verbose)   a1[[nb.x]]=do.call(aa1, as)
       else a1[[nb.x]]=suppressWarnings(do.call(aa1, as))
    }
   else  a1[[nb.x]]=create.bspline.basis(rangeval=rtt,nbasis=basis.x[[nb.x]])
   }
   basis.x=a1
  }
  else   lenbasis.x=length(nbasis1)
  if (!is.list(basis.b))  {
    lenbasis.y=length(basis.b)
    maxbasis.y<-which.max(basis.b)
    for (nb.y in 1:lenbasis.y) {
      if (!is.null(type.basis))  {
        if (length(type.basis)>1)       aa1 <- paste("create.", type.basis[2], ".basis", sep = "")
        else     aa1 <- paste("create.", type.basis[1], ".basis", sep = "")
        as <- list()
        as[[1]] <- rtt
        names(as)[[1]] <- "rangeval"
        as[[2]] <- basis.b[nb.y]
        names(as)[[2]] <- "nbasis"
        
      if (verbose) a2[[nb.y]]=do.call(aa1, as)
      else  a2[[nb.y]]=suppressWarnings(do.call(aa1, as))
        
        }
    else  a2[[nb.y]]=create.bspline.basis(rangeval=rtt,nbasis=basis.b[[nb.y]])
   }
   basis.b=a2
  }
  else  {
        lenbasis.y=length(nbasis2)
        maxbasis.y<-which.max(nbasis2)
  #      nb.y<-
        }
   pr=Inf
   i.lambda.opt=1;i.nb.y.opt=1;i.nb.x.opt=1
   xx<-fdata.cen(fdataobj)
   xmean=xx[[2]]
   xcen=xx[[1]]
   ymean=mean(y)
   ycen=y-ymean
   if (fou)   	 {
     basis.x<-basis.b
     nbasis1<-nbasis2
     if (verbose)   warning("Same number of basis elements in the basis.x and basis.b")
           lenbasis.x<-lenbasis.y
     nbasis12<-rep(NA,lenbasis.x)
     nbasis22<-rep(NA,lenbasis.y)
     x.fdfou<-Cfou<-Cmfou<- list()
     for (nb.x in 1:lenbasis.x) {
           xaux <- fdata2basis(fdataobj,basis.x[[nb.x]])
           # poner eso en fdata2basis!!!!!!!!!!!!!!!!!
           colnames(xaux$coefs) <- paste(fdataobj$names$main,".",colnames(xaux$coefs),sep="")
           # x.fd <- x.fdfou[[nb.x]] <- Data2fd(argvals=tt,y=t(xcen$data),basisobj=basis.x[[nb.x]])         
           Cfou[[nb.x]] <- xaux$coefs
           # Cmfou[[nb.x]] <- matrix(t(mean.fd(x.fdfou[[nb.x]])$coefs))
           nbasis12[nb.x] <- basis.x[[nb.x]]$type
           # x.fd <- x.fdfou[[nb.x]] <- Data2fd(argvals=tt,y=t(xcen$data),basisobj=basis.x[[nb.x]])         
           # Cfou[[nb.x]]=t(x.fdfou[[nb.x]]$coefs)
           # Cmfou[[nb.x]]=matrix(t(mean.fd(x.fdfou[[nb.x]])$coefs))
           # nbasis12[nb.x]<-basis.x[[nb.x]]$type
           }
    }
    else {
     nbasis12<-rep(NA,lenbasis.x)
     nbasis22<-rep(NA,lenbasis.y)
     x.fdfou<-Cfou<-Cmfou<- list()
     for (nb.x in 1:lenbasis.x) {
           # x.fd=x.fdfou[[nb.x]] =Data2fd(argvals=tt,y=t(xcen$data),basisobj=basis.x[[nb.x]])         
            xaux <- fdata2basis(fdataobj,basis.x[[nb.x]])
            colnames(xaux$coefs) <- paste(fdataobj$names$main,".",colnames(xaux$coefs),sep="")
            Cfou[[nb.x]] <- xaux$coefs
           # Cfou[[nb.x]]=t(x.fdfou[[nb.x]]$coefs)
           # Cmfou[[nb.x]]=matrix(t(mean.fd(x.fdfou[[nb.x]])$coefs))
           nbasis12[nb.x]<-basis.x[[nb.x]]$type
           }
    }
    gcv=array(NA,dim=c(lenbasis.x,lenbasis.y,lenlambda))
    for (nb.x in 1:lenbasis.x) {
     if (fou) {ifou<-nb.x;iifou<-nb.x}
     else {ifou<-1;iifou<-lenbasis.y  }
     for (nb.y in ifou:iifou) {
       nbasis22[nb.y]<-basis.b[[nb.y]]$type
       C<-Cfou[[nb.x]]
       # Cm<-Cmfou[[nb.x]] # viejo
       J<-inprod(basis.x[[nb.x]],basis.b[[nb.y]])        
  	   Z=C%*%J
       # Z=cbind(rep(1,len=n),Z)
  	   if (any(lambda!=0)) {
         # R=diag(0,ncol= basis.b[[nb.y]]$nbasis+1,nrow=basis.b[[nb.y]]$nbasis+1)
         # R[-1,-1]<-eval.penalty(basis.b[[nb.y]],Lfdobj)
  	     R=diag(0,ncol= basis.b[[nb.y]]$nbasis,nrow=basis.b[[nb.y]]$nbasis)
         }
       else R=0
       for (k in 1:lenlambda) {
      
# norm(R)/norm(t(Z)%*%W%*%Z
         Sb=t(Z)%*%W%*%Z+lambda[k]*R
        #Cinv<-solve(Sb)
        Cinv<-Minverse(Sb)
         Sb2=Cinv%*%t(Z)%*%W
         par.CV$S <-Z%*%Sb2
         par.CV$y<-y
         gcv[nb.x,nb.y,k]<- do.call(type.CV,par.CV)
         if ((gcv[nb.x,nb.y,k]+tol) < pr) {
            pr=gcv[nb.x,nb.y,k]
            lambda.opt=lambda[k]
            basis.b.opt=basis.b[[nb.y]]
            basis.x.opt=basis.x[[nb.x]]
            Sb.opt=Sb2
            Z.opt=Z
           # Cm.opt=Cm
            J.opt=J
            Cinv.opt=Cinv
            R.opt=R
          } 
         }
      }  }
  if (all(is.na(gcv))) stop("System is computationally singular. Try to reduce the number of basis elements")
         l = which.min(gcv)[1]
         
      gcv.opt=min(gcv,na.rm=TRUE)
      
      S=Z.opt%*%Sb.opt
      DD<-t(Z.opt)%*%y
      yp=S%*%y
      b.est=Sb.opt%*%y
      bet<-Cinv.opt%*%DD
      rownames(b.est)<-1:nrow(b.est)
      rownames(b.est)[1]<- "(Intercept)"
       # codigo viejo
      #beta.est2=fd(b.est2[-1,1]*diff(rtt),basis.b)
      beta.est=fd(b.est[,1],basis.b.opt) # nuevo
      
      beta.est=fd(b.est[,1],basis.b.opt)
      a.est=b.est[1,1]
      e=drop(y)-drop(yp)
      names(e)<-rownames(x)
      df=sum(diag(S))
#      df=basis.b.opt$nbasis+1
      sr2=sum(e^2)/(n-df)
      Vp<-sr2*Cinv.opt    
      r2=1-sum(e^2)/sum(ycen^2)
      object.lm=list()
      object.lm$coefficients<-drop(b.est)
      object.lm$residuals<-drop(e)
      object.lm$fitted.values<-yp
      object.lm$y<-y
      object.lm$x<-Z.opt
      object.lm$rank <- df
      object.lm$df.residual<-n-df
      vfunc=call[[2]]
      colnames(Z.opt) <-  paste(vfunc,".",basis.b.opt$names, sep = "")
      #colnames(Z.opt)[1]="(Intercept)"
      vcov2=sr2*Cinv.opt
      std.error=sqrt(diag(vcov2))
      t.value=b.est/std.error
      p.value= 2 * pt(abs(t.value),n-df, lower.tail = FALSE)
      coefficients<-cbind(b.est,std.error,t.value,p.value)
      colnames(coefficients) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
      rownames(coefficients)<-colnames(Z.opt)
      class(object.lm)<-"lm"
  b.est=b.est[-1]
  names(b.est)<-rownames(coefficients)[-1]
  lambda2<-paste("lambda=",lambda,sep="")
  #if (lenlambda==1)  gcv<-gcv[,,1]
  if (fou)  {
   nbasis12<-paste(nbasis12,nbasis2,sep="")
   gcv<-data.matrix(apply(gcv,3,diag))
   rownames(gcv)<-nbasis12
   colnames(gcv)<-lambda2
      }
  else{
  nbasis12<-paste(nbasis12,nbasis1,sep="")
  nbasis22<-paste(nbasis22,nbasis2,sep="")
  dimnames(gcv)<-list(nbasis12,nbasis22,lambda2)
  }
  # x.fd=Data2fd(argvals=tt,y=t(xcen$data),basisobj=basis.x.opt)          
  model<-fregre.basis(fdataobj = fdataobj, y =y, basis.x=basis.x.opt,
                    basis.b= basis.b.opt, lambda=lambda.opt, 
                    Lfdobj= Lfdobj, weights=weights,...)
  model$coefficients <- coefficients
  class(model) <- "fregre.fd"
  out <- list("fregre.basis"=model,    
            basis.x.opt=basis.x.opt$nbasis,
            basis.b.opt=basis.b.opt$nbasis,
           "lambda.opt"=lambda.opt,
           "gcv.opt"=gcv.opt,"gcv"=gcv)
  return(invisible(out))
}
