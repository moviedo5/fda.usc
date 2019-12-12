#' @rdname  optim.basis
#' @aliases optim.basis min.basis
#' @title Select the number of basis using GCV method.
#' @note min.basis deprecated.
#' @description Functional data estimation via basis representation using cross-validation
#' (CV) or generalized cross-validation (GCV) method with a roughness penalty.
#' 
#' @details Provides the least GCV for functional data for a list of number of basis
#' \code{num.basis} and lambda values \code{lambda}. You can define the type of
#' CV to use with the \code{type.CV}, the default is used \code{GCV.S}. 
#' 
#' Smoothing matrix is performed by \code{\link{S.basis}}. \code{W} is the
#' matrix of weights of the discretization points.
#' 
#' @param fdataobj \code{\link{fdata}} class object.
#' @param type.CV Type of cross-validation. By default generalized
#' cross-validation (GCV) method.
#' @param W Matrix of weights.
#' @param lambda A roughness penalty. By default, no penalty \code{lambda=0}.
#' @param numbasis Number of basis to use.
#' @param type.basis Character string which determines type of basis. By
#' default \emph{"bspline"}.
#' @param par.CV List of parameters for type.CV: trim, the alpha of the
#' trimming and \code{draw=TRUE}.
#' @param verbose If \code{TRUE} information about GCV values and input
#' parameters is printed. Default is \code{FALSE}.
#' @param \dots Further arguments passed to or from other methods. Arguments to
#' be passed by default to \link[fda]{create.basis}.
#' @return
#' \itemize{
#' \item \code{gcv}{ Returns GCV values calculated for input parameters.}
#' \item \code{fdataobj}{ Matrix of set cases with dimension (\code{n} x \code{m}),
#' where \code{n} is the number of curves and \code{m} are the points observed
#' in each curve.} 
#' \item \code{fdata.est}{ Estimated \code{fdata} class object.}
#' \item \code{numbasis.opt}{ \code{numbasis} value that minimizes CV or GCV method.}
#' \item \code{lambda.opt}{ \code{lambda} value that minimizes CV or GCV method.}
#' \item \code{basis.opt}{ \code{basis} for the minimum CV or GCV method.}
#' \item \code{S.opt}{ Smoothing matrix for the minimum CV or GCV method.}
#' \item \code{gcv.opt}{ Minimum of CV or GCV method.} 
#' \item \code{lambda}{ A roughness penalty. By default, no penalty \code{lambda=0}.}
#' \item \code{numbasis}{ Number of basis to use.} 
#' \item \code{verbose}{ If \code{TRUE} information about GCV values
#' and input parameters is printed. Default is \code{FALSE}.}
#'  }
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@usc.es}
#' @seealso See Also as \code{\link{S.basis}}. \cr Alternative method:
#' \code{\link{optim.np}}
#' @references Ramsay, James O., and Silverman, Bernard W. (2006), \emph{
#' Functional Data Analysis}, 2nd ed., Springer, New York.
#' 
#' Wasserman, L. \emph{All of Nonparametric Statistics}. Springer Texts in
#' Statistics, 2006.
#' 
#' Hardle, W. \emph{Applied Nonparametric Regression}. Cambridge University
#' Press, 1994.
#' 
#' Febrero-Bande, M., Oviedo de la Fuente, M. (2012).  \emph{Statistical
#' Computing in Functional Data Analysis: The R Package fda.usc.} Journal of
#' Statistical Software, 51(4), 1-28. \url{http://www.jstatsoft.org/v51/i04/}
#' @keywords nonparametric
#' @examples
#' \dontrun{ 
#' a1<-seq(0,1,by=.01)
#' a2=rnorm(length(a1),sd=0.2)
#' f1<-(sin(2*pi*a1))+rnorm(length(a1),sd=0.2)
#' nc<-50
#' np<-length(f1)
#' tt=1:101
#' S<-S.NW(tt,2)
#' mdata<-matrix(NA,ncol=np,nrow=50)
#' for (i in 1:50) mdata[i,]<- (sin(2*pi*a1))+rnorm(length(a1),sd=0.2)
#' mdata<-fdata(mdata)
#' nb<-floor(seq(5,29,len=5))
#' l<-2^(-5:15)
#' out<-optim.basis(mdata,lambda=l,numbasis=nb,type.basis="fourier")
#' matplot(t(out$gcv),type="l",main="GCV with fourier basis")
#' 
#' # out1<-optim.basis(mdata,type.CV = CV.S,lambda=l,numbasis=nb)
#' # out2<-optim.basis(mdata,lambda=l,numbasis=nb)
#' 
#' # variance calculations
#' y<-mdata
#' i<-3
#' z=qnorm(0.025/np)
#' fdata.est<-out$fdata.est
#' var.e<-Var.e(mdata,out$S.opt)
#' var.y<-Var.y(mdata,out$S.opt)
#' var.y2<-Var.y(mdata,out$S.opt,var.e)
#' 
#' # estimated fdata and point confidence interval
#' upper.var.e<-out$fdata.est[["data"]][i,]-z*sqrt(diag(var.e))
#' lower.var.e<-out$fdata.est[["data"]][i,]+z*sqrt(diag(var.e))
#' dev.new()
#' plot(y[i,],lwd=1,ylim=c(min(lower.var.e),max(upper.var.e)))
#' lines(out$fdata.est[["data"]][i,],col=gray(.1),lwd=1)
#' lines(out$fdata.est[["data"]][i,]+z*sqrt(diag(var.y)),col=gray(0.7),lwd=2)
#' lines(out$fdata.est[["data"]][i,]-z*sqrt(diag(var.y)),col=gray(0.7),lwd=2)
#' lines(upper.var.e,col=gray(.3),lwd=2,lty=2)
#' lines(lower.var.e,col=gray(.3),lwd=2,lty=2)
#' legend("top",legend=c("Var.y","Var.error"), col = c(gray(0.7),
#' gray(0.3)),lty=c(1,2))
#' }
#' 
#' @export 
optim.basis<-function(fdataobj,type.CV=GCV.S,W=NULL,lambda=0,
numbasis=floor(seq(ncol(fdataobj)/16,ncol(fdataobj)/2,len=10)),
type.basis="bspline",par.CV=list(trim=0,draw=FALSE), verbose = FALSE,...){
 if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
  nas1<-is.na.fdata(fdataobj)
 if (any(nas1))  stop("fdataobj contain ",sum(nas1)," curves with some NA value \n")

x<-fdataobj[["data"]]
tt<-fdataobj[["argvals"]]
rtt<-fdataobj[["rangeval"]]
nam<-fdataobj[["nam"]]
   lenlambda<-length(lambda)
   lenbasis<-length(numbasis)
   nc<-nrow(fdataobj)
   np<-ncol(fdataobj)
   gcv<-array(Inf,dim=c(lenbasis,lenlambda))
   GCV.basis.min=Inf 
   as <- list()
   as[[1]] <- rtt
   names(as)[[1]]<-'rangeval'
   C <- match.call()
   mf <- match.call(expand.dots = FALSE)
   m <- match(c("fdataobj","tt","type.CV","W","lambda","numbasis","type.basis","verbose"),
   names(mf),0L)
   imetric <- m[7]
   if (imetric == 0) {
        a1 <- create.bspline.basis
        len.metricc <- length(formals(a1))
        vv <- array(0, dim = c(len.metricc))     }
   else {
         a1 <- paste('create.',type.basis,'.basis',sep="")
        len.metricc <- length(formals(a1))
        vv <- array(0, dim = c(len.metricc))     }
   ii <- imetric + 1
   if (C[ii] != "NULL()") {
        ind.m <- 3
        while (C[ii] != "NULL()" && ind.m <= len.metricc) {
            aa <- any(names(C) == names(formals(a1))[ind.m])
            if (aa) {
                vv[ind.m] <- which(names(C) == names(formals(a1)[ind.m]))
                ii <- ii + 1
                as[[ind.m]] <- C[[vv[ind.m]]]
                names(as)[[ind.m]]<-names(formals(a1)[ind.m])           }
            else {             as[[ind.m]] <- formals(a1)[[ind.m]]    }
            ind.m <- ind.m + 1
        }
    }
    for (i in 1:lenbasis) {
      as[[2]] <- numbasis[i]
      names(as)[[2]]<-'nbasis'
      base <- do.call(a1, as)
      for (k in 1:lenlambda) {
          S2<-S.basis(tt,base,lambda[k])
          gcv[i,k]<-type.CV(fdataobj,S=S2,W=W,trim=par.CV$trim,draw=par.CV$draw,...) }
    }
l=which.min(gcv)
i= (l %% lenbasis)
k= (l %/% lenbasis)+1
if (i==0) {i=lenbasis;k=k-1}
lambda.opt<-lambda[k]
numbasis.opt<-numbasis[i]
gcv.opt<-gcv[l]
as[[2]]=numbasis[i]
names(as)[[2]]<-'nbasis'
base.opt=do.call(a1,as)
S.opt<-S.basis(tt,base.opt,lambda[k])
fdata.est<-S.opt%*%t(x) ### 
if (length(numbasis)>1) dimnames(gcv)[[1]]<-numbasis
if (length(lambda)>1) dimnames(gcv)[[2]]<-lambda
if (verbose) {
  cat("\n The minimum GCV (GCV.OPT=",round(gcv.opt,4),sep="",") is achieved with
 the number of basis (numbasis.opt=",numbasis.opt,")\n and lambda value    (lambda.opt=",lambda.opt,")\n\n")
if (lenbasis>1) {
  if (numbasis.opt==min(numbasis))  cat(" Warning: numbasis.opt is the minimum number of basis provided, range(numbasis)=",range(numbasis),"\n")
  else if (numbasis.opt==max(numbasis)) cat(" Warning: numbasis.opt is the maximum number of basis provided, range(numbasis)=",range(numbasis),"\n")
}
if (lenlambda > 1) {
  if (lambda.opt==min(lambda))  cat(" Warning: lambda.opt is the minimum lambda value provided, range(lambda)=",range(lambda),"\n")
  else   if (lambda.opt==max(lambda))  cat(" Warning: lambda.opt is the maximum lambda value provided, range(lambda)=",range(lambda),"\n")
}
}
fdata.est=fdata(t(fdata.est),tt,rtt,nam)
output<-list(gcv=gcv,numbasis=numbasis,lambda=lambda,fdataobj=fdataobj,fdata.est=fdata.est,gcv.opt=gcv.opt,numbasis.opt=numbasis.opt,lambda.opt=lambda.opt,S.opt=S.opt,base.opt=base.opt)
}

