#' Create Basis Set for Functional Data of fdata class
#' 
#' @description Compute basis for functional data.
#' 
#' @aliases create.fdata.basis create.pc.basis create.pls.basis
#' create.raw.fdata 
#' @param fdataobj \code{\link{fdata}} class object.
#' @param y Vector of response (scalar).
#' @param l Vector of basis index.
#' @param maxl maximum number of basis
#' @param type.basis Type of basis (see create.basis function).
#' @param rangeval A vector of length 2 giving the lower and upper limits of
#' the range of permissible values for the function argument.
#' @param norm If \code{TRUE} the norm of eigenvectors \code{basis} is 1.
#' @param class.out =="fd" basisfd class, =="fdata" fdata class.
#' @param basis "fd" basis object.
#' @param lambda Amount of penalization. Default value is 0, i.e. no
#' penalization is used.
#' @param P If P is a vector: coefficients to define the penalty matrix object.
#' By default P=c(0,0,1) penalize the second derivative (curvature) or
#' acceleration.  If P is a matrix: the penalty matrix object.
#' @param \dots Further arguments passed to or from other methods.
#' @return 
#' \itemize{
#' \item \code{basis}{ basis} 
#' \item \code{x}{ if \code{TRUE} the value of the rotated data (the centred data multiplied by the rotation matrix) is returned}
#' \item \code{mean}{ functional mean of \code{fdataobj}} 
#' \item \code{df}{ degree of freedom} 
#' \item \code{type}{ type of basis}
#' }
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@usc.es}
#' @seealso See Also as \link[fda]{create.basis} and \code{\link{fdata2pc}}.
#' @references Ramsay, James O. and Silverman, Bernard W. (2006),
#' \emph{Functional Data Analysis}, 2nd ed., Springer, New York.
#' 
#' N. Kraemer, A.-L. Boulsteix, and G. Tutz (2008). Penalized Partial Least
#' Squares with Applications to B-Spline Transformations and Functional Data.
#' Chemometrics and Intelligent Laboratory Systems, 94, 60 - 69.
#' \url{http://dx.doi.org/10.1016/j.chemolab.2008.06.009}
#' @keywords multivariate
#' @examples
#' \dontrun{
#' data(tecator)
#' basis.pc<-create.pc.basis(tecator$absorp.fdata,c(1,4,5))
#' plot(basis.pc$basis,col=1)
#' basis.pls<-create.pls.basis(tecator$absorp.fdata,y=tecator$y[,1],c(1,4,5))
#' lines(basis.pls$basis,col=2)
#' 
#' basis.fd<-create.fdata.basis(tecator$absorp.fdata,c(1,4,5),
#' type.basis="fourier")
#' plot(basis.pc$basis)
#' basis.fdata<-create.fdata.basis(tecator$absorp.fdata,c(1,4,5),
#' type.basis="fourier",class.out="fdata")
#' plot(basis.fd,col=2,lty=1)
#' lines(basis.fdata,col=3,lty=1)
#' }

#' @export
create.fdata.basis <-function(fdataobj, l = 1:5, maxl = max(l), type.basis = "bspline", 
          rangeval = fdataobj$rangeval, class.out = "fd") 
{
  aa1 <- paste("create.", type.basis, ".basis", sep = "")
  if (type.basis == "pc") {
    as <- list()
    as$fdataobj <- fdataobj
    as$l <- l
    basis = do.call(aa1, as)
  }
  if (type.basis %in% c("bspline", "fourier", "constant", "exponential", 
                        "polygonal", "power")) {
    as <- list()
    as$rangeval <- rangeval
    as$nbasis <- maxl
    basis = do.call(aa1, as)
    #basis$params <- diff(rangeval)
    basis$dropind <- setdiff(1:maxl, l)
    if (class.out == "fdata") {
      nam <- basis$names[intersect(1:maxl, l)]
      basis = fdata(t(eval.basis(fdataobj$argvals, basis)), 
                    fdataobj$argvals, fdataobj$rangeval)
      rownames(basis$data) <- nam
      basis$type <- type.basis
      basis$nbasis <- maxl
      basis$dropind <- as$dropind
    }
  }
  basis
}


scores.basis <-function(fdataobj,l=1:5,maxl=max(l),type.basis="bspline", lambda ){
      aa1 <- paste("create.",type.basis,".basis", sep = "")
      if (type.basis=="pc"){
        
        as <- list() 
        as$fdataobj <- fdataobj
        as$l <- l
        if (missing(lambda)) lambda = 0
        as$lambda <- lambda
        basis=do.call(aa1,as)
      }
      if (type.basis %in% c("bspline","fourier","constant","exponential"
                            ,"polygonal","power")){
        
################
#        fdataobj<-tecator$absorp.fdata
#        l=1:5
#        type.basis="bspline"
#        lambda = NULL
################        
        # as <- list()
        maxl=max(l)
        if (missing(lambda)) lambda = NULL
        #as$lambda <- lambda
        #basis=do.call(aa1,as)
        basis0 = fdata2fd(fdataobj, type.basis, nbasis = maxl, lambda = lambda)
        basis<-list()
        basis$basis <- basis0$basis
        basis$basis.fdata=fdata(t(eval.basis(fdataobj$argvals,basis0$basis)),
                          fdataobj$argvals,fdataobj$rangeval)
        basis$x <- t(basis0$coefs)  
        basis$fdataobj <- fdataobj
        basis$l <- l
        basis$lambda <- as$lambda
        basis$type <- type.basis
        basis$nbasis <-maxl
        basis$dropind<-setdiff(1:maxl,l)
        basis$rangeval <- fdataobj$rangeval
        basis$mean <- func.mean(fdataobj)
        
        #x.fd = Data2fd(argvals = tt,
         #              y = t(fdata.cen(fdataobj,object$mean[[vfunc[i]]])[[1]]$data), 
          #             basisobj = object$basis.x[[vfunc[i]]]$basis, 
           #            fdnames = fdnames)
        
        
        #as$nbasis <-maxl
        #as$dropind<-setdiff(1:maxl,l)
        #as$rangeval <- fdataobj$rangeval
        #nam<-basis$names[intersect(1:maxl,l)]
        #rownames(basis$data)<-nam
        #basis$type<-type.basis
        #basis$nbasis<-maxl
        #basis$dropind<-as$dropind
      }
  basis
} 
#######################
     
#' @rdname create.fdata.basis
#' @export
create.pc.basis<-function(fdataobj,l=1:5,norm=TRUE,basis=NULL,
                          lambda=0,P = c(0, 0, 1),...){
 tt<-fdataobj$argvals
 rtt<-fdataobj$rangeval
 dropind=NULL
 if (lambda>0) pc<-fdata2pc(fdataobj,norm=norm,ncomp=max(l),lambda=lambda,P=P,...)
 else  pc<-fdata2pc(fdataobj,norm=norm,ncomp=max(l),...)
 vs<-pc$rotation$data    
 lenl<-length(l)  
 pc.fdata<-pc$u[,l,drop=FALSE]%*%(diag(lenl)*pc$d[l])%*%vs[l,,drop=FALSE]
 pc.fdata<-sweep(pc.fdata,2,matrix(pc$mean$data,ncol=1),"+")
 basis.pc = pc$rotation[l, ,drop=FALSE]
 rownames(basis.pc$data) <- paste("PC", l, sep = "")
 basisobj<-pc
 fdnames<- colnames(pc$x[,l,drop=FALSE])
 if (is.null(basis)) {
   pc.fdata<-fdata(pc.fdata,tt,rtt,fdataobj$names)
   out <- list(fdataobj.pc=pc.fdata,basis = basis.pc, x = pc$x, mean = pc$mean,
   fdataobj.cen = pc$fdataobj.cen,fdataobj = fdataobj,l = l,norm=norm,
   lambda=lambda,P=P,type = "pc")
   class(out) <- "fdata.comp"
   }
 else {
      fdobj<- Data2fd(argvals = tt, y = t(pc.fdata),basisobj = basis)
      out<-list()
      out$harmonics<-fdobj
      colnames(out$harmonics$coefs)<-rownames(fdataobj$data)
      out$values<-pc$newd^2
      out$scores<-pc$x[,l,drop=FALSE]
      rownames(out$scores)<-rownames(fdataobj$data)
      out$varprop<-out$values[l]/sum(out$values)
      out$meanfd<- Data2fd(argvals = tt, y = pc$mean$data[1,],basisobj = basis)
      class(out) <- "pca.fd"
      }
 return(out) 
}

#' @rdname create.fdata.basis
#' @export
create.pls.basis<-function(fdataobj,y,l=1:5,norm=TRUE,lambda=0,P = c(0, 0, 1),...){
if (lambda>0) pls<-fdata2pls(fdataobj,y,norm=norm,ncomp=max(l),lambda=lambda,P=P,...)
 else  pls<-fdata2pls(fdataobj,y,norm=norm,ncomp=max(l),...)
     basis=pls$rotation[l,,drop=FALSE]
     rownames(basis$data)<-paste("PLS",l,sep="")
out<-list("basis"=basis,"x"=pls$x,"mean"=pls$mean,"df"=pls$df,
"fdataobj.cen"=pls$fdataobj.cen,"fdataobj"=fdataobj,norm=norm,
"l"=l,"type"="pls","y"=y)
class(out) <- "fdata.comp"
return(out)
} 

#' @rdname create.fdata.basis
#' @export
create.raw.fdata=function (fdataobj, l = 1:ncol(fdataobj))
{
    return(list(basis =fdataobj[,l] , type = "raw"))
}

#########################
create.mfdata.basis <- function(mfdata, l = 1:5
                                , type.basis = "bspline"
                                , class.out = "fd") 
{
  nvar <- length(mfdata)
  nam <- names(mfdata)
  aa1 <- paste("create.", type.basis, ".basis", sep = "")
  basis.x <- NULL
  maxl <- max(l)
  if (type.basis == "pc") {
    for (i in 1:nvar) {
      as <- list()
      as$fdataobj <- mfdata[[nam[i]]]
      as$l <- l
      basis = do.call(aa1, as)
      basis.x[[nam[i]]] <- basis
  }}
  
  if (type.basis %in% c("bspline", "fourier", "constant", "exponential", 
                        "polygonal", "power")) {
    for (i in 1:nvar) {
    as <- list()
    fdataobj <- mfdata[[nam[i]]]
    rangeval = fdataobj$rangeval
    as$rangeval <- rangeval
    as$nbasis <- maxl
    as$dropind <- setdiff(1:maxl, l)
    basis = do.call(aa1, as)
    
    if (class.out == "fdata") {
      nam <- basis$names[intersect(1:maxl, l)]
      basis = fdata(t(eval.basis(fdataobj$argvals, basis)), 
                    fdataobj$argvals, fdataobj$rangeval)
      rownames(basis$data) <- nam
      basis$type <- type.basis
      basis$nbasis <- maxl
      basis$dropind <- as$dropind
    }
   basis.x[[nam[i]]] <- basis
   }
  }
  return(invisible(basis.x))
}

#########################
create.ldata.basis <- function(ldata, l = 1:5
                               , type.basis = "bspline"
                               , class.out = "fd") 
{
  clases <- sapply(ldata,class)
  ifdata <- which(clases == "fdata")
  basis.x <- create.mfdata.basis(ldata[ifdata], l = l, 
             type.basis = type.basis, class.out = "fd") 
  return(basis.x)
}  
