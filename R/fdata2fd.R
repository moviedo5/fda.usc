#' Converts fdata class object into fd class object
#' 
#' Converts \code{fdata} class object into \code{fd} class object using
#' \code{Data2fd} function.
#' 
#' @param fdataobj \code{\link{fdata}} class object.
#' @param type.basis Type of basis. A function \code{create."type.basis".basis}
#' must exists. By default, \code{bspline} basis is used.
#' @param nbasis Number of basis which is used in \code{create.basis} function.
#' @param nderiv Order of derivation which is used in \code{deriv.fd} function
#' (optional).
#' @param lambda Weight on the smoothing operator specified by \code{nderiv}.
#' @param \dots Further arguments passed to or from other methods.
#' @return Return an object of the \code{fd} class.
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' @seealso See Also as \code{\link{fdata}} and \link[fda]{Data2fd}
#' @references Febrero-Bande, M., Oviedo de la Fuente, M. (2012).
#' \emph{Statistical Computing in Functional Data Analysis: The R Package
#' fda.usc.} Journal of Statistical Software, 51(4), 1-28.
#' \url{https://www.jstatsoft.org/v51/i04/}
#' @keywords manip
#' @examples
#' \dontrun{
#' data(phoneme)
#' mlearn<-phoneme$learn[1]
#' fdata2=fdata2fd(mlearn)
#' class(mlearn)
#' class(fdata2)
#' fdata3=fdata2fd(mlearn,type.basis="fourier",nbasis=7)
#' plot(mlearn)
#' lines(fdata2,col=2)
#' lines(fdata3,col=3)
#' fdata5=fdata2fd(mlearn,nderiv=1)
#' }
#' @export
fdata2fd=function(fdataobj,type.basis=NULL,nbasis=NULL,nderiv=0,lambda=NULL,...) {
if (is.fdata(fdataobj)) DATA=fdataobj[["data"]]
else stop("No fdata object")
np=ncol(DATA)
tt =fdataobj[["argvals"]]
rtt<-fdataobj[["rangeval"]]
if (is.null(lambda)) lambda=3e-08/diff(rtt)
if (is.null(nbasis)) {
       nbasis=ifelse(floor(np/3) > floor((np - nderiv - 4)/2),
       floor((np - nderiv - 4)/2), floor(np/3))
       }
   as <- list()
   as[[1]] <-rtt
   names(as)[[1]]<-"rangeval"
   as[[2]] <- nbasis
   names(as)[[2]]<-"nbasis"
   C <- match.call()
   mf <- match.call(expand.dots = FALSE)
   m<-match(c("DATA","type.basis","nbasis","nderiv"),names(mf),0L)
   imetri <- m[2]
   if (imetri == 0) {
        type.basis1="bspline"
        a1 <- create.bspline.basis
        len.metric <- length(formals(a1))
        vv <- array(0, dim = c(len.metric))    }
   else {  a1 <- paste('create.',type.basis,'.basis',sep="")
        len.metric <- length(formals(a1))
        vv <- array(0, dim = c(len.metric)) }
  ii <- imetri + 1
  if (C[ii] != "NULL()") {
        ind.m <- 3
        while (C[ii] != "NULL()" && ind.m <= len.metric) {
            aa <- any(names(C) == names(formals(a1))[ind.m])
            if (aa) {
                vv[ind.m] <- which(names(C) == names(formals(a1)[ind.m]))
                ii <- ii + 1
                as[[ind.m]] <- C[[vv[ind.m]]]
                names(as)[[ind.m]]<-names(formals(a1)[ind.m])            }
#            else {                 as[[ind.m]] <- formals(a1)[[ind.m]]   }
            ind.m <- ind.m + 1            }
  }
 b1.1<- do.call(a1, as)
 class(DATA)="matrix"

 fd1.1 <- Data2fd(argvals=tt,y=t(DATA),basisobj=b1.1,lambda=lambda,...) ######
 if (nderiv>0) fd1.1=deriv.fd(fd1.1,int2Lfd(nderiv)) #######
 fd1.1
}

