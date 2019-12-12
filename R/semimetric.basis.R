#' Proximities between functional data
#' 
#' Approximates semi-metric distances for functional data of class \code{fdata}
#' or \code{fd}.
#' 
#' @details Approximates semi-metric distances for functional data of two \code{fd}
#' class objects.  If functional data are not functional \code{fd} class, the
#' \code{semimetric.basis} function creates a basis to represent the functional
#' data, by default is used \code{\link{create.bspline.basis}} and the
#' \code{fdata} class object is converted to \code{fd} class using the
#' \code{\link{Data2fd}}.\cr The function calculates distances between the
#' derivative of order \code{nderiv} of curves using \code{\link{deriv.fd}}
#' function.
#' 
#' @param fdata1 Functional data 1 or curve 1.
#' @param fdata2 Functional data 2 or curve 2.
#' @param nderiv Order of derivation, used in \code{deriv.fd}
#' @param type.basis1 Type of Basis for \code{fdata1}.
#' @param nbasis1 Number of Basis for \code{fdata1}.
#' @param type.basis2 Type of Basis for \code{fdata2}.
#' @param nbasis2 Number of Basis for \code{fdata2.}
#' @param \dots Further arguments passed to or from other methods.
#' @return Returns a proximities matrix between functional data.
#' @seealso See also \code{\link{metric.lp}}, \code{\link{semimetric.NPFDA}}
#' and \code{\link{deriv.fd}}
#' @references Ferraty, F. and Vieu, P. (2006). \emph{Nonparametric functional
#' data analysis.} Springer Series in Statistics, New York.
#' @keywords cluster
#' @examples
#' \dontrun{
#' data(phoneme)
#' DATA1<-phoneme$learn[c(30:50,210:230)]
#' DATA2<-phoneme$test[231:250]
#' a1=semimetric.basis(DATA1,DATA2)
#' a2=semimetric.basis(DATA1,DATA2,type.basis1="fourier",
#' nbasis1=11, type.basis2="fourier",nbasis2=11)
#' fd1 <- fdata2fd(DATA1)
#' fd2 <- fdata2fd(DATA2)
#' a3=semimetric.basis(fd1,fd2)
#' a4=semimetric.basis(fd1,fd2,nderiv=1)
#' }
#' 
#' @export
semimetric.basis=function(fdata1, fdata2 = fdata1,nderiv=0,type.basis1=NULL,
nbasis1=NULL,type.basis2=type.basis1,nbasis2=NULL,...) {
 if (any(class(fdata1)=="fd")) {
   r=fdata1$basis[[3]]
   tt=seq(r[1],r[2],len=length(fdata1$fdnames$time))
   df1=deriv.fd(fdata1,nderiv)
   df2=deriv.fd(fdata2,nderiv)
   fd1=fdata(t(eval.fd(tt,df1)),tt,r)
   fd2=fdata(t(eval.fd(tt,df2)),tt,r)
   mdist=metric.lp(fd1,fd2,...)
  }
 else {
 if (!is.fdata(fdata1)) fdata1<-fdata(fdata1)
 tt<-fdata1[["argvals"]]
 rtt<-fdata1[["rangeval"]]
 nas1<-is.na.fdata(fdata1)
 if (any(nas1))  stop("fdata1 contain ",sum(nas1)," curves with some NA value \n")
 else  if (!is.fdata(fdata2))  {fdata2<-fdata(fdata2,tt,rtt) }
 nas2<-is.na.fdata(fdata2)
 if (any(nas2))  stop("fdata2 contain ",sum(nas2)," curves with some NA value \n")
 rtt<-fdata1[["rangeval"]]
#   print("Raw class")
   np=ncol(fdata1)
   if (is.null(nbasis1)) {
       nbasis1=ifelse(floor(np/3) > floor((np - nderiv - 4)/2),
       floor((np - nderiv - 4)/2), floor(np/3))
       }
   if (is.null(nbasis2)) nbasis2=nbasis1
   as <- list()
   bs <- list()
   as[[1]] <- rtt
   bs[[1]] <- rtt
   names(as)[[1]]<-"rangeval"
   names(bs)[[1]]<-"rangeval"
   as[[2]] <- nbasis1
   names(as)[[2]]<-"nbasis"
   C <- match.call()
   mf <- match.call(expand.dots = FALSE)
   m<-match(c("fdata1", "fdata2","nderiv","type.basis1","nbasis1","type.basis2","nbasis2"),names(mf),0L)
   imetric <- m[4]
   imetric2 <- m[6]
   if (imetric == 0) {
        type.basis1="bspline"
        a1 <- create.bspline.basis
        len.metric <- length(formals(a1))
        vv <- array(0, dim = c(len.metric))    }
   else {  a1 <- paste('create.',type.basis1,'.basis',sep="")
        len.metric <- length(formals(a1))
        vv <- array(0, dim = c(len.metric)) }
  ii <- imetric + 1
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
 if (imetric2 == 0) {
         b1 <-  paste('create.',type.basis1,'.basis',sep="")
        len.metric <- length(formals(b1))
        vv <- array(0, dim = c(len.metric))
}
else {  b1 <- paste('create.',type.basis2,'.basis',sep="")
        len.metric <- length(formals(b1))
        vv <- array(0, dim = c(len.metric)) }
 ii <- imetric2 + 1
 if (C[ii] != "NULL()") {
        ind.m <- 3
        while (C[ii] != "NULL()" && ind.m <= len.metric) {
            aa <- any(names(C) == names(formals(b1))[ind.m])
            if (aa) {
                vv[ind.m] <- which(names(C) == names(formals(b1)[ind.m]))
                ii <- ii + 1
                bs[[ind.m]] <- C[[vv[ind.m]]]
                names(bs)[[ind.m]]<-names(formals(b1)[ind.m])            }
            else { as[[ind.m]] <- formals(b1)[[ind.m]]}
            ind.m <- ind.m + 1}
  }

   bs[[2]] <- nbasis2
   names(bs)[[2]]<-'nbasis'
   b1.2<- do.call(b1, bs)
   fd1.1 <- Data2fd(argvals=tt,y=t(fdata1$data),basisobj=b1.1)
   fd1.2 <- Data2fd(argvals=tt,y=t(fdata2$data),basisobj=b1.2)
   df1=deriv.fd(fd1.1,nderiv)
   df2=deriv.fd(fd1.2,nderiv)
   fd1=fdata(t(eval.fd(tt,df1)),tt)
   fd2=fdata(t(eval.fd(tt,df2)),tt)
   mdist=metric.lp(fd1,fd2,...)
}
attr(mdist,"call")<-"semimetric.basis"
attr(mdist,"par.metric")<-list("nderiv"=nderiv,"type.basis1"=type.basis1,
"nbasis1"=nbasis1,"type.basis2"=type.basis2,"nbasis2"=nbasis2)
mdist
}


