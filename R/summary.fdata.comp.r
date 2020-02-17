#' @title Correlation for functional data by Principal Component Analysis
#' 
#' @description Summary of functional principal components
#' 
#' @param object fdata.comp class object calculated by: \code{fdata2pc},
#' \code{fdata2pls}, \code{fregre.pc} or \code{fregre.pls}.
#' @param biplot =TRUE draw the biplot and PC (or PLS) components.
#' @param \dots Further arguments passed to or from other methods.
#' @return If \code{corplot}=TRUE, are displaying the biplot between the PC (or PLS) components.\cr
#'  If \code{ask}=TRUE, draw each graph in a window, waiting to confirm the change
#' of page with a click of the mouse or pressing ENTER.  If \code{ask}=FALSE draw graphs in one window.
#' @author Manuel Febrero-Bande and Manuel Oviedo de la Fuente \email{manuel.oviedo@@usc.es}
#' @seealso See Also as \code{\link{fdata2pc}}, \code{\link{fdata2pls}} and \link[stats]{cor}
#' @references Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988).
#' \emph{The New S Language}. Wadsworth \& Brooks/Cole.
#' 
#' Venables, W. N. and B. D. Ripley (2002). \emph{Modern Applied Statistics
#' with S}. Springer-Verlag.
#' @keywords multivariate
#' @examples
#' \dontrun{
#' n <- 200
#' tt <- seq(0,1,len=101)
#' x0 <- rproc2fdata(n,tt,sigma="wiener")
#' x1 <- rproc2fdata(n,tt,sigma=0.1)
#' x <- x0*3+x1
#' beta <- tt*sin(2*pi*tt)^2
#' fbeta <- fdata(beta,tt)
#' y <- inprod.fdata(x,fbeta) + rnorm(n,sd=0.1)
#' pc1 <- fdata2pc(x,3)
#' summary(pc1)
#' pls1 <- fdata2pls(x,y)
#' summary(pls1)
#' }
#' 
#' @export
summary.fdata.comp=function(object,biplot=TRUE,...) {
  if (inherits(object, "fdata.comp"))         {
     a1=TRUE
     pr.com<-object
     if (is.null(y)) {
        if (object$call[[1]]=="fdata2pls" | object$call[[1]]=="fdata2ppls")
         y<-object$y                  
     }
  } else if (inherits(object, "fregre.fd"))     {
     a1=FALSE
     pr.com<-object$fdata.comp
     y<-object$y
  } else stop("Error in input data")
 out2=pr.com$x
 l<-object$l
 le=length(l)
 rotation=aperm(pr.com$rotation$data)
 p=pr.com$x
 if (object$call[[1]]=="fdata2pls" | object$call[[1]]=="fdata2ppls") {
  var.1<-apply(p, 2, var)
  pr.x2= var.1/sum(var.1)
 } else {
 d<-object$d
 if (is.null(object$rn)) rn<-0
 pr.x2<-(d^2+rn)/sum(d^2+rn)
 names(pr.x2)<-colnames(out2)[l]
 }
 C<-match.call()
 lenC=length(C)
 cor.y.pc=rep(NA,le)
 xxx=cbind(y,pr.com$x)
 cor.y.pc=round(cor(xxx[,c(1,l+1)]),3)[1,-1]
 types<-colnames(pr.com$x)
 cat("\n     - SUMMARY:  ",object$call[[1]]," object   -\n")
 #if (object$call[[1]]=="fdata2pc" | object$call[[1]]=="fdata2ppc") 
   {
   cat("\n-With",le," components are explained ",round(sum(pr.x2[l])*100
 ,2),"%\n of the variability of explicative variables.\n \n-Variability for each component (%):\n")
  # print(round(pr.x[l] * 100, 2))
  print(round(pr.x2[l] * 100, 2))
  }
 if (biplot){
  j=1
  while (j<=lenC) {
        if (names(C)[j]=="ask") {
           ask=C[[j]]
           j=lenC +1
           }
        else { j=j+1; ask=FALSE  }
  }
   #dev.new()
  
   if (ask) {
          par(mfrow=c(1,1))
          dev.interactive()
          oask <- devAskNewPage(TRUE)
          on.exit(devAskNewPage(oask))
          for (i in 1:le) {
          ts.plot(rotation[,l[i]],ylab=c("loadings",l[i],sep=""),
      main=c(paste("components",l[i],"- Expl. Var. ",round(pr.x2[l[i]] * 100, 2),"%",sep="")))
      if (i<le)
      for (j in (i+1):le) {

            if (nrow(out2)<50)   {
                        plot(p[,c(i,l[j])],main="BIPLOT",type="n")
                        text(p[,c(i,l[j])])#,rownames(out2))
                        }
             else                         plot(p[,c(i,l[j])],main="BIPLOT")
            if (nrow(out2)<50)      {
                           plot(p[,c(l[j],i)],main="BIPLOT",type="n")
                           text(p[,c(l[j],i)])#,rownames(out2))
               }
           else  plot(p[,c(l[j],i)],main="BIPLOT")
      } }  
    }
    else   {
      par(mfrow=c(le,le))
      for (kk in 1:(le*le)){
        plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
      }
    for (i in 1:le) {
      par(mfg=c(i,i))
      plot(rotation[,l[i]],ylab=c("loadings",l[i],sep=""),type="l",
       main=c(paste("Component",l[i],"- Expl. Var. ",round(pr.x2[l[i]] * 100, 2),"%",sep="")))
      if (i<le)
      for (j in (i+1):le) {
            par(mfg=c(i,j))
            if (nrow(out2)<50)     {
                      plot(p[,c(i,l[j])],main="BIPLOT")
                      text(p[,c(i,l[j])])
                      }
            else plot(p[,c(i,l[j])],main="BIPLOT")
            par(mfg=c(j,i))
            if (nrow(out2)<50)            {
               plot(p[,c(l[j],i)],main="BIPLOT",type="n")
               text(p[,c(l[j],i)])#,rownames(out2))
               }
            else plot(p[,c(l[j],i)],main="BIPLOT")
     }  }    }  }
  #return(invisible(pr.com))
  return(invisible(object))
}

