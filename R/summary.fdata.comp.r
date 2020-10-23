#' @title Correlation for functional data by Principal Component Analysis
#' 
#' @description Summary of functional principal components
#' 
#' @param object fdata.comp class object calculated by: \code{fdata2pc},
#' \code{fdata2pls}, \code{fregre.pc} or \code{fregre.pls}.
#' @param biplot =TRUE draw the biplot and PC (or PLS) components.
#' @param \dots Further arguments to be passed to plot
#' @return If \code{biplot}=TRUE, are displaying the biplot between the PC (or PLS) components.\cr
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
#' summary(pls1,pch=19,col="red")
#' }
#' 
#' @export
summary.fdata.comp=function(object,biplot=TRUE,...) {
  if (inherits(object, "fdata.comp"))         {
     pr.com<-object
	 d<-object$d
	if (is.null(object$rn)) rn<-0
	 ccall=object$call[[1]]
     if (ccall=="fdata2pls" | ccall=="fdata2ppls"){
         y<-object$y
     }
  } else if (inherits(object, "fregre.fd") & !is.null(object$fdata.comp))     {
     pr.com<-object$fdata.comp
     y<-object$y	
	 d<-object$fdata.comp$d
	if (is.null(object$rn)) rn<-0
	 ccall<-object$fdata.comp$call[[1]]
  } else stop("Error in input data")
 out2=pr.com$x
 types<-colnames(out2)
 #print(types)
 l<-object$l
 le=length(l)
 if (is.null(object$rn)) rn<-0
 if (ccall=="fdata2pc"){

	pr.x2<-round((d^2+rn)/sum(d^2+rn),4)
	names(pr.x2)<-colnames(out2)
 }
 if (ccall=="fdata2pls") {
  pr.xy= round(cor(cbind(y,out2))[1,-1]^2,4)
  names(pr.xy)<-colnames(out2)
 } 


 cat("\n     - SUMMARY:  ",object$call[[1]]," object   -\n")
 if (ccall=="fdata2pc"){
 cat("\n-With",le," components are explained ",round(sum(pr.x2[l])*100
 ,2),"%\n of the variability of explicative variables.\n \n-Variability for each component (%):\n")
  print(round(pr.x2[l] * 100, 2))
 }
if (ccall=="fdata2pls"){
 cat("\n   R squared with y:\n")
  print(round(pr.xy[l], 4))
 }
 
  C<-match.call()
  lenC=length(C)
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
titlediag=if(ccall=="fdata2pc"){
		 paste(types,"- Expl. Var. ",round(pr.x2* 100, 2),"%",sep="")
		 }else{
		 paste(types,"- R Squared. ",round(pr.xy* 100, 2),"%",sep="")}

   if (ask) {
          par(mfrow=c(1,1))
          dev.interactive()
          oask <- devAskNewPage(TRUE)
          on.exit(devAskNewPage(oask))
          for (i in 1:le) {
		  for (j in 1:le){
		  if (i==j) {
          plot(pr.com$rotation[l[i]],ylab=paste("Comp.",l[i]),main=titlediag[l[i]])
		} else {
            if (nrow(out2)<50)   {
                        plot(out2[,c(l[i],l[j])],main="BIPLOT",type="n",...)
                        text(out2[,c(l[i],l[j])])#,rownames(out2))
                        }
             else                         plot(out2[,c(l[i],l[j])],main="BIPLOT",...)
            if (nrow(out2)<50)      {
                           plot(out2[,c(l[j],l[i])],main="BIPLOT",type="n",...)
                           text(out2[,c(l[j],l[i])])#,rownames(out2))
               }
           else  plot(out2[,c(l[j],l[i])],main="BIPLOT",...)
      } 
    }}}
    else   {
      par(mfrow=c(le,le))
      for (kk in 1:(le*le)){
        plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
      }
    for (i in 1:le) {
      par(mfg=c(i,i))
      plot(pr.com$rotation[l[i]],ylab=paste("Comp.",l[i]),main=titlediag[l[i]])
      if (i<le)
      for (j in (i+1):le) {
            par(mfg=c(i,j))
            if (nrow(out2)<50)     {
                      plot(out2[,c(l[i],l[j])],main="BIPLOT",...)
                      text(out2[,c(l[i],l[j])])
                      }
            else plot(out2[,c(l[i],l[j])],main="BIPLOT",...)
            par(mfg=c(j,i))
            if (nrow(out2)<50)            {
               plot(out2[,c(l[j],l[i])],main="BIPLOT",type="n",...)
               text(out2[,c(l[j],l[i])])#,rownames(out2))
               }
            else plot(out2[,c(l[j],l[i])],main="BIPLOT",...)
     }  }    }  }
  #return(invisible(pr.com))
  return(invisible(object))
}

