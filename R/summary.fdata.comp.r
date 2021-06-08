#' @title Correlation for functional data by Principal Component Analysis
#' 
#' @description Summary of functional principal components
#' 
#' @param object fdata.comp class object calculated by: \code{fdata2pc},
#' \code{fdata2pls}, \code{fregre.pc} or \code{fregre.pls}.
#' @param biplot =TRUE draw the biplot and PC (or PLS) components.
#' @param draw =TRUE original curves  and their basis representation are plotted
#' @param index NULL, by default the first n curves are plotted, where n = min(4, length(fdataobj)). 
#' Otherwise, the curves indicated in the \code{index} vector are plotted.
#' @param \dots Further arguments passed to or from other methods.
#' @return If \code{corplot}=TRUE, are displaying the biplot between the PC (or PLS) components.\cr
#'  If \code{ask}=TRUE, draw each graph in a window, waiting to confirm the change
#' of page with a click of the mouse or pressing ENTER.  If \code{ask}=FALSE draw graphs in one window.
#' @author Manuel Febrero-Bande and Manuel Oviedo de la Fuente \email{manuel.oviedo@@udc.es}
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
#' pc1 <- fdata2pc(x,3)
#' summary(pc1)
#' y <- inprod.fdata(x,fbeta) #+ rnorm(n,sd=0.1)
#' pls1 <- fdata2pls(x,y,2)
#' summary(pls1)
#' res <- fdata2basis(x0,pc1)
#' summary(res)

#' }
#' 
#' @export
summary.fdata.comp=function(object,biplot=TRUE,...) {
  pc <- TRUE
  if (inherits(object, "fdata.comp"))         {
     a1=TRUE
     pr.com<-object
     #if (is.null(object$y)) {
        if (object$call[[1]]=="fdata2pls" | object$call[[1]]=="create.pls.basis"
            | object$call[[1]]=="fdata2ppls" |  object$call[[1]]=="fdata2ppls")
         y<-object$y                  
        pc <- FALSE
     
     #}
  } else if (inherits(object, "fregre.fd"))     {
     a1=FALSE
     pr.com<-object$fdata.comp
     y<-object$y
  } else stop("Error in input data")
 out2 <- pr.com$coefs
 l <- object$l
 le <- length(l)
 rotation=aperm(pr.com$basis$data)
 p <- pr.com$coefs
 if (!pc) {
  var.1<-apply(p, 2, cor,y)
  pr.x2 <- object$Xvar / object$Xtotvar
  

 } else {
   if (is.null(object$rn)) rn<-0
     d<-object$d
     if (is.null(object$rn)) rn<-0
     pr.x2<-(d^2+rn)/sum(d^2+rn)
     #colnames(pr.x2)<-colname
     names(pr.x2)<-names(d)
       }
 C <- match.call()
 lenC=length(C)
 # cor.y.pc=rep(NA,le)
 # xxx=cbind(y,pr.com$x)
 # cor.y.pc=round(cor(xxx[,c(1,l+1)]),3)[1,-1]
 types<-colnames(pr.com$x)
 if (pc){
 cat("\n     - SUMMARY:  ",object$call[[1]]," object   -\n")
# if (object$call[[1]]=="create.pc.basis")    {   d<-1:length(d) }
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
          #ts.plot(rotation[,l[i]],ylab=c("loadings",l[i],sep=""),
            ts.plot(rotation[,i],ylab=c("loadings",l[i],sep=""),
      main=c(paste("components",l[i],"- Expl. Var. ",round(pr.x2[l[i]] * 100, 2),"%",sep="")))
      if (i<le)
      for (j in (i+1):le) {

            if (nrow(out2)<50)   {
                        plot(p[,c(i,j)],main="BIPLOT",type="n")
                        text(p[,c(i,j)])#,rownames(out2))
                        }
             else                         plot(p[,c(i,j)],main="BIPLOT")
            if (nrow(out2)<50)      {
                           plot(p[,c(j,i)],main="BIPLOT",type="n")
                           text(p[,c(j,i)])#,rownames(out2))
               }
           else  plot(p[,c(j,i)],main="BIPLOT")
      } }  
    }
    else   {
      par(mfrow=c(le,le))
      for (kk in 1:(le*le)){
        plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
      }
    for (i in 1:le) {
      par(mfg=c(i,i))
      plot(rotation[,i],ylab=c("loadings",l[i],sep=""),type="l",
       main=c(paste("Component",l[i],"- Expl. Var. ",round(pr.x2[l[i]] * 100, 2),"%",sep="")))
      if (i<le)
      for (j in (i+1):le) {
            par(mfg=c(i,j))
            if (nrow(out2)<50)     {
                      plot(p[,c(i,j)],main="BIPLOT")
                      text(p[,c(i,j)])
                      }
            else plot(p[,c(i,j)],main="BIPLOT")
            par(mfg=c(j,i))
            if (nrow(out2)<50)            {
               plot(p[,c(j,i)],main="BIPLOT",type="n")
               text(p[,c(j,i)])#,rownames(out2))
               }
            else plot(p[,c(j,i)],main="BIPLOT")
     }  }    }  }
  #return(invisible(pr.com))
  return(invisible(object))
}

#' @export
summary.basis.fdata=function(object, draw=TRUE, index=NULL) {
  #object <- bb
  cat("\n     - SUMMARY -\n")
  le <- length(object$basis)
  R <- numeric(le)
  #print(R)
  n <- nrow(object$fdataobj)
  R2 <- matrix(NA,n,le)
  colnames(R2)<-paste0(substr(object$type,1,3),"(1:",1:le,")")
  rownames(R2)<-rownames(object$fdataobj)
  Xcen <- fdata.cen(object$fdataobj)$Xcen
  type=FALSE
  cat("Type of basis:",object$type,"\nNum. of basis:",le,"\nRangeval:",object$basis$rangeval,"\n")
  for (l in 1:le){
    if (object$type =="pc" | object$type == "pls"){
      type=TRUE
      #xmean <- func.mean(object$fdataobj)
      xmean <- object$mean
      fdata.est <- basis2fdata(object$coefs[,1:l,drop=F],
                               object$basis[1:l])
    } else{ fdata.est <- basis2fdata(object$coefs[,1:l,drop=F],
                                     object$basis[1:l])
    #R2[l] <- 1 - sum(norm.fdata(fdata.est-object$fdataobj)^2)/
    #   sum(norm.fdata(fdata.cen(object$fdataobj)$Xcen)^2)
    }
    R[l] <- 1-sum(norm.fdata(fdata.est-object$fdataob)^2)/sum(norm.fdata(Xcen)^2)   
    R2[,l] <- 1 - norm.fdata(fdata.est-object$fdataobj)^2/norm.fdata(Xcen)^2    
  }
  if (is.null(index)) {
    ymin <- min(4,n);     index<- 1:ymin
  }
  txt <- paste(length(index),"curves (in grey) and their basis representation (in red) are plotted")
  if (draw) {
    # yl <- c(min(object$fdataobj,fdata.est),max(object$fdataobj,fdata.est))
    yl <- c(min(object$fdataobj[index],fdata.est[index]),max(object$fdataobj[index],fdata.est[index]))
    # mn <- expression( 1 - frac(paste( "||X - ",hat(X),"||"),paste( "||X - ",bar(X),"||")))
    mn <-  expression( paste("X(t) vs ",hat(X),"(t)"))
    plot(object$fdataobj[index],main=mn,col="grey",ylim=yl)
    lines(fdata.est[index],lty=2,col=2)
  }
  #  En el texto pon 1- Sum(||---||^2)/Sum(||---||^2)
  
  #  message(paste( "||X - ",expression(bar(X)),"||"))
  cat("\nMeasure of fit: 1 - Sum||X - hat(X)||^2  / Sum||X - bar(X)||^2:\n")
  names(R) <- paste0("basis(1:",1:le,")",sep="")
  if (!type)    R <- R[le]
  #print(expression( 1 - paste( "||X - ",hat(X),"|| / ||X - ",bar(X),"||")))
  print(R)
  if (draw) cat("\n",txt)
  
  return(invisible(R2))
}