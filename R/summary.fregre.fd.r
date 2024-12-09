#' @title Summarizes information from fregre.fd objects.
#' 
#' @description Summary function for \code{\link{fregre.pc}}, \code{\link{fregre.basis}},
#' \code{\link{fregre.pls}}, \code{\link{fregre.np}}\cr and
#' \code{\link{fregre.plm}} functions.
#' 
#' Shows:\cr \tabular{ll}{ \tab -Call.\cr \tab -R squared.\cr \tab -Residual
#' variance.\cr \tab -Index of possible atypical curves or possible
#' outliers.\cr \tab -Index of possible influence curves.\cr } If the
#' \code{fregre.fd} object comes from the \code{\link{fregre.pc}} then shows:
#' \tabular{ll}{ \tab -Variability of explicative variables explained by
#' Principal Components.\cr \tab -Variability for each principal components
#' -PC-.\cr }
#' 
#' If draw=TRUE plot: \cr \tabular{ll}{ \tab -y vs y fitted values.\cr \tab
#' -Residuals vs fitted values.\cr \tab -Standarized residuals vs fitted
#' values.\cr \tab -Levarage.\cr \tab -Residual boxplot.\cr \tab
#' -Quantile-Quantile Plot (qqnorm).\cr } If \code{ask}=FALSE draw graphs in
#' one window, by default. If \code{ask}=TRUE, draw each graph in a window,
#' waiting to confirm.
#' 
#' @aliases summary.fregre.fd  summary.fregre.lm plot.summary.lm
#' summary.fregre.igls print.fregre.igls  print.fregre.plm print.fregre.fd
#' @param object Estimated by functional regression, \code{fregre.fd} object.
#' @param times.influ Limit for detect possible infuence curves.
#' @param times.sigma Limit for detect possible oultiers or atypical curves.
#' @param draw =TRUE draw estimation and residuals graphics.
#' @param \dots Further arguments passed to or from other methods.
#' @return 
#' \itemize{
#' \item \code{Influence}: Vector of influence measures.
#' \item \code{i.influence}: Index of possible influence curves.
#' \item \code{i.atypical}: Index of possible atypical curves or possible outliers.
#' }
#' @author Manuel Febrero-Bande and Manuel Oviedo de la Fuente \email{manuel.oviedo@@udc.es}
#' @seealso Summary function for \code{\link{fregre.pc}},
#' \code{\link{fregre.basis}}, \code{\link{fregre.pls}}, \cr
#' \code{\link{fregre.np}} and \code{\link{fregre.plm}}.
#' @keywords print
#' @examples
#' \dontrun{
#' # Ex 1. Simulated data
#' n= 200;tt= seq(0,1,len=101)
#' x0<-rproc2fdata(n,tt,sigma="wiener")
#' x1<-rproc2fdata(n,tt,sigma=0.1)
#' x<-x0*3+x1
#' beta = tt*sin(2*pi*tt)^2
#' fbeta = fdata(beta,tt)
#' y<-inprod.fdata(x,fbeta)+rnorm(n,sd=0.1)
#' 
#' # Functional regression
#' res=fregre.pc(x,y,l=c(1:5))
#' summary(res,3,ask=TRUE)
#' 
#' res2=fregre.pls(x,y,l=c(1:4))
#' summary(res2)
#' 
#' res3=fregre.pls(x,y)
#' summary(res3)
#' }
#' 
#' @export 
summary.fregre.fd<-function(object,times.influ=3,times.sigma=3,draw=TRUE,...){
    x<-object$fdataobj$data
    t=object$fdataobj$argvals
    y<-object$y
    isfdata<-is.fdata(y)
    n=nrow(x)
    if (!isfdata) {
     up=mean(object$residuals)+times.sigma*sqrt(object$sr2)
     lo=mean(object$residuals)-times.sigma*sqrt(object$sr2)
     i.atypical=which(object$residuals>up|object$residuals<lo)
     lim.influ=fdata.trace(object$H)/n
     influence=diag(object$H)
     i.influence=which(influence>times.influ*lim.influ)
     if (length(i.influence) == 0) i.influence=NA
     if (length(i.atypical) == 0) i.atypical=NA
     }
     if (object$call[[1]]=="fregre.pc") {
     if (object$lambda==0)     {
     cat(" *** Summary Functional Data Regression with Principal Components ***\n")
      object$lm$call<-object$call
      print(summary(object$lm))}
      else  {
     cat(" *** Summary Functional Ridge Regression with Principal Components*** \n\n")
            cat("-Call: ");    print(object$call)
            cat("\n")
            print(object$coefs)
            cat("\n-R squared: ",object$r2)
#            cat("\n-Residual variance: ",object$sr2,"\n")
     cat("\n-Residual variance: ",
            object$sr2,"on ",n-object$df.residual," degrees of freedom\n")
       cat("-Penalization parameter (lambda): ",object$lambda,"\n")
            }
            
#     object$lm$call<-object$call
#     print(summary(object$lm))
     var.1<-apply(object$fdata.comp$x, 2, var)
     pr.x= var.1/sum(var.1)
 cat("\n-With",length(object$l),"Principal Components is  explained ",round(sum(pr.x[object$l])*100
 ,2),"%\n of the variability of explicative variables. \n
-Variability for each  principal components -PC- (%):\n")
    print(round(pr.x[object$l] * 100, 2))
    }
     if (object$call[[1]]=="fregre.ppc") {
     cat(" *** Summary Functional Regression with Penalized Principal Components ***\n")
      object$lm$call<-object$call
      print(summary(object$lm))
#            cat("\n-R squared: ",object$r2)
#      cat("\n-Residual variance: ",
#            object$sr2,"on ",n-object$df.residual," degrees of freedom\n")
       cat("-Lambda penalty: ",object$lambda)
       #     object$lm$call<-object$call
#     print(summary(object$lm))
     var.1<-apply(object$fdata.comp$x, 2, var)
     pr.x= var.1/sum(var.1)
 cat("\n-With",length(object$l),"Principal Components is explained ",round(sum(pr.x[object$l])*100
 ,2),"%\n of the variability of explicative variables. \n -Variability for each  principal components -PC- (%):\n")
    print(round(pr.x[object$l] * 100, 2))
    }
     if (object$call[[1]]=="fregre.pls") {
##     cat(" *** Summary Functional Data Regression with Partial Least Squares ***\n")
##      object$lm$call<-object$call
##      print(summary(object$lm))
     cat(" *** Summary Functional Regression with Partial Least Squares*** \n\n")
            cat("-Call: ");    print(object$call)
            cat("\n")
            print(object$coefficients)
              cat("\n-R squared: ",object$r2)
#            cat("\n-Residual variance: ",object$sr2,"\n")
              cat("\n-Residual variance: ",
            object$sr2,"on ",n-object$df.residual," degrees of freedom\n")

#     object$lm$call<-object$call
#     print(summary(object$lm))
#     var.1<-apply(object$fdata.comp$x, 2, var)
#     pr.x= var.1/sum(var.1)
# cat("\n-With",length(object$l),"Partial Least Squares is  explained ",round(sum(pr.x[object$l])*100
# ,2),"%\n of the variability of explicative variables. \n
#-Variability for each  Partial Least Squares -PLS- (%):\n")
#    print(round(pr.x[object$l] * 100, 2))
    }
     if (object$call[[1]]=="fregre.ppls") {
     cat(" *** Summary Functional Regression with Penalized Partial Least Squares ***\n")
            cat("-Call: ");    print(object$call)
            cat("\n")
            print(object$coefs)
            cat("\n-R squared: ",object$r2)
#            cat("\n-Residual variance: ",object$sr2,"\n")
              cat("\n-Residual variance: ",
            object$sr2,"on ",n-object$df.residual," degrees of freedom\n")
       cat("-Lambda penalty: ",object$lambda)
      #     object$lm$call<-object$call
#     print(summary(object$lm))
#     var.1<-apply(object$fdata.comp$x, 2, var)
#     pr.x= var.1/sum(var.1)
# cat("\n-With",length(object$l),"Partial Least Squares is  explained ",round(sum(pr.x[object$l])*100
# ,2),"%\n of the variability of explicative variables. \n -Variability for each Partial Least Squares -PLS- (%):\n")
#    print(round(pr.x[object$l] * 100, 2))
    }
     if (object$call[[1]]=="fregre.basis") {
     cat(" *** Summary Functional Data Regression with representation in Basis *** \n")
     if (object$lambda==0)     {object$lm$call<-object$call
                                print(summary(object$lm))}
      else  {
            cat("-Call: ");    print(object$call)
            cat("\n")
            print(object$coefs)
            cat("\n-R squared: ",object$r2)
            cat("\n-Residual variance: ",
            object$sr2,"on ",n-object$df.residual," degrees of freedom\n")
            }
    }
     if (object$call[[1]]=="fregre.basis.cv") {
      cat(" *** Summary Functional Data Regression with representation in Basis *** \n\n")
      cat("-Call: ");    print(object$call)
            cat("\n")
            print(object$coefficients)
            cat("\n-R squared: ",object$r2)
#            cat("\n-Residual variance: ",object$sr2,"\n")
cat("\n-Residual variance: ",
            object$sr2,"on ",n-object$df.residual," degrees of freedom\n")
      cat("-Optimal Beta Basis: \n")
      print(object$basis.b.opt)
      cat("\n-Optimal lambda penalty=",object$lambda.opt,"\n")
#      print(object$Lfdobj)
    }
     if (object$call[[1]]=="fregre.np") {
       cat(" *** Summary Functional Non-linear Model *** \n\n")
       cat("-Call: ");    print(object$call)
       cat("\n-Bandwidth (h): ",object$h.opt)
       cat("\n-R squared: ",object$r2)
   #   cat("\n-Residual variance: ",object$sr2,"\n")
       cat("\n-Residual variance: ",
              object$sr2,"on ",n-object$df.residual," degrees of freedom\n")
      }
     if (object$call[[1]]=="fregre.np.cv") {
      cat(" *** Summary Functional Non-linear Model *** \n\n")
      cat("-Call: ");    print(object$call)
      cat("\n-Bandwidth (h): ",object$h.opt)
      cat("\n-R squared: ",object$r2)
    #    cat("\n-Residual variance: ",object$sr2,"\n")
      cat("\n-Residual variance: ",
            object$sr2,"on ",n-object$df.residual," degrees of freedom\n")
    }
     if (object$call[[1]]=="fregre.plm") {
       cat(" *** Summary Functional Semi-linear Model *** \n\n")
       cat("-Call: ");    print(object$call)
       cat("\n-Coefficients:  non functional covariates\n")
       print(object$coefficients)
       cat("\n-Bandwidth (h): ",object$h.opt)
       cat("\n-R squared: ",object$r2)
      #    cat("\n-Residual variance: ",object$sr2,"\n")
       cat("\n-Residual variance: ",
              object$sr2,"on ",n-object$df.residual," degrees of freedom\n")
       }
    if (!isfdata) {
     cat("-Names of possible atypical curves: ");
     if (is.na(i.atypical[1]))     cat("No atypical curves \n")
     else   if (length(i.atypical)<11)  cat(rownames(x)[i.atypical],"\n")
           else cat(rownames(x)[i.atypical[1:10]],
           "\n It prints only the 10 most atypical curves. \n")
     cat("-Names of possible influence curves: ");
     if (is.na(i.influence[1])) cat("No influence curves \n")
     else  if (length(i.influence)<11) cat(rownames(x)[i.influence],"\n")
     else cat(rownames(x)[i.influence[1:10]],
     "\n It prints only the 10 most influence curves \n")
   }
   else  return(invisible(object)) #draw=FALSE
    if (draw) {
      oldpar <- par()
      C<-match.call()
      lenC=length(C)
      j=1
      while (j<=lenC) {
        if (names(C)[j]=="ask") {
           ask=C[[j]]
           j=lenC +1             }
        else {      j=j+1
                    ask=FALSE             }
       }
       if (ask) {
          par(mfrow=c(1,1))
          dev.interactive()
          oask <- devAskNewPage(TRUE)
          on.exit(devAskNewPage(oask))
          }
       else   par(mfrow=c(2,3))
 plot(object$fitted.values,y,xlab="Fitted values",main=paste("R-squared=",
     round(object$r2,2)))
 plot(object$fitted.values,object$residuals,ylab="Residuals",
    xlab="Fitted values",main="Residuals vs fitted.values")
    text(object$fitted.values[i.atypical],object$residuals[i.atypical],
    rownames(x)[i.atypical],cex=0.7)
    abline(h=mean(object$residuals),lwd=1,lty=2)
    abline(h=up,col=2,lwd=2,lty=2)
    abline(h=lo,col=2,lwd=2,lty=2)
#############
resid.sd=sqrt(abs(object$residuals/sd(object$residuals)))
main= "Scale-Location"
ylab23<-"Standardized residuals"
ylim <- c(0, max(resid.sd, na.rm = TRUE))
yl <- as.expression(substitute(sqrt(abs(YL)), list(YL = as.name(ylab23))))
plot(object$fitted.values,resid.sd, xlab = "Fitted values",
 ylab = yl, main = main,ylim = ylim)
 text(object$fitted.values[i.atypical],resid.sd[i.atypical],
 rownames(x)[i.atypical],cex=0.7)
 plot(diag(object$H),1:nrow(x),xlab="Leverage",ylab="Index.curves",
    main="Leverage")
text(diag(object$H)[i.influence],i.influence,
rownames(x)[i.influence],cex=0.7)
abline(v=times.influ*lim.influ,col=2,lwd=2,lty=2)
#  plot(density(object$residuals),main="Residuals")
    qqnorm(object$residuals,main="Residuals")
    boxplot(object$residuals,main="Residuals")
    par(mfrow=c(1,1))
    }
#    cat("\n")
return(invisible(list("Influence"=influence,"i.influence"=i.influence,
"i.atypical"=i.atypical)))
}


#' @export 
print.fregre.fd<-function (x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\n-Call: ", deparse(x$call), "\n", sep = "")
  if (length(coef(x))) {
    cat("\n-Coefficients:\n")
    print.default(format(coef(x), digits = digits), print.gap = 2,
                  quote = FALSE)
    if (x$call[[1]]=="fregre.lm")      print(x$beta.est[[2]])
  }
  else {
    if (x$call[[1]]=="fregre.np" || x$call[[1]]=="fregre.np.cv") {
      cat("\n-Bandwidth (h): ",x$h.opt)
    }
  }
  cat("\n-R squared: ",x$r2)
  cat("\n-Residual variance: ",x$sr2,"\n")
  invisible(x)
}

#' @export 
print.fregre.plm<-function (x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("-Call: ");    print(x$call)
  cat("-Coefficients:  non functional covariates\n")
  print(x$coefficients)
  cat("-Bandwidth (h): ",x$h.opt)
  cat("\n-R squared: ",x$r2," -Residual variance: ",x$sr2)
}
