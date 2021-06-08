# @S3method summary fregre.lm


#' @export 
summary.fregre.lm<-function (object, correlation = FALSE, 
                             symbolic.cor = FALSE, times.influ=3,
                             times.sigma=3,...) {
if (object$lambda)  {
  x<-object$fdataobj$data
  n<-length(object$residuals)
  cat(" *** Summary Functional Linear Regression with Penalization*** \n\n")
  cat("-Call: fregre.lm")#;    print(object$call)
  cat("\n")
  cat("-Residuals:\n ")
  print(summary(object$residuals))#;  cat("\n")
  cat("-Coefficients: \n")
  if (is.null(object$coefs)) print(object$coefficients)
  else  print(object$coefs)
  cat("\n-R squared: ",object$r2)
  cat("\n-Residual variance: ",
      object$sr2,"on ",round(n-object$df)," degrees of freedom\n")
  cat("-Penalization parameter (lambda): \n")
  print(unlist(object$lambda.opt))
  

  up=mean(object$residuals)+times.sigma*sqrt(object$sr2)
  lo=mean(object$residuals)-times.sigma*sqrt(object$sr2)
  i.atypical=which(object$residuals>up|object$residuals<lo)
  lim.influ=fdata.trace(object$H)/n
  influence=diag(object$H)
  i.influence=which(influence>times.influ*lim.influ)
  if (length(i.influence) == 0) i.influence=NA
  if (length(i.atypical) == 0) i.atypical=NA
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
  
  }  else {
      summary.lm(object,correlation=correlation, symbolic.cor=symbolic.cor,...)
      }
}  

#' @export 
plot.summary.lm<-function(x,times.influ=3,times.sigma=3,...){
  object<-x
  if (object$lambda)  {
    x<-object$fdataobj$data
    t=object$fdataobj$argvals
    y<-object$y
    isfdata<-is.fdata(y)
    n=length(object$residuals)
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

   else  return(invisible(object)) #draw=FALSE
  #  if (draw) {
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
  plot(object$fitted.values,y,xlab="Fitted values",main=paste("R-squared=", round(object$r2,2)))
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
 dd<-diag(object$H)
 plot(dd,1:n,xlab="Leverage",ylab="Index.curves",
    main="Leverage")
 text(diag(object$H)[i.influence],i.influence,
  rownames(x)[i.influence],cex=0.7)
  abline(v=times.influ*lim.influ,col=2,lwd=2,lty=2)
#  plot(density(object$residuals),main="Residuals")
    qqnorm(object$residuals,main="Residuals")
    boxplot(object$residuals,main="Residuals")
    par(mfrow=c(1,1))
    #}
#    cat("\n")
return(invisible(list("Influence"=influence,"i.influence"=i.influence,
"i.atypical"=i.atypical)))
  }
  else {
    plot(object)
  }
}  

