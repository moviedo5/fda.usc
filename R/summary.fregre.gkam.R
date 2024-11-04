#' @title Summarizes information from fregre.gkam objects.
#' 
#' @description Summary function for \code{\link{fregre.gkam}} function.
#' 
#' \tabular{ll}{ \tab -Family used.\cr \tab -Number or iteration of algorithm
#' and if it has converged. \cr \tab -Residual and null deviance.\cr \tab
#' -Number of data.\cr } 
#' 
#' Produces a list of summary information for a fitted
#' fregre.np object for each functional covariate.  
#' 
#' \tabular{ll}{ \tab
#' -Call.\cr \tab -R squared.\cr \tab -Residual variance.\cr \tab -Index of
#' possible atypical curves or possible outliers.\cr \tab -Index of possible
#' influence curves.\cr } If draw=TRUE plot: \cr \tabular{ll}{ \tab -y vs y
#' fitted values.\cr \tab -Residuals vs fitted values.\cr \tab -Residual
#' boxplot.\cr \tab -Quantile-Quantile Plot (qqnorm).\cr \tab -Plot for a each
#' single model term.\cr }
#'  If \code{ask}=FALSE draw graphs in one window, by
#' default. If \code{ask}=TRUE, draw each graph in a window, waiting to
#' confirm.
#' 
#' @aliases summary.fregre.gkam print.fregre.gkam
#' @param object Estimated by functional regression, \code{fregre.fd} object.
#' @param draw =TRUE draw estimation and residuals graphics.
#' @param selec Allows the plot for a single model term to be selected for
#' printing. e.g. if you just want the plot for the second smooth term set
#' selec=2.  .
#' @param times.influ Limit for detect possible infuence curves.
#' @param \dots Further arguments passed to or from other methods.
#' 
#' @author Manuel Febrero-Bande and Manuel Oviedo de la Fuente \email{manuel.oviedo@@udc.es}
#' @seealso Summary function for \code{\link{fregre.gkam}}.
#' @keywords print
#' @examples
#' \dontrun{
#' # Time consuming
#' data(tecator)
#' ind<-1:129
#' ab=tecator$absorp.fdata[ind]
#' ab2=fdata.deriv(ab,2)
#' yfat=as.integer(cut(tecator$y[ind,"Fat"],c(0,15,100)))-1
#' xlist=list("df"=data.frame(yfat),"ab2"=ab2,"ab"=ab)
#' f<-yfat~ab+ab2
#' res=fregre.gkam(f, data = xlist, family = binomial("logit"),
#'                 control = list(maxit = 2))
#' summary(res)
#' 
#' }
#' 
#' @export 
summary.fregre.gkam<-function(object,draw=TRUE,selec=NULL
                              ,times.influ=3,...){
  #effects=TRUE
  cat(" *** Summary Functional Data Regression with backfiting algorithm *** \n")
  print(object$family,"\n")
  nobs <- n<-length(object$y)
  cat("alpha=",as.numeric(format(object$alpha,digits=3)),"  n= ",n,"\n")
  cat("Algorithm converged?",ifelse(object$converged,"Yes","No")," Number  of iterations ",object$iter,"\n\n")
  res<- object$result
  namesres<-names(res)
  lenres<-length(res)
  eta<-object$linear.predictors
  x<-object$fdataobj$data
  t=object$fdataobj$argvals
  y<-object$y
  if (object$family$family=="binomial") y<-as.numeric(factor(y,labels=c(0,1)))-1
  lenres<-length(object$result)
  tab<-matrix(NA,nrow=(lenres),ncol=3)
  rownames(tab)<-paste("f(",namesres,")",sep="")
  colnames(tab)<-c("h","cor(f(X),eta)","edf")
  cat("****    ****    ****    ****    ****    ****\n")
  for (i in 1:lenres) {
    hh<-res[[i]]$h.opt
    tab[i,1]<-as.numeric(format(hh,digits=3))
    tab[i,2]<-round(cor(object$effects[,i],object$fitted.values),3)
    tab[i,3]<-round(object$eqrank[[i]],1)
  }
  print(tab)
  cat("****    ****    ****    ****    ****    ****\n\n")
  cat("edf: Equivalent degrees of freedom\n")
  rowname<-rownames(res[[1]]$fdataobj)
  if (is.null(rowname)) rowname<-1:n
  else rowname<-rownames(x)
  H<-object$H
  influence=diag(H)
  lim.influ=sum(influence)/n
  i.influence=which(influence>times.influ*lim.influ)
  if (length(i.influence) == 0) i.influence=NA
  residual.df <- nobs - sum(object$eqrank)
  w<-object$prior.weights
  ycen = y - mean(y)
  r.sq <- round(1 - var(w * (y - object$fitted.values))*(nobs - 1)/(var(w * y)*residual.df),3)
  r.sq2 <- round( 1 - sum((object$fitted.values-y)^2)/sum(ycen^2),3)
  cat("Residual deviance=",round(object$deviance,3)," Null deviance=",round(object$null.deviance,3),"\n")
  dev.expl<-100*round((object$null.deviance - object$deviance)/object$null.deviance,3)
  #     cat("R-sq. = ",1 - sum((object$residuals)^2)/sum(ycen^2),"  Deviance explained = ", dev.expl,"%","\n")
  cat("AIC= ",round(object$aic,3),"Deviance explained=",dev.expl,"%","\n")
  cat("R-sq.=",r.sq2," R-sq.(adj)=",r.sq,"\n")
  cat("Names of possible influence curves: ");
  if (is.na(i.influence[1])) cat("No influence curves \n")
  else  if (length(i.influence)<11) cat(rowname[i.influence],"\n")
  else cat(rowname[i.influence[1:10]],
           "\n It prints only the 10 most influence curves \n")
  cat("\n")
  if (draw) {
    le<-lenres
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
    
    ############
    if (ask) {          par(mfrow=c(1,1))
    }
    else                par(mfrow = c(2,3))
    # plot(object$fitted.values,y,xlab="Fitted values",main=paste("R-squared=",
    #     round(r.sq,2)))
    col1<-rep(1,n)
    #if (is.factor(y)) col1=as.numeric(y)
    #    if (family$family=="binomial" & !is.factor(y)) y<-factor(y)
    if (object$family$family=="binomial") col1=y+2
    
    plot(object$linear.predictors,y,xlab="Linear predictors",main=paste("R-sq=",
                                                                        round(r.sq2,2)))
    points(object$linear.predictors,object$fitted.values,col=col1)
    legend(min(object$linear.predictors),0.95,legend=c("y","fitted values"),col=1:2,lty=1:2,
           box.col="white",xjust=0)
    plot(object$fitted.values,object$residuals,ylab="Residuals",
         xlab="Fitted values",main="Residuals vs fitted.values",col=col1)
    abline(h=mean(object$residuals),lwd=1,lty=2)
    #############
    resid.sd=sqrt(abs(object$residuals/sd(object$residuals)))
    main= "Scale-Location"
    ylab23<-"Standardized residuals"
    ylim <- c(0, max(resid.sd, na.rm = TRUE))
    yl <- as.expression(substitute(sqrt(abs(YL)), list(YL = as.name(ylab23))))
    plot(object$fitted.values,resid.sd, xlab = "Fitted values",
         ylab = yl, main = main,ylim = ylim,col=col1)
    # text(object$fitted.values[i.atypical],resid.sd[i.atypical],rowname[i.atypical],cex=0.7)
    plot(diag(H),1:n,xlab="Leverage",ylab="Index.curves",main="Leverage",col=col1)
    text(influence[i.influence],i.influence,rowname[i.influence],cex=0.7)
    abline(v=times.influ*lim.influ,col=4,lwd=2,lty=2)
    qqnorm(object$residuals,main="Residuals",col=col1)
    boxplot(object$residuals,main="Residuals")
    X<-object$effects
    if (length(table(y) > 6)) {colores = 1}
    else { colores = y + 1 }
    dev.interactive()
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
    if (!ask) {
      #dev.new()
      cc<-if(le>6) {cc<-c(ceiling(le/3), 3)}
      else {
        if (le>=2) cc<-c(ceiling(le/2), 2)
        else cc<-c(1,1)
      }
      par(mfrow = cc)
      #for (kk in 1:(cc[1]*cc[2])){
      # plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
      #}
    }
    else {
      par(mfrow=c(1,1))
      dev.interactive()
      oask <- devAskNewPage(TRUE)
      on.exit(devAskNewPage(oask))
    }
    if (!is.null(selec)) {
      par(mfrow=c(1,1))
      plot(X[, selec], eta, col = col1, ylab = "Linear predictors",
           xlab = paste("f(", namesres[selec], ")", sep = ""),
           main = paste(namesres[selec], "edf:", round(object$eqrank[selec],1)))
      if (length(table(y)) == 2) {
        abline(h = 0); abline(v = 0)
      }
    }
    else {
      for (i in 1:(ncol(X)-1)) {
        plot(X[, i], eta, col = col1, ylab = "Linear predictors",
             xlab = paste("f(", namesres[i], ")", sep = ""),
             main = paste(namesres[i], "edf:", round(object$eqrank[i],1)))
        if (length(table(y)) == 2) {
          abline(h = 0);abline(v = 0)
        }
      }
    }
    par(mfrow=c(1,1))
  }
  cat("\n")
  return(invisible(list("Influence"=influence,"object"=object)))
}

#' @export 
print.fregre.gkam<-function(x,digits = max(3, getOption("digits") - 3),...){
  object<-x
  print(object$family,"\n")
  # cat("Algorithm converged?",ifelse(object$converged,"Yes","No")," Number  of iterations ",object$iter,"\n")
  res<- object$result
  namesres<-names(res)
  lenres<-length(res)
  x<-object$fdataobj$data
  t=object$fdataobj$argvals
  y<-object$y
  if (object$family$family=="binomial") y<-as.numeric(factor(y,labels=c(0,1)))-1
  
  nobs <- n<-length(y)
  lenres<-length(object$result)
  tab<-matrix(NA,nrow=lenres,ncol=3)
  rownames(tab)<-paste("f(",namesres,")",sep="")
  colnames(tab)<-c("h","cor(f(X),eta)","edf")
  cat("alpha=",as.numeric(format(object$alpha,digits=3)),"  n= ",n,"\n")
  cat("****    ****    ****    ****    ****    ****\n")
  for (i in 1:lenres) {
    hh<-res[[i]]$h.opt
    tab[i,1]<-as.numeric(format(hh,digits=3))
    tab[i,2]<-round(cor(object$effects[,i],object$fitted.values),3)
    tab[i,3]<-round(object$eqrank[i],1)
  }
  print(tab)
  cat("****    ****    ****    ****    ****    ****\n")
  cat("edf: Equivalent degrees of freedom\n")
  residual.df <- nobs - sum(object$eqrank)
  w<-object$prior.weights
  ycen = y - mean(y)
  r.sq <- round(1 - var(w * (y - object$fitted.values))*(nobs - 1)/(var(w * y)*residual.df),3)
  r.sq2 <- round( 1 - sum((object$fitted.values-y)^2)/sum(ycen^2),3)
  #    cat("Residual deviance= ",round(object$deviance,3),"  Null deviance= ",round(object$null.deviance,3),"\n")
  dev.expl<-100*round((object$null.deviance - object$deviance)/object$null.deviance,3)
  #     cat("R-sq. = ",1 - sum((object$residuals)^2)/sum(ycen^2),"  Deviance explained = ", dev.expl,"%","\n")
  cat("AIC=",round(object$aic,3),"  Deviance explained =", dev.expl,"%","\n")
  cat("R-sq.=",r.sq2,"  R-sq.(adj)=",r.sq)
  cat("\n")
  return(invisible(object))
}



kgam.H<-function(object,inverse="svd") {
  #print("kgam.H")
  lenH<-length(object)
  if (lenH==1) return(object[[1]]$H)
  else {
    SS.list<-SS2<-list()
    n<-ncol(object[[1]]$H)
    SS<-matrix(NA,ncol=(lenH)*n,nrow=(lenH)*n)
    II=diag(n)
    unos<-matrix(1,ncol=n,nrow=n)/n
    M<-object[[1]]$H
    # MM=sweep(M,1,apply(M,1,mean),"-")
    MM=sweep(M,1,rowMeans(M),"-")
    if (lenH>1) {
      for (i in 2:lenH) {
        #   MMaux=sweep(object[[i]]$H,1,apply(object[[i]]$H,1,mean),"-")
        MMaux=sweep(object[[i]]$H,1,rowMeans(object[[i]]$H),"-")
        MM<-rbind(MM,MMaux)
      }
    }
    MM1<-matrix(rep(MM,lenH),nrow=(n*lenH))
    DD<-kronecker(diag(lenH),outer(rep(1,n),rep(1,n)))
    D1<-abs(DD-1)
    SS<-MM1*D1+diag(n*lenH)
    #  cat("kappa=",kappa(SS),"\n")    
    slv <- try(solve(SS),silent=TRUE)
    #    print(slv)
    if (is(slv,"try-error")) {
      sv<-svd(SS)
      slv<-drop((sv$v%*%diag(1/sv$d)%*%t(sv$u)))
      warning("Inverse of sigma computed by SVD")
    }    
    SSinv<-switch (inverse,solve=slv,
                   svd={
                     res.X=svd(SS)
                     (res.X$v)%*%diag(res.X$d^(-1))%*%t(res.X$u)})
    H2=SSinv%*%MM
    HH<-unos+H2[1:n,1:n]
    if (lenH>1) {
      for (i in 2:lenH) {HH=HH+H2[(n*(i-1)+1):(n*i),]  }
    }
    HH
  }
}