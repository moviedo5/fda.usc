#' @title Predict method for functional linear model (fregre.fd class)
#' 
#' @description Computes predictions for regression between functional explanatory variables
#' and scalar response using: basis representation, Principal Components
#' Analysis, Partial least squares or nonparametric kernel estimation.
#' 
#' Predicts from a fitted \code{fregre.basis} object,see
#' \code{\link{fregre.basis}} or \code{\link{fregre.basis.cv}}\cr Predicts from
#' a fitted \code{fregre.pc} object,see \code{\link{fregre.pc}} or
#' \code{\link{fregre.pc.cv}}\cr Predicts from a fitted \code{fregre.pls}
#' object,see \code{\link{fregre.pls}} or \code{\link{fregre.pls.cv}}\cr
#' Predicts from a fitted \code{fregre.np} object, see \code{\link{fregre.np}}
#' or \code{\link{fregre.np.cv}}.
#' 
#' @param object \code{fregre.fd} object.
#' @param new.fdataobj New functional explanatory data of \code{fdata} class.
#' @param se.fit =TRUE (not default) standard error estimates are returned for
#' each prediction.
#' @param scale Scale parameter for std.err. calculation.
#' @param df Degrees of freedom for scale.
#' @param interval Type of interval calculation.
#' @param level Tolerance/confidence level.
#' @param pred.var the variance(s) for future observations to be assumed for
#' prediction intervals. See \code{link{predict.lm}} for more details.
#' @param weights variance weights for prediction. This can be a numeric vector
#' or a one-sided model formula. In the latter case, it is interpreted as an
#' expression evaluated in newdata
#' @param \dots Further arguments passed to or from other methods.
#' @return If \code{se.fit = FALSE}, a vector of predictions of scalar response
#' is returned or a matrix of predictions and bounds with column names fit,
#' lwr, and upr if interval is set.
#' If \code{se.fit =TRUE} a list with the following components is returned: 
#' \itemize{
#' \item  \code{fit}: A vector of predictions or a matrix of predictions and bounds as above.
#' \item  \code{se.fit}: Associated standard error estimates of predictions.
#' \item  \code{residual.scale}: Residual standard deviations.
#' \item  \code{df}: Degrees of freedom for residual.
#' }
# If \code{se.fit} is TRUE then a 2 item list is returned with items (both
# arrays) fit and se.fit containing predictions and associated standard error
# estimates, otherwise an array of predictions of scalar response is returned.

#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' @seealso See Also as: \code{\link{fregre.basis}},
#' \code{\link{fregre.basis.cv}}, \code{\link{fregre.np}},
#' \code{\link{fregre.np.cv}}, \cr \code{\link{fregre.pc}},
#' \code{\link{fregre.pc.cv}}, \code{\link{fregre.pls}},
#' \code{\link{fregre.pls.cv}} \cr and \code{\link{summary.fregre.fd}}.\cr
#' @references Cai TT, Hall P. 2006. \emph{Prediction in functional linear
#' regression}. Annals of Statistics 34: 2159-2179.
#' 
#' Cardot H, Ferraty F, Sarda P. 1999. \emph{Functional linear model}.
#' Statistics and Probability Letters 45: 11-22.
#' 
#' Ferraty, F. and Vieu, P. (2006). \emph{Nonparametric functional data
#' analysis.} Springer Series in Statistics, New York.
#' 
#' Hall P, Hosseini-Nasab M. 2006. \emph{On properties of functional principal
#' components analysis}. Journal of the Royal Statistical Society B 68:
#' 109-126.
#' 
#' Hardle, W. \emph{Applied Nonparametric Regression}. Cambridge University
#' Press, 1994.
#' 
#' Ramsay, James O., and Silverman, Bernard W. (2006), \emph{ Functional Data
#' Analysis}, 2nd ed., Springer, New York.
#' 
#' Febrero-Bande, M., Oviedo de la Fuente, M. (2012).  \emph{Statistical
#' Computing in Functional Data Analysis: The R Package fda.usc.} Journal of
#' Statistical Software, 51(4), 1-28. \url{https://www.jstatsoft.org/v51/i04/}
#' @keywords regression
#' @examples 
#' \dontrun{
#' data(tecator)
#' absorp=tecator$absorp.fdata
#' ind=1:129
#' x=absorp[ind,]
#' y=tecator$y$Fat[ind]
#' newx=absorp[-ind,]
#' newy=matrix(tecator$y$Fat[-ind],ncol=1)
#' ## Functional PC regression
#' res.pc=fregre.pc(x,y,1:6)
#' pred.pc=predict(res.pc,newx)
#' # Functional PLS regression
#' res.pls=fregre.pls(x,y,1:6)
#' pred.pls=predict(res.pls,newx)
#' # Functional nonparametric regression
#' res.np=fregre.np(x,y,Ker=AKer.tri,metric=semimetric.deriv)
#' pred.np=predict(res.np,newx)
#' # Functional regression with basis representation
#' res.basis=fregre.basis.cv(x,y)
#' pred.basis=predict(res.basis[[1]],newx)
#'  
#' dev.new()
#' plot(pred.pc-newy)
#' points(pred.pls-newy,col=2,pch=2)
#' points(pred.np-newy,col=3,pch=3)
#' points(pred.basis-newy,col=4,pch=4)
#' sum((pred.pc-newy)^2,na.rm=TRUE)/sum((newy-mean(newy))^2,na.rm=TRUE)
#' sum((pred.pls-newy)^2,na.rm=TRUE)/sum((newy-mean(newy))^2,na.rm=TRUE)
#' sum((pred.np-newy)^2,na.rm=TRUE)/sum((newy-mean(newy))^2,na.rm=TRUE)
#' sum((pred.basis-newy)^2,na.rm=TRUE)/sum((newy-mean(newy))^2,na.rm=TRUE)
#' }
#' 
#' @export
predict.fregre.fd<-function(object,new.fdataobj=NULL,se.fit=FALSE,scale = NULL,df=df,
    interval = "none", level = 0.95,weights = 1, pred.var = res.var/weights, ...){
if (is.null(object)) stop("No fregre.fd object entered")
if (is.null(new.fdataobj)) return(object$fitted.values)
if (object$call[[1]]=="gam")  return(predict(object,new.fdataobj,...))
if (object$call[[1]]=="lm")  return(predict(object,new.fdataobj,...)  )
if (object$call[[1]]=="glm")  return(predict(object,new.fdataobj,...))
if (!is.fdata(new.fdataobj)) 
  new.fdataobj=fdata(new.fdataobj,object$fdataobj[["argvals"]],object$fdataobj[["rangeval"]],object$fdataobj[["names"]])
y=object$y
isfdata<-is.fdata(y)
if (!isfdata) {
   if (is.vector(y)) y.mat<-matrix(y,ncol=1)
   }
else y.mat<-y$data
gg<-1:nrow(new.fdataobj)
nas<- is.na.fdata(new.fdataobj)
if (any(nas)) {
   bb<-!nas
   cat("Warning: ",sum(nas)," curves with NA are omited\n")
   new.fdataobj$data<-new.fdataobj$data[bb,]
   gg<-gg[bb]
   }
newx<-new.fdataobj[["data"]]
tt<-new.fdataobj[["argvals"]]
rtt<-new.fdataobj[["rangeval"]]
nn <- nrow(new.fdataobj)
np <- ncol(new.fdataobj)
 if (is.null(rownames(newx)))         rownames(newx) <- 1:nn
 if (object$call[[1]]=="fregre.pc" || object$call[[1]]=="fregre.ppc") {
  Z<- inprod.fdata(fdata.cen(new.fdataobj,object$fdata.comp$mean)[[1]],object$fdata.comp$rotation[object$l])
  colnames(Z)<-names(object$lm$coefficients[-1])
  XX<-data.frame(Z)   
  res.var<-object$sr2 
   if (object$lambda==0) return(predict.lm(object=object$lm,newdata=XX,se.fit=se.fit,
     interval=interval,level=level,scale=scale,df=df,weights =weights,pred.var=pred.var,...))
   else {  
    a1<-object$coefficients[1]*rep(1,len=nrow(newx))
    object$beta.est$data<-matrix(object$beta.est$data,nrow=1)
    b1<-inprod.fdata(fdata.cen(new.fdataobj,object$fdata.comp$mean)[[1]],object$beta.est)#/(ncol(newx)-1)
    yp<- a1+b1      
    yp<-drop(yp)
    XX2<-cbind(rep(1,len=nn),Z)
  names(yp) <-gg                      
  predictor<-yp 
  if (se.fit || interval != "none") {  
    XX2<-cbind(rep(1,len=nn),Z) 
    ip<-rowSums((XX2 %*%object$Vp*XX2))   
#    res.var<-sum(fit$residuals^2)/fit$df.residual      
    df<-object$lm$df.residual
    if (interval != "none") {
        tfrac <- qt((1 - level)/2, df)
        hwid <- tfrac * switch(interval, confidence = sqrt(ip),prediction = sqrt(ip + pred.var))
        predictor <- cbind(predictor, predictor + hwid %o%c(1, -1))
        colnames(predictor) <- c("fit", "lwr", "upr")
    }
   }
  if (se.fit)   {
        se <- sqrt(ip)
        return(list(fit = predictor, se.fit = se, df = df, residual.scale = sqrt(res.var)))
        }
else return(predictor)
    }
 }
else{
  if (object$call[[1]]=="fregre.pls" || object$call[[1]]=="fregre.ppls") {
  newXcen<-fdata.cen(new.fdataobj,object$fdata.comp$mean)[[1]]
  if (object$fdata.comp$norm)  {
    sd.X <- sqrt(apply(object$fdataobj$data, 2, var))
    newXcen$data<- newXcen$data/(rep(1, nn) %*% t(sd.X))
  }
  Z<- inprod.fdata(newXcen,object$fdata.comp$rotation) 
  colnames(Z)<-names(object$lm$coefficients[-1])  
  XX<-data.frame(Z)  
  return(predict.lm(object=object$lm,newdata=XX,se.fit=se.fit,...))
 } 
else {
 if (object$call[[1]]=="fregre.basis" || object$call[[1]]=="fregre.basis.cv" ){
  x=newx
  basis.x=object$basis.x.opt            
  xcen<-fdata.cen(new.fdataobj,object$mean)[[1]]
	x.fd=Data2fd(argvals=tt,y=t(xcen$data),basisobj=basis.x)
  C=t(x.fd$coefs)
  if (is.vector(object$b.est)) object$b.est<-matrix(object$b.est,ncol=1,nrow=length(object$b.est))
  Z<-C%*%object$J 
  yp=object$a.est* rep(1,len=nn) + Z%*%object$b.est
   if (isfdata) {    return(fdata(yp,y$argvals,y$rtt,y$names))    }
  predictor<-yp
  if (se.fit || interval != "none") {    
    XX2<-cbind(rep(1,len=nn),Z) 
    ip<-(rowSums((XX2 %*%object$Vp*XX2)))
    res.var<-object$sr2    
    df<-object$lm$df.residual   
    if (interval != "none") {
        tfrac <- qt((1 - level)/2, df)
        hwid <- tfrac * switch(interval, confidence = sqrt(ip),prediction = sqrt(ip + pred.var))
        predictor <- cbind(predictor, predictor + hwid %o%c(1, -1))
        colnames(predictor) <- c("fit", "lwr", "upr")     }
   }
  if (se.fit)   {
        se <- sqrt(ip)
        return(list(fit = predictor, se.fit = se, df = df, residual.scale = sqrt(res.var)))
        }
else return(predictor)  
 }
 else {
 if (object$call[[1]]=="fregre.np" || object$call[[1]]=="fregre.np.cv"){
 x=object$fdataobj
 h=object$h.opt
 n = nrow(x)
 nn = nrow(newx)
 np <- ncol(x)
 if (is.null(rownames(newx)))         rownames(newx) <- 1:nn
 par.S<-object$par.S
 bs=as<-list() 
 Ker=object$Ker
   par.metric<-attr(object$mdist,"par.metric")
   par.metric[["fdata1"]]<-new.fdataobj
   par.metric[["fdata2"]]<-x
  a1<-attr(object$mdist,"call")
  a2<-attr(object$par.S,"call")
  nmdist <- do.call(a1,par.metric)
  par.S$tt<-nmdist
  par.S$cv=FALSE
  H<-do.call(a2,par.S)
  yp=H%*%y.mat    
 if (isfdata) {    return(fdata(yp,y$argvals,y$rtt,y$names))       }
 else{ 
  predictor<-drop(yp)
  if (se.fit || interval != "none") {    
    ip<- drop(H%*%y.mat^2-yp^2  )
    df<-fdata.trace(object$H)
    edf<-n-df
    res.var<-object$sr2    
    se<-sqrt(ip/edf)  
    if (interval != "none") {
        tfrac <- qt((1 - level)/2, df)
        hwid <- tfrac * switch(interval, confidence = se,prediction = sqrt((ip + pred.var)/edf))
        predictor <- cbind(predictor, predictor + hwid %o%c(1, -1))
        colnames(predictor) <- c("fit", "lwr", "upr")
      }
   }
  if (se.fit)   {           
        return(list(fit = predictor, se.fit = se, df = df, residual.scale = sqrt(res.var)))
        }
else return(predictor)      
 }   
 }     }
}     }
return(yp)
}

