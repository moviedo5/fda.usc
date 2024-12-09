predict_mregre <- function(object, newdata=NULL, se.fit=FALSE,
                           scale = NULL, df=df, interval = "none",
                           level = 0.95,weights = 1,
                          #  pred.var = res.var/weights,
                          ...){
if (is.null(object)) stop("No fregre.fd object entered")
if (is.null(newdata)) return(object$fitted.values)
if (object$call[[1]]=="gam")  return(predict(object,newdata,...))
if (object$call[[1]]=="lm")  return(predict(object,newdata,...)  )
if (object$call[[1]]=="glm")  return(predict(object,newdata,...))
#if (!is.fdata(newdata)) 
#  newdata=fdata(newdata,object$fdataobj[["argvals"]],object$fdataobj[["rangeval"]],object$fdataobj[["names"]])
y <- object$y
if (is.vector(y)) y.mat<-matrix(y,ncol=1)
if (is.vector(newdata)) newdata<-matrix(newdata,ncol=1)
gg<-1:nrow(newdata)
# nas<-apply(newdata,1,count.na)
nas<-apply(newdata,1,anyNA)
if (any(nas)) {
   bb<-!nas
   cat("Warning: ",sum(nas)," curves with NA are omited\n")
   newdata<-newdata[bb,]
   gg<-gg[bb]
   }
newx<-newdata
nn <- NROW(newdata)
np <- NCOL(newdata)
 if (is.null(rownames(newx)))         rownames(newx) <- 1:nn
 if (object$call[[1]]=="mregre.np" || object$call[[1]]=="mregre.np.cv"){
 x=object$x
 h=object$h.opt
 n = nrow(x)
 nn = nrow(newx)
 np <- ncol(x)
 if (is.null(rownames(newx)))         rownames(newx) <- 1:nn
 par.S<-object$par.S
 Ker <- object$Ker
 par.metric<-attr(object$mdist,"par.metric")
 par.metric[["x"]] <- newdata
 par.metric[["y"]] <- x
 a1 <- attr(object$mdist,"call")
 a2 <- attr(object$par.S,"call")
 par.S$tt <- do.call(a1,par.metric)
 par.S$cv <- FALSE
 H <- do.call(a2,par.S)
 yp <- H %*% y.mat    
 }   
return(yp)
}
