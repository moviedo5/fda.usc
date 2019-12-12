################################################################################
predict.mregre <- function(object, newdata=NULL, se.fit=FALSE,
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
# print("******************entra preict mgregre***********************")
# print("newdata")
# print(newdata[1:3])
  y=object$y
if (is.vector(y)) y.mat<-matrix(y,ncol=1)
if (is.vector(newdata)) newdata<-matrix(newdata,ncol=1)
gg<-1:nrow(newdata)
#nas<-apply(newdata,1,count.na)
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
 bs=as<-list()
 Ker=object$Ker
   par.metric<-attr(object$mdist,"par.metric")
   par.metric[["x"]]<-newdata
# print(dim(x))
# print(dim(newdata))   
# print("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa")
# print(class(x))
# print(class(newdata))   
   par.metric[["y"]]<-x
  a1<-attr(object$mdist,"call")
  a2<-attr(object$par.S,"call")
# print(a1)  
  ########$par.metric$x<-object$x
  nmdist <- do.call(a1,par.metric)
  # print(object$mdist[1:2,1:3])
  # print(nmdist[1:2,1:3])
  # print(a1)
  nmdist <- do.call(a1,par.metric)
  par.S$tt<-nmdist
  # print(dim(nmdist))  
  # print(attributes(nmdist))  
  par.S$cv=FALSE
  #do.call(ty,par.S)
  H<-do.call(a2,par.S)
  yp=H%*%y.mat    
  # print(names(par.metric))
 }   
# print("saaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaleeeeeeeeeeeeeeeeee")
return(yp)
}
