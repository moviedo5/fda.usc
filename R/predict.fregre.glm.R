# @S3method predict fregre.glm
#' @rdname predict.fregre.lm
#' @export 
predict.fregre.glm<-function(object, newx = NULL, type = "response",...){
 if (is.null(object)) stop("No fregre.glm object entered")
 if (is.null(newx)) {
    if (type == "effects"){
      fake  = predict.glm(object, type = "terms", ...) 
      yp <- effect.fake(object,fake)
    } else{
      yp  = predict.glm(object, type = type, ...)
    }
  return(yp)
 } else {
# data=newx
 basis.x=object$basis.x
 basis.b=object$basis.b
 formula=object$formula
 tf <- terms.formula(formula)
 terms <- attr(tf, "term.labels")
# nt <- length(terms)
# pf <- rf <- "~"  	 
	 
##########
 vtab<-rownames(attr(tf,"factors"))
 vnf=intersect(terms,names(newx$df))
# vnf2=intersect(vtab[-1],names(data$df)[-1])
 vfunc2=setdiff(terms,vnf)
 vint=setdiff(terms,vtab)
 vfunc=setdiff(vfunc2,vint)
# vnf=c(vnf2,vint)
 off<-attr(tf,"offset")
 beta.l=list()
 kterms=1

if (length(vnf)>0) {
 first=FALSE
 XX=data.frame(newx[["df"]][,c(vnf)])
 names(XX)=vnf } else  first=TRUE

if (length(vfunc)>0)  {
#   yp2<-a1 <- object$coefficients[1] * rep(1, len = nrow(data[[vfunc[1]]]))
   for (i in 1:length(vfunc)) {
   if(inherits(newx[[vfunc[i]]],"fdata"))  {
      fdataobj<-newx[[vfunc[i]]]
      if (nrow(newx[[vfunc[i]]])==1) rwn<-NULL
      else rwn<-rownames(newx[[vfunc[i]]]$data)

	  xaux<-fdata2basis(newx[[vfunc[i]]],basis.x[[vfunc[i]]])
	  Z <- xaux$coefs%*%object$vs.list[[vfunc[i]]]
	  colnames(Z)<-paste(vfunc[i],".",colnames(object$vs.list[[vfunc[i]]]),sep="")

	  if (first) {XX=Z; first=FALSE} else {XX=cbind(XX,Z)}
	  } else {
	  if (inherits(newx[[vfunc[i]]],"fd")) {
            if (!inherits(object$basis.x[[vfunc[i]]], "pca.fd")) {
              x.fd <- newx[[vfunc[i]]]
              r = x.fd[["basis"]][["rangeval"]]
              J <- object$vs.list[[vfunc[i]]]
              x.fd$coefs <- x.fd$coefs - object$mean[[vfunc[i]]]$coefs[,1]
              Z = t(x.fd$coefs) %*% J
              colnames(Z) = colnames(J)
            }
            else {
              name.coef[[vfunc[i]]] = paste(vfunc[i], 
                                            ".", colnames(object$basis.x[[vfunc[i]]]$harmonics$coefs), 
                                            sep = "")
              newx[[vfunc[i]]]$coefs <- sweep(data[[vfunc[i]]]$coefs, 
                                              1, (object$basis.x[[vfunc[i]]]$meanfd$coefs), 
                                              FUN = "-")
              fd.cen <- newx[[vfunc[i]]]
              Z <- inprod(fd.cen, object$basis.x[[vfunc[i]]]$harmonics)
              colnames(Z) <- paste(vfunc[i],".",colnames(object$vs.list[[vfunc[i]]]),sep="")
            }
            if (first) {
              XX = Z
              first = FALSE
            }
            else XX = cbind(XX, Z)
          }
          else stop("Please, enter functional covariate")
        
	  }
	  }
	  }
	  
 if (first) return(rep(object$coefficients,length=nrow(newx[[1]])) )        
 if (!is.data.frame(XX)) XX=data.frame(XX)    
 if (type == "effects"){
   fake  = predict.glm(object, newdata = XX, type = "terms",x=TRUE,y=TRUE, ...)
   yp <- effect.fake(object,fake)
 } else{
   yp <- predict.glm(object = object, newdata = XX, type = type, x=TRUE,y=TRUE,...)
 }
 return(yp)
}
}
