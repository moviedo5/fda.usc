#' @rdname predict.fregre.lm
#' @export 
predict.fregre.gsam <- function(object, newx = NULL, type = "response",...){
    if (is.null(object)) stop("No fregre.gsam.fr object entered")
    if (is.null(newx)) {
      if (type == "effects"){
        fake  = predict.gam(object, type = "terms", ...) 
        yp <- effect.gam(object,fake)
      } else{
        yp  = predict.gam(object, type = type, ...)
      }
      return(yp)
    } else {
      #print(2)
      data=newx
      basis.x=object$basis.x
      formula=object$formula
      tf <- terms.formula(formula, specials = c("s", "te", "t2"))
      terms <- attr(tf, "term.labels")
      if (length(terms)==0) return(rep(object$coefficient,length=nrow(newx[[1]])) ) 
      special <- attr(tf, "specials")
      nt <- length(terms)
      if (attr(tf, "response") > 0) {
        response <- as.character(attr(tf, "variables")[2])
      }
      vtab <- rownames(attr(tf,"factors"))
      gp <- interpret.gam(formula)
      #print(gp$smooth.spec)
      len.smo <- length(gp$smooth.spec)
      name.coef <- NULL
      vfunc <- object$vfunc
      vnf <- object$vnf
      nnf<-length(vnf)
      if (!is.null(vnf)) {
        first=FALSE
        #XX=NULL
        XX=data.frame(data[["df"]][,c(vnf)])
        names(XX) <- vnf
      } else {  
        first=TRUE
      }
      lenfunc <- length(vfunc)
      bsp1 <- object$bsp
      if (lenfunc>0) {
        k <- 1
        mean.list=vs.list=JJ=list()
        for (i in 1:lenfunc) {
          if(class(newx[[vfunc[i]]])[1]=="fdata"){
            tt <- data[[vfunc[i]]][["argvals"]]
            rtt <- data[[vfunc[i]]][["rangeval"]]
            fdataobj <- data[[vfunc[i]]]
            fdat <- data[[vfunc[i]]]
            dat <- fdataobj$data
            if (nrow(dat) == 1)  rwn <- NULL         else rwn<-rownames(dat)
            #  if (basis.x[[vfunc[i]]]$type=="pc" 
            #      | basis.x[[vfunc[i]]]$type=="pls")
            #    bsp1=FALSE      else bsp1 <- TRUE
            xaux <- fdata2basis(data[[vfunc[i]]], basis.x[[vfunc[i]]],method="inprod")
            name.coef[[vfunc[i]]] <- colnames(xaux$coefs) <- paste(vfunc[i],".",colnames(xaux$coefs),sep="")
            Z <- xaux$coefs
            if ( object$bsp ){
              Z = Z %*% inprod(basis.x[[vfunc[i]]], object$basis.b[[vfunc[i]]])  
             }
            if (first) {
              XX=Z
              first=FALSE
            }   else {
              XX = cbind(XX,Z)} }
        }
      }
    if (!is.data.frame(XX)) 
      XX=data.frame(XX)
      names(XX)<-colnames(object$XX)[-1]
    yp=predict.gam(object=object, newdata=XX, type=type,...)
    return(yp)
    }
  }
  
effect.gam <- function(object,terms){
  #fake<-predict(object,type = "terms")
  fake <- terms
  ff<-colnames(fake)
  nc <- nchar(ff)
  ismo <- which(substr(ff,1,2)=="s(")
  ilin <- which(substr(ff,1,2)!="s(")
  nam.smo <- substr(ff[ismo],3,nc[ismo]-1)
  nam.lin <- ff[-ismo]
  #nam <- c(nam.lin,nam.smo)
  vfunc <- names(object$JJ)
  nfunc <- length(vfunc)
  effects <- NULL
  vf <- NULL
  for (i in 1:nfunc){
    ifunc <- colnames(object$JJ[[i]])
    #dfnames <- intersect(nam,ifunc)
    iss<-intersect(nam.smo,ifunc)
    if (length(iss)>0)  dfnames <- paste0("s(",ifunc,")")
    else  dfnames <- ifunc
    vf <- c(vf,dfnames)
    effects <- cbind(effects,rowSums(fake[,dfnames,drop=F]))
  }
  colnames(effects)<-vfunc
  effects <- cbind(fake[,setdiff(ff,vf),drop=F],effects)
  # colnames(effects) <- c(effects.df,vfunc )
  effects  
}
###############################
# Modification 20221007 (adapted from  fregre.gsam.fr)
# New code included:
# if ( object$bsp ){
#    Z = Z %*% inprod(basis.x[[vfunc[i]]], object$basis.b[[vfunc[i]]])  
#   colnames(Z)<-colnames(object$XX)[-1]
# }

