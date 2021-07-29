#' @title Compute fucntional coefficients from functional data represented in a base of functions
#' 
#' @aliases fdata2basis summary.basis.fdata
#' @description Compute fucntional coefficients from functional data (\code{\link{fdata}} class object) 
#' represented in a basis (fixed of data-driven basis).
#' 
#' @param fdataobj \code{\link{fdata}} class object.
#' @param basis  	a functional basis object defining the basis
#' @param method character string, if it is "grid" the fdata object is evaluated in the grid (\code{argvals} of fdata),
#'  if it is "inprod" the basis representation of functional data is computed by inner product 
#'  (\code{\link{inprod.fdata}(fdataobj,basis)}).
#' @param object \code{basis.fdata} class object calculated by: \code{\link{fdata2basis}}
#' @param draw logical, original curves  and their basis representation are plotted
#' @param index vector, by default (if NULL) the first n curves are plotted, where n = min(4, length(fdataobj)). 
#' Otherwise, index vector indicates taht curvesare plotted.
#' @param \dots Further arguments passed to or from other methods.
#' @return fdata2basis \code{fdata2bais} function return: 
#' \itemize{
#' \item {coef}{a matrix or two-dimensional array of coefficients.}
#' \item {basis}{basis of \code{\link{fdata}} class evaluated in the same grid of \code{fdataobj}. } 
#' }
#' 
#' summary function return: 
#' \itemize{
#' \item {R}{a matrix with a measure similar to R-sq for each curve aproximation (by row) and number of basis elements (by column).}
#' }
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente  \email{manuel.oviedo@@usc.es}
#' @seealso  Inverse function: \code{\link{gridfdata}}.
#' Alternative method: \code{\link{fdata2pc}}, \code{\link{fdata2pls}}.
# @references
#' @keywords multivariate
#' @examples 
#' \dontrun{
#' T <- 71
#' S <- 51
#' tj <- round(seq(0,1,len=T),3)
#' si <- round(seq(0,1,len=S),3)
#' beta1 <- outer(si,tj,function(si,tj){exp(-5*abs((tj-si)/5))})
#' nbasis.s =7
#' nbasis.t=11
#' base.s <- create.fourier.basis(c(0,1),nbasis=nbasis.s)
#' base.t <- create.fourier.basis(c(0,1),nbasis=nbasis.t)
#' y1 <- fdata(rbind(log(1+tj),1-5*(tj-0.5)^2),argvals=tj,rangeval=c(0,1))
#' aa <- fdata2basis(y1,base.t,method="inprod")
#' summary(aa)
#' plot(gridfdata(aa$coefs,aa$basis))
#' lines(y1,lwd=2,col=c(3,4),lty=2)
#' }
#' @export
fdata2basis <- function(fdataobj, basis, method=c("grid","inprod")){
  xmean <- NULL
  
  
  if (is.basis(basis)){
   # print(1)
      bb=fdata(t(eval.basis(fdataobj$argvals,basis)),
               argvals=fdataobj$argvals,rangeval=fdataobj$rangeval)
  } else{
    if (class(basis) %in% c("fdata")){
      bb=basis
      xcen <- fdata.cen(fdataobj)
      fdataobj <- Xcen$Xcen
      xmean <- Xcen$meanX
      #xmean <- func.mean(fdataobj)
  } else {
    bb=  basis$basis   
    xmean <- basis$mean
    fdataobj$data <- sweep(fdataobj$data,2,xmean$data,"-")
  }
  }
  #print(2)
  if (method[1]=="grid"){
    #print(3)
  	A=t(fdataobj$data)
	  B=t(bb$data)
  	coefs=t(Minverse(t(B)%*%B)%*%t(B)%*%A)
  } else {
   # print(4)
  	A=t(inprod.fdata(fdataobj,bb))
  	B=inprod.fdata(bb)
  	coefs=t(solve(B)%*%A)
  }
  #print(5)
  if (!is.null(basis$type)) type <- basis$type
  #if (!is.null(bb$names$main)) type <- bb$names$main
  out <- list(coefs=coefs,basis=bb, 
              fdataobj=fdataobj,
              type=type,
              mean = xmean
              )
  class(out) <- "basis.fdata"
  return(out)
}

#' @rdname fdata2basis
#' @export
summary.basis.fdata=function(object, draw=TRUE, index=NULL,...) {
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

#wrapper of gridfdata
basis2fdata<-function (coef, fdataobj, mu) 
{
  return(gridfdata(coef, fdataobj, mu) )
}

# Devolver la matriz entera R2 para todas las curvas en un summary 
# parece excesivamente prolijo. Como texto habría que devolver un R2 
# conjunto  con el tamaño global de la base:
#   1-sum(norm.fdata(.,.)^2)/sum(norm.fdata(Xcen)^2). 
# La función debe devolver la matriz R2 como invisible por si 
# se quiere usar posteriormente pero devolver el mismo objeto 
# que entra no tiene mucha utilidad. Respecto al texto que se 
# incluye al principio (solo la llamada de la función) no sé 
# si sería mejor poner un texto con el object$type, la longitud 
# y alguna info sobre los argvals y el rangeval.  


# no ayuda por el momento
fdata2basis2d=function(fdata2d,basis.s,basis.t){
  n=dim(fdata2d$data)[1]
  if (is.basis(basis.s)) {
    basis.s=fdata(t(eval.basis(fdata2d$argvals[[1]],basis.s)),
                  argvals=fdata2d$argvals[[1]],rangeval=fdata2d$rangeval[[1]])
  }  
  if (is.basis(basis.t)) {
    basis.t=fdata(t(eval.basis(fdata2d$argvals[[2]],basis.t)),
                  argvals=fdata2d$argvals[[2]],rangeval=fdata2d$rangeval[[2]])
    }
  n.s=nrow(basis.s)
  n.t=nrow(basis.t)
  eta=basis.s$data;theta=basis.t$data
  Jeta=tcrossprod(eta);Jtheta=tcrossprod(theta)
  AA=Minverse(Jeta)%*%eta
  BB=t(theta)%*% Minverse(Jtheta)
  B=array(NA,dim=c(n,n.s,n.t))
  #B=solve(Jeta)%*%eta%*%fdata2d%*%t(theta)%*%solve(Jtheta)
  for (i in 1:n){B[i,,]=AA%*%fdata2d$data[i,,]%*%BB}
  return(list(coefs=B,b.s=basis.s,b.t=basis.t))
}

# no ayuda por el momento
basis2d2fdata=function(coefs,b.s,b.t){
  if (!is.array(coefs)) 
    stop("coefs must be an array")
  if (dim(coefs)[2]!=length(b.s) | dim(coefs)[3]!=length(b.t)) 
    stop("Dimensions of coefs and basis.s or basis.t do not match")
  n=dim(coefs)[1]
  n.s=length(b.s$argvals)
  n.t=length(b.t$argvals)
  Beta=array(NA,dim=c(n,n.s,n.t))
  for (i in 1:n){
   Beta[i,,]=t(b.s$data)%*%coefs[i,,]%*%b.t$data
  }
  Beta=fdata(Beta,argvals=list(b.s$argvals,b.t$argvals),
             rangeval=list(b.s$rangeval,b.t$rangeval),
  		       names=list(main="B(s,t)",xlab=list("s","t")),
  		       fdata2d=TRUE)
}

# no ayuda
intXbeta=function(X,Beta,arg.t=1:ncol(Beta),
                  rangeval=range(arg.t),
                  equi=argvals.equi(X$argvals),
                  method=NULL){
  if (is.null(method)) {
    par.fda.usc=eval(parse(text="fda.usc:::par.fda.usc"), envir=.GlobalEnv)
    method=par.fda.usc$int.method
  }
  if (inherits(X,"fdata")) {
  	nX=NROW(X)
  	arg.s=X$argvals
  	equi=argvals.equi(arg.s)
  	n.s=length(X$argvals)
  	if (!is.matrix(Beta)) { stop("Only implemented for Beta being a matrix")} else {
  		if (nrow(Beta)!=n.s) stop ("Dimensions of X and Beta do not coincide")
  		if (ncol(Beta)!=length(arg.t)) stop("Beta and argvals.t do not have the same length")
  			n.t=ncol(Beta);arg.t=arg.t
  						} 
  	if (equi) {
  	  res=X$data%*%Beta
  	} else{
    	res=matrix(NA,nrow=nX,ncol=n.t)
    	for (i in 1:nX){
    		for (j in 1:n.t){
    	    res[i,j]=int.simpson2(arg.s,X$data[i,]*Beta[,j],equi=equi,method=method)
    	}}
  	}
  	nms <- list(main="X(s)Beta(s,t)",xlab="t",ylab="Y(t)")
  	res=fdata(res,argvals=arg.t,rangeval=rangeval,names=nms)
  }
}

# gridfdata2d=function(coefs,fdataobj,mu){
# if (!is.fdata(fdataobj) | length(fdataobj$argvals)!=2) stop("fdataobj is not an fdata2d object")

# if (is.matrix(coefs)) {
		# n=1
		# if (nrow(coefs)!=length(fdataobj$argvals[[1]]) | ncol(coefs)!=length(fdataobj$argvals[[2]])) stop("Dimensions of coefs and fdataobj do not coincide")
		# n.s=nrow(coefs);n.t=ncol(coefs)
		# coefs=array(coefs,dim=c(1,nrow=n.s,ncol=n.t))
	# } else if (is.array(coefs)){
		# n=dim(coefs)[1]
		# if (dim(coefs)[2]!=length(fdataobj$argvals[[1]]) | dim(coefs)[3]!=length(fdataobj$argvals[[2]])) stop("Dimensions of coefs and fdataobj do not coincide")
		# n.s=dim(coefs)[2]
		# n.t=dim(coefs)[3]
# } else stop("coefs must be a matrix or an array")
# res=array(NA,dim=c(


# if (!missing(mu)) {mu=fdata(array(rep(0,n.s*n.t),dim=c(1,n.s,n.t)),argvals=fdataobj$argvals, 
						# rangeval=fdataobj$rangeval, names=fdataobj$names,fdata2d=TRUE)}
		
# }

#' Minverse <- fda.usc.devel:::Minverse
#' T <- 71
#' S <- 51
#' tj <- round(seq(0,1,len=T),3)
#' si <- round(seq(0,1,len=S),3)
#' beta1 <- outer(si,tj,function(si,tj){exp(-5*abs((tj-si)/5))})
#' nbasis.s =31
#' nbasis.t=51
#' base.s <- create.fourier.basis(c(0,1),nbasis=nbasis.s)
#' base.t <- create.fourier.basis(c(0,1),nbasis=nbasis.t)
#----------------------------------------------
#' y1 <- fdata(rbind(log(1+tj),1-5*(tj-0.5)^2),argvals=tj,rangeval=c(0,1))
#' aa <- fdata2basis(y1,base.t,method="inprod")
#' plot(gridfdata(aa$coefs,aa$basis))
#' lines(y1,lwd=2,col=c(3,4),lty=2)
#----------------------------------------------
#' z <- array(NA,dim=c(2,S,T))
#' for (i in 1:2) {z[i,,] <- outer(si, tj, function(u, v) {(i*u)*(v-0.5)^i})}
#' z.fdata <- fdata(z,argvals=list(si,tj),rangeval=list(c(0,1),c(0,1)),fdata2d=TRUE)
#' aa <- fdata2basis2d(z.fdata,base.s,base.t)
#' z2.fdata <- basis2d2fdata(aa$coefs,aa$b.s,aa$b.t)
#' plot(z.fdata)   # Original
#' plot(z2.fdata)  #Reconstruida
#' par(mfrow=c(1,2))
#' plot(z.fdata,type="persp")
#' plot(z2.fdata,type="persp")
#-------------------------------------------



#-------------------------------------------
#' xx=fdata(rbind(-1/((si-0.5)^2+0.25),si,si^2),argvals=si,rangeval=c(0,1))
#' miBeta=outer(si,tj,function(u,v){-(v-0.25)^2})
#' miBeta=fdata(miBeta,argvals=list(si,tj),rangeval=list(c(0,1),c(0,1)),fdata2d=TRUE)
#' yy=intXbeta(xx,miBeta$data[1,,],tj,rangeval=c(0,1))
#' par(mfrow=c(1,2))
#' plot(xx);plot(yy)

# eta=t(eval.basis(si,base.s))
# theta=t(eval.basis(tj,base.t))

# Jeta=tcrossprod(eta)
# Jtheta=tcrossprod(theta)

# B=solve(Jeta)%*%eta%*%beta1%*%t(theta)%*%solve(Jtheta)

# betatilde=t(eta)%*%B%*%theta
# par(mfrow=c(1,2))

# image(si,tj,beta1,main="beta")
# image(si,tj,betatilde,main="betatilde")
#--------------------------
# npc.x=5
# npc.y=5

# Sigma.X=outer(si,si,function(u,v){exp(-abs(u-v))})
# Sigma.Y=outer(tj,tj,function(u,v){exp(-abs(u-v))})
# res.X=eigen(Sigma.X)
# bpc.X=fdata(t(res.X$vector[,1:npc.x]),argvals=si)
# bpc.X[["data"]]=sweep(bpc.X[["data"]],1,norm.fdata(bpc.X),"/")

# res.Y=eigen(Sigma.Y)
# bpc.Y=fdata(t(res.Y$vector[,1:npc.y]),argvals=tj)
# bpc.Y[["data"]]=sweep(bpc.Y[["data"]],1,norm.fdata(bpc.Y),"/")


# X=rproc2fdata(500,t=si,sigma="vexponential")
# Y=rproc2fdata(500,t=tj,sigma="vexponential")
# bpc.X=create.pc.basis(X,1:npc.x,lambda=10)$basis
# bpc.Y=create.pc.basis(Y,1:npc.y,lambda=10)$basis
# eta.pc=bpc.X$data
# theta.pc=bpc.Y$data
# Jeta.pc=tcrossprod(eta.pc)
# Jtheta.pc=tcrossprod(theta.pc)

# B2=solve(Jeta.pc)%*%eta.pc%*%beta1%*%t(theta.pc)%*%solve(Jtheta.pc)
# betatilde2=t(eta.pc)%*%B2%*%theta.pc
# par(mfrow=c(1,2))
# image(si,tj,beta1,main="beta")
# image(si,tj,betatilde2,main="betatilde2")