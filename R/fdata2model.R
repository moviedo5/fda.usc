inprodbasis<-function(basis1,basis2){
if (class(basis1)=="pca.fd") basis1$type="pca.fd"
if (class(basis2)=="pca.fd") basis2$type="pca.fd"
if (is.null(basis1$type) | is.null(basis2$type)) stop("Basis type not recognized")
if (basis1$type=="pc" | basis1$type=="pls"){
	nms1=rownames(basis1$basis$data)
	if (basis2$type=="pc" | basis2$type=="pls"){
	nms2=rownames(basis2$basis$data)
	J=inprod.fdata(basis1$basis,basis2$basis)
	}else{
	if (class(basis2)=="pca.fd") {
	baux2=fdata(t(eval.fd(basis1$basis$argvals,basis2$harmonics)),argvals=basis1$basis$argvals,rangeval=basis1$basis$rangeval)
	nms2=basis2$harmonics$fdnames[[2]]
		} else { 
	baux2=fdata(t(eval.basis(basis1$basis$argvals,basis2)),argvals=basis1$basis$argvals,rangeval=basis1$basis$rangeval)
	nms2=basis2$names
		}
	J=inprod.fdata(basis1$basis,baux2)
	}
} else if (basis2$type=="pc" | basis2$type=="pls"){
	nms2=rownames(basis2$basis$data)
	if (class(basis1)=="pca.fd") {
	baux1=fdata(t(eval.fd(basis2$basis$argvals,basis1$harmonics)),argvals=basis2$basis$argvals,rangeval=basis2$basis$rangeval)
	nms1=basis1$harmonics$fdnames[[2]]
		} else { 
	baux1=fdata(t(eval.basis(basis2$basis$argvals,basis1)),argvals=basis2$basis$argvals,rangeval=basis2$basis$rangeval)
	nms1=basis1$names
	}
	J=inprod.fdata(baux1,basis2$basis)
} else {
	if (class(basis1)=="pca.fd") {baux1=basis1$harmonics;nms1=basis1$harmonics$fdnames[[2]]} else {baux1=basis1;nms1=basis1$names}
	if (class(basis2)=="pca.fd") {baux2=basis2$harmonics;nms2=basis2$harmonics$fdnames[[2]]} else {baux2=basis2;nms2=basis2$names}
	J=inprod(baux1,baux2)
}
rownames(J)=nms1
colnames(J)=nms2
return(J)
}


fdata2model <- function(vfunc, vnf, response, data,
                        basis.x=NULL, basis.b = NULL, pf,tf){
return(fdata2model.penalty(vfunc=vfunc,vnf=vnf,response=response,data=data,
							basis.x=basis.x,basis.b=basis.b,pf=pf,tf=tf))
	}

fdata2model.penalty <- function(vfunc, vnf, response, data,
                        basis.x=NULL, basis.b = NULL, pf,tf
                        , lambda=NULL, P=NULL){
  
  kterms <- 1
  vs.list =name.coef=nam=beta.l=list()
  mean.list = basis.list = list()
  if (is.null(basis.x)) {basis.x=vector("list",length(vfunc));names(basis.x)=vfunc}
  if (length(vnf) > 0) {
    XX=data[["df"]][,c(response,vnf)] #data.frame el 1er elemento de la lista
    for ( i in 1:length(vnf)){
      #     print(paste("Non functional covariate:",vnf[i]))
      if (kterms > 1)   pf <- paste(pf, "+", vnf[i], sep = "")
      else pf <- paste(pf, vnf[i], sep = "")
      kterms <- kterms + 1
    }
    if   (attr(tf,"intercept")==0) {
      pf<- paste(pf,-1,sep="")
    }
  } else {
    XX=data$df[,response,drop=F]
    names(XX)=response
  }
  lpenalty <- ipenalty <- list()
  if (is.null(lambda)) {
    lambda0 <- FALSE}  else {
    lambda0 <- TRUE    }
  if (is.null(P)) {
    lambda0 <- FALSE}  else {
      lambda0 <- TRUE    }
  
  lenfunc <- length(vfunc)>0

  ipenal <- NCOL(XX)
  if (lenfunc) {
    for (i in 1:length(vfunc)) {
      # print(2)
	  fdat<-data[[vfunc[i]]]

      if (inherits(fdat, "fdata")) {
	  	  if (is.null(basis.x[[vfunc[i]]])) {
          if (is.null(basis.b[[vfunc[i]]])) { 
		  basis.b[[vfunc[i]]] <- basis.x[[vfunc[i]]] <- create.fdata.basis(fdat, l = 1:7)
          } else  {
		  basis.x[[vfunc[i]]] <-  basis.b[[vfunc[i]]]
		  }
        } else {
		  if (is.null(basis.b[[vfunc[i]]])) basis.b[[vfunc[i]]]<-basis.x[[vfunc[i]]] 
		}
		
	    nms <- fdat$names
		xaux<-fdata2basis(fdat,basis.x[[vfunc[i]]])

    Z <- xaux$coefs
		if ((basis.x[[vfunc[i]]]$type=="pc" | basis.x[[vfunc[i]]]$type=="pls") & identical(basis.b[[vfunc[i]]],basis.x[[vfunc[i]]])){
		J=diag(ncol(Z));colnames(J)=colnames(xaux$coefs)} else { 
		J <- inprodbasis(basis.x[[vfunc[i]]],basis.b[[vfunc[i]]])
		}
		Z <- Z %*% J
		colnames(Z)=paste0(vfunc[i],".",colnames(J))
    name.coef[[vfunc[i]]] <- paste(vfunc[i],".",colnames(J),sep="")
        lencoef <- length(colnames(Z))
        XX = cbind(XX, Z)
        for (j in 1:lencoef) {
          pf <- paste(pf, "+", name.coef[[vfunc[i]]][j], sep = "")
          kterms <- kterms + 1
        }       
        basis.list[[vfunc[i]]] <- xaux$basis
        vs.list[[vfunc[i]]] = J
		    mean.list[[vfunc[i]]]=xaux$mean

        if (lambda0) {
        #lpenalty[[vfunc[i]]] <- createMatrixPenalty(tt,lambda[[vfunc[i]]],P[[vfunc[i]]],vs=NULL)
        lpenalty[[vfunc[i]]] <- createMatrixPenalty(1:lencoef,lambda[[vfunc[i]]],P[[vfunc[i]]],vs=NULL)
        ipenalty[[vfunc[i]]] <- (ipenal+1):(ipenal+lencoef)
        # print(lpenalty[[vfunc[i]]])
        }
      } 
      else {
	    if (inherits(fdat,"fd")){
          if (is.null(basis.x[[vfunc[i]]]))  basis.x[[vfunc[i]]] <- fdat$basis
		      if (is.null(basis.b[[vfunc[i]]]))  basis.b[[vfunc[i]]] <- basis.x[[vfunc[i]]] 

          if (inherits(basis.x[[vfunc[i]]],"basisfd")) {
            r=fdat[["basis"]][["rangeval"]]
            if (!is.null( basis.x[[vfunc[i]]]$dropind)) {
              int<-setdiff(1:basis.x[[vfunc[i]]]$nbasis,basis.x[[vfunc[i]]]$dropind)
              basis.x[[vfunc[i]]]$nbasis<-length(int)
              basis.x[[vfunc[i]]]$dropind<-NULL
              basis.x[[vfunc[i]]]$names<-basis.x[[vfunc[i]]]$names[int]
            }

            J = inprodbasis(basis.x[[vfunc[i]]],basis.b[[vfunc[i]]]) 
            mean.list[[vfunc[i]]]<-mean.fd(fdat)
            x.fd<-center.fd(fdat)
            Z =t(fdat$coefs) %*% J
            colnames(J)=colnames(Z) = name.coef[[vfunc[i]]]=paste(vfunc[i],".",basis.x[[vfunc[i]]]$names,sep="")
            XX = cbind(XX,Z)
            for ( j in 1:length(colnames(Z))){
              if (kterms >= 1)  pf <- paste(pf, "+", colnames(Z)[j], sep = "")
              else pf <- paste(pf, colnames(Z)[j], sep = "")
              kterms <- kterms + 1
            }
            vs.list[[vfunc[i]]]<- J
          
          } else {  # basis.x is a pca.fd object
            l<-ncol(basis.x[[vfunc[i]]]$scores)
            vs <- diag(l) # Now matrix J is diagonal because basis.b is ignored.

            Z<-basis.x[[vfunc[i]]]$scores
            colnames(Z) <- colnames(vs)<-name.coef[[vfunc[i]]]=paste(vfunc[i], ".",colnames(basis.x[[vfunc[i]]]$harmonics$coefs),sep ="")
            XX = cbind(XX,Z)
            vs.list[[vfunc[i]]] = vs
            mean.list[[vfunc[i]]] = basis.x[[vfunc[i]]]$meanfd
            for ( j in 1:length(colnames(Z))){
              if (kterms >= 1)  pf <- paste(pf, "+", name.coef[[vfunc[i]]][j], sep = "")
              else pf <- paste(pf, name.coef[[vfunc[i]]][j], sep = "")
              kterms <- kterms + 1
            }
          }
        }
        else stop(paste(vfunc[i],"seems not to be a functional covariate"))
      }
    }  }
  else pf <- tf
  pf <- as.formula(pf)
  if (!is.data.frame(XX)) XX=data.frame(XX)
  return(list(pf=pf, mean.list=mean.list, basis.list = basis.list,XX=XX,
              basis.x=basis.x,basis.b=basis.b, vs.list=vs.list,name.coef=name.coef,
			  lpenalty=lpenalty, ipenalty=ipenalty, penalty=lambda0))
}
#####################################################

#####################################################
# # createMatrixPenalty <- function(tt,lambda,P,vs=NULL){
createMatrixPenalty <- function(tt,lambda,P,vs=NULL){
  np <- length(tt)
  # x <- Z
  for (i in 1:length(P))   {      if (P[i]!=0)         order.deriv<-i}
  P <- P.penalty(tt,P)
  if (!is.null(vs)){
     # Se escala para cuando son PCs
     P<-t(vs) %*% P %*% vs  # see fregre.pc
     rtt <- range(tt)
     P<-P*(diff(rtt) / (np -1))^(order.deriv*2-1)
     }
  #P <- P.penalty(1:3,c(1))
  return(lambda*P)
}


# #####################################################
# createMatrixPenalty(1:5,2,c(0,0,1))
#####################################################

  # incluir escalado para cuando son PC's
  # lpenalty[[vfunc[i]]] <- createMatrixPenalty(tt,lambda[[vfunc[i]]],P[[vfunc[i]]],vs=NULL){

