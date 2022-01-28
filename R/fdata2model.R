#####################################################
fdata2model <- function(vfunc, vnf, response, data, basis.x=NULL,pf,tf){
  kterms = 1
  vs.list=mean.list=name.coef=nam=beta.l=list()
  bsp1=TRUE
  if (length(vnf)>0) {
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
    XX=data$df[response]
    names(XX)=response
  }
  lenfunc<-length(vfunc)>0
  if (lenfunc) {
    for (i in 1:length(vfunc)) {
      if (class(data[[vfunc[i]]])[1]=="fdata"){
        tt<-data[[vfunc[i]]][["argvals"]]
        rtt<-data[[vfunc[i]]][["rangeval"]]
        fdat<-data[[vfunc[i]]];      dat<-data[[vfunc[i]]]$data
        if (is.null(basis.x[[vfunc[i]]]))  basis.x[[vfunc[i]]]<-create.fdata.basis(fdat,l=1:7)
        else   if (basis.x[[vfunc[i]]]$type=="pc" | basis.x[[vfunc[i]]]$type=="pls") bsp1=FALSE
        if (bsp1) {
          if (is.null(rownames(dat)))    rownames(fdat$data)<-1:nrow(dat)
          fdnames=list("time"=tt,"reps"=rownames(fdat[["data"]]),"values"="values")
          xcc<-fdata.cen(data[[vfunc[i]]])
          mean.list[[vfunc[i]]]=xcc[[2]]
          if (!is.null( basis.x[[vfunc[i]]]$dropind)) {
            int<-setdiff(1:basis.x[[vfunc[i]]]$nbasis,basis.x[[vfunc[i]]]$dropind)
            basis.x[[vfunc[i]]]$nbasis<-length(int)
            basis.x[[vfunc[i]]]$dropind<-NULL
            basis.x[[vfunc[i]]]$names<-basis.x[[vfunc[i]]]$names[int]
          }
          x.fd = Data2fd(argvals = tt, y = t(xcc[[1]]$data),basisobj = basis.x[[vfunc[i]]],fdnames=fdnames)
          r=x.fd[[2]][[3]]
          J=inprod(basis.x[[vfunc[i]]],basis.x[[vfunc[i]]])
          J12=inprod(basis.x[[vfunc[i]]])
# print(J[1:2,1:3]); print(J12[1:2,1:3])          
          Z =t(x.fd$coefs) %*% J
          colnames(J)=colnames(Z) = name.coef[[vfunc[i]]]=paste(vfunc[i],".",basis.x[[vfunc[i]]]$names,sep="")
          XX = cbind(XX,Z)
          for ( j in 1:length(colnames(Z))){
            if (kterms >= 1)  pf <- paste(pf, "+", colnames(Z)[j], sep = "")
            else pf <- paste(pf, colnames(Z)[j], sep = "")
            kterms <- kterms + 1
          }
          vs.list[[vfunc[i]]]<-J
        }        else {    #PC o PLS
          l<-basis.x[[vfunc[i]]]$l
          lenl<-length(l)
          vs <- t(basis.x[[vfunc[i]]]$basis$data)
          Z<-basis.x[[vfunc[i]]]$x[,l,drop=FALSE]     
          response = "y"
          colnames(Z) = name.coef[[vfunc[i]]]=paste(vfunc[i], ".",colnames(Z),sep ="")      
          name.coef[[vfunc[i]]]<-colnames(Z)
          XX = cbind(XX,Z)
          vs.list[[vfunc[i]]]=basis.x[[vfunc[i]]]$basis
          mean.list[[vfunc[i]]]=basis.x[[vfunc[i]]]$mean
          for ( j in 1:length(colnames(Z))){
            if (kterms >= 1)  pf <- paste(pf, "+", name.coef[[vfunc[i]]][j], sep = "")
            else pf <- paste(pf, name.coef[[vfunc[i]]][j], sep = "")
            kterms <- kterms + 1
          }       
        }
      }
      # else {
      #   if(class(data[[vfunc[i]]])[1]=="fd"){
      #     fdat<-data[[vfunc[i]]]
      #     if (is.null(basis.x[[vfunc[i]]]))  basis.x[[vfunc[i]]]<-fdat$basis
      #     else   if (class(basis.x[[vfunc[i]]])=="pca.fd") bsp1=FALSE
      #     if (bsp1) {
      #       r=fdat[[2]][[3]]
      #       if (!is.null( basis.x[[vfunc[i]]]$dropind)) {
      #         int<-setdiff(1:basis.x[[vfunc[i]]]$nbasis,basis.x[[vfunc[i]]]$dropind)
      #         basis.x[[vfunc[i]]]$nbasis<-length(int)
      #         basis.x[[vfunc[i]]]$dropind<-NULL
      #         basis.x[[vfunc[i]]]$names<-basis.x[[vfunc[i]]]$names[int]
      #       }
      #       J=inprod(basis.x[[vfunc[i]]],basis.x[[vfunc[i]]])
      #       mean.list[[vfunc[i]]]<-mean.fd(x.fd)
      #       x.fd<-center.fd(x.fd)
      #       Z =t(x.fd$coefs) %*% J
      #       colnames(J)=colnames(Z) = name.coef[[vfunc[i]]]=paste(vfunc[i],".",basis.x[[vfunc[i]]]$names,sep="")
      #       XX = cbind(XX,Z)
      #       for ( j in 1:length(colnames(Z))){
      #         if (kterms >= 1)  pf <- paste(pf, "+", colnames(Z)[j], sep = "")
      #         else pf <- paste(pf, colnames(Z)[j], sep = "")
      #         kterms <- kterms + 1
      #       }
      #       vs.list[[vfunc[i]]]<-J
      #     }          else {
      #       l<-ncol(basis.x[[vfunc[i]]]$scores)
      #       vs <- basis.x[[vfunc[i]]]$harmonics$coefs
      #       Z<-basis.x[[vfunc[i]]]$scores
      #       response = "y"
      #       colnames(Z) = name.coef[[vfunc[i]]]=paste(vfunc[i], ".",colnames(basis.x[[vfunc[i]]]$harmonics$coefs),sep ="")
      #       XX = cbind(XX,Z)
      #       vs.list[[vfunc[i]]]=vs
      #       mean.list[[vfunc[i]]]=basis.x[[vfunc[i]]]$meanfd
      #       for ( j in 1:length(colnames(Z))){
      #         if (kterms >= 1)  pf <- paste(pf, "+", name.coef[[vfunc[i]]][j], sep = "")
      #         else pf <- paste(pf, name.coef[[vfunc[i]]][j], sep = "")
      #         kterms <- kterms + 1
      #       }
      #     }
      #   }
      #   else stop("Please, enter functional covariate")
      # }
    }  }
   else pf <- tf
  pf<-as.formula(pf)
  if (!is.data.frame(XX)) XX=data.frame(XX)
  return(list(pf=pf,vs.list=vs.list,mean.list=mean.list,XX=XX,basis.x=basis.x))
}
#####################################################
