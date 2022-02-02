# Internal function used in fregre.glm and fregre.gsam


# fdata2model.penalty <- function(vfunc, vnf, response, data,
#                                 basis.x=NULL, basis.b = NULL, pf,tf
#                                 , lambda=NULL, P=NULL){
#   # print("entra fdata2model.penalty")
#   kterms <- 1
#   vs.list = mean.list=name.coef=nam=beta.l=list()
#   bsp1 <- TRUE
#   if (length(vnf) > 0) {
#     XX=data[["df"]][,c(response,vnf)] #data.frame el 1er elemento de la lista
#     for ( i in 1:length(vnf)){
#       #     print(paste("Non functional covariate:",vnf[i]))
#       if (kterms > 1)   pf <- paste(pf, "+", vnf[i], sep = "")
#       else pf <- paste(pf, vnf[i], sep = "")
#       kterms <- kterms + 1
#     }
#     if   (attr(tf,"intercept")==0) {
#       pf<- paste(pf,-1,sep="")
#     }
#   } else {
#     XX=data$df[,response,drop=F]
#     names(XX)=response
#   }
#   lpenalty <- ipenalty <- list()
#   if (is.null(lambda)) {
#     lambda0 <- FALSE}  else {
#       lambda0 <- TRUE    }
#   if (is.null(P)) {
#     lambda0 <- FALSE}  else {
#       lambda0 <- TRUE    }
#   
#   lenfunc <- length(vfunc)>0
#   mean.list = basis.list = list()
#   ipenal <- NCOL(XX)
#   if (lenfunc) {
#     for (i in 1:length(vfunc)) {
#       # print(2)
#       if (is(data[[vfunc[i]]], "fdata")) {
#         tt <- data[[vfunc[i]]][["argvals"]]
#         rtt <- data[[vfunc[i]]][["rangeval"]]
#         fdat <- data[[vfunc[i]]]
#         nms <- data[[vfunc[i]]]$names
#         #dat <- data[[vfunc[i]]]$data
#         if (is.null(basis.x[[vfunc[i]]])) {
#           if (is.null(basis.b[[vfunc[i]]])) { basis.b[[vfunc[i]]] <- basis.x[[vfunc[i]]] <- 
#             create.fdata.basis(fdat, l = 1:7)
#           } else  {basis.x[[vfunc[i]]] <-  basis.b[[vfunc[i]]] }
#         } else
#           if (basis.x[[vfunc[i]]]$type == "pc" | basis.x[[vfunc[i]]]$type == "pls") 
#             bsp1 = FALSE
#         #         if (bsp1) {
#         xaux <- fdata2basis(data[[vfunc[i]]],basis.x[[vfunc[i]]])
#         name.coef[[vfunc[i]]] <- colnames(xaux$coefs) <- paste(vfunc[i],".",colnames(xaux$coefs),sep="")
#         Z <- xaux$coefs
#         if (!is.null(basis.b)){
#           J = inprod(basis.x[[vfunc[i]]], basis.b[[vfunc[i]]])
#           #mean.list[[vfunc[i]]] <- mean.fd(x.fd);          x.fd <- center.fd(x.fd)
#           colnam <- colnames(Z)
#           Z <- Z %*% J
#           name.coef[[vfunc[i]]] <- colnames(Z) <- colnam[1:NCOL(Z)]
#         }
#         lencoef <- length(colnames(Z))
#         XX = cbind(XX, Z)
#         for (j in 1:lencoef) {
#           pf <- paste(pf, "+", name.coef[[vfunc[i]]][j], sep = "")
#           kterms <- kterms + 1
#         }       
#         basis.list[[vfunc[i]]] <- xaux$basis
#         # J=inprod(basis.x[[vfunc[i]]],basis.b[[vfunc[i]]])
#         #   vs.list[[vfunc[i]]] = basis.x[[vfunc[i]]]$basis
#         if (!bsp1) 
#           mean.list[[vfunc[i]]] = basis.x[[vfunc[i]]]$mean
#         else {
#           xcc <- fdata.cen(data[[vfunc[i]]])
#           mean.list[[vfunc[i]]] = xcc[[2]]
#         }
#         if (lambda0) {
#           #lpenalty[[vfunc[i]]] <- createMatrixPenalty(tt,lambda[[vfunc[i]]],P[[vfunc[i]]],vs=NULL)
#           lpenalty[[vfunc[i]]] <- createMatrixPenalty(1:lencoef,lambda[[vfunc[i]]],P[[vfunc[i]]],vs=NULL)
#           ipenalty[[vfunc[i]]] <- (ipenal+1):(ipenal+lencoef)
#           # print(lpenalty[[vfunc[i]]])
#         }
#       } 
#       else {
#         # print("fd")
#         if(class(data[[vfunc[i]]])[1]=="fd"){
#           fdat<-data[[vfunc[i]]]
#           if (is.null(basis.x[[vfunc[i]]]))  basis.x[[vfunc[i]]] <- fdat$basis
#           else   if (class(basis.x[[vfunc[i]]])=="pca.fd") bsp1=FALSE
#           if (bsp1) {
#             r=fdat[[2]][[3]]
#             if (!is.null( basis.x[[vfunc[i]]]$dropind)) {
#               int<-setdiff(1:basis.x[[vfunc[i]]]$nbasis,basis.x[[vfunc[i]]]$dropind)
#               basis.x[[vfunc[i]]]$nbasis<-length(int)
#               basis.x[[vfunc[i]]]$dropind<-NULL
#               basis.x[[vfunc[i]]]$names<-basis.x[[vfunc[i]]]$names[int]
#             }
#             if (is.null(basis.b[[vfunc[i]]]))
#               basis.b[[vfunc[i]]] <- basis.x[[vfunc[i]]]
#             J = inprod(basis.x[[vfunc[i]]],basis.b[[vfunc[i]]])
#             mean.list[[vfunc[i]]]<-mean.fd(x.fd)
#             x.fd<-center.fd(x.fd)
#             Z =t(x.fd$coefs) %*% J
#             colnames(J)=colnames(Z) = name.coef[[vfunc[i]]]=paste(vfunc[i],".",basis.x[[vfunc[i]]]$names,sep="")
#             XX = cbind(XX,Z)
#             for ( j in 1:length(colnames(Z))){
#               if (kterms >= 1)  pf <- paste(pf, "+", colnames(Z)[j], sep = "")
#               else pf <- paste(pf, colnames(Z)[j], sep = "")
#               kterms <- kterms + 1
#             }
#             vs.list[[vfunc[i]]]<- J
#           }          else {
#             l<-ncol(basis.x[[vfunc[i]]]$scores)
#             vs <- basis.x[[vfunc[i]]]$harmonics$coefs
#             Z<-basis.x[[vfunc[i]]]$scores
#             #response = "y"
#             colnames(Z) = name.coef[[vfunc[i]]]=paste(vfunc[i], ".",colnames(basis.x[[vfunc[i]]]$harmonics$coefs),sep ="")
#             XX = cbind(XX,Z)
#             vs.list[[vfunc[i]]] = vs
#             mean.list[[vfunc[i]]] = basis.x[[vfunc[i]]]$meanfd
#             for ( j in 1:length(colnames(Z))){
#               if (kterms >= 1)  pf <- paste(pf, "+", name.coef[[vfunc[i]]][j], sep = "")
#               else pf <- paste(pf, name.coef[[vfunc[i]]][j], sep = "")
#               kterms <- kterms + 1
#             }
#           }
#         }
#         else stop("Please, enter functional covariate")
#       }
#     }  }
#   else pf <- tf
#   pf <- as.formula(pf)
#   if (!is.data.frame(XX)) XX=data.frame(XX)
#   #  print("sale fdata2model")  
#   colnames(XX)[1] <-  response
#   # print("sale fdata2model.penalty")
#   return(list(pf=pf, vs.list=vs.list, mean.list=mean.list,
#               basis.list = basis.list,XX=XX,
#               basis.x=basis.x,basis.b=basis.b, name.coef=name.coef, bsp1=bsp1
#               ,lpenalty=lpenalty, ipenalty=ipenalty, penalty=lambda0))
# }

#####################################################
# Internal function used in .....
#
# If basis.b is null, returns the design (or model) matrix, 
# otherwise the functional data coefficients.
# 
fdata2model <- function(vfunc, vnf, response, data,
                        basis.x=NULL, basis.b = NULL, pf,tf){
# print("fdata2model")
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
# print(1)  
  lenfunc<-length(vfunc)>0
  if (lenfunc) {
    for (i in 1:length(vfunc)) {
      # print(2)
      if (class(data[[vfunc[i]]])[1]=="fdata"){
        tt<-data[[vfunc[i]]][["argvals"]]
        rtt<-data[[vfunc[i]]][["rangeval"]]
        fdat<-data[[vfunc[i]]];      dat<-data[[vfunc[i]]]$data
        if (is.null(basis.x[[vfunc[i]]]))  basis.x[[vfunc[i]]]<-create.fdata.basis(fdat,l=1:7)
        else   if (basis.x[[vfunc[i]]]$type=="pc" | basis.x[[vfunc[i]]]$type=="pls") bsp1=FALSE
        if (bsp1) {
          # print(bsp1)
          # print(23)
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
#print(3)
          x.fd = Data2fd(argvals = tt, y = t(xcc[[1]]$data),basisobj = basis.x[[vfunc[i]]],fdnames=fdnames)
          r=x.fd[[2]][[3]]
          if (is.null(basis.b[[vfunc[i]]]))
            basis.b[[vfunc[i]]] <- basis.x[[vfunc[i]]]
          J = inprod(basis.x[[vfunc[i]]],basis.b[[vfunc[i]]]) 
  #        print(4)          
          Z = t(x.fd$coefs) %*% J
   #       print(5)          
colnames(J) = colnames(Z) = name.coef[[vfunc[i]]] = paste(vfunc[i],".",basis.b[[vfunc[i]]]$names,sep="")
# print(6)          
          XX = cbind(XX,Z)
    #      print(6)          
          for ( j in 1:length(colnames(Z))){
            if (kterms >= 1)  pf <- paste(pf, "+", colnames(Z)[j], sep = "")
            else pf <- paste(pf, colnames(Z)[j], sep = "")
            kterms <- kterms + 1
          }
          vs.list[[vfunc[i]]] <- J
        }        else {    #PC o PLS
#          print(PC)
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
      else {
     #   print("fd")
        if(class(data[[vfunc[i]]])[1]=="fd"){
          fdat<-data[[vfunc[i]]]
          if (is.null(basis.x[[vfunc[i]]]))  basis.x[[vfunc[i]]]<-fdat$basis
          else   if (class(basis.x[[vfunc[i]]])=="pca.fd") bsp1=FALSE
          if (bsp1) {
            r=fdat[[2]][[3]]
            if (!is.null( basis.x[[vfunc[i]]]$dropind)) {
              int<-setdiff(1:basis.x[[vfunc[i]]]$nbasis,basis.x[[vfunc[i]]]$dropind)
              basis.x[[vfunc[i]]]$nbasis<-length(int)
              basis.x[[vfunc[i]]]$dropind<-NULL
              basis.x[[vfunc[i]]]$names<-basis.x[[vfunc[i]]]$names[int]
            }
            if (is.null(basis.b[[vfunc[i]]]))
              basis.b[[vfunc[i]]] <- basis.x[[vfunc[i]]]
            J = inprod(basis.x[[vfunc[i]]],basis.b[[vfunc[i]]])
            mean.list[[vfunc[i]]]<-mean.fd(x.fd)
            x.fd<-center.fd(x.fd)
            Z =t(x.fd$coefs) %*% J
            colnames(J)=colnames(Z) = name.coef[[vfunc[i]]]=paste(vfunc[i],".",basis.x[[vfunc[i]]]$names,sep="")
            XX = cbind(XX,Z)
            for ( j in 1:length(colnames(Z))){
              if (kterms >= 1)  pf <- paste(pf, "+", colnames(Z)[j], sep = "")
              else pf <- paste(pf, colnames(Z)[j], sep = "")
              kterms <- kterms + 1
            }
            vs.list[[vfunc[i]]]<- J
          }          else {
            l<-ncol(basis.x[[vfunc[i]]]$scores)
            vs <- basis.x[[vfunc[i]]]$harmonics$coefs
            Z<-basis.x[[vfunc[i]]]$scores
            response = "y"
            colnames(Z) = name.coef[[vfunc[i]]]=paste(vfunc[i], ".",colnames(basis.x[[vfunc[i]]]$harmonics$coefs),sep ="")
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
        else stop("Please, enter functional covariate")
      }
    }  }
   else pf <- tf
# print(2)  ;#print(pf)
  pf<-as.formula(pf)
  
  if (!is.data.frame(XX)) XX=data.frame(XX)
#  print("sale fdata2model")  
  return(list(pf=pf,vs.list=vs.list,mean.list=mean.list,XX=XX,
              basis.x=basis.x,name.coef=name.coef,bsp1=bsp1))
}
#####################################################

#####################################################
# Internal function used in .....
#
# If basis.b is null, returns the design (or model) matrix, 
# otherwise the functional data coefficients.
# 
fdata2model.penalty <- function(vfunc, vnf, response, data,
                        basis.x=NULL, basis.b = NULL, pf,tf
                        , lambda=NULL, P=NULL){
  #print("entra fdata2model.penalty")
  #out <- fdata2model.penalty(vfunc, vnf, response, data, basis.x,basis.b,pf,tf,lambda,P)
  
  # Para borrar
  #######################
  # print(lambda)
  # print(P)
  kterms <- 1
  vs.list = mean.list=name.coef=nam=beta.l=list()
  bsp1 <- TRUE
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
  mean.list = basis.list = list()
  ipenal <- NCOL(XX)
  if (lenfunc) {
    for (i in 1:length(vfunc)) {
      # print(2)
      if (is(data[[vfunc[i]]], "fdata")) {
        tt <- data[[vfunc[i]]][["argvals"]]
        rtt <- data[[vfunc[i]]][["rangeval"]]
        fdat <- data[[vfunc[i]]]
        nms <- data[[vfunc[i]]]$names
        #dat <- data[[vfunc[i]]]$data
        if (is.null(basis.x[[vfunc[i]]])) {
          if (is.null(basis.b[[vfunc[i]]])) { basis.b[[vfunc[i]]] <- basis.x[[vfunc[i]]] <- 
          create.fdata.basis(fdat, l = 1:7)
          } else  {basis.x[[vfunc[i]]] <-  basis.b[[vfunc[i]]] }
        } else
            if (basis.x[[vfunc[i]]]$type == "pc" | basis.x[[vfunc[i]]]$type == "pls") 
              bsp1 = FALSE
        #         if (bsp1) {
        xaux <- fdata2basis(data[[vfunc[i]]],basis.x[[vfunc[i]]])
        name.coef[[vfunc[i]]] <- colnames(xaux$coefs) <- paste(vfunc[i],".",colnames(xaux$coefs),sep="")
        Z <- xaux$coefs
        if (!is.null(basis.b)){
          J = inprod(basis.x[[vfunc[i]]], basis.b[[vfunc[i]]])
          #mean.list[[vfunc[i]]] <- mean.fd(x.fd);          x.fd <- center.fd(x.fd)
          colnam <- colnames(Z)
          Z <- Z %*% J
          name.coef[[vfunc[i]]] <- colnames(Z) <- colnam[1:NCOL(Z)]
        }
        lencoef <- length(colnames(Z))
        XX = cbind(XX, Z)
        for (j in 1:lencoef) {
          pf <- paste(pf, "+", name.coef[[vfunc[i]]][j], sep = "")
          kterms <- kterms + 1
        }       
        basis.list[[vfunc[i]]] <- xaux$basis
        # J=inprod(basis.x[[vfunc[i]]],basis.b[[vfunc[i]]])
        #   vs.list[[vfunc[i]]] = basis.x[[vfunc[i]]]$basis
        if (!bsp1) 
          mean.list[[vfunc[i]]] = basis.x[[vfunc[i]]]$mean
        else {
          xcc <- fdata.cen(data[[vfunc[i]]])
          mean.list[[vfunc[i]]] = xcc[[2]]
        }
        if (lambda0) {
        #lpenalty[[vfunc[i]]] <- createMatrixPenalty(tt,lambda[[vfunc[i]]],P[[vfunc[i]]],vs=NULL)
        lpenalty[[vfunc[i]]] <- createMatrixPenalty(1:lencoef,lambda[[vfunc[i]]],P[[vfunc[i]]],vs=NULL)
        ipenalty[[vfunc[i]]] <- (ipenal+1):(ipenal+lencoef)
        # print(lpenalty[[vfunc[i]]])
        }
      } 
      else {
        # print("fd")
        if(class(data[[vfunc[i]]])[1]=="fd"){
          fdat<-data[[vfunc[i]]]
          if (is.null(basis.x[[vfunc[i]]]))  basis.x[[vfunc[i]]] <- fdat$basis
          else   if (class(basis.x[[vfunc[i]]])=="pca.fd") bsp1=FALSE
          if (bsp1) {
            r=fdat[[2]][[3]]
            if (!is.null( basis.x[[vfunc[i]]]$dropind)) {
              int<-setdiff(1:basis.x[[vfunc[i]]]$nbasis,basis.x[[vfunc[i]]]$dropind)
              basis.x[[vfunc[i]]]$nbasis<-length(int)
              basis.x[[vfunc[i]]]$dropind<-NULL
              basis.x[[vfunc[i]]]$names<-basis.x[[vfunc[i]]]$names[int]
            }
            if (is.null(basis.b[[vfunc[i]]]))
              basis.b[[vfunc[i]]] <- basis.x[[vfunc[i]]]
            J = inprod(basis.x[[vfunc[i]]],basis.b[[vfunc[i]]])
            mean.list[[vfunc[i]]]<-mean.fd(x.fd)
            x.fd<-center.fd(x.fd)
            Z =t(x.fd$coefs) %*% J
            colnames(J)=colnames(Z) = name.coef[[vfunc[i]]]=paste(vfunc[i],".",basis.x[[vfunc[i]]]$names,sep="")
            XX = cbind(XX,Z)
            for ( j in 1:length(colnames(Z))){
              if (kterms >= 1)  pf <- paste(pf, "+", colnames(Z)[j], sep = "")
              else pf <- paste(pf, colnames(Z)[j], sep = "")
              kterms <- kterms + 1
            }
            vs.list[[vfunc[i]]]<- J
          }          else {
            l<-ncol(basis.x[[vfunc[i]]]$scores)
            vs <- basis.x[[vfunc[i]]]$harmonics$coefs
            Z<-basis.x[[vfunc[i]]]$scores
#response = "y"
            colnames(Z) = name.coef[[vfunc[i]]]=paste(vfunc[i], ".",colnames(basis.x[[vfunc[i]]]$harmonics$coefs),sep ="")
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
        else stop("Please, enter functional covariate")
      }
    }  }
  else pf <- tf
  pf <- as.formula(pf)
  if (!is.data.frame(XX)) XX=data.frame(XX)
  #  print("sale fdata2model")  
  colnames(XX)[1] <-  response
  # print("sale fdata2model.penalty")
  return(list(pf=pf, vs.list=vs.list, mean.list=mean.list,
              basis.list = basis.list,XX=XX,
              basis.x=basis.x,basis.b=basis.b, name.coef=name.coef, bsp1=bsp1
              ,lpenalty=lpenalty, ipenalty=ipenalty, penalty=lambda0))
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