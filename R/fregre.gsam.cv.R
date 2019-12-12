##############################################
# Internal function (without predict)
##############################################
fregre.gsam.cv<-function(data, y, x, family = gaussian()
                         , alpha = 0.05, type.basis = "pc"
                         , ncomp = 8, ncomp.fix = TRUE
                         , par.basis=list(), kbs, par.model){
  #0-
#  print("eeentra CV")
  if (missing(y)) {stop("Please, specify the name of the response variable")  }
  #if (missing(alpha)) alpha<-.05
  if (missing(par.model)) par.model<-list() 
  if (missing(kbs)) kbs<--1
  fact <- 1.001
  #criterio <- "sp"
# print("entra gsam CV")
  namdata<-names(data)
  idf<-which("df"==namdata )
  ydat<-data$df[,y]
  namfunc<-names(data[-idf])
  namnfunc<-names(data$df)
  infunc<-setdiff(namnfunc,y)
  basis1<-NULL
  if (missing(x)) {  
    #    ldata<-as.list(data[idf][,y,drop=F])
    #xdatos<-data[-idf][1]    
     #x<-names(data[-idf])[1]
    ifunc<-namfunc
    #infunc<-intersect(c(y,x),namnfunc)
    xdatos<-as.list(data$df[,,drop=F])
    xdatos<-c(xdatos,data[namfunc])
    x<-c(namnfunc,namfunc)
   # warning("La primera variables funcional es considerada")
  }   else {
    ifunc<-intersect(x,namfunc)
    infunc<-intersect(c(x),namnfunc)
    xdatos<-as.list(data$df[,c(y,infunc),drop=F])
    xdatos<-c(xdatos,data[ifunc])
  } 
  ldata0<-xdatos
  resp<-y
  xentran<-NULL
  tipoentran<-NULL
  parar<-FALSE
  it<-1
  xmodelo<-basisx<-NULL
  npredictors<-length(ldata0)-1
  ipredictors<-numeric(npredictors)
  names(ipredictors)<- setdiff(names(ldata0),resp)
  xentra <- x
  tbasis <-  type.basis
  iposibles<-xposibles<-NULL
  xdata<-NULL
  xposibles <- c(infunc)
  xdata<-cbind(data$df[,c(y,infunc),drop=F])
  lenfunc <- ifelse(length(ifunc)==0,FALSE,TRUE)
  if (lenfunc){
    for (i in 1:length(ifunc)){
    switch(tbasis,"pc"={
          basis1<-create.pc.basis(ldata0[[ifunc[i]]],1:ncomp)
          icomp <- colnames(basis1$x)
        },
        "pls"={
          basis1<-create.pls.basis(ldata0[[ifunc[i]]],ldata0[[y]],1:ncomp)
          icomp <- colnames(basis1$x)
        },"basis"={
          ncomp.fix=TRUE
          basis1<-create.fdata.basis(ldata0[[ifunc[i]]],
                                     1:ncomp)
          xx <- fdata.cen(ldata0[[ifunc[i]]])
          xmean = xx[[2]]
          xcen = xx[[1]]
          x.fd = Data2fd(argvals = ldata0[[ifunc[i]]]$argvals,
                         y = t(xcen$data), 
                         basisobj = basis1)
          
          basis1$x = t(x.fd$coefs)          
        #  print(dim(basis1$x ))
          icomp <-  basis1$names
          print(icomp)
      }
    )
#print(paste(ifunc[i],".",colnames(basis1$x),sep="") )
    colnames(basis1$x)<-paste(ifunc[i],".",icomp,sep="") 
    xposibles <- c(xposibles,colnames(basis1$x)[1:ncomp])
      xdata<-cbind(xdata,basis1$x[,1:ncomp,drop=F])
    }
  }
  iposibles<-1:length(xposibles)
  maxvar <- length(xposibles)
  iposibles <- 1:maxvar
  form.nl <- paste(y,"~")
  nam.model<-"gam"    
  par.model<-list()
  ibest <- 1
  k <- 1
  form.ini <-  paste(y,"~ 1")
  # print("best.model")
  best.model <- gam(as.formula(form.ini),data=xdata,family=family)#,weights=weights)
  #print("best.model");  print(best.model)
  nvar <- 1
  res<-NULL
  best.aic<-Inf
  parar <- FALSE
  imodelo<-xmodelo<-NULL
  while (!parar) {
    # print("while");  print(nvar)
    imax <- 1:(maxvar-nvar+1)
    form.nl <- form.ini
    lmax<-length(imax)
    aic<-numeric(lmax)
    for (k in 1:lmax) {  
    #   print(k)
    #  print(xposibles[k])
    if (ncomp.fix) ik<-1:k
    else ik<-k
       if (is.factor(xdata[[xposibles[k]]]))
         fpredictors.nl <- paste("~.+",xposibles[k],sep="",collapse="+")
       else {
         fpredictors.nl <- paste("s(",xposibles[ik],",k=",kbs,")",sep="",collapse="+")
         fpredictors.nl <- paste0("~.+", fpredictors.nl)
         }
       form.nl <-update.formula(as.formula(form.ini), as.formula(fpredictors.nl))
       #as.formula(form.ini,fpredictors.nl)
       par.model$formula<-form.nl #as.formula(form.nl)
       par.model$data<-xdata
       res[[k]]<-do.call(nam.model,par.model)             
       
       aic[k] <- res[[k]]$gcv.ubre
      } 
  posible.aic <- which.min(aic)
  pmejora <- alpha
  if ((aic[posible.aic]*fact)< best.aic) {
    #print("entra mejora");   print(aic[posible.aic]);   print(best.aic)
    #print(best.model$formula)
    #print(res[[posible.aic]]$formula)
    pmejora <- pvalue.anova(best.model,res[[posible.aic]])
    pmejora<- pmejora[["Pr(>F)"]][2]#na.omit(pmejora[["Pr(>F)"]])
  #  print(pmejora);    print(alpha)  
  if (pmejora < alpha & !is.na(pmejora)) {
        #if (criterio=="p.signif") best.aic<-alpha # best.aic deberia llamarse best.criteria
        #else best.aic<-aic[posible.aic]  
      best.aic <- aic[posible.aic]
      best.model <- res[[posible.aic]]
      if (ncomp.fix){
        nvar <- posible.aic
        xmodelo <- c(xmodelo,xposibles[1:posible.aic]) 
        imodelo <- c(imodelo,iposibles[1:posible.aic]) 
        xposibles <- xposibles[-c(1:posible.aic)]
        iposibles <- iposibles[-c(1:posible.aic)]
      } else {
        xmodelo <- c(xmodelo,xposibles[posible.aic]) 
        imodelo <-c(imodelo,iposibles[posible.aic]) 
        xposibles <- xposibles[-posible.aic]
        iposibles <- iposibles[-posible.aic]
      }
      form.ini <-   res[[posible.aic]]$formula
      }    else {      parar=TRUE    }
    #cat("nvar");print(nvar)
    if (nvar==ncomp) {      parar=TRUE    }    else { nvar<- nvar+1}
  } else {     parar=TRUE    }
} #fin while
 #para PC
 vfunc<-x
 i <- 1
 z <- best.model
 z$basis.x=z$basis.b=vs.list=mean.list=list()
 #basis1$basis<- basis1$basis#[imodelo,,drop=F]
 vs.list[[vfunc[i]]] = basis1$basis
 mean.list[[vfunc[i]]] = basis1$mean
 JJ = NULL
 z <- best.model
 if (lenfunc) {
   z$mean = mean.list
   z$basis.x[[vfunc[i]]] = basis1
   z$basis.b[[vfunc[i]]] = basis1
   z$JJ <- JJ
   z$vs.list = vs.list
   z$vfunc <- vfunc
   colnames(z$basis.x[[vfunc[i]]]$x)<-paste("PC",1:ncol(z$basis.x[[vfunc[i]]]$x),sep="")
 }
 #print("finnn")
 if (is.null(imodelo)) {
   warning("No component of the variable is significant") 
 }
 #else{
 
  z$ncomp.opt<-imodelo
  z$formula <- as.formula(z$formula) 
  z$formula.ini =    as.formula( paste(y,"~+s(",x,",k=",kbs,")",sep=""))
# print(summary(z)); print(z$formula) ; print(z$formula.ini) 
  z$data = data
  z$XX = best.model$model
#print(colnames(z$XX)) ;print("optimas");print(z$pc.opt)
# z$nnf <- nnf #var no  funcionales
 class(z) <- c(class(z), "fregre.gsam")
 #print("sale fregre.gsam.cv") 
 return(z)
}
###################################################
#sp<-fda.usc:::sp
#pvalue.anova<-fda.usc:::pvalue.anova
 # res.gsam.cv<-fregre.gsam.cv(la,resp,alpha=0.05,
 #                             type.basis="basis",
 #                             criterio="sp",ncomp=5)
 # summary(res.gsam.cv)
 # 
 # pr<-predict.fregre.gsam(res.gsam.cv,la)
 # 
 # res.gsam.cv<-fregre.gsam.cv(la,resp,alpha=0.05,
 #                             type.basis="pc",
 #                             criterio="sp",ncomp=8)
 # summary(res.gsam.cv)
 #  class(res.gsam.cv) 
 #  names(res.gsam.cv)
 #  res.gsam.cv$
 #  pr<-predict.fregre.gsam(res.gsam.cv,la)
 # 
# res.gsam.cv<-fregre.gsam.cv(la,resp,
#                             c("x1","x2","x"),ncomp=2)
# summary(res.gsam.cv)
# 
# 
# res.gsam.cv<-fregre.gsam.cv(la,resp,c("xcat0","xcat2" ,"xcat1"),alpha=0.05,
#                             type.basis="pc",
#                             criterio="sp",ncomp=2)
# summary(res.gsam.cv)
# 
# 
# 
# resp<-"y"
# res.gsam.cv<-fregre.gsam.cv(la,resp,c("x1","x2","x"),alpha=0.05,
#                             type.basis="pc",
#                             criterio="sp",ncomp=8)
# summary(res.gsam.cv)
# 
# 
# res.gsam.cv<-fregre.gsam.cv(la,resp,alpha=0.05,
#                             type.basis="pc",
#                             criterio="sp",ncomp=8)
# summary(res.gsam.cv)
# 
