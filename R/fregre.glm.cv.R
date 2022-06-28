##############################################
# Internal function (without predict)
# Pensada para selecionar por AIC u otro criterios las bases/coef significativos 
# de una covariable funcional
##############################################
fregre.glm.cv<-function(data, y, x, family = gaussian()
                         , alpha = 0.05, type.basis = "pc"
                         , numbasis = 3, numbasis.opt = FALSE
                        , par.model){
  #0-
#print("eeentra glm.cv")
  
  if (missing(y)) {stop("Please, specify the name of the response variable")  }
  #if (missing(alpha)) alpha<-.05
  if (missing(par.model)) par.model<-list() 
  fact <- 1.001
  #criterio <- "sp"
# print("entra glm CV")
  namdata<-names(data)
  idf<-which("df"==namdata )
  ydat<-data$df[,y]
  namfunc<-names(data[-idf])
  namnfunc<-names(data$df)
  infunc<-setdiff(namnfunc,y)
  basis1<-NULL
  if (missing(x)) {  
    ifunc<-namfunc
    xdatos<-as.list(data$df[,,drop=F])
    xdatos<-c(xdatos,data[namfunc])
    x<-c(namnfunc,namfunc)
  }   else {
    ifunc<-intersect(x,namfunc)
    infunc<-intersect(c(x),namnfunc)
    xdatos<-as.list(data$df[,c(y,infunc),drop=F])
    xdatos<-c(xdatos,data[ifunc])
  } 
  #print(names(xdatos))
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
  iposibles<-xposibles<-NULL
  xdata<-NULL
  xposibles <- c(infunc)
  xdata<-cbind(data$df[,c(y,infunc),drop=F])
  lenfunc <- ifelse(length(ifunc)==0,FALSE,TRUE)
  tbasis <- type.basis
  if (type.basis %in% c("fourier","bspline")) 
    {
    tbasis="basis"
    k <- min(5,numbasis)
    if (numbasis<5) warning("Numbasis has been increased to 5")
  }
  if (lenfunc){
    for (i in 1:length(ifunc)){
      switch(tbasis,"pc"={
         
          basis1<-create.pc.basis(ldata0[[ifunc[i]]],1:numbasis)
          icomp <- colnames(basis1$coefs)
        },
        "pls"={
          basis1<-create.pls.basis(ldata0[[ifunc[i]]],ldata0[[y]],1:numbasis)
          icomp <- colnames(basis1$coefs)
        },"basis"={
          #numbasis.opt=TRUE
          basis1<-create.fdata.basis(ldata0[[ifunc[i]]],
                                     1:numbasis,type.basis=type.basis)
          xx <- fdata.cen(ldata0[[ifunc[i]]])
          xmean = xx[[2]]
          xcen = xx[[1]]
          x.fd = Data2fd(argvals = ldata0[[ifunc[i]]]$argvals,
                         y = t(xcen$data), 
                         basisobj = basis1)
          
          basis1$coefs = t(x.fd$coefs)          
          icomp <-  basis1$names
          }
    )
      # print(4)
     
      colnames(basis1$coefs)<-paste(ifunc[i],".",icomp,sep="") 
      xposibles <- c(xposibles,colnames(basis1$coefs)[1:numbasis])
      xdata<-cbind(xdata,basis1$coefs[,1:numbasis,drop=F])
    }
  }
  #print(8)
  #print(xposibles)
  iposibles<-1:length(xposibles)
  maxvar <- length(xposibles)
  iposibles <- 1:maxvar
  form.nl <- paste(y,"~")
  nam.model<-"glm"    
  par.model<-list()
  ibest <- 1

  form.ini <-  paste(y,"~ 1")
  # print("best.model")
  best.model <- glm(as.formula(form.ini),data=xdata,family=family)#,weights=weights)
  # print("best.model");  print(best.model)
  nvar <- 1
  res<-NULL
  best.aic<-Inf
  parar <- FALSE
  imodelo<-xmodelo<-NULL
  while (!parar) {
     # print("while");  print(nvar)
    imax <- 1:(maxvar-nvar+1)
    form.nl <- form.ini
    lmax <- length(imax)
    aic <- numeric(lmax)
    for (k in 1:lmax) {  
    # print("k")  ; print(k)
    
   #  print(xposibles[k])
    #if (!numbasis.opt) 
    ik<-1:k
    #else ik<-k
    if (is.factor(xdata[[xposibles[k]]]))   
      fpredictors.nl <- paste("~.+",xposibles[k],sep="",collapse="+")
    #     fpredictors.nl <- paste("~.+",xposibles[ik],sep="",collapse="+")
    else {
      fpredictors.nl <- paste(xposibles[ik],sep="",collapse="+")
      #     fpredictors.nl <- paste("s(",xposibles[ik],",k=",kbs,")",sep="",collapse="+")
         fpredictors.nl <- paste0("~.+", fpredictors.nl)
         }
         # print(form.ini)
         # print(ik)
         # print(fpredictors.nl)
         
       # print("fpredictors.nl")
    #if (!is.na(fpredictors.nl))
      form.nl <-update.formula(as.formula(form.ini), as.formula(fpredictors.nl))
    par.model$formula<-form.nl #as.formula(form.nl)
    par.model$data<-xdata
    par.model$family<- family
    res[[k]]<-do.call(nam.model,par.model)             
    aic[k] <- res[[k]]$aic
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
      if (numbasis.opt){
        nvar <- posible.aic
        xmodelo <- c(xmodelo,xposibles[1:posible.aic]) 
        imodelo <- c(imodelo,iposibles[1:posible.aic]) 
        xposibles <- xposibles[-c(1:posible.aic)]
        iposibles <- iposibles[-c(1:posible.aic)]
      } else {
        xmodelo <- c(xmodelo,xposibles[posible.aic]) 
        imodelo <- c(imodelo,iposibles[posible.aic]) 
        xposibles <- xposibles[-posible.aic]
        iposibles <- iposibles[-posible.aic]
      }
      form.ini <-   res[[posible.aic]]$formula
      }    else {      parar=TRUE    }
    #cat("nvar");print(nvar)
    if (nvar==numbasis) {      parar=TRUE    }    else { nvar<- nvar+1}
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
 #JJ = NULL
 z <- best.model
 if (lenfunc) {
   z$mean = mean.list
   z$basis.x[[vfunc[i]]] = basis1
   #z$basis.b[[vfunc[i]]] = basis1
   #z$JJ <- JJ
   z$vs.list = vs.list
   z$vfunc <- vfunc
   colnames(z$basis.x[[vfunc[i]]]$coefs)<-paste("PC",1:ncol(z$basis.x[[vfunc[i]]]$coefs),sep="")
 }
 # print("finnn")
 if (is.null(imodelo)) {
   warning("No component of the variable is significant") 
 }
 #else{
 
  z$numbasis.opt<-imodelo
  z$formula <- as.formula(z$formula) 
  z$formula.ini =    as.formula( paste(y,"~+",x,sep=""))
# print(summary(z)); print(z$formula) ; print(z$formula.ini) 
  z$data = data
  z$XX = best.model$model
#print(colnames(z$XX)) ;print("optimas");print(z$pc.opt)
# z$nnf <- nnf #var no  funcionales
 class(z) <- c(class(z), "fregre.glm")
 #print("sale fregre.glm.cv") 
 z$call <- z$call[1:2]
 return(z)
}
###################################################
# is.na.fdata <- fda.usc:::is.na.fdata
# pvalue.anova <- fda.usc:::pvalue.anova
# res2=fregre.glm.cv(data=ldata,x="x",y="Fat")
# summary(res2)
# res2=fregre.glm.cv(data=ldata,x="x",y="Fat",numbasis=28)
# summary(res2)
# 
# res2=fregre.glm.cv(data=ldata,y="Fat",numbasis.opt = T,numbasis=28)
# summary(res2)
# 
# 
# res2=fregre.glm.cv(data=ldata,y="Fat",numbasis.opt = T,numbasis=28,type.basis="fourier")
# summary(res2)
# 
# res2=fregre.glm.cv(data=ldata,y="Fat",numbasis.opt = F,numbasis=28,type.basis="fourier")
# summary(res2)
# res2$call[1:2]
# 
# data(tecator)
# x=tecator$absorp.fdata
# y=tecator$y$Fat
# tt=x[["argvals"]]
# dataf=as.data.frame(tecator$y)
# nbasis.x=11
# nbasis.b=7
# basis1=create.bspline.basis(rangeval=range(tt),nbasis=nbasis.x)
# basis2=create.bspline.basis(rangeval=range(tt),nbasis=nbasis.b)
# f=Fat~Protein+x
# basis.x=list("x"=basis1)
# basis.b=list("x"=basis2)
# ldata=list("df"=dataf,"x"=x)
# res=fregre.glm(f,family=gaussian(),data=ldata,basis.x=basis.x,
#                basis.b=basis.b)
# summary(res)
