classif.ML.vs = function(data = list(), y, x, family = binomial(),
                         classif="classif.svm"
                          , weights = "equal"
                          , basis.x = NULL, basis.b = NULL
                          ,type = "1vsall", prob = 0.5, alpha=0.05
                          ,dcor.min =0.01#, smooth=TRUE
                          ,measure= "accuracy",xydist,...){
  #print("ini")
  if (missing(y)) {stop("The name of the response must be specified in the 'y' argument")  }
  resp <- y
  n=length(data$df[[resp]])
  nesc=length(names(data$df))-1 #Number of scalar variables
  namesc=names(data$df)[which(names(data$df)!=resp)]
  namfunc=names(data)[which(names(data)!="df")]
  
  if (missing(x)) {  
    x<-c(namesc,namfunc)
    ifunc<-namfunc
    infunc<-namesc
  }  else {
    ifunc<-intersect(x,namfunc)
    infunc<-intersect(c(y,x),namesc)
  } 
  
  # xdatos<-as.list(data$df[,y,drop=F])
  # xdatos<-c(xdatos,as.list(data$df[,infunc,drop=F]),data[ifunc])
  
  #nvar=length(names(data$df))-1+length(data)-1 #Global number of variates
  nvar =length(x)
  
  dcor = matrix(0,nrow=nvar,ncol=nvar)
  ind = numeric(nvar)
  colnames(dcor) = x #c(namesc,namfunc)
  covar0 <- covar <- x #c(namesc,namfunc)
  names(ind) = colnames(dcor)
  rownames(dcor) = 1:nvar
  
  form = as.formula(paste0(resp,"~ 1"))#################################
  
  j=1
  Mset=c()
  # print(form)  
  #print("1")   
  #modelo <-classif.gsam(formula=form, family = family, data=data,basis.x=basis.x,...)
  #par.classif=list(formula=form, data=data,basis.x=basis.x)     	 
  #print(names(par.classif))  
  #modelo <- do.call(classif,par.classif)
  
  #res<-NULL
  #for (i in 1:length(modelo$fit))
  #    res<-cbind(res,modelo$fit[[i]]$residuals)
  res0<-res<-model.matrix(as.formula(paste0("~-1+",resp)),data$df)
  # print("3")  
  #print("aa")
  #print(head(res0)  )
  #res=data$df[,resp]-mean(data$df[,resp])
  modfinish=FALSE
  len.esc <- length(namesc)
  len.fun <- length(namfunc)
  i.predictor<-NULL
  ipredictor<-numeric(nvar)
  names(ipredictor)<-x
  #print("a1")
  if (missing(xydist)) {
# print("entra xydist")
    calc.dist<-TRUE
    xydist<-list()
    xydist[[resp]]<-metric.dist(res)  
    if (len.esc>0){
      #  cat(len.esc,"Escalar j=",j,"\n") 
      for (i in 1:len.esc){
        aux<-NA
        if (is.factor(data$df[[namesc[i]]])) {
          y=model.matrix(as.formula(paste0("~-1+",namesc[i])),data$df)
        } else {
          y=data$df[,namesc[i],drop=F]
        }
        #      cat(i,namesc[i],dim(y),"\n")
        xydist[[namesc[i]]]<-metric.dist(y) 
      }
      #    print(names(xydist))
    }
    if (len.fun>0){
      #    cat("Funcional j=",j,"\n")        
      for (i in 1:len.fun){
        xydist[[namfunc[i]]]<-metric.lp(data[[namfunc[i]]])
      }}
    
  }# else  calc.dist<-FALSE
  
  
  #print("a2")
  len.var<-length(covar)
  #print(calc.dist)
  iiii<-0#borrraar
  while(!modfinish){
    iiii<-iiii+1 #borrraar
    # cat("iiiii ");print(iiii)
    iname<-NULL
    # print(res)
    xydist[[resp]]<-metric.dist(res) ##############
    # print(xydist[[resp]][1:4,1:4])
    len.var<-length(covar)
    for (i in 1:len.var){
      icovar<-covar[i]
# print(resp);      print(icovar)
#print(dim(xydist[[resp]]));      print(dim(xydist[[icovar]]))
      tt=dcor.test(xydist[[resp]],xydist[[icovar]],n=n)
      #dcor[j,i]=tt$estimate*(tt$p.value<alpha)
      dcor[j,icovar]=tt$estimate*(tt$p.value<alpha)
    }
    
    #print(max(dcor[j,]))
    #print(which.max(dcor[j,]))
    #dcor[j,covar %in% Mset]=-dcor[j,covar %in% Mset]
    jj=which.max(dcor[j,])
    aux<-jj
    iname<-covar0[jj]
    #print("iname") ;print(iname)
    if  (iname %in% namfunc) esfuncional<-TRUE
    else esfuncional<-FALSE
    # print(esfuncional)
    # print(dcor[j,])
    #print(Mset) ;print(jj)
    #print(colnames(dcor)) 
    
  # cat(paste0("Iter:",j," dcor(",covar0[jj],")= ",round(dcor[j,jj],4),"\n"))  
    # print(iname)
    if (max(dcor[j,],na.rm=TRUE) > dcor.min ) {
      if (j>1)      modant=modelo
      else modant<-NULL
      formant=form
      Mset=c(Mset,iname)
      covar<-setdiff(covar,iname)
      if (!esfuncional){
        # print("escalar")		  
        #namesc[jj]
        esfactor<-is.factor(data$df[,iname])
        namesc<-setdiff(namesc,iname)
        len.esc<-len.esc-1
      } 		else { 
        # print("funcional")
        esfactor<-FALSE
        namfunc<-setdiff(namfunc,iname)
        len.fun<-len.fun-1
      }
      a1<- Inf
      # cat("esfactor " ,esfactor,"\n")
      form.lin=update.formula(form, paste0(".~.+",iname))
      #old		  modelo=fregre.gsam(formula=form.lin,family = family, data=data,basis.x=basis.x,...)
      par.classif<-list(formula=form.lin,data=data,basis.x=basis.x)
#print(form.lin)
      modelo=do.call(classif,par.classif)
      #modelo=classif.gsam(formula=form.lin,family = family, data=data,basis.x=basis.x,...)
      
      #a1=do.call(msc,list("model"=modelo))	
      a1 <-  cat2meas(modelo$group,modelo$group.est, measure = measure)
      form <- form.lin
      #print(form)
      j=j+1
      
      if (j>nvar){modfinish=TRUE} # entran todas o ninguna
    } else {modfinish=TRUE}
    #!modfinish
    # aa=do.call(msc,list("model"=modant))
    # ab=do.call(msc,list("model"=modelo))
    ab <-  cat2meas(modelo$group,modelo$group.est, measure = measure)
    if (j==2){
      form.ant<-form
      modant<-modelo
      aa <-0
    }    else{
    aa <-  cat2meas(modant$group,modant$group.est, measure = measure)    }
    if ((ab/aa) <= 1){
      form=formant
      modelo=modant	
      #	modfinish=TRUE
    } else {
      if (!modfinish) {
      #   print("Entra y actualiza")    
      #  print(iname)
        i.predictor<-c(i.predictor,iname)
        ipredictor[aux] <- 1
      }
      if (ncol(res0)==1)  res<-res0-modelo$prob.group[,1]
      else res<-res0-modelo$prob.group
    }
       if (ab==1) {modfinish=TRUE} # prob.classif=1!
  }
  #return(list(form=form, data=data, basis.x=basis.x, model=modelo,
  #            dcor=dcor,i.predictor=ipredictor,ipredictor=i.predictor))
  modelo$dcor=dcor
  modelo$i.predictor=ipredictor
  modelo$ipredictor=i.predictor
  #devolver DCOR en cada paso
  
  return(modelo)
}

# res.vs <- classif.ML.vs(data = ltrain, "class", 
#                        covar100, classif="classif.svm")
# res.vs$ipredictor
# res.vs$max.prob
# res1 <- classif.ML.vs(ltrain,"class",nam.f[1:14],
#                       classif=classif, xydist=ldist,
#                       ,dcor.min = dc)
# summary(res1)
# res1$i.predictor

