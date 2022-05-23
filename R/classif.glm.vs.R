classif.glm.vs = function(data = list(), y, x, family = binomial(),classif="classif.glm"
                         , weights = "equal"
                         , basis.x = NULL, basis.b = NULL
                         ,type = "1vsall", prob = 0.5, alpha=0.05
                         ,dcor.min =0.01, smooth=TRUE
                        ,measure= "accuracy",xydist,...){
# print("ini")
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
# print("1")   
  #modelo <-classif.gsam(formula=form, family = family, data=data,basis.x=basis.x,...)
  par.classif=list(formula=form, data=data,basis.x=basis.x)     	 
# print(names(par.classif))  
  
  modelo <- do.call(classif,par.classif)
# print("2")  
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
# print("a1")
if (missing(xydist)) {
#  print("entra xydist")
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
#cat("iiiii ");print(iiii)
  iname<-NULL
  xydist[[resp]]<-metric.dist(res)
  len.var<-length(covar)
  for (i in 1:len.var){
      icovar<-covar[i]
  # print(i)
  # print(n)
  # print(dim(xydist[[resp]]))
  # print(dim(xydist[[icovar]]))
  # print(resp)
  # print(icovar)
      tt=dcor.test(xydist[[resp]],xydist[[icovar]],n=n)
# print("peta")      
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
#print(esfuncional)
#print(dcor[j,])
#print(Mset) ;print(jj)
#print(colnames(dcor)) 
 
# cat(paste0("Iter:",j," dcor(",covar0[jj],")= ",round(dcor[j,jj],4),"\n"))  
 if (max(dcor[j,],na.rm=TRUE) > dcor.min ) {
		modant=modelo
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
    
# print("es factor o !smooth")
		form.lin=update.formula(form, paste0(".~.+",iname))
#old		  modelo=fregre.gsam(formula=form.lin,family = family, data=data,basis.x=basis.x,...)
		par.classif<-list(formula=form.lin,data=data,basis.x=basis.x)
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
# print("5")
#	#print(names(modelo))
##########	 
# aa=do.call(msc,list("model"=modant))
# ab=do.call(msc,list("model"=modelo))

aa <-  cat2meas(modant$group,modant$group.est, measure = measure)
ab <-  cat2meas(modelo$group,modelo$group.est, measure = measure)

#print(msc)
#cat(aa,ab,(ab/aa)<=1,"\n")
#	#print(names(modelo))
 #print(" #################      modelo anterior")
# print(summary(modant))
# print(" #################      modelo actual")
# print(summary(modelo))
if ((ab/aa) <= 1){
  # print("ni entra ni actualiza")
	form=formant
	modelo=modant	
#	modfinish=TRUE
} else {
  if (!modfinish) {
    # print("Entra y actualiza")    
#print(iname)
    i.predictor<-c(i.predictor,iname)
    ipredictor[aux] <- 1
  }
#  res<-NULL
  #for (i in 1:length(modelo$fit))
   # res<-cbind(res,modelo$fit[[i]]$residuals)
# print(head(res0))  
# print(head(modelo$prob.group))
# si 2 grupos solo cogemos 1 grupo
  if (ncol(res0)==1)  res<-res0-modelo$prob.group[,1]
  else res<-res0-modelo$prob.group
  }
# print(head(res))
# print("aaa")
# print(7)
	 if (ab==1) {modfinish=TRUE} # prob.classif=1!
}

#return(list(form=form, data=data, basis.x=basis.x, model=modelo,
#            dcor=dcor,i.predictor=ipredictor,ipredictor=i.predictor))
modelo$dcor=dcor
modelo$i.predictor=ipredictor
modelo$ipredictor=i.predictor


return(modelo)
}

# res.vs <- classif.glm.vs(data = ltrain, "class", 
#                         covar100, classif="classif.glm")
# res.vs$ipredictor
# summary(res.vs)
