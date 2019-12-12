##############################################
# Internal function (without predict)
##############################################
fregre.gsam.ncomp<-function(form.ini,data,y,covar,alpha=0.05,
                        type.basis="pc",ncomp=8,
                        par.basis=list(),
                        kbs=-1,#por omision, family=gaussian()
                        criterio="sp",fac=1.01
                        ){
#print("entra ncomp")  
  pc.opt<-0
##print(form.ini)
  iposibles<-  iposibles <- 1:ncomp
    xposibles<-covar
  ##print(xposibles)
  nam.model <- "gam"
  par.model <- list("data"=data)
  aic <- numeric(ncomp)
  res <- NULL
  best.aic <- Inf
  imodelo <- xmodelo <- NULL
  best.model <- gam(as.formula(form.ini), data = data)
  
for (k in 1:ncomp) { 
#print("kkkk")  
#print(k)
     fpredictors.nl <- paste("s(",xposibles[1:k],",k=",kbs,")",sep="",collapse="+")
     fpredictors.nl <- paste0("~.+", fpredictors.nl)
     ##print(fpredictors.nl)
     form.nl <-update.formula(as.formula(form.ini), as.formula(fpredictors.nl))
       #as.formula(form.ini,fpredictors.nl)
     ##print(form.nl)
     #print(fpredictors.nl)
     #print(par.model$formula)
     par.model$formula<-form.nl #as.formula(form.nl)
     
        ##print(nam.model)       
        ##print(names(data))
        #print(names(par.model))
       res[[k]]<-do.call(nam.model,par.model)             
##print("res")       
       aic[k] <- do.call(criterio,list(model=res[[k]]))
 
##print(aic)
  if ((aic[k] *fac)< best.aic) {
    ##print("entra mejora");
    ##print(aic[k] );   ##print(best.aic)
    ##print(best.model$formula)
    pc.opt<-1:k

######################  debe hacerse para el modelo inical vs el final con las PC's seleccionadas #######################
   # 
   pmejora <- pvalue.anova(best.model,res[[k]])
  #  ##print(pmejora)
    pmejora<- pmejora[["Pr(>F)"]][2]#na.omit(pmejora[["Pr(>F)"]])
    if (pmejora<alpha & !is.na(pmejora)) {
##print("aaa")      

      best.aic<-aic[k]
      best.model <- res[[k]]
      xmodelo <- c(xmodelo,xposibles[k]) 
      form.ini <-   res[[k]]$formula
##print(form.ini)      
      }
    else {
      parar=TRUE
      ##print(k)
      ##print("no se mejora el modelo anova")
    }
  }
  } 
 #para PC
 vfunc<-xposibles##########error aqui y en  z$formula.ini =   mas abajo
 i<-1
 z<-best.model
 z$basis.x=z$basis.b=vs.list=mean.list=list()
 #basis1$basis<- basis1$basis#[imodelo,,drop=F]
 
 z <- best.model
 ###print("finnn")
  z$pc.opt<-pc.opt
  z$aic<-aic
  
  z$formula <- as.formula(z$formula) 
  z$formula.ini =    as.formula( paste(y,"~+s(",xposibles,",k=",kbs,")",sep=""))
  z$data = data
  z$XX = best.model$model
  class(z) <- c(class(z), "fregre.gsam")
  #print("sale ncomp")    
 return(z)
}
###################################################
# ncomp<-14
# xx<-create.pc.basis(ldata$x.d2,ncomp)$x[,1:ncomp]
# dim(xx)
# dat<-data.frame(Fat=ldata$df$Fat,zz=ldata$df$zz,xx)
# dim(dat)
# dat$zz<-ldata$df$zz
# covar<-names(dat)[-c(1:2)]
# res1=fregre.gsam.ncomp(Fat~zz,dat,"Fat",covar,ncomp=ncomp,fac=1.05)
# 
# summary(res1)
# res1$pc.opt
# res1$aic
# names(res1$XX)
# plot(res1$aic)
