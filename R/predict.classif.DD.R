
#' Predicts from a fitted classif.DD object.
#' 
#' Classifier of functional (and multivariate) data by DD--classifier.
#' 
#' @details Returns the groups or classes predicted using a previously trained model.
#' 
#' @param object Object \code{object} estimated by \code{classif.DD}.
#' @param new.fdataobj By default, new p functional explanatory dataset or new
#' mulitvariate data of \code{data.frame} class
#' @param type !=''predictive'', for each row of data shows the probability of
#' each group membership.
#' @param \dots Further arguments passed to or from other methods.
#' @return 
#' \itemize{
#' \item \code{group.pred}:Vector of groups or classes predicted
#' \item \code{prob.group}: For each functional data shows the probability of each group membership.
#' }
#' @author Febrero-Bande, M., and Oviedo de la Fuente, M.
#' @seealso See also \code{\link{classif.DD}} .
#' @references Li, J., P.C., Cuesta-Albertos, J.A. and Liu, R.
#' \emph{DD--Classifier: Nonparametric Classification Procedure Based on
#' DD-plot}.  Journal of the American Statistical Association (2012), Vol. 107,
#' 737--753.
#' @keywords classif
#' @examples
#' \dontrun{
#' # DD-classif for multivariate data 
#' data(iris)
#' iris<-iris[1:100,]
#' ii<-sample(1:100,80)
#' group.train<-factor(iris[ii,5])
#' x.train<-iris[ii,1:4]
#' out1=classif.DD(group.train,x.train,depth="MhD",classif="lda")
#' out2=classif.DD(group.train,x.train,depth="MhD",classif="glm")
#' summary(out1)
#' summary(out2)
#' x.test<-iris[-ii,1:4]
#' pred1=predict(out1,x.test)
#' pred2=predict(out2,x.test)
#' group.test<-iris[-ii,5]
#' table(pred1,group.test)
#' table(pred2,group.test)
#' 
#' # DD-classif for Functional data
#' data(phoneme)
#' mlearn<-phoneme[["learn"]]
#' glearn<-phoneme[["classlearn"]]
#' 
#' #	ESTIMATION 
#' out1=classif.DD(glearn,mlearn,depth="FM",classif="glm")
#' summary(out1)
#' #	PREDICTION 
#' mtest<-phoneme[["test"]]
#' gtest<-phoneme[["classtest"]]
#' pred1=predict(out1,mtest)
#' table(pred1,gtest)
#' }
#' @export
predict.classif.DD<-function(object,new.fdataobj=NULL,type="class",...){
#print("pred classif DD")  
  if (is.null(new.fdataobj)) return(object$group.est)
  depth<-object$depth
  par.depth<-object$par.depth
  multi<-object$multi                           
  if (is.fdata(new.fdataobj)| is.matrix(new.fdataobj) | is.data.frame(new.fdataobj)) {  ldata<-list(new.fdataobj)}
  else {   ldata<-new.fdataobj   }
  lenlista<-length(ldata)
  lendepth<-length(depth)
  par.ldata<-object$par.depth
  w<-object$w
  lenpesos<-length(w)
  depth<-depthl<-object$depth
  ng2<-ncol(object$dep)
  par.Df<-list()  
  group<-object$group
  nn<-length(group)
  par.depth<-object$par.depth  
  lev<-levels(group)
  ng<-length(lev)
  ind<-matrix(NA,nrow=nn,ncol=ng)
  integrated<-FALSE   
  
  if (missing(w)) {   
    model<-TRUE
    lenl<-lenlista
    w<-rep(1/lenlista,len=lenlista)
    ng2<-ng*lenlista
  }
  else {
    if (depth[1]=="modep") {
      integrated<-TRUE
      n<-nrow(ldata[[1]])
      Df<-matrix(NA,ncol=ng,nrow=n)
      hq<-par.depth$h
      par.depth$lfdata<-new.fdataobj
      for (i in 1:ng) {
        ind[,i]<-group==lev[i]
        par.depth$lfdataref<-c.ldata(object$par.depth$lfdata,ind[,i])
        par.depth$h<-hq[i] #misma ventana para toodos los datos
        aa<-do.call("depth.modep",par.depth)$dep        #call a la funcion depth
        Df[,i]<-aa
      }    
    }
    if (depth[1]=="FMp" | depth[1]=="RPp"){ 
      integrated<-TRUE
      n<-nrow(new.fdataobj[[1]])
      Df<-matrix(NA,ncol=ng,nrow=n)       
      par.depth$lfdata<-new.fdataobj
      nam.dep<-paste("depth.",depth[1],sep="")
      for (i in 1:ng) {
        ind[,i]<-group==lev[i]
        par.depth$lfdataref<-c.ldata(object$par.depth$lfdata,ind[,i])
        Df[,i]<-do.call(nam.dep,par.depth)$dep        
      }    
    }  
  }     
  
  if (!integrated){
#print("no integrated")    
    for (idat in 1:lenlista) {
      depth<-depthl[idat]
      par.depth<-par.ldata[[idat]]
      fdataobj<-par.depth$fdataobj
      new.fdataobj<-ldata[[idat]]
      nc<-ncol(fdataobj)
      if (multi) {
        depth.long<-paste("mdepth.",depth,sep="")
        if (is.vector(new.fdataobj)) new.fdataobj<-par.depth$x<-matrix(new.fdataobj,nrow=1)
        else par.depth$x<-new.fdataobj
        fdataobj<-object$par.depth[[idat]]$x
      }
      else    {
        depth.long<-paste("depth.",depth,sep="")
        if (is.vector(new.fdataobj[["data"]])) new.fdataobj[["data"]]<-matrix(new.fdataobj,nrow=1)
        par.depth[["fdataobj"]]<-new.fdataobj
        fdataobj<-object$par.depth[[idat]]$fdataobj
      }    
      n<-nrow(new.fdataobj)
      nc<-ncol(new.fdataobj)
      
      nvec<-table(group)
      p<-nvec[1]/n    
      Df<-matrix(NA,ncol=ng,nrow=n) 
      ismdist<-is.matrix(par.depth$metric)
      if (ismdist) {
        mdist<-par.depth$metric
        par.depth$metric<-attr(mdist, "call")
      }
      dmode<-c(depth.long=="depth.mode" | depth.long=="mdepth.mode")
      if (dmode)  hq<-par.depth$h
      for (i in 1:ng) {
        ind[,i]<-group==lev[i]
        if (multi)  par.depth$xx<-fdataobj[ind[,i],]
        else   par.depth$fdataori<-fdataobj[ind[,i],]
        if (dmode)      par.depth$h<-hq[i]
        Df[,i]<-do.call(depth.long,par.depth)$dep
      } 
      if (object$model) {
        if (idat==1) Df2<-Df
        else  Df2<-cbind(Df2,Df)
        par.ldata[[idat]]<-par.depth
        par.Df[[paste("dep",idat,sep="")]]<-fdata(Df)
      }
      else{
        if (idat==1) Df2<-w[idat]*Df
        else  Df2<-Df2+w[idat]*Df
        par.ldata[[idat]]<-par.depth
      }
    } 
    Df<-Df2
  }
  colnames(Df)<-colnames(object$dep)
  #  print(object$classif)
  if (object$classif=="DDk"){
    if (length(object$par.classif$pol)==2) object$classif="DD2"
    if (length(object$par.classif$pol)==3) object$classif="DD3"
  }
  #print(object$classif)
  group.est<-switch(object$classif,
                    # MD={gest<-apply(Df,1,which.max)},
                    DD1={
                      if (ng>2) {#majority voting option 
                        #  stop("DD-plot for more than 2 levels not available")
                        warning("Majority voting classification")
                        #ojo ng es num de grupos y ng2 ng*ndepth
                        cvot<-combn(ng2,2)
                        nvot<-ncol(cvot)
                        votos<-matrix(0,ng,n)
                        eps<-1e-10
                        Df[Df<eps]<-eps
                        for (ivot in 1:nvot) {      
                          #       cat("votando",ivot)
                          #       i2a2<-which(group==lev[cvot[1,ivot]] | group==lev[cvot[2,ivot]] )
                          #Df0<-Df[i2a2,cvot[,ivot]]
                          Df0<-Df[,cvot[,ivot]]
                          #       ind0<-ind[i2a2,cvot[,ivot]]
                          #       b<-unique(Df0[,1]/Df0[,2])       
                          #       mis <- sapply(b,MCR0,Df0,ind0)
                          #       b0 <- min(b[which.min(mis)])
                          b0<-object$par.classif$pol[[ivot]]
                          group.log<-b0*Df0[,1]<Df0[,2]
                          votos[cvot[1,ivot],]<-votos[cvot[1,ivot],]+as.numeric(!group.log)
                          votos[cvot[2,ivot],]<-votos[cvot[2,ivot],]+as.numeric(group.log)    
                        }
                        maj.voto<-apply(votos,2,which.max)
                        group.est<-maj.voto
                        for (ii in 1:n) {
                          l = seq_along(votos[,ii])[votos[,ii] == max(votos[,ii], na.rm = TRUE)]      
                          if (length(l) > 1) {
                            abc<-which(Df[ii,]== max(Df[ii,l ], na.rm = TRUE))  
                            group.est[ii] =group[abc]
                          }
                          group.est <-  factor(group.est,levels = lev)
                        }
                        group.est
                      } #######  fin   voting en prediccion #####################   
                      else group.est<- factor(ifelse(object$par.classif$pol*Df[,1]>Df[,2],lev[1],lev[2]),levels=lev)},   
                    DD2={
                      if (ng>2) {#majority voting option 
                        #  stop("DD-plot for more than 2 levels not available")
                        warning("Majority voting classification")     
                        cvot<-combn(ng2,2)
                        nvot<-ncol(cvot)
                        votos<-matrix(0,ng,n)
                        eps<-1e-10
                        Df[Df<eps]<-eps
                        for (ivot in 1:nvot) {      
                          #       cat("votando",ivot)
                          Df0<-Df[,cvot[,ivot]]
                          b0<-object$par.classif$pol[[ivot]]
                          #group.log<-b0*Df0[,1]<Df0[,2]
                          group.log<-sapply(Df0[,1],RR,a=b0)<Df0[,2]   
                          votos[cvot[1,ivot],]<-votos[cvot[1,ivot],]+as.numeric(!group.log)
                          votos[cvot[2,ivot],]<-votos[cvot[2,ivot],]+as.numeric(group.log)    
                        }
                        maj.voto<-apply(votos,2,which.max)
                        group.est<-maj.voto
                        for (ii in 1:n) {
                          l = seq_along(votos[,ii])[votos[,ii] == max(votos[,ii], na.rm = TRUE)]      
                          if (length(l) > 1) {
                            abc<-which(Df[ii,]== max(Df[ii,l ], na.rm = TRUE))  
                            group.est[ii] =group[abc]
                          }
                          group.est <-  factor(group.est,levels = lev)
                        }
                        group.est
                      } #######  fin   voting en prediccion #####################   
                      else {
                        if (object$par.classif$rotate) group.est<- factor(ifelse(sapply(Df[,2],RR,a=object$par.classif$pol)>Df[,1],lev[2],lev[1]),levels=lev)
                        else group.est<- factor(ifelse(sapply(Df[,1],RR,a=object$par.classif$pol)>Df[,2],lev[1],lev[2]),levels=lev)
                      }
                    },
                    DD3={
                      if (ng>2) {#majority voting option 
                        #  stop("DD-plot for more than 2 levels not available")
                        warning("Majority voting classification")
                        cvot<-combn(ng,2)
                        nvot<-ncol(cvot)
                        votos<-matrix(0,ng,n)
                        eps<-1e-10
                        Df[Df<eps]<-eps 
                        for (ivot in 1:nvot) {      
                          #cat("votando",ivot)
                          Df0<-Df[,cvot[,ivot]]
                          b0<-object$par.classif$pol[[ivot]]
                          #group.log<-b0*Df0[,1]<Df0[,2]
                          group.log<-sapply(Df0[,1],RR,a=b0)<Df0[,2]       
                          votos[cvot[1,ivot],]<-votos[cvot[1,ivot],]+as.numeric(!group.log)
                          votos[cvot[2,ivot],]<-votos[cvot[2,ivot],]+as.numeric(group.log)    
                        }
                        maj.voto<-apply(votos,2,which.max)
                        group.est<-maj.voto
                        for (ii in 1:n) {
                          l = seq_along(votos[,ii])[votos[,ii] == max(votos[,ii], na.rm = TRUE)]      
                          if (length(l) > 1) {
                            abc<-which(Df[ii,]== max(Df[ii,l ], na.rm = TRUE))  
                            group.est[ii] =group[abc]
                          }
                          group.est <-  factor(group.est,levels = lev)
                        }
                        group.est
                      }
                      else{
                        if (object$par.classif$rotate) group.est<- factor(ifelse(sapply(Df[,2],RR,a=object$par.classif$pol)>Df[,1],lev[2],lev[1]),levels=lev)
                        else group.est<- factor(ifelse(sapply(Df[,1],RR,a=object$par.classif$pol)>Df[,2],lev[1],lev[2]),levels=lev)
                      }
                    },  
                    lda={group.es<-predict(object$fit,Df)$class},
                    qda={group.es<-predict(object$fit,Df)$class},
                    glm={
                      dat<-data.frame(Df)
                      #group.est<-pred2glm2boost(object$fit,list("df"=dat))$gro
                      # print(group.est)
                      group.est<-predict(object$fit,list("df"=dat),type = "class")
                      group.est
                      },
                    gam={
                      dat<-data.frame(Df) 
                      
                      #group.est<-pred2gsam2boost(object,list("df"=dat))
                      group.est<-predict(object$fit,list("df"=dat),type = "class")
                    },
                    rpart={
                      dat<-data.frame(Df)
                      names(dat)<-colnames(object$dep)
                      group.est<-predict(object$fit,dat,type = "class")
                    }, 
                    knn={
                      dat<-data.frame(Df)
                      group.est<-predict.classif(object$fit,dat,type = "class")
                    },
                    np={
                      dat<-data.frame(Df)
                      group.est<-predict.classif(object$fit,dat)
                    },
                    grm={
                      dat<-data.frame(Df)
                      group.est<-predict.classif(object$fit,par.Df,type = "class")
                    }
  )
  
  if (type=="dep") return(list("group.pred"=group.est,"dep"=Df))
  else   return(group.est)
}
###########################################
