#######################
classif.gsam2boost=function(group,fdataobj,family=binomial(),
                            weights=NULL, basis.x=NULL,basis.b=NULL,
                            par.gsam=NULL,CV=FALSE,...){
 C <- match.call()
 a <- list()
 mf <- match.call(expand.dots = FALSE)
 m <- match(c( "group","fdataobj","family","basis.x","basis.b","par.gsam","CV"), names(mf), 0L)
 numg <- nlevels(as.factor(group)) 
 if (is.fdata(fdataobj))   {
   gsam<-C[[3]];
   if (is.null(par.gsam$func)&is.null(par.gsam$k))   {gsam<-paste("s(X)",sep="")}
   else  {
    if (!is.null(par.gsam$func)&!is.null(par.gsam$k))
        gsam<-paste(par.gsam$func,"(X,k=",par.gsam$k,")",sep="")
    else if (is.null(par.gsam$func)&!is.null(par.gsam$k))
      gsam<-paste("s(X,k=",par.gsam$k,")",sep="")
         else         if (!is.null(par.gsam$func)&is.null(par.gsam$k))
         gsam<-paste(par.gsam$func,"(X)",sep="")
  }
  pf2<-as.formula(paste("y~", gsam,sep = ""))
  dataf<-data.frame("y"=group)
  ldata<-list("df"=dataf,"X"=fdataobj)
  X<-C[[3]]       
  }   else {   
      if (is.null(par.gsam$func)) par.gsam$func<-"s"
      if (is.null(par.gsam$k)) par.gsam$k<--1        
      dataf<-data.frame("y"=group,fdataobj)
      ldata<-list("df"=dataf)
      if (is.matrix(fdataobj)) nms<-colnames(fdataobj)
      else nms<-names(fdataobj)               
      gsam<-paste("+",par.gsam$func,"(",nms,",k=",par.gsam$k,")",sep="",collapse="")
      pf2<-as.formula(paste("y~", gsam,sep = ""))
      X<-C[[3]]              
  }   
 newy <- y <- ldata$df$y
 if (!is.factor(y)) y <- as.factor(y)
 n <- length(y);
 newdata <- ldata
 ngroup<-nlevels(y)
 lev<-levels(y)
 prob<-ngroup<-length(table(y))
 if (!is.null(basis.x)) basis.x=list("X"=basis.x)
 if (!is.null(basis.b)) basis.b=list("X"=basis.b)
 formula<-pf2
 lev <-levels(y)
 if (ngroup==2) {
      #lev<-(as.numeric(names(table(y)))
      newy<-ifelse(y==lev[1],0,1)
      newdata$df$y<-newy
      a[[1]]<-fregre.gsam(formula,data=newdata,family=family,weights=weights,basis.x=basis.x,
      basis.b=basis.b,CV=CV,...)
      yest<-ifelse(a[[1]]$fitted.values<.5,lev[1],lev[2])
      yest<-factor(yest,levels=lev)
      tab<-table(yest,y)
      prob[1]=tab[1,1]/sum(tab[,1])
      dtab<-dim(tab)
      if (dtab[1]==dtab[2])    {
           prob[2]=tab[2,2]/sum(tab[,2])
           names(prob)<-lev
         }
      else prob[2]<-0
      prob.group<-a[[1]]$fitted.values
 } else {
  prob.group<-array(NA,dim=c(n,ngroup))
  colnames(prob.group)<-lev
  for (i in 1:ngroup) {
      newy<-ifelse(y==lev[i],0,1)
      newdata$df$y<-newy
      a[[i]]<-fregre.gsam(formula,data=newdata,family=family, 
                          weights=weights,basis.x=basis.x,
                          basis.b=basis.b,CV=CV,...)
      prob.group[,i]<-a[[i]]$fitted.values
  }
  yest<-lev[apply(prob.group,1,which.min)]######no sera which.max
  yest<-factor(yest,levels=lev)
  tab<-table(yest,y)
  for (i in 1:ngroup) {     prob[i]=tab[i,i]/sum(tab[,i])     }
  names(prob)<-lev
 }
 max.prob=sum(diag(tab))/sum(tab)
 output<-list(fdataobj=fdataobj,group=y,group.est=as.factor(yest),
 prob.classification=prob,prob.group=prob.group,C=C,m=m,max.prob=max.prob,
 formula=formula,data=newdata)
 output$fit <- a
 class(output) <- "classif"
 return(output)
}


#############################
classif.glm2boost=function(group,fdataobj,family=binomial(),basis.x=NULL,
                           basis.b=NULL,CV=FALSE,...){
  C<-match.call()
  a<-list()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("group","fdataobj","family","basis.x","basis.b","CV"), names(mf), 0L)
  if (is.fdata(fdataobj))   {
    dataf<-data.frame("y"=group)
    ldata<-list("df"=dataf,"X"=fdataobj)
    formula<-as.formula("y~X")       
  }
  else {
    if (is.matrix(fdataobj)) nms<-colnames(fdataobj)
    else nms<-names(fdataobj)        
    #quita [ o $ del names
    #for (i in 1:length(nms)) nms[i]<-unlist(strsplit(nms[i], "[$]"))[[1]]
    #for (i in 1:length(nms)) nms[i]<-unlist(strsplit(nms[i], "[[]"))[[1]]
    dataf<-data.frame("y"=group,fdataobj)
    ldata<-list("df"=dataf)
    aaa<-paste(nms,collapse="+")
    formula<-as.formula(paste("y~",aaa,sep=""))      
  }
  newy<-y<-ldata$df$y
  if (!is.factor(y)) y<-as.factor(y)
  n<-length(y)
  newdata<-ldata
  ngroup<-nlevels(y)
  lev<-levels(y)
  prob<-rep(NA,ngroup)
  if (!is.null(basis.x)) basis.x=list("X"=basis.x)
  if (!is.null(basis.b)) basis.b=list("X"=basis.b)  
  #nlevels(y) cambiar todos
  a<-list()
  if (ngroup==2) {
    #      lev<-as.numeric(names(table(y)))
    newy<-ifelse(y==lev[1],0,1)
    newdata$df$y<-newy
     a[[1]]<-fregre.glm(formula,data=newdata,family=binomial,basis.x=basis.x,
                       basis.b=basis.b,CV=CV)                

#    a[[1]]<-fregre.gsam(formula,data=newdata,family=binomial,basis.x=basis.x,
#                       basis.b=basis.b,CV=CV)                
    yest<-ifelse(a[[1]]$fitted.values<.5,lev[1],lev[2])
    tab<-table(yest,y)
    prob[1]=tab[1,1]/sum(tab[,1])
    dtab<-dim(tab)
    if (dtab[1]==dtab[2])    {
      prob[2]=tab[2,2]/sum(tab[,2])
      names(prob)<-lev
    }
    else prob[2]<-0
    prob.group<- a[[1]]$fitted.values
    yest<-factor(yest,levels=lev)
  }
  else {
    #   lev<-as.numeric(names(table(y)))
    prob.group<-array(NA,dim=c(n,ngroup))
    colnames(prob.group)<-lev
    for (i in 1:ngroup) {
      newy<-ifelse(y==lev[i],0,1)
      newdata$df$y<-newy
      a[[i]]<-fregre.glm(formula,data=newdata,family=family,basis.x=basis.x,
                         basis.b=basis.b)
#      a[[i]]<-fregre.gsam(formula,data=newdata,family=family,basis.x=basis.x,
#                         basis.b=basis.b)
      prob.group[,i]<-a[[i]]$fitted.values
    }
    yest<-lev[apply(prob.group,1,which.min)]
    yest<-factor(yest,levels=lev)
    tab<-table(yest,y)
    for (ii in 1:ngroup) {prob[ii]=tab[ii,ii]/sum(tab[,ii])}
    names(prob)<-lev
  }
  max.prob=sum(diag(tab))/sum(tab)
  output<-list(fdataobj=fdataobj,group=y,group.est=as.factor(yest),
               prob.classification=prob,prob.group=prob.group,C=C,
               m=m,max.prob=max.prob,formula=formula,data=newdata)
  output$fit<-a  
  class(output) <- "classif"
  return(output)
}

#######################
classif.rpart2boost=function(group,fdataobj,basis.x=NULL,basis.b=NULL,...){   
if (!is.factor(group)) group<-as.factor(group)
   C<-match.call()
   mf <- match.call(expand.dots = FALSE)
   m <- match(c("group","fdataobj","basis.x"), names(mf), 0L)

#  if (is.null(basis.x))  {
#   basis.x2<-create.pc.basis(fdataobj,l=1:5)
#   X<-basis.x<- basis.x2$x
#   }
#else {
#    basis.x2<-basis.x
#    X<-switch(basis.x$type,
#    "pc"= basis.x$x,
#    "raw"=basis.x$basis$data,
#    "bspline"=t(Data2fd(argvals = argvals(fdataobj), y = t(fdataobj$data),basisobj = basis.x)$coef)
#    )
#}
   if (!is.null(basis.x))  basis.x<-list("X"=basis.x)
   if (!is.null(basis.b))  basis.b<-list("X"=basis.b)     
   dataf<-data.frame("y"=group)
   ldata<-list("df"=dataf,"X"=fdataobj)
#   dataf<-list("y"=group,"X"=X)
   formula<-formula(y~X)     
#   basis.x<-list("X"=basis.x2)
   fit<-classif.rpart(formula,data=ldata,basis.x=basis.x,...)       
  output<-list(fit=fit,formula=formula,fdataobj=dataf,basis.x=basis.x,
  basis.b=basis.b,group=group,C=C,max.prob=fit$max.prob,"prob.group"=fit$prob.group,
  "prob.classification"=fit$prob.classification,"group.est"=fit$group.est)
class(output) <- "classif"
return(output)
}



################################################################################
classif.gkam2boost=function(group,fdataobj,family = binomial(),weights=rep(1,n),
par.metric = NULL,par.np=NULL, offset=NULL,
control = list(maxit = 100,epsilon = 0.001, trace = FALSE,inverse="solve"),...)  {
formula<-formula(y~X)
C<-match.call()
a<-list()
mf <- match.call(expand.dots = FALSE)
m <- match(c("group","fdataobj","family","weights","par.metric","par.np","offset",
"control"), names(mf),0L)
   dataf<-data.frame("y"=group)
   ldata<-list("df"=dataf,"X"=fdataobj)
   newy<-y<-ldata$df$y
   if (!is.factor(y)) y<-as.factor(y)
   n<-length(y)
   newdata<-ldata
   ngroup<-nlevels(y)
   prob<-rep(NA,ngroup)
   lev<-levels(y)
#newdata<-data
if (!is.null(par.np)) {
   par.np=list("X"=par.np)
   if (is.null(par.np[["X"]]$Ker)) par.np[["X"]]$Ker=AKer.norm
   if (is.null(par.np[["X"]]$type.S)) par.np[["X"]]$type.S="S.NW"
   }
else          par.np =list("X"=list(Ker=AKer.norm,type.S="S.NW"))

if (ngroup==2) {
      newy<-ifelse(y==lev[1],0,1)
      newdata$df$y<-newy
      a[[1]]<-fregre.gkam(formula,data=newdata,family=family,weights=weights,
      par.metric=par.metric,par.np=par.np,offset=offset,control=control,...)
      yest<-ifelse(a[[1]]$fitted.values<.5,lev[1],lev[2])
            tab<-table(yest,y)
      prob[1]=tab[1,1]/sum(tab[,1])
      dtab<-dim(tab)
      if (dtab[1]==dtab[2])    {
           prob[2]=tab[2,2]/sum(tab[,2])
           names(prob)<-lev
         }
      else prob[2]<-0
      prob.group<-a$fitted.values
      yest<-factor(yest,levels=lev)
   }
else {
#   lev<-levels(y)
   prob.group<-array(NA,dim=c(n,ngroup))
   colnames(prob.group)<-lev
   for (i in 1:ngroup) {
              newy<-ifelse(y==lev[i],0,1)
              newdata$df$y<-newy
              a[[i]]<-fregre.gkam(formula,data=newdata,family=family,weights=weights,
              par.metric=par.metric,par.np=par.np,offset=offset,control=control,...)
              prob.group[,i]<-a[[i]]$fitted.values
              }
   yest<-lev[apply(prob.group,1,which.min)]#no sera which.max
   yest<-factor(yest,levels=lev)
   tab<-table(yest,y)
   for (ii in 1:ngroup) {
       prob[ii]=tab[ii,ii]/sum(tab[,ii])
       }
  names(prob)<-lev
}
max.prob=sum(diag(tab))/sum(tab)
output<-list(fdataobj=fdataobj,group=y,group.est=yest,
prob.classification=prob,prob.group=prob.group,C=C,m=m,max.prob=max.prob,fit=a
,formula=formula,data=newdata)
class(output) <- "classif"
return(output)
}

#######################
#######################
classif.glm.cv=function(formula,data,w=NULL,family = binomial(),
basis.x=NULL,basis.b=NULL,CV=FALSE,...){
C<-match.call()
a<-list()
mf <- match.call(expand.dots = FALSE)
m <- match(c( "formula","data","w","family","basis.x","basis.b","CV"), names(mf), 0L)
 tf <- terms.formula(formula)
 terms <- attr(tf, "term.labels")
 nt <- length(terms)
 if (attr(tf, "response") > 0) {
        response <- as.character(attr(tf, "variables")[2])
        pf <- rf <- paste(response, "~", sep = "")
    } else pf <- rf <- "~"
newy<-y<-data$df[[response]]
if (!is.factor(y)) y<-as.factor(y)
n<-length(y)
lev<-levels(y)
newdata<-data
prob2<-prob<-ngroup<-length(table(y))
if (ngroup==2) {
      lev<-as.numeric(names(table(y)))
      newy<-ifelse(y==lev[1],0,1)
      newdata$df$y<-newy
      a[[1]]<-fregre.glm(formula,data,family=family,basis.x=basis.x,basis.b=basis.b,CV=CV,...)

      if (CV) prediction<-a[[1]]$pred.cv
      else prediction<-a[[1]]$fitted.values
      yest<-ifelse(prediction<.5,0,1)
      yest<-factor(yest,levels=lev)
      tab<-table(yest,y)
   prob[1]=tab[1,1]/sum(tab[,1])
      dtab<-dim(tab)
      if (dtab[1]==dtab[2])    {
           prob[2]=tab[2,2]/sum(tab[,2])
           names(prob)<-lev
         }
      else prob[2]<-0
      prob.group<-prediction

      #devolver a mayores y estimada!
   }
else {
   lev<-as.numeric(names(table(y)))
   prob.group<-array(NA,dim=c(n,ngroup))
   colnames(prob.group)<-lev
   for (i in 1:ngroup) {
              newy<-ifelse(y==lev[i],0,1)
              newdata$df$y<-newy
              a[[i]]<-fregre.glm(formula,data,family=family,basis.x=basis.x,basis.b=basis.b,CV=CV,...)
      if (CV) prediction<-a[[1]]$pred.cv
      else prediction<-a[[1]]$fitted.values
            }
   yest<-lev[apply(prob.group,1,which.min)]
   yest<-factor(yest,levels=lev)
   tab<-table(yest,y)
   for (i in 1:ngroup) { prob[i]=tab[i,i]/sum(tab[,i])     }
           names(prob)<-lev
}
lev<-levels(y)
max.prob=sum(diag(tab))/sum(tab)
output<-list(formula=formula,data=data,group=y,group.est=factor(yest,levels=lev),
prob.classification=prob,prob.group=prob.group,C=C,fit=a,m=m,max.prob=max.prob)
class(output)="classif"
return(output)
}
