### R code from vignette source 'FDA2.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: setup
###################################################
opts_chunk$set(fig.path='revclass/beamer-',fig.align='center',
				fig.show='hold',size='scriptsize',fig.height=4,
				comment=">",tidy=TRUE,tidy.opts=list(blank=FALSE))

# Dividir salidas grandes, incluir outlines=inicio:fin
hook_output = knit_hooks$get('output')
knit_hooks$set(output = function(x, options) {
  if (!is.null(n <- options$out.lines)) {
    if(length(n)==1) n<-c(n,n)
    x = unlist(stringr::str_split(x, '\n'))
    len <- length(x)
    if (length(x) > n[1]) {
      # truncate the output
      x = x[intersect(seq_along(x), n)]
      if(n[1] >1) x = c('....\n', x)
      if(range(n)[2] < len) x = c(x, '....\n')
    }
    x = paste(x, collapse = '\n') # paste first n lines together
  }
hook_output(x, options)
})


###################################################
### code chunk number 2: 'poblenou'
###################################################
options(width=60) # make the printing fit on the page
library(fda.usc)
#library(fda.classif)
#Mientras no se actuliza el paquete
#source("../../FebreroOviedo/fda.classif/R/predict.classif.R") 
#source("../../FebreroOviedo/fda.classif/R/predict.classif.DD.R") 

set.seed(27031967)
data(poblenou)
dayw=ifelse(poblenou$df$day.week==7 | poblenou$df$day.festive==1,3,ifelse(poblenou$df$day.week==6,2,1))
mc=c("black","red","blue")
plot(poblenou$nox,col=mc[dayw],lwd=1,lty=1)
legend("topleft",c("Workday","Sat","Sun & Festive"),col=mc,lwd=1)


###################################################
### code chunk number 3: 'phoneme'
###################################################
data(phoneme)
mc=c("black","green","cyan","red","blue")
mc2=rgb(t(col2rgb(mc)/255),alpha=0.4)
plot(phoneme$learn,col=mc2[phoneme$classlearn],lwd=0.75,lty=2,ylim=c(0,30))
legend("topright",c("sh","iy","dcl","aa","ao"),fill=mc,cex=0.8)
ldata=list(df=data.frame(g=phoneme$classlearn),learn=phoneme$learn)
lines(func.mean.formula(learn~g,data=ldata),lwd=3,lty=1,col=mc)


###################################################
### code chunk number 4: intro
###################################################
tt=seq(0,1,len=51)
x1=rproc2fdata(100,tt,mu=10*tt*(1-tt)^2,sigma="vexponential")
x2=rproc2fdata(100,tt,mu=-10*tt^2*(1-tt),sigma="vexponential")
x=c(x1,x2)
g=factor(c(rep(1,100),rep(2,100)))
ldata=list(df=data.frame(gr=g),X=x)
b.x=list(X=create.pc.basis(x,1:4))
xnew1=rproc2fdata(5,tt,mu=10*tt*(1-tt)^2,sigma="vexponential")
xnew2=rproc2fdata(5,tt,mu=-10*tt^2*(1-tt),sigma="vexponential")
zero=fdata(rep(0,length(tt)),argvals=tt)
xnew=c(xnew1,xnew2,zero)
ldatanew=list(X=xnew)
mc=c("red","blue")
par(mfrow=c(1,3))
plot(x,lwd=0.5,lty=2,col=mc[g])
lines(c(func.mean(x1),func.mean(x2)),col=mc[1:2],lty=1,lwd=2)
plot(b.x$X$x[,1:2],col=mc[g],pch=19,main="PC1 vs PC2 by group")
plot(xnew,lty=1,col=c(rep(mc[1],5),rep(mc[2],5),"gray"),lwd=2,main="New curves")


###################################################
### code chunk number 5: mv
###################################################
res.svm=classif.svm(gr~X,data=ldata,basis.x=b.x)
res.nnet=classif.nnet(gr~X,data=ldata,basis.x=b.x,trace=FALSE)
res.rpart=classif.rpart(gr~X,data=ldata,basis.x=b.x)
print(paste0("Est. SVM=",round(mean(res.svm$prob.classification),3),", Est. NNet=",round(mean(res.nnet$prob.classification),3),", Est. Tree=",round(mean(res.rpart$prob.classification),3)))
pr.svm=predict(res.svm,ldatanew,type="probs")
pr.nnet=predict(res.nnet,ldatanew, type="probs")
pr.rpart=predict(res.rpart,ldatanew,type="probs")

dfaux0=data.frame(svm.g=pr.svm$group.pred,nnet.g=pr.nnet$group.pred,tree.g=pr.rpart$group.pred,
	svm1=round(pr.svm$prob.group[,1],3),svm2=round(pr.svm$prob.group[,2],3),
	nnet1=round(pr.nnet$prob.group[,1],3),nnet2=round(pr.nnet$prob.group[,2],3),
	tree1=round(pr.rpart$prob.group[,1],3),tree2=round(pr.rpart$prob.group[,2],3))
dfaux0


###################################################
### code chunk number 6: glm
###################################################
res.glm=classif.glm(gr~X,data=ldata,basis.x=b.x)
res.gsam=classif.gsam(gr~s(X),data=ldata,basis.x=b.x)
res.gkam=classif.gkam(gr~X,data=ldata,par.metric=list(X=list(metric=metric.lp,lp=2)))
print(paste0("Est. GLM=",round(mean(res.glm$prob.classification),3),", Est. GSAM=",round(mean(res.gsam$prob.classification),3),", Est. GKAM=",round(mean(res.gkam$prob.classification),3)))
pr.glm=predict(res.glm,ldatanew,type="probs")
pr.gsam=predict(res.gsam,ldatanew, type="probs")
pr.gkam=predict.classif(res.gkam,ldatanew,type="probs")

dfaux=data.frame(glm.g=pr.glm$group.pred,gsam.g=pr.gsam$group.pred,gkam.g=pr.gkam$group.pred,
	glm1=round(pr.glm$prob.group[,1],3),glm2=round(pr.glm$prob.group[,2],3),
	gsam1=round(pr.gsam$prob.group[,1],3),gsam2=round(pr.gsam$prob.group[,2],3),
	gkam1=round(pr.gkam$prob.group[,1],3),gkam2=round(pr.gkam$prob.group[,2],3))
dfaux


###################################################
### code chunk number 7: np
###################################################
res.np=classif.np(g,x,h=seq(.1,.4,len=7))
res.knn=classif.knn(g,x,knn=seq(1,9,by=2))
print(paste0("Optim. h=",res.np$h.opt,", Optim. knn=",res.knn$h.opt))
print(paste0("Est. NP=",round(mean(res.np$prob.classification),3),", Est. KNN=",round(mean(res.knn$prob.classification),3)))
pr.np=predict.classif(res.np,xnew,type="probs")
pr.knn=predict.classif(res.knn,xnew,type="probs")
dfaux2=data.frame(np.g=pr.np$group.pred,knn.g=pr.knn$group.pred,
	np1=round(pr.np$prob.group[,1],3),np2=round(pr.np$prob.group[,2],3),
	knn1=round(pr.knn$prob.group[,1],3),knn2=round(pr.knn$prob.group[,2],3))
dfaux2


###################################################
### code chunk number 8: DD
###################################################
par(mfrow=c(1,3))
res.DD1=classif.DD(g,x,depth="mode",classif="qda")
res.DD2=classif.DD(g,x,depth="FM",classif="gam")
res.DD3=classif.DD(g,x,depth="mode",classif="np")
print(paste0("Est. DD1=",round(mean(res.DD1$prob.classification),3),", Est. DD2=",round(mean(res.DD2$prob.classification),3),", Est. DD3=",round(mean(res.DD3$prob.classification),3)))
pr.DD1=predict(res.DD1,xnew)
pr.DD2=predict(res.DD2,xnew)
pr.DD3=predict(res.DD3,xnew)
dfaux3=data.frame(DD1.g=pr.DD1,DD2.g=pr.DD2,DD3.g=pr.DD3)
dfaux3


###################################################
### code chunk number 9: boost
###################################################
classif.adaboost<-fda.usc:::classif.adaboost
predict.classif.adaboost<-fda.usc:::predict.classif.adaboost

res.boost=classif.adaboost(gr~X,data=ldata,classif="classif.glm",par.classif=list(basis.x=b.x))
#pr.boost=predict.classif.boosting(res.boost,newdata=ldatanew)
pr.boost=predict.classif.adaboost(res.boost,newdata=ldatanew)
print(paste0("Est. GLM=",round(mean(res.glm$prob.classification),3),", Est. Boost=",round(mean(res.boost$prob.classification),3)))
dfaux4=data.frame(glm.g=pr.glm$group.pred,Boost.g=pr.boost)
dfaux4

#fda.usc:::meas2accuracy(ldata$df$gr,pred1)
#fda.usc:::meas2accuracy(ldata$df$gr,pred2)

###################################################
### code chunk number 10: convex
###################################################
#setwd("C:/Users/moviedo/OneDrive - Universidade de Santiago de Compostela/GitHub/fda.usc.beta_2019_06_21/inst/script/potential")
#??fblocconvex
#source("remuestreoconvexo.R")
library(fda.usc)
bconvex<-fda.usc:::bconvex
blocconvex<-fda.usc:::blocconvex
fblocconvex<-fda.usc:::fblocconvex
bvar<-fda.usc:::bvar
bblocvar<-fda.usc:::bblocvar
blocvar<-fda.usc:::blocvar

traceback()
aa=cbind(rnorm(100),rnorm(100))
faa<-fdata(aa)
xx=bconvex(aa,J=4,B=5000)
xx2=blocconvex(aa,k=25,J=4,B=5000)
fxx2=fblocconvex(faa,k=25,J=4,B=5000)

xx3=bvar(aa,B=5000)
xx4=blocvar(aa,k=25,B=5000)
par(mfrow=c(1,2))

plot(rbind(xx2,aa),col=2,main="Local Convex",xlab="X1",ylab="X2")
points(aa,pch=19)
dd=kde2d(xx2[,1],xx2[,2],h=c(bw.SJ(xx2[,1]),bw.SJ(xx2[,2]))*3,n=51)
contour(dd,add=TRUE,col="blue")



plot(rbind(xx3,aa),col=2,main="Gaussian",xlab="X1",ylab="X2")
points(aa,pch=19)
ee=kde2d(xx3[,1],xx3[,2],h=c(bw.SJ(xx3[,1]),bw.SJ(xx3[,2]))*3,n=51)
contour(ee,add=TRUE,col="blue")

plot(rbind(xx4,aa),col=2,main="Local Gaussian",xlab="X1",ylab="X2")
points(aa,pch=19)
ff=kde2d(xx4[,1],xx4[,2],h=c(bw.SJ(xx3[,1]),bw.SJ(xx3[,2]))*3,n=51)
contour(ff,add=TRUE,col="blue")


###################################################
### code chunk number 11: tec
###################################################
classif.adaboost<-fda.usc:::classif.adaboost
predict.classif.adaboost<-fda.usc:::predict.classif.adaboost

data(tecator)
ab=tecator$absorp.fdata
ab2=fdata.deriv(ab,2)
ifat=factor(ifelse(tecator$y$Fat<8,1,ifelse(tecator$y$Fat>16,3,2)),label=c("F0","F8","F16"))
ltec=list(df=data.frame(ifat=ifat),ab=ab,ab2=ab2);class(ltec)="ldata"
btec=list(ab=create.pc.basis(ab,1:4),ab2=create.pc.basis(ab2,1:4))
print(tt<-table(ifat))
wei=1/tt[ifat]
wei=wei/sum(wei)
ww=(1/tt)/sum(1/tt)
cat("Frequencies:",tt,"\n")
cat("Weights per group:",round(ww,3),"\n")


# error    en el gsam (glm ok!)
ctrl2<- list(trace = FALSE, draw = TRUE,gray.scale=FALSE,fine=51)
colores=c("green","orange","red")

tec3 = classif.DD(ifat, ab, depth = "mode", classif = "gam",
                  control = ctrl2)

tec.glm1 = classif.glm(ifat ~ ab, ltec, basis.x = btec, type = "1vsall")
tec.glm2 = classif.glm(ifat ~ ab, ltec, basis.x = btec, type = "majority")
tec.glmw = classif.glm(ifat ~ ab, ltec, basis.x = btec, type = "majority",
                       weights = wei)
tec.ada = classif.adaboost(ifat ~ ab, ltec, classif = "classif.glm",
                           par.classif = list(basis.x = btec, type = "majority"))


ctrl2<- list(trace = FALSE, draw = TRUE,gray.scale=FALSE,fine=51)
colores=c("green","orange","red")

par(mfrow=c(1,2))
plot(ab,lty=1,col=colores[ifat])
legend("topleft",legend=c("F[0,8]","F[8,16]","F[16-]"),col=colores,lwd=2)
plot(ab2,lty=1,col=colores[ifat],main="2nd derivative")
legend("topleft",legend=c("F[0,8]","F[8,16]","F[16-]"),col=colores,lwd=2)

  
tec3=classif.DD(ifat,ab,depth="mode",classif="gam",control=ctrl2)
tec.glm1=classif.glm(ifat~ab,ltec,basis.x=btec,type="1vsall")
tec.glm2=classif.glm(ifat~ab,ltec,basis.x=btec,type="majority")
tec.glmw=classif.glm(ifat~ab,ltec,basis.x=btec,type="majority",weights=wei)
tec.ada=classif.adaboost(ifat~ab,ltec,classif="classif.glm",par.classif=list(basis.x=btec,type="majority"))
M=rbind(tec3$prob.classification,tec.glm1$prob.classification,tec.glm2$prob.classification,tec.glmw$prob.classification,tec.ada$prob.classification)%*%ww
M=cbind(rbind(tec3$prob.classification,tec.glm1$prob.classification,tec.glm2$prob.classification,tec.glmw$prob.classification,tec.ada$prob.classification),M)
rownames(M)=c("DD","GLM/OVA","GLM/OVO","GLM/OVO-W","Adaboost")
colnames(M)=c("F0","F8","F16","Weig. Prob. Class.")
print(round(M,4))
  
  
tec3=classif.DD(ifat,ab2,depth="mode",classif="gam",control=ctrl2)
tec.glm1=classif.glm(ifat~ab2,ltec,basis.x=btec,type="1vsall")
tec.glm2=classif.glm(ifat~ab2,ltec,basis.x=btec,type="majority")
tec.glmw=classif.glm(ifat~ab2,ltec,basis.x=btec,type="majority",weights=wei)
tec.ada=classif.adaboost(ifat~ab2,ltec,classif="classif.glm",par.classif=list(basis.x=btec,type="majority"))
M=rbind(tec3$prob.classification,tec.glm1$prob.classification,tec.glm2$prob.classification,tec.glmw$prob.classification,tec.ada$prob.classification)%*%ww
M=cbind(rbind(tec3$prob.classification,tec.glm1$prob.classification,tec.glm2$prob.classification,tec.glmw$prob.classification,tec.ada$prob.classification),M)
rownames(M)=c("DD","GLM/OVA","GLM/OVO","GLM/OVO-W","Adaboost")
colnames(M)=c("F0","F8","F16","Weig. Prob. Class.")
print(round(M,4))
 
# ArrowHead
# cargar 
#library(fda.classif)
data(ArrowHead,package="fda.classif")
#library(fdatasets)

data(ArrowHead)
ah=fdata(ArrowHead[,1:83],argvals=seq(2,248,by=3),rangeval=c(0,250),names=list(main="Arrow Heads",xlab="t",ylab="X(t)"))
gr=factor(ArrowHead[,84],labels=c("Avonlea","Clovis","Mix"))
ah1=fdata.deriv(ah,1,nbasis=21)
tt=table(gr)
wei=1/tt[gr]
wei=wei/sum(wei)
ww=1/tt/sum(1/tt)
ah.np=classif.np(gr,ah,par.S=list(w=wei))
ah.DD=classif.DD(gr,ah,depth="mode",classif="gam")
ah.np2=classif.np(gr,ah1,par.S=list(w=wei))

par(mfrow=c(1,2))
plot(ah,col=c("blue","red","green")[gr],lty=1)
legend("topright",levels(gr),lwd=2,col=c("blue","red","green"))
plot(ah1,col=c("blue","red","green")[gr],lty=1)
legend("topright",levels(gr),lwd=2,col=c("blue","red","green"))
cat("Frequencies:",tt,"\n")
cat("Weights per group",ww,"\n")
M=rbind(ah.np$prob.classification,ah.np2$prob.classification,ah.DD$prob.classification)%*%ww
M=cbind(rbind(ah.np$prob.classification,ah.np2$prob.classification,ah.DD$prob.classification),M)
rownames(M)=c("NP","NP.Der.","DD(mode,gam)")
colnames(M)=c("Avonlea","Clovis","Mix","Weig. Prob. Class.")
print(round(M,4))

# pruebas mias
# ldata<-list("df"=data.frame(gr),"ah1"=ah1)
# arr.glm1=classif.glm(gr~ah1,type="1vsall",data=ldata)
# lpc<-list("ah1"=create.pc.basis(ah1,1:8))
# arr.glm1=classif.glm(gr~ah1,type="majority",data=ldata,basis.x=lpc)
# arr.glm1

#BeetleFly
data(BeetleFly,package="fda.classif")
bfly=fdata(BeetleFly[,1:512],argvals=0:511,rangeval=c(0,511),names=list(main="Beetle/Fly",xlab="t",ylab="X(t)"))
gr=factor(BeetleFly[,513],labels=c("Beetle","Fly"))
bfly1=fdata.deriv(bfly,1,nbasis=21)
tt=table(gr)
wei=1/tt[gr]
wei=wei/sum(wei)
ww=1/tt/sum(1/tt)
bf.np=classif.np(gr,bfly,par.S=list(w=wei))
bf.DD=classif.DD(gr,bfly,depth="mode",classif="gam")
bf.np2=classif.np(gr,bfly,metric=metric.DTW,par.S=list(w=wei))

plot(bfly,col=c("blue","red")[gr],lty=1)
legend("topright",levels(gr),lwd=2,col=c("blue","red"))
#plot(bfly1,col=gr,lty=1)
#legend("topright",levels(gr),lwd=2,col=c("blue","red"))
cat("Frequencies:",tt,"\n")
cat("Weights per group",ww,"\n")
M=rbind(bf.np$prob.classification,bf.np2$prob.classification,bf.DD$prob.classification)%*%ww
M=cbind(rbind(bf.np$prob.classification,bf.np2$prob.classification,bf.DD$prob.classification),M)
rownames(M)=c("NP","NP/DTW","DD-mode")
colnames(M)=c("Beetle","Fly","Weig. Prob. Class.")
print(round(M,4))
# pruebas mias
ldata<-list("df"=data.frame(gr),"bfly1"=bfly1)
classif.glm(gr~bfly1,type="1vsall",data=ldata)
lpc<-list("bfly1"=create.pc.basis(bfly1,1:8))
classif.glm(gr~bfly1,type="majority",data=ldata,basis.x=lpc)
 
#ECG
data(ECG,package="fda.classif")
ecg=ECG$learn
gr=ECG$classlearn

tt=table(gr)
wei=1/tt[gr]
wei=wei/sum(wei)
ww=1/tt/sum(1/tt)
ecg.np=classif.np(gr,ecg,par.S=list(w=wei))
ecg.DD=classif.DD(gr,ecg,depth="mode",classif="gam")

plot(ecg,col=c("blue","red")[gr],lty=1)
legend("topright",levels(gr),lwd=2,col=c("blue","red"))
#plot(bfly1,col=gr,lty=1)
#legend("topright",levels(gr),lwd=2,col=c("blue","red"))
cat("Frequencies:",tt,"\n")
cat("Weights per group",ww,"\n")
M=rbind(ecg.np$prob.classification,ecg.DD$prob.classification)%*%ww
M=cbind(rbind(ecg.np$prob.classification,ecg.DD$prob.classification),M)
rownames(M)=c("NP","DD-mode")
colnames(M)=c(levels(gr),"Weig. Prob. Class.")
print(round(M,4))
pr=predict(ecg.np,ECG$test)
pr2=predict(ecg.DD,ECG$test)
M2=cbind(table(ECG$classtest,pr),table(ECG$classtest,pr2))
colnames(M2)=c(paste0("NP-",levels(gr)),paste0("DD-",levels(gr)))
cat("Test\n")
print(M2)

# GUN
data(gun,package="fda.classif")
fgun=gun$learn
gr=gun$classlearn
lgun=list(df=data.frame(gr=gr),fgun=fgun)
#class(lgun)="ldata"
bgun=list(fgun=create.pc.basis(fgun,l=1:5))
tt=table(gr)
wei=1/tt[gr]
wei=wei/sum(wei)
ww=1/tt/sum(1/tt)
gun.np=classif.np(gr,fgun,par.S=list(w=wei))
gun.glm=classif.gsam(gr~s(fgun),data=lgun,basis.x=bgun)
gun.gsam=classif.gsam(gr~s(fgun),data=lgun,basis.x=bgun,weights=wei)
gun.gkam=classif.gkam(gr~fgun,data=lgun,weights=wei)
plot(fgun,col=c("blue","red")[gr],lty=1)
legend("bottom",levels(gr),lwd=2,col=c("blue","red"))
cat("Frequencies:",tt,"\n")
cat("Weights per group",ww,"\n")
M=rbind(gun.np$prob.classification,gun.glm$prob.classification,gun.gsam$prob.classification,gun.gkam$prob.classification)%*%ww
M=cbind(rbind(gun.np$prob.classification,gun.glm$prob.classification,gun.gsam$prob.classification,gun.gkam$prob.classification),M)
rownames(M)=c("NP","GSAM/E.","GSAM/W.","GKAM/W")
colnames(M)=c(strtrim(levels(gr),3),"Weig. Prob. Class.")
print(round(M,4))
ltest=list(fgun=gun$test)
pr=predict(gun.np,gun$test)
pr2=predict(gun.glm,ltest)
pr3=predict(gun.gsam,ltest)
pr4=predict(gun.gkam,ltest)
M2=cbind(table(gun$classtest,pr),table(gun$classtest,pr2),table(gun$classtest,pr3),table(gun$classtest,pr4))
colnames(M2)=c(paste0("NP-",strtrim(levels(gr),3)),paste0("GAM/E-",strtrim(levels(gr),3)),paste0("GAM/W-",strtrim(levels(gr),3)),paste0("GKAM/W-",strtrim(levels(gr),3)))
rownames(M2)=strtrim(levels(gr),3)
cat("Test\n")
print(M2)

# WAFER
data(wafer,package="fda.classif")
waf=log(wafer$wafer+3)
waf1=fdata.deriv(waf,nbasis=31)
gr=wafer$class
lwaf=list(df=data.frame(gr=gr),waf=waf,waf1=waf1)
#class(lwaf)="ldata"
bwaf=list(waf=create.pc.basis(waf,l=1:8),waf1=create.pc.basis(waf1,l=1:8))
tt=table(gr)
wei=1/tt[gr]
wei=wei/sum(wei)
ww=1/tt/sum(1/tt)
waf.np=classif.np(gr,waf,par.S=list(w=wei))
waf.glm=classif.glm(gr~waf+waf1,data=lwaf,basis.x=bwaf,weights=wei)
waf.gsam=classif.gsam(gr~s(waf)+s(waf1),data=lwaf,basis.x=bwaf,weights=wei)
# Time consuming 
ctrl1<-list(maxit = 2) 

waf.gkam=classif.gkam(gr~waf+waf1,data=lwaf,weights=wei,control=ctrl1)
par(mfrow=c(1,2))
plot(waf,col=c("blue","red")[gr],lty=1)
legend("topright",levels(gr),lwd=2,col=c("blue","red"))
plot(waf1,col=c("blue","red")[gr],lty=1)
legend("topright",levels(gr),lwd=2,col=c("blue","red"))
cat("Frequencies:",tt,"\n")
cat("Weights per group",ww,"\n")
M=rbind(waf.np$prob.classification,waf.glm$prob.classification,waf.gsam$prob.classification,waf.gkam$prob.classification)%*%ww
M=cbind(rbind(waf.np$prob.classification,waf.glm$prob.classification,waf.gsam$prob.classification,waf.gkam$prob.classification),M)
rownames(M)=c("NP","GLM","GSAM","GKAM")
colnames(M)=c(strtrim(levels(gr),3),"Weig. Prob. Class.")
print(round(M,4))

# phoneme
predict.adaboost=function(obj,datanew){
    n=nrow(datanew[[1]])
    nlev=nlevels(obj$group)
    m=length(obj$list.fit)
    Mpr=matrix(0,nrow=n,ncol=nlev)
    colnames(Mpr)=levels(obj$group)
    for (i in 1:m){
      pr=predict(obj$list.fit[[i]],datanew)
      for ( j in 1:n){
        Mpr[j,as.numeric(pr[j])]=Mpr[j,as.numeric(pr[j])]+obj$alpha.boost[i]}
    }
    final=factor(levels(obj$group)[apply(Mpr,1,which.max)],levels=levels(obj$group))
    return(final)
  }



data(phoneme)
pho=phoneme$learn
gr=factor(phoneme$classlearn,labels=c("sh","iy","dcl","aa","ao"))
mc=c("black","green","cyan","red","blue")
lpho=list(df=data.frame(gr=gr),pho=pho)
class(lpho)="ldata"
bpho=list(pho=create.pc.basis(pho,1:5))
tt=table(gr)
wei=1/tt[gr]
wei=wei/sum(wei)
ww=1/tt/sum(1/tt)
pho.np=classif.np(gr,pho,par.S=list(w=wei))
pho.glm=classif.glm(gr~pho,data=lpho,weights=wei,basis.x=bpho,type="majority")
pho.ada=classif.adaboost(gr~pho,lpho,classif="classif.glm",par.classif=list(basis.x=bpho,type="majority"))

mc2=rgb(t(col2rgb(mc)/255),alpha=0.4)
plot(phoneme$learn,col=mc2[phoneme$classlearn],lwd=0.75,lty=2,ylim=c(0,30))
legend("topright",c("sh","iy","dcl","aa","ao"),fill=mc,cex=0.8)
lines(func.mean.formula(pho~gr,data=lpho),lwd=3,lty=1,col=mc)

cat("Frequencies:",tt,"\n")
cat("Weights per group",ww,"\n")
M=rbind(pho.np$prob.classification,pho.glm$prob.classification,pho.ada$prob.classification)%*%ww
M=cbind(rbind(pho.np$prob.classification,pho.glm$prob.classification,pho.ada$prob.classification),M)
rownames(M)=c("NP","GLM","Adaboost")
colnames(M)=c(levels(gr),"Weig. Prob. Class.")
print(round(M,4))
lphonew=list(pho=phoneme$test)
pr=predict(pho.np,phoneme$test)
pr2=predict(pho.glm,lphonew)
#pr3=predict.adaboost(pho.ada,lphonew)
pr3=predict.classif.adaboost(pho.ada,lphonew)
to1=table(phoneme$classtest,pr)
to2=table(phoneme$classtest,pr2)
to3=table(phoneme$classtest,pr3)
M2=cbind(diag(prop.table(to1,1)),diag(prop.table(to2,1)),diag(prop.table(to3,1)))
colnames(M2)=c("NP","GLM","AdaBoost")
rownames(M2)=levels(gr)
cat("Test\n")
print(M2)
