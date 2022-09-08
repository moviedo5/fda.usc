
mean.test.fdata=function(X.fdata,Y.fdata,method=c("X2","Boot"),p=5,B=1000,draw=FALSE){
  if ("X2" %in% method) ix2=TRUE else {ix2=FALSE;Unm1=NA;Unm2=NA}
  if ("Boot" %in% method) iboot=TRUE else {iboot=FALSE;Unm=NA;Uboot=NA;qboot=NA;draw=FALSE}
  n=nrow(X.fdata)
  m=nrow(Y.fdata)
  Xcen=fdata.cen(X.fdata)
  Ycen=fdata.cen(Y.fdata)
  mX=Xcen$meanX
  mY=Ycen$meanX
  difm=mX-mY
  Sig=var(rbind(Xcen$Xcen$data,Ycen$Xcen$data))
  if (iboot) {
    Unm=as.numeric((n*m/(n+m))*norm.fdata(difm)^2)
    Xb=rproc2fdata(B,t=argvals(mX),sigma=Sig)
    Uboot=drop(norm.fdata(Xb)^2)
    qboot=quantile(Uboot,prob=.95)
  }
  if (ix2){
    if (p<0 & abs(p)<1) {
      pcall=fdata2pc(c(Xcen$Xcen,Ycen$Xcen),ncomp=min(ncol(X.fdata),nrow(X.fdata)))
      sumvar=cumsum(pcall$d^2)/sum(pcall$d^2)
      p=min(which(sumvar>=abs(p)))
    } else {
      pcall=fdata2pc(c(Xcen$Xcen,Ycen$Xcen),ncomp=p)
    }
    ai=inprod.fdata(difm,pcall$rotation[1:p])
    Unm1=(n*m/(n+m))*sum(ai[1:p]^2)
    Unm2=(n*m/(n+m))*sum(ai[1:p]^2/(pcall$d[1:p]^2/(n+m)))
  }
  if (draw) {
    plot(density(Uboot),main=paste0("Density of ",B," replicas"))
    abline(v=Unm,col=2,lwd=2)
    abline(v=Unm1,col=3,lwd=2)
  }
  stat=c(Unm2,Unm1,Unm)
  pv=c(1-pchisq(Unm2,p),sum(Unm<Uboot)/B)
  vcrit=c(qchisq(.95,p),qboot)
  names(stat)=c("X2","Unm(1)","Unm")
  names(pv)=c("X2","Boot")
  names(vcrit)=c("X2","Boot")
  return(list(stat=stat,pvalue=pv,
              vcrit=vcrit,p=p,B=B))
}



cov.test.fdata=function(X.fdata,Y.fdata,method=c("X2","Boot"),p=5,B=1000,draw=FALSE){
  if ("X2" %in% method) ix2=TRUE else {ix2=FALSE;Unm1=NA;Unm2=NA}
  if ("Boot" %in% method) iboot=TRUE else {iboot=FALSE;Unm=NA;Uboot=NA;qboot=NA;draw=FALSE}
  n=nrow(X.fdata)
  m=nrow(Y.fdata)
  Xcen=fdata.cen(X.fdata)
  Ycen=fdata.cen(Y.fdata)
  
  Sig=var(rbind(Xcen$Xcen$data,Ycen$Xcen$data))
  
  if (ix2){
    if (p<0 & abs(p)<1) {
      pcall=fdata2pc(c(Xcen$Xcen,Ycen$Xcen),ncomp=min(ncol(X.fdata),nrow(X.fdata)))
      sumvar=cumsum(pcall$d^2)/sum(pcall$d^2)
      p=min(which(sumvar>=abs(p)))
    } else {
      pcall=fdata2pc(c(Xcen$Xcen,Ycen$Xcen),ncomp=p)
    }
    propvar=sum(pcall$d[1:p]^2)/sum(pcall$d^2)
    ax=inprod.fdata(Xcen$Xcen,pcall$rotation[1:p])
    ay=inprod.fdata(Ycen$Xcen,pcall$rotation[1:p])
    lambdax=t(ax)%*%ax/n
    lambday=t(ay)%*%ay/m
    dlambda=(n*diag(lambdax)+m*diag(lambday))/(n+m)
    ldiag=outer(dlambda,dlambda)
    Tn=sum((lambdax-lambday)^2/ldiag)*n*m/(2*(n+m))
    T1=(.5*(n*m)/(n+m))*sum((log(diag(lambdax))-log(diag(lambday)))^2)
    stat=c(Tn,T1)
    vcrit=c(qchisq(.95,p*(p+1)/2),qchisq(.95,p))
    pv=c(1-pchisq(Tn,p*(p+1)/2),1-pchisq(T1,p))
    
  }
  
  if (iboot) {
    Sigx=var(Xcen$Xcen$data)
    Sigy=var(Ycen$Xcen$data)
    statHS=sum((Sigx-Sigy)^2)
    Tboot=numeric(B)
    for (i in 1:B){
      Xb=rproc2fdata(n,t=argvals(X.fdata),sigma=Sig)
      Yb=rproc2fdata(m,t=argvals(Y.fdata),sigma=Sig)
      SbX=var(Xb$data)
      SbY=var(Yb$data)
      Tboot[i]=sum((SbX-SbY)^2)
    }
    qboot=quantile(Tboot,prob=.95)
  }
  if (draw) {
    plot(density(Tboot),main=paste0("Density of ",B," replicas"))
    abline(v=statHS,col=2,lwd=2)
    abline(v=sum((lambdax-lambday)^2),col=3,lwd=2)
  }
  stat=c(stat,statHS)
  pv=c(pv,sum(statHS<Tboot)/B)
  vcrit=c(vcrit,qboot)
  names(stat)=c("Tn","T1","HS")
  names(pv)=c("X2(p(p+1)/2)","X2(p)","Boot")
  names(vcrit)=c("X2(p(p+1)/2)","X2(p)","Boot")
  return(list(stat=stat,pvalue=pv,
              vcrit=vcrit,p=c(p*(p+1)/2,p),B=B))
}



XY.RP=function(X,Y,nproj=10,alpha=0.95,npc=5,test=c("KS","AD")){
#require(twosamples)
require(kSamples)
pvalues=matrix(NA,nrow=nproj,ncol=length(test))
colnames(pvalues)=test
rownames(pvalues)=paste0("h",1:nproj)
aa=create.pc.basis(c(X,Y),1:npc)
h=rcombfdata(nproj,aa$basis)
for (i in 1:nproj){
Xpr=inprod.fdata(h[i],X)
Ypr=inprod.fdata(h[i],Y)
for (k in 1:length(test)){
if (test[k]=="KS") {
pvalues[i,k]=ks.test(Xpr,Ypr,exact=TRUE)$p.value } else {
pvalues[i,k]=ad.test(Xpr,Ypr,Nsim=1000)$ad[1,3]
					}
			}	
}
return(pvalues)
}

MMD.test=function(X.fdata,Y.fdata,metric="metric.lp",nMC=1000,ops.metric=list(lp=2),draw=FALSE){
n=nrow(X.fdata)
m=nrow(Y.fdata)
if (is.null(n) | is.null(m)) stop("One of the objects X.fdata, Y.fdata has no rows")
if (is.function(metric)) metric=deparse(substitute(metric))
DX=do.call(metric,c(list(X.fdata),ops.metric))
DY=do.call(metric,c(list(Y.fdata),ops.metric))
DXY=do.call(metric,c(list(X.fdata,Y.fdata),ops.metric))
D=rbind(cbind(DX,DXY),cbind(t(DXY),DY))

#hb=quantile(D[D>0],prob=.25)
#MKp=exp(-0.5*(D/hb)^2)   # RBF kernel

vnorm=drop(do.call("norm.fdata",c(list(fdataobj=c(X.fdata,Y.fdata),metric=get(metric)),ops.metric)))
onorm=outer(vnorm,vnorm,"+")
#D=do.call(metric,c(list(fdata1=c(X.fdata,Y.fdata)),ops.metric))
MKp=0.5*(onorm-D)

#H=diag(n+m)-matrix(1/(n+m),nrow=nrow(MKp),ncol=ncol(MKp))
#Ktilde=H%*%MKp%*%H
#MKp=Ktilde

#lambda=eigen(Ktilde)$values
#k=min(which(cumsum(lambda^2)/sum(lambda^2)>0.98))
#vcrit=apply(sweep(matrix(rnorm(nMC*k)^2,nrow=nMC,ncol=k),2,lambda[1:k],"*"),1,sum)
MMD2b=mean(MKp[1:n,1:n])+mean(MKp[(n+1):(n+m),(n+1):(n+m)])-2*mean(MKp[(n+1):(n+m),1:n])
#Permutations
MMD2H0=numeric(nMC)
for (i in 1:nMC){
pp=sample(1:(n+m))
MMD2H0[i]=mean(MKp[pp[1:n],pp[1:n]])+mean(MKp[pp[(n+1):(n+m)],pp[(n+1):(n+m)]])-2*mean(MKp[pp[(n+1):(n+m)],pp[1:n]])
}
#K=max(MKp)
#mp=2*n*m/(n+m)

#par(mfrow=c(1,2))
#plot(density(vcrit*2/mp))
#abline(v=MMD2b)
if (draw){
plot(density(MMD2H0),main="Density of MMD(H0) by Shuffling")
abline(v=MMD2b)
}
#pvasym=1/exp((sqrt(MMD2b/(2*K/mp))-1)^2/2)
#print(paste0("N.PC's:",k," Thr.As.:",round((2*K/mp)*(1+sqrt(2*log(1/0.05)))^2,4), " Thr.MC:",round(quantile(2*vcrit/mp,0.95),3), 
#	" Thr.Per:",round(quantile(MMD2H0,0.95),3)))
#pvnum=mean(MMD2b>2*vcrit/mp)
pvnum2=mean(MMD2b<=MMD2H0)

#result=data.frame(Stat=c(sqrt(MMD2b),mp*MMD2b/2,MMD2b),p.values=c(pvasym,pvnum,pvnum2))
result=list(stat=MMD2b,p.value=pvnum2,thresh=quantile(MMD2H0,.95)) 
#rownames(result)=c("Asymp","MC","Perm")
return(result)
}


fEqDistrib.test=function(X.fdata,Y.fdata,metric="metric.lp",method=c("Exch","WildB"),nboot=5000,ops.metric=list(lp=2),iboot=FALSE){
n=nrow(X.fdata)
m=nrow(Y.fdata)

if (is.null(n) | is.null(m)) stop("One of the objects X.fdata, Y.fdata has no rows")
B=nboot
DX=do.call(metric,c(list(X.fdata),ops.metric))
DY=do.call(metric,c(list(Y.fdata),ops.metric))
DXY=do.call(metric,c(list(X.fdata,Y.fdata),ops.metric))
Tn=n*(sum(DXY)*2/(n*m)-sum(DX)/n^2-sum(DY)/m^2)
VX=(1/(2*n^2))*sum(DX)
VY=(1/(2*m^2))*sum(DY)
TnS=Tn/(VX+(n/m)*VY)


iexch=FALSE;iwildB=FALSE
if ("Exch" %in% method){iexch=TRUE}
if ("WildB" %in% method) {iwildB=TRUE}
if (!iexch & !iwildB) iexch=TRUE
if (iexch){
BootETN=numeric(B)
BootESTN=numeric(B)
}
if (iwildB){
BootWTN=numeric(B)
BootWSTN=numeric(B)
}

for (i in 1:B){
if (iexch) {
W=runif(n)
V=runif(m)
W=n*W/sum(W)
V=m*V/sum(V)
w=(W-1)/n
v=(V-1)/m
cW=sqrt(mean((W-1)^2))
cV=sqrt(mean((V-1)^2))
VEX=(1/2)*sum(outer(w+1/n,w+1/n)*DX)
VEY=(1/2)*sum(outer(v+1/m,v+1/m)*DY)
#BootETN[i]=n*(sum(outer(w,v)*DXY)*2-sum(outer(w,w)*DX)-sum(outer(v,v)*DY))/(cV*cW)		#mean((c(W,V)-1)^2)
BootETN[i]=n*(sum(outer(w,v)*DXY)*2/(cV*cW)-sum(outer(w,w)*DX)/(cW^2)-sum(outer(v,v)*DY)/(cV^2))		#mean((c(W,V)-1)^2)
BootESTN[i]=BootETN[i]/(VEX+(n/m)*VEY)
}
if (iwildB){
psi=rnorm(n)
eta=rnorm(m)
psie=(psi-mean(psi))/n
etae=(eta-mean(eta))/m
VWX=(1/2)*sum(outer((psie+1/n),(psie+1/n))*DX)
VWY=(1/2)*sum(outer((etae+1/m),(etae+1/m))*DY)
BootWTN[i]=n*(2*sum(outer(psie,etae)*DXY)-sum(outer(psie,psie)*DX)-sum(outer(etae,etae)*DY))
BootWSTN[i]=BootWTN[i]/(VWX+(n/m)*VWY)
}
}
if (sum(iexch,iwildB)==2){
#	result=vector("list",4)
	result=data.frame(Stat=rep(c(Tn,TnS),2),p.value=c(mean(BootETN>Tn),mean(BootESTN>TnS),mean(BootWTN>Tn),mean(BootWSTN>TnS)))
	rownames(result)=c(as.vector(outer(c("","Stud. "),c("Exch.","Wild.Boot"), paste0)))
} else if (iexch){
#	result=vector("list",2)
	result=data.frame(Stat=c(Tn,TnS),p.value=c(mean(BootETN>Tn),mean(BootESTN>TnS)))
	rownames(result)=c(paste0(c("","Stud. "),c("Exchangeable")))
	} else {
#	result=vector("list",2)
	result=data.frame(Stat=c(Tn,TnS),p.value=c(mean(BootWTN>Tn),mean(BootWSTN>TnS)))
	rownames(result)=c(paste0(c("","Stud. "),c("Wild Bootstrap")))
	}
#	result=list(statistic=statistic,p.value=p.value,method=method)
#	class(result)="htest"
	if (iboot) {
	return(list(result=result,Boot=data.frame(BExch=BootETN,BSExch=BootESTN,BWB=BootWTN,BSWB=BootWSTN)))} else {
	return(result)
	}
}

MMDA.test=function(X.fdata,Y.fdata,kern="RBF",metric="metric.lp",nMC=1000,ops.metric=list(lp=2),draw=FALSE){
#kern="RBF" Kernel o "metric"
n=nrow(X.fdata)
m=nrow(Y.fdata)
if (is.null(n) | is.null(m)) stop("One of the objects X.fdata, Y.fdata has no rows")
if (is.function(metric)) metric=deparse(substitute(metric))
DX=do.call(metric,c(list(X.fdata),ops.metric))
DY=do.call(metric,c(list(Y.fdata),ops.metric))
DXY=do.call(metric,c(list(X.fdata,Y.fdata),ops.metric))
D=rbind(cbind(DX,DXY),cbind(t(DXY),DY))
#D=do.call(metric,c(list(fdata1=c(X.fdata,Y.fdata)),ops.metric))
if (kern=="RBF"){
hb=quantile(D[D>0],prob=.25)
MKp=exp(-0.5*(D/hb)^2)   # RBF kernel
} else {
# metrica
vnorm=drop(do.call("norm.fdata",c(list(fdataobj=c(X.fdata,Y.fdata),metric=get(metric)),ops.metric)))
onorm=outer(vnorm,vnorm,"+")
MKp=0.5*(onorm-D)
}

H=diag(n+m)-matrix(1/(n+m),nrow=nrow(MKp),ncol=ncol(MKp))
Ktilde=H%*%MKp%*%H
MMD2b=n*(mean(MKp[1:n,1:n])+mean(MKp[(n+1):(n+m),(n+1):(n+m)])-2*mean(MKp[(n+1):(n+m),1:n]))

lambda=eigen(Ktilde)$values
ik=min(which(cumsum(lambda)/sum(lambda)>0.98))
lambda=lambda/(n+m)
vcrit=2*apply(sweep(matrix(rnorm(nMC*ik)^2,nrow=nMC,ncol=ik),2,lambda[1:ik],"*"),1,sum)
pvnum=mean(vcrit>MMD2b)
#print(paste0(round(MMD2b,3),"-",round(quantile(vcrit,0.95),3)))
result=list(stat=MMD2b,p.value=pvnum,thresh=quantile(vcrit,0.95)) 
return(result)
}



# fEqDistrib2.test=function(X.fdata,Y.fdata,metric="metric.lp",method=c("Exch","WildB"),B=5000,ops.metric=list(lp=2)){
# n=nrow(X.fdata)
# m=nrow(Y.fdata)

# if (is.null(n) | is.null(m)) stop("One of the objects X.fdata, Y.fdata has no rows")

# DX=do.call(metric,c(list(X.fdata),ops.metric))
# DY=do.call(metric,c(list(Y.fdata),ops.metric))
# DXY=do.call(metric,c(list(X.fdata,Y.fdata),ops.metric))
# Tn=n*(sum(DXY)*2/(n*m)-sum(DX)/n^2-sum(DY)/m^2)

# iexch=FALSE;iwildB=FALSE
# if ("Exch" %in% method){iexch=TRUE}
# if ("WildB" %in% method) {iwildB=TRUE}
# if (!iexch & !iwildB) iexch=TRUE
# if (iexch){
# BootETN=numeric(B)
# }
# if (iwildB){
# BootWTN=numeric(B)
# }

# for (i in 1:B){
# if (iexch) {
# W=runif(n)
# V=runif(m)
# W=n*W/sum(W)
# V=m*V/sum(V)
# w=(W-1)/n
# v=(V-1)/m
# cW=sqrt(mean((W-1)^2))
# cV=sqrt(mean((V-1)^2))
# BootETN[i]=n*(sum(outer(w,v)*DXY)*2/(cV*cW)-sum(outer(w,w)*DX)/(cW^2)-sum(outer(v,v)*DY)/(cV^2))		#mean((c(W,V)-1)^2)
# }
# if (iwildB){
# psi=rnorm(n)
# eta=rnorm(m)
# psie=(psi-mean(psi))/n
# etae=(eta-mean(eta))/m
# BootWTN[i]=n*(2*sum(outer(psie,etae)*DXY)-sum(outer(psie,psie)*DX)-sum(outer(etae,etae)*DY))
# }
# }

# if (sum(iexch,iwildB)==2){
	# result=list(statistic=Tn,p.value=c(mean(BootETN>Tn),mean(BootWTN>Tn)))
	# names(result$p.value)=c("Exch.Boot.","Wild Boot.")
# } else if (iexch){
	# result=list(statistic=Tn,p.value=mean(BootETN>Tn))
	# names(result$p.value)="Exchangeable"
	# } else {
# #	result=vector("list",2)
	# result=list(statistic=Tn,p.value=mean(BootWTN>Tn))
	# names(result$p.value)="Wild Bootstrap"
	# }
	# return(result)
# }




