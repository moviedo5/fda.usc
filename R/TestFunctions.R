#' @name fEqMoments.test
#' 
#' @title Tests for checking the equality of means and/or covariance between two populations under gaussianity. 
#' 
#' @description Two tests for the equality of means and covariances of two populations are provided.
#'  Both tests are constructed under gaussianity following Horvath & Kokoszka, 2012, Chapter 5.
#' 
#' @details \code{\link{fmean.test.fdata}} computes the test for equality of means. 
#' \code{\link{cov.test.fdata}} computes the test for equality of covariance operators.
#' Both tests have asymptotic distributions under the null related with chi-square distribution. Also, a 
#' parametric bootstrap procedure is implemented in both cases. 
#' 
#' @aliases fmean.test.fdata cov.test.fdata
#' @param X.fdata \code{fdata} object containing the curves from the first population.
#' @param Y.fdata \code{fdata} object containing the curves from the second population.
#' @param method c("X2","Boot"). "X2" includes the asymptotic distribution. "Boot" computes the bootstrap approximation.
#' @param npc The number of principal components employed. If \code{npc} is negative and 0<\code{abs(npc)}<1, the number of components 
#' are determined for explaining, at least, \code{abs(p)}\% of variability.
#' @param alpha Confidence level. By default =0.95.
#' @param B Number of bootstrap replicas when method="Boot".
#' @param draw By default, \code{FALSE}. Plots the density of the bootstrap replicas jointly with the statistic. 

#' @return Return a list with:
#' \itemize{
#' \item \code{stat}: Value of the statistic.
#' \item \code{pvalue}: P-values for the test.
#' \item \code{vcrit}: Critical cutoff for rejecting the null hypothesis.
#' \item \code{p}: Degrees of freedom for X2 statistic.
#' \item \code{B}: Number of bootstrap replicas.
#' }
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.febrero@@usc.es}
#' @seealso See Also as \code{\link{fanova.RPm}, \link{fanova.onefactor}}.
#' @references Inference for Functional Data with Applications. Horvath, L and Kokoszka, P. (2012). Springer. 
#' @keywords htest
#' @examples 
#' \dontrun{
#' tt=seq(0,1,len=51)
#' bet=0
#' mu1=fdata(10*tt*(1-tt)^(1+bet),tt)
#' mu2=fdata(10*tt^(1+bet)*(1-tt),tt) 
#' fsig=1
#' X=rproc2fdata(100,tt,mu1,sigma="vexponential",par.list=list(scale=0.2,theta=0.35))
#' Y=rproc2fdata(100,tt,mu2,sigma="vexponential",par.list=list(scale=0.2*fsig,theta=0.35))
#' fmean.test.fdata(X,Y,npc=-.98,draw=TRUE)
#' cov.test.fdata(X,Y,npc=5,draw=TRUE)
#' bet=0.1
#' mu1=fdata(10*tt*(1-tt)^(1+bet),tt)
#' mu2=fdata(10*tt^(1+bet)*(1-tt),tt) 
#' fsig=1.5
#' X=rproc2fdata(100,tt,mu1,sigma="vexponential",par.list=list(scale=0.2,theta=0.35))
#' Y=rproc2fdata(100,tt,mu2,sigma="vexponential",par.list=list(scale=0.2*fsig,theta=0.35))
#' fmean.test.fdata(X,Y,npc=-.98,draw=TRUE)
#' cov.test.fdata(X,Y,npc=5,draw=TRUE)
#' }
#' 
#' @rdname fEqMoments.test
#' @export fmean.test.fdata

fmean.test.fdata=function(X.fdata,Y.fdata,method=c("X2","Boot"),npc=5,alpha=0.95,B=1000,draw=FALSE){
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
    qboot=quantile(Uboot,prob=alpha)
  }
  if (ix2){
    if (npc<0 & abs(npc)<1) {
      pcall=fdata2pc(c(Xcen$Xcen,Ycen$Xcen),ncomp=min(ncol(X.fdata),nrow(X.fdata)))
      sumvar=cumsum(pcall$d^2)/sum(pcall$d^2)
      npc=min(which(sumvar>=abs(npc)))
    } else {
      pcall=fdata2pc(c(Xcen$Xcen,Ycen$Xcen),ncomp=npc)
    }
    ai=inprod.fdata(difm,pcall$rotation[1:npc])
    Unm1=(n*m/(n+m))*sum(ai[1:npc]^2)
    Unm2=(n*m/(n+m))*sum(ai[1:npc]^2/(pcall$d[1:npc]^2/(n+m)))
  }
  if (draw) {
    plot(density(Uboot),main=paste0("Density of ",B," replicas"))
    abline(v=Unm,col=2,lwd=2)
    abline(v=Unm1,col=3,lwd=2)
  }
  stat=c(Unm2,Unm1,Unm)
  pv=c(1-pchisq(Unm2,npc),sum(Unm<Uboot)/B)
  vcrit=c(qchisq(alpha,npc),qboot)
  names(stat)=c("X2","Unm(1)","Unm")
  names(pv)=c("X2","Boot")
  names(vcrit)=paste0(c("X2:","Boot:"),alpha)
  return(list(stat=stat,pvalue=pv,
              vcrit=vcrit,p=npc,B=B))
}


#' @rdname fEqMoments.test
#' @export cov.test.fdata
cov.test.fdata=function(X.fdata,Y.fdata,method=c("X2","Boot"),npc=5,alpha=0.95,B=1000,draw=FALSE){
  if ("X2" %in% method) ix2=TRUE else {ix2=FALSE;Unm1=NA;Unm2=NA}
  if ("Boot" %in% method) iboot=TRUE else {iboot=FALSE;Unm=NA;Uboot=NA;qboot=NA;draw=FALSE}
  n=nrow(X.fdata)
  m=nrow(Y.fdata)
  Xcen=fdata.cen(X.fdata)
  Ycen=fdata.cen(Y.fdata)
  
  Sig=var(rbind(Xcen$Xcen$data,Ycen$Xcen$data))
  
  if (ix2){
    if (npc<0 & abs(npc)<1) {
      pcall=fdata2pc(c(Xcen$Xcen,Ycen$Xcen),ncomp=min(ncol(X.fdata),nrow(X.fdata)))
      sumvar=cumsum(pcall$d^2)/sum(pcall$d^2)
      npc=min(which(sumvar>=abs(npc)))
    } else {
      pcall=fdata2pc(c(Xcen$Xcen,Ycen$Xcen),ncomp=npc)
    }
    propvar=sum(pcall$d[1:npc]^2)/sum(pcall$d^2)
    ax=inprod.fdata(Xcen$Xcen,pcall$rotation[1:npc])
    ay=inprod.fdata(Ycen$Xcen,pcall$rotation[1:npc])
    lambdax=t(ax)%*%ax/n
    lambday=t(ay)%*%ay/m
    dlambda=(n*diag(lambdax)+m*diag(lambday))/(n+m)
    ldiag=outer(dlambda,dlambda)
    Tn=sum((lambdax-lambday)^2/ldiag)*n*m/(2*(n+m))
    T1=(.5*(n*m)/(n+m))*sum((log(diag(lambdax))-log(diag(lambday)))^2)
    stat=c(Tn,T1)
    vcrit=c(qchisq(alpha,npc*(npc+1)/2),qchisq(alpha,npc))
    pv=c(1-pchisq(Tn,npc*(npc+1)/2),1-pchisq(T1,npc))
    
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
    qboot=quantile(Tboot,prob=alpha)
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
  names(vcrit)=paste0(c("X2(p(p+1)/2)","X2(p)","Boot"),":",alpha)
  return(list(stat=stat,pvalue=pv,
              vcrit=vcrit,p=c(npc*(npc+1)/2,npc),B=B))
}

#' @name fEqDistrib.test
#' 
#' @title Tests for checking the equality of distributions between two functional populations. 
#' 
#' @description Three tests for the equality of distributions of two populations are provided. The null hypothesis is that the two populations are the same
#' 
#' @details \code{\link{XYRP.test}} computes the p-values using random projections. Requires \code{kSamples} library. 
#' \code{\link{MMD.test}} computes Maximum Mean Discrepancy p-values using permutations (see Sejdinovic et al, (2013)) and \code{\link{MMDA.test}} 
#' does the same using an asymptotic approximation. 
#' \code{\link{fEqDistrib.test}} checks the equality of distributions using an embedding in a RKHS and two bootstrap approximations for 
#' calibration. 
#' 
#' @aliases XYRP.test MMD.test MMDA.test fEqDistrib.test
#' @param X.fdata \code{fdata} object containing the curves from the first population.
#' @param Y.fdata \code{fdata} object containing the curves from the second population.
#' @param nproj Number of projections for \code{XYRP.test}.
#' @param npc The number of principal components employed for generating the random projections.
#' @param test For \code{XYRP.test} "KS" and/or "AD" for computing Kolmogorov-Smirnov or Anderson-Darling p-values in the projections.
#' @param kern For \code{MMDA.test} "RBF" or "metric" for indicating the use of Radial Basis Function or directly, the distances.
#' @param metric Character with the metric function for computing distances among curves. 
#' @param ops.metric List of parameters to be used with \code{metric}.
#' @param method In \code{fEqDistrib.test} a character indicating the bootstrap method for computing the distribution under H0.
#'  "Exch" for Exchangeable bootstrap and "WildB" for Wild Bootstrap. By default, both are provided. 
#' @param B Number of bootstrap or Monte Carlo replicas.
#' @param alpha Confidence level for computing the threshold. By default =0.95.
#' @param iboot In \code{fEqDistrib.test} returns the bootstrap replicas.
#' @param draw By default, FALSE. Plots the density of the bootstrap replicas jointly with the statistic. 

#' @return A list with the following components by function:
#' \itemize{
#' \item \code{XYRP.test}, \code{FDR.pv}: p-value using FDR, \code{proj.pv}: Matrix of p-values obtained for projections.
#' \item \code{MMD.test}, \code{MMDA.test}: \code{stat}: Statistic, \code{p.value}: p-value, \code{thresh}: Threshold at level \code{alpha}.
#' \item \code{fEqDistrib.test}, \code{result}: \code{data.frame} with columns \code{Stat} and \code{p.value}, 
#' \code{Boot}: \code{data.frame} with bootstrap replicas if \code{iboot=TRUE}.
#' }
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.febrero@@usc.es}
#' @seealso  \code{\link{fmean.test.fdata}, \link{cov.test.fdata}}.
#' @references Sejdinovic, D., Sriperumbudur, B., Gretton, A., Fukumizu, K. \emph{Equivalence of distance-based and RKHS-based statistics in Hypothesis Testing} The Annals of Statistics, 2013. 
#' DOI \bold{10.1214/13-AOS1140}. 
#' @keywords htest
#' @examples 
#' \dontrun{
#' tt=seq(0,1,len=51)
#' bet=0
#' mu1=fdata(10*tt*(1-tt)^(1+bet),tt)
#' mu2=fdata(10*tt^(1+bet)*(1-tt),tt) 
#' fsig=1
#' X=rproc2fdata(100,tt,mu1,sigma="vexponential",par.list=list(scale=0.2,theta=0.35))
#' Y=rproc2fdata(100,tt,mu2,sigma="vexponential",par.list=list(scale=0.2*fsig,theta=0.35))
#' fmean.test.fdata(X,Y,npc=-.98,draw=TRUE)
#' cov.test.fdata(X,Y,npc=5,draw=TRUE)
#' bet=0.1
#' mu1=fdata(10*tt*(1-tt)^(1+bet),tt)
#' mu2=fdata(10*tt^(1+bet)*(1-tt),tt) 
#' fsig=1.5
#' X=rproc2fdata(100,tt,mu1,sigma="vexponential",par.list=list(scale=0.2,theta=0.35))
#' Y=rproc2fdata(100,tt,mu2,sigma="vexponential",par.list=list(scale=0.2*fsig,theta=0.35))
#' fmean.test.fdata(X,Y,npc=-.98,draw=TRUE)
#' cov.test.fdata(X,Y,npc=5,draw=TRUE)
#' XYRP.test(X,Y,nproj=15)
#' MMD.test(X,Y,B=1000)
#' fEqDistrib.test(X,Y,B=1000)
#' }
#' 
#' @rdname fEqDistrib.test
#' @export XYRP.test

XYRP.test=function(X.fdata,Y.fdata,nproj=10,npc=5,test=c("KS","AD")){
pvalues=matrix(NA,nrow=nproj,ncol=length(test))
colnames(pvalues)=test
rownames(pvalues)=paste0("h",1:nproj)
aa=create.pc.basis(c(X.fdata,Y.fdata),1:npc)
h=rcombfdata(nproj,aa$basis)
for (i in 1:nproj){
Xpr=inprod.fdata(h[i],X.fdata)
Ypr=inprod.fdata(h[i],Y.fdata)
for (k in 1:length(test)){
if (test[k]=="KS") {
pvalues[i,k]=ks.test(Xpr,Ypr,exact=TRUE)$p.value } else {
pvalues[i,k]=ad.test(Xpr,Ypr,Nsim=1000)$ad[1,3]
					}
			}
 	
}
FDR.pv=apply(pvalues,2,pvalue.FDR)
names(FDR.pv)=paste0("FDR-",test)
return(list(FDR.pv=FDR.pv,proj.pv=pvalues))
}

#' @rdname fEqDistrib.test
#' @export MMD.test

MMD.test=function(X.fdata,Y.fdata,metric="metric.lp",B=1000,alpha=.95,kern="RBF",ops.metric=list(lp=2),draw=FALSE){
n=nrow(X.fdata)
m=nrow(Y.fdata)
if (is.null(n) | is.null(m)) stop("One of the objects X.fdata, Y.fdata has no rows")
if (is.function(metric)) metric=deparse(substitute(metric))
DX=do.call(metric,c(list(X.fdata),ops.metric))
DY=do.call(metric,c(list(Y.fdata),ops.metric))
DXY=do.call(metric,c(list(X.fdata,Y.fdata),ops.metric))
D=rbind(cbind(DX,DXY),cbind(t(DXY),DY))
if (kern=="RBF"){
hb=quantile(D[D>0],prob=.25)
MKp=exp(-0.5*(D/hb)^2)   # RBF kernel
} else {
# metrica
vnorm=drop(do.call("norm.fdata",c(list(fdataobj=c(X.fdata,Y.fdata),metric=get(metric)),ops.metric)))
onorm=outer(vnorm,vnorm,"+")
MKp=0.5*(onorm-D)
}

MMD2b=mean(MKp[1:n,1:n])+mean(MKp[(n+1):(n+m),(n+1):(n+m)])-2*mean(MKp[(n+1):(n+m),1:n])
#Permutations
MMD2H0=numeric(B)
for (i in 1:B){
pp=sample(1:(n+m))
MMD2H0[i]=mean(MKp[pp[1:n],pp[1:n]])+mean(MKp[pp[(n+1):(n+m)],pp[(n+1):(n+m)]])-2*mean(MKp[pp[(n+1):(n+m)],pp[1:n]])
}

if (draw){
plot(density(MMD2H0),main="Density of MMD(H0) by Shuffling",xlim=range(MMD2b,MMD2H0))
abline(v=MMD2b,col=2)
}
pvnum2=mean(MMD2b<=MMD2H0)

result=list(stat=MMD2b,p.value=pvnum2,thresh=quantile(MMD2H0,alpha,na.rm=TRUE)) 
return(result)
}

#' @rdname fEqDistrib.test
#' @export MMDA.test

MMDA.test=function(X.fdata,Y.fdata,metric="metric.lp",B=1000,alpha=.95,kern="RBF",ops.metric=list(lp=2),draw=FALSE){
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
vcrit=2*apply(sweep(matrix(rnorm(B*ik)^2,nrow=B,ncol=ik),2,lambda[1:ik],"*"),1,sum)
pvnum=mean(vcrit>MMD2b)
if (draw){
plot(density(vcrit),main="Density of MMD(H0) by Asymptotic Approx.",xlim=range(MMD2b,vcrit))
abline(v=MMD2b,col=2)
}
#print(paste0(round(MMD2b,3),"-",round(quantile(vcrit,0.95),3)))
result=list(stat=MMD2b,p.value=pvnum,thresh=quantile(vcrit,alpha,na.rm=TRUE)) 
return(result)
}


#' @rdname fEqDistrib.test
#' @export fEqDistrib.test


fEqDistrib.test=function(X.fdata,Y.fdata,metric="metric.lp",method=c("Exch","WildB"),B=5000,ops.metric=list(lp=2),iboot=FALSE){
n=nrow(X.fdata)
m=nrow(Y.fdata)

if (is.null(n) | is.null(m)) stop("One of the objects X.fdata, Y.fdata has no rows")
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




