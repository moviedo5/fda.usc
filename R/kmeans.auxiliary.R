kmeans.assig.groups=function(out,draw=TRUE,...){
if (!is.null(out$lcenters))  lxm=out$lcenters
else  lxm=NULL
mdist=out$z.dist
nr=nrow(out$fdataobj)
nc=ncol(out$fdataobj)
xm=out$centers[["data"]]
ncl=nrow(xm)
grupo=rep(0,nr)
d=out$d
par(mfrow=c(1,2))
ngroups=nrow(d)-nrow(out$fdataobj[["data"]])
for (i in 1:nr){    grupo[i]=which.min(d[(nr+1):(nr+ngroups),i])      }
if (draw){
 if (nr==2){
  plot(out$fdataobj,main="Assigning groups")
  for (i in 1:ngroups){points(xm[i,1],xm[i,2],col=i+1,pch=8,cex=1.5)}
  }
 else{
  plot(out$fdataobj,col=grupo+1,main="Assigning groups",lwd=.3,lty=2)
  lines(out$centers,col=2:(length(grupo+1)),lwd=3,lty=1)     #new
 }
}
if (nc==2) { for (j in 1:nc){points(xm[j,1],xm[j,2],col=j+1,pch=7,cex=1.5)}}
return(list("centers"=out$centers,"cluster"=grupo))
}
##

# @rdname kmeans.fd
# @export
kmeans.centers.update=function(out,group,dfunc=func.trim.FM,draw=TRUE,par.dfunc=list(trim=0.05),...){
  if (class(out)!="kmeans.fd") stop("Error: incorrect input data")
  z=out$fdataobj[["data"]]
  tt=out$fdataobj[["argvals"]]
  rtt<-out$fdataobj[["rangeval"]]
  names=out$fdataobj[["names"]]
  mdist=out$z.dist
  xm=out$centers[["data"]]
  centers=out$centers
  nr=nrow(z)
  nc=ncol(z)
  grupo=group
  ngroups=length(unique(group))
  d=out$d
  ncl=nrow(xm)
  for (j in 1:ngroups){
    if (sum((grupo==j))>0) {
      dm=z[grupo==j,]
      ind=which(grupo==j)
      if (is.vector(dm) || nrow(dm)<3) {k=j}#revisar pq  k no hace nada!!
      else   {
        par.dfunc$fdataobj<-centers
        par.dfunc$fdataobj$data<-dm
        
        stat=do.call(dfunc,par.dfunc)
      }
      
      if (is.fdata(stat)) xm[j,]=stat[["data"]]
      else  xm[j,]=stat
    }
  }
  centers$data=xm
  rownames(centers$data)<-paste("center ",1:ngroups,sep="")
  if (draw){
    if (nr==2){
      plot(out$fdataobj,main="Center update")
      for (i in 1:ngroups){points(xm[i,1],xm[i,2],col=i+1,pch=8,cex=1.5)}}
    else{
      plot(out$fdataobj,col="grey",lty=grupo+1,lwd=0.15,cex=0.2,main="Update centers")
      lines(centers,col=2:(length(grupo+1)),lwd=3,lty=1)
    }}
  return(list("centers"=centers,"cluster"=grupo))
}


#' @rdname kmeans.fd
#' @export
kmeans.center.ini=function(fdataobj,ncl=2,metric=metric.lp,draw=TRUE,method="sample",iter=100,par.metric=NULL,...){
if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
if (is.null(par.metric)) par.metric=list()
par.metric$fdata1<-fdataobj
#mdist=metric(fdataobj,...)
mdist=do.call(metric,par.metric)
z<-fdataobj[["data"]]
tt<-fdataobj[["argvals"]]
rtt<-fdataobj[["rangeval"]]
names<-fdataobj[["names"]]
nr=nrow(fdataobj)
nc=ncol(fdataobj)
#else stop("The argument metric is not properly defined")
if (is.vector(ncl)) {
  len.ncl=length(ncl)
  if    (len.ncl==1) {
    ind=1
    ngroups=ncl
    if (method=="sample")  {
      #            lxm=sample(nr,ngroups,replace=FALSE) #max.iter=1
      max.combn<-choose(nr,ngroups)
      max.iter<-min(iter,max.combn)
      vec<-array(NA,dim=c(ngroups,max.iter))
      vec.d<-rep(NA,max.iter)
      for (i in 1:max.iter)  {
        vec[,i]<-sample(1:nr,ngroups,replace=FALSE)
        vec.d[i]<-sum(mdist[vec[,i],vec[,i]])
      }
      ind.max<-which.max(vec.d)
      lxm<-vec[,ind.max]
    }
    else if (method=="exact") {
      co<-combn(1:nr,ngroups)
      vec<-rep(NA,ncol(co))
      for (i in 1:ncol(co)) {vec[i]<-sum(mdist[co[,i],co[,i]])}
      max.vec<-which.max(vec)
      lxm<-co[,max.vec]
    }
    else stop("Center initialization method unknown")
    xm=z[lxm,]
  }
  else      stop("the argument ncl is expected the number of groups to detect")
}         else      stop("the argument ncl is expected the number of groups to detect")
d=rbind(mdist,mdist[lxm,])
centers=fdata(xm,tt,rtt,names)
if (draw){
  if (nr==2){
    plot(fdataobj)
    for (i in 1:ngroups){points(xm[i,1],xm[i,2],col=i+1,pch=8,cex=1.5)}}
  else{
    plot(fdataobj,col="grey",lty=2)
    lines(centers,col=2:(ngroups+1),lwd=3,lty=1)
    #  for (i in 1:ngroups){    points(tt,xm[i,],col=i+1,lwd=3)}
  }
}
out=list("centers"=centers,"lcenters"=lxm,"z.dist"=mdist,"fdataobj"=fdataobj)
class(out)="kmeans.fd"
return(invisible(out))
}


