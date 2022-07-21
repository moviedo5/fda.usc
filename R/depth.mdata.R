#' @name depth.mdata
#' 
#' @title Provides the depth measure for multivariate data
#' 
#' @description Compute measure of centrality of the multivariate data. Type of depth
#' function: simplicial depth (SD), Mahalanobis depth (MhD), Random Half--Space
#' depth (HS), random projection depth (RP) and Likelihood Depth (LD).
#' 
#' @details Type of depth measures:
#' \itemize{ 
#' \item The \code{\link{mdepth.SD}} calculates the simplicial depth (HD) of the points in \code{x} w.r.t.
#' \code{xx} (for bivariate data).  
#' \item The \code{\link{mdepth.HS}} function calculates the random half--space depth (HS)
#'  of the points in \code{x} w.r.t. \code{xx} based on random projections \code{proj}.  
#' \item The \code{\link{mdepth.MhD}} function calculates the Mahalanobis depth (MhD)
#'  of the points in \code{x} w.r.t. \code{xx}.  
#' \item The \code{\link{mdepth.RP}} calculates the random' projection depth (RP) 
#'  of the points in \code{x} w.r.t. \code{xx} based on random projections \code{proj}.
#' \item The \code{\link{mdepth.LD}} calculates the Likelihood depth (LD) of the points 
#' in \code{x} w.r.t. \code{xx}.  
#' \item The \code{\link{mdepth.TD}} function provides the Tukey depth measure for multivariate data.
# \item The \code{\link{mdepth.MB}} calculates the Modified Band
#depth (MB) of the points in \code{x} w.r.t. \code{xx}.  }
#' }
#' @aliases Depth.Multivariate mdepth.FM mdepth.HS mdepth.MhD mdepth.SD mdepth.LD
#' mdepth.TD mdepth.RP mdepth.FSD mdepth.KFSD
#' @param x is a set of points, a d-column matrix.
#' @param xx is a d-dimension multivariate reference sample (a d-column matrix)
#' where \code{x} points are evaluated.
#' @param proj are the directions for random projections, by default 500 random
#' projections generated from a scaled \code{runif(500,-1,1)}.
#' @param scale =TRUE, scale the depth, see \link[base]{scale}.
#' @param metric Metric function, by default \code{\link{metric.dist}}.
#' Distance matrix between \code{x} and \code{xx} is computed.
#' @param xeps Accuracy. The left limit of the empirical distribution function.
#' @param random =TRUE for random projections. =FALSE for deterministic
#' projections.
#' @param h Bandwidth, \code{h>0}. Default argument values are provided as the
#' 15\%--quantile of the distance between \code{x} and \code{xx}.
#' @param trim The alpha of the trimming.
#' @param draw =TRUE, draw the curves, the sample median and trimmed mean.
#' @param dfunc type of univariate depth function used inside depth function:
#' "FM1" refers to the original Fraiman and Muniz univariate depth (default),
#' "TD1" Tukey (Halfspace),"Liu1" for simplical depth, "LD1" for Likelihood
#' depth and "MhD1" for Mahalanobis 1D depth. Also, any user function
#' fulfilling the following pattern \code{FUN.USER(x,xx,...)} and returning a
#' \code{dep} component can be included.
#' @param \dots Further arguments passed to or from other methods.
#' @return 
#' \itemize{
#' \item {lmed}{ Index of deepest element \code{median} of \code{xx}.}
#' \item {ltrim}{ Index of set of points \code{x} with trimmed mean
#' \code{mtrim}. } 
#' \item {dep}{ Depth of each point \code{x} w.r.t. \code{xx}.}
#' \item {proj}{ The projection value of each point on set of points. }
#' \item {x}{is a set of points to be evaluated.}
#' \item {xx}{ a reference sample}
#' \item {name}{ Name of depth method }
#' }
#' @author \code{\link{mdepth.RP}}, \code{\link{mdepth.MhD}} and
#' \code{\link{mdepth.HS}} are versions created by Manuel Febrero Bande and
#' Manuel Oviedo de la Fuente of the original version created by Jun Li, Juan
#' A. Cuesta Albertos and Regina Y. Liu for polynomial classifier.
#' 
#' @seealso Functional depth functions: \code{\link{depth.FM}},
#' \code{\link{depth.mode}}, \code{\link{depth.RP}}, \code{\link{depth.RPD}}
#' and \code{\link{depth.RT}}.
#' 
#' @references Liu, R. Y., Parelius, J. M., and Singh, K. (1999). Multivariate
#' analysis by data depth: descriptive statistics, graphics and inference,(with
#' discussion and a rejoinder by Liu and Singh). \emph{The Annals of
#' Statistics}, 27(3), 783-858.
#' 
#' @keywords descriptive
#' @examples
#' \dontrun{
#' data(iris)
#' group<-iris[,5]
#' x<-iris[,1:2]
#'                                   
#' MhD<-mdepth.MhD(x)
#' PD<-mdepth.RP(x)
#' HD<-mdepth.HS(x)
#' SD<-mdepth.SD(x)
#' 
#' x.setosa<-x[group=="setosa",]
#' x.versicolor<-x[group=="versicolor",] 
#' x.virginica<-x[group=="virginica",]
#' d1<-mdepth.SD(x,x.setosa)$dep
#' d2<-mdepth.SD(x,x.versicolor)$dep
#' d3<-mdepth.SD(x,x.virginica)$dep
#' }
#' @rdname depth.mdata
#' @export 
mdepth.LD=function(x,xx=x,metric=metric.dist,
                   h=NULL,scale=FALSE,...){    
  if (is.vector(x)){
    if (all(xx==x)) x<-xx<-matrix(x,ncol=1)#stop("One of x or xx must be a matrix")
    else   {
      m2=ncol(xx)
      if (length(x)!=m2) stop("Length of x does not match with dimension of xx")
      x=matrix(x,ncol=m2)
    }
  }    
  m2=ncol(xx)
  n2=nrow(xx)
  n<-nrow(x)
  m<-ncol(x)
  if (is.null(rownames(x)))  rownames(x)<-1:nrow(x)
  nms<-rownames(x)
  x<-na.omit(x)  
  xx<-na.omit(xx)
  nas<-na.action(x)
  nullans<-!is.null(nas) 
  d <- ncol(x)
  if (is.null(n) && is.null(m)) stop("ERROR IN THE DATA DIMENSIONS")
  if (is.null(m) && is.null(m2)) stop("ERROR IN THE DATA DIMENSIONS")
  if (is.matrix(metric)) {mdist=metric}
  else {  mdist=metric(xx,xx,...)  }
  class(mdist)<-"matrix"
  
  if (n==n2 & m==m2) {
    equal<-all(x==xx)
    if (equal) mdist2<-mdist
    else mdist2<-metric(x,xx,...)}
  else  mdist2<-metric(x,xx,...)
  if (is.null(h))   {
    h<-0.15
    hq2=quantile(mdist,probs=h,na.rm=TRUE)
  }
  else {
    if (is.numeric(h))    hq2<-h  
    else hq2=quantile(mdist,probs=as.numeric(h),na.rm=TRUE)  
  }
  class(mdist)<-class(mdist2)<-c("matrix")
  ans<-Ker.norm(mdist2/hq2)    ####
  ans<-apply(ans,1,sum,na.rm=TRUE)                                    
  mx = scale
  if (scale)   {
    #  mn<-min(ans,na.rm=TRUE)
    ans2=Ker.norm(mdist/hq2)
    ans2=apply(ans2,1,sum,na.rm=TRUE)
    mx<-max(ans2,na.rm=TRUE)
    ans=as.vector(ans/mx)   
  }                                
  if  (nullans){
    ans1<-rep(NA,len=n)
    ans1[-nas] <-ans 
    ans<-ans1      
  }
  names(ans)<-nms   
  out<-list("dep" = ans,"hq"=hq2,dscale=mx,
            x=x,xx=xx,name="LD")
  class(out)<-"mdepth"
  return(invisible(out))
}


#' @rdname depth.mdata
#' @export 
mdepth.HS <-function(x, xx=x,proj=50,scale=FALSE,xeps=1e-15,random=FALSE)
{
#mdepth.HS:  calculates the half-space depth (HS) of the points in x w.r.t. xx based on projections 
#xx is a d-dimension multivariate sample, a d-column matrix
#x is a set of points, a d-column matrix, x can be missing
#proj are the directions for projections    
      if (is.vector(x)) {
         x<-matrix(x,nrow=1)
#         ans<-pmin(sum(x<=xx)/m,sum(x>=xx)/m)
         Fn=ecdf(xx)
         ans=pmin(Fn(x),(1-Fn(x-xeps)))
         if (scale) ans<-ans*2   
         return(invisible(list(dep = ans, Fn=Fn)))   
        }
       if ( is.null(rownames(x)))  rownames(x)<-1:nrow(x)
       nms<-rownames(x)
       m0<-nrow(x)
       xx<-na.omit(xx)
       x<-na.omit(x)
       nas<-na.action(x)
       nullans<-!is.null(nas)        
        n <- nrow(x)
        d <- dim(xx)[2]
        lenn<-length(proj)
        if (lenn==1) {
          mm<-proj[1]
          if (d==2 & !random) {
              sq<-seq(0,2*pi,len=mm)
#              proj<-data.matrix(expand.grid(cos(sq),sin(sq)))             
               proj2d<-function(angl){matrix(c(cos(angl),sin(angl)),2)}
               proj<-t(sapply(sq,proj2d   ))
             }
          else {
          if (d==3 & !random) {
			      mmr=floor(sqrt(mm))+1
      		  phi=seq(0,2*pi,len=mmr)
			      theta=seq(0,pi,len=mmr)
            exgrid=expand.grid(phi=phi,theta=theta)
			      proj=cbind(sin(exgrid$theta)*cos(exgrid$phi),sin(exgrid$theta)*sin(exgrid$phi),cos(exgrid$theta))
            mm<-nrow(proj)          
            }
          else{             
          warning("Method based on Random Projections")
          u <- matrix(runif(d*mm,-1,1),mm,d)
          norm <- sqrt(rowSums(u*u))
          proj <- (u/norm) 
         }
        }            
        }
        else   mm<-nrow(proj)
      out <- matrix(0, mm,n) 
      if (d==3 & !random) {
               for(i in 1:mm) {
#                 z<-t(a[,,i]%*%t(x))# este calculo esta malm  debe dar nproj*ndatos
#                 z2<-t(a[,,i]%*%t(xx))
        z  <- proj %*% t(x)            
        z2 <- proj %*% t(xx)
                  
                 out[i,] <- sapply(z[i,], function(y) min(sum(y<=z2[i,])/n,sum(y>=z2[i,])/n))
               }       
      }
      else{     
        z  <- proj %*% t(x)            
        z2 <- proj %*% t(xx)
        for(i in 1:mm) {
          out[i,] <- sapply(z[i,], function(y) min(sum(y<=z2[i,])/n,sum(y>=z2[i,])/n))
        }        
        }
    ans = as.vector(apply(out,2,min))        
   if (scale)          ans<-ans*2
  if  (nullans){
        ans1<-rep(NA,len=m0)
        ans1[-nas] <-ans 
        ans<-ans1      
        }
   names(ans)<-nms   
   out <- list(dep = ans, proj = proj, x=x,xx=xx,name="HS")
   return(invisible(out))
}


#' @rdname depth.mdata
#' @export   
mdepth.RP<-function(x, xx=x,proj=50,scale=FALSE){
  #################################################################################################################
  # depth.RD: calculates the projection depth (PD) of the points in x w.r.t. xx based on random projections proj
  #xx is a d-dimension multivariate sample, a d-column matrix
  #x is a set of points, a d-column matrix, x can be missing
  #proj are the directions fo random projections
  #trim the alpha of the trimming
  #draw=TRUE, draw the points in a gray scale of its depth, the sample median (in red) and trimmed mean (in yellow)
  #################################################################################################################
  if (is.vector(x)){
    if (all(xx==x)) x<-xx<-matrix(x,ncol=1)#stop("One of x or xx must be a matrix")
    else   {
      m2=ncol(xx)
      if (length(x)!=m2) stop("Length of x does not match with dimension of xx")
      x=matrix(x,ncol=m2)
    }
  }
  if ( is.null(rownames(x)))  rownames(x)<-1:nrow(x)
  nms<-rownames(x)
  m0<-nrow(x)
  xx<-na.omit(xx)
  x<-na.omit(x)
  nas<-na.action(x)
  nullans<-!is.null(nas) 
  n <- nrow(x)
  d <- ncol(x)
  lenn<-length(proj)
  if (lenn==1) {
    u <- matrix(runif(d*proj,-1,1),proj,d)
    norm <- sqrt(rowSums(u*u))
    proj <- u/norm
  }
  z <- proj %*% t(xx)
  z1 <- proj %*% t(x)        
  mm<-nrow(proj)
  pdep=matrix(NA,nrow=n,ncol=mm)
  for (i in 1:mm){
    pdep[,i]=mdepth.TD1(z1[i,],z[i,],scale=scale)$dep
    #        m1 <- m2 <- rep(0, mm)
    #        for(i in 1:mm) {
    #				 m1[i]= min(z[i,])
    #				 m2[i]= max(z[i,])
    #                m1[i] <- median(z[i,  ])
    #                m2[i] <- median(abs(z[i,  ] - m1[i]))
  }
  out1 <- apply(pdep,1,mean)  
  #        for(j in 1:n) {  out1[j] <- max(abs(z1[, j] - m1)/m2)       }             
  #        for(j in 1:n) {  out1[j] <- mean(abs(z1[, j] - m1)/(m2-m1))       }             
  ans=out1
  #        ans = 1/(1 + out1)        
  #        if (scale){ 
  #          ans<-ans/max(ans) #*2 =/.5
  #        }   
  if  (nullans){
    ans1<-rep(NA,len=m0)
    ans1[-nas] <-ans 
    ans<-ans1
  }
  names(ans)<-nms    
  out <- list( dep = ans,  proj = proj,x=x,xx=xx,name="RP")
  return(invisible(out))
}






















#' @rdname depth.mdata
#' @export 
mdepth.MhD <- function(x,xx=x,scale=FALSE){
#mdepth.MhD:  calculates the Mahalanobis depth (MhD) of the points in x w.r.t. xx
#xx is a d-dimension multivariate sample, a d-column matrix
#x is a set of points, a d-column matrix, x can be missing
#trim the alpha of the trimming
#draw=TRUE, draw the points in a gray scale of its depth, the sample median (in red) and trimmed mean (in yellow)
 m0 <- nrow(x)
 if (!is.vector(x) & is.null(rownames(x)))  rownames(x)<-1:nrow(x)
 nms<-rownames(x)
 x<-na.omit(x)    
 xx<-na.omit(xx)
 nas<-na.action(x)
 nullans<-!is.null(nas) 
 if (is.vector(x)) {
              D<- (1+(x-(mean(xx)))^2/sd(xx)^2)
              }
 else{

  n <- nrow(xx)
	m <- nrow(x)
	d<-ncol(x)
	mu <-colMeans(xx)
	sigma <- cov(xx)
	D <-rep(0,m)
  sigma.inv <- try(solve(sigma),silent=TRUE)#new
  
  if (!is.matrix(sigma.inv)) {
     sv<-svd(sigma)    
     sigma.inv<-sv$v%*%diag(1/sv$d)%*%t(sv$u)
     warning("Inverse of sigma computed by SVD")
    }
   D <- 1+apply(t(x)-mu,2, function(x) t(x)%*%sigma.inv%*%x)
   }
ans<-1/D      
if  (nullans){
        ans1<-rep(NA,len=m0)
        ans1[-nas] <-1/D 
        ans<-ans1
        }
  names(ans)<-nms   
#   if (scale) {        ans <- ans/max(ans)    }
  out <- list( dep = ans,  x=x,xx=xx,name="MhD")
   return(invisible(out))
}  

#' @rdname depth.mdata
#' @export 
mdepth.KFSD=function (x, xx = x, trim = 0.25,
                      h=NULL,scale = FALSE, draw = FALSE){
  if (is.matrix(x) & is.matrix(xx)){
    m0=nrow(x)
    rownames(x) <- 1:nrow(x)
    nms <- rownames(x)
    x <- na.omit(x)
    xx <- na.omit(xx)
    m2 <- ncol(xx)
    n2 <- nrow(xx)	
    n <- nrow(x)
    m <- ncol(x)		
    if (m2!=m) stop ("Error in dimensions of x and xx")
    nas <- na.action(x)
    nullans <- !is.null(nas)
  }
  else {
    stop("no matrix object")
  }
  
  if (is.null(n) && is.null(m)) 
    stop("ERROR IN THE DATA DIMENSIONS")
  if (is.null(m) && is.null(m2)) 
    stop("ERROR IN THE DATA DIMENSIONS")
  mdist=matrix(NA,ncol=n2,nrow=n2)
  for (i in 1:(n2-1)){for (j in (i+1):n2){
    mdist[i,j]<-mdist[j,i]<-sqrt(sum((xx[i,]-xx[j,])^2))
  }}
  if (is.null(h))   {
    h<-0.15
    hq2=quantile(mdist,probs=h,na.rm=TRUE)
    #print("es nulo h")  
  }
  else {
    #cat("no es nulo h ",h,"\n")    
    if (is.numeric(h))    hq2<-h  
    else hq2=quantile(mdist,probs=as.numeric(h),na.rm=TRUE)
  }	
  kern=function(x,y,h=hq2){exp(-sum((x-y)^2)/h^2)}
  K02=rep(1,nrow(x))
  K01=rep(1,nrow(xx))
  M1=matrix(NA,nrow=n2,ncol=n2)
  M2=matrix(NA,nrow=n,ncol=n2)
  M=array(NA,dim=c(n,n2,n2))
  for (i in 1:n2){for (j in i:n2){
    if (i==j) M1[i,i]=K01[i] else M1[i,j]<-M1[j,i]<-kern(xx[i,],xx[j,],h=hq2)
  }}
  same.dim <- FALSE
  if (n==n2 & m==m2){ same.dim <- TRUE}
  if (same.dim)
    if (all(x==xx)) {
      M2=M1
    } else {same.dim <- FALSE}
  if (!same.dim){	
    for (i in 1:n){for (j in 1:n2){
      if (all(x[i,] == xx[j,])) M2[i,j]=K02[i] else M2[i,j]=kern(x[i,],xx[j,],h=hq2)
    }}}  
  for (i in 1:n){for (j in 1:n2){for (k in 1:n2){	
    M[i,j,k]<-(K02[i]+M1[j,k]-M2[i,j]-M2[i,k])/(sqrt(K02[i]+K01[j]-2*M2[i,j])*sqrt(K02[i]+K01[k]-2*M2[i,k]))
  }}}
  l=which(!is.finite(M))
  M[l]=NA
  dep=1-sqrt(apply(M,1,sum,na.rm=TRUE))/n
  if (scale) {
    MO=array(NA,dim=c(n2,n2,n2))
    for (i in 1:n2){for (j in 1:n2){for (k in 1:n2){	
      MO[i,j,k]<-(K02[i]+M1[j,k]-M1[i,j]-M1[i,k])/(sqrt(K02[i]+K01[j]-2*M1[i,j])*sqrt(K02[i]+K01[k]-2*M1[i,k]))
    }}}
    l=which(!is.finite(MO))
    MO[l]=NA
    dep2=1-sqrt(apply(MO,1,sum,na.rm=TRUE))/n2
    mn <- min(dep2, na.rm = TRUE)
    mx <- max(dep2, na.rm = TRUE)
    dep = as.vector(dep/mx)
  }
  if (nullans) {
    ans1 <- rep(NA, len = m0)
    ans1[-nas] <- dep
    dep <- ans1
  }
  names(dep)=nms
  k = which.max(dep)
  med = matrix(x[k,],nrow=1)
  nl = length(trim)
  lista=vector("list",nl)
  tr <- paste("KFSD.tr", round(trim * 100,2), "%", sep = "")
  if (nl>1) names(lista)=paste0("tr",round(trim*100,2))
  mtrim = matrix(NA, nrow = nl, ncol = m)
  for (j in 1:length(trim)) {
    lista[[j]] = which(dep >= quantile(dep, probs = trim[j], na.rm = TRUE))
    if (length(lista[[j]])==1) {
      mtrim[j,]<-x[lista[[j]],]
      if (draw) {draw=FALSE;warning("Too few curves in mtrim. The plot is not drawn")}
    }
    else mtrim[j,]=apply(x[lista[[j]],],2,mean)       
  }
  
  rownames(med) <- "KFSD.med"
  rownames(mtrim) <- tr
  out <- list(median = med, lmed = k, mtrim = mtrim, ltrim = if (nl==1) unlist(lista) else lista, 
              dep = dep, h = h, hq=hq2,name="KFSD")
  if (scale) 
    out$dscale = mx
  class(out) <- "mdepth"
  if (draw) {
    plot.mdepth(x,xx,out)
  }
  return(invisible(out))	
}

#' @rdname depth.mdata
#' @export 
mdepth.FSD=function (x, xx = x, 
                     trim = 0.25, scale = FALSE, draw = FALSE){
  if (is.matrix(x) & is.matrix(xx)) {
    if (is.null(rownames(x))) 
      rownames(x) <- 1:nrow(x)
    nms <- rownames(x)
    m0 <- nrow(x)
    x <- na.omit(x)
    xx <- na.omit(xx)
    nas <- na.action(x)
    nullans <- !is.null(nas)
    n <- nrow(x)
    m <- ncol(x)
    m2 <- ncol(xx)
    n2 <- nrow(xx)
    if (m2!=m) stop("Error in data dimensions")
  }
  else {
    stop("no matrix object")
  }
  if (is.null(n) && is.null(m)) 
    stop("ERROR IN THE DATA DIMENSIONS")
  if (is.null(m) && is.null(m2)) 
    stop("ERROR IN THE DATA DIMENSIONS")
  kern=function(x,y){sum(x*y)}
  K02=diag(x%*%t(x))
  K01=diag(xx%*%t(xx))
  M1=matrix(NA,nrow=n2,ncol=n2)
  M2=matrix(NA,nrow=n,ncol=n2)
  M=array(NA,dim=c(n,n2,n2))
  if (scale) MO=array(NA,dim=c(n2,n2,n2))
  for (i in 1:n2){for (j in i:n2){
    if (i==j) M1[i,i]=K01[i] else M1[i,j]<-M1[j,i]<-kern(xx[i,],xx[j,])
  }}
  same.dim <- FALSE
  if (n==n2 & m==m2){ same.dim <- TRUE}
  if (same.dim)
    if (all(x==xx)) {
      M2=M1
    } else {same.dim <- FALSE}
  if (!same.dim){
    for (i in 1:n){for (j in 1:n2){
      if (all(x[i,] == xx[j,])) M2[i,j]=K02[i] else M2[i,j]=kern(x[i,],xx[j,])
    }}
  }
  for (i in 1:n){for (j in 1:n2){for (k in 1:n2){
    if (all(x[i,] == xx[j,]) | all(x[i,] == xx[k,])) M[i,j,k]=NA else M[i,j,k]<-(K02[i]+M1[j,k]-M2[i,j]-M2[i,k])/(sqrt(K02[i]+K01[j]-2*M2[i,j])*sqrt(K02[i]+K01[k]-2*M2[i,k]))
  }}}
  l=which(!is.finite(M))
  M[l]=NA
  dep=1-sqrt(apply(M,1,sum,na.rm=TRUE))/n
  if (scale) {
    for (i in 1:n2){for (j in 1:n2){for (k in 1:n2){
      if (all(xx[i,] == xx[j,]) | all(xx[i,] == xx[k,])) MO[i,j,k]=NA else MO[i,j,k]<-(K01[i]+M1[j,k]-M1[i,j]-M1[i,k])/(sqrt(K01[i]+K01[j]-2*M1[i,j])*sqrt(K01[i]+K01[k]-2*M1[i,k]))
    }}}
    l=which(!is.finite(MO))
    MO[l]=NA
    dep2=1-sqrt(apply(MO,1,sum,na.rm=TRUE))/n2
    mn <- min(dep2, na.rm = TRUE)
    mx <- max(dep2, na.rm = TRUE)
    dep = as.vector(dep/mx)
  }
  if (nullans) {
    ans1 <- rep(NA, len = m0)
    ans1[-nas] <- dep
    dep <- ans1
  }
  names(dep)=nms
  k = which.max(dep)
  med = matrix(x[k,],nrow=1)
  nl = length(trim)
  lista=vector("list",nl)
  tr <- paste("FSD.tr", round(trim * 100,2), "%", sep = "")
  if (nl>1) names(lista)=paste0("tr",round(trim*100,2))
  mtrim = matrix(NA, nrow = nl, ncol = m)
  for (j in 1:length(trim)) {
    lista[[j]] = which(dep >= quantile(dep, probs = trim[j], na.rm = TRUE))
    if (length(lista[[j]])==1) {
      mtrim[j,]<-x[lista[[j]],]
      if (draw) {draw=FALSE;warning("Too data in mtrim. The plot is not drawn")}
    }
    else mtrim[j,]=apply(x[lista[[j]],],2,mean)       
  }
  
  rownames(med) <- "FSD.med"
  rownames(mtrim) <- tr
  out <- list(median = med, lmed = k, mtrim = mtrim, ltrim = if (nl==1) unlist(lista) else lista, 
              dep = dep,name="FSD")
  if (scale) 
    out$dscale <- mx
  class(out) <- "mdepth"
  
  if (draw) {
    plot.mdepth(out)
  }
  return(invisible(out))
}

mdepth.MB <-function (x, xx = NULL,trim=0.25, scale=FALSE, draw =FALSE, 
                      grayscale = FALSE, band = FALSE, band.limits = NULL, lty = 1, lwd = 2, 
                      col = NULL, cold = NULL, colRef = NULL, ylim = NULL, cex = 1, ...)
{
  x <- data.matrix(x)
  n <- nrow(x)
  d <- ncol(x)
  depth.ori<-NULL
  fdataori<-xx
  fdataobj<-x
  tt<-1:d
  #    dtt<-c(0,dtt,0)/sum(dtt)*d
  if (length(xx) == 0) {
    if (ncol(x) == 1) {
      x <- t(x)
    }
    depth <- matrix(0, 1,n)
    ordered.matrix <- x
    if (n > 1) {
      for (columns in 1:d) {
        ordered.matrix[, columns] <- sort(x[, columns])
        for (element in 1:n) {
          index.1 <- length(which(ordered.matrix[, columns] <
                                    (x[element, columns])))
          index.2 <- length(which(ordered.matrix[, columns] <=
                                    (x[element, columns])))
          multiplicity <- index.2 - index.1
          depth[element] <- depth[element]+ index.1 *
            (n - (index.2)) + multiplicity * (n - index.2 +
                                                index.1) + choose(multiplicity, 2)
        }
      }
      depth <- depth/(d * choose(n, 2))            
    }
    if (n == 1) {
      deepest <- x
      depth <- 0
    }
    ordering <- order(depth, decreasing = TRUE)
    #########################################
    if (draw) {
      par(mar = c(4, 5, 3, 3), xpd = FALSE)
      fobj <- t(x[ordering[n:2], ,drop=FALSE ])
      lwdd <- lwd[1]
      if (is.null(cold)) {
        cold <- 2
      }
      if (is.null(ylim)) {
        ylim <- range(x)
      }
      if (band) {
        lty <- lty[1]
        if (is.null(band.limits)) {
          band.limits <- c(0.5, 1)
        }
        else {
          band.limits <- unique(sort(band.limits))
        }
        if (floor(n * (band.limits[1])) < 2) {
          stop("Check the limits. The band must contain at least 2 curves...")
        }
        no.poly <- length(band.limits)
        if (is.null(col)) {
          if (grayscale) {
            color <- rev(gray((1:no.poly)/(no.poly +
                                             1)))
          }
          else {
            color <- 3:(2 + no.poly)
          }
        }
        else {
          if (length(col) < no.poly) {
            color <- rep(col, length.out = no.poly)[1:no.poly]
          }
          else {
            color <- col[1:no.poly]
          }
        }
        par(mar = c(4, 5, 5, 3), xpd = FALSE)
        #                matplot(Gene.Expression, ylim = ylim, type = "l",lty = 0, ...)
        matplot(fobj, ylim = ylim, type = "l",lty = 0, ...)
        for (poly in no.poly:1) {
          limit <- band.limits[poly]
          no.points <- floor(limit * n)
          xx2 <- t(x[ordering[no.points:1],])
          upper <- apply(xx2, 1, max)
          lower <- apply(xx2, 1, min)
          polygon(c(tt, rev(tt)), c(upper, rev(lower)),
                  col = color[poly])
        }
      }
      else {
        if (is.null(col)) {
          if (grayscale) {
            color <- rev(gray((1:n)/(n + 1)))
          }
          else {
            color <- rep(8, n)
          }
        }
        else {
          if (length(col) < n) {
            color <- rep(col, length.out = n)[n:1]
          }
          else {
            color <- col[n:1]
          }
        }
        matplot(fobj, type = "l", ylim = ylim,lty = lty, lwd = lwd, col = color[n:2], ...)
      }
      lines(x[ordering[1], ], lty = lty, lwd = lwdd, col = cold)
      par(xpd = TRUE)
      legend("top", legend = "deepest sample", col = cold,
             lty = lty, lwd = lwdd, cex = cex)
      if (band) {
        legend("top", inset = -0.1 * cex, horiz = TRUE,
               legend = band.limits, col = color, pch = 15,
               title = "Proportion of central samples", cex = cex,
               bty = "n")
      }
    }
    #########################################
  }
  else {
    xRef<-xx
    xRef <- data.matrix(xRef)
    if (ncol(xRef) != d) {
      stop("Dimensions of x and xRef do not match")
    }
    n0 <- nrow(xRef)
    if (ncol(x) == 1) {
      x <- t(x)
    }
    depth <- matrix(0, 1,n)
    depth.ori <- matrix(0, 1,n0)
    ordered.matrix <- xRef
    if (n0 > 1) {
      for (columns in 1:d) {
        ordered.matrix[, columns] <- sort(xRef[, columns])
        for (element in 1:n) {
          index.1 <- length(which(ordered.matrix[, columns] <
                                    x[element, columns]))
          index.2 <- length(which(ordered.matrix[, columns] <=
                                    x[element, columns]))
          multiplicity <- index.2 - index.1
          depth[element] <- depth[element]+(index.1 +
                                              multiplicity) * (n0 - index.1 - multiplicity) +
            multiplicity * (index.1 + (multiplicity -
                                         1)/2)
        }
        for (element in 1:n0) {
          index.1 <- length(which(ordered.matrix[, columns] <
                                    xRef[element, columns]))
          index.2 <- length(which(ordered.matrix[, columns] <=
                                    xRef[element, columns]))
          multiplicity <- index.2 - index.1
          depth.ori[element] <- depth.ori[element]+(index.1 +
                                                      multiplicity) * (n0 - index.1 - multiplicity) +
            multiplicity * (index.1 + (multiplicity -
                                         1)/2)
        }
      }
      depth <- depth/(d * choose(n0, 2))
      depth.ori <- depth.ori/(d * choose(n0, 2))
    }
    if (n == 1) {
      deepest <- x
      depth <- 0
    }
    ordering <- order(depth, decreasing = TRUE)
    if (draw) {
      par(mar = c(4, 5, 3, 3), xpd = FALSE)
      fdobj<- t(xRef)
      if (is.null(colRef)) {
        colRef <- 4
      }
      else {
        colRef <- colRef[1]
      }
      if (is.null(cold)) {
        cold <- 2
      }
      if (is.null(ylim)) {
        ylim <- range(x, xRef)
      }
      lwdd <- lwd[1]
      if (band) {
        lty <- lty[1]
        if (is.null(band.limits)) {
          band.limits <- c(0.5, 1)
        }
        band.limits <- unique(sort(band.limits))
        if (floor(n * (band.limits[1])) < 2) {
          stop("Check the limits. The band must contain at least 2 curves...")
        }
        no.poly <- length(band.limits)
        par(mar = c(4, 5, 5, 3), xpd = FALSE)
        matplot(fdobj, type = "l", ylim = ylim,
                lty = lty, lwd = lwd/2, col = colRef, ...)
        if (is.null(col)) {
          if (grayscale) {
            color <- rev(gray((1:no.poly)/(no.poly +
                                             1)))
          }
          else {
            color <- 5:(no.poly + 4)
          }
        }
        else {
          if (length(col) < no.poly) {
            color <- rep(col, length.out = no.poly)[1:no.poly]
          }
          else {
            color <- col[1:no.poly]
          }
        }
        for (poly in no.poly:1) {
          limit <- band.limits[poly]
          no.points <- floor(limit * n)
          xx2<- t(x[ordering[no.points:1],])
          upper <- apply(xx2, 1, max)
          lower <- apply(xx2, 1, min)
          polygon(c(tt, rev(tt)), c(upper, rev(lower)),
                  col = color[poly])
        }
        lines(x[ordering[1], ], lty = lty, lwd = lwdd,
              col = cold)
      }
      else {
        matplot(fdobj, type = "l", ylim = ylim,
                lty = lty, lwd = lwd/2, col = colRef, ...)
        if (is.null(col)) {
          if (grayscale) {
            color <- rev(gray(1/(n + 1):(n/(n + 1))))
          }
          else {
            color <- rep(8, n)
          }
        }
        else {
          if (length(col) < n) {
            color <- rep(col, length.out = n)[n:1]
          }
          else {
            color <- col[n:1]
          }
        }
        matlines(t(x[ordering[n:2], ]), lwd = lwd, col = color[n:2],
                 lty = lty)
        lines(x[ordering == 1, ], lwd = lwdd, lty = lty[1],
              col = cold)
      }
      par(xpd = TRUE)
      legend("top", legend = c("deepest sample", "reference set"),
             col = c(cold, colRef), lty = lty, lwd = c(lwdd,
                                                       lwd/2), cex = cex)
      if (band) {
        legend("top", inset = -0.1 * cex, horiz = TRUE,
               legend = band.limits, col = color, pch = 15,
               title = "Proportion of central samples", cex = cex,
               bty = "n")
      }
    }
    if (scale)    depth.ori<-depth.ori/max(depth.ori)
  }
  med<-fdataobj[ordering[1]]
  lista = which(depth >= quantile(depth[1,], probs = trim, na.rm = TRUE))
  mtrim= apply(x[lista, ], 2, mean)
  #  tr<-paste("band.tr",trim*100,"\u0025",sep="")
  #  names2<-fdataobj[["names"]]
  #  med$names$main<-"depth.mode median"
  #  names2$main<-paste("depth.mbandtrim ",trim*100,"\u0025",sep="")
  #  mtrim<-fdata(mtrim,tt,rtt,names2)
  #  rownames(med)<-"mode.med"
  #  rownames(mtrim)<-tr
  if (scale)    depth[1,]<-depth[1,]/max(depth[1,])
  return(invisible(list("median"=med,"lmed"=ordering[1],ordering=ordering,
                        "mtrim"=mtrim,"ltrim"=lista,"dep"=depth[1,],"dep.ori"=depth.ori)))
}


# multivariate version
#' @rdname depth.mdata
#' @export
mdepth.FM=function(x,xx=x,scale=FALSE,dfunc="TD1"){
  n<-nrow(x)
  m<-ncol(xx)
  d<-matrix(NA,nrow=n,ncol=m)
  Fn<-list()
  print(m)
  par.dfunc<-list()
  for (i in 1:m)   {
    if (dfunc %in% c("TD1","Liu1","FM1")){
      Fn[[i]]=ecdf(xx[,i])
      par.dfunc$x=x[,i]
      par.dfunc$Fn=Fn[[i]]
      d[,i]=do.call(dfunc,par.dfunc)
    }
    else     {
      par.dfunc$x=x[,i]
      par.dfunc$xx=xx[,i]
      d[,i]=do.call(dfunc,par.dfunc)
    }
  }
  d<-apply(d,1,mean,na.rm=TRUE)[1:n]
  return(list("dep"=d,"Fn"=Fn))  #
}

mdepth.Liu=function(x,xx=x,xeps=1e-15,scale=FALSE){
  n<-nrow(x)
  m<-ncol(x)
  d<-matrix(NA,nrow=n,ncol=m)
  Fn<-list() 
  for (i in 1:m){
    Fn[[i]]=ecdf(xx[,i])
    d[,i]=Fn[[i]](x[,i])*(1-Fn[[i]](x[,i]-xeps))
  }
  if (scale){   d<-d*4   }
  d<-apply(d,1,mean,na.rm=TRUE)[1:n]
  return(list("dep"=d,"Fn"=Fn))  
}     

# multivariate version
#' @rdname depth.mdata
#' @export
mdepth.TD=function(x,xx=x,xeps=1e-15,scale=FALSE){
  #print("entra TD")
  n<-nrow(x)
  m<-ncol(x)
  d<-matrix(NA,nrow=n,ncol=m)
  Fn<-list() 
  for (i in 1:m){
    Fn[[i]]=ecdf(xx[,i])
    d[,i]=pmin(Fn[[i]](x[,i]),(1-Fn[[i]](x[,i]-xeps)))
  }
  if (scale){ d<-d*2      }
  d<-apply(d,1,mean,na.rm=TRUE)[1:n]
  return(list("dep"=d,"Fn"=Fn)) #"dep.ori"=d2,
}


isiindata<-function(i,x){
  n<-nrow(x)
  ind<-rep(NA,len=n)
  for (j in 1:n){
    ind[j]<-all(i==x[j,])
  }
  ind
}
triang=function(Q,P1,P2,P3){
  q1=Q-P1;a1=atan2(q1[2],q1[1])
  q2=Q-P2;a2=atan2(q2[2],q2[1])
  q3=Q-P3;a3=atan2(q3[2],q3[1])
  a=sort(c(a1,a2,a3))
  return(((a[2]-a[1])<=pi) && ((a[3]-a[2])<=pi) && ((a[3]-a[1])>=pi))
}

ntriang<-function(q,puntos,plot=FALSE){
  #print("entra ntriang")
  #   if (dim(x)[2]!=2) stop("The dimension of the data points is not a 2-column matrix")
  #   if (dim(x)[1]<3) stop("The dimension of the data points must be at least 3 rows.")
  #   if (length(q)!=2) stop("El punto q debe tener dos componentes")
  n=dim(puntos)[1]
  #nq=q-puntos
  nq<-sweep(puntos,2,q,"-")
  alfa=atan2(nq[,2],nq[,1])
  a=sort(alfa)
  di<-ci<-numeric(n);di=numeric(n)
  #ci=di=numeric(n)
  ind.a<-a>0
  k<-ifelse(sum(ind.a)>0,min(which(ind.a)),n)
  #k=min(min(which(a>0)),n)
  #if (is.infinite(k)) k=n
  lista=1:(k-1)
  for (i in 1:n){
    ci[i]=max(which((a-a[i])<=pi))
    di[i]=min(c(which((a-a[i])>=pi)),n+1)
  }
  e1=ci[k]*(ci[k]+1-2*k)+k*(3-k)-2
  d1=n*(n-k+1)+0.5*k*(k-1)-0.5*n*(n+1)
  e=(lista+n-1)*ci[lista]
  f=(ci[lista]+1)*ci[lista]+(di[lista]-1)*(di[lista]-2*lista-2)
  sume=sum(e)+0.5*n*e1
  sumf=sum(f/2)+n*d1
  ttriang=n*(n-1)*(n-2)/6 #n sobre 3
  return(list("triang"=sume-sumf,"ttriang"=ttriang))
}

#' @rdname depth.mdata
#' @export
mdepth.SD=function(x,xx=NULL,scale=FALSE){
  if (dim(x)[2]!=2) stop("The dimension of the data points is not a 2-column matrix")
  if (dim(x)[1]<3) stop("The dimension of the data points must be at least 3 rows.")
  if (is.data.frame(x)) {
    x <- data.matrix(x)
  }
  if (is.null(xx)) {
    xx<-x
    al<-TRUE
  }
  else {
    al<-FALSE
    xx <- data.matrix(xx)
  }
  nam<-colnames(xx)
  if ( is.null(rownames(x)))  rownames(x)<-1:nrow(x)
  nms<-rownames(x)
  m0<-nrow(x)
  xx<-na.omit(xx)
  x<-na.omit(x)
  nas<-na.action(x)
  nullans<-!is.null(nas) 
  
  n=nrow(x)
  nn=nrow(xx)
  rownames(xx)<-NULL#1:nn
  d=ncol(xx)
  dep=numeric(n)
  if (al){
    # print("entra al l")
    ttriang=nn*(nn-1)*(nn-2)/6    
    #    print(al);    print("all(x==xx)")
    fac<-(nn-1)*(nn-2)/2
    for (i in 1:n){
      #        dep[i]=(ntriang(x[i,],x[-i,],TRUE)$triang+(n-1)*(n-2)/2)/ttriang
      dep[i]<-ntriang(x[i,],xx[-i,])$triang      
      #        cat("i ",i);        print(dep[i])
    }
    dep<-(dep+fac)/ttriang
  }
  else{
    ttriang=(nn+1)*(nn)*(nn-1)/6
    fac<-(nn)*(nn-1)/2
    for (i in 1:n){
      dep[i]=ntriang(x[i,],xx)$triang
    }      
    dep<-(dep+fac)/ttriang 
  }
  if (scale){ 
    dep<-dep/max(dep)
  } 
  if  (nullans){
    ans1<-rep(NA,len=m0)
    ans1[-nas] <-dep
    dep<-ans1      
  }
  names(dep)<-nms         
  return(invisible(list(dep = dep,x=x,xx=xx)))
}




FM1<-function(x,Fn,scale=FALSE) {
  d=1-abs(0.5-Fn(x))
  if (scale){     d<-(d-.5)*2    }  
  d
}
Liu1<-function(x,Fn,xeps=1e-15,scale=FALSE) {
  d=Fn(x)*(1-Fn(x-xeps))
  if (scale){   d<-d*4   }
  d
}
TD1<-function(x,Fn,xeps=1e-15,scale=FALSE) {
  d=pmin(Fn(x),(1-Fn(x-xeps)))      #HS1D
  if (scale){ d<-d*2      }
  d
}
LD1<-function(x,xx=x,scale=TRUE) {
  mdepth.LD(matrix(x,ncol=1),matrix(xx,ncol=1),scale=scale)$dep
}      
MhD1 <- function(x,xx=x,scale=FALSE){
  if (is.vector(x)) {
    D<- (1+(x-(mean(xx)))^2/sd(xx)^2)
  }
  else{stop("x is no a vector") }
  ans <- 1/D
  ans
}   



# Univariate version
mdepth.FM1=function(x,xx=x,scale=FALSE){
  isx<-is.vector(x)
  isxx<-is.vector(xx)
  if (!isx)  stop("x object is not a vector")
  if (!isxx)  stop("xx object is not a vector")
  n<-length(x)
  m<-length(xx)
  if (n!=m)  stop("The length of x is not the same of xx")
  Fn<-list()
  Fn=ecdf(xx)
  d=1-abs(0.5-Fn(x))
  if (scale){     d<-(d-.5)*2    }   
  return(list("dep"=d,"Fn"=Fn))  #
}

# univariate version
mdepth.TD1=function(x,xx=x,xeps=1e-15,scale=FALSE){
  #print("entra TD")
  isx<-is.vector(x)
  isxx<-is.vector(xx)
  if (!isx)  stop("x object is not a vector")
  if (!isxx)  stop("xx object is not a vector")  
  n<-length(x)
  m<-length(xx)
  if (n!=m)  stop("The length of x is not the same of xx")
  Fn=ecdf(xx)
  d=pmin(Fn(x),(1-Fn(x-xeps)))
  if (scale){ d<-d*2      }
  return(list("dep"=d,"Fn"=Fn))  #
}

# univariate version
mdepth.Liu1=function(x,xx=x,xeps=1e-15,scale=FALSE){
  isx<-is.vector(x)
  isxx<-is.vector(xx)
  if (!isx)  stop("x object is not a vector")
  if (!isxx)  stop("xx object is not a vector")
  n<-length(x)
  m<-length(xx)
  if (n!=m)  stop("The length of x is not the same of xx")
  Fn=ecdf(xx)
  d=Fn(x)*(1-Fn(x-xeps))
  if (scale){  d<-d*4      }
  return(list("dep"=d,"Fn"=Fn))  
}
