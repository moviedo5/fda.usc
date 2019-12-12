#' @title DTW: Dynamic time warping
#'
#' @description Computes distances time warping for functional data
#'
#' @param fdata1 Functional data 1 or curve 1. If \code{fdata} class, the dimension of \code{fdata1$data} object is (\code{n1} x \code{m}), where \code{n1} is the number of curves and \code{m} are the points observed in each curve.
#' @param fdata2 Functional data 2 or curve 2. If \code{fdata} class, the dimension of \code{fdata2$data} object is (\code{n2} x \code{m}), where \code{n2} is the number of curves and \code{m} are the points observed in each curve.
#' @param p Lp norm, by default it uses \code{p = 2}
#' @param w Vector of weights with length \code{m}, If \code{w = 1} approximates the metric Lp by Simpson's rule. By default it uses \code{w = 1} 
#' @param lambda \code{numeric} lambda value (0 by default)
#' @param wmax \code{numeric} maximum value of weight, (1 by default)
#' @param g  \code{numeric} \code{g=0} (constant), \code{0.05} (linear) by default, 0.25 \code{sigmoid}, 3 two weight values
#' @param nu  \code{numeric} constant value, (0 by default)
#' 
#' @aliases metric.DTW metric.WDTW metric.TWED
#' @details Three optins:
#' \itemize{
#'   \item DTW: Dynamic time warping
#'   \item WDTW: Weight Dynamic time warping
#'   \item TWED: twed   
#' }
# Futher information in \url{http://timeseriesclassification.com/}
#' 
#' @return DTW matrix 
#'
#' @examples
#' \dontrun{
#' data(tecator)
#' metric.DTW(tecator$absorp.fdata[1:4,])
#' ab=tecator[[1]]
#' D1=fda.usc:::DTW(ab$data[1,],ab$data[2,],p=2)
#' aa1=fda.usc:::findPath(D1$D)
#' D2=fda.usc:::DTW(ab$data[1,],ab$data[2,],p=2,w=5)
#' aa2=fda.usc:::findPath(D2$D)
#' D3=fda.usc:::WDTW(ab$data[1,],ab$data[2,],p=2,g=0.05) 
#' aa3=fda.usc:::findPath(D3$D)
#' D4=fda.usc:::TWED(ab$data[1,],ab$data[2,],p=2,lambda=0,nu=0)
#' aa4=fda.usc:::findPath(D4$D)
#' par(mfrow=c(2,2))
#' plot(c(ab[1:2]))
#' segments(ab$argvals[aa1[,1]],ab[1]$data[aa1[,1]],ab$argvals[aa1[,2]],ab[2]$data[aa1[,2]])
#' plot(c(ab[1:2]))
#' segments(ab$argvals[aa2[,1]],ab[1]$data[aa2[,1]],ab$argvals[aa2[,2]],ab[2]$data[aa2[,2]],col=2)
#' plot(c(ab[1:2]))
#' segments(ab$argvals[aa3[,1]],ab[1]$data[aa3[,1]],ab$argvals[aa3[,2]],ab[2]$data[aa3[,2]],col=3)
#' plot(c(ab[1:2]))
#' segments(ab$argvals[aa4[,1]],ab[1]$data[aa4[,1]],ab$argvals[aa4[,2]],ab[2]$data[aa4[,2]],col=4)
#' }
#' 
#' @references
#' Jeong, Y. S., Jeong, M. K., & Omitaomu, O. A. (2011). Weighted dynamic time warping for 
#' time series classification. \emph{Pattern Recognition}, 44(9), 2231-2240
#' 
#' @seealso See also  \code{\link{semimetric.basis}} and \code{\link{semimetric.NPFDA}}
#' 
#' @author
#'   Manuel Febrero-Bande, Manuel Oviedo de la Fuente \email{manuel.oviedo@@usc.es}
#'   
#'@keywords cluster 

#' @rdname metric.DTW 
#' @export 
metric.DTW=function(fdata1, fdata2=NULL, p=2,
                    w=min(ncol(fdata1),ncol(fdata2))){
C1=match.call()
same=FALSE
if (is.null(fdata2)) {
  fdata2 <- fdata1
  same = TRUE
}
if (is.fdata(fdata1)){
	nas1=is.na.fdata(fdata1)
	if (any(nas1)) warning("fdata1 contains ",sum(nas1), " curves with some NA value \n")
	if (is.null(fdata2)) {fdata2=fdata1;same=TRUE} else if (!is.fdata(fdata2)) {fdata2=fdata(fdata2,fdata1$argvals,fdata1$rangeval,fdata1$names)}
	nas2=is.na.fdata(fdata2) 
	if (any(nas2)) warning("fdata2 contains ",sum(nas1), " curves with some NA value \n")
	data1=fdata1[["data"]]
	data2=fdata2[["data"]]
} else {
	data1=as.matrix(fdata1)
	if (is.null(fdata2)) {data2=data1;same=TRUE} else {data2=as.matrix(fdata2)}
}
#if (all(data1==data2)) same=TRUE
D=matrix(0,nrow=nrow(data1),ncol=nrow(data2))
if (same){
for (i in 1:(nrow(data1)-1)){
	for (j in (i+1):nrow(data2)){
	D[i,j]=DTW(data1[i,],data2[j,],p=p,w=w)$dist
	}}
D=(D+t(D))	
} else {
for (i in 1:nrow(data1)){
	for (j in 1:nrow(data2)){
	D[i,j]=DTW(data1[i,],data2[j,],p=p,w=w)$dist
	}}
}
attr(D,"call")="metric.DTW"
attr(D,"par.metric")=list(p=p,w=w)
return(D)
}

#' @rdname metric.DTW
#' @export 
metric.WDTW=function(fdata1,fdata2=NULL,p=2,
                       w=min(ncol(fdata1),ncol(fdata2)),wmax=1,g=0.05){
C1=match.call()
same=FALSE
if (is.null(fdata2)) {
  fdata2 <- fdata1
  same = TRUE
}
if (is.fdata(fdata1)){
	nas1=is.na.fdata(fdata1)
	if (any(nas1)) warning("fdata1 contains ",sum(nas1), " curves with some NA value \n")
	if (is.null(fdata2)) {fdata2=fdata1;same=TRUE} else if (!is.fdata(fdata2)) {fdata2=fdata(fdata2,fdata1$argvals,fdata1$rangeval,fdata1$names)}
	nas2=is.na.fdata(fdata2) 
	if (any(nas2)) warning("fdata2 contains ",sum(nas1), " curves with some NA value \n")
	data1=fdata1[["data"]]
	data2=fdata2[["data"]]
} else {
	data1=as.matrix(fdata1)
	if (is.null(fdata2)) {data2=data1;same=TRUE} else {data2=as.matrix(fdata2)}
}
#if (all(data1==data2)) same=TRUE
D=matrix(0,nrow=nrow(data1),ncol=nrow(data2))
if (same){
for (i in 1:(nrow(data1)-1)){
	for (j in (i+1):nrow(data2)){
	D[i,j]=WDTW(data1[i,],data2[j,],p=p,w=w,wmax=wmax,g=g)$dist
	}}
D=(D+t(D))	
} else {
for (i in 1:nrow(data1)){
	for (j in 1:nrow(data2)){
	D[i,j]=WDTW(data1[i,],data2[j,],p=p,w=w,wmax=wmax,g=g)$dist
	}}
}
attr(D,"call")="metric.WDTW"
attr(D,"par.metric")=list(p=p,w=w,wmax=wmax,g=g)
return(D)
}

#' @rdname metric.DTW
#' @export 
metric.TWED=function(fdata1, fdata2=NULL,p=2,lambda=1,nu=0.05){
C1=match.call()
same=FALSE
if (is.null(fdata2)) {
  fdata2 <- fdata1
  same = TRUE
}
if (is.fdata(fdata1)){
	nas1=is.na.fdata(fdata1)
	if (any(nas1)) warning("fdata1 contains ",sum(nas1), " curves with some NA value \n")
	if (is.null(fdata2)) {fdata2=fdata1;same=TRUE} else if (!is.fdata(fdata2)) {fdata2=fdata(fdata2,fdata1$argvals,fdata1$rangeval,fdata1$names)}
	nas2=is.na.fdata(fdata2) 
	if (any(nas2)) warning("fdata2 contains ",sum(nas1), " curves with some NA value \n")
	data1=fdata1[["data"]]
	data2=fdata2[["data"]]
} else {
	data1=as.matrix(fdata1)
	if (is.null(fdata2)) {data2=data1;same=TRUE} else {data2=as.matrix(fdata2)}
}
#if (all(data1==data2)) same=TRUE
D=matrix(0,nrow=nrow(data1),ncol=nrow(data2))
if (same){
for (i in 1:(nrow(data1)-1)){
	for (j in (i+1):nrow(data2)){
	D[i,j]=TWED(data1[i,],data2[j,],p=p,lambda=lambda,nu=nu)$dist
	}}
D=(D+t(D))	
} else {
for (i in 1:nrow(data1)){
	for (j in 1:nrow(data2)){
	D[i,j]=TWED(data1[i,],data2[j,],p=p,lambda=lambda,nu=nu)$dist
	}}
}
attr(D,"call")="metric.TWED"
attr(D,"par.metric")=list(p=p,lambda=lambda,nu=nu)
return(D)
}




############################
cumgamma=function(D,w=min(nrow(D),ncol(D))){
  n=nrow(D)
  m=ncol(D)
  DTW=matrix(0,nrow=n+1,ncol=m+1)
  DTW[1,]=Inf
  DTW[,1]=Inf
  DTW[1,1]=0
  for (i in 2:(n+1)){for (j in max(2,i-w):min(m+1,i+w)){DTW[i,j]=0}}
  for (i in 2:(n+1)){
    for (j in max(2,i-w):min(m+1,i+w)){
      DTW[i,j]=D[i-1,j-1]+min(DTW[i-1,j],DTW[i,j-1],DTW[i-1,j-1])
    }}
  return(DTW)
}

############################
findPath=function(D,i=nrow(D),j=ncol(D)){
  ca=c(i)
  cb=c(j)
  while(i>1 & j>1) {
    if (i>1 & j>1){
      aa=min(D[i-1,j-1],D[i-1,j],D[i,j-1])
      if (aa==D[i-1,j-1]) {i=i-1;j=j-1} else if (aa==D[i-1,j]){i=i-1} else {j=j-1}
    } else if (j>1){j=j-1} else {i=i-1}
    ca=c(ca,i)
    cb=c(cb,j)
  }
  aa=cbind(ca,cb)
  aa=aa[-nrow(aa),]-1
  return(aa)
}

############################
DTW=function(a,b,w=min(length(a),length(b)),p=2){
  D=outer(a,b,"-")
  D=abs(D)^p
  DTW=cumgamma(D,w=w)
  return(list(dist=DTW[length(a)+1,length(b)+1],D=DTW))
}

############################
WDTW=function(a,b,w=length(a),p=2,wmax=1,g=0.05){
  #g=0 (constante), 0.05 (linear), 0.25 (sigmoid), 3 (dos pesos)
  n=length(b)
  wei=outer(1:n,1:n,"-")
  wei=wmax/(1+exp(-g*(abs(w)-n/2)))
  M=outer(a,b,"-")
  D=abs(wei*M)^p
  DTW=cumgamma(D,w=w)
  return(list(dist=DTW[length(a)+1,length(b)+1],D=DTW))
}

############################
TWED=function(a,b,p=2,lambda=1,nu=1){
  n=length(a)
  m=length(b)
  D=matrix(0,nrow=n+1,ncol=m+1)
  D[1,1]=0
  D[2,1]=abs(a[1])^p
  D[1,2]=abs(b[1])^p
  for (i in 3:(n+1)){ D[i,1]=D[i-1,1]+abs(a[i-2]-a[i-1])^p}
  for (j in 3:(m+1)){ D[1,j]=D[1,j-1]+abs(b[i-2]-b[i-1])^p}
  for (i in 2:(n+1)){
    for (j in 2:(m+1)){
      if (i>2 & j>2) {d1=D[i-1,j-1]+2*nu*abs(i-j)+abs(a[i-1]-b[j-1])^p+abs(a[i-2]-b[j-2])^p} else {d1=D[i-1,j-1]+nu*abs(i-j)+abs(a[i-1]-b[j-1])^p}
      if (i>2) {d2=D[i-1,j]+abs(a[i-1]-a[i-2])^p+lambda+nu} else {d2=D[i-1,j]+abs(a[i-1])^p+lambda}
      if (j>2) {d3=D[i,j-1]+abs(b[j-1]-b[j-2])^p+lambda+nu} else {d3=D[i,j-1]+abs(b[j-1])^p+lambda}
      D[i,j]=min(d1,d2,d3)
    }}
  return(list(dist=D[length(a)+1,length(b)+1],D=D))
}
