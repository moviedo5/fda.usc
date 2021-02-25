predict.kmeans.fd <- function(object, newdata, ...){
  x = object$centers
  n = NROW(x)
  np <- NCOL(x)
  nn = NROW(newx)
  pp = NCOL(newx)
  if (pp != np) stop("newdata have wrong dimension")
  if (is.null(rownames(newx))) 
    rownames(newx) <- 1:nn
  par.metric <- attr(object$mdist, "par.metric")
  par.metric[["fdata1"]] <- newdata
  par.metric[["fdata2"]] <- x
  a1 <- attr(object$mdist, "call")
  nmdist <- do.call(a1, par.metric)
  #for (i in 1:nr) {  grupo[i] = which.min(d[(nr + 1):(nr + ngroups), i])}  
  grupo <- apply(nmdist,1,which.min)
  return(grupo)
}


# @name kmeans.fd
# @title K-Means Clustering for functional data
# 
# @description Perform k-means clustering on functional data.
# 
# @details The method searches the locations around which are grouped data (for a
# predetermined number of groups).\cr
# 
# If \code{ncl=NULL}, randomizes the initial centers, \code{ncl=2} using
# \code{kmeans.center.ini} function.\cr If \code{ncl} is an integer,
# indicating the number of groups to classify,\cr are selected \code{ncl}
# initial centers using \code{kmeans.center.ini} function.\cr If \code{ncl} is
# a vector of integers, indicating the position of the initial centers with
# \code{length(ncl)} equal to number of groups.\cr If \code{ncl} is a
# \code{fdata} class objecct, \code{ncl} are the initial centers curves with
# \code{nrow(ncl)} number of groups.\cr
# 
# @aliases kmeans.fd kmeans.center.ini kmeans.centers.update
# 
# @param fdataobj \code{\link{fdata}} class object.
# @param ncl See details section.
# @param metric Metric function, by default \code{\link{metric.lp}}.
# @param dfunc Type of depth measure, by default FM depth.
# @param max.iter Maximum number of iterations for the detection of centers.
# @param draw =TRUE, draw the curves in the color of the centers.
# @param par.dfunc List of arguments to pass to the \code{dfunc} function .
# @param par.ini List of arguments to pass to the \code{kmeans.center.ini} function .
# @param method Method for selecting initial centers. If
# \code{method}=\emph{"Sample"} (by default) takes \code{n} times a random
# selection by the \code{ncl} centers. The \code{ncl} curves with greater
# distance are the initial centers. If \code{method}=\emph{"Exact"} calculated
# all combinations (if < 1e+6) of \code{ncl} centers. The \code{ncl} curves with greater
# distance are the initial centers (this method may be too slow).
# @param cluster.size Minimum cluster size (by default is 5). If a cluster has fewer curves,
# it is eliminated and the process is continued with a less cluster.
# @param max.comb Maximum number of initial selection of
# centers (only used when \code{method="exact"}).
# @param par.metric List of arguments to pass to the \code{metric} function.
# @param group groups or classes
# @param \dots Further arguments passed to or from other methods.
# 
# @return Return:
# \itemize{
#  \item \code{cluster}{ Indexes of groups assigned.}  
# \item \code{centers}{ Curves centers.}
#  %\item{lcenters}{ Indexes of initial curves centers.}
# }
# 
# @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
# \email{manuel.oviedo@@udc.es}
# 
# @seealso See Also generic \link[stats]{kmeans} function.
# 
# @references Hartigan, J. A. and Wong, M. A. (1979). \emph{A K-means
# clustering algorithm}. Applied Statistics 28, 100 \-108.
# 
# @keywords cluster
# @examples
# \dontrun{
# data(phoneme)
# mlearn<-phoneme$learn[c(1:50,101:150,201:250),]
# 
# # Unsupervised classification
# out.fd1=kmeans.fd(mlearn,ncl=3,draw=TRUE)
# out.fd2=kmeans.fd(mlearn,ncl=3,draw=TRUE,method="exact")
# # Different Depth function
# ind=c(17,77,126)
# out.fd3=kmeans.fd(mlearn,ncl=mlearn[ind,],draw=FALSE,
# dfunc=func.trim.FM,par.dfunc=list(trim=0.1))
# out.fd4=kmeans.fd(mlearn,ncl=mlearn[ind,],draw=FALSE,
# dfunc=func.med.FM)
# group=c(rep(1,50),rep(2,50),rep(3,50))
# table(out.fd4$cluster,group)
# }
# 
# @rdname kmeans.fd
# @export
# 
 
# kmeans.fd=function(fdataobj,ncl=2,metric=metric.lp
#                    ,dfunc=func.trim.FM,max.iter=100
#                    ,par.metric=NULL,par.dfunc=list(trim=0.05)
#                    ,method="sample", cluster.size=5,draw=TRUE,...) {
# #if (is.data.frame(z)) z=data.matrix(z)
# #else if (is.vector(z))     z <- data.matrix(t(z))
# if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
#  # nas1<-is.na.fdata(fdataobj)
# nas1<-is.na(fdataobj)
#  if (any(nas1))  stop("fdataobj contain ",sum(nas1)," curves with some NA value \n")
# z<-fdataobj[["data"]]
# tt<-fdataobj[["argvals"]]
# rtt<-fdataobj[["rangeval"]]
# nr=nrow(z)
# nc=ncol(z)
# if (is.vector(ncl)) {
#   len.ncl=length(ncl)
#   if (len.ncl==1) {
#     par.ini<- list()
#     par.ini$fdataobj=fdataobj
#     par.ini$method=method
#     par.ini$ncl=ncl
#     par.ini$metric=metric
#     par.ini$draw=draw
#     par.ini$max.comb = 1e6
#     par.ini$max.iter = max.iter
#     if (!is.null(par.metric)) 
#        par.ini$par.metric<-par.metric
#     par.ini$... <- par.metric
#     out1=do.call(kmeans.center.ini,par.ini)
#     lxm<-out1$lcenters
#     out1$d=rbind(out1$z.dist,out1$z.dist[lxm,])
#     }  else {
#      ngroups=length(ncl)
#      lxm=ncl
#      xm=z[lxm,]
#      out1=list()
#      out1$fdataobj<-fdataobj
#      out1$ncl=len.ncl
#      if (is.null(par.metric)) 
#        par.metric=list("p"=2,"w"=1)
#      par.metric$fdata1<-fdataobj
#      mdist=do.call(metric,par.metric)
#      out1$z.dist<-mdist
#      out1$d=rbind(mdist,mdist[lxm,])
#      out1$centers<-fdataobj[ncl,]
#      out1$lcenters<-ncl
#      class(out1)="kmeans.fd"
#     }
#   } else if (is.fdata(ncl)) {   # fdata centers
#    lxm=NULL
#    xm=ncl[["data"]]
#    if (is.null(par.metric)) par.metric=list("p"=2,"w"=1)
#    par.metric$fdata1<-fdataobj
#    #mdist=metric(fdataobj,...)
#    mdist=do.call(metric,par.metric)
#    par.metric2<-par.metric
#    par.metric2$fdata2<-ncl
#    mdist2=do.call(metric,par.metric2)
#    out1 = list()
#    out1$fdataobj<-fdataobj
#    out1$centers = ncl
#    out1$lcenters <- NULL
#    ngroups=nrow(ncl)
#    ncl = nrow(ncl)
#    out1$d = rbind(mdist,t(mdist2))
#    class(out1) = "kmeans.fd"
# }
#  ngroups=nrow(out1$centers[["data"]])
#  a=0;aa<-i<-1
# # C <- match.call()
# # mf <- match.call(expand.dots = FALSE)
#  same_centers=FALSE
# #while ((i<max.iter) && (a!=aa)) {
# if (is.null(colnames(out1$d))) 
#   colnames(out1$d)<-1:NCOL(out1$d)
#  cnames<-colnames(out1$d)
# while ((i<max.iter) && (!same_centers)) {
#   iterar<-FALSE
# # print(i)
# # print("1. assign")
#   out3=kmeans.assig.groups(out1,draw=draw)
#   names(out3$cluster) <- cnames
#   # print("tala out 3")
#   # print(table(out3$cluster))
#   # print("tala out 4")
#   tab <- table(out3$cluster)
#   imin <- which.min(tab)[1]
# print(tab)  
# # print(imin)
#   
#   if (cluster.size > tab[imin] ) {
# # print("entra tab")
#     # ngroups <- ngroups-1
#     if (fda.usc:::par.fda.usc$warning) warning(paste0(" One of the clusters only has ",tab[imin]
#     , " curves and the minimum cluster size is ", cluster.size
#     ,".\n The cluster is completed with the closest curves of the other clusters."))
#     # names(sort(aa)[1:5])
#     iclust <- out3$cluster == imin
#  #print(   range(out1$d[imin,]))
#  #print(sum(iclust))
#     out1$d[nr+imin,iclust] <- 0
#     out1$d[nr+imin,!iclust] <- out1$d[nr+imin,!iclust]*.99
#   #  print(dim(out1$d))
#   #  print(dim(out1$z.dist))
#     #  print(range(out1$d[imin,]))
#     #print(range(out1$z.dist))
#     #out1$z.dist<-out1$d
#     iterar <-TRUE
#     #i=i+1
#  }
# else{
#      # print("2. update")
#      out2=kmeans.centers.update(out1, group=out3$cluster
#                              , dfunc=dfunc, draw=draw
#                              , par.dfunc=par.dfunc
#                              #,cluster.size=cluster.size
#                              ,...)
#     if (!iterar){
#     same_centers <- out2$centers$data == out3$centers$data
#     out1$centers <- out2$centers
#     i=i+1
#     }
#   }
# cat("iterations: ",i)
#   }
# 
# out<-list("cluster"=out2$cluster,"centers"=out2$centers
#           ,mdist=mdist,metric = metric)
# return(out)
# }
# 
# 
