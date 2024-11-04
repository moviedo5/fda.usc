#' @name depth.mfdata
#' 
#' @title Provides the depth measure for a list of p--functional data objects
#' 
#' @description This function computes the depth measure for a list of p--functional data
#' objects.  The procedure extends the Fraiman and Muniz (FM), modal, and
#' random project depth functions from 1 functional dataset to p functional
#' datasets.
#' @details
#' \itemize{ 
#' \item \code{\link{depth.FMp}}, this procedure suposes that each
#' curve of the mfdataobj have the same support [0,T] (same argvals and
#' rangeval).  The FMp depth is defined as: \eqn{FM_i^p
#' =\int_{0}^{T}Z_i^p(t)dt}{FM_i^p =\int Z_i^p(t)dt} where
#' \eqn{Z_i^p(t)}{Z_i^p(t)} is a \eqn{p}--variate depth of the vector
#' \eqn{(x_i^1(t),\ldots,x_i^p(t))}{(x_i^1(t),...,x_i^p(t))} w.r.t. the sample
#' at \eqn{t}. % This happens, for instance, when using the curves and its
#' derivatives. In this case,note solo un dato funcional se reduce
#' depth.FM=depth.FM1 
#' \item The \code{\link{depth.RPp}} function calculates the depth in two
#' steps. It builds random projections for the each curve of the \code{mfdata}
#' w.r.t. each curve of the \code{mfdataref} object. Then it applyes a
#' multivariate depth function specified in \code{dfunc} argument to the set of
#' random projections.  This procedure is a generalization of Random Projection
#' with derivatives (RPD) implemented in \code{\link{depth.RPD}} function.
#' Now, the procedure computes a p-variate depth with the projections using the
#' \eqn{p} functional dataset.
#' \item The modal depth \code{\link{depth.modep}} function calculates the
#' depth in three steps.  First, the function calculates a suitable metrics or
#' semi--metrics \eqn{m_1+\cdots+m_p}{m1,...,mp} for each curve of the
#' \code{mfdata} w.r.t. each curve in the \code{mfdataref} object using the
#' \code{metric} and \code{par.metric} arguments, see \code{\link{metric.lp}}
#' or \code{\link{semimetric.NPFDA}} for more details.  Second, the function
#' uses the \eqn{p}--dimensional metrics to construct a new metric, specified
#' in \code{method} argument, by default if \code{method="euclidean"}, i.e.
#' \eqn{m:=\sqrt{m_1^2+\cdots+m_p^2}}{m:=\sqrt(m1^2+...+mp^2)}. Finally, the
#' empirical \emph{h}--depth is computed as:
#' \deqn{\hat{f}_h(x_0)=N^{-1}\sum_{i=1}^{N}{K(m/h)}}{\hat{f}(x)=1/N\sum
#' K(m/h)} where \eqn{x} is dataset with p observed fucntional data, \eqn{m} is
#' a suitable metric or semi--metric, \eqn{K(t)} is an asymmetric kernel
#' function and \code{h} is the bandwidth parameter. 
#' }
#' 
#' @aliases depth.FMp depth.modep depth.RPp
#' @param mfdata A list of new curves (list of fdata ojects) to evaluate the
#' depth.
#' @param mfdataref A set of reference curves (list of fdata ojects) w.r.t. the
#' depth of mfdata is computed.
#' @param trim The alpha of the trimming.
#' @param dfunc Type of multivariate depth (of order p) function used in
#' Framiman and Muniz depth, \code{depth.FMp} or in Random Projection
#' depth,\code{depth.FMp}: : \itemize{ \item The \code{\link{mdepth.SD}}
#' function provides the simplicial depth measure for bivariate data.  \item
#' The \code{\link{mdepth.LD}} function provides the Likelihood depth measure
#' based on Nadaraya--Watson estimator of empirical density function.  \item
#' The \code{\link{mdepth.HS}} function implements a half-space depth measure
#' based on random projections.  \item The \code{\link{mdepth.TD}} function
#' implements a Tukey depth measure.  \item The
#' \code{\link{mdepth.MhD}}function implements a Mahalanobis depth measure.
#' \item The \code{\link{mdepth.RP}} function provides the depth measure using
#' random projections for multivariate data.  }
#' @param par.dfunc list of parameters for the \code{dfunc} depth function, see
#' \code{\link{Depth.Multivariate}}.
#' @param nproj The number projection.
#' @param proj if is a character: create the random projection using a
#' covariance matrix by process indicated in the argument (by default, proj=1,
#' sigma=diag(ncol(fdataobj))), else if is a matrix of random projection
#' provided by the user.
#' @param h Bandwidth, \code{h>0}. Default argument values are provided as the
#' 15\%--quantile of the distance between \code{fdataobj} and \code{fdataori}.
#' @param metric Metric or semi--metric function used for compute the distance
#' between each element in \code{ldata} w.r.t. \code{ldataref}, by default
#' \code{\link{metric.lp}}.
#' @param par.metric list of parameters for the metric function.
#' @param method Type of the distance measure (by default \code{euclidean}) to
#' compute the metric between the p-distance matrix computed from the p
#' functional data elements.
#' @param scale =TRUE, scale the depth.
#' @param draw =TRUE, draw the curves, the sample median and trimmed mean.
#' @param ask logical. If \code{TRUE} (and the R session is interactive) the
#' user is asked for input, before a new figure is drawn.
#' @param \dots Further arguments passed to or from other methods.
#' @return
#' \itemize{
#' \item \code{lmed}: Index deepest element \code{median}.
#' \item \code{ltrim}: Index of curves with trimmed mean \code{mtrim}.
#' \item \code{dep}: Depth of each curve of \code{fdataobj} w.r.t. \code{fdataori}.
#' \item \code{dfunc}: Second depth function used as multivariate depth, see details section.
#' \item \code{par.dfunc}: List of parameters for the \code{dfunc} depth function.
#' \item \code{proj}: The projection value of each point on the curves.
#' \item \code{dist}: Distance matrix between curves or functional data.
#' }
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' @seealso See Also as \code{\link{Descriptive}}.
#' @references Cuevas, A., Febrero-Bande, M. and Fraiman, R. (2007).
#' \emph{Robust estimation and classification for functional data via
#' projection-based depth notions.} Computational Statistics 22, 3, 481-496.
#' 
#' @keywords descriptive
#' @examples
#' \dontrun{
#' data(tecator)
#' xx<-tecator$absorp
#' xx1<-fdata.deriv(xx,1)
#' lx<-list(xx=xx,xx=xx1)
#' # Fraiman-Muniz Depth
#' par.df<-list(scale =TRUE)
#' out.FM1p=depth.FMp(lx,trim=0.1,draw=TRUE, par.dfunc = par.df)
#' out.FM2p=depth.FMp(lx,trim=0.1,dfunc="mdepth.LD",
#' par.dfunc = par.df, draw=TRUE)
#' 
#' # Random Project Depth
#' out.RP1p=depth.RPp(lx,trim=0.1,dfunc="mdepth.TD",
#' draw=TRUE,par.dfunc = par.df)
#' out.RP2p=depth.RPp(lx,trim=0.1,dfunc="mdepth.LD",
#' draw=TRUE,par.dfunc = par.df)
#' 
#' #Modal Depth
#' out.mode1p=depth.modep(lx,trim=0.1,draw=T,scale=T)
#' out.mode2p=depth.modep(lx,trim=0.1,method="manhattan",
#' draw=T,scale=T)
#' 
#' par(mfrow=c(2,3))
#' plot(out.FM1p$dep,out.FM2p$dep)
#' plot(out.RP1p$dep,out.RP2p$dep)
#' plot(out.mode1p$dep,out.mode2p$dep)
#' plot(out.FM1p$dep,out.RP1p$dep)
#' plot(out.RP1p$dep,out.mode1p$dep)
#' plot(out.FM1p$dep,out.mode1p$dep)
#' }
#' 
#' @rdname depth.mfdata
#' @export
depth.modep=function(mfdata,mfdataref=mfdata,h=NULL,metric,
                     par.metric=list(),method="euclidean",
                     scale=FALSE,trim=0.25,draw=FALSE,ask=FALSE) 
{
  equal<-identical(mfdata,mfdataref)
  if (inherits(mfdata,"list")){
    lenl <- length(mfdata)
    lenl2 <- length(mfdataref)
    if (lenl!=lenl2) stop("Length of mfdata!= length of mfdataref")
    m0 <- nrow(mfdata[[1]])
    if (is.null(rownames(mfdata[[1]]$data))) 
      rownames(mfdata[[1]]$data) <- 1:m0
    nms <- rownames(mfdata[[1]]$data)
    nas <- NULL
    for (i in 1:lenl) {
      nas <- c(nas, na.action(na.omit(mfdata[[i]])))
    }
    nas <- unique(nas)
    nullans <- !is.null(nas)
    nam1<-names(mfdata)
    if (missing(metric)){
      if (is.fdata(mfdata[[nam1[1]]])) metric<-rep("metric.lp",len=lenl)
      else    metric<-rep("metric.dist",len=lenl)   
    }
    mdist2<-metric.ldata(ldata1=mfdata,ldata2=mfdataref,metric=metric,par.metric=par.metric,method=method)
    if (equal) mdist<-mdist2
    else mdist<-metric.ldata(ldata1=mfdataref,ldata2=mfdataref,metric=metric,par.metric=par.metric,method=method)
    m<-m0
    n<-nrow(mfdata[[1]])
  }
  else stop("Error in mfdata argument")
  #attr(mdist, "par.metric") <-attr(atr, "par.metric")                                  
  #attr(mdist, "call") <-attr(atr, "call")                                  
  #attr(mdist, "method") <-method                                    
  if (is.null(n) && is.null(m)) stop("ERROR IN THE DATA DIMENSIONS")
  #if (is.null(m) && is.null(m2)) stop("ERROR IN THE DATA DIMENSIONS")
  #h="h=0.15"
  if (is.null(h))   {
    h<-0.1
    hq2=quantile(mdist,probs=h,na.rm=TRUE)
    #print("es nulo h")  
  }
  else {
    #cat("no es nulo h ",h,"\n")    
    if (is.numeric(h))    hq2<-h  
    else hq2=quantile(mdist,probs=as.numeric(h),na.rm=TRUE)
  }
  #print(hq2)
  class(mdist2) <- class(mdist)<-c("matrix","fdist")
  dep <- Ker.norm(mdist2/hq2)    #### en (nxm) matrix (new X ref)
  dep <- apply(dep,1,sum,na.rm=TRUE)                                      
  if (scale)   {
    dep2<-Ker.norm(mdist/hq2)
    dep2<-apply(dep2,1,sum,na.rm=TRUE)
    # mn<-min(dep2,na.rm=TRUE)
    mx<-max(dep2,na.rm=TRUE)
    #   scl<-mx-mn
    #   ans=as.vector(scale(ans,center=mn,scale=scl))
    #   ans2=as.vector(scale(ans2,center=mn,scale=scl))
    dep=as.vector(dep/mx)   
  }
  if (nullans) {
    ans1<-rep(NA,len=m0)
    ans1[-nas] <-dep
    dep<-ans1      
  }
  names(dep)=nms
  k = which.max(dep) 
  nl=length(trim)
  tr<-paste("modep.tr",round(trim*100,2),"\u0025",sep="")
  lista=vector("list",nl)
  med=vector("list",lenl)
  mtrim=vector("list",lenl)
  for (ik in 1:lenl){
    med[[ik]]=mfdata[[ik]][k]
    mtrim[[ik]]=fdata(matrix(NA,ncol=length(mfdata[[ik]]$argvals),nrow=length(trim)),mfdata[[ik]]$argvals,mfdata[[ik]]$rangeval,mfdata[[ik]]$names)
    for (j in 1:nl){
      lista[[j]]=which(dep>=quantile(dep,probs=trim[j],na.rm=TRUE))
      if (length(lista[[j]])==1) {
        mtrim[[ik]]$data[j,]<-mfdata[[ik]][lista[[j]]]$data
        if (draw) {draw=FALSE;warning("Too few curves in mtrim. The plot is not drawn")}
      }
      else   mtrim[[ik]]$data[j,]=apply(mfdata[[ik]]$data[lista[[j]],,drop=FALSE],2,mean,na.rm=TRUE)
    }
    rownames(med[[ik]]$data)="modep.med"
    rownames(mtrim[[ik]]$data)=tr
  }
  if (nl>1) names(lista)=paste0("tr",trim*100)	
  if (draw) {
    mf=5
    if (lenl>4) ask=TRUE
    if (ask) {par(mfrow = c(1, 1))
      dev.interactive()
      oask <- devAskNewPage(TRUE)
      on.exit(devAskNewPage(oask))}
    else{    mf<-switch(lenl,
                        "1"={c(1,1)},
                        "2"={c(1,2)},
                        "3"={c(1,3)},
                        "4"={c(2,2)})            
    par(mfrow =mf)                    
    }
    ans <- dep
    ind1 <- !is.nan(ans)
    ans[is.nan(ans)] = NA
    color=colorRampPalette(c("red","blue"))(nl+1)
    cgray = 1 - (ans - min(ans, na.rm = TRUE))/(max(ans,
                                                    na.rm = TRUE) - min(ans, na.rm = TRUE))
    for (ik in 1:lenl) {
      plot(mfdata[[ik]][ind1, ], col =gray(cgray[ind1]),lty=1, main = paste(nam1[ik],
                                                                            " modep Depth",sep=""))
      lines(mtrim[[ik]],lwd=3,col=color[-1],lty=1)
      lines(med[[ik]],col=color[1],lwd=3)
      legend("topleft",legend=c("Median",tr),lwd=3,col=color,box.col=0)
    }
  }
  
  out<- list(median = med, lmed = k, mtrim = mtrim,
             ltrim = if (nl==1) unlist(lista) else lista, dep = dep,
             "metric"=metric,"par.metric"=par.metric,"mdist"=mdist,"hq"=hq2)
  class(out)="depth"
  return(invisible(out))
} 

#' @rdname depth.mfdata
#' @export
depth.RPp<-function (mfdata,mfdataref=mfdata, nproj = 50, proj="vexponential",trim = 0.25,
dfunc="mdepth.TD",par.dfunc=list(scale=TRUE),draw = FALSE,ask=FALSE){  
verbose<-TRUE
if (!is(mfdata,"list"))     stop("mfdata input must be a list of fdata objects")
if (!is(mfdataref,"list")) stop("mfdataref input must be a list of fdata objects")
 lenl <- length(mfdata)
 lenl2 <- length(mfdataref)
 n <- nrow(mfdata[[1]])
 m <- nrow(mfdataref[[1]]) 
 # mdist2 <- matrix(0,n,n)
 # amdist <- array(NA,dim=c(n,m,lenl))
 # mdist <- list()
 nam1 <- names(mfdata)
 nam2 <- names(mfdataref) 
 if (is.null(nam1)) {names(mfdata)<-nam1<-paste("var",1:lenl,sep="")}
 if (is.null(nam2)) {names(mfdataref)<-nam2<-paste("var",1:lenl2,sep="")} 
# depth<-paste("mdepth.",dfunc,sep="")
# if (missing(dfunc))  dfunc<-depth.mode#dfunc<-rep("depth.mode",len=lenl)
# else if (length(dfunc)==1)   dfunc<-rep("depth.mode",len=lenl)
  if (is.fdata(mfdata[[nam1[1]]]))  {
    fdataobj<-mfdata[[nam1[1]]]
    fdataori<-mfdataref[[nam1[1]]]}
  else stop("No fdata object in the mfdata argument")
  data <- fdataobj[["data"]]
  data2 <- fdataori[["data"]]    
  names1 <- names2 <- names <- fdataobj[["names"]]
  names1$main <- "depth.RPD median"
  names2$main <- paste("RPD trim ", trim * 100, "%", sep = "")
  n <- nrow(data)
  m <- ncol(data)
  m2 <- ncol(data2)
  n2 <- nrow(data2) 
  if (is.null(m) && is.null(m2)) stop("ERROR IN THE DATA DIMENSIONS")
  if (is.null(n) || is.null(m))       stop("Input must be a matrix")
  tt = fdataobj[["argvals"]]
  rtt <- fdataobj[["rangeval"]]
  newfunc = array(NA, dim = c(n, m, lenl))
  newfunc2 = array(NA, dim = c(n2, m2, lenl))    
##################
  for (ider in 1:lenl) {
          newfunc[, , ider] = mfdata[[nam1[ider]]]$data     
          newfunc2[, , ider] = mfdataref[[nam1[ider]]]$data 
  }
  dep = rep(0, n)
#  dep2 = rep(0, n2)    
  vproject = matrix(0, nrow = n, ncol = lenl)
  vproject2 = matrix(0, nrow = n2, ncol = lenl)
#   modulo = function(z) {        sqrt(sum(z^2))    }
#  z = matrix(rnorm(m * nproj), nrow = nproj, ncol = m)
#  modu = apply(z, 1, modulo)
#  z = z/modu
  if (is.fdata(proj)) {
     if (any(fdataobj$argvals!=proj$argvals) || ncol(fdataobj)!=ncol(proj)) stop("Error en proj dimension")
     z<-proj
     nproj<-nrow(z)
    }
  else {	 z<-rproc2fdata(nproj,tt,sigma=proj,norm=TRUE)	}
  if (verbose) 
  pb = txtProgressBar(min = 0, max = nproj, style = 3)
  for (j in 1:nproj) {
      if (verbose)         setTxtProgressBar(pb, j - 0.5)
        for (ider in 1:lenl) {
            matriz = newfunc[, , ider]
            vproject[, ider] = matriz %*%z$data[j, ]
            vproject2[, ider] = newfunc2[, , ider] %*%z$data[j, ]              
        }
#        resul = dfunc2(vproject, ...)       
#        par.dfunc = list()
     if (substr(dfunc,1,1)=="m") {
        par.dfunc$x<- vproject       
        par.dfunc$xx <- vproject2
        }
     else{ 
        par.dfunc$fdataobj <- fdata(vproject)        
        par.dfunc$fdataori <- fdata(vproject2)
        }
#        par.dfunc$scale<-TRUE
        resul = do.call(dfunc, par.dfunc)
        if (dfunc=="depth.RP" || dfunc=="mdepth.HS" || dfunc=="mdepth.RP") par.dfunc$proj<-resul$proj        
        dep = dep + resul$dep
     if (verbose)         setTxtProgressBar(pb, j)
  }                                                               
  if (verbose)     close(pb)
    dep = dep/nproj
    k = which.max(dep)
    med = data[k,,drop=FALSE ]
    lista = which(dep >= quantile(dep, probs = trim,na.rm=TRUE))
    mtrim=apply(data[lista,,drop=FALSE],2,mean,na.rm=TRUE)
#else mtrim<-data     
    tr <- paste("RPD.tr", trim * 100, "%", sep = "")
    med <- fdata(med, tt, rtt, names1)
    mtrim <- fdata(mtrim, tt, rtt, names2)
    rownames(med$data) <- "RPD.med"
    rownames(mtrim$data) <- tr
if (draw){
  mf=5
  if (lenl>4) ask=TRUE
  if (ask) {par(mfrow = c(1, 1))
            dev.interactive()
            oask <- devAskNewPage(TRUE)
            on.exit(devAskNewPage(oask))}
   else{    mf<-switch(lenl,
   "1"={c(1,1)},
   "2"={c(1,2)},
   "3"={c(1,3)},
   "4"={c(2,2)})            
            par(mfrow =mf)                    }          
 for (idat in 1:lenl) {
   data<-mfdata[[idat]]$data
   med<-data[k,,drop=FALSE]      
  mtrim=apply(data[lista,,drop=FALSE],2,mean,na.rm=TRUE)
# mtrim=data[lista,]
   med<-fdata(med,tt,rtt,names1)
   mtrim<-fdata(mtrim,tt,rtt,names2)
   rownames(med$data)<-"RPd.med"
   rownames(mtrim$data)<-tr
   ind1<-!is.nan(dep)
   dep[is.nan(dep)]=NA
   cgray=1-(dep-min(dep,na.rm=TRUE))/(max(dep,na.rm=TRUE)-min(dep,na.rm=TRUE))
   plot(mfdata[[idat]][ind1, ], col =  gray(cgray[ind1]),lty=1, main = paste(nam1[idat],"curves with RPp Depth",sep=""))
   lines(mtrim,lwd=2,col="yellow")
   lines(med,col="red",lwd=2)
   legend("topleft",legend=c(tr,"Median"),lwd=2,col=c("yellow","red"),box.col=0)
 }
}
    out<- list(median = med, lmed = k, mtrim = mtrim,
               ltrim = lista, dep = dep,proj = z,dfunc=dfunc,par.dfunc=par.dfunc)
    out$trim= out$trim
    out$fdataobj=fdataobj
    out$fdataori=fdataori
    class(out)<-"depth"
return(invisible(out))
}  

#' @rdname depth.mfdata
#' @export 
depth.FMp=function(mfdata,mfdataref=mfdata,trim=0.25,dfunc="mdepth.MhD",
                   par.dfunc=list(scale=FALSE),draw=FALSE,ask=FALSE,...){
  #   print("FMp")
  #if (!is.list(mfdata)) stop("mfdata1 must be a list")
  #if (!is.list(mfdataref)) stop("mfdata1 must be a list")
  nam1<-names(mfdata)
  nam2<-names(mfdata)
  len1<-length(mfdata)
  len2<-length(mfdataref) 
  if (len1!=len2) stop("mfdata and mfdataref have no the same length")
  m0<-nrow(mfdata[[1]])
  if (is.null(rownames(mfdata[[1]]$data)))  rownames(mfdata[[1]]$data)<-1:m0
  nms<-rownames(mfdata[[1]]$data)
  nas<-NULL
  for (i in 1:len1) {
    nas<-c(nas,na.action(na.omit(mfdata[[i]])))
  } 
  nas<-  (nas)
  #nullans<-!is.null(nas) 
  #for (i in 1:len1) {
    #  if (nullans) mfdata[[i]]<-mfdata[[i]][-nas]
    #  if (nullans)  mfdataref[[i]]<-mfdataref[[i]][-nas]  
  #}
  #comprobar is.fdatalist
  if (is.null(nam1)) nam1<-paste("var",1:len1,sep="")
  if (is.null(nam2)) nam2<-nam1
  fdataobj<-mfdata[[1]]
  fdataori<-mfdataref[[1]] 
  if (is.null(rownames(fdataobj$data)))  rownames(fdataobj$data)<-seq_len(nrow(fdataobj$data))
  nms<-rownames(fdataobj$data)
  n<-nrow(fdataobj)
  m<-nrow(fdataori)
  p<-ncol(fdataori)
  tt<-fdataori$argvals
  rtt<-fdataori$rangeval
  xref<-matrix(NA,m,len1)
  x0<-matrix(NA,n,len1)
  #depth<-paste("mdepth.",dfunc,sep="")
  d<-array(NA,dim=c(n,p)) 
  for (idat in 1:len1) {
    if (any(mfdata[[idat]]$argvals!=tt))  stop("Incorrect argvals in the fdata objects")}
  # verificar dimensiones mfdata y mismos argvals/rangeval
  for (iti in 1:p) {
    #  mfdata[[idat]]
    for (idat in 1:len1) { 
      x0[,idat]<-mfdata[[idat]]$data[,iti]  
      xref[,idat]<-mfdataref[[idat]]$data[,iti]    
    } 
    par.dfunc$x<-x0
    par.dfunc$xx<-xref
    d[,iti]<-do.call(dfunc,par.dfunc)$dep 
  }   
  #print(d)         
  data<-fdataobj[["data"]]
  data2<-fdataori[["data"]]
  n<-nrow(data)
  m<-ncol(data)
  m2<-ncol(data2)
  n2<-nrow(data2)
  names1<-names2<-names<-mfdata[[1]][["names"]]
  names1$main<-"depth.FM median"
  
  tr<-paste("FM.tr",trim*100,"\u0025",sep="")
  names2$main<-tr
  dtt <- diff(tt)
  eps <- as.double(.Machine[[1]] * 100)
  inf <- dtt - eps
  sup <- dtt + eps
# Creo que no se usa equi
#  if (all(dtt > inf) && all(dtt < sup)) {
#    equi = TRUE
#  }
#  else equi = FALSE
  #d2<-int.simpson(fdata(d,tt,rtt),equi=equi)/n
  ans<-apply(d,1,mean)[1:n]
  #if (nullans) {
  #  ans1<-rep(NA,len=m0)
  #ans1[-nas] <-ans 
  # ans<-ans1    
  # }       
  names(ans)<-nms    
  k=which.max(ans)    
  lista=which(ans>=quantile(ans,probs=trim,na.rm=TRUE))
  if (draw){
    mf=5
    if (len1>4) ask=TRUE
    if (ask) {par(mfrow = c(1, 1))
      dev.interactive()
      oask <- devAskNewPage(TRUE)
      on.exit(devAskNewPage(oask))}
    else{    mf<-switch(len1,
                        "1"={c(1,1)},
                        "2"={c(1,2)},
                        "3"={c(1,3)},
                        "4"={c(2,2)})            
    par(mfrow =mf)                    }          
    for (idat in 1:len1) {
      data<-mfdata[[idat]]$data
      med<-data[k,]    
      
      mtrim=apply(data[lista,,drop=FALSE],2,mean,na.rm=TRUE)
      med<-fdata(med,tt,rtt,names1)
      mtrim<-fdata(mtrim,tt,rtt,names2)
      rownames(med$data)<-"FMp.med"
      rownames(mtrim$data)<-tr
      ind1<-!is.nan(ans)
      ans[is.nan(ans)]=NA
      cgray=1-(ans-min(ans,na.rm=TRUE))/(max(ans,na.rm=TRUE)-min(ans,na.rm=TRUE))
      plot(mfdata[[idat]][ind1, ], col =  gray(cgray[ind1]),lty=1, 
                        main = paste(nam1[idat]," FMp Depth",sep=""))
      lines(mtrim,lwd=2,col="yellow")
      lines(med,col="red",lwd=2)
      legend("topleft",legend=c(tr,"Median"),lwd=2,col=c("yellow","red"),box.col=0)
    }
  }
  out <- list("lmed"=k,"ltrim"=lista,"dep"=ans,"par.dfunc"=par.dfunc)
  class(out) <- "depth"
  return(invisible(out))
}