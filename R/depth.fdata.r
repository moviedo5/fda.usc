#' @name depth.fdata
#' 
#' @title Computation of depth measures for functional data
#' 
#' @description Several depth measures can be computed for functional data for descriptive
#' or classification purposes. 
#' 
#' @details Type of depth functions: Fraiman and Muniz (FM)
#' depth, modal depth, random Tukey (RT), random projection (RP) depth and
#' double random projection depth (RPD).
#' \itemize{ 
#' \item \code{\link{depth.FM}} computes the integration of an univariate depth
#' along the axis x (see Fraiman and Muniz 2001). It is also known as
#' Integrated Depth.
#' 
#' \item \code{\link{depth.mode}} implements the modal depth (see Cuevas et al
#' 2007).
#' 
#' \item \code{\link{depth.RT}} implements the Random Tukey depth (see
#' Cuesta--Albertos and Nieto--Reyes 2008).
#' 
#' \item \code{\link{depth.RP}} computes the Random Projection depth (see
#' Cuevas et al. 2007).
#' 
#' \item \code{\link{depth.RPD}} implements a depth measure based on random
#' projections possibly using several derivatives (see Cuevas et al. 2007).
#' 
#' \item \code{\link{depth.FSD}} computes the Functional Spatial Depth (see
#' Sguera et al. 2014).
#' 
#' \item \code{\link{depth.KFSD}} implements the Kernelized Functional Spatial
#' Depth (see Sguera et al. 2014).  
#' } 
#' \itemize{ 
#' \item The \code{\link{depth.mode}} function calculates the depth of a datum
#' accounting the number of curves in its neighbourhood. By default, the
#' distance is calculated using \code{\link{metric.lp}} function although any
#' other distance could be employed through argument \code{metric} (with the
#' general pattern \code{USER.DIST(fdataobj,fdataori)}).
#' \item The \code{\link{depth.RP}} function summarizes the random projections
#' through averages whereas the \code{\link{depth.RT}} function uses the
#' minimum of all projections.
#' \item The \code{\link{depth.RPD}} function involves the original
#' trajectories and the derivatives of each curve in two steps. It builds
#' random projections for the function and their derivatives (indicated in the
#' parameter \code{deriv}) and then applies a depth function (by default
#' \code{\link{depth.mode}}) to this set of random projections (by default the
#' Tukey one).
#' \item The \code{\link{depth.FSD}} and \code{\link{depth.KFSD}} are the
#' implementations of the default versions of the functional spatial depths
#' proposed in Sguera et al 2014. At this moment, it is not possible to change
#' the kernel in the second one.#' 
#' }
#' 
#' @aliases Depth depth.FM depth.mode depth.RP depth.RPD depth.RT depth.FSD
#' depth.KFSD
#' @param fdataobj The set of new curves to evaluate the depth.
#' \code{\link{fdata}} class object.
#' @param fdataori The set of reference curves respect to which the depth is
#' computed.  \code{\link{fdata}} class object.
#' @param trim The alpha of the trimming.
#' @param nproj The number of projections. Ignored if a \code{fdata} class
#' object is provided in \code{proj}
#' @param proj if a \code{fdata} class, projections provided by the user.
#' Otherwise, it is the \code{sigma} parameter of \code{\link{rproc2fdata}}
#' function.
#' @param dfunc type of univariate depth function used inside depth function:
#' "FM1" refers to the original Fraiman and Muniz univariate depth (default),
#' "TD1" Tukey (Halfspace),"Liu1" for simplical depth, "LD1" for Likelihood
#' depth and "MhD1" for Mahalanobis 1D depth. Also, any user function
#' fulfilling the following pattern \code{FUN.USER(x,xx,...)} and returning a
#' \code{dep} component can be included.f
#' @param par.dfunc List of parameters for \emph{dfunc}.
#' @param dfunc2 Multivariate depth function (second step depth function) in
#' RPD depth, by default \code{\link{mdepth.LD}}. Any user function with the
#' pattern \code{FUN.USER(x,xx,...)} can be employed.
#' @param deriv Number of derivatives described in integer vector \code{deriv}.
#' \code{=0} means no derivative.
#' @param method Type of derivative method. See \code{\link{fdata.deriv}} for
#' more details.
#' @param h Bandwidth parameter. 
#'  \itemize{ 
#' \item If \code{h} is a numerical value, the procedure considers the argument
#' value as the bandwidth.  
#' \item If is \code{NULL} (by default) the bandwidth is provided as the 
#'  15\%--quantile of the distance among curves of
#' \code{fdataori}.
#' \item If \code{h} is a character string (like \code{"0.15"}), the procedure
#' reads the numeric value and consider it as the quantile of the distance in
#' \code{fdataori} (as in the second case).  
#' }
#' @param metric Metric function, by default \code{\link{metric.lp}}. Distance
#' matrix between \code{fdataobj} and \code{fdataori}.
#' @param scale =TRUE, the depth is scaled respect to depths in
#' \code{fdataori}.
#' @param xeps Accuracy. The left limit of the empirical distribution function.
#' @param draw =TRUE, draw the curves, the sample median and trimmed mean.
#' @param \dots Further arguments passed to or from other methods. For
#' \code{depth.mode} parameters for \code{metric}. For random projection
#' depths, parameters to be included in \code{rproc2fdata} not included before.
#' @return  Return a list with:
#' \itemize{
#' \item {median}{ Deepest curve.} 
#' \item {lmed}{ Index deepest element \code{median}.}
#' \item {mtrim}{ \code{fdata} class object with the average from the \code{(1-trim)\%} deepest curves. }
#' \item {ltrim}{ Indexes of curves that conform the trimmed mean \code{mtrim}. }
#' \item {dep}{ Depth of each curve of fdataobj w.r.t. fdataori.}
#' \item {dep.ori}{ Depth of each curve of fdataori w.r.t. fdataori.}
#' \item {proj}{ The projection value of each point on the curves. } 
#' \item {dist}{ Distance matrix between curves or functional data.}
#' }
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' @seealso See Also as \code{\link{Descriptive}}.
#' @references Cuevas, A., Febrero-Bande, M., Fraiman, R. (2007). Robust
#' estimation and classification for functional data via projection-based depth
#' notions. \emph{Computational Statistics} 22, 3, 481-496.
#' 
#' Fraiman R, Muniz G. (2001). Trimmed means for functional data. \emph{Test}
#' 10: 419-440.
#' 
#' Cuesta--Albertos, JA, Nieto--Reyes, A. (2008) The Random Tukey Depth.
#' \emph{Computational Statistics and Data Analysis} Vol. 52, Issue 11,
#' 4979-4988.
#' 
#' Febrero-Bande, M, Oviedo de la Fuente, M. (2012).  Statistical Computing in
#' Functional Data Analysis: The R Package fda.usc. \emph{Journal of
#' Statistical Software}, 51(4), 1-28. \url{https://www.jstatsoft.org/v51/i04/}
#' 
#' Sguera C, Galeano P, Lillo R (2014). Spatial depth based classification for
#' functional data. \emph{TEST} 23(4):725--750.
#' @keywords descriptive
#' @examples 
#' \dontrun{
#' #Ex: CanadianWeather data
#' tt=1:365
#' fdataobj<-fdata(t(CanadianWeather$dailyAv[,,1]),tt)
#' # Fraiman-Muniz Depth
#' out.FM=depth.FM(fdataobj,trim=0.1,draw=TRUE)
#' #Modal Depth
#' out.mode=depth.mode(fdataobj,trim=0.1,draw=TRUE)
#' out.RP=depth.RP(fdataobj,trim=0.1,draw=TRUE)
#' out.RT=depth.RT(fdataobj,trim=0.1,draw=TRUE)
#' out.FSD=depth.FSD(fdataobj,trim=0.1,draw=TRUE)
#' out.KFSD=depth.KFSD(fdataobj,trim=0.1,draw=TRUE)
#' ## Double Random Projections
#' out.RPD=depth.RPD(fdataobj,deriv=c(0,1),dfunc2=mdepth.LD,
#' trim=0.1,draw=TRUE)
#' out<-c(out.FM$mtrim,out.mode$mtrim,out.RP$mtrim,out.RPD$mtrim)
#' plot(fdataobj,col="grey")
#' lines(out)
#' cdep<-cbind(out.FM$dep,out.mode$dep,out.RP$dep,out.RT$dep,out.FSD$dep,out.KFSD$dep)
#' colnames(cdep)<-c("FM","mode","RP","RT","FSD","KFSD")
#' pairs(cdep)
#' round(cor(cdep),2)
#' }
#' 
#' @rdname depth.fdata
#' @export depth.mode
depth.mode=function(fdataobj,fdataori=fdataobj,trim=0.25,
                    metric=metric.lp,h=NULL,scale=FALSE,
                    draw=FALSE,...){    
#if (is.fdata(fdataobj)) {
# fdat<-TRUE
# nas<-is.na.fdata(fdataobj)
#if (any(nas))  {
#   fdataobj$data<-fdataobj$data[!nas,]
#   cat("Warning: ",sum(nas)," curves with NA are not used in the calculations \n")
#   }
 if (is.null(rownames(fdataobj$data)))  rownames(fdataobj$data)<-1:nrow(fdataobj$data)
 nms<-rownames(fdataobj$data)
 m0<-nrow(fdataobj)
 fdataobj2=fdataobj
 fdataobj<-na.omit.fdata(fdataobj)
 fdataori<-na.omit.fdata(fdataori) 
 nas<-na.action(fdataobj)
 nullans<-!is.null(nas) 
# data<-fdataobj[["data"]]
# data2<-fdataori[["data"]]
 names1<-names2<-fdataobj[["names"]]
 names1$main<-"depth.mode median"
 names2$main<-paste("depth.mode trim ",trim*100,"\u0025",sep="")
 tt=fdataobj[["argvals"]]
 rtt<-fdataobj[["rangeval"]]
#print("is fdata") 
#}
#else { stop("no fdata class object")
#        data<-fdataobj
#        data2<-fdataori
#        fdat<-FALSE   
#  }
n<-nrow(fdataobj)
m<-ncol(fdataobj)
m2<-ncol(fdataori)
n2<-nrow(fdataori)
if (is.null(n) && is.null(m)) stop("ERROR IN THE DATA DIMENSIONS")
if (is.null(m) && is.null(m2)) stop("ERROR IN THE DATA DIMENSIONS")
if (is.matrix(metric)) mdist=metric                    
else  mdist=metric(fdataori,fdataori,...)
class(mdist) <- "matrix"

if (n==n2 & m==m2) {
  equal<-all(fdataobj$data==fdataori$data)
  if (equal) mdist2<-mdist
  else mdist2<-metric(fdataobj,fdataori,...)}
else  mdist2<-metric(fdataobj,fdataori,...)
if (is.null(h))   {
  h<-0.15
  hq2=quantile(mdist,probs=h,na.rm=TRUE)
  }
else {
  if (is.numeric(h))    hq2<-h  
  else hq2=quantile(mdist,probs=as.numeric(h),na.rm=TRUE)
}
class(mdist) <- class(mdist2) <- c("matrix","fdist")
dep<-Ker.norm(mdist2/hq2)    ####
dep<-apply(dep,1,sum,na.rm=TRUE)                                    
if (scale)   {
	dep2=Ker.norm(mdist/hq2)
	dep2=apply(dep2,1,sum,na.rm=TRUE)
   mn<-min(dep2,na.rm=TRUE)
   mx<-max(dep2,na.rm=TRUE)
   dep=as.vector(dep/mx)   
}
if (nullans) {
        ans1<-rep(NA,len=m0)
        ans1[-nas] <-dep
        dep<-ans1      
        }
names(dep)<-nms		
k=which.max(dep)
med=fdataobj[k]
   nl=length(trim)
   lista=vector("list",nl)
   tr<-paste("mode.tr",round(trim*100,2),"\u0025",sep="")
   if (nl>1) names(lista)=paste0("tr",round(trim*100,2))
   mtrim=fdata(matrix(NA,ncol=length(tt),nrow=length(trim)),tt,rtt,names2)
   for (j in 1:nl){
   lista[[j]]=which(dep>=quantile(dep,probs=trim[j],na.rm=TRUE))
   	if (length(lista[[j]])==1) {
			mtrim$data[j,]<-fdataobj2[lista[[j]]]$data
	if (draw) {draw=FALSE;warning("Too few curves in mtrim. The plot is not drawn")}
	}
	else   mtrim$data[j,]=apply(fdataobj2$data[lista[[j]],,drop=FALSE],2,mean,na.rm=TRUE)
	}
   rownames(med$data)<-"mode.med"
   rownames(mtrim$data)<-tr
   out<-list("median"=med,"lmed"=k,"mtrim"=mtrim,"ltrim"=if (nl==1) unlist(lista) else lista,
"dep"=dep,"hq"=hq2,name="Mode") 
if (scale) out$dscale=mx
   out$trim <- trim
   out$name <- "mode"
   out$fdataobj <- fdataobj
   out$fdataori <- fdataori
  class(out) <- "depth"
  if (draw){
		plot.depth(out)
  }
  return(invisible(out))
}

#' @rdname depth.fdata
#' @export depth.RP
depth.RP<-function(fdataobj,fdataori=fdataobj,trim=0.25,nproj=50,proj="vexponential",
                   dfunc="TD1",par.dfunc=list(),scale=FALSE,draw=FALSE,...){
  if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
  if (!is.fdata(fdataori)) fdataori=fdata(fdataori)
  
  # nas<-is.na.fdata(fdataobj)
  # if (any(nas))  {
  #    fdataobj$data<-fdataobj$data[!nas,]
  #    cat("Warning: ",sum(nas)," curves with NA are not used in the calculations \n")
  #    }
  
  if (is.null(rownames(fdataobj$data)))  rownames(fdataobj$data)<-1:nrow(fdataobj$data)
  nms<-rownames(fdataobj$data)
  m0<-nrow(fdataobj)
  fdataobj2<-fdataobj
  fdataobj<-na.omit.fdata(fdataobj)
  fdataori<-na.omit.fdata(fdataori) 
  nas<-na.action(fdataobj)
  nullans<-!is.null(nas)
  #data<-fdataobj[["data"]]
  #data2<-fdataori[["data"]]
  n<-nrow(fdataobj)
  m<-ncol(fdataobj)
  m2<-ncol(fdataori)
  n2<-nrow(fdataori)
  if (is.null(n) && is.null(m))  stop("ERROR IN THE DATA DIMENSIONS")
  if (is.null(m) && is.null(m2)) stop("ERROR IN THE DATA DIMENSIONS")
  tt=fdataobj[["argvals"]]
  rtt<-fdataobj[["rangeval"]]
  names1<-names2<-fdataobj[["names"]]
  names1$main<-"depth.RP median"
  names2$main<-paste("RP trim",trim*100,"\u0025",sep="")
  
  #### new
  if (is.fdata(proj)) {
    nproj<-nrow(proj)
    if (fdataobj$argvals!=proj$argvals || ncol(fdataobj)!=ncol(proj)) stop("Error in proj dimension")
    z<-proj
  }
  else {
    z<-rproc2fdata(nproj,tt,sigma=proj,norm=TRUE,...)
  }
  ##
  dep=matrix(0.0,ncol=nproj,nrow=n)
  dep2=matrix(0.0,ncol=nproj,nrow=n2)
  if (dfunc %in% c("Liu1","TD1","FM1")) {
    Fn<-list()
    for (j in 1:nproj){
      valor=inprod.fdata(fdataobj,z[j])
      valor2=inprod.fdata(fdataori,z[j])	
      Fn[[j]]=ecdf(valor2)
      par.dfunc$x<-valor
      par.dfunc$Fn<-Fn[[j]]
      dep[,j]<-do.call(dfunc,par.dfunc)
      if (scale) {par.dfunc$x=valor2;dep2[,j]<-do.call(dfunc,par.dfunc)}
    }
    dep=apply(dep,1,mean)   #dep/nproj
  }
  else if (dfunc %in% c("MhD1","LD1")) {
    for (j in 1:nproj){
      valor=inprod.fdata(fdataobj,z[j])
      valor2=inprod.fdata(fdataori,z[j])
      par.dfunc$x<-drop(valor)
      par.dfunc$xx<-drop(valor2)
      dep[,j]<-do.call(dfunc,par.dfunc)
      if (scale) {par.dfunc$x=valor2;dep2[,j]<-do.call(dfunc,par.dfunc)}		 
      #         dep<-dep+dp
    }
    dep=apply(dep,1,mean)   #dep/nproj
  }
  if (scale){    dep2=apply(dep2,1,mean);dep<-dep/max(dep2)    }
  names(dep)<-nms
  k=which.max(dep)
  med=fdataobj[k]
  if (nullans) {
    ans1<-rep(NA,len=m0)
    ans1[-nas] <-dep
    dep<-ans1
  }
  nl=length(trim)
  lista=vector("list",nl)
  tr<-paste("RP.tr",round(trim*100,2),"\u0025",sep="")
  if (nl>1) names(lista)=paste0("tr",round(trim*100,2))
  mtrim=fdata(matrix(NA,ncol=length(tt),nrow=length(trim)),tt,rtt,names1)
  for (j in 1:nl){
    lista[[j]]=which(dep>=quantile(dep,probs=trim[j],na.rm=TRUE))
    if (length(lista[[j]])==1) {
      mtrim$data[j,]<-fdataobj2[lista[[j]]]$data
      if (draw) {draw=FALSE;warning("Too few curves in mtrim. The plot is not drawn")}
    }
    else   mtrim$data[j,]=apply(fdataobj2$data[lista[[j]],,drop=FALSE],2,mean,na.rm=TRUE)
  }
  rownames(med$data)<-"RP.med"
  rownames(mtrim$data)<-tr
  out=list("median"=med,"lmed"=k,"mtrim"=mtrim,"ltrim"=if (nl==1) unlist(lista) else lista, "dep"=dep,"proj"=z,dfunc=dfunc,par.dfunc=par.dfunc)
  out$trim <- trim
  out$name <- "RP"
  out$fdataobj <- fdataobj
  out$fdataori <- fdataori
  class(out) <- "depth"
  if (draw){
    plot.depth(out)
  }
  return(invisible(out))
}

#' @rdname depth.fdata
#' @export depth.RPD
depth.RPD<-function (fdataobj,fdataori=fdataobj, nproj = 20, proj=1,deriv = c(0, 1), trim = 0.25,
                     dfunc2 = mdepth.LD, method = "fmm", draw = FALSE, ...)
{
  if (!is.fdata(fdataobj))         fdataobj = fdata(fdataobj)
  if (!is.fdata(fdataori)) fdataobj=fdata(fdataori)    
  if (is.null(rownames(fdataobj$data)))  rownames(fdataobj$data)<-1:nrow(fdataobj$data)
  nms<-rownames(fdataobj$data)
  m0<-nrow(fdataobj)
  fdataobj2<- fdataobj
  fdataobj<-na.omit.fdata(fdataobj)
  fdataori<-na.omit.fdata(fdataori) 
  nas<-na.action(fdataobj)
  nullans<-!is.null(nas) 
  names1 <- names2 <- fdataobj[["names"]]
  names1$main <- "depth.RPD median"
  names2$main <- paste("RPD trim ", trim * 100, "%", sep = "")
  n <- nrow(fdataobj)
  m <- ncol(fdataobj)
  m2<-ncol(fdataori)
  n2<-nrow(fdataori)    
  if (is.null(m) && is.null(m2)) stop("ERROR IN THE DATA DIMENSIONS")
  if (is.null(n) || is.null(m))       stop("Input must be a matrix")
  tt = fdataobj[["argvals"]]
  rtt <- fdataobj[["rangeval"]]
  newfunc=vector("list",length(deriv))
  newfunc2=vector("list",length(deriv))
  for (ider in 1:length(deriv)) {
    if (deriv[ider] == 0) {   
      newfunc[[ider]]=fdataobj
      newfunc2[[ider]]=fdataori
    } 
    else {
      newfunc[[ider]] = fdata.deriv(fdataobj, nderiv = deriv[ider], method = method, ...)
      newfunc2[[ider]] = fdata.deriv(fdataori, nderiv = deriv[ider], method = method, ...) 
    }
  }
  dep = rep(0, n)
  dep2 = rep(0, n2)    
  vproject = matrix(0, nrow = n, ncol = length(deriv))
  vproject2 = matrix(0, nrow = n2, ncol = length(deriv))    
  if (is.fdata(proj)) {
    if (fdataobj$argvals!=proj$argvals || m!=ncol(proj)) stop("Error in proj dimension")
    z<-proj
    nproj<-nrow(z)
  }
  else {	 z<-rproc2fdata(nproj,tt,sigma=proj,norm=TRUE,...)	}
  pb = txtProgressBar(min = 0, max = nproj, style = 3)
  for (j in 1:nproj) {
    setTxtProgressBar(pb, j - 0.5)
    for (ider in 1:length(deriv)) {
      vproject[, ider] = inprod.fdata(newfunc[[ider]],z[j])
      vproject2[, ider] = inprod.fdata(newfunc2[[ider]],z[j]) 
    }
    par.dfunc = list()
    par.dfunc$x <- vproject
    par.dfunc$xx <- vproject2
    #        par.dfunc$trim <- trim              
    par.dfunc$scale<-TRUE
    resul = do.call(dfunc2, par.dfunc)
    dep = dep + resul$dep
    setTxtProgressBar(pb, j)
  }                                                               
  close(pb)
  names(dep)<-nms       
  dep = dep/nproj
  if (nullans) {
    ans1<-rep(NA,len=m0)
    ans1[-nas] <-dep
    dep<-ans1      
  }
  k = which.max(dep)    
  med = fdataobj[k]
  nl=length(trim)
  lista=vector("list",nl)
  tr<-paste("RPD.tr",round(trim*100,2),"\u0025",sep="")
  if (nl>1) names(lista)=paste0("tr",round(trim*100,2))
  mtrim=fdata(matrix(NA,ncol=length(tt),nrow=length(trim)),tt,rtt,names2)
  for (j in 1:nl){
    lista[[j]]=which(dep>=quantile(dep,probs=trim[j],na.rm=TRUE))
    if (length(lista[[j]])==1) {
      mtrim$data[j,]<-fdataobj2[lista[[j]]]$data
      if (draw) {draw=FALSE;warning("Too few curves in mtrim. The plot is not drawn")}
    }
    else   mtrim$data[j,]=apply(fdataobj2$data[lista[[j]],,drop=FALSE],2,mean,na.rm=TRUE)
  }
  #else mtrim=data[lista,,drop=FALSE]
  #   med<-fdata(med,tt,rtt,names1)
  #   mtrim<-fdata(mtrim,tt,rtt,names2)
  rownames(med$data)<-"RPD.med"
  rownames(mtrim$data)<-tr
  out=list(median = med, lmed = k, mtrim = mtrim,
           ltrim = if (nl==1) unlist(lista) else lista, dep = dep,deriv=deriv,proj = z,name="RPD")
  out$trim <- trim
  out$name <- "RPD"
  out$fdataobj <- fdataobj
  out$fdataori <- fdataori
  class(out)="depth"
  if (draw) {
    plot.depth(out)
  }
  return(invisible(out))
}  

#' @rdname depth.fdata
#' @export 
depth.RT  <-function (fdataobj,fdataori=fdataobj, trim = 0.25, nproj = 10, proj = 1, xeps = 1e-07, 
                      draw = FALSE, ...) 
{
  if (!is.fdata(fdataobj)) fdataobj = fdata(fdataobj)
  if (!is.fdata(fdataori)) fdataobj=fdata(fdataori)    
  if (is.null(rownames(fdataobj$data)))  rownames(fdataobj$data)<-1:nrow(fdataobj$data)
  nms<-rownames(fdataobj$data)
  m0<-nrow(fdataobj)
  fdataobj2<-fdataobj
  fdataobj<-na.omit.fdata(fdataobj)
  nas<-na.action(fdataobj)
  nullans<-!is.null(nas) 
  
  fdataori=na.omit.fdata(fdataori)
  data  <- fdataobj[["data"]]
  data2 <- fdataori[["data"]]    
  n <- nrow(fdataobj)
  m <- ncol(fdataobj)
  m2<-ncol(fdataori)
  n2<-nrow(fdataori)
  if (is.null(n) && is.null(m)) stop("ERROR IN THE DATA DIMENSIONS OF fdataobj")
  if (is.null(n2) && is.null(m2)) stop("ERROR IN THE DATA DIMENSIONS fdataori")
  tt = fdataobj[["argvals"]]
  rtt <- fdataobj[["rangeval"]]
  names1 <- names2 <- names <- fdataobj[["names"]]
  names1$main <- "depth.RT median"
  names2$main <- paste("RT trim", trim * 100, "%", sep = "")
  
  if (is.fdata(proj)) {
    nproj <- nrow(proj)
    if (fdataobj$argvals != proj$argvals || ncol(fdataobj) != 
        ncol(proj)) 
      stop("Error en proj dimension")
    z <- proj
  }
  else {
    z <- rproc2fdata(nproj, tt, sigma = proj, norm = TRUE,...)
  }
  Fn <- list()
  Prod=inprod.fdata(fdataobj,z)
  Prod2=inprod.fdata(fdataori,z)
  dep = array(NA, dim = c(nproj, n))
  
  for (j in 1:nproj) {
    Fn[[j]] = ecdf(Prod2[,j])  
    dep[j, ] = pmin(Fn[[j]](Prod[,j]-xeps) ,(1 - Fn[[j]](Prod[,j]-xeps)))       
  }
  #    print(dep)
  dep = apply(dep, 2, min)
  if (nullans) {
    ans1<-rep(NA,len=m0)
    ans1[-nas] <-dep
    dep<-ans1      
  }
  names(dep)<-nms      
  k = which.max(dep)
  med = fdataobj2[k]
  nl = length(trim)
  lista=vector("list",nl)
  tr <- paste("RT.tr", round(trim * 100,2), "%", sep = "")	
  if (nl>1) names(lista)=paste0("tr",round(trim*100,2))
  mtrim = fdata(matrix(NA, nrow = nl, ncol = m),tt,rtt,names)
  for (j in 1:length(trim)) {
    lista[[j]] = which(dep >= quantile(dep, probs = trim[j], na.rm = TRUE))
    if (length(lista[[j]])==1) {
      mtrim$data[j,]<-fdataobj2[lista[[j]]]$data
      if (draw) {draw=FALSE;warning("Too few curves in mtrim. The plot is not drawn")}
    }
    else mtrim$data[j,]=func.mean(fdataobj2[lista[[j]]])$data       
  }
  rownames(med$data) <- "RT.med"
  rownames(mtrim$data) <- tr
  out=list(median = med, lmed = k, mtrim = mtrim, 
           ltrim = if (nl==1) unlist(lista) else lista, dep = dep, proj = z, Fn = Fn,name="RT")
  out$trim <- trim
  out$name <- "RT"
  out$fdataobj <- fdataobj
  out$fdataori <- fdataori
  class(out)="depth"
  if (draw) {
    plot.depth(out)
  }
  return(invisible(out))
}

#' @rdname depth.fdata
#' @export 
depth.KFSD=function (fdataobj, fdataori = fdataobj, trim = 0.25,
                     h=NULL, scale = FALSE, draw = FALSE){
  if (is.fdata(fdataobj) & is.fdata(fdataori)) {
    if (is.null(rownames(fdataobj$data))) 
      rownames(fdataobj$data) <- 1:nrow(fdataobj$data)
    nms <- rownames(fdataobj$data)
    m0 <- nrow(fdataobj)
    fdataobj <- na.omit.fdata(fdataobj)
    fdataori <- na.omit.fdata(fdataori)
    nas <- na.action(fdataobj)
    nullans <- !is.null(nas)
    names1 <- names2 <- names <- fdataobj[["names"]]
    names1$main <- "depth.KFSD median"
    names2$main <- paste("depth.KFSD trim ", trim * 100, 
                         "%", sep = "")
    tt=fdataobj$argvals
    rtt=fdataobj$rangval
    n <- nrow(fdataobj)
    m <- ncol(fdataobj)
    m2 <- ncol(fdataori)
    n2 <- nrow(fdataori)
  }
  else {
    stop("no fdata class object on input")
  }
  if (is.null(n) && is.null(m)) 
    stop("ERROR IN THE DATA DIMENSIONS")
  if (is.null(m) && is.null(m2)) 
    stop("ERROR IN THE DATA DIMENSIONS")
  mdist=matrix(NA,ncol=n2,nrow=n2)
  
  for (i in 1:(n2-1)){for (j in (i+1):n2){
    mdist[i,j]<-mdist[j,i]<-norm.fdata(fdataori[i]-fdataori[j])
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
  kern=function(x,y,h=hq2){exp(-norm.fdata(x-y)^2/h^2)}
  K02=rep(1,nrow(fdataobj))
  K01=rep(1,nrow(fdataori))
  M1=matrix(NA,nrow=n2,ncol=n2)
  M2=matrix(NA,nrow=n,ncol=n2)
  M=array(NA,dim=c(n,n2,n2))
  for (i in 1:n2){for (j in i:n2){
    if (i==j) M1[i,i]=K01[i] else M1[i,j]<-M1[j,i]<-kern(fdataori[i],fdataori[j],h=hq2)
  }}
  same.dim <- FALSE
  if (n==n2 & m==m2){ same.dim <- TRUE}
  if (same.dim)
    if (all(fdataobj$data==fdataori$data)) {
      M2=M1
    } else {same.dim <- FALSE}
  if (!same.dim){	
    for (i in 1:n){for (j in 1:n2){
      if (all(fdataobj[i]$data == fdataori[j]$data)) M2[i,j]=K02[i] else M2[i,j]=kern(fdataobj[i],fdataori[j],h=hq2)
    }}
  }  
  for (i in 1:n){for (j in 1:n2){for (k in 1:n2){
    M[i,j,k]<-(K02[i]+M1[j,k]-M2[i,j]-M2[i,k])/(sqrt(K02[i]+K01[j]-2*M2[i,j])*sqrt(K02[i]+K01[k]-2*M2[i,k]))
  }}}
  l=which(!is.finite(M))
  M[l]=NA
  # Bug 2019/08/30  (change n by n2)
  # dep=1-sqrt(apply(M,1,sum,na.rm=TRUE))/n
  dep=1-sqrt(apply(M,1,sum,na.rm=TRUE))/n2
  if (scale) {
    MO=array(NA,dim=c(n2,n2,n2))
    for (i in 1:n2){for (j in 1:n2){for (k in 1:n2){
      MO[i,j,k]<-(K01[i]+M1[j,k]-M1[i,j]-M1[i,k])/(sqrt(K01[i]+K01[j]-2*M1[i,j])*sqrt(K01[i]+K01[k]-2*M1[i,k]))
    }}}
    l=which(!is.finite(MO))
    MO[l]=NA
    dep2=1-sqrt(apply(MO,1,sum,na.rm=TRUE))/n2
    #mn <- min(dep2, na.rm = TRUE)
    #
    #print(range(dep));   print(range(dep2))
    mx <- max(dep, na.rm = TRUE)
    dep = as.vector(dep/mx)
  }
  if (nullans) {
    ans1 <- rep(NA, len = m0)
    ans1[-nas] <- dep
    dep <- ans1
  }
  names(dep)=nms
  k = which.max(dep)
  med = fdataobj[k]
  nl = length(trim)
  lista = vector("list",nl)
  tr <- paste("KFSD.tr", round(trim * 100,2), "%", sep = "")
  if (nl>1) names(lista)=paste0("tr",round(trim*100,2))
  mtrim = fdata(matrix(NA, nrow = nl, ncol = m),tt,rtt,names)
  for (j in 1:length(trim)) {
    lista[[j]] = which(dep >= quantile(dep, probs = trim[j], na.rm = TRUE))
    if (length(lista[[j]])==1) {
      mtrim$data[j,]<-fdataobj[lista[[j]]]$data
      if (draw) {draw=FALSE;warning("Too few curves in mtrim. The plot is not drawn")}
    }
    else mtrim$data[j,]=func.mean(fdataobj[lista[[j]]])$data       
  }
  rownames(med$data) <- "KFSD.med"
  rownames(mtrim$data) <- tr
  out <- list(median = med, lmed = k, mtrim = mtrim, ltrim = if (nl==1) unlist(lista) else lista, 
              dep = dep, h = h, hq=hq2,name="KFSD")
  if (scale) out$dscale = mx
  class(out)="depth"
  out$trim <- trim
  out$name <- "KFSD"
  out$fdataobj=fdataobj
  out$fdataori=fdataori
  if (draw) {
    plot.depth(out)
  }
  return(invisible(out))	
}

#' @rdname depth.fdata
#' @export 
depth.FSD=function (fdataobj, fdataori = fdataobj, 
                    trim = 0.25, scale = FALSE, draw = FALSE){
  if (is.fdata(fdataobj)) {
    if (is.null(rownames(fdataobj$data))) 
      rownames(fdataobj$data) <- 1:nrow(fdataobj$data)
    m0=nrow(fdataobj)
    nms <- rownames(fdataobj$data)
    fdataobj <- na.omit.fdata(fdataobj)
    fdataori <- na.omit.fdata(fdataori)
    nas <- na.action(fdataobj)
    n <- nrow(fdataobj)
    m <- ncol(fdataobj)
    m2 <- ncol(fdataori)
    n2 <- nrow(fdataori)
    nullans <- !is.null(nas)
    names1 <- names2 <- names <- fdataobj[["names"]]
    names1$main <- "depth.FSD median"
    names2$main <- paste("depth.FSD trim ", trim * 100, 
                         "%", sep = "")
    tt = fdataobj[["argvals"]]
    rtt <- fdataobj[["rangeval"]]
  } else {
    stop("no fdata class object")
  }
  if (is.null(n) && is.null(m)) 
    stop("ERROR IN THE DATA DIMENSIONS")
  if (is.null(m) && is.null(m2)) 
    stop("ERROR IN THE DATA DIMENSIONS")
  kern=function(x,y){inprod.fdata(x,y)}
  K02=norm.fdata(fdataobj)^2
  K01=norm.fdata(fdataori)^2
  M1=matrix(NA,nrow=n2,ncol=n2)
  M2=matrix(NA,nrow=n,ncol=n2)
  M=array(NA,dim=c(n,n2,n2))
#  if (scale) MO=array(NA,dim=c(n2,n2,n2))
  for (i in 1:n2){
    for (j in i:n2){
      if (i==j) 
        M1[i,i]=K01[i] else M1[i,j]<-M1[j,i]<-kern(fdataori[i],fdataori[j])
    }
  }
  same.dim <- FALSE
  if (n==n2 & m==m2){ same.dim <- TRUE}
  if (same.dim)
    if (all(fdataobj$data==fdataori$data)) {
      M2=M1
    } else {same.dim <- FALSE}
  if (!same.dim){
  #  print("entra")
    for (i in 1:n){
      for (j in 1:n2){
      if (all(fdataobj[i]$data == fdataori[j]$data)) 
        M2[i,j]=K02[i] else {
          M2[i,j]=kern(fdataobj[i],fdataori[j])
        }
    }
    }
  }
  #print(M2)
  for (i in 1:n){
    for (j in 1:n2){
      for (k in 1:n2){
       if (all(fdataobj[i]$data == fdataori[j]$data) | 
              all(fdataobj[i]$data == fdataori[k]$data)) 
               M[i,j,k]=NA else M[i,j,k]<-(K02[i]+M1[j,k]-M2[i,j]-M2[i,k])/(sqrt(K02[i]+K01[j]-2*M2[i,j])*sqrt(K02[i]+K01[k]-2*M2[i,k]))
      }}}
  l=which(!is.finite(M))
  M[l]=NA
  # Bug 2019/08/30  (change n by n2)
  # dep=1-sqrt(apply(M,1,sum,na.rm=TRUE))/n
  dep=1-sqrt(apply(M,1,sum,na.rm=TRUE))/n2
  if (scale) {
    MO=array(NA,dim=c(n2,n2,n2))
    for (i in 1:n2){for (j in 1:n2){for (k in 1:n2){
      if (all(fdataori[i]$data == fdataori[j]$data) | all(fdataori[i]$data == fdataori[k]$data)) MO[i,j,k]=NA else MO[i,j,k]<-(K01[i]+M1[j,k]-M1[i,j]-M1[i,k])/(sqrt(K01[i]+K01[j]-2*M1[i,j])*sqrt(K01[i]+K01[k]-2*M1[i,k]))
    }}}
    l=which(!is.finite(MO))
    MO[l]=NA
    dep2=1-sqrt(apply(MO,1,sum,na.rm=TRUE))/n2
 #   mn <- min(dep2, na.rm = TRUE) 
    mx <- max(dep2, na.rm = TRUE)
    dep = as.vector(dep/mx)
  }
 # print("M M2 M")
#  print(range(M,na.rm=T))
 #       print(range(M2,na.rm=T))
  #      print(range(M0,na.rm=T)
  if (nullans) {
    ans1 <- rep(NA, len = m0)
    ans1[-nas] <- dep
    dep <- ans1
  }
  names(dep)=nms
  k = which.max(dep)
  med = fdataobj[k]
  nl = length(trim)
  lista=vector("list",nl)
  tr <- paste("FSD.tr", round(trim * 100,2), "%", sep = "")
  if (nl>1) names(lista)=paste0("tr",round(trim*100,2))	
  mtrim = fdata(matrix(NA, nrow = nl, ncol = m),tt,rtt,names2)
  for (j in 1:length(trim)) {
    lista[[j]] = which(dep >= quantile(dep, probs = trim[j], na.rm = TRUE))
    if (length(lista[[j]])==1) {
      mtrim$data[j,]<-fdataobj[lista[[j]]]$data
      if (draw) {draw=FALSE;warning("Too few curves in mtrim. The plot is not drawn")}
    }
    else mtrim$data[j,]=func.mean(fdataobj[lista[[j]]])$data       
  }
  rownames(med$data) <- "FSD.med"
  rownames(mtrim$data) <- tr
  out <- list(median = med, lmed = k, mtrim = mtrim, ltrim = if (nl==1) unlist(lista) else lista, 
              dep = dep,name="FSD")
  if (scale) 
    out$dscale = mx
  class(out)="depth"
  out$trim <- trim
  out$name <- "FSD"
  out$fdataobj=fdataobj
  out$fdataori=fdataori
  if (draw) {
    plot.depth(out)
  }
  return(invisible(out))
}

#' @rdname depth.fdata
#' @export
depth.FM=function(fdataobj,fdataori=fdataobj,trim=0.25,scale=FALSE,dfunc="FM1",par.dfunc=list(scale=TRUE),draw=FALSE){
  if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
  if (!is.fdata(fdataori)) fdataori=fdata(fdataori)
  if (is.null(fdataobj))  rownames(fdataobj$data)<-1:nrow(fdataobj$data)
  nms<-rownames(fdataobj$data)
  
  #nas<-apply(fd-ataobj$data,1,count.na)
  #if (any(nas))  {
  #   fdataobj$data<-fdataobj$data[!nas,]
  #   cat("Warning: ",sum(nas)," curves with NA are not used in the calculations \n")
  #   }             
  data<-fdataobj[["data"]]
  data2<-fdataori[["data"]]
  n<-nrow(data)
  m<-ncol(data)
  m2<-ncol(data2)
  n2<-nrow(data2)
  d<-matrix(NA,nrow=n,ncol=m)
  if (is.null(n) && is.null(m)) stop("ERROR IN THE DATA DIMENSIONS")
  if (is.null(m) && is.null(m2)) stop("ERROR IN THE DATA DIMENSIONS")
  if (is.list(dfunc)) dfunc<-dfunc$dfunc
  t=fdataobj[["argvals"]]
  rtt<-fdataobj[["rangeval"]]
  names1<-names2<-names<-fdataobj[["names"]]
  names1$main<-"depth.FM median"
  Fn<-list()
  tr<-paste("FM.tr",trim*100,"\u0025",sep="")
  names2$main<- tr
  dtt <- diff(t)
  eps <- as.double(.Machine[[1]] * 100)
  inf <- dtt - eps
  sup <- dtt + eps
  if (all(dtt > inf) & all(dtt < sup)) {  equi = TRUE  }
  else equi = FALSE
  for (i in 1:m)   {
    if (dfunc %in% c("TD1","Liu1","FM1")){
      Fn[[i]]=ecdf(data2[,i])
      par.dfunc$x=data[,i]
      par.dfunc$Fn=Fn[[i]]
      d[,i]=do.call(dfunc,par.dfunc)
    }
    else     {
      par.dfunc$x=data[,i]
      par.dfunc$xx=data2[,i]
      d[,i]=do.call(dfunc,par.dfunc)
    }
  }
  #d<-int.simpson(fdata(d,t,rtt),equi=equi)
  d<-apply(d,1,mean)[1:n]
  ans<-d
  if (scale) { ans=d/max(d)}
  names(ans)<-nms
  k=which.max(ans)    
  med=data[k,]                                                      
  lista=which(ans>=quantile(ans,probs=trim,na.rm=TRUE))
  if (n>1)  mtrim=apply(data[lista,,drop=FALSE],2,mean,na.rm=TRUE)
  else mtrim<-data                                     
  med<-fdata(med,t,rtt,names1)
  mtrim<-fdata(mtrim,t,rtt,names2)
  rownames(med$data)<-"FM.med"
  rownames(mtrim$data)<-tr
  out=list("median"=med,"lmed"=k,"mtrim"=mtrim,"ltrim"=lista,
           "dep"=ans,"Fn"=Fn)
  class(out)="depth"
  out$trim= out$trim
  out$fdataobj=fdataobj
  out$fdataori=fdataori
  if (draw){
    plot.depth(out)
  }
  return(invisible(out))
} 



depth.MB<-function (fdataobj, fdataori = NULL,trim=0.25, scale=FALSE,
                    draw =FALSE, grayscale = FALSE, band = FALSE, band.limits = NULL, lty = 1,
                    lwd = 2, col = NULL, cold = NULL, colRef = NULL, ylim = NULL, cex = 1, ...)
{
  ################################################################################
  # Wrapper version of modified band depth (original code from depthTools::MBD) 
  # adapted from mulitvariate to functioinal case by Manuel Oviedo de la Fuente 
  ################################################################################
  if (!is.fdata(fdataobj)) {fdataobj=fdata(fdataobj)}
  x<-fdataobj$data
  n <- nrow(x)
  d <- ncol(x)
  x <- as.matrix(x)
  tt<-fdataobj$argvals
  rtt<-diff(fdataobj$rangeval)
  dtt<-diff(tt/rtt)
  depth.ori<-NULL
  #    dtt<-c(0,dtt,0)/sum(dtt)*d
  if (length(fdataori) == 0) {
    if (ncol(x) == 1) {
      x <- t(x)
    }
    depth <- matrix(0, n,d)
    ordered.matrix <- x
    if (n > 1) {
      for (columns in 1:d) {
        if (columns>1 & columns<d)
          wei <- .5*dtt[columns]+.5*dtt[columns-1]
        else {
          if (columns==1)     {    wei <-.5*dtt[columns]}
          else {if (columns==d)    wei <- .5*dtt[columns-1]}
        }
        ordered.matrix[, columns] <- sort(x[, columns])
        for (element in 1:n) {
          index.1 <- length(which(ordered.matrix[, columns] <
                                    (x[element, columns])))
          index.2 <- length(which(ordered.matrix[, columns] <=
                                    (x[element, columns])))
          multiplicity <- index.2 - index.1
          depth[element,columns] <- index.1 *
            (n - (index.2)) + multiplicity * (n - index.2 +
                                                index.1) + choose(multiplicity, 2)
        }
        depth[, columns] <- depth[, columns]* wei *d
      }
      depth <- rowSums(depth/(d * choose(n, 2)))
    }
    if (n == 1) {
      deepest <- x
      depth <- 0
    }
    ordering <- order(depth, decreasing = TRUE)
    #########################################
    if (draw) {
      par(mar = c(4, 5, 3, 3), xpd = FALSE)
      fobj <- fdataobj[ordering[n:2], ]
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
        plot(fobj, ylim = ylim, type = "l",lty = 0, ...)
        for (poly in no.poly:1) {
          limit <- band.limits[poly]
          no.points <- floor(limit * n)
          xx <- t(x[ordering[no.points:1],])
          upper <- apply(xx, 1, max)
          lower <- apply(xx, 1, min)
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
        plot(fobj, type = "l", ylim = ylim,lty = lty, lwd = lwd, col = color[n:2], ...)
      }
      lines(fdataobj[ordering[1]], lty = lty, lwd = lwdd, col = cold)
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
    if (!is.fdata(fdataori)) fdataori=fdata(fdataori)
    xRef<-fdataori$data
    xRef <- as.matrix(xRef)
    if (ncol(xRef) != d) {
      stop("Dimensions of x and xRef do not match")
    }
    n0 <- nrow(xRef)
    if (ncol(x) == 1) {
      x <- t(x)
    }
    depth <- matrix(0, n,d)
    depth.ori <- matrix(0, n0,d)
    ordered.matrix <- xRef
    if (n0 > 1) {
      for (columns in 1:d) {
        if (columns>1 & columns<d)
          wei <- .5*dtt[columns]+.5*dtt[columns-1]
        else {
          if (columns==1)     {    wei <-.5*dtt[columns]}
          else {if (columns==d)    wei <- .5*dtt[columns-1]}
        }
        ordered.matrix[, columns] <- sort(xRef[, columns])
        for (element in 1:n) {
          index.1 <- length(which(ordered.matrix[, columns] <
                                    x[element, columns]))
          index.2 <- length(which(ordered.matrix[, columns] <=
                                    x[element, columns]))
          multiplicity <- index.2 - index.1
          depth[element,columns] <- (index.1 +
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
          depth.ori[element,columns] <- (index.1 +
                                           multiplicity) * (n0 - index.1 - multiplicity) +
            multiplicity * (index.1 + (multiplicity -
                                         1)/2)
        }
        depth[, columns] <- depth[, columns]* wei *d
        depth.ori[, columns] <- depth.ori[, columns]* wei *d
      }
      depth <- rowSums(depth/(d * choose(n0, 2)))
      depth.ori <- rowSums(depth.ori/(d * choose(n0, 2)))
      
      #            depth <- depth/(d * choose(n0, 2))
    }
    if (n == 1) {
      deepest <- x
      depth <- 0
    }
    ordering <- order(depth, decreasing = TRUE)
    if (draw) {
      par(mar = c(4, 5, 3, 3), xpd = FALSE)
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
        plot(fdataori, type = "l", ylim = ylim,
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
          xx<- t(x[ordering[no.points:1],])
          upper <- apply(xx, 1, max)
          lower <- apply(xx, 1, min)
          polygon(c(tt, rev(tt)), c(upper, rev(lower)),
                  col = color[poly])
        }
        lines(fdataobj[ordering[1], ], lty = lty, lwd = lwdd,
              col = cold)
      }
      else {
        plot(fdataori, type = "l", ylim = ylim,
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
        lines(fdataobj[ordering[n:2], ], lwd = lwd, col = color[n:2],
              lty = lty)
        lines(fdataobj[ordering [1], ], lwd = lwdd, lty = lty[1],
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
    if (scale)     depth.ori<-depth.ori/max(depth.ori)           
  }
  med<-fdataobj[ordering[1]]
  lista = which(depth >= quantile(depth, probs = trim, na.rm = TRUE))
  mtrim= apply(x[lista, ], 2, mean)
  tr<-paste("band.tr",trim*100,"\u0025",sep="")
  names2<-fdataobj[["names"]]
  med$names$main<-"depth.mode median"
  names2$main<-paste("depth.mbandtrim ",trim*100,"\u0025",sep="")
  mtrim<-fdata(mtrim,tt,rtt,names2)
  rownames(med$data)<-"mode.med"
  rownames(mtrim$data)<-tr
  if (scale)    depth<-depth/max(depth)
  return(invisible(list("median"=med,"lmed"=ordering[1],ordering=ordering,
                        "mtrim"=mtrim,"ltrim"=lista,"dep"=depth,"dep.ori"=depth.ori)))
}

