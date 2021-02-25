#' @title Bootstrap samples of a functional statistic
#' 
#' @description provides bootstrap samples for functional data.
#' 
#' @param fdataobj \code{fdata} class object.
#' @param statistic Sample statistic. It must be a function that returns an
#' object of class \code{fdata}. By default, it uses sample mean
#' \code{\link{func.mean}}.  See \code{\link{Descriptive}} for other
#' statistics.
#' @param alpha Significance value.
#' @param nb Number of bootstrap resamples.
#' @param smo The smoothing parameter for the bootstrap samples as a proportion
#' of the sample variance matrix.
#' @param draw If \code{TRUE}, plot the bootstrap samples and the statistic.
#' @param draw.control List that it specifies the \code{col}, \code{lty} and
#' \code{lwd} for objects: \code{fdataobj}, \code{statistic}, \code{IN} and
#' \code{OUT}.
#' @param \dots Further arguments passed to or from other methods.
#' 
#' @details The \code{fdata.bootstrap()} computes a confidence ball using bootstrap in
#' the following way: 
#' \itemize{ 
#' \item Let \eqn{X_1(t),\ldots,X_n(t)}{X_1(t),...,X_n(t)} the original data and
#' \eqn{T=T(X_1(t),\ldots,X_n(t))}{T=T(X_1(t),...,X_n(t))} the sample' statistic.
#' 
#' \item Calculate the \code{nb} bootstrap resamples
#' \eqn{\left\{X_{1}^{*}{(t)},\cdots,X_n^*(t)\right\}}{(X*_1(t),...,X*_n(t))},
#' using the following scheme \eqn{X_i^*(t)=X_i(t)+Z(t)}{X*_i(t)=X_i(t)+Z(t)} 
#' where \eqn{Z(t)}{Z(t)} is normally distributed with mean 0 and covariance matrix
#' \eqn{\gamma\Sigma_x}{\gamma\Sigma_x}, where \eqn{\Sigma_x}{\Sigma_x} is the
#' covariance matrix of' \eqn{\left\{X_1(t),\ldots,X_n(t)\right\}}{(X*_1(t),...,X*_n(t))} 
#' and \eqn{\gamma}{\gamma} is the smoothing parameter.  
#' 
#' \item Let \eqn{T^{*j}=T(X^{*j}_1(t),...,X^{*j}_n(t))}{T^{*j}=T(X^{*j}_1(t),...,X^{*j}_n(t))}
#' the estimate using the \eqn{j} resample.  
#' 
#' \item Compute \eqn{d(T,T^{*j})}, \eqn{j=1,\ldots,nb}. Define the bootstrap 
#' confidence ball of level \eqn{1-\alpha}{1-\alpha} as \eqn{CB(\alpha)=X\in E}{CB(\alpha)=X \in E} 
#' such that \eqn{d(T,X)\leq d_{\alpha}}{d(T,X)<= d\alpha} being
#' \eqn{d_{\alpha}}{d\alpha} the quantile \eqn{(1-\alpha)}{(1-\alpha)} of the
#' distances between the bootstrap resamples and the sample estimate. 
#' }
#' 
#' The \code{fdata.bootstrap} function allows us to define a statistic
#' calculated on the \code{nb} resamples, control the degree of smoothing by
#' \code{smo} argument and represent the confidence ball with level
#' \eqn{1-\alpha}{1-\alpha} as those resamples that fulfill the condition of
#' belonging to \eqn{CB(\alpha)}{CB(\alpha)}.  The \code{statistic} used by
#' default is the mean (\code{\link{func.mean}}) but also other depth-based
#' functions can be used (see \code{help(Descriptive)}).
#' 
#' @return 
#' \itemize{
#' \item{statistic}{ \code{fdata} class object with the statistic
#' estimate from \code{nb} bootstrap samples.} 
#' \item{dband}{ Bootstrap estimate of \code{(1-alpha)\%} distance.} 
#' \item{rep.dist}{ Distance from every replicate.} 
#' \item{resamples}{ \code{fdata} class object with the bootstrap resamples.} 
#' \item{fdataobj}{ \code{fdata} class object.}
#' }
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' @seealso See Also as \code{\link{Descriptive}}
#' 
#' @references Cuevas A., Febrero-Bande, M. and Fraiman, R. (2007).
#' \emph{Robust estimation and classification for functional data via
#' projection-based depth notions.} Computational Statistics 22, 3: 481-496.
#' 
#' Cuevas A., Febrero-Bande, M., Fraiman R. 2006.  \emph{On the use of
#' bootstrap for estimating functions with functional data.} Computational
#' Statistics and Data Analysis 51: 1063-1074.
#' 
#' Febrero-Bande, M., Oviedo de la Fuente, M. (2012).  \emph{Statistical
#' Computing in Functional Data Analysis: The R Package fda.usc.} Journal of
#' Statistical Software, 51(4), 1-28. \url{http://www.jstatsoft.org/v51/i04/}
#' @keywords datagen
#' @examples
#' \dontrun{
#' data(tecator)
#' absorp<-tecator$absorp.fdata
#' # Time consuming
#' #Bootstrap for Trimmed Mean with depth mode
#' out.boot=fdata.bootstrap(absorp,statistic=func.trim.FM,nb=200,draw=TRUE)
#' names(out.boot)
#' #Bootstrap for Median with with depth mode
#' control=list("col"=c("grey","blue","cyan"),"lty"=c(2,1,1),"lwd"=c(1,3,1))
#' out.boot=fdata.bootstrap(absorp,statistic=func.med.mode,
#' draw=TRUE,draw.control=control)
#' }
#' @aliases fdata.bootstrap fdata.bootstrap2
#' 
#' @export 
fdata.bootstrap <- function(fdataobj, statistic  =func.mean
                            ,alpha = 0.05, nb= 200, smo = 0
                            ,draw = FALSE,draw.control = NULL,...){
  if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
  data<-fdataobj[["data"]]
  if (smo > 0) varX <- var(data) * smo
  estmues<-statistic(fdataobj,...)
  nr<-nrow(fdataobj)
  nc<-ncol(fdataobj)
  tt =fdataobj[["argvals"]]
  rtt=fdataobj[["rangeval"]]
  names=fdataobj[["names"]]
  distboot<-matrix(NA,nrow=nb)
  #estboot<-matrix(NA,nrow=nb,ncol=nc)
  #pb=txtProgressBar(min=0,max=nb,style=3)
  #for (i in 1:nb){
 
  par.fda.usc <- eval(parse(text="fda.usc:::par.fda.usc"), envir=.GlobalEnv)
  ncores <- par.fda.usc$ncores
  #group.bootstrap<-fda.usc:::group.bootstrap
  #ncores <- ops.fda.usc()$ncores
  inumgr <- icount(nb)

  estboot <-foreach(i = inumgr, .packages='fda.usc', .inorder = FALSE, .combine="rbind") %dopar% {
    # set.seed(i)
    bmuestra <- fdataobj[sample(1:nr,size=nr,replace=TRUE),]
    if ( smo > 0) {
      bmuestra[["data"]] <- bmuestra[["data"]] + mvrnorm(n=nr, rep(0,nc), varX )
    }
    stat <- statistic(bmuestra,...)
    #estboot[i,] <- 
    stat[["data"]]
  }
  # pb=txtProgressBar(min=0,max=nb,style=3)
  # for (i in 1:nb){
  #   setTxtProgressBar(pb,i-0.5)
  #   bmuestra<-fdataobj[sample(1:nr,size=nr,replace=TRUE),]
  #   if (smo>0) {
  #     bmuestra[["data"]]<-bmuestra[["data"]]+mvrnorm(n=nr,rep(0,nc),var(data)*smo)
  #   }
  #   stat<-statistic(bmuestra,...)
  #   estboot[i,]<-stat[["data"]]
  #   setTxtProgressBar(pb,i)
  # }
  # close(pb)
  
  center<-estmues
  #for (i in 1:nb){  aux<-fdata(estboot[i,],tt,rtt)
  #  distboot[i]<-metric.lp(center,aux,...)  }
  #print(dim(estboot))  print(class(estboot))
  
  distboot <- metric.lp(fdata(estboot,tt,rtt),center,...)
  dist <- max(distboot[rank(distboot)<=floor((1-alpha)*nb)])
  resample <- fdata(estboot,tt,rtt,names)
  if (draw){
    if (is.null(draw.control)) draw.control=list("col"=c("grey","blue","pink"),"lty"=c(2,1,1),"lwd"=c(1,2,1))
    if (is.null(draw.control$lwd)) draw.control$lwd=c(1,2,1)
    if (is.null(draw.control$lty)) draw.control$lty=c(2,1,1)
    if (is.null(draw.control$col)) draw.control$col=c("grey","blue","pink")
    plot(fdataobj,lwd=draw.control$lwd[1],lty=draw.control$lty[1],col=draw.control$col[1])
    # for(i in 1:nb){
    #if (distboot[i]<=dist) lines(tt,estboot[i,],lwd=draw.control$lwd[3],lty=draw.control$lty[3],col=draw.control$col[3])
    #else lines(tt,estboot[i,],lwd=draw.control$lwd[4],lty=draw.control$lty[4],
    #col=draw.control$col[4])  }
    lines(resample[distboot<=dist],lwd=draw.control$lwd[3],lty=draw.control$lty[3],col=draw.control$col[3])
    lines(estmues,lwd=draw.control$lwd[2],lty=draw.control$lty[2],col=draw.control$col[2])
    legend(x=min(tt),y=0.99*max(data),legend=c("original curves", estmues$names$main,"bootstrap curves IN"),
           lty=draw.control$lty,lwd=draw.control$lwd,col=draw.control$col,cex=0.9,box.col=0)
  }
  return(list("statistic" = estmues, "dband" = dist, "rep.dist"=distboot,
              "resample"=resample,fdataobj=fdataobj))
}

# old version (until 2.0.1)
# fdata.bootstrap <- function(fdataobj, statistic  =func.mean
#                             ,alpha = 0.05, nb= 200, smo = 0
#                             ,draw = FALSE,draw.control = NULL,...){
#   if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
#   data<-fdataobj[["data"]]
#   if (smo > 0) varX <- var(data)
#   estmues<-statistic(fdataobj,...)
#   nr<-nrow(fdataobj)
#   nc<-ncol(fdataobj)
#   tt =fdataobj[["argvals"]]
#   rtt=fdataobj[["rangeval"]]
#   names=fdataobj[["names"]]
#   distboot<-matrix(NA,nrow=nb)
#   estboot<-matrix(NA,nrow=nb,ncol=nc)
#   pb=txtProgressBar(min=0,max=nb,style=3)
#   for (i in 1:nb){
#     setTxtProgressBar(pb,i-0.5)
#     bmuestra<-fdataobj[sample(1:nr,size=nr,replace=TRUE),]
#     if (smo>0) {
#       bmuestra[["data"]]<-bmuestra[["data"]]+mvrnorm(n=nr,rep(0,nc),varX*smo)
#     }
#     stat<-statistic(bmuestra,...)
#     estboot[i,]<-stat[["data"]]
#     setTxtProgressBar(pb,i)
#   }
#   close(pb)
#   center<-estmues
#   #for (i in 1:nb){  aux<-fdata(estboot[i,],tt,rtt)
#   #  distboot[i]<-metric.lp(center,aux,...)  }
#   #print(dim(distboot))
#   distboot<-metric.lp(fdata(estboot,tt,rtt),center,...)
#   dist<-max(distboot[rank(distboot)<=floor((1-alpha)*nb)])
#   resample<-fdata(estboot,tt,rtt,names)
#   if (draw){
#     if (is.null(draw.control)) draw.control=list("col"=c("grey","blue","pink"),"lty"=c(2,1,1),"lwd"=c(1,2,1))
#     if (is.null(draw.control$lwd)) draw.control$lwd=c(1,2,1)
#     if (is.null(draw.control$lty)) draw.control$lty=c(2,1,1)
#     if (is.null(draw.control$col)) draw.control$col=c("grey","blue","pink")
#     plot(fdataobj,lwd=draw.control$lwd[1],lty=draw.control$lty[1],col=draw.control$col[1])
#     # for(i in 1:nb){
#     #if (distboot[i]<=dist) lines(tt,estboot[i,],lwd=draw.control$lwd[3],lty=draw.control$lty[3],col=draw.control$col[3])
#     #else lines(tt,estboot[i,],lwd=draw.control$lwd[4],lty=draw.control$lty[4],
#     #col=draw.control$col[4])  }
#     lines(resample[distboot<=dist],lwd=draw.control$lwd[3],lty=draw.control$lty[3],col=draw.control$col[3])
#     lines(estmues,lwd=draw.control$lwd[2],lty=draw.control$lty[2],col=draw.control$col[2])
#     legend(x=min(tt),y=0.99*max(data),legend=c("original curves",stat$names$main,"bootstrap curves IN"),
#            lty=draw.control$lty,lwd=draw.control$lwd,col=draw.control$col,cex=0.9,box.col=0)
#   }
#   return(list("statistic"=estmues,"dband"= dist,"rep.dist"=distboot,
#               "resample"=resample,fdataobj=fdataobj))
# }
