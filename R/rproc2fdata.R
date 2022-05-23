#' @title Simulate several random processes.
#' 
#' @description Simulate Functional Data from different processes: Ornstein Uhlenbeck,
#' Brownian, Fractional Brownian, Gaussian or Exponential variogram.
#' 
#' @param n Number of functional curves to be generated.
#' @param t Discretization points.
#' @param mu \code{vector} which specifies the trend values at the
#' discretization points, by default \code{mu}=\eqn{\mu(t)=0}.  If \code{mu} is
#' a \code{fdata} class object, \code{t}\eqn{=}\code{argvals(mu)}.
#' @param sigma A positive-definite symmetric matrix,
#' \eqn{\Sigma_{s,t}}{\Sigma}, specifying the covariance matrix among grid
#' points.  If \code{sigma} is a \code{scalar}, creates a random Gaussian
#' process with \eqn{\Sigma_{s,t}=}{\Sigma=}\code{sigmaI} (by default
#' \code{sigma=1}).\cr If \code{sigma} is a \code{vector}, creates a random
#' Gaussian process with \eqn{\Sigma_{s,t}=}{\Sigma=}\code{diag(sigma)}.\cr If
#' \code{sigma} is a character: create a random process using the covariance
#' matrix \eqn{\Sigma_{s,t}}{\Sigma} indicated in the argument, \itemize{ \item
#' \code{"OU"} or \code{"OrnsteinUhlenbeck"}, creates a random Ornstein
#' Uhlenbeck process with
#' \eqn{\Sigma_{s,t}=\frac{\sigma^2}{2\theta}e^{-\theta\left(s+t\right)}
#' \left(e^{2\theta\left(s+t\right)}-1\right)}{\Sigma(s,t)=\sigma^2/\theta
#' exp(-\theta(s+t))(exp(2\theta(s+t)-1))}, by default
#' \eqn{\theta=1/(3range(t))}{\theta=1/(3range(t))},
#' \eqn{\sigma^2={1}}{\sigma^2=1}.  \item \code{"brownian"} or \code{"wiener"},
#' creates a random Wiener process with \eqn{\Sigma_{s,t}=\sigma^2
#' min(s,t)}{\Sigma(s,t)=\sigma^2 min(s,t)}, by default
#' \eqn{\sigma^2=1}{\sigma^2=1}. \item \code{"fbrownian"}, creates a random
#' fractional brownian process with
#' \eqn{\Sigma_{s,t}=\sigma^{2H}/2{|s|^{2H}+|t|^{2H}-|s-t|^{2H}}}{\Sigma(s,t)=\sigma^{2H}/2{|s|^{2H}+|t|^{2H}-|s-t|^{2H}}},
#' by default \eqn{\sigma^2=1}{\sigma^2=1} and \eqn{H=0.5}{H=0.5} (brownian
#' process).  \item \code{"vexponential"}, creates a random gaussian process
#' with exponential variogram \eqn{\Sigma_{s,t}=\sigma^2
#' e^{\left(-\frac{\left|s-t\right|}{\theta}\right)}}{\Sigma=\sigma^2
#' exp(-|s-t|/\theta )}, by default \eqn{\theta={0.2
#' range(t)}}{\theta=.2*range(t)}, \eqn{\sigma^2={1}}{\sigma^2=1}. }
#' @param par.list List of parameter to process, by default \code{"scale"}
#' \eqn{\sigma^2=1}{\sigma^2=1}, \code{"theta"} \eqn{\theta=0.2
#' range(t)}{\theta=0.2 range(t)} and \code{"H"}=0.5.
#' @param norm If \code{TRUE} the norm of random projection is 1. Default is
#' \code{FALSE}
#' @param verbose If \code{TRUE}, information about procedure is printed.
#' Default is \code{FALSE}.
#' @param \dots Further arguments passed to or from other methods.
#' 
#' @return Return the functional random processes as a \code{fdata} class
#' object.
#' 
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' 
#' @keywords datagen
#' 
#' @examples
#' \dontrun{
#' par(mfrow=c(3,2))
#' lent<-30
#' tt<-seq(0,1,len=lent)
#' mu<-fdata(rep(0,lent),tt)
#' plot(rproc2fdata(200,t=tt,sigma="OU",par.list=list("scale"=1)))
#' plot(rproc2fdata(200,mu=mu,sigma="OU",par.list=list("scale"=1)))
#' plot(rproc2fdata(200,t=tt,sigma="vexponential"))
#' plot(rproc2fdata(200,t=tt,sigma=1:lent))
#' plot(rproc2fdata(200,t=tt,sigma="brownian"))
#' plot(rproc2fdata(200,t=tt,sigma="wiener"))
#' #plot(rproc2fdata(200,seq(0,1,len=30),sigma="oo")) # this is an error 
#' }
#' 
#' @rdname rproc2fdata
#' @export
rproc2fdata=function(n, t=NULL, mu=rep(0,length(t)), sigma=1,
                     par.list=list("scale"=1,"theta"=.2*diff(rtt),"H"=.5),
                     norm=FALSE, verbose=FALSE,...) {
sigma2<-sigma
if (is.fdata(mu)){   
  if (!is.null(t) & verbose) warnings("The argvals of argument t are ignored, the function uses the argvals(mu)")
  t<-mu$argvals
  p=length(t)   
  rtt<-mu$rangeval
  mu<-drop(mu$data[1,])
 } 
else {
  if (is.null(t)) stop("At least t or fdata class mu must be provided")
  p=length(t)
  if (p!= length(mu)) stop("t and mu must have the same length")
  rtt=range(t)
  }
if (is.null(par.list$mu0)) par.list$mu0<-0
if (is.null(par.list$theta)) par.list$theta<-1/(3*(diff(rtt)))
#if (is.null(par.list$theta)) par.list$theta<-1/3
if (is.null(par.list$scale)) par.list$scale<-1
if (is.character(sigma)) {
 type.proc<-c("brownian","wiener","OU","OrnsteinUhlenbeck","vexponential","fbrownian")
 if (!is.element(sigma,type.proc)) stop("Error in sigma argument label")
 ss<-which(sigma==type.proc)
   sigma2<-sigma
  sigma=switch(sigma,"brownian"=par.list$scale*outer(t,t,function(u,v){pmin(u,v)}),
         "wiener"= par.list$scale*outer(t,t,function(u,v){pmin(u,v)}),
         "OU"= {#m<-100;t2<-t+m
         par.list$scale/(2*par.list$theta)*outer(t,t,function(u,v){
            exp(-par.list$theta*(u+v))*(exp(2*par.list$theta*pmin(u,v))-1)})},
"OrnsteinUhlenbeck"={#m<-100;#t2<-t+m
            par.list$scale/(2*par.list$theta)*outer(t,t,function(u,v){
            exp(-par.list$theta*(u+v))*(exp(2*par.list$theta*pmin(u,v))-1)})},
             "vexponential"=par.list$scale*outer(t,t,function(u,v){
             exp(-abs(u-v)/par.list$theta)}),
			"fbrownian"=0.5*abs(par.list$scale)^par.list$H*outer(t,t,function(u,v){abs(u)^(2*par.list$H)+abs(v)^(2*par.list$H)-abs(u-v)^(2*par.list$H)})
             )
  sigma<-t(sigma)

}
else {
 sigma2<-"Gaussian"
 if   (is.matrix(sigma)) {
  if (dim(sigma)[2]!=p) stop("Error in sigma argument")
 }
 else if (length(sigma)==1 | length(sigma)==p) sigma<-diag(p)*sigma
 else stop("Error in sigma argument")
}
C=svd(t(sigma))
L=C$u%*%diag(sqrt(C$d))
X=matrix(rnorm(n*p),ncol=p)
X=t(L%*%t(X))
X=sweep(X,2,mu,"+")
X=fdata(X,t,rtt,names=list("main"=paste(sigma2," process",sep="")))
if (norm) {
if (sigma2[1]=="brownian") print("The normalization is not done")
else{
        no <- norm.fdata(X,...)
        X$data <- sweep(X$data, 1, no, "/")
        }
}
return(X)
}

