#' Integrate Smoothing Kernels.
#' 
#' @description Represent integrate kernels: normal, cosine, triweight, quartic and uniform.
#' 
#' @details Type of integrate kernel: \tabular{ll}{ \tab Integrate Normal Kernel:
#' \code{IKer.norm}\cr \tab Integrate Cosine Kernel: \code{IKer.cos}\cr \tab
#' Integrate Epanechnikov Kernel: \code{IKer.epa}\cr \tab Integrate Triweight
#' Kernel: \code{IKer.tri}\cr \tab Integrate Quartic Kernel:
#' \code{IKer.quar}\cr \tab Integrate Uniform Kernel: \code{IKer.unif}\cr }
#' 
#' @aliases Kernel.integrate IKer.norm IKer.cos IKer.epa IKer.tri IKer.tri
#' IKer.quar IKer.unif
#' @param u data
#' @param Ker Type of Kernel. By default normal kernel.
#' @param a Lower limit of integration.
#' @return Returns integrate kernel.
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' @seealso See Also as: \code{\link{Kernel}} and \link[stats]{integrate}.
#' @references Ferraty, F. and Vieu, P. (2006). \emph{Nonparametric functional
#' data analysis.} Springer Series in Statistics, New York.
#' 
#' Hardle, W. \emph{Applied Nonparametric Regression}. Cambridge University
#' Press, 1994.
#' @keywords kernel
#' @examples
#' y=qnorm(seq(.1,.9,len=100))
#' d=IKer.tri(y)
#' e=IKer.cos(y)
#' e2=Kernel.integrate(u=y,Ker=Ker.cos)
#' e-e2
#' f=IKer.epa(y)
#' f2=Kernel.integrate(u=y,Ker=Ker.epa)
#' f-f2
#' plot(d,type="l",ylab="Integrate Kernel")
#' lines(e,col=2,type="l")
#' lines(f,col=4,type="l") 
#' 
#' @export Kernel.integrate
Kernel.integrate=function(u,Ker=Ker.norm,a=-1){
    kaux=function(u){integrate(Ker,a,u)$value}
    mapply(kaux,u)
}


#' @export AKer.norm
AKer.norm=function(u)  {ifelse(u>=0,2*dnorm(u),0)}
#' @export AKer.cos
AKer.cos=function(u)  {ifelse(u>=0,pi/2*(cos(pi*u/2)),0)}
#' @export AKer.epa
AKer.epa=function(u) {ifelse(u>=0 & u<=1,1.5*(1-u^2),0)}
#' @export AKer.tri
AKer.tri=function(u) {ifelse(u>=0 & u<=1,35/16*(1-u^2)^3,0)}
#' @export AKer.quar
AKer.quar=function(u) {ifelse(u>=0 & u<=1,15/8*(1-u^2)^2,0)}
#' @export AKer.unif
AKer.unif=function(u) {ifelse(u>=0 & u<=1,1,0)}
#' @export IKer.norm
IKer.norm=function(u){
	kaux=function(u){integrate(Ker.norm,-1,u)$value}
	mapply(kaux,u)
}
#' @export IKer.cos
IKer.cos=function(u){
	kaux=function(u){integrate(Ker.cos,-1,u)$value}
	mapply(kaux,u)
}
#' @export IKer.tri
IKer.tri=function(u){
	kaux=function(u){integrate(Ker.tri,-1,u)$value}
	mapply(kaux,u)
}
#' @export IKer.quar
IKer.quar=function(u){
	kaux=function(u){integrate(Ker.quar,-1,u)$value}
	mapply(kaux,u)
}
#' @export IKer.unif
IKer.unif=function(u){
	kaux=function(u){integrate(Ker.unif,-1,u)$value}
	mapply(kaux,u)
}
#' @export IKer.epa
IKer.epa=function(u){
	kaux=function(u){integrate(Ker.epa,-1,u)$value}
	mapply(kaux,u)
}



#
rkernel<-function(x,y,Ker=Ker.norm){
  s1=sum(Ker(x)*y,na.rm=TRUE)
  ss=sum(Ker(x),na.rm=TRUE)
  return(s1/ss)}
