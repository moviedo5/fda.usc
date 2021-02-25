#' @title Symmetric Smoothing Kernels.
#' 
#' @description Represent symmetric smoothing kernels:: normal, cosine, triweight, quartic
#' and uniform.
#' 
#' @details 
#' \tabular{ll}{ \tab Ker.norm=dnorm(u)\cr \tab
#' Ker.cos=ifelse(abs(u)<=1,pi/4*(cos(pi*u/2)),0)\cr \tab
#' Ker.epa=ifelse(abs(u)<=1,3/4*(1-u^2),0)\cr \tab
#' Ker.tri=ifelse(abs(u)<=1,35/32*(1-u^2)^3,0)\cr \tab
#' Ker.quar=ifelse(abs(u)<=1,15/16*(1-u^2)^2,0)\cr \tab
#' Ker.unif=ifelse(abs(u)<=1,1/2,0)\cr }
#' 
#' Type of kernel: \tabular{ll}{ \tab Normal Kernel: \code{Ker.norm}\cr \tab
#' Cosine Kernel: \code{Ker.cos}\cr \tab Epanechnikov Kernel: \code{Ker.epa}\cr
#' \tab Triweight Kernel: \code{Ker.tri}\cr \tab Quartic Kernel:
#' \code{Ker.quar}\cr \tab Uniform Kernel: \code{Ker.unif}\cr }
#' 
#' @aliases Kernel Ker.norm Ker.cos Ker.epa Ker.tri Ker.tri Ker.quar Ker.unif
#' @param type.Ker Type of Kernel. By default normal kernel.
#' @param u Data.
#' @return  Returns symmetric kernel.
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' @references Ferraty, F. and Vieu, P. (2006). \emph{Nonparametric functional
#' data analysis.} Springer Series in Statistics, New York. \cr
#' 
#' Hardle, W. \emph{Applied Nonparametric Regression}. Cambridge University
#' Press, 1994.
#' @keywords kernel
#' @examples
#' y=qnorm(seq(.1,.9,len=100))
#' a<-Kernel(u=y)
#' b<-Kernel(type.Ker="Ker.tri",u=y)
#' c=Ker.cos(y)
#' @export Kernel
Kernel=function(u,type.Ker="Ker.norm"){
  tab=list("Ker.norm","Ker.cos","Ker.epa","Ker.tri","Ker.quar","Ker.unif")
  type.i=pmatch(type.Ker,tab)
  if (is.na(type.i))   Ker=dnorm(u)
  else {
    if (type.i==1)   Ker=dnorm(u)
    if (type.i==2)   Ker=ifelse(abs(u)<=1,pi/4*(cos(pi*u/2)),0)
    if (type.i==3)   Ker=ifelse(abs(u)<=1,0.75*(1-u^2),0)
    if (type.i==4)   Ker=ifelse(abs(u)<=1,35/32*(1-u^2)^3,0)
    if (type.i==5)   Ker=ifelse(abs(u)<=1,15/16*(1-u^2)^2,0)
    if (type.i==6)   Ker=ifelse(abs(u)<=1,0.5,0)
  }
  return(Ker)
}



#' @export Ker.norm
Ker.norm=function(u){dnorm(u)}
#' @export Ker.cos
Ker.cos=function(u){ifelse(abs(u)<=1,pi/4*(cos(pi*u/2)),0)}
#' @export Ker.epa
Ker.epa=function(u){ifelse(abs(u)<=1,0.75*(1-u^2),0)}
#' @export Ker.tri
Ker.tri=function(u){ifelse(abs(u)<=1,35/32*(1-u^2)^3,0)}
#' @export Ker.quar
Ker.quar=function(u){ifelse(abs(u)<=1,15/16*(1-u^2)^2,0)}
#' @export Ker.unif
Ker.unif=function(u){ifelse(abs(u)<=1,0.5,0)}