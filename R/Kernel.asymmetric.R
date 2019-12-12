#' Asymmetric Smoothing Kernel
#' 
#' @description Represent Asymmetric Smoothing Kernels: normal, cosine, triweight, quartic
#' and uniform. \tabular{ll}{ \tab AKer.norm=ifelse(u>=0,2*dnorm(u),0)\cr \tab
#' AKer.cos=ifelse(u>=0,pi/2*(cos(pi*u/2)),0)\cr \tab AKer.epa=ifelse(u>=0 &
#' u<=1,3/2*(1-u^2),0)\cr \tab AKer.tri=ifelse(u>=0 &
#' u<=1,35/16*(1-u^2)^3,0)\cr \tab AKer.quar=ifelse(u>=0 &
#' u<=1,15/8*(1-u^2)^2,0)\cr \tab AKer.unif=ifelse(u>=0 & u<=1,1,0)\cr }
#' 
#' @details Type of Asymmetric kernel: \tabular{ll}{ \tab Asymmetric Normal Kernel:
#' \code{AKer.norm}\cr \tab Asymmetric Cosine Kernel: \code{AKer.cos}\cr \tab
#' Asymmetric Epanechnikov Kernel: \code{AKer.epa}\cr \tab Asymmetric Triweight
#' Kernel: \code{AKer.tri}\cr \tab Asymmetric Quartic Kernel:
#' \code{AKer.quar}\cr \tab Asymmetric Uniform Kernel: \code{AKer.unif}\cr }
#' 
#' @aliases AKer.norm AKer.cos AKer.epa AKer.tri AKer.tri AKer.quar AKer.unif
#' Kernel.asymmetric
#' @param type.Ker Type of asymmetric metric kernel, by default asymmetric
#' normal kernel.
#' @param u Data.
#' @return Returns asymmetric kernel.
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@usc.es}
#' @references Ferraty, F. and Vieu, P. (2006). \emph{Nonparametric functional
#' data analysis.} Springer Series in Statistics, New York.
#' 
#' Hardle, W. \emph{Applied Nonparametric Regression}. Cambridge University
#' Press, 1994.
#' @keywords kernel
#' @examples
#' y=qnorm(seq(.1,.9,len=100))
#' a<-Kernel.asymmetric(u=y)
#' b<-Kernel.asymmetric(type.Ker="AKer.tri",u=y)
#' c=AKer.cos(y)
#' @export Kernel.asymmetric
Kernel.asymmetric=function(u,type.Ker="AKer.norm"){
   tab=list("AKer.norm","AKer.cos","AKer.epa","AKer.tri","AKer.quar","AKer.unif")
   type.i=pmatch(type.Ker,tab)
   if (is.na(type.i))   Ker=ifelse(u>=0,2*dnorm(u),0)
     else {
     if (type.i==1)   Ker=ifelse(u>=0,2*dnorm(u),0)
		 if (type.i==2)   Ker=ifelse(u>=0,pi/2*(cos(pi*u/2)),0)
     if (type.i==3)   Ker=ifelse(u>=0 & u<=1,1.5*(1-u^2),0)
     if (type.i==4)   Ker=ifelse(u>=0 & u<=1,35/16*(1-u^2)^3,0)
     if (type.i==5)   Ker=ifelse(u>=0 & u<=1,15/8*(1-u^2)^2,0)
     if (type.i==6)   Ker=ifelse(u>=0 & u<=1,1,0)
}
  return(Ker)
}

