#' @title Utils for generate functional data
#' 
#' @description \code{gridfdata} generates \code{n} curves as lineal combination of the
#' original curves \code{fdataobj} plus a functional trend \code{mu}.
#'  
#' @details \code{rcombfdata} generates \code{n} random linear combinations of the
#' \code{fdataobj} curves plus a functional trend \code{mu}. The coefficients
#' of the combinations follows a normal distribution with zero mean and
#' standard deviation \code{sdarg}. 
#' 
#' @aliases gridfdata rcombfdata
#' @param coef Coefficients of the combination. A matrix with number of columns
#' equal to number of curves in \code{fdataobj}
#' @param fdataobj \code{\link{fdata}} class object.
#' @param mu Functional trend, by default \code{mu}=\eqn{\mu(t)=0}. An object
#' of class \code{\link{fdata}}.  %If \code{mu} is a \code{fdata} class object,
#' \code{t}\eqn{=}\code{argvals(mu)}.
#' @param n Number of curves to be generated
#' @param sdarg Standard deviation of the coefficients.
#' @param norm Norm of the coefficients. The norm is adjusted before the
#' transformation for \code{sdarg} is performed.
#' @return Return the functional trajectories as a \code{fdata} class object.
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' @seealso See Also as \code{\link{rproc2fdata}}
#' @keywords generation
#' @examples
#' \dontrun{
#' tt=seq(0,1,len=51)
#' fou3=create.fourier.basis(c(0,1),nbasis=3)
#' fdataobj=fdata(t(eval.basis(tt,fou3)),argvals=tt)
#' 
#' coef=expand.grid(0,seq(-1,1,len=11),seq(-1,1,len=11))
#' grid=gridfdata(coef,fdataobj)
#' plot(grid,lty=1)
#' 
#' rcomb=rcombfdata(n=51,fdataobj,mu=fdata(30*tt*(1-tt),tt))
#' plot(rcomb,lty=1)
#' }
#' 
#' @rdname rcombfdata
#' @export
rcombfdata=function(n = 10, fdataobj, mu,
                    sdarg = rep(1,nrow(fdataobj)),
                    norm = 1){
  if (class(fdataobj)!="fdata") 
    stop("Argument fdataobj must be of class fdata")
  tt <- argvals(fdataobj)
  if (missing(mu)) 
    mu=fdata(rep(0,ncol(fdataobj)),argvals=tt)
  nr=nrow(fdataobj)
  xx=matrix(rnorm(n*nr),ncol=nr)
  xx=sweep(xx,1,apply(xx,1,function(v){sqrt(sum(v^2))})/norm,"/")
  xx=sweep(xx,2,sdarg,"*")
  res=xx%*%fdataobj[["data"]]
  res=fdata(sweep(res,2,mu[["data"]],"+"),argvals=tt)
  return(res)
}

#' @rdname rcombfdata
#' @export
gridfdata=function(coef,fdataobj,mu){
  if (class(fdataobj)!="fdata") stop("Argument fdataobj must be of class fdata")
  nr=nrow(fdataobj)
  tt <- argvals(fdataobj)
  if (missing(mu)) 
    mu=fdata(rep(0,ncol(fdataobj)),argvals=tt)  
  coef=data.matrix(coef)
  if (ncol(coef)!=nr) stop("Argument coef must be a matrix with ncol(coef)==nrow(fdataobj)")
  res=coef%*%fdataobj[["data"]]
  res=fdata(sweep(res,2,mu[["data"]],"+"),argvals=tt)
  return(res)
}

