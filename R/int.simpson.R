###################################
# 20190406 equi argument is deleted
# 20190406 argvals.equi function is created and used
###################################


#' @title Simpson integration
#' 
#' @description Computes the integral of \code{fdataobj$data} with respect to
#' \code{fdataobj$argvals} using simpson or trapezoid rule integration.
#' 
#' @details Posible values for \code{method} are: \itemize{ \item \code{"TRAPZ"}:
#' Trapezoid rule integration. \item \code{"CSR"}: Composite Simpson's rule
#' integration.  \item \code{"ESR"}: Extended Simpson's rule integration. } If
#' \code{method=NULL} (default), the value of \code{par.fda.usc$int.method} is
#' used.
#' 
#' @aliases int.simpson int.simpson2
#' @param fdataobj fdata objtect.
#' @param x Sorted vector of x-axis values: \code{argvals}.
#' @param y Vector of y-axis values.
#' @param equi =TRUE, the observed points on each curve are equally spaced (by
#' default).
#' @param method Method for numerical integration, see details.
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' @seealso See also \code{\link{integrate}}.
#' @keywords cluster
#' @examples
#' \dontrun{
#' x<-seq(0,2*pi,length=1001)
#' fx<-fdata(sin(x)/sqrt(pi),x)
#' fx0<-fdata(rep(0,length(x)),x)
#' int.simpson(fx0)
#' int.simpson(fx)
#' }
#' @export
int.simpson=function(fdataobj,method=NULL){
 par.fda.usc <-      eval(parse(text="fda.usc:::par.fda.usc"), envir=.GlobalEnv)
 method <- par.fda.usc$int.method
 if (!inherits(fdataobj, "fdata"))  {
   fdataobj = fdata(fdataobj)
   equi = TRUE
   tt = fdataobj$argvals
 } else {
   tt = fdataobj$argvals
   equi = argvals.equi(tt)
 }
 n = NROW(fdataobj)
 out = rep(NA,n)
 for (i in 1:n) {
   out[i] <- int.simpson2(tt,fdataobj$data[i,],equi=equi,method=method)
 }
	return(out)
}

#' @rdname int.simpson
#' @export int.simpson2
int.simpson2=function(x,y,equi=TRUE,
                      method=NULL){
  n=length(x);ny=length(y)
  if (is.null(method)) {
    par.fda.usc <-      eval(parse(text="fda.usc:::par.fda.usc"), envir=.GlobalEnv)
    method <- par.fda.usc$int.method
  } 
  if (n!=ny) stop("Different length in the input data")
  if (n==2 || ny==2) method="TRAPZ"
  out <- switch(method,
                "TRAPZ" = {
                  if (!equi){
                    idx=2:n
                    value<-as.double((x[idx]-x[idx-1])%*%(y[idx]+y[idx-1]))/2
                  } else {
                    h=(max(x)-min(x))/(n-1)
                    y[c(1,n)]=y[c(1,n)]/2
                    value<-h*sum(y)
                  }
                  value
                },"CSR" = {
                  if (!equi){
                    n=2*n-1
                    app=approx(x,y,n=n);x=app$x;y=app$y}
                  h=(max(x)-min(x))/(n-1)
                  value=(h/3)*(y[n]+y[1]+2*sum(y[2*(1:((n-1)/2))+1])+4*sum(y[2*(1:((n-1)/2))]))
                },
                "ESR" = {
                  if (!equi){
                    n = 2*n-1
                    app = approx(x,y,n=n)
                    x=app$x
                    y=app$y
                  }
                  h = (max(x)-min(x))/(n-1)
                  if (n<=4) stop("This method needs n>4")
                  value=17*(y[1]+y[n])+59*(y[2]+y[n-1])+43*(y[3]+y[n-2])+49*(y[4]+y[n-3])
                  value=value+48*sum(y[5:(n-4)])
                  value=(h/48)*value
                }
  )
  return(out)
}