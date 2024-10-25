#' @title Conditional mode
#' 
#' @description Computes the mode for conditional distribution function.
#' 
#' @details The conditional mode is calculated as the maximum argument of the derivative
#' of the conditional distribution function (density function \code{f}).
#' 
#' @param Fc Object estimated by \code{cond.F} function.
#' @param method Specifies the type of spline to be used. Possible values are
#' \emph{"diff"}, \emph{"fmm"}, \emph{"natural"}, \emph{"periodic"} and
#' \emph{"monoH.FC"}.
#' @param draw =TRUE, plots the conditional distribution and density function.
#' 
#' @return Return the mode for conditional distribution function.
#' \itemize{
#' \item {\code{mode.cond}}{ Conditional mode.} 
#' \item {\code{x}}{ Grid of length \code{n} where the the conditional density function is evaluated.} 
#' \item {\code{f}}{ The conditional density function evaluated in \code{x}.}
#' }
#' 
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' 
#' @seealso See Also as: \code{\link{cond.F}}, \code{\link{cond.quantile}} and
#' \link[stats]{splinefun} .
#' 
#' @references Ferraty, F. and Vieu, P. (2006). \emph{Nonparametric functional
#' data analysis.} Springer Series in Statistics, New York.
#' 
#' @keywords distribution
#' @examples
#' \dontrun{
#' n= 500
#' t= seq(0,1,len=101)
#' beta = t*sin(2*pi*t)^2
#' x = matrix(NA, ncol=101, nrow=n)
#' y=numeric(n)
#' x0<-rproc2fdata(n,seq(0,1,len=101),sigma="wiener")
#' x1<-rproc2fdata(n,seq(0,1,len=101),sigma=0.1)
#' x<-x0*3+x1
#' fbeta = fdata(beta,t)
#' y<-inprod.fdata(x,fbeta)+rnorm(n,sd=0.1)
#' prx=x[1:100];pry=y[1:100]
#' ind=101;ind2=101:110
#' pr0=x[ind];pr10=x[ind2]
#' ndist=161
#' gridy=seq(-1.598069,1.598069, len=ndist)
#' # Conditional Function
#' I=5
#' # Time consuming
#' res = cond.F(pr10[I], gridy, prx, pry, h=1)
#' mcond=cond.mode(res)
#' mcond2=cond.mode(res,method="diff")
#' }
#' 
#' @export
cond.mode=function (Fc, method = "monoH.FC",draw=TRUE) 
{
    x = Fc$y0
    Fc = Fc$Fc
    ndist = length(x)
    if (draw) par(mfrow = c(2, 1))
    if (method == "diff") {
		if (draw) {
        plot(x, Fc, type = "l", main = "Conditional Distribution Function", 
            ylab = "Fc")}
        if (is.vector(Fc)) 
            Fc = matrix(Fc, ncol = 1)
        ind3 = 1
        der = (Fc[2:ndist, ind3] - Fc[1:(ndist - 1), ind3])/(x[2:ndist] - 
            x[1:(ndist - 1)])
        x<-x[2:(ndist)]    
        names(der) = round(x, 2)
    }
    else {
        f <- splinefun(x, Fc, method = method)
        if (draw) {curve(f(x), x[1], x[length(Fc)], main = "Conditional Distribution Function", 
            ylab = "Fc")}
        der = f(x, deriv = 1)
    }
    max.der = max(der)
    ind.mode.cond = which.max(der)
    mode.cond = x[ind.mode.cond]
#    print((x));    print((der))
	if (draw) {plot(x,der, type = "l", main = c("Conditional mode=", round(mode.cond, 
        3)), xlab = "y",ylab = "Derivative")
    abline(v = x[ind.mode.cond], col = 2, lty = 2)}
    return(list(mode.cond = mode.cond, x = x,f=der))
}

