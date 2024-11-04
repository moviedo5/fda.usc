#' @name dfv.test
#' @aliases  dfv.test dfv.statistic
#' @title Delsol, Ferraty and Vieu test for no functional-scalar interaction
#' 
#' @description  The function \code{dfv.test} tests the null hypothesis of no interaction between a functional covariate and a scalar response in a general framework. The null hypothesis is \deqn{H_0:\,m(X)=0,}{H_0: m(X)=0,} where \eqn{m(\cdot)}{m(.)} denotes the regression function of the functional variate \eqn{X}{X} over the centred scalar response \eqn{Y}{Y} (\eqn{E[Y]=0}{E[Y]=0}). The null hypothesis is tested by the smoothed integrated square error of the response (see Details). 
#' 
#' @param X.fdata Functional covariate. The object must be in the class \code{\link{fdata}}.
#' @param Y Scalar response. Must be a vector with the same number of elements as functions are in \code{X.fdata}.
#' @param h Bandwidth parameter for the kernel smoothing. This is a crucial parameter that affects the power performance of the test. One possibility to choose it is considering the Cross-validatory bandwidth of the nonparametric functional regression, given by the function \code{\link{fregre.np}} (see Examples). 
#'     Other possibility is to consider a grid of bandwidths. This is the default option, considering the grid given by the quantiles 0.05, 0.10, 0.15, 0.25 and 0.50 of the functional \eqn{L^2}{L^2} distances of the data.
#' @param B Number of bootstrap replicates to calibrate the distribution of the test statistic. \code{B=5000} replicates are the recommended for carry out the test, although for exploratory analysis (\bold{not inferential}), an acceptable less time-consuming option is \code{B=500}.
#' @param K Kernel function. If no specified it is taken to be the rescaled right part of the normal density.
#' @param weights A vector of weights for the sample data. The default is the uniform weights \code{rep(1,dim(X.fdata$data)[1])}.
#' @param d Semimetric to use in the kernel smoothers. By default is the \eqn{L^2}{L^2} distance given by \code{\link{metric.lp}}.
#' @param dist Matrix of distances of the functional data, used to save time in the bootstrap calibration. If not given, the matrix is automatically computed using the semimetric \code{d}.
#' @param verbose Either to show or not information about computing progress.
#' @details  The Delsol, Ferraty and Vieu statistic is defined as 
#' \deqn{T_n=\int\bigg(\sum_{i=1}^n(Y_i-m(X_i))K\bigg(\frac{d(X,X_i)}{h}\bigg)\bigg)^2\omega(X)dP_X(X)}{T_n=\int(\sum_{i=1}^n(Y_i-m(X_i))K(\frac{d(X,X_i)}{h}))^2\omega(X)dP_X(X)}
#' and in the case of no interaction with centred scalar response (when \eqn{H_0:\,m(X)=0}{H_0: m(X)=0} 
#' holds), its sample version is computed from 
#' \deqn{T_n=\frac{1}{n}\sum_{j=1}^n\bigg(\sum_{i=1}^n Y_iK\bigg(\frac{d(X_j,X_i)}{h}\bigg)\bigg)^2\omega(X_j).}{T_n=\frac{1}{n}\sum_{j=1}^n(\sum_{i=1}^n Y_iK(\frac{d(X_j,X_i)}{h}))^2\omega(X_j).}
#' The sample version implemented here does not consider a splitting of the sample, as the authors 
#' comment in their paper. The statistic is computed by the function \code{dfv.statistic} and, before
#' applying the test, the response \eqn{Y}{Y} is centred. The distribution of the test statistic
#' is approximated by a wild bootstrap on the residuals, using the \emph{golden section bootstrap}.
#' 
#' Please note that if a grid of bandwidths is passed, a harmless warning message will prompt at the
#' end of the test (it comes from returning several p-values in the \code{htest} class).
#' 
#' @return  The value of \code{dfv.statistic} is a vector of length \code{length(h)} with the values
#'  of the statistic for each bandwidth. The value of \code{dfv.test} is an object with class 
#'  \code{"htest"} whose underlying structure is a list containing the following components:
#' \itemize{
#' \item \code{statistic}: The value of the Delsol, Ferraty and Vieu test statistic.
#' \item \code{boot.statistics}: A vector of length \code{B} with the values of the bootstrap test statistics.
#' \item \code{p.value}: The p-value of the test.
#' \item \code{method}: The character string "Delsol, Ferraty and Vieu test for no functional-scalar interaction".
#' \item \code{B}: The number of bootstrap replicates used.
#' \item \code{h}: Bandwidth parameters for the test.
#' \item \code{K}: Kernel function used.
#' \item \code{weights}: The weights considered.
#' \item \code{d}: Matrix of distances of the functional data.
#' \item \code{data.name}: The character string "Y=0+e".
#' }
#' 
#' @references
#' Delsol, L., Ferraty, F. and Vieu, P. (2011). Structural test in regression on functional variables. 
#' Journal of Multivariate Analysis, 102, 422-447. \doi{10.1016/j.jmva.2010.10.003}
#' 
#' Delsol, L. (2013). No effect tests in regression on functional variable and some applications to spectrometric studies. Computational Statistics, 28(4), 1775-1811. \doi{10.1007/s00180-012-0378-1}
#'  
#' @note No NA's are allowed neither in the functional covariate nor in the scalar response.
#' @author Eduardo Garcia-Portugues. Please, report bugs and suggestions
#'  to \if{latex}{\cr}\email{eduardo.garcia.portugues@uc3m.es}
#'  
#' @seealso \code{\link{rwild}}, \code{\link{flm.test}}, \code{\link{flm.Ftest}}, \code{\link{fregre.np}}
#' @examples
#' \dontrun{
#' ## Simulated example ##
#' X=rproc2fdata(n=50,t=seq(0,1,l=101),sigma="OU")
#' 
#' beta0=fdata(mdata=rep(0,length=101)+rnorm(101,sd=0.05),
#' argvals=seq(0,1,l=101),rangeval=c(0,1))
#' beta1=fdata(mdata=cos(2*pi*seq(0,1,l=101))-(seq(0,1,l=101)-0.5)^2+
#' rnorm(101,sd=0.05),argvals=seq(0,1,l=101),rangeval=c(0,1))
#' 
#' # Null hypothesis holds
#' Y0=drop(inprod.fdata(X,beta0)+rnorm(50,sd=0.1))
#' 
#' # Null hypothesis does not hold
#' Y1=drop(inprod.fdata(X,beta1)+rnorm(50,sd=0.1))
#' 
#' # We use the CV bandwidth given by fregre.np
#' # Do not reject H0
#' dfv.test(X,Y0,h=fregre.np(X,Y0)$h.opt,B=100)
#' # dfv.test(X,Y0,B=5000)
#' 
#' # Reject H0
#' dfv.test(X,Y1,B=100)
#' # dfv.test(X,Y1,B=5000)
#' }
#' 
#' @keywords htest models regression

##############################################
##   Delsol, Ferraty and Vieu test for no   ##
##      functional-scalar interaction       ##
##############################################

##############################################
## File created by Eduardo Garcia-Portugues ##
## using code from library fda.usc          ##
##############################################


# Delsol, Ferraty and Vieu test statisitc for the simple hypothesis of no interaction
#' @rdname dfv.test
#' @export 
dfv.statistic=function(X.fdata,Y,h=quantile(x=metric.lp(X.fdata),probs=c(0.05,0.10,0.15,0.25,0.50)),K=function(x)2*dnorm(abs(x)),weights=rep(1,dim(X.fdata$data)[1]),d=metric.lp,dist=NULL){
	
	# Check if it is necesary to compute the distances matrix
	if(is.null(dist)) dist=d(X.fdata,X.fdata)
	
	# Statistic
	res=sapply(h,function(hh) mean(drop(K(dist/hh)%*%Y)^2*weights))

	return(res)

}

# Delsol, Ferraty and Vieu test for the simple hypothesis of no interaction with bootstrap calibration
#' @rdname dfv.test
#' @export 
dfv.test=function(X.fdata,Y,B=5000,h=quantile(x=metric.lp(X.fdata),probs=c(0.05,0.10,0.15,0.25,0.50)),K=function(x)2*dnorm(abs(x)),weights=rep(1,dim(X.fdata$data)[1]),d=metric.lp,verbose=TRUE){
	
	# REAL WORLD
	dist=d(X.fdata,X.fdata)
	Tn=dfv.statistic(X.fdata=X.fdata,Y=Y,h=h,K=K,weights=weights,dist=dist)
	
	# BOOTSTRAP WORLD
	Tn.star=matrix(nrow=B,ncol=length(h))
	colnames(Tn.star)=paste("Tn.boot(",sprintf("h=%.3f",h),")",sep="")
	
	if(verbose) pb=txtProgressBar(style=3)
	for(i in 1:B){
		
		# Bootsrtap version of Tn: perturbation with a centred and unit variance noise
		Tn.star[i,]=dfv.statistic(X.fdata=X.fdata,Y=rwild(Y,"golden"),h=h,K=K,weights=weights,dist=dist)
		
		# Progress bar
		if(verbose) setTxtProgressBar(pb,i/B)
		
	}
	
	# P-value
	pvalue=sapply(1:length(h),function(i) sum(Tn.star[,i]>Tn[i]))/B
	names(pvalue)=paste("p-value(",sprintf("h=%.3f",h),")",sep="")

	# Result: class htest
	names(Tn)=paste("Tn(",sprintf("h=%.3f",h),")",sep="")
	result=list(statistic=Tn,boot.statistics=Tn.star,p.value=pvalue,method="Delsol, Ferraty and Vieu test for no functional-scalar interaction",B=B,h=h,K=K,weights=weights,d=d,data.name="Y=0+e")
	class(result) <- "htest"
	return(result)
	
}
