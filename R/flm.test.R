#' @aliases  flm.test
#' 
#' @title Goodness-of-fit test for the Functional Linear Model with scalar response
#' 
#' @description The function \code{flm.test} tests the composite null hypothesis of
#' a Functional Linear Model with scalar response (FLM),
#' \deqn{H_0:\,Y=\big<X,\beta\big>+\epsilon,}{H_0: Y=<X,\beta>+\epsilon,}  versus
#' a general alternative. If \eqn{\beta=\beta_0}{\beta=\beta_0} is provided, then the 
#' simple hypothesis \eqn{H_0:\,Y=\big<X,\beta_0\big>+\epsilon}{H_0: Y=<X,\beta_0>+\epsilon} is tested.
#' The testing of the null hypothesis is done by a Projected Cramer-von Mises statistic (see Details). 
#' 
#' @param X.fdata Functional covariate for the FLM. The object must be in the class 
#' \code{\link{fdata}}.
#' @param Y Scalar response for the FLM. Must be a vector with the same number of elements
#'  as functions are in \code{X.fdata}.
#' @param beta0.fdata Functional parameter for the simple null hypothesis, in the \code{\link{fdata}} class. 
#' Recall that the \code{argvals} and \code{rangeval} arguments of \code{beta0.fdata} must be the same
#' of \code{X.fdata}. A possibility to do this is to consider, for example for \eqn{\beta_0=0}{\beta_0=0} 
#' (the simple null hypothesis of no interaction),\if{latex}{\cr}
#'  \code{beta0.fdata=fdata(mdata=rep(0,length(X.fdata$argvals)),}\if{latex}{\cr}\code{argvals=X.fdata$argvals,rangeval=X.fdata$rangeval)}.\if{latex}{\cr}
#'  If \code{beta0.fdata=NULL} (default), the function will test for the composite null hypothesis.
#' 
#' @param B Number of bootstrap replicates to calibrate the distribution of the test statistic.
#'  \code{B=5000} replicates are the recommended for carry out the test, although for exploratory analysis
#'   (\bold{not inferential}), an acceptable less time-consuming option is \code{B=500}.
#' @param est.method Estimation method for the unknown parameter \eqn{\beta}{\beta}, 
#' only used in the composite case. Mainly, there are two options: specify the number of basis 
#' elements for the estimated \eqn{\beta}{\beta} by \code{p} or optimally select \code{p} by a
#' data-driven criteria (see Details section for discussion). Then, it must be one of the following 
#' methods:
#' \itemize{
#'  \item \code{"pc"}: If \code{p}, the number of basis elements, is given, then \eqn{\beta}{\beta} is estimated by \code{\link{fregre.pc}}. Otherwise, an optimum \code{p} is chosen using \code{\link{fregre.pc.cv}} and the \code{"SICc"} criteria.
#'  \item \code{"pls"}: If \code{p} is given, \eqn{\beta}{\beta} is estimated by \code{\link{fregre.pls}}. Otherwise, an optimum \code{p} is chosen using \code{\link{fregre.pls.cv}} and the \code{"SICc"} criteria. 
#'  This is the default argument as it has been checked empirically that provides a good balance between the performance of the test and the estimation of \eqn{\beta}{\beta}.
#'  \item \code{"basis"}: If \code{p} is given, \eqn{\beta}{\beta} is estimated by \code{\link{fregre.basis}}. Otherwise, an optimum \code{p} is chosen using \code{\link{fregre.basis.cv}} and the \code{"GCV.S"} criteria. In these functions, the same basis for the arguments \code{basis.x} and \code{basis.b} is considered.
#'  The type of basis used will be the given by the argument \code{type.basis} and must be one of the class of \code{create.basis}. Further arguments passed to \code{\link{create.basis}} (not \code{rangeval} that is taken as the \code{rangeval} of \code{X.fdata}), can be passed throughout \code{\dots}.
#' }
#' @param p Number of elements of the basis considered. If it is not given, an optimal \code{p} will be chosen using a specific criteria (see \code{est.method} and \code{type.basis} arguments). 
#' @param type.basis Type of basis used to represent the functional process. Depending on the hypothesis, it will have a different interpretation:
#' \itemize{
#'  \item Simple hypothesis. One of these options:
#'   \itemize{
#'   \item \code{"bspline"}: If \code{p} is given, the functional process is expressed in a basis of \code{p} B-splines. If not, an optimal \code{p} will be chosen by \code{\link{optim.basis}}, using the \code{"GCV.S"} criteria.
#'   \item \code{"fourier"}: If \code{p} is given, the functional process is expressed in a basis of \code{p} Fourier functions. If not, an optimal \code{p} will be chosen by \code{\link{optim.basis}}, using the \code{"GCV.S"} criteria.
#'   \item \code{"pc"}: \code{p} must be given. Expresses the functional process in a basis of \code{p} principal components.
#'   \item \code{"pls"}: \code{p} must be given. Expresses the functional process in a basis of \code{p} partial least squares.
#'   }
#'  Although other basis types supported by \code{\link{create.basis}} are possible, \code{"bspline"} and \code{"fourier"} are recommended. Other basis types may cause incompatibilities.
#'  \item Composite hypothesis. This argument is only used when \if{latex}{\cr}\code{est.method="basis"} and, in this case, it specifies the type of basis used in the basis estimation method of the functional parameter. Again, basis
#'  \code{"bspline"} and \code{"fourier"} are recommended, as other basis types may cause incompatibilities.
#'  }
#'  
#' @param verbose Either to show or not information about computing progress.
#' @param plot.it Either to show or not a graph of the observed trajectory, 
#'  and the bootstrap trajectories under the null composite hypothesis, of the 
#'  process \eqn{R_n(\cdot)}{R_n(.)} (see Details). Note that if \code{plot.it=TRUE}, 
#'  the function takes more time to run. 
#' @param B.plot Number of bootstrap trajectories to show in the resulting plot of the test.
#'  As the trajectories shown are the first \code{B.plot} of \code{B}, \code{B.plot} must be 
#'  lower or equal to \code{B}.
#' @param G Number of projections used to compute the trajectories of the process
#'  \eqn{R_n(\cdot)}{R_n(.)} by Monte Carlo.
#' @param \dots Further arguments passed to \code{\link{create.basis}}.
#'  
#' @details The Functional Linear Model with scalar response (FLM), is defined as 
#' \eqn{Y=\big<X,\beta\big>+\epsilon}{Y=<X,\beta>+\epsilon}, for a functional process 
#' \eqn{X}{X} such that \eqn{E[X(t)]=0}{E[X(t)]=0}, \eqn{E[X(t)\epsilon]=0}{E[X(t)\epsilon]=0}
#'  for all \eqn{t}{t} and for a scalar variable \eqn{Y}{Y} such that \eqn{E[Y]=0}{E[Y]=0}.
#'  Then, the test assumes that \code{Y} and \code{X.fdata} are \bold{centred} and will automatically 
#'  center them. So, bear in mind that when you apply the test for \code{Y} and \code{X.fdata}, 
#'  actually,  you are applying it to \code{Y-mean(Y)} and \code{fdata.cen(X.fdata)$Xcen}.
#'  The test statistic corresponds to the Cramer-von Mises norm of the \emph{Residual Marked 
#'  empirical Process based on Projections} \eqn{R_n(u,\gamma)}{R_n(u,\gamma)} defined in 
#'  Garcia-Portugues \emph{et al.} (2014). 
#'  The expression of this process in a \eqn{p}{p}-truncated basis of the space \eqn{L^2[0,T]}{L^2[0,T]}
#'  leads to the \eqn{p}{p}-multivariate process \eqn{R_{n,p}\big(u,\gamma^{(p)}\big)}{R_{n,p}(u,\gamma^{(p)})}, 
#'  whose Cramer-von Mises norm is computed.
#'  The choice of an appropriate \eqn{p}{p} to represent the functional process \eqn{X}{X}, 
#'  in case that is not provided, is done via the estimation of \eqn{\beta}{\beta} for the composite 
#'  hypothesis. For the simple hypothesis, as no estimation of \eqn{\beta}{\beta} is done, the choice 
#'  of \eqn{p}{p} depends only on the functional process \eqn{X}{X}. As the result of the test may 
#'  change for different \eqn{p}{p}'s, we recommend to use an automatic criterion to select \eqn{p}{p} 
#'  instead of provide a fixed one.
#'  The distribution of the test statistic is approximated by a wild bootstrap resampling on the 
#'  residuals, using the \emph{golden section bootstrap}.
#'  Finally, the graph shown if \code{plot.it=TRUE} represents the observed trajectory, and the 
#'  bootstrap trajectories under the null, of the process RMPP \emph{integrated on the projections}:
#'  \deqn{R_n(u)\approx\frac{1}{G}\sum_{g=1}^G R_n(u,\gamma_g),}{R_n(u) \approx \frac{1}{G} \sum_{g=1}^G R_n(u,\gamma_g),}
#'   where \eqn{\gamma_g}{\gamma_g} are simulated as Gaussians processes. This gives a graphical idea of
#'    how \emph{distant} is the observed trajectory from the null hypothesis.
#'  
#' @return An object with class \code{"htest"} whose underlying structure is a list containing 
#' the following components:
#' \itemize{
#'  \item \code{statistic}: The value of the test statistic.
#'  \item \code{boot.statistics}: A vector of length \code{B} with the values of the bootstrap test statistics.
#'  \item \code{p.value}: The p-value of the test.
#'  \item \code{method}: The method used.
#'  \item \code{B}: The number of bootstrap replicates used.
#'  \item \code{type.basis}: The type of basis used.
#'  \item \code{beta.est}: The estimated functional parameter \eqn{\beta}{\beta} in the composite 
#'  hypothesis. For the simple hypothesis, the given \code{beta0.fdata}.
#'  \item \code{p}: The number of basis elements passed or automatically chosen.
#'  \item \code{ord}: The optimal order for PC and PLS given by \code{\link{fregre.pc.cv}} and \if{latex}{\cr} \code{\link{fregre.pls.cv}}. For other methods, it is set to \code{1:p}.
#'  \item \code{data.name}: The character string "Y=<X,b>+e".
#' }
#'  
#' @references 
#' Escanciano, J. C. (2006). A consistent diagnostic test for regression models using projections. Econometric Theory, 22, 1030-1051. \doi{10.1017/S0266466606060506}
#' 
#' Garcia-Portugues, E., Gonzalez-Manteiga, W. and Febrero-Bande, M. (2014). A goodness--of--fit test for the functional linear model with scalar response. Journal of Computational and Graphical Statistics, 23(3), 761-778. \doi{10.1080/10618600.2013.812519}
#' 
#' @note No NA's are allowed neither in the functional covariate nor in the scalar response.
#' 
#' @author  Eduardo Garcia-Portugues. Please, report bugs and suggestions to
#'  \if{latex}{\cr}\email{edgarcia@@est-econ.uc3m.es}
#'  
#' @seealso \code{\link{Adot}}, \code{\link{PCvM.statistic}}, \code{\link{rwild}}, 
#'  \code{\link{flm.Ftest}}, \code{\link{dfv.test}},
#'  \code{\link{fregre.pc}}, \code{\link{fregre.pls}},\if{latex}{\cr} \code{\link{fregre.basis}}, 
#'  \code{\link{fregre.pc.cv}}, \code{\link{fregre.pls.cv}},
#'  \code{\link{fregre.basis.cv}}, \code{\link{optim.basis}}, \if{latex}{\cr}
#'  \link[fda]{create.basis}
#' @examples
#' # Simulated example #
#' X=rproc2fdata(n=100,t=seq(0,1,l=101),sigma="OU")
#' beta0=fdata(mdata=cos(2*pi*seq(0,1,l=101))-(seq(0,1,l=101)-0.5)^2+
#'             rnorm(101,sd=0.05),argvals=seq(0,1,l=101),rangeval=c(0,1))
#' Y=inprod.fdata(X,beta0)+rnorm(100,sd=0.1)
#' 
#' dev.new(width=21,height=7)
#' par(mfrow=c(1,3))
#' plot(X,main="X")
#' plot(beta0,main="beta0")
#' plot(density(Y),main="Density of Y",xlab="Y",ylab="Density")
#' rug(Y)
#' 
#' \dontrun{
#' # Composite hypothesis: do not reject FLM
#' pcvm.sim=flm.test(X,Y,B=50,B.plot=50,G=100,plot.it=TRUE)
#' pcvm.sim
#' flm.test(X,Y,B=5000)
#'  
#' # Estimated beta
#' dev.new()
#' plot(pcvm.sim$beta.est)
#' 
#' # Simple hypothesis: do not reject beta=beta0
#' flm.test(X,Y,beta0.fdata=beta0,B=50,B.plot=50,G=100)
#' flm.test(X,Y,beta0.fdata=beta0,B=5000) 
#' 
#' # AEMET dataset #
#' data(aemet)
#' # Remove the 5\% of the curves with less depth (i.e. 4 curves)
#' dev.new()
#' res.FM=depth.FM(aemet$temp,draw=TRUE)
#' qu=quantile(res.FM$dep,prob=0.05)
#' l=which(res.FM$dep<=qu)
#' lines(aemet$temp[l],col=3)
#' aemet$df$name[l]
#' 
#' # Data without outliers 
#' wind.speed=apply(aemet$wind.speed$data,1,mean)[-l]
#' temp=aemet$temp[-l]
#' # Exploratory analysis: accept the FLM
#' pcvm.aemet=flm.test(temp,wind.speed,est.method="pls",B=100,B.plot=50,G=100)
#' pcvm.aemet
#' 
#' # Estimated beta
#' dev.new()
#' plot(pcvm.aemet$beta.est,lwd=2,col=2)
#' # B=5000 for more precision on calibration of the test: also accept the FLM
#' flm.test(temp,wind.speed,est.method="pls",B=5000) 
#' 
#' # Simple hypothesis: rejection of beta0=0? Limiting p-value...
#' dat=rep(0,length(temp$argvals))
#' flm.test(temp,wind.speed, beta0.fdata=fdata(mdata=dat,argvals=temp$argvals,
#'                                             rangeval=temp$rangeval),B=100)
#' flm.test(temp,wind.speed, beta0.fdata=fdata(mdata=dat,argvals=temp$argvals,
#'                                             rangeval=temp$rangeval),B=5000) 
#'                                             
#' # Tecator dataset #
#' data(tecator)
#' names(tecator)
#' absorp=tecator$absorp.fdata
#' ind=1:129 # or ind=1:215
#' x=absorp[ind,]
#' y=tecator$y$Fat[ind]
#' tt=absorp[["argvals"]]
#' 
#' # Exploratory analysis for composite hypothesis with automatic choose of p
#' pcvm.tecat=flm.test(x,y,B=100,B.plot=50,G=100)
#' pcvm.tecat
#' 
#' # B=5000 for more precision on calibration of the test: also reject the FLM
#' flm.test(x,y,B=5000) 
#' 
#' # Distribution of the PCvM statistic
#' plot(density(pcvm.tecat$boot.statistics),lwd=2,xlim=c(0,10),
#'               main="PCvM distribution", xlab="PCvM*",ylab="Density")
#' rug(pcvm.tecat$boot.statistics)
#' abline(v=pcvm.tecat$statistic,col=2,lwd=2)
#' legend("top",legend=c("PCvM observed"),lwd=2,col=2)
#' 
#' # Simple hypothesis: fixed p
#' dat=rep(0,length(x$argvals))
#' flm.test(x,y,beta0.fdata=fdata(mdata=dat,argvals=x$argvals,
#'                                rangeval=x$rangeval),B=100,p=11)
#'                                
#' # Simple hypothesis, automatic choose of p
#' flm.test(x,y,beta0.fdata=fdata(mdata=dat,argvals=x$argvals,
#'                                rangeval=x$rangeval),B=100)
#' flm.test(x,y,beta0.fdata=fdata(mdata=dat,argvals=x$argvals,
#'                                rangeval=x$rangeval),B=5000)
#' }
#' @keywords htest models regression

# PCvM test for the composite hypothesis with bootstrap calibration
#' @rdname flm.test
#' @export
flm.test=function(X.fdata,Y,beta0.fdata=NULL,B=5000,est.method="pls",
                  p=NULL,type.basis="bspline",verbose=TRUE,plot.it=TRUE,
                  B.plot=100,G=200,...){
	
	# Check B.plot
	if(plot.it & B.plot>B) stop("B.plot must be less or equal than B")
	
	# Number of functions
	n=dim(X.fdata)[1]
	
	if(verbose) cat("Computing estimation of beta... ")	
	
	## COMPOSITE HYPOTHESIS ##
	if(is.null(beta0.fdata)){

		# Center the data first
		X.fdata=fdata.cen(X.fdata)$Xcen
		Y=Y-mean(Y)
		
		## 1. Optimal estimation of beta and the basis order ##
		
		if(est.method=="pc"){
			
			if(is.null(p)){

				# Method
				meth="PCvM test for the functional linear model using optimal PC basis representation"

				# Choose the number of basis elements: SICc is probably the best criteria
				mod.pc=fregre.pc.cv(fdataobj=X.fdata,y=Y,kmax=1:10,criteria="SICc")
				p.opt=length(mod.pc$pc.opt)
				ord.opt=mod.pc$pc.opt
				
				# PC components to be passed to the bootstrap
				pc.comp=mod.pc$fregre.pc$fdata.comp # pc.comp=mod.pc$fregre.pc$pc
				pc.comp$l=mod.pc$pc.opt
				
				# Express X.fdata and beta.est in the PC basis
				basis.pc=mod.pc$fregre.pc$fdata.comp$basis
				if(length(pc.comp$l)!=1){
					X.est=fdata(mdata=mod.pc$fregre.pc$fdata.comp$coefs[,mod.pc$fregre.pc$l]%*%mod.pc$fregre.pc$fdata.comp$basis$data[mod.pc$fregre.pc$l,],argvals=X.fdata$argvals,rangeval=X.fdata$rangeval) # X.est=fdata(mdata=mod.pc$fregre.pc$pc$coefs[,mod.pc$fregre.pc$l]%*%mod.pc$fregre.pc$pc$basis$data[mod.pc$fregre.pc$l,],argvals=X.fdata$argvals,rangeval=X.fdata$rangeval)
				}else{
					X.est=fdata(mdata=mod.pc$fregre.pc$fdata.comp$coefs[,mod.pc$fregre.pc$l]%*%t(mod.pc$fregre.pc$fdata.comp$basis$data[mod.pc$fregre.pc$l,]),argvals=X.fdata$argvals,rangeval=X.fdata$rangeval) # X.est=fdata(mdata=mod.pc$fregre.pc$pc$coefs[,mod.pc$fregre.pc$l]%*%t(mod.pc$fregre.pc$pc$basis$data[mod.pc$fregre.pc$l,]),argvals=X.fdata$argvals,rangeval=X.fdata$rangeval)
				}
				beta.est=mod.pc$fregre.pc$beta.est
				norm.beta.est=norm.fdata(beta.est)
				
				# Compute the residuals
				e=mod.pc$fregre.pc$residuals
				
			}else{
			
				# Method
				meth=paste("PCvM test for the functional linear model using a representation in a PC basis of ",p,"elements") 
				
				# Estimation of beta on the given fixed basis
				mod.pc=fregre.pc(fdataobj=X.fdata,y=Y,l=1:p)
				p.opt=p
				ord.opt=mod.pc$l

				# PC components to be passed to the bootstrap
				pc.comp=mod.pc$pc
				pc.comp$l=mod.pc$l
				
				# Express X.fdata and beta.est in the basis
				if(p!=1){
					X.est=fdata(mdata=mod.pc$fdata.comp$coefs[,mod.pc$l]%*%mod.pc$fdata.comp$basis$data[mod.pc$l,],argvals=X.fdata$argvals,rangeval=X.fdata$rangeval)
				}else{
					X.est=fdata(mdata=mod.pc$fdata.comp$coefs[,mod.pc$l]%*%t(mod.pc$fdata.comp$basis$data[mod.pc$l,]),argvals=X.fdata$argvals,rangeval=X.fdata$rangeval)
				}
				beta.est=mod.pc$beta.est
				norm.beta.est=norm.fdata(beta.est)
				
				# Compute the residuals
				e=mod.pc$residuals
			
			}
			
		}else if(est.method=="pls"){
			
			if(is.null(p)){
        # print("PLS1")
				# Method
				meth="PCvM test for the functional linear model using optimal PLS basis representation"
			
				# Choose the number of the basis: SICc is probably the best criteria
				mod.pls=fregre.pls.cv(fdataobj=X.fdata,y=Y,kmax=10,criteria="SICc") 
        #   print("PLS2")				
				p.opt=length(mod.pls$pls.opt)
				ord.opt=mod.pls$pls.opt
				
				# PLS components to be passed to the bootstrap
				pls.comp=mod.pls$fregre.pls$fdata.comp
				pls.comp$l=mod.pls$pls.opt
						
				# Express X.fdata and beta.est in the PLS basis
				basis.pls=mod.pls$fregre.pls$fdata.comp$basis
				if(length(pls.comp$l)!=1){
					X.est=fdata(mdata=mod.pls$fregre.pls$fdata.comp$coefs[,mod.pls$fregre.pls$l]%*%mod.pls$fregre.pls$fdata.comp$basis$data[mod.pls$fregre.pls$l,],argvals=X.fdata$argvals,rangeval=X.fdata$rangeval)
				}else{
					X.est=fdata(mdata=mod.pls$fregre.pls$fdata.comp$coefs[,mod.pls$fregre.pls$l]%*%t(mod.pls$fregre.pls$fdata.comp$basis$data[mod.pls$fregre.pls$l,]),argvals=X.fdata$argvals,rangeval=X.fdata$rangeval)
				}
				beta.est=mod.pls$fregre.pls$beta.est
				norm.beta.est=norm.fdata(beta.est)
				
				# Compute the residuals
				e=mod.pls$fregre.pls$residuals
				
			}else{
			
				# Method
				meth=paste("PCvM test for the functional linear model using a representation in a PLS basis of ",p,"elements") 
			
				# Estimation of beta on the given fixed basis
				mod.pls=fregre.pc(fdataobj=X.fdata,y=Y,l=1:p)
				p.opt=p
				ord.opt=mod.pls$l

				# PLS components to be passed to the bootstrap
				pls.comp=mod.pls$fdata.comp
				pls.comp$l=mod.pls$l
				
				# Express X.fdata and beta.est in the basis
				if(p!=1){
					X.est=fdata(mdata=mod.pls$fdata.comp$coefs[,mod.pls$l]%*%mod.pls$fdata.comp$basis$data[mod.pls$l,],argvals=X.fdata$argvals,rangeval=X.fdata$rangeval)
				}else{
					X.est=fdata(mdata=mod.pls$fdata.comp$coefs[,mod.pls$l]%*%t(mod.pls$fdata.comp$basis$data[mod.pls$l,]),argvals=X.fdata$argvals,rangeval=X.fdata$rangeval)
				}
				beta.est=mod.pls$beta.est
				norm.beta.est=norm.fdata(beta.est)
				
				# Compute the residuals
				e=mod.pls$residuals
			
			}
			
		}else if(est.method=="basis"){
			
			if(is.null(p)){
			
				# Method
				meth=paste("PCvM test for the functional linear model using optimal",type.basis,"basis representation")
			
				# Choose the number of the bspline basis with GCV.S
				mod.basis=fregre.basis.cv(fdataobj=X.fdata,y=Y,basis.x=seq(5,30,by=1),basis.b=NULL,type.basis=type.basis,type.CV=GCV.S,verbose=FALSE,...)
				p.opt=mod.basis$basis.x.opt$nbasis
				ord.opt=1:p.opt
				
				# Express X.fdata and beta.est in the optimal basis
				basis.opt=mod.basis$basis.x.opt
				X.est=mod.basis$x.fd
				beta.est=mod.basis$beta.est
				norm.beta.est=norm.fd(beta.est)
				
				# Compute the residuals
				e=mod.basis$residuals
			
			}else{
				
				# Method
				meth=paste("PCvM test for the functional linear model using a representation in a",type.basis,"basis of ",p,"elements") 

				# Estimation of beta on the given fixed basis
				basis.opt=do.call(what=paste("create.",type.basis,".basis",sep=""),args=list(rangeval=X.fdata$rangeval,nbasis=p,...))
				mod.basis=fregre.basis(fdataobj=X.fdata,y=Y,basis.x=basis.opt,basis.b=basis.opt)
				p.opt=p
				ord.opt=1:p.opt

				# Express X.fdata and beta.est in the basis
				X.est=mod.basis$x.fd
				beta.est=mod.basis$beta.est
				norm.beta.est=norm.fd(beta.est)
				
				# Compute the residuals
				e=mod.basis$residuals
			
			}
			
		}else{
		
			stop(paste("Estimation method",est.method,"not implemented."))
			
		}
	
	## SIMPLE HYPOTHESIS ##
	}else{
	
		## 1. Optimal representation of X and beta0 ##
		
		# Choose the number of basis elements
		if(type.basis!="pc" & type.basis!="pls"){
			
			# Basis method: select the number of elements of the basis is it is not given
			if(is.null(p)){
				
				# Method
				meth=paste("PCvM test for the simple hypothesis in a functional linear model, using a representation in an optimal ",type.basis,"basis") 
				
				p.opt=optim.basis(X.fdata,type.basis=type.basis,numbasis=seq(31,71,by=2),verbose=FALSE)$numbasis.opt
				if(p.opt==31){
					cat("Readapting interval...\n")
					p.opt=optim.basis(X.fdata,type.basis=type.basis,numbasis=seq(5,31,by=2),verbose=FALSE)$numbasis.opt
				}else if(p.opt==71){
					cat("Readapting interval...\n")
					p.opt=optim.basis(X.fdata,type.basis=type.basis,numbasis=seq(71,101,by=2),verbose=FALSE)$numbasis.opt
				}
				
				ord.opt=1:p.opt
				
			}else{
				
				# Method
				meth=paste("PCvM test for the simple hypothesis in a functional linear model, using a representation in a ",type.basis," basis of ",p,"elements") 
				
				p.opt=p
				ord.opt=1:p.opt
			}
			
			# Express X.fdata in the basis
			X.est=fdata2fd(X.fdata,type.basis=type.basis,nbasis=p)
			beta.est=fdata2fd(beta0.fdata,type.basis=type.basis,nbasis=p)
			
			# Compute the residuals
			e=Y-inprod(X.est,beta.est)

		}else if (type.basis=="pc"){
			
			if(is.null(p)) stop("Simple hypothesis with type.basis=\"pc\" need the number of components p") 
			# Method
			meth=paste("PCvM test for the simple hypothesis in a functional linear model, using a representation in a PC basis of ",p,"elements") 
			
			# Express X.fdata in a PC basis
			fd2pc=fdata2pc(X.fdata,ncomp=p)
			if(length(fd2pc$l)!=1){
				X.est=fdata(mdata=fd2pc$coefs[,fd2pc$l]%*%fd2pc$basis$data[fd2pc$l,],argvals=X.fdata$argvals,rangeval=X.fdata$rangeval)
			}else{
				X.est=fdata(mdata=fd2pc$coefs[,fd2pc$l]%*%t(fd2pc$basis$data[fd2pc$l,]),argvals=X.fdata$argvals,rangeval=X.fdata$rangeval)
			}
			beta.est=beta0.fdata
			p.opt=p
			ord.opt=1:p.opt
			
			# Compute the residuals
			e=Y-inprod.fdata(X.est,beta.est)
			
		}else if(type.basis=="pls"){
		
			if(is.null(p)) stop("Simple hypothesis with type.basis=\"pls\" need the number of components p")

			# Method
			meth=paste("PCvM test for the simple hypothesis in a functional linear model, using a representation in a PLS basis of ",p,"elements") 
			
			# Express X.fdata in a PLS basis
			fd2pls=fdata2pls(X.fdata,Y,ncomp=p)
			if(length(fd2pls$l)!=1){
				X.est=fdata(mdata=fd2pls$coefs[,fd2pls$l]%*%fd2pls$basis$data[fd2pls$l,],argvals=X.fdata$argvals,rangeval=X.fdata$rangeval)
			}else{
				X.est=fdata(mdata=fd2pls$coefs[,fd2pls$l]%*%t(fd2pls$basis$data[fd2pls$l,]),argvals=X.fdata$argvals,rangeval=X.fdata$rangeval)
			}		
			beta.est=beta0.fdata
			p.opt=p
			ord.opt=1:p.opt
			
			# Compute the residuals
			e=Y-inprod.fdata(X.est,beta.est)
			
		}else{
		
			stop(paste("Type of basis method",type.basis,"not implemented."))

		}
	
	}
	
	## 2. Bootstrap calibration ##
	
	# Start up
	pcvm.star=numeric(B)
	e.hat.star=matrix(ncol=n,nrow=B)
	
	# Calculus of the Adot.vec
	Adot.vec=Adot(X.est)

	# REAL WORLD
	pcvm=PCvM.statistic(X=X.est,residuals=e,p=p.opt,Adot.vec=Adot.vec)
	
	# BOOTSTRAP WORLD 
	if(verbose) cat("Done.\nBootstrap calibration...\n ")
	if(verbose) pb=txtProgressBar(style=3)

	## COMPOSITE HYPOTHESIS ##
	if(is.null(beta0.fdata)){
					
		# Calculate the design matrix of the linear model
		# This allows to resample efficiently the residuals without estimating again the beta
		if(est.method=="pc"){
			
			if(is.null(p)){
				# Design matrix for the PC estimation
				X.matrix=mod.pc$fregre.pc$lm$x
			}else{
				# Design matrix for the PC estimation
				X.matrix=mod.pc$lm$x
			}
			
		}else if(est.method=="pls"){

			if(is.null(p)){	
				# Design matrix for the PLS estimation
				X.matrix=mod.pls$fregre.pls$lm$x
			}else{
				# Design matrix for the PLS estimation
				X.matrix=mod.pls$lm$x
			}		
		}else if(est.method=="basis"){

			if(is.null(p)){	
				# Design matrix for the basis estimation
				X.matrix=mod.basis$lm$x
			}else{
				# Design matrix for the basis estimation
				X.matrix=mod.basis$lm$x
			}		
		}
	  
	  # Projection matrix
	  P=(diag(rep(1,n))-X.matrix%*%solve(t(X.matrix)%*%X.matrix)%*%t(X.matrix))
	  
		# Bootstrap resampling
		for(i in 1:B){
		
			# Generate bootstrap errors
			e.hat=rwild(e,"golden")
	
			# Calculate Y.star
			Y.star=Y-e+e.hat
			
			# Residuals from the bootstrap estimated model
			e.hat.star[i,]=P%*%Y.star
		
			# Calculate PCVM.star
			pcvm.star[i]=PCvM.statistic(X=X.est,residuals=e.hat.star[i,],p=p.opt,Adot.vec=Adot.vec)
			
			if(verbose) setTxtProgressBar(pb,i/B)
		
		}
					
	## SIMPLE HYPOTHEIS ##
	}else{
		
		# Bootstrap resampling
		for(i in 1:B){
				
			# Generate bootstrap errors
			e.hat.star[i,]=rwild(e,"golden")
		
			# Calculate PCVM.star
			pcvm.star[i]=PCvM.statistic(X=X.est,residuals=e.hat.star[i,],p=p.opt,Adot.vec=Adot.vec)
				
			if(verbose) setTxtProgressBar(pb,i/B)
		
		}
					
	}

	## 3. MC estimation of the p-value and order the result ##
	
	# Compute the p-value
	pvalue=sum(pcvm.star>pcvm)/B
	
	## 4. Graphical representation of the integrated process ##
	
	if(verbose) cat("\nDone.\nComputing graphical representation... ")
	if(is.null(beta0.fdata) & plot.it){
	
		gamma=rproc2fdata(n=G,t=X.fdata$argvals,sigma="OU",par.list=list(theta=2/diff(range(X.fdata$argvals))))
		gamma=gamma/drop(norm.fdata(gamma))
		ind=drop(inprod.fdata(X.fdata,gamma))
		
		r=0.9*max(max(ind),-min(ind))
		u=seq(-r,r,l=200)
		mean.proc=numeric(length(u))
		mean.boot.proc=matrix(ncol=length(u),nrow=B.plot)
		
		res=apply(ind,2,function(iind){
			iind.sort=sort(iind,index.return=TRUE)
			stepfun(x=iind.sort$x,y=c(0,cumsum(e[iind.sort$ix])))(u)/sqrt(n)
		}
		)
		mean.proc=apply(res,1,mean)
	
		for(i in 1:B.plot){
		
			res=apply(ind,2,function(iind){
				iind.sort=sort(iind,index.return=TRUE)
				stepfun(x=iind.sort$x,y=c(0,cumsum(e.hat.star[i,iind.sort$ix])))(u)/sqrt(n)
			}
			)
			mean.boot.proc[i,]=apply(res,1,mean)
		
		}
		
		# Plot
		dev.new()
		plot(u,mean.proc,ylim=c(min(mean.proc,mean.boot.proc),max(mean.proc,mean.boot.proc))*1.05,type="l",xlab=expression(paste(symbol("\341"),list(X, gamma),symbol("\361"))),ylab=expression(R[n](u)),main="")
		for(i in 1:B.plot) lines(u,mean.boot.proc[i,],lty=2,col=gray(0.8))
		lines(u,mean.proc)
		text(x=0.75*u[1],y=0.75*min(mean.proc,mean.boot.proc),labels=sprintf("p-value=%.3f",pvalue))
		
	}
	if(verbose) cat("Done.\n")
	# Result: class htest
	names(pcvm)="PCvM statistic"
	result=structure(list(statistic=pcvm,boot.statistics=pcvm.star,p.value=pvalue,method=meth,B=B,type.basis=type.basis,beta.est=beta.est,p=p.opt,ord=ord.opt,data.name="Y=<X,b>+e"))
							
	class(result) <- "htest"
	return(result)

}

               