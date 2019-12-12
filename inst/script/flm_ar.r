
################################################################################
# Estimation functional beta parameter  using simulated data
# where the theoretical  beta is known
#
# For more details see: Febrero-Bande, M., Galeano, P. and Gonzalez-Manteiga, W.
# (2010). Measures of influence for the functional linear model with scalar
# response. Journal of Multivariate Analysis 101, 327-339.
################################################################################


#####################################################################################
# The function fdata.brown() generates functional data independent of a movement
# Brownian defined in the interval [0, Tend]. The data appear in a matrix
# (J x n) size where n is the number of curves and J is the number of discretized
# points, which is supposed to be the same for all curves.
#####################################################################################
fdata.brown <- function(n=1000,J=101,tt=seq(0,1,len=J),rtt=range(tt),names=NULL){
  X <- matrix(NA,ncol=J,nrow=n)
  for (i in 1:n){
    X[i,1] <- 0
    for (j in 2:J){X[i,j] <- X[i,j-1] + rnorm(1,0,sqrt(tt[j]-tt[j-1]))}
  }
  X<-fdata(X,tt,rtt,names=names)
  return(X)
}
#####################################################################################
# The function fdata.cfs.2003.a()   generates simulated data of the functional
# linear linear model with scalar response, Cardot, Ferraty y Sarda (2003).
#####################################################################################
fdata.cfs.2003.a <- function(x=NULL,snr,rho=0,p=1,s2eps,sd){
  if  (is.null(x)) x<- fdata.brown()
  tj <- x$argvals
  rtt<-x$rangeval
  n<-nrow(x)
  X <- x$data
  J<-ncol(x)
  n<-nrow(x)
  v1 <- sqrt(2) * sin(.5 * pi * tj)
  v2 <- sqrt(2) * sin(1.5 * pi * tj)
  v3 <- sqrt(2) * sin(2.5 * pi * tj)
  lamb1 <- 1 / (.5 * pi)^2
  lamb2 <- 1 / (1.5 * pi)^2
  lamb3 <- 1 / (2.5 * pi)^2
  b1 <- 2 / sqrt(2)
  b2 <- 4 / sqrt(2)
  b3 <- 5 / sqrt(2)
  bet <- b1 * v1 + b2 * v2 + b3 * v3
  varXg <- b1^2 * lamb1 + b2^2 * lamb2 + b3^2 * lamb3
  Xmean <- fdata.cen(x)
  bet<-fdata(bet,tj,rtt)                        	
  yp<-p*drop(inprod.fdata(x,bet))+(1-p)*as.numeric(norm.fdata(x))
  s2eps<-snr*var(yp)#/(1-snr)
  if ( any(rho==0)) {
    if (missing(sd))     sd<-sqrt(s2eps)
    e <- rnorm(n,0,sd=sd)
    ar.model=NULL
  }
  else{
    if (length(rho)==1) 
      sd<-sqrt(s2eps*(1-rho^2))     
    if (length(rho)==2) {
      f1<-rho[1]
      f2<-rho[2]
      aux<-(1-f1*(f1/(1-f2))-f2*(f1^2/(1-f2)+f2))
      sd<-sqrt(s2eps*aux)
    }
    ar.model<-arima.sim(n =n,list(order = c(length(rho),0,0), ar =rho),
                        n.start=10000,mean=0,sd=sd)
    e<-as.vector(ar.model)
  }                      
  #print(sd)
  yp <- yp +e
  return(list("x"=x,"y"=yp,"v1"=fdata(v1,tj,rtt),"v2"=fdata(v2,tj,rtt),
              "v3"=fdata(v3,tj,rtt),"bet"=bet,"s2eps"=s2eps,"snr"=snr,"rho"=rho,"e"=e,"sd"=sd,"ar.model"=ar.model))
}

#####################################################################################
# The function data.cfs.2003.b()  generates simulated data of the functional
# linear linear model with scalar response, Cardot, Ferraty y Sarda (2003).
# One has to notice that case (a) favors the SPCR estimator since ?1 is a linear
# combination of the first three eigenfunctions of ?. Case (b) is more general since
# ?2 combines log and periodic components
#####################################################################################
fdata.cfs.2003.b <- function(x=NULL,snr,rho,p=1,s2eps,sd){
  if  (is.null(x)) x<- fdata.brown()
  tj <- x$argvals
  rtt<-x$rangeval
  X <- x$data
  J<-ncol(x)
  n<-nrow(x)
  maxK <- qr(X)$rank
  varXg <- rep(NA,len=maxK)
  bet <- log(15*tj^2 + 10) + cos(4*pi*tj)
  for (k in 1:maxK){
    integrand <- function(t) {sqrt(2)*(log(15*t^2+10)+cos(4*pi*t))*sin((k-.5)*pi*t)}
    bk <- integrate(integrand, lower = rtt[1], upper = rtt[2])$value
    varXg[k] <- bk^2 / ((k - .5) * pi)^2
  }                  	
  Xmean <- fdata.cen(x)
  bet<-fdata(bet,tj,rtt)
  yp<-p*drop(inprod.fdata(x,bet))+(1-p)*as.numeric(norm.fdata(x))
  s2eps<-snr*var(yp)
  if ( any(rho==0)) {
    if (missing(sd))     sd<-sqrt(s2eps)
    e <- rnorm(n,0,sd=sd)
    ar.model=NULL
  }
  else{
    
    if (length(rho)==1) 
      sd<-sqrt(s2eps*(1-rho^2))     
    if (length(rho)==2) {
      f1<-rho[1]
      f2<-rho[2]
      aux<-(1-f1*(f1/(1-f2))-f2*(f1^2/(1-f2)+f2))
      sd<-sqrt(s2eps*aux)
    }
    ar.model<-arima.sim(n =n,list(order = c(length(rho),0,0), ar =rho),
                        n.start=10000,mean=0,sd=sd)
    e<-as.vector(ar.model)
  }                      
  yp <- yp +e
  return(list("x"=x,"y"=yp,"v1"=fdata(v1,tj,rtt),"v2"=fdata(v2,tj,rtt),
              "v3"=fdata(v3,tj,rtt),"bet"=bet,"s2eps"=s2eps,"snr"=snr,
              "rho"=rho,"e"=e,"sd"=sd,"ar.model"=ar.model))
}
#####################################################################################
# The function data.hh.2006 generates simulated data of the functional linear
# model with scalar response, de Hall and Housseini-Nasab (2006).
#####################################################################################
fdata.hh.2006 <- function(n,J,R2,rho,nls=1){
  tj <- seq(0,1,length.out=J)
  maxK <- min(c(n,J-1))
  v <- matrix(NA,ncol=maxK,nrow=J)
  bet <- pi^2 * (tj^2 - 1/3)
  varXg <- rep(NA,len=maxK)
  for (k in 1 : maxK){
    v[,k] <- sqrt(2) * cos(k * pi * tj)
    bk <- (-1)^k * k^(-2) * 2^(3/2)
    varXg[k] <-  bk^2 / k^2
  }
  s2eps <- sum(varXg) * (1 / R2 - 1)
  X <- matrix(0,ncol=J,nrow=n)
  for (i in 1 : n){for (k in 1:maxK){X[i,] <- X[i,] + rnorm(1,0,k^(-1)) * t(v[,k])}}
  x<-fdata(X,tj)
  
  Xmean <- fdata.cen(x)
  yp <- (Xmean$Xcen$data %*% bet)^nls * (1 / J)
  if (missing(rho)) {
    e <- matrix(rnorm(n,0,sqrt(s2eps)),ncol=1)
    rho=0             }
  else{e<-simul.ar1(n, rho=rho,burn=1000,mu=0,sd=sqrt(s2eps*(1-rho^2)))} 
  y <- yp + e
  R2 <- varXg/(1 + varXg)
  return(list("x"=x,"y"=y,"bet"=fdata(bet,tj),"s2eps"=s2eps,"rho"=rho,"varX"=varXg,"e"=e))
}

#####################################################################################
# The function data.cfs.2003.c()generates simulated data of the functional linear
# model with scalar response, Cardot, Ferraty y Sarda (2003), but with the
# inclusion of different eigenfunctions.
#####################################################################################
fdata.cfs.2003.c <- function(x=NULL,snr,rho,nls=1){
  if  (is.null(x)) x<- fdata.brown()
  tj <- x$argvals
  rtt<-x$rangeval
  X <- x$data
  J<-ncol(x)
  n<-nrow(x)
  v3 <- sqrt(2) * sin(2.5 * pi * tj)
  v5 <- sqrt(2) * sin(4.5 * pi * tj)
  v7 <- sqrt(2) * sin(6.5 * pi * tj)
  lamb3 <- 1 / (2.5 * pi)^2
  lamb5 <- 1 / (4.5 * pi)^2
  lamb7 <- 1 / (6.5 * pi)^2
  b3 <- 2 / sqrt(2)
  b5 <- 4 / sqrt(2)
  b7 <- 5 / sqrt(2)
  bet <- b3 * v3 + b5 * v5 + b7 * v7
  varXg <- b3^2 * lamb3 + b5^2 * lamb5 + b7^2 * lamb7
  #s2eps <-  * (1 / R2 - 1)
  s2eps<-snr*varXg/(1-snr)  
  Xmean <- fdata.cen(X)
  yp <- (Xmean$Xcen$data %*% bet^nls) * (1 / J)
  if (missing(rho)) {
    e <- matrix(rnorm(n,0,sqrt(s2eps)),ncol=1)
    rho=0             }
  else{e<-simul.ar1(n, rho=rho,burn=1000,mu=0,sd=sqrt(s2eps*(1-rho^2)))} 
  y <- yp + e
  return(list("x"=x,"y"=y,"v3"=v3,"v5"=v5,"v7"=v5,"bet"=fdata(bet,tj,rtt),
              "s2eps"=s2eps,"rho"=rho))
}


################################################################################
# Case studied
#n=100;J=50;r2=0.3
#fx=fdata.brown(n,J)
#fdatos.a=fdata.cfs.2003.a(fx,R2=r2)
#fdatos.b=fdata.cfs.2003.b(fx,R2=r2)
#fdatos.c=fdata.cfs.2003.c(fx,R2=r2)
#fdatos.d<-fdata.hh.2006(n,J,R2=r2)

################################################################################
# Example for the data generated in Figure 1 Febrero et al. (2010)
#par(mfrow=c(2,2))
#plot(fdatos.a$bet,main="cfs.2003.a")
#plot(fdatos.b$bet,col=2,main="cfs.2003.b")
#plot(fdatos.c$bet,col=3,main="cfs.2003.c")
#plot(fdatos.d$bet,col=4,main="hh.2006")
#par(mfrow=c(2,2))
#plot(fdatos.a$x,main="cfs.2003.a")
#plot(fdatos.b$x,col=2,main="cfs.2003.b")
#plot(fdatos.c$x,col=3,main="cfs.2003.c")
#plot(fdatos.d$x,col=4,main="hh.2006")

################################################################################
# Example estimation functional beta parameter  using  simulated data  generates
# via fdata.cfs.2003.a()
# obj<-fdatos.a
# x<-obj$x
# y<-obj$y
# dev.off()
# plot(obj$bet,main="Beta theoretical",ylab="Beta(t)")

################################################################################
# Functional linear model with principal components
# res=fregre.pc(x,y,1:3)
# Displaying Theoretical and estimated  beta
# plot(c(obj$bet,res$beta.est),lwd=1:2,main="Beta(t) - Beta.est(t)",ylab="beta(t)")

##########
#####################################################################################
# The function simul.st  spatio-temporal errors, AR(1)+Exp(-s/rho)
#####################################################################################
simul.st<-function(n, phi=0.8,rho=.8,burn=1000,mu=0,sd=1) {
n<-100
J=101
tt=seq(0,1,len=J)
rho=.8;burn=1000;mu=0;sd=1
k<-0:(J-1)
t<-s<-tj <- tt
tjs<-expand.grid(t,s)         #sample(x,y)
#sigmas<-matrix(.5*exp(-.8*abs(tjs[,1]-tjs[,2])),ncol=J)
#sigmat<-matrix((rho^k)*100*abs(tjs[,1]-tjs[,2]),ncol=J)
#xs<-fdata(mvrnorm(n=n,rep(0,len=J), Sigma=sigmas),tj,rtt)  
#xt<-fdata(mvrnorm(n=n,rep(0,len=J), Sigma=sigmat),tj,rtt)
rho=.95
rho2=.5
 W0<-toeplitz(rho^k)*matrix(exp(-(tjs[,1]-tjs[,2])/rho2),ncol=J)
 #W0<-matrix(exp(-(tjs[,1]-tjs[,2])/rho2),ncol=J)
# W0<-toeplitz(rho^k)*toeplitz(rho2^k)
} 
#ss<-simul.st(100)
#round(ss[50,],2)
#image(simul.st(100))

#####################################################################################
# The function fdata.lm()   generates simulated data of the functional
# linear linear model with scalar response, Febrero-Bande, Gonzalez-Manteiga (2012).
#####################################################################################
 fdata.lm <- function(n=1000,J=101,tt=seq(0,1,len=J),rtt=range(tt),rho,snr,rate=5){
	t<-s<-tj <- tt
  tjs<-expand.grid(t,s)
  sigma1<-matrix(.5*exp(-.8*abs(tjs[1]-tjs[2]))[,1],ncol=J)
  sigma2<-matrix(.4*exp(-.6*abs(tjs[1]-tjs[2]))[,1],ncol=J)
  x1<-fdata(mvrnorm(n=n,sin(2*pi * tj), Sigma=sigma1),tj,rtt)
  x2<-fdata(mvrnorm(n=n,rep(0,len=J), Sigma=sigma1),tj,rtt)  #mov browniano!
#  plot(x1);  plot(x2)
	v1 <- sqrt(2) * sin(pi * tj)
	b1 <- 2 / sqrt(2)
	bet<-bet1 <- fdata(b1 * v1 *(t*(1-t)),tj,rtt)
	bet2 <- fdata(rep(1,len=J),tj,rtt)	
	S1<-inprod.fdata(x1,bet1)
  S2<-inprod.fdata(x2,bet2)
  varS<-var(drop(S1+S2))*snr
  rate<-1
	yp <-   rate*S1	- S2 
  yp <-   rate*S1#	- S2 
  s2eps<-snr*var(yp)/(1-snr)   
  if (missing(rho)) {
   	e <- rnorm(n,sd=sqrt(varS))
    rho=0             }
  else{
  #e<-simul.ar1(n, rho=rho,burn=1000,mu=0,sd=sqrt(varS))
      e<-as.vector(arima.sim(n =n,list(order = c(length(rho),0,0), ar =rho),n.start=1000,mean=0,sd=sqrt(s2eps)))
      }
  yp <- yp +e
	return(list("x"=x1,bet=bet,"x1"=x1,"x2"=x2,"y"=yp,"bet1"=bet1,"bet2"=bet2,"snr"=snr,rho=rho))
}
#####################################################################################
# The function fdata.nlm()   generates simulated data of the functional
# non linear model with scalar response,  modifid version (only X1) Febrero-Bande, Gonzalez-Manteiga (2012).
#####################################################################################
 fdata.nlm <- function(n=1000,J=101,tt=seq(0,1,len=J),rtt=range(tt),rho,snr,nls=1){
	t<-s<-tj <- tt
  tjs<-expand.grid(t,s) 
  sigma1<-matrix(.5*exp(-.8*abs(tjs[1]-tjs[2]))[,1],ncol=J)
  sigma2<-matrix(.4*exp(-.6*abs(tjs[1]-tjs[2]))[,1],ncol=J)
  x1<-fdata(mvrnorm(n=n,sin(2*pi * tj), Sigma=sigma1),tj,rtt)
  x2<-fdata(mvrnorm(n=n,rep(0,len=J), Sigma=sigma1),tj,rtt)  #mov browniano!
#  plot(x1);  plot(x2)
	v1 <- sqrt(2) * sin(pi * tj)
	b1 <- 2 / sqrt(2)
	bet<-bet1 <- fdata(b1 * v1 *(t*(1-t)),tj,rtt)
	S1<-5*norm.fdata(x)[,1]^nls####################inprod.fdata(x1,bet1)^nls  #as.numeric(norm.fdata(x)^nls#
  varS<-var(drop(S1))*snr
  yp <-   S1#	- S2 
  s2eps<-snr*var(yp)#/(1-snr)
  if (missing(rho)) {
    sd<-sqrt(varS)
   	e <- rnorm(n,sd=sd)
    rho=0             }
  else{
  #e<-simul.ar1(n, rho=rho,burn=1000,mu=0,sd=sqrt(varS))
      sd<-sqrt(s2eps)
      e<-as.vector(arima.sim(n =n,list(order = c(length(rho),0,0), ar =rho),n.start=1000,mean=0,sd=sd))
      }
  yp <- yp +e
	return(list("x"=x1,bet=bet,"x1"=x1,"x2"=x2,"y"=yp,"e"=e,"bet1"=bet1,"bet2"=bet2,"snr"=snr,rho=rho,sd=sd))
}
#    fdatos.a<- fdata.nlm(n,J,tt=tt,rtt=range(tt),rho[j],snr=r2[i],nls=2)
 #     fx<-fdatos.a$x     

# a
# ar(fregre.np(aa$x,aa$y)$residuals,order.max=1)
# ar(fregre.pc(aa$x,aa$y)$residuals,order.max=1)
# plot(fregre.pc(aa$x,aa$y)$beta.est)

 #####################################################################################
# The function fdata.nls()   generates simulated data of the functional
# non-linear model with scalar response.
#####################################################################################
 fdata.nls<- function(n=1000,J=101,tt=seq(0,1,len=J),rtt=range(tt),rho=0,snr,nls=1){
	t<-s<-tj <- tt
  tjs<-expand.grid(t,s) 
  k3<-5
  x1<-k3*rproc2fdata(n,tt,sigma="OU")      
	S1<-norm.fdata(x1)[,1]^nls
  varS<-var(drop(S1))
  yp <-  S1
#  s2eps<-snr*var(yp)/(1-snr)  
  s2eps<-var(yp)
  sd<-snr*sqrt(s2eps)
  if (rho==0) {   	e <- rnorm(n,sd=sd)      }
  else{
#  e<-simul.ar1(n, rho=rho,burn=1000,mu=0,sd=sqrt(s2eps))
      e<-as.vector(arima.sim(n =n,list(order = c(length(rho),0,0), ar =rho),n.start=1000,mean=0,sd=sd))
      }
  yp <- yp +e
	return(list("x"=x1,"y"=yp,"snr"=snr,rho=rho,"residuals"=e,sd=))
}
#  fdatos.a<- fdata.nls(n,J,tt=tt,rtt=range(tt),rho[j],snr=r2[i],nls=1)
#  fregre.np(fdatos.a$x,fdatos.a$y) 

#	yp <- 3+5*inprod.fdata(x1,bet1)	- inprod.fdata(x1,bet)
#	plot(density(yp))
# res<-fregre.pc(x1,y)
# summary(res)

# An ARIMA simulation
#ts.sim <- arima.sim(list(order = c(1,1,0), ar = 0.7), n = 200)
#ts.plot(ts.sim)
#####################################################################################
# The function fdata.carmack()   generates simulated data of the functional
# linear linear model with scalar response, Carmack, Spence and Schucany (2012).
#####################################################################################
 fdata.carmack<- function(n=1000,J=101,tt=seq(0,1,len=J),rtt=range(tt),rho=0,snr,model="M1"){
	t<-s<-tj <- tt
  tjs<-expand.grid(t,s)
  
#  e<-outer(tt,tt, function(u,v) {(theta*exp(-abs(u-v)/scala))})
#e2<-outer(tt,tt, function(u,v) {(theta/2*exp(-abs(u-v)/scala))}) 
#x.test<-c(rproc2fdata(ntest,tt,mu=mu0,sigma=e),rproc2fdata(ntest,tt,mu=mu1,sigma=e2))
  theta=1/10
  scala=.8
  k0<-1.2
  k1<-1.2
  k2<-50
  k3<-5
  sigma1<-outer(tt,tt, function(u,v) {(theta*exp(-abs(u-v)/scala))}) 
#  sigma1<-.0001
  #matrix(.00001*exp(-20*abs(tjs[1]-tjs[2]))[,1],ncol=J)
# sigma1<-matrix(.5*exp(-.8*(tjs[1]-tjs[2])^2)[,1],ncol=J)  
  x1<-switch(model,
  M1={
  print("aa")
  bet<-k2*tt^k0*(1-tt)^k1
  x1<-rproc2fdata(n,tt,mu=bet,sigma=sigma1)},
  M2={                    
  tt<-tt/2   
  bet<-tt^k0*(1-tt)^k1
  x1<-rproc2fdata(n,tt,mu=bet,sigma=sigma1)},
  M3={
   bet<-k2*(2*tt^10*(1-tt)^2+tt^2*(1-tt)^10)
  x1<-rproc2fdata(n,tt,mu=bet,sigma=sigma1)},
  M4={
    bet<-k2*ifelse(tt<1/3,.0212*exp(tt-1/3),.0212*exp(-2*(tt-1/3)))       
  x1<-5*rproc2fdata(n,tt,mu=bet,sigma=sigma1)},
  M0={
      x1<-k3*rproc2fdata(n,tt,sigma="OU")      
      })
	v1 <- sqrt(2) * sin(pi * tj)
	b1 <- 2 / sqrt(2)
  bet1 <- fdata(bet,tj,rtt)		
	S1<-norm.fdata(x1)[,1]^nls
  varS<-var(drop(S1))
  yp <-  5+S1
  s2eps<-snr*var(yp)/(1-snr)  
 s2eps<-var(yp)
  if (rho==0) {
   	e <- rnorm(n,sd=sqrt(s2eps))
       }
  else{
     e<-as.vector(arima.sim(n =n,list(order = c(length(rho),0,0), ar =rho),n.start=1000,mean=0,sd=snr*sqrt(s2eps)))
      }
  yp <- yp +e
	return(list("x"=x1,bet=bet1,"y"=yp,"bet"=bet,"snr"=snr,rho=rho,"residuals"=e))
}
   
################################################################################

