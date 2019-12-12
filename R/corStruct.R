
#library(geoR)
###########################################################################
# Simula una ARMA que permite pasarle model y si el mu=noise y sd=0 
# puede utilizarse para la prediccion
simul.arma<-function(n, model,n.start=1000,mu=0,sd=1) {
    if (!is.list(model)) stop("'model' must be list")
    if (n <= 0L)         stop("'n' must be strictly positive")
    p <- length(model$ar)    
    if (p) {
        minroots <- min(Mod(polyroot(c(1, -model$ar))))
        if (minroots <= 1) 
            stop("'ar' part of model is not stationary")
    }
    q <- length(model$ma)
    if (is.na(n.start)) 
        n.start <- p + q + ifelse(p > 0, ceiling(6/log(minroots)), 
            0)
    if (n.start < p + q) 
        stop("burn-in 'n.start' must be as long as 'ar + ma'")        
    nmax<-n.start+n
    noise <- rnorm(nmax,mu,sd)
    x <- numeric(nmax)  
    maxpq<-max(p,q)
    x[1:(maxpq-1)] <- noise[1:(maxpq-1)]
    for (i in (maxpq+1):nmax)   x[i] <- sum(model$ar * x[(i-1):(i-p)]) + sum(c(1,model$ma)*noise[i:(i-q)])    
    x[(n.start+1):nmax]
    }
# xx<-simul.arma(10,list(ar = c(rho),ma=phi))
# xx<-simul.arma(10,list(ar = c(rho)),1,1,0)
# xx<-simul.arma(10,list(ar = c(rho),ma=.75),2,1,0)

###########################################################################
# Recupera los residuos de un modelo ARIMA
filter.arima<-function(model,x){
##print("entra filter arima")
    p <- length(model$ar)    
    if (p)  x <- filter(x, c(1,-model$ar),method = "convolution",sides = 2L)
    q <- length(model$ma)
    if (q) {
            x[seq_along(model$ma)] <- 0
            x <- filter(x,-model$ma, method = "recursive")
            }
#     x[seq_along(model$ar)] <- 0
##print(x)   
   x #eps[t]<-x[t-1]
}  



################################################################################
################################################################################
corUnstruc<-function(x){return(x)}
################################################################################
################################################################################
cor.AR <- function (x, order.max = 8, p=1, method = "lm")   {
    if (is.vector(x)) {
      n <- length(x)
      W0 <- diag(n)
      x <- x - mean(x)
      aa <- ar(x, order.max = order.max, intercept = FALSE, 
               demean = FALSE, method = "mle")
      if (aa$order != 0) {
        aa <- arima(x, order = c(aa$order, 0, 0), transform.pars = TRUE, 
                    include.mean = FALSE, method = "CSS")
        W0 <- toeplitz(ARMAacf(aa$coef, lag.max = c(n - 1)))
      }
      else aa$ar <- 0
      return(list(W = W0, ar = aa))
    }
    if (method == "pc") {
      x <- t(x)
      if (!is.matrix(x)) 
        x <- as.matrix(x)
      Xcen <- x
      eigenres <- svd(Xcen)
      v <- eigenres$v
      u <- eigenres$u
      d <- eigenres$d
      l <- 1
      lenl <- length(l)
      pc.fdata <- u[, l, drop = FALSE] %*% (diag(lenl) * d[l]) %*% 
        v[l, , drop = FALSE]
      aa <- ar(pc.fdata[, 1], order.max = order.max)
      n <- length(pc.fdata[, 1])
      W0 <- diag(n)
      if (aa$order != 0) 
        W0 <- toeplitz(ARMAacf(aa$ar, lag.max = n - 1))
      else aa$ar <- 0
      return(list(W = W0, ar = aa, pc = eigenres))
    }
    if (method == "lm") {
      x <- t(x)
      n <- nrow(x)
      np <- ncol(x)
      if (missing(p)) 
        p <- 1:order.max
      else {
        order.max <- p
        p <- 1:p
      }
      lenp <- order.max
      maxp <- max(p)
      aabb <- list()
      nam <- paste("x.", p, sep = "")
      aabb <- list()
      aabb[["y"]] <- as.vector(x[(order.max + 1):n, ])
      best <- lm(y ~ 1, data = aabb)
      aic <- AIC(best)
      bestp <- 0
      coeff <- 0
      for (i in 1:lenp) {
        aabb[[nam[i]]] <- as.vector(x[(order.max + 1 - i):(n - 
                                                             i), ])
        f <- as.formula(paste("y~-1+", paste(nam[1:i], collapse = "+"), 
                              sep = ""))
        aa <- lm(f, data = aabb)
        aic2 <- AIC(aa)
        if (aic2 < aic) {
          bestp <- i
          best <- aa
          aic <- aic2
        }
      }
      if (bestp != 0) {
        coeff <- best$coefficients
        best$res.lm <- aabb[1:(bestp + 1)]
        best$res.x <- x[n:(n - bestp + 1), , drop = FALSE]
        W0 <- toeplitz(ARMAacf(coeff, lag.max = n - 1))
      }
      else {
        W0 <- diag(np)
        coeff <- 0
        best$res.x <- matrix(0, ncol = ncol(x))
      }
      return(list(W = W0, lm = best))
    }
  }
################################################################################
################################################################################
cor.ARMA<-function (x, p, d = 0, q = 0, method = "lm", order.max = 1) 
{
  if (is.vector(x)) {
    n <- length(x)
    W0 <- diag(n)
    if (missing(p)) 
      p <- 1
    aa <- arima(x, order = c(p, d, q), include.mean = FALSE, 
                transform.pars = TRUE)
    W0 <- toeplitz(ARMAacf(aa$coef, lag.max = c(n - 1)))
    return(list(W = W0, ar = aa))
  }
  if (method == "pc") {
    x <- t(x)
    if (!is.matrix(x)) 
      x <- as.matrix(x)
    Xcen <- x
    eigenres <- svd(Xcen)
    v <- eigenres$v
    u <- eigenres$u
    d <- eigenres$d
    l <- 1
    lenl <- length(l)
    pc.fdata <- u[, l, drop = FALSE] %*% (diag(lenl) * d[l]) %*% 
      v[l, , drop = FALSE]
    aa <- ar(pc.fdata[, 1], order.max = order.max)
    n <- length(pc.fdata[, 1])
    W0 <- diag(n)
    if (aa$order != 0) 
      W0 <- toeplitz(ARMAacf(aa$ar, lag.max = n - 1))
    else aa$ar <- 0
    return(list(W = W0, ar = aa, pc = eigenres))
  }
  if (method == "lm") {
    x <- (x)
    n <- nrow(x)
    np <- ncol(x)
    if (missing(p)) 
      p <- 1:order.max
    else {
      order.max <- p
      p <- 1:p
    }
    lenp <- order.max
    maxp <- max(p)
    aabb <- list()
    nam <- paste("x.", p, sep = "")
    aabb[["y"]] <- as.vector(x[(maxp + 1):n, ])
    for (i in 1:maxp) aabb[[nam[i]]] <- as.vector(x[(maxp + 
                                                       1 - i):(n - i), ])
    f <- as.formula(paste("y~-1+", paste(nam[p], collapse = "+"), 
                          sep = ""))
    aa <- lm(f, data = aabb)
    coeff <- aa$coefficients
    aa$res.lm <- aabb[1:(maxp + 1)]
    aa$res.x <- x[n:(n - maxp + 1), , drop = FALSE]
    W0 <- toeplitz(ARMAacf(coeff, lag.max = n - 1))
    return(list(W = W0, lm = aa))
  }
}
#######################################################################


###############################################################################
###############################################################################
corExpo<-function(xy,range, method = "euclidean",p=2){
 if (is.data.frame(xy)) xy<-as.matrix(xy)
 if (class(xy)=="dist") dxy<-as.matrix(xy)
 else dxy<- as.matrix(dist(xy,method=method,p=p,diag =TRUE, upper = TRUE))
 vdxy<-as.vector(dxy)
if (missing(range)) range<-quantile(vdxy,.9)/3
else range/3
  W0=exp(-dxy/range)
  return(W0)
 }


# 2016/10/31 se elimina cor.Exp y corCloud 
###############################################################################
