quantile_outliers_pond <- function(x, dfunc = depth.mode,
                                   nb = 200, smo = 0.05,
                                   ns = 0.01,...){
 if (!is.fdata(x)) 
   x=fdata(x)
 dat <- x[["data"]]
 tt <- x[["argvals"]]
 rtt <- x[["rangeval"]]
 n <- NROW(dat)
 m <- NCOL(dat)
 if (is.null(n) && is.null(m)) stop("ERROR IN THE DATA DIMENSIONS")
 d <- dfunc(x,...)$dep
 cuantiles <- numeric(nb)
 vv <- var(dat)
 pr <- d/sum(d)
 rr <- rep(0,m)
 for (i in 1:nb){
   bmuestra <- x[sample(1:n,size=n,replace=TRUE,prob=pr),]
   if (smo>0) {
     bmuestra[["data"]] <- bmuestra[["data"]] + mvrnorm(n=n, rr, vv*smo)
     }
   d <- dfunc(bmuestra,...)$dep
   cuantiles[i] <- quantile(d, probs=ns, type=8)
 }
return(cuantiles)
}


