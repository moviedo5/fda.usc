#' @name metric.kl
#' @title Kullback--Leibler distance
#' 
#' @description Measures the proximity between two groups of densities (of class
#' \code{fdata}) by computing the Kullback--Leibler distance.
#' 
#' @details Kullback--Leibler distance between \eqn{f(t)}{f(t)} and \eqn{g(t)}{g(t)} is
#' \deqn{metric.kl(f(t),g(t))= \int_{a}^{b} {f(t) log\left(\frac{f(t)}{g(t)}\right)dt}}{dist(f(t),g(t))= \int_{a}^{b} f(t)*
#' log(f(t)/g(t)) dt} where \eqn{t} are the \code{m} coordinates of the points
#' where the density is observed (the \code{argvals} of the \code{fdata} object).
#' 
#' The Kullback--Leibler distance is asymmetric, 
#' \deqn{metric.kl(f(t),g(t))\neq metric.kl(g(t),f(t))}{dist(f(t),g(t))!=dist(g(t),f(t))}
#'  A symmetry version of K--L distance (by default) can be obtained by
#' \deqn{0.5\left(metric.kl(f(t),g(t))+metric.kl(g(t),f(t))\right)}{0.5*(dist(f(t),g(t))+dist(g(t),f(t)))}
#' 
#' If \eqn{\left(f_i(t)=0\ \& \ g_j(t)=0\right) \Longrightarrow
#' metric.kl(f(t),g(t))=0}{(f(t)=0 and g(t)=0), then metric.kl(f(t),g(t))=0}.
#' 
#' If \eqn{\left|f_i(t)\-g_i(t) \right|\leq \epsilon \Longrightarrow
#' f_i(t)=f_i(t)+\epsilon}{abs(f(t)-g(t))<\epsilon, then f(t)=f(t)+\epsilon},
#' where \eqn{\epsilon} is the tolerance value (by default \code{eps=1e-10}).
#' 
#' The coordinates of the points where the density is observed (discretization
#' points \eqn{t}) can be equally spaced (by default) or not. 
#' 
#' @param fdata1 Functional data 1 (\code{fdata} class) with the densities. The
#' dimension of \code{fdata1} object is (\code{n1} x \code{m}), where \code{n1}
#' is the number of densities and \code{m} is the number of coordinates of the
#' points where the density is observed.
#' @param fdata2 Functional data 2 (\code{fdata} class) with the densities. The
#' dimension of \code{fdata2} object is (\code{n2} x \code{m}).
#' @param symm If \code{TRUE} the symmetric K--L distance is computed, see
#' details section.
#' @param base The logarithm base used to compute the distance.
#' @param eps Tolerance value.
#' @param \dots Further arguments passed to or from other methods.
#' 
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@usc.es}
#' 
#' @seealso See also \code{\link{metric.lp}} and \code{\link{fdata}}
#' 
#' @references Kullback, S., Leibler, R.A. (1951). \emph{On information and
#' sufficiency.} Annals of Mathematical Statistics, 22: 79-86
#' 
#' @keywords cluster
#' 
#' @examples
#' \dontrun{   
#' n<-201                                                                                       
#' tt01<-seq(0,1,len=n)                                                                         
#' rtt01<-c(0,1)  
#' x1<-dbeta(tt01,20,5)                                                                           
#' x2<-dbeta(tt01,21,5)                                                                           
#' y1<-dbeta(tt01,5,20)                                                                           
#' y2<-dbeta(tt01,5,21)                                                                           
#' xy<-fdata(rbind(x1,x2,y1,y2),tt01,rtt01)
#' plot(xy)
#' round(metric.kl(xy,xy,eps=1e-5),6)  
#' round(metric.kl(xy,eps=1e-5),6)
#' round(metric.kl(xy,eps=1e-6),6)
#' round(metric.kl(xy,xy,symm=FALSE,eps=1e-5),6)  
#' round(metric.kl(xy,symm=FALSE,eps=1e-5),6)
#' 
#' plot(c(fdata(y1[1:101]),fdata(y2[1:101])))                       
#' metric.kl(fdata(x1))  
#' metric.kl(fdata(x1),fdata(x2),eps=1e-5,symm=F)       
#' metric.kl(fdata(x1),fdata(x2),eps=1e-6,symm=F)       
#' metric.kl(fdata(y1[1:101]),fdata(y2[1:101]),eps=1e-13,symm=F)  
#' metric.kl(fdata(y1[1:101]),fdata(y2[1:101]),eps=1e-14,symm=F)  
#' }   
#' 
#' @rdname metric.kl
#' @export
metric.kl=function(fdata1, fdata2 = NULL,symm=TRUE, 
                   base=exp(1),eps=1e-10,...) {                                                                                
#verbose=TRUE
    C1 <- match.call()
    same <- FALSE
    if (is.fdata(fdata1)) {
        fdat <- TRUE
        tt <- tt1 <- fdata1[["argvals"]]
        rtt <- fdata1[["rangeval"]]
        nas1 <- is.na.fdata(fdata1)
        if (any(nas1)) 
            warning("fdata1 contain ", sum(nas1), " curves with some NA value \n")
        if (is.null(fdata2)) {
            fdata2 <- fdata1
            same <- TRUE
        }
        else if (!is.fdata(fdata2)) {
            fdata2 <- fdata(fdata2, tt1, rtt, fdata1$names)
        }
        nas2 <- is.na.fdata(fdata2)
        if (any(nas2)) {
            warning("fdata2 contain ", sum(nas2), " curves with some NA value \n")
        }
        DATA1 <- fdata1[["data"]]
        DATA2 <- fdata2[["data"]]
#print(DATA1[,1:3 ])        
#print(DATA2[,1:3 ])
        tt2 <- fdata2[["argvals"]]
        if (!same) {
            if (sum(tt1 != tt2) != 0) {
                stop("Error: different discretization points in the input data.\n")
            }
        }
    }
    else {
        fdat <- FALSE
        DATA1 <- fdata1
        if (is.null(fdata2)) {
            fdata2 <- fdata1
            same <- TRUE
        }
        DATA2 <- fdata2
    }
    numgr = nrow(DATA1)
    numgr2 = nrow(DATA2)
    np <- ncol(DATA1)   
#eps2<-eps
testfordim <- sum(dim(DATA1) == dim(DATA2)) == 2
	etiq1=rownames(DATA1)
	etiq2=rownames(DATA2)
twodatasets <- TRUE
if (testfordim)   twodatasets <- (all(DATA1== DATA2)) 
#if (verbose) {print(testfordim);print(twodatasets);print(symm)}
        eps2 <- as.double(.Machine[[1]] * 100)
        inf <- 1 - eps2                                      
        sup <- 1 + eps2                                      
        zz<-apply(DATA1,1,sum)
        if (all(zz > inf) & all(zz < sup)) {                
            densi=TRUE
        }                                                     
        else {
            densi = FALSE                                     
          warnings("Ponderation fdata1/sum(fdata1) is done")
          DATA1<-DATA1/apply(DATA1,1,sum)  #calcular integral
          }
        zz<-apply(DATA2,1,sum)          
         if (all(zz > inf) & all(zz < sup)) {   densi=TRUE        }                                                                       
         else {                                                                  
           densi = FALSE                                                       
           warnings("Ponderation fdata2/sum(fdata2) is done")       
           DATA2<-DATA2/apply(DATA2,1,sum)  #calcular integral                                     
           }                                                                      
     #& symm
#if (verbose) {print(testfordim)        ;print(twodatasets)}

   if (fdat) {
        dtt <- diff(tt)
        inf <- dtt - eps2
        sup <- dtt + eps2
        if (all(dtt > inf) & all(dtt < sup)) {
            equi = TRUE
        }
        else equi = FALSE
        mdist = array(0, dim = c(numgr, numgr2))
        predi <- TRUE
        if (testfordim) {
            if (twodatasets) {  
                predi <- FALSE
#if (verbose) {print(" no predi")  ;print("symm")  ;print(symm)  }
if (symm){              
                for (i in 1:numgr) {
                  ii = i + 1           
                  for (ii in i:numgr2) {
a4<-DATA1[i, ]  
a44<-DATA2[ii, ]
#      a4[a4<eps & a44<eps]<-min(eps2,min(a4[a4>0]))              
#                    a44[a4<eps & a44<eps]<-min(eps2,min(a44[a44>0]))  
if (any((a4<eps & a44>(eps+a4))|(a4>(eps+a44) & a44<eps)))   {
#          print("1")
          mdist[i, ii]<-Inf }
else  {
                    a4[a4<eps & a44<eps]<-min(eps2,min(a4[a4>0]))              
                    a44[a4<eps & a44<eps]<-min(eps2,min(a44[a44>0]))              
                    a4<-a4/sum(a4)                  
                    a44<-a44/sum(a44)
                    f = a4*log(a4/a44,base=base)
                    d1 = (int.simpson2(tt, f, equi))
                    f = a44*log(a44/a4,base=base)
                    d2= (int.simpson2(tt, f, equi))
                    mdist[i, ii] <-(d1+d2)/2    
  }                      
                  }
#print(                   mdist[i,] )
                }
              mdist = t(mdist) + mdist
                }
            else{
                for (i in 1:numgr) {
                  ii = i + 1
                  for (ii in 1:numgr2) {
a4<-DATA1[i, ]  
a44<-DATA2[ii, ] 
if (any((a4<eps & a44>(eps+a4))|(a4>(eps+a44) & a44<eps)))   {
#print("2")
mdist[i, ii]<-Inf }
else  { 
             a4[a4<eps & a44<eps]<-min(eps2,min(a4[a4>0]))              
             a44[a4<eps & a44<eps]<-min(eps2,min(a44[a44>0]))  
             a4<-a4/sum(a4)                  
             a44<-a44/sum(a44)                    
                    f <- a4*log(a4/a44,base=base)
                    mdist[i, ii]=(int.simpson2(tt, f, equi))                                     
#print( mdist[i, ii])                    
                  
                  }
#print(                   mdist[i,] )                  
                }
              }   }             

                
            }
        }
        if (predi) {
           if (symm){
            for (i in 1:numgr) {
                for (ii in 1:numgr2) {
a4<-DATA1[i, ]  
a44<-DATA2[ii, ]     
if (any((a4<eps & a44>(eps+a4))|(a4>(eps+a44) & a44<eps)))   {
#print(3)
mdist[i, ii]<-Inf }
else  { 
             a4[a4<eps & a44<eps]<-min(eps2,min(a4[a4>0]))              
             a44[a4<eps & a44<eps]<-min(eps2,min(a44[a44>0]))  
             a4<-a4/sum(a4)                  
             a44<-a44/sum(a44)                    
                    f = a4*log(a4/a44,base=base)
                    d1 = (int.simpson2(tt, f, equi))
                     f = a44*log(a44/a4,base=base)
                    d2= (int.simpson2(tt, f, equi))
                    mdist[i, ii] <-(d1+d2)/2
                }
            }
                        }}
           else{     
            for (i in 1:numgr) {
                for (ii in 1:numgr2) {
a4<-DATA1[i, ]  
a44<-DATA2[ii, ] 
if (any((a4<eps & a44>(eps+a4))|(a4>(eps+a44) & a44<eps)))   {
#print(4)
mdist[i, ii]<-Inf }
else  {    
                    a4[a4<eps & a44<eps]<-min(eps2,min(a4[a4>0]))              
                    a44[a4<eps & a44<eps]<-min(eps2,min(a44[a44>0]))  
 f = a4*log(a4/a44,base=base)
                    mdist[i, ii] = (int.simpson2(tt, f, equi))
}                     
                }
            }
          }
        }       
    }
    else {
   stop("No fdata class object")
     }
	rownames(mdist)<-etiq1
	colnames(mdist)<-etiq2     
    attr(mdist, "call") <- "metric.kl"
    attr(mdist, "par.metric") <- list( base=base,symm=symm)
    return(mdist)
}

