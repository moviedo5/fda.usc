#' Distance Correlation Statistic and t-Test  
#' 
#' @description Distance correlation t-test of multivariate and functional 
#' independence (wrapper functions of energy package).
#' 
#' @details  These wrapper functions extend the functions of the \code{energy} package
#' for multivariate data to functional data.  Distance correlation is a measure
#' of dependence between random vectors introduced by Szekely, Rizzo, and
#' Bakirov (2007).
#' \code{dcor.xy} performs a nonparametric t-test of multivariate or functional
#' independence in high dimension. The distribution of the test statistic is
#' approximately Student t with \eqn{n(n-3)/2-1} degrees of freedom and for
#' \eqn{n \geq 10} the statistic is approximately distributed as standard
#' normal.  Wrapper function of \code{energy:::dcor.ttest}.  The t statistic is
#' a transformation of a bias corrected version of distance correlation (see SR
#' 2013 for details). Large values (upper tail) of the t statistic are
#' significant.\cr 
#' \code{dcor.test} similar to \code{dcor.xy} but only for distance matrix.
#' \code{dcor.dist} compute distance correlation statistic.  Wrapper function
#' of \code{energy::dcor} but only for distance matrix
#' \code{bcdcor.dist} compute bias corrected distance correlation statistic.
#' Wrapper function of \code{energy:::bcdcor} but only for distance matrix.
#' 
#' @aliases dcor.test dcor.dist bcdcor.dist dcor.xy 
#' 
#' @param x data (fdata, matrix or data.frame class) of first sample.
#' @param y data (fdata, matrix or data.frame class) of second sample.
#' @param test if TRUE, compute bias corrected distance correlation statistic
#' and the corresponding t-test, else compute distance correlation statistic.
#' @param metric.x,metric.y Name of metric or semi-metric function used for
#' compute the distances of \code{x} and \code{y} object respectively. By
#' default, \code{\link{metric.lp}} for functional data and
#' \code{\link{metric.dist}} for multivariate data.
#' @param par.metric.x,par.metric.y List of parameters for the corresponding
#' metric function.
#' @param n The sample size used in bias corrected version of distance
#' correlation, by default is the number of rows of \code{x}.
#' @param D1 Distances of first sample (x data).
#' @param D2 Distances of second sample (y data).
#' 
#' @return 
#' \code{dcor.test} returns a list with class \code{htest} containing
#' \itemize{
#' \item \code{method}{ description of test} 
#' \item \code{statistic}{ observed value of the test statistic} 
#' \item \code{parameter}{ degrees of freedom} 
#' \item \code{estimate}{ bias corrected distance correlation \code{bcdcor(x,y)}} 
#' \item \code{p.value}{ p-value of the t-test} 
#' \item \code{data.name}{ description of data}
#' }
#' 
#' \code{dcor.xy} returns the previous list with class \code{htest} and 
#' \itemize{
#' \item \code{D1}{ the distance matrix of \code{x}} 
#' \item \code{D2}{ the distance matrix of \code{y}}
#' }
#' 
#' \code{dcor.dist} returns the distance correlation statistic.
#' 
#' \code{bcdcor.dist} returns the bias corrected distance correlation
#' statistic.
#' 
#' @author Manuel Oviedo de la Fuente \email{manuel.oviedo@@udc.es} and Manuel
#' Febrero Bande
#' 
#' @seealso \code{\link{metric.lp}} amd \code{\link{metric.dist}}.
#' @references Szekely, G.J. and Rizzo, M.L. (2013). The distance correlation
#' t-test of independence in high dimension. \emph{Journal of Multivariate
#' Analysis}, Volume 117, pp. 193-213. 
#' 
#' Szekely, G.J., Rizzo, M.L., and Bakirov, N.K. (2007), Measuring and Testing
#' Dependence by Correlation of Distances, \emph{Annals of Statistics}, Vol. 35
#' No. 6, pp. 2769-2794.  
#' 
#' @keywords htest multivariate nonparametric
#' @examples
#' \dontrun{ 
#' x<-rproc2fdata(100,1:50)
#' y<-rproc2fdata(100,1:50)
#' dcor.xy(x, y,test=TRUE)
#' dx <- metric.lp(x)
#' dy <- metric.lp(y)
#' dcor.test(dx, dy)
#' bcdcor.dist(dx, dy)
#' dcor.xy(x, y,test=FALSE)
#' dcor.dist(dx, dy)
#' }
#' 
#' @rdname dcor.xy
#' @export 
dcor.xy<-function (x, y,test=TRUE,metric.x,metric.y,par.metric.x,par.metric.y,n)
{
  # print("entra dcor.xy")
  dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  if (is.vector(x)) x<-matrix(x,ncol=1)
  if (is.vector(y)) y<-matrix(y,ncol=1)
  if (is.factor(y)) y<-model.matrix(~y-1)
  if (is.factor(x)) x<-model.matrix(~x-1)

  #   if (anyNA(x)) {
  #   nas <- !apply(x,1,is.na)
  #   x<-x[nas,,drop=F]
  #   y<-y[nas,,drop=F]
  #   warning("NA's in x data")
  # }
  # if (anyNA(y)) {
  #   nas <- !apply(y,1,is.na)
  #   x<-x[nas,,drop=F]
  #   y<-y[nas,,drop=F]
  #   warning("NA's in y data")
  # }

  if (missing(n)) n <- nrow(x)
  isfdata1<-is.fdata(x)
  true.dist<-FALSE
  if (missing(metric.x)){
    if (isfdata1) metric.x="metric.lp"
    else metric.x="metric.dist"
   }
  if (missing(par.metric.x)) par.metric.x<-list()
  if (missing(n)) n <- nrow(x)
  isfdata2<-is.fdata(y)
  if (isfdata1)         par.metric.x$fdata1<-x
  else                  par.metric.x$x<-as.matrix(x)
  D1=do.call(metric.x,par.metric.x)
  if (missing(metric.y)){
          if (isfdata2) metric.y="metric.lp"
          else metric.y="metric.dist"
  }
  if (missing(par.metric.y)) par.metric.y<-list()
  if (isfdata2)        par.metric.y$fdata1<-y
  else                 par.metric.y$x<-as.matrix(y)
  D2=do.call(metric.y,par.metric.y)
  D1<-as.matrix(D1)
  D2<-as.matrix(D2)
  if (test) {
   rval<-dcor.test(D1, D2,n)
   rval$D1<-D1
   rval$D2<-D2  
   }
  else    {
   rval<- dcor.dist(D1=D1,D2=D2) 
   names(rval)< "dcor" 
   }
  return(rval)
} 



#' @rdname dcor.xy
#' @export 
dcor.dist=function(D1,D2){
  # Description
  # Computes distance correlation for functional data.
  #
  # Arguments
  # D1: distances of 1st sample
  # D2: distances of 2nd sample 
  # Returns the sample distance correlation
  # @rdname dcor.xy
  m1row=rowMeans(D1)
  m1col=colMeans(D1)
  m2row=rowMeans(D2)
  m2col=colMeans(D2)
  n=nrow(D1)
  ones=rep(1,n)
  pD1=D1-outer(ones,m1row)-outer(m1col,ones)+mean(D1)
  pD2=D2-outer(ones,m2row)-outer(m2col,ones)+mean(D2)
  res=sqrt(sum(pD1*pD2))/n
  res1=sqrt(sum(pD1*pD1))/n        
  res2=sqrt(sum(pD2*pD2))/n        
  out<-res/sqrt(res1*res2)
  return(out)
}

#' @rdname dcor.xy
#' @export
bcdcor.dist=function(D1,D2,n){
  # Description
  # Computes the bias corrected distance correlation for functional data.
  #
  # Arguments
  # D: distances of 1st sample
  # D2: distances of 2nd sample 
  # n:  sample dimension, by default nrow(D1)
  # Returns the bias corrected dcor statistic
  # print("entra bcdcor.dist")
  # print(n)
  # print(dim(D1))
  # print(dim(D2))
  if (missing(n)) n <- nrow(D1)
  AA <- Astar2(D1,n)
  BB <- Astar2(D2,n)
  res <- sum(AA * BB) - (n/(n - 2)) * sum(diag(AA * BB))
  res1<- sum(AA * AA) - (n/(n - 2)) * sum(diag(AA * AA))
  res2<- sum(BB * BB) - (n/(n - 2)) * sum(diag(BB * BB))
  out<-res/sqrt(res1*res2)          
  #print("sale bcdcor.dist")
  return(out)
}

#' @rdname dcor.xy
#' @export 
dcor.test <- function (D1, D2,n)
{
  dname <- paste(deparse(substitute(D1)), "and", deparse(substitute(D2)))
  if (missing(n)) n <- nrow(D1)
  R<-bcdcor.dist(D1=D1,D2=D2,n=n) 
  M <- n * (n - 3)/2
  df <- M - 1
  names(df) <- "df"
  tstat <- sqrt(M - 1) * R/sqrt(1 - R^2)
  names(tstat) <- "T"
  names(R) <- "Bias corrected dcor"
  pval <- 1 - pt(tstat, df = df)
  method <- "dcor t-test of independence"
  rval <- list(statistic = tstat, parameter = df, p.value = pval,
               estimate = R, method = method, data.name = dname)
  class(rval) <- "htest"
  return(rval)
}


Astar2<-function(d,n){
  #  d <- as.matrix(d)
  m <- rowMeans(d)
  M <- mean(d)
  a <- sweep(d, 1, m)
  b <- sweep(a, 2, m)
  A <- b + M
  A <- A - d/n
  diag(A) <- m - M
  (n/(n - 1)) * A
}

dcor.y <- function(ldist,response,n,bcdcor=TRUE){
  lenldist<-length(ldist)
  namldist<-names(ldist)
  if (missing(n)) {
    if (is.fdata(ldist[[1]]) | is.matrix(ldist[[1]]|is.data.frame(ldist[[1]])) )
      n<-nrow(ldist[[1]])
    if (is.vector(ldist[[1]]))    n <- length(response)
  }
  
  if (missing(response)) {#print("se calculan todas las distancia")
    dst<-diag(lenldist)
    for (i in 1:(lenldist-1)) {
      for (j in i:(lenldist)) {
        if (bcdcor)     dst[i,j]<-dst[j,i]<- bcdcor.dist(ldist[[i]],ldist[[j]],n=n)
        else            dst[i,j]<-dst[j,i]<-cor.fdist(ldist[[i]],ldist[[j]])
      }}
    
    rownames(dst)<-colnames(dst)<-namldist
  }
  else{                     #se calculan todas las distancias respecto la respuesta
    if (length(response)>1) stop("Incorrect response label")
    ii<-response==namldist
    if (sum(ii)==0) stop("Incorrect response label")
    iresponse<-which(ii)
    dst<-numeric(lenldist-1)
    predictors<-setdiff(namldist,response)
    
    for (i in 1:(lenldist-1)) {
      #           dst[i]<-cor.fdist(ldist[[response]],ldist[[predictors[i]]])
      if (bcdcor)  dst[i]<-bcdcor.dist(ldist[[response]],ldist[[predictors[i]]],n=n)
      else dst[i]<-dcor.dist(ldist[[response]],ldist[[predictors[i]]])
    }
    names(dst)<-predictors
  }
  dst
}


cor.fdist=function(D1,D2=NULL,...){
  # Description
  # Computes distance correlation for functional data.
  #
  # Arguments
  # D1: class(D1)=matrix: distances of first sample,  class(D1)=fdata: first fdata data sample
  # D2: class(D2)=matrix: distances of second sample, class(D2)=fdata: second fdata data sample
  # Returns the sample distance correlation
  if (is.null(D2)) {
    if (is.fdata(D1)) D1=metric.dist(D1,...)
    m1row=rowMeans(D1)#apply(D1,1,mean)
    m1col=colMeans(D1)#apply(D1,2,mean)
    n=nrow(D1)
    ones=rep(1,n)
    pD1=D1-outer(ones,m1row)-outer(m1col,ones)+mean(D1)
    res=sqrt(sum(pD1*pD1))/n
    out<-res/sqrt(res*res) 
  }
  else {
    if (is.fdata(D1)) D1=metric.lp(D1,...)      
    if (is.fdata(D2)) D2=metric.lp(D2,...)                    
    m1row=rowMeans(D1) 
    m2row=rowMeans(D2)
    m1col=colMeans(D1)
    m2col=colMeans(D2)
    n=nrow(D1)
    ones=rep(1,n)
    pD1=D1-outer(ones,m1row)-outer(m1col,ones)+mean(D1)
    pD2=D2-outer(ones,m2row)-outer(m2col,ones)+mean(D2)
    res=sqrt(sum(pD1*pD2))/n
    res1=sqrt(sum(pD1*pD1))/n        
    res2=sqrt(sum(pD2*pD2))/n        
    out<-res/sqrt(res1*res2)
  }
  return(out)
}
