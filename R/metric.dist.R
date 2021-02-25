#' @title Distance Matrix Computation
#' 
#' @description This function computes the distances between the rows of a data matrix by
#' using the specified distance measure.
#' 
#' This function returns a distance matrix by using \code{\link{dist}}
#' function. \cr The matrix dimension is (\code{n1} x \code{n1}) if
#' \code{y=NULL}, (\code{n1} x \code{n2}) otherwise.
#' 
#' @param x Data frame 1. The dimension is (\code{n1} x \code{m}).
#' @param y Data frame 2. The dimension is (\code{n2} x \code{m}).
#' @param method The distance measure to be used. This must be one of
#' "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski".
#' @param p The power of the Minkowski distance.
#' @param dscale If scale is a numeric, the distance matrix is divided by the
#' scale value. If scale is a function (as the mean for example) the distance
#' matrix is divided by the corresponding value from the output of the
#' function.
#' @param \dots Further arguments passed to \code{\link{dist}} function.
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' @seealso See also \code{\link{dist}} for multivariate date case and
#' \code{\link{metric.lp} for functional data case}
#' @keywords cluster
#' @examples 
#' \dontrun{
#' data(iris)
#' d<-metric.dist(iris[,1:4])
#' matplot(d,type="l",col=as.numeric(iris[,5]))
#' }
#' @export
metric.dist <- function(x,y=NULL,method="euclidean",p=2,dscale=1,...){
if (is.vector(x)) x<-matrix(x,nrow=1)
else x<-data.matrix(x)
ynull<-is.null(y)
if (method=="mahalanobis"){
    if (ynull)   {
        y <- x
        vc <- var(x)
        }
    else {
    y <- data.matrix(y)
    vc <- var(rbind(x, y))
    }
    n <- nrow(x)
    m <- nrow(y)
    mdist<- matrix(0, n, m)
    for (i in 1:m) {
        mdist[, i] <- mahalanobis(x, y[i, ], cov = vc)
    } 
    mdist<-sqrt(mdist)
}
else{
 if (!ynull) {    
    if (is.vector(y)) y<-matrix(y,nrow=1) 
    n<-nrow(y)
    nn<-nrow(x)
    mdist<-as.matrix(dist(rbind(x,y) , method = method, diag = TRUE, upper = TRUE,p=p))[1:nn,(nn+1):(nn+n)] 
    }
 else   mdist<-as.matrix(dist(x, method = method, diag = TRUE, upper = TRUE,p=p))  
 }
 if (is.vector(mdist)) mdist<-matrix(mdist,nrow=nn)   
 	etiq1=rownames(x)
	etiq2=rownames(y)
# namesx<-rownames(x)
# if (ynull) dimnames(mdist) <- list(namesx,namesx)
# else dimnames(mdist) <- list(namesx, rownames(y))
 if (is.function(dscale)) {
   if (nrow(mdist)==ncol(mdist)) diag(mdist)<-NA ################# ojjjooo solo para matrices cuadradas
   dscale<-dscale(as.dist(mdist))
   if (nrow(mdist)==ncol(mdist)) diag(mdist)<-0   
 } 
 mdist<-mdist/dscale
 attr(mdist, "call") <- "metric.dist"
 attr(mdist, "par.metric") <- list(method =method,p=p,dscale=dscale) 
 rownames(mdist)<-etiq1
 colnames(mdist)<-etiq2
 return(mdist)
}
#####################
