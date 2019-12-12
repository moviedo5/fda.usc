#' @title Compute the Hausdorff distances between two curves.
#' 
#' @description Hausdorff distance is the greatest of all the distances from a point in one
#' curve to the closest point in the other curve (been closest the euclidean
#' distance).
#' 
#' @details Let \eqn{G(X)=\left\{ (t,X(t))\in R^2 \right\}}{G(X)={(t,X(t)) \in R^2}} and
#' \eqn{G(Y)=\left\{(t,Y(t))\in R^2\right\}}{G(Y)={(s,Y(s)) \in R^2}} be two
#' graphs of the considered curves \eqn{X} and \eqn{Y} respectively, the
#' Hausdorff distance \eqn{d_H(X, Y)} is defined as,
#' 
#' \deqn{ d_H(X,Y)=max\left\{ sup_{x\in G(X)} inf_{y\in G(Y)} d_2(x,y),
#' sup_{y\in G(Y)} inf_{x\in G(X)}d_2(x,y)\right\},}{d_H(X,Y)=max{sup_{x\in
#' G(X)} inf_{y\in G(Y)}d_2(x,y),sup_{y\in G(Y)} inf_{x\in G(X)}d_2(x,y)},}
#' where \eqn{d_2(x,y)} is the euclidean distance, see \code{\link{metric.lp}.}
#' 
#' @param fdata1 Curves 1 of \code{fdata} class. The dimension of \code{fdata1}
#' object is (\code{n1} x \code{m}), where \code{n1} is the number of points
#' observed in \code{t} coordinates with lenght \code{m}.
#' @param fdata2 Curves 2 of \code{fdata} class. The dimension of \code{fdata2}
#' object is (\code{n2} x \code{m}), where \code{n2} is the number of points
#' observed in \code{t} coordinates with lenght \code{m}.
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@usc.es}
#' @keywords cluster
#' @examples
#' \dontrun{   
#' data(poblenou)
#' nox<-poblenou$nox[1:6]
#' # Hausdorff vs maximum distance
#' out1<-metric.hausdorff(nox)       
#' out2<-metric.lp(nox,lp=0) 
#' out1
#' out2
#' par(mfrow=c(1,3))
#' plot(nox)
#' plot(hclust(as.dist(out1)))
#' plot(hclust(as.dist(out2)))
#' }   
#'
#' @rdname metric.hausdorff
#' @export
metric.hausdorff=function (fdata1, fdata2 = fdata1) 
{
    if (is.fdata(fdata1)) {
        tt <- fdata1[["argvals"]]
        rtt <- fdata1[["rangeval"]]
        nas1 <- is.na.fdata(fdata1)
        if (any(nas1)) 
            stop("fdata1 contain ", sum(nas1), " curves with some NA value \n")
        else if (!is.fdata(fdata2)) {
            fdata2 <- fdata(fdata2, tt, rtt)
        }
        nas2 <- is.na.fdata(fdata2)
        if (any(nas2)) 
            stop("fdata2 contain ", sum(nas2), " curves with some NA value \n")
        DATA1 <- fdata1[["data"]]
        DATA2 <- fdata2[["data"]]
        range.t <- rtt
    }
    else {
        if (is.vector(fdata1)) 
            fdata1 <- as.matrix(t(fdata1))
        if (is.vector(fdata2)) 
            fdata2 <- as.matrix(t(fdata2))
        DATA1 <- fdata1
        DATA2 <- fdata2
        range.t <- c(1, ncol(DATA1))
    }

    testfordim <- sum(dim(DATA1) == dim(DATA2)) == 2
    twodatasets <- TRUE
    if (testfordim) 
        twodatasets <- sum(DATA1 == DATA2) != prod(dim(DATA1))
    numgr1 <- nrow(DATA1)
	numgr2 <- nrow(DATA2)
	Mtt=outer(tt,tt,"-")^2
	mdist=array(0,dim=c(numgr1,numgr2))
	etiq1=rownames(DATA1)
	etiq2=rownames(DATA2)
	predi=TRUE
	if (testfordim) {
                if (twodatasets) {
                  predi <- FALSE
                  for (i in 1:numgr1) {
                    ii = i + 1
                    for (ii in i:numgr2) {
					  Mxy=sqrt(outer(DATA1[i,],DATA2[ii,],"-")^2+Mtt)
					  mdist[i, ii] = max(max(apply(Mxy,1,min)),max(apply(Mxy,2,min)))
					  }
                  }
                  mdist = t(mdist) + mdist
                }
            }
            if (predi) {
                for (i in 1:numgr1) {
                  for (ii in 1:numgr2) {
					  Mxy=sqrt(outer(DATA1[i,],DATA2[ii,],"-")^2+Mtt)
					  mdist[i, ii] = max(max(apply(Mxy,1,min)),max(apply(Mxy,2,min)))
                    }
                }
            }
			
	attr(mdist, "call") <- "metric.hausdorff"
	rownames(mdist)<-etiq1
	colnames(mdist)<-etiq2
    return(mdist)
}
