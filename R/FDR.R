#' @title False Discorvery Rate (FDR)
#' 
#' @description Compute the False Discovery Rate for a vector of p-values and alpha value.
#' 
#' @details \code{FDR} method is used for multiple hypothesis testing to correct
#' problems of multiple contrasts.\cr If \code{dep = 1}, the tests are
#' positively correlated, for example when many tests are the same contrast.
#' \cr If \code{dep < 1} the tests are negatively correlated.
#' 
#' @aliases FDR pvalue.FDR
#' 
#' @param pvalues Vector of p-values
#' @param alpha Alpha value (level of significance).
#' @param dep Parameter dependence test. By default \code{dep = 1}, direct
#' dependence between tests.
#' @return Return:
#' \itemize{
#' \item \code{out.FDR=TRUE}: If there are significative differences.
#' \item \code{pv.FDR}: p-value for False Discovery Rate test.
#' }
#' @author Febrero-Bande, M.  and Oviedo de la Fuente, M.
#' @seealso Function used in \code{\link{fanova.RPm}}
#' @references Benjamini, Y., Yekutieli, D. (2001). \emph{The control of the
#' false discovery rate in multiple testing under dependency}. Annals of
#' Statistics. 29 (4): 1165-1188. DOI:10.1214/aos/1013699998.
#' @examples
#'  p=seq(1:50)/1000
#'  FDR(p)
#'  pvalue.FDR(p)
#'  FDR(p,alpha=0.9999)
#'  FDR(p,alpha=0.9)
#'  FDR(p,alpha=0.9,dep=-1)
#' @name FDR
#' @rdname FDR
#' @export
FDR <- function(pvalues=NULL,alpha=0.95,dep=1){
 if (is.null(pvalues)) stop("No p-values entered")
	m=length(pvalues)
	if (dep<0) {const.m=sum(1/(1:m))} else {const.m=1}
	spvalues=sort(pvalues)
	FDR=(1-alpha)*(1:m)/(m*const.m)
	return(any(spvalues<FDR))
}

#' @rdname FDR
#' @export
pvalue.FDR <- function(pvalues=NULL,dep=1){
 if (is.null(pvalues)) stop("No p-values entered")
	m=length(pvalues)
	if (dep<0) {const.m=sum(1/(1:m))} else {const.m=1}
	spvalues=sort(pvalues)
	pv.FDR=min(spvalues*(m*const.m)/(1:m))
	return(pv.FDR)
}
###

