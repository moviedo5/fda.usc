#' The cross-validation (CV) score
#' 
#' Compute the leave-one-out cross-validation score.
#' 
#' @details A.-If \code{trim=0}:\cr
#' \deqn{CV(h)=\frac{1}{n}
#' \sum_{i=1}^{n}{\Bigg(\frac{y_i-r_{i}(x_i)}{(1-S_{ii})}\Bigg)^{2}w(x_{i})}}{CV(h)=1/n\,
#' \sum_i ((y_i\, -\, r_{i}(x_i))\, /\, (1\, -\, S_ii))^2\, w(x_i),\,
#' i=1,...,n} \eqn{S_{ii}}{Sii} is the ith diagonal element of the smoothing
#' matrix \eqn{S}{S}.
#' 
#' B.-If \code{trim>0}:\cr \deqn{CV(h)=\frac{1}{l}
#' \sum_{i=1}^{l}{\Bigg(\frac{y_i-r_{i}(x_i)}{(1-S_{ii})}\Bigg)^{2}w(x_{i})}}{CV(h)=1/n\
#' \sum_i ((y_i-r_{i}(x_i))/(1-S_ii))^2 w(x_i),\, i=1,...,l} \eqn{S_{ii}}{Sii}
#' is the ith diagonal element of the smoothing matrix \eqn{S}{S} and l the
#' index of \code{(1-trim)} curves with less error.
#' 
#' @param y Matrix of set cases with dimension (\code{n} x \code{m}), where
#' \code{n} is the number of curves and \code{m} are the points observed in
#' each curve.
#' @param S Smoothing matrix, see \code{\link{S.NW}}, \code{\link{S.LLR}} or
#' \eqn{S.KNN}.
#' @param W Matrix of weights.
#' @param trim The alpha of the trimming.
#' @param draw =TRUE, draw the curves, the sample median and trimmed mean.
#' @param metric Metric function, by default \code{\link{metric.lp}}.
#' @param \dots Further arguments passed to or from other methods.
#' @return { Returns CV score calculated for input parameters.  }
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' @seealso See Also as \code{\link{optim.np}} \cr Alternative method:
#' \code{\link{GCV.S}}
#' @references Wasserman, L. \emph{All of Nonparametric Statistics}. Springer
#' Texts in Statistics, 2006.
#' @keywords utilities
#' @examples
#' \dontrun{
#' data(tecator)
#' x<-tecator$absorp.fdata
#' np<-ncol(x)
#' tt<-1:np
#' S1 <- S.NW(tt,3,Ker.epa)
#' S2 <- S.LLR(tt,3,Ker.epa)
#' S3 <- S.NW(tt,5,Ker.epa)
#' S4 <- S.LLR(tt,5,Ker.epa)
#' cv1 <- CV.S(x, S1)
#' cv2 <- CV.S(x, S2)
#' cv3 <- CV.S(x, S3)
#' cv4 <- CV.S(x, S4)
#' cv5 <- CV.S(x, S4,trim=0.1,draw=TRUE)
#' cv1;cv2;cv3;cv4;cv5
#' S6 <- S.KNN(tt,1,Ker.unif,cv=TRUE)
#' S7 <- S.KNN(tt,5,Ker.unif,cv=TRUE)
#' cv6 <- CV.S(x, S6)
#' cv7 <- CV.S(x, S7)
#' cv6;cv7
#' }
#'  
#' @export CV.S
CV.S=function (y, S, W = NULL, trim = 0, draw = FALSE, metric = metric.lp, ...) {
    n = ncol(S)
    isfdata <- is.fdata(y)
    if (isfdata) {
        nn<-nrow(y)
        if (is.null(W)) W<-diag(nn)
        y2 = t(y$data)
        y.est = t(S %*% y2)
        y.est <- fdata(y.est, y$argvals, y$rangeval, y$names)
        #e <- (      y - y.est)/(1-diag(S))
        e <- fdata(sweep((y - y.est)$data,2,1-diag(S),"%/%"), y$argvals, y$rangeval, y$names)
        ee <- drop(norm.fdata(e, metric = metric, ...)[, 1]^2)
        if (trim > 0) {
            e.trunc = quantile(ee, probs = (1 - trim), na.rm = TRUE, type = 4)
            ind <- ee <= e.trunc
            if (draw)   plot(y, col = (2 - ind))
            res = mean(ee[ind])
        }
        else res = mean(ee)
    }
    else {
       if (is.null(W)) W<-diag(n)
        y2 <- y
        y.est = S %*% y2
        I = diag(n)/(1 - diag(S))^2
        W = W * I
        e <- y2 - y.est
        if (trim > 0) {
            ee = t(e)
            e.trunc = quantile(abs(ee), probs = (1 - trim), na.rm = TRUE, 
                type = 4)
            l <- which(abs(ee) <= e.trunc)
            res = mean(diag(W)[l] * e[l]^2)
        }
        res = mean(diag(W) * e^2)
    }
    if (is.nan(res))    res = Inf
    return(res)
}

