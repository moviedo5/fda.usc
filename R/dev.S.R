#' @title The deviance score
#' 
#' @description Returns the deviance of a fitted model object by GCV score.
#' 
#' @details Up to a constant, minus twice the maximized log-likelihood. Where sensible,
#' the constant is chosen so that a saturated model has deviance zero.
#' 
#' \deqn{GCV(h)=p(h) \Xi(n^{-1}h^{-1})}{} 
#' 
#' Where \cr \deqn{p(h)=\frac{1}{n}
#' \sum_{i=1}^{n}{\Big(y_i-r_{i}(x_i)\Big)^{2}w(x_i)}}{} \cr and penalty
#' function \deqn{\Xi()}{} can be selected from the following criteria: \cr
#' 
#' Generalized Cross-validation (GCV):
#' 
#' \deqn{\Xi_{GCV}(n^{-1}h^{-1})=(1-n^{-1}S_{ii})^{-2}}{} \cr Akaike's
#' Information Criterion (AIC): 
#' 
#' \deqn{\Xi_{AIC}(n^{-1}h^{-1})=exp(2n^{-1}S_{ii})}{} 
#' Finite Prediction Error (FPE) 
#' 
#' \deqn{\Xi_{FPE}(n^{-1}h^{-1})=\frac{(1+n^{-1}S_{ii})}{(1-n^{-1}S_{ii})}}{}
#' 
#' Shibata's model selector (Shibata): 
#' 
#' \deqn{\Xi_{Shibata}(n^{-1}h^{-1})=(1+2n^{-1}S_{ii})}{}
#' Rice's bandwidth selector (Rice):
#' \deqn{\Xi_{Rice}(n^{-1}h^{-1})=(1-2n^{-1}S_{ii})^{-1}}{}
#' 
#' @param y Matrix of set cases with dimension (\code{n} x \code{m}), where
#' \code{n} is the number of curves and \code{m} are the points observed in
#' each curve.
#' @param obs observed response.
#' @param S Smoothing matrix.
#' @param family a description of the error distribution and link function to
#' be used in the model. This can be a character string naming a family
#' function, a family function or the result of a call to a family function.
#' (See \code{\link{family}} for details of family functions.)
#' @param off off
#' @param offdf off, degrees of freedom
#' @param criteria The penalizing function. By default \emph{"Rice"} criteria.
#' Possible values are \emph{"GCV"}, \emph{"AIC"}, \emph{"FPE"},
#' \emph{"Shibata"}, \emph{"Rice"}.
#' @param W Matrix of weights.
#' @param trim The alpha of the trimming.
#' @param draw =TRUE, draw the curves, the sample median and trimmed mean.
#' @param \dots Further arguments passed to or from other methods.
#' 
#' @return Returns GCV score calculated for input parameters.  
#' 
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' 
#' @seealso See Also as \code{\link{GCV.S}}. \cr Alternative method:
#' \code{\link{CV.S}}
#' 
#' @references 
#' Wasserman, L. \emph{All of Nonparametric Statistics}. Springer
#' Texts in Statistics, 2006.
#' 
#' Hardle, W. \emph{Applied Nonparametric Regression}. Cambridge University
#' Press, 1994.
#' 
#' Febrero-Bande, M., Oviedo de la Fuente, M. (2012).  \emph{Statistical
#' Computing in Functional Data Analysis: The R Package fda.usc.} Journal of
#' Statistical Software, 51(4), 1-28. \url{http://www.jstatsoft.org/v51/i04/}
#' 
#' @keywords utilities
#' 
#' @examples
#' data(phoneme)
#' mlearn<-phoneme$learn
#' np<-ncol(mlearn)
#' tt<-mlearn[["argvals"]]
#' S1 <- S.NW(tt,2.5)
#' gcv1 <- dev.S(mlearn$data[1,],obs=(sample(150)), 
#' S1,off=rep(1,150),offdf=3)
#' gcv2 <- dev.S(mlearn$data[1,],obs=sort(sample(150)), 
#' S1,off=rep(1,150),offdf=3)
#'
#' @aliases dev.S
#' @export   
dev.S=function (y, S, obs,family = gaussian(),off,offdf,criteria="GCV",
	W = diag(1, ncol = ncol(S), nrow = nrow(S)), trim = 0, draw = FALSE, ...)
{
    n = ncol(S)
    mu = family$linkinv(S%*%y+off)
    sigma=sqrt(family$variance(mu))
    tab = list("GCV", "AIC", "FPE", "Shibata", "Rice")
    type.i = pmatch(criteria, tab)
    e = (obs - mu)/sigma
    if (trim > 0) {
            e.trunc = quantile(abs(e), probs = (1 - trim), na.rm = TRUE,
                type = 4)
            l <- which(abs(e) <= e.trunc)
    }
    else {  l = 1:n  }
    res = fdata.trace(t(e[l]) %*% W[l,l] %*% e[l])
    ndf <- sum(diag(S)[l],na.rm=TRUE)+offdf
    if (is.na(type.i)) {
        if (ndf>0.8*n)  vv = Inf
        else vv = 1/(1 - 2 * ndf/n)
    }
    else {
        vv<-switch(type.i,
                   "1"=if (ndf>0.8*n){vv=Inf} else {vv = (1 - ndf/n)^(-2)},
                   "2"=if (ndf>0.8*n){vv=Inf} else {vv = exp(2 * ndf/n)},
                   "3"=if (ndf>0.8*n){vv=Inf} else {vv = (1 + ndf/n)/(1 - ndf/n)},
                   "4"=if (ndf>0.8*n){vv=Inf} else {vv = (1 + ndf/n)/(1 - ndf/n)},
                   "5"=if (ndf>0.8*n){vv=Inf} else { vv = 1/(1 - 2 * ndf/n)})
         }
    return(res * vv/n)
}

