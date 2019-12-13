#' @rdname rwild
#' @aliases rwild
#' @title Wild bootstrap residuals
#' @description   The wild bootstrap residuals are computed as \eqn{residuals*V}{residualsV}, where \eqn{V} is a sampling from a random variable (see details section). 

#' @param residuals residuals
#' @param type Type of distribution of V.

#' @details  For the construction of wild bootstrap residuals, sampling from a random variable \eqn{V}{V} such that \eqn{E[V^2]=0}{E[V^2]=0} and \eqn{E[V]=0}{E[V]=0} is needed. 
#'  A simple and suitable \eqn{V}{V} is obtained with a discrete variable of the form:
#'  \itemize{
#'  \item ``golden'', Sampling from golden section bootstrap values suggested by Mammen (1993).
#'  \deqn{P\Bigg\{ V=\frac{1-\sqrt{5}}{2} \Bigg\} = \frac{5+\sqrt{5}}{10} \, and \, P\Bigg\{ V=\frac{1+\sqrt{5}}{2} \Bigg\} = \frac{5-\sqrt{5}}{10},}{P\{ V=(1-\sqrt 5)/2 \} = (5+\sqrt 5)/10 and P\{ V=(1+\sqrt 5)/2 \} = (5-\sqrt 5)/10,}
#'  which leads to the \emph{golden section bootstrap}.  
#'  \item ``Rademacher'', Sampling from Rademacher distribution values
#'  \eqn{\big\{-1,\,1\big\}}{-1,1} with probabilities \eqn{\big\{\frac{1}{2},\,\frac{1}{2}\big\}}{\{1/2, 1/2\}}, respectively.  
#'  \item ``normal'', Sampling from a standard normal distribution.
#'  }

#' @return 
#'   The wild bootstrap residuals computed using a sample of the random variable \eqn{V}{V}.
#' @author 
#'   Eduardo Garcia-Portugues, Manuel Febrero-Bande and  Manuel Oviedo de la Fuente \email{manuel.oviedo@@usc.es}.  
#' @seealso 
#'  \code{\link{flm.test}}, \code{\link{flm.Ftest}}, \code{\link{dfv.test}}, \code{\link{fregre.bootstrap}}
#' @references 
#'   Mammen, E. (1993). \emph{Bootstrap and wild bootstrap for high dimensional linear models}.
#'     Annals of Statistics 21, 255 285.
#'       Davidson, R. and E. Flachaire (2001). \emph{The wild bootstrap, tamed at last}. working paper IER1000, Queens University.
#' @examples 
#' n<-100
#' # For golden wild bootstrap variable
#' e.boot0=rwild(rep(1,len=n),"golden")
#' # Construction of wild bootstrap residuals
#' e=rnorm(n)
#' e.boot1=rwild(e,"golden")
#' e.boot2=rwild(e,"Rademacher")
#' e.boot3=rwild(e,"normal")
#' summary(e.boot1)
#' summary(e.boot2)
#' summary(e.boot3)
#'              
#' @keywords {distribution}    

##############################################
##          Functional Linear Model         ##
##             Goodness-of-fit              ##
##############################################

##############################################
## File created by Eduardo Garcia-Portugues ##
## using code from library fda.usc          ##
##############################################

# Generation of centred random variables with unit variance for the wild bootstrap
#' @export
rwild=function(residuals,type="golden"){
  n=length(residuals)
  res=switch(type,
             golden=sample(c((1-sqrt(5))/2,(1+sqrt(5))/2),size=n,prob=c((5+sqrt(5))/10,(5-sqrt(5))/10),replace=TRUE), 
             Rademacher={	
               sample(c(-1,1),size=n,prob=c(.5,.5),replace=TRUE)
             },
             normal=rnorm(n)
  )
  return(residuals*res)
}

