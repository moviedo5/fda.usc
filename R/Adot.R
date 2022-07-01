
#' @name PCvM.statistic
#' @aliases PCvM.statistic Adot
#' @title PCvM statistic for the Functional Linear Model with scalar response
#' @description Projected Cramer-von Mises statistic (PCvM) for the Functional Linear Model with scalar response (FLM): 
#' \eqn{Y=\big<X,\beta\big>+\varepsilon}{Y=<X,\beta>+\epsilon}.
#' 
#' @param  X Functional covariate for the FLM. 
#' The object must be either in the class \code{\link{fdata}} or in the class \code{\link{fd}}.
#' It is used to compute the matrix of inner products.
#' @param residuals Residuals of the estimated FLM.
#' @param p Number of elements of the functional basis where the functional covariate is represented.
#' @param Adot.vec Output from the \code{Adot} function (see Details). Computed if not given.
#' @param inpr Matrix of inner products of \code{X}. Computed if not given.
#' 
#' @details In order to optimize the computation of the statistic, the critical parts
#' of these two functions are coded in FORTRAN. The hardest part corresponds to the 
#' function \code{Adot}, which involves the computation of a symmetric matrix of dimension 
#' \eqn{n\times n}{n x n} where each entry is a sum of \eqn{n}{n} elements.
#' As this matrix is symmetric, the order of the method can be reduced from \eqn{O(n^3)}{O(n^3)}
#' to \eqn{O\big(\frac{n^3-n^2}{2}\big)}{O((n^3-n^2)/2)}. The memory requirement can also be reduced
#' to \eqn{O\big(\frac{n^2-n+2}{2}\big)}{O((n^2-n+2)/2)}. The value of \code{Adot} is a vector of 
#' length \eqn{\frac{n^2-n+2}{2}}{(n^2-n+2)/2} where the first element is the common diagonal 
#' element and the rest are the lower triangle entries of the matrix, sorted by rows (see Examples).
#'  
#' @return   For \code{PCvM.statistic}, the value of the statistic. For \code{Adot},
#' a suitable output to be used in the argument \code{Adot.vec}.
#' 
#' @references
#' Escanciano, J. C. (2006). A consistent diagnostic test for regression models using projections. 
#' Econometric  Theory, 22, 1030-1051. \doi{10.1017/S0266466606060506}
#'  
#' Garcia-Portugues, E., Gonzalez-Manteiga, W. and Febrero-Bande, M. (2014). A goodness--of--fit
#' test for the functional linear model with scalar response. Journal of Computational and 
#' Graphical Statistics, 23(3), 761-778. \doi{10.1080/10618600.2013.812519}
#'  
#' @note No NA's are allowed in the functional covariate.
#'  
#' @author Eduardo Garcia-Portugues. Please, report bugs and suggestions
#' to \if{latex}{\cr}\email{eduardo.garcia.portugues@@uc3m.es}
#'  
#' @seealso \code{\link{flm.test}}
#' @examples
#' # Functional process
#' X=rproc2fdata(n=10,t=seq(0,1,l=101))
#' # Adot
#' Adot.vec=Adot(X)
#' 
#' # Obtain the entire matrix Adot
#' Ad=diag(rep(Adot.vec[1],dim(X$data)[1]))
#' Ad[upper.tri(Ad,diag=FALSE)]=Adot.vec[-1]
#' Ad=t(Ad)
#' Ad=Ad+t(Ad)-diag(diag(Ad))
#' Ad
#' # Statistic
#' PCvM.statistic(X,residuals=rnorm(10),p=5)
#' @keywords  htest
#' 
# Adot function to call Fortran code
#' @rdname  PCvM.statistic
#' @export 
Adot=function(X,inpr){
  
  # Check if the Fortran function is correctly loaded
  if(!is.loaded("adot")) stop("Fortran function adot is not loaded")
  
  # Compute the inproduct
  if(missing(inpr)){
    if(is.fd(X)){
      inpr=t(X$coefs)%*%inprod(X$basis,X$basis)%*%X$coefs
    }else if(is.fdata(X)){
      inpr=inprod.fdata(X,X)
    }else{
      stop("X is not a fd or fdata object")
    }
  }
  
  # Number of functions in X
  n=dim(inpr)[1]
  
  # Check for repeated interior products (will cause error in the adot subroutine)
  if(anyDuplicated(diag(inpr))){
    diag(inpr)=diag(inpr)*(1+rexp(n,rate=1000))
  }
  
  # Call Fortran function
  res=.Fortran("adot",n=as.integer(n),inprod=t(inpr)[upper.tri(inpr,diag=TRUE)],Adot_vec=as.double(rep(0,(n*n-n+2)/2)))$Adot_vec
  
  # Return result
  return(res)
  
}

# PCvM statistic
#' @rdname  PCvM.statistic
#' @export 
PCvM.statistic=function(X,residuals,p,Adot.vec){
  
  # Check if the Fortran function is correctly loaded
  if(!is.loaded("pcvm_statistic")) stop("Fortran function pcvm_statistic is not loaded")
  
  # Obtain the sample size
  n=length(residuals)
  
  # Check if the lengths of e and X are the same
  if(is.fd(X)){
    if(n!=dim(X$coefs)[2]) stop("Incompatible lenghts of X.fd and residuals")
  }else if(is.fdata(X)){
    if(n!=dim(X$data)[1]) stop("Incompatible lenghts of X.fdata and residuals")
  }
  
  # Compute Adot.vec if missing
  if(missing(Adot.vec)) Adot.vec=Adot(X)
  
  # Call Fortran function
  res=n^(-2)*pi^((p/2)-1)/gamma(p/2)*.Fortran("pcvm_statistic",n=as.integer(n),Adot_vec=Adot.vec,residuals=residuals,statistic=0)$statistic
  
  return(res)
  
}

