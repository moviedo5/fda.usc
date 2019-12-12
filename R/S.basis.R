#' Smoothing matrix with roughness penalties by basis representation.
#' 
#' @description Provides the smoothing matrix \code{S} with roughness penalties.
#' 
#' @details Provides the smoothing matrix S for the discretization points \code{tt} and
#' b\code{basis} with roughness penalties. If \code{lambda=0} is not used
#' penalty, else a basis roughness penalty matrix is caluclated using
#' \link[fda]{getbasispenalty}.
#' 
#' @param tt Discretization points.
#' @param basis Basis to use. See \link[fda]{create.basis}.
#' @param lambda A roughness penalty. By default, no penalty \code{lambda}=0.
#' @param Lfdobj See \link[fda]{eval.penalty}.
#' @param w Optional case weights.
#' @param \dots Further arguments passed to or from other methods. Arguments to
#' be passed by default to \link[fda]{create.basis}
#' @return Return the smoothing matrix \code{S}.
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@usc.es}
#' @seealso See Also as \code{\link{S.np}}
#' @references
#' 
#' Ramsay, James O. and Silverman, Bernard W. (2006). \emph{ Functional Data
#' Analysis}, 2nd ed., Springer, New York.
#' 
#' Wasserman, L. \emph{All of Nonparametric Statistics}. Springer Texts in
#' Statistics, 2006.
#' @keywords smooth
#' @examples
#' \dontrun{
#' np=101
#' tt=seq(0,1,len=np)
#' 
#' nbasis=11
#' base1 <- create.bspline.basis(c(0, np), nbasis)
#' base2 <- create.fourier.basis(c(0, np), nbasis)
#' 
#' S1<-S.basis(tt,basis=base1,lambda=3)
#' image(S1) 
#' S2<-S.basis(tt,basis=base2,lambda=3)
#' image(S2)
#' }
#' @export  S.basis
S.basis=function(tt,basis,lambda=0,Lfdobj=vec2Lfd(c(0,0)),w=NULL,...){
    phi=getbasismatrix(tt,basis,...)
    np<-length(tt)
    if (is.null(w)) w<-rep(1,np)
    if (!is.matrix(w)) w<-diag(w)
    if (lambda!=0) {
       R=getbasispenalty(basis,Lfdobj)
       S=phi%*%solve(t(phi)%*%w%*%phi+lambda*R)%*%t(phi)%*%w}
   else {S=phi%*%solve(t(phi)%*%w%*%phi)%*%t(phi)%*%w}
    return(S)
}
