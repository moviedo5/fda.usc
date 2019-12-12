#' The generalized correlated cross-validation (GCCV) score.
#' 
#' The generalized correlated cross-validation (GCV) score.
#' 
#' @details \deqn{ }{\sum(y-y.fit)^2 / (1-tr(C)/n)^2}\deqn{GCCV=\frac{\sum_{i=1}^n
#' {y_{i}-\hat{y}_{i,b}}^2}{1-\frac{tr(C)}{n}^2} }{\sum(y-y.fit)^2 /
#' (1-tr(C)/n)^2} %where \eqn{C=2{S\Sigma(\theta)}-{S\Sigma(\theta)S'}} \cr
#' %and \eqn{\Sigma=\sigma C} is the n x n covariance matrix with
#' \eqn{cor(\epsilon_i,\epsilon_j ) =\sigma}\cr
#' 
#' where \eqn{S} is the smoothing matrix \eqn{S} and:\cr A.-If \eqn{C=2S\Sigma
#' - S\Sigma S} \cr B.-If \eqn{C=S\Sigma} \cr C.-If \eqn{C=S\Sigma S'} \cr with
#' \eqn{\Sigma} is the n x n covariance matrix with
#' \eqn{cor(\epsilon_i,\epsilon_j ) =\sigma}
#' 
#' @note Provided that \eqn{C = I} and the smoother matrix S is symmetric and
#' idempotent, as is the case for many linear fitting techniques, the trace
#' term reduces to \eqn{n - tr[S]}, which is proportional to the familiar
#' denominator in GCV.
#' 
#' @param y Response vectorith length \code{n} or Matrix of set cases with
#' dimension (\code{n} x \code{m}), where \code{n} is the number of curves and
#' \code{m} are the points observed in each curve.
#' @param S Smoothing matrix, see \code{\link{S.NW}}, \code{\link{S.LLR}} or
#' \eqn{S.KNN}.
#' @param criteria The penalizing function. By default \emph{"Rice"} criteria.
#' "GCCV1","GCCV2","GCCV3","GCV") Possible values are \emph{"GCCV1"},
#' \emph{"GCCV2"}, \emph{"GCCV3"}, \emph{"GCV"}.
#' @param W Matrix of weights.
#' @param trim The alpha of the trimming.
#' @param draw =TRUE, draw the curves, the sample median and trimmed mean.
#' @param metric Metric function, by default \code{\link{metric.lp}}.
#' @param \dots Further arguments passed to or from other methods.
#' @return  Returns GCCV score calculated for input parameters.  
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@usc.es}
#' @seealso See Also as \code{\link{optim.np}}. \cr Alternative method
#' (independent case): \code{\link{GCV.S}}
#' @references Carmack, P. S., Spence, J. S., and Schucany, W. R. (2012).
#' Generalised correlated cross-validation. Journal of Nonparametric
#' Statistics, 24(2):269--282.
#' 
#' Oviedo de la Fuente, M., Febrero-Bande, M., Pilar Munoz, and Dominguez, A.
#' Predicting seasonal influenza transmission using Functional Regression
#' Models with Temporal Dependence. arXiv:1610.08718.
#' \url{https://arxiv.org/abs/1610.08718}
#' @keywords utilities
#' @examples 
#' \dontrun{
#' data(tecator)
#' x=tecator$absorp.fdata
#' x.d2<-fdata.deriv(x,nderiv=)
#' tt<-x[["argvals"]]
#' dataf=as.data.frame(tecator$y)
#' y=tecator$y$Fat
#' # plot the response
#' plot(ts(tecator$y$Fat))
#' 
#' nbasis.x=11;nbasis.b=7
#' basis1=create.bspline.basis(rangeval=range(tt),nbasis=nbasis.x)
#' basis2=create.bspline.basis(rangeval=range(tt),nbasis=nbasis.b)
#' basis.x=list("x.d2"=basis1)
#' basis.b=list("x.d2"=basis2)
#' ldata=list("df"=dataf,"x.d2"=x.d2)
#' # No correlation
#' res.gls=fregre.gls(Fat~x.d2,data=ldata, 
#'                    basis.x=basis.x,basis.b=basis.b)
#' # AR1 correlation                   
#' res.gls=fregre.gls(Fat~x.d2,data=ldata, correlation=corAR1(),
#'                    basis.x=basis.x,basis.b=basis.b)
#' GCCV.S(y,res.gls$H,"GCCV1",W=res.gls$W)
#' res.gls$gcv
#' }
#' @export GCCV.S
GCCV.S=function(y,S,criteria="GCCV1",W=NULL,trim=0,draw=FALSE,metric=metric.lp,...){
    isfdata<-is.fdata(y)
#    tab=list("GCV","AIC","FPE","Shibata","Rice")
    tab=list("GCCV1","GCCV2","GCCV3","GCV")
    type.i=pmatch(criteria,tab)
    n=ncol(S)
    l=1:n
    if (isfdata) #to smooth function (optim.np,optim.basis)
    { 
         nn<-ncol(y)
         if (is.null(W)) W<-diag(nn)
         y2=t(y$data)
         y.est=t(S%*%y2)
         y.est<-fdata(y.est,y$argvals, y$rangeval, y$names)
         e <- y - y.est
#         e$data<-sqrt(W)%*%(e$data)   
         ee <- drop(norm.fdata(e,metric=metric,...)[,1]^2)
         if (trim>0) {
            e.trunc=quantile(ee,probs=(1-trim),na.rm=TRUE,type=4)
            ind<-ee<=e.trunc
            if (draw)  plot(y,col=(2-ind))
            l<-which(abs(ee)<=e.trunc)
            res = mean(ee[ind],na.rm=TRUE)
            }
        else  res = mean(ee, na.rm = TRUE)
        ee<-e$data
   }
   else   #to regression function
   {
        if (is.null(W)) W<-diag(n)         
        if (is.matrix(y)&&(ncol(y)==1) ){y2<-y;draw<-FALSE}
        else if (is.vector(y)){y2<-y;draw<-FALSE}
        else stop("y is not a fdata,  vector or matrix")
     y.est=S%*%y2
     e=y2-y.est
     if (trim>0) {
             ee = t(e)
             e.trunc=quantile(abs(ee),probs=(1-trim),na.rm=TRUE,type=4)
             l<-which(abs(ee)<=e.trunc)
             e<-e[l]
             }
             ee<-e
    }
    d<-diag(S)[l] 
    vv<-switch(type.i,
                   "1"={
                   Sigma<-solve(W)
                   SC<-S%*%Sigma
                   df<-trace.matrix(2*SC-SC%*%t(S))
                   mean(ee^2)/(1-df/n)^2                                                    
                   },
                   "2"={
                   Sigma<-solve(W)
                   df<-trace.matrix(S%*%Sigma)
                   mean(ee^2)/(1-df/n)^2
                   },                                                     
                   "3"={
                   Sigma<-solve(W)
                   df<-trace.matrix(S%*%Sigma%*%t(S))
                   mean(ee^2)/(1-df/n)^2
                   }
                   ,"4"={
#                   Sigma<-solve(W)
                   df<-trace.matrix(S)
                   mean(ee^2)/(1-df/n)^2
                   })
 attr(vv, "df") <- df    
 return(vv)
 }
  

################################################################################
################################################################################
GCCV.S.version.simple=function(y,S,criteria="GCV",W=diag(ncol(S)),Sigma=W,trim = 0, draw = FALSE, metric = metric.lp, 
                ...) {
  # print("entra GCV.GLSS")
  #    tab=list("GCV","AIC","FPE","Shibata","Rice","GCCV1","GCCV2","GCCV3")
  #    type.i=pmatch(criteria,tab)
  # print(type.i)    
  n=ncol(S);l=1:n
  y.est=(S%*%y)
  e <- y - y.est
  # print("entra 2")  
  vv<-switch(criteria,                  
             "GCV"=n*trace.matrix(t(e) %*% W %*% e)/(n-trace.matrix(S))^2,
             "GCCV1"={
               SC<-S%*%Sigma
               mean(e^2)/(1-trace.matrix(2*SC-SC%*%t(S))/n)^2                                                    
             },
             "GCCV2"=mean(e^2)/(1-trace.matrix(S%*%Sigma)/n)^2,                                                     
             "GCCV3"=mean(e^2)/(1-trace.matrix(S%*%Sigma%*%t(S))/n)^2                                                     
  )
  return(vv)
}
################################################################################
