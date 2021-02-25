#' The generalized correlated cross-validation (GCCV) score
#' 
#' Compute the  generalized correlated cross-validation (GCV) score. 
#' 
#' @details A.-If \code{trim=0}:\cr
#' \deqn{GCCV=\frac{\sum_{i=1}^n {y_{i}-\hat{y}_{i,b}}^2}{1-\frac{tr(C)}{n}^2}}{\sum(y-y.fit)^2 / (1-tr(C)/n)^2} 
#' where \eqn{S} is the smoothing matrix \eqn{S} and:\cr
#' A.-If \eqn{C=2S\Sigma - S\Sigma S} \cr
#' B.-If  \eqn{C=S\Sigma} \cr
#' C.-If  \eqn{C=S\Sigma S'} \cr
#' with \eqn{\Sigma} is the n x n covariance matrix with
#' \eqn{cor(\epsilon_i,\epsilon_j ) =\sigma}
#' Note: Provided that \eqn{C = I} and the smoother matrix S is symmetric and idempotent, as is the case for many linear fitting techniques, the trace term reduces to  \eqn{n - tr[S]}, which is proportional to the familiar denominator in GCV. 

#' @param y Matrix of set cases with dimension (\code{n} x \code{m}), where
#' \code{n} is the number of curves and \code{m} are the points observed in
#' each curve.
#' @param S Smoothing matrix, see \code{\link{S.NW}}, \code{\link{S.LLR}} or
#' @param criteria The penalizing function. By default \emph{"Rice"} criteria. 
#'  Possible values are \emph{"GCCV1"}, \emph{"GCCV2"}, \emph{"GCCV3"}, \emph{"GCV"}. 
#' @param W Matrix of weights.
#' @param trim The alpha of the trimming.
#' @param draw =TRUE, draw the curves, the sample median and trimmed mean.
#' @param metric Metric function, by default \code{\link{metric.lp}}.
#' @param \dots Further arguments passed to or from other methods.
#' 
#' @return { Returns GCV score calculated for input parameters.  }
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' @seealso See Also as \code{\link{optim.np}} \cr Alternative method:
#' \code{\link{CV.S}}
#' @references 
#' Wasserman, L. \emph{All of Nonparametric Statistics}. Springer Texts in Statistics, 2006. 
#' Hardle, W. \emph{Applied Nonparametric Regression}. Cambridge University Press, 1994.
#' Febrero-Bande,  M., Oviedo de la Fuente, M. (2012).  \emph{Statistical Computing in Functional Data Analysis: The R Package fda.usc.}
#' Journal of Statistical Software, 51(4), 1-28. \url{http://www.jstatsoft.org/v51/i04/}
#' @keywords utilities
#' @examples
#' \dontrun{
#' data(phoneme)
#' mlearn<-phoneme$learn
#' tt<-1:ncol(mlearn)
#' S1 <- S.NW(tt,2.5)
#' S2 <- S.LLR(tt,2.5)
#' gcv1 <- GCV.S(mlearn, S1)
#' gcv2 <- GCV.S(mlearn, S2)
#' gcv3 <- GCV.S(mlearn, S1,criteria="AIC")
#' gcv4 <- GCV.S(mlearn, S2,criteria="AIC")
#' gcv1; gcv2; gcv3; gcv4
#' }
#'  
#' @export
GCV.S=function(y,S,criteria="GCV",W=NULL,trim=0,draw=FALSE,metric=metric.lp,...){
    isfdata<-is.fdata(y)
    tab=list("GCV","AIC","FPE","Shibata","Rice")
    type.i=pmatch(criteria,tab)
    n=ncol(S);l=1:n
    if (isfdata)  {
         nn<-nrow(y)
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
         }
   else {
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
             res=fdata.trace(t(e[l])%*%W[l,l]%*%e[l])
             }
     else    res=fdata.trace(t(e)%*%W%*%e)
    }
    d<-diag(S)[l] 
    df<-sum(d)
    if (is.na(type.i))   {
                   if (mean(d,na.rm=TRUE)>0.5) vv=Inf
                   else   vv= 1/(1-2*mean(d,na.rm=TRUE))        }
    else {
        vv<-switch(type.i,
                   "1"=if (type.i==1)  vv=(1-mean(d,na.rm=TRUE))^(-2),
                   "2"=if (type.i==2)  vv= exp(2*mean(d)),
                   "3"=if (type.i==3)  vv=(1+mean(d,na.rm=TRUE))/(1-mean(d,na.rm=TRUE)),
                   "4"=if (type.i==4)  vv=(1+2*mean(d,na.rm=TRUE)),
                   "5"=if (type.i==5)  {
                           if (mean(d,na.rm=TRUE)>0.5) vv=Inf
                           else   vv= 1/(1-2*mean(d,na.rm=TRUE))   }
                  )
         }
out<-res*vv/n
 attr(out, "df") <- df           
return(out) }
  
