#' Calculation of the smoothing parameter (h) for a functional data
#' 
#' Calculation of the smoothing parameter (h) for a functional data using
#' nonparametric kernel estimation.
#' 
#' @param fdataobj \code{\link{fdata}} class object.
#' @param prob Vector of probabilities for extracting the quantiles of the distance matrix. If \code{length(prob)=2} 
#' a sequence between \code{prob[1]} and \code{prob[2]} of length \code{len}.
#' @param len Vector length of smoothing parameter \code{h} to return only used when \code{length(prob)=2}.
#' @param metric If is a function: name of the function to calculate the
#' distance matrix between the curves, by default \code{\link{metric.lp}}.  If
#' is a matrix: distance matrix between the curves.
# @param Ker Type of asymmetric kernel used, by default asymmetric normal
#' kernel.
#' @param type.S Type of smothing matrix \code{S}.  Possible values are:
#' Nadaraya-Watson estimator \emph{"S.NW"} and K nearest neighbors estimator
#' \emph{"S.KNN"}
#' @param Ker Kernel function. By default, \emph{Ker.norm}. Useful for scaling the bandwidth values
#' according to Kernel
#' @param \dots Arguments to be passed for metric argument.
#' @return Returns the vector of smoothing parameter or bandwidth \code{h}.
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' @seealso See Also as \code{\link{metric.lp}}, \code{\link{Kernel}} and
#' \code{\link{S.NW}}. \cr Function used in \code{\link{fregre.np}} and
#' \code{\link{fregre.np.cv}} function.
#' @keywords nonparametric
#' @examples
#' \dontrun{
#' data(aemet)
#' h1<-h.default(aemet$temp,prob=c(0.025, 0.25),len=2)
#' mdist<-metric.lp(aemet$temp)
#' h2<-h.default(aemet$temp,len=2,metric=mdist)
#' h3<-h.default(aemet$temp,len=2,metric=semimetric.pca,q=2)
#' h4<-h.default(aemet$temp,len=2,metric=semimetric.pca,q=4)
#' h5<-h.default(aemet$temp,prob=c(.2),type.S="S.KNN")
#' h1;h2;h3;h4;h5
#' }
#' @export
h.default=function (fdataobj, prob=c(0.025,0.25),len=51,
                    metric = metric.lp,type.S ="S.NW",Ker=Ker.norm,...)
{
    #type.S<-deparse(substitute(type.S))
 # print(type.S)
    if (length(prob)==1) {
      prob<-c(prob,prob)
      len=1
    }
    else {if (len==1) len=2}
    if  (prob[1]==prob[2]) len=1
    
    if (!is.fdata(fdataobj))   fdataobj = fdata(fdataobj)
    n=nrow(fdataobj)
    
    if (is.character(Ker)){
      nker <- function(u,mik=Ker){get(mik)(u)}
    } else {
      nker <- function(u,mik=Ker){mik(u)} 
    }
    
    ck=integrate(function(u){nker(u)^2},-4,4)$value
    dk=integrate(function(u){u^2*nker(u)},-4,4)$value
    c0=(ck/dk^2)^(1/5)/1.3510   #1.3510 from unif      
      if (type.S %in% c("S.NW","S.LLR","S.LPR","S.LCR")) {
        if (is.matrix(metric)) {mdist=metric}
        else {mdist=metric(fdataobj,fdataobj,...)}
        diag(mdist) = Inf
        class(mdist)<-"matrix"
		  if (length(prob)>2) {
		  qua=quantile(mdist,probs=prob,na.rm=TRUE,type=4)
		  h=sort(unique(qua))*c0
			} else {
        q.min = quantile(mdist, probs = min(prob), na.rm = TRUE,type=4)
        q.max = quantile(mdist, probs = max(prob), na.rm = TRUE,type=4)
        h.min = q.min
        h.max = max(h.min,q.max)
        h = unique(seq(h.min, h.max, len = len))*c0
		  }
        #h<-h+h*0.0001
        #if (any(h==0)) h<-h[(which(h==0)+1):len]
        h <- h[h>0]
        }
        else if (type.S=="S.KNN") {
         # print("S.knnn")
           if (len>n) {len=n-2;print("len=nrow-2")}
           h.min = min(1,max(1,floor(min(prob)*n)))
           h.max = max(h.min,floor(max(prob)*n))
           if (len==1)           h = h.max
           else  h = seq_len(h.max)# h.min:h.max
         # print(h)
        }        #else {        #           stop("Error in type.S argument                }
        return(h)
 }

