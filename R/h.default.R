#' Calculation of the smoothing parameter (h) for a functional data
#' 
#' Calculation of the smoothing parameter (h) for a functional data using
#' nonparametric kernel estimation.
#' 
#' @param fdataobj \code{\link{fdata}} class object.
#' @param prob Range of probabilities for the quantiles of the distance matrix.
#' @param len Vector length of smoothing parameter \code{h} to return.
#' @param metric If is a function: name of the function to calculate the
#' distance matrix between the curves, by default \code{\link{metric.lp}}.  If
#' is a matrix: distance matrix between the curves.
# @param Ker Type of asymmetric kernel used, by default asymmetric normal
#' kernel.
#' @param type.S Type of smothing matrix \code{S}.  Possible values are:
#' Nadaraya-Watson estimator \emph{"S.NW"} and K nearest neighbors estimator
#' \emph{"S.KNN"}
#' @param \dots Arguments to be passed for metric argument.
#' @return Returns the vector of smoothing parameter or bandwidth \code{h}.
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@usc.es}
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
                    metric = metric.lp,type.S ="S.NW",...)
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
  #  if (is.function(Ker))    if (body(Ker)==body(AKer.norm))  Ker= "AKer.norm"
  # S.LLR<-function (tt, h, Ker = Ker.norm,w=NULL,cv=FALSE)
  # S.LPR<-function (tt, h, p=1, Ker = Ker.norm, w = NULL, cv = FALSE) 
  # S.LCR<-function(tt, h, Ker=Ker.norm, w=NULL, cv=FALSE
  # S.NW<-function (tt, h=NULL, Ker = Ker.norm,w=NULL,cv=FALSE)
  # S.KNN<-function(tt,h=NULL,Ker=Ker.unif,w=NULL,cv=FALSE)
      
      if (type.S %in% c("S.NW","S.LLR","S.LPR","S.LCR")) {
        if (is.matrix(metric)) {mdist=metric}
        else {mdist=metric(fdataobj,fdataobj,...)}
        diag(mdist) = Inf
        #h0 <- apply(mdist, 1, min, na.rm = TRUE)
        #h.max = max(h0)
        #h.med = median(h0)
        class(mdist)<-"matrix"
        q.min = quantile(mdist, probs = prob[1], na.rm = TRUE,type=4)
        q.max = quantile(mdist, probs = prob[2], na.rm = TRUE,type=4)
        #h.min = max(c(drop(q.min),h.max))#h.med
        #h.max = max(c(drop(q.max), h.max))
        #if (Ker== "AKer.norm" ) {
        #        h.max = min(q.max, h.max)
        #        h.min = min(q.min, h.med)
        #}
        h.min = q.min
        h.max = max(h.min,q.max)
        h = unique(seq(h.min, h.max, len = len))
        #h<-h+h*0.0001
        #if (any(h==0)) h<-h[(which(h==0)+1):len]
        h <- h[h>0]
        }
        else if (type.S=="S.KNN") {
         # print("S.knnn")
           if (len>n) {len=n-2;print("len=nrow-2")}
           h.min = min(1,max(1,floor(prob[1]*n)))
           h.max = max(h.min,floor(prob[2]*n))
           if (len==1)           h = h.max
           else  h = seq_len(h.max)# h.min:h.max
         # print(h)
        }        #else {        #           stop("Error in type.S argument                }
        return(h)
 }

