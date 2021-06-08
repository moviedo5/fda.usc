#' @name Outliers.fdata
#' @aliases outliers.depth.pond outliers.depth.trim outliers.thres.lrt outliers.lrt quantile.outliers.trim quantile.outliers.pond
#' @title outliers for functional dataset
#' @description Procedure for detecting funcitonal outliers.
#' 
#' @param fdataobj \code{\link{fdata}} class object.
#' @param nb The number of bootstrap samples.
#' @param smo The smoothing parameter for the bootstrap samples.
#' @param trim The alpha of the trimming.
#' @param quan Quantile to determine the cutoff from the Bootstrap procedure (by default=0.5)
#' @param dfunc Type of depth measure, by default \code{depth.mode}.
# @param ns Significance level, by defaul 1\%. 
#' @param \dots Further arguments passed to or from other methods.
#' @details
#' Outlier detection in functional data by likelihood ratio test (\code{outliers.lrt}). The threshold for outlier detection is given by the 
#'   \code{outliers.thres.lrt}.
#'  Outlier detection in functional data by depth measures: 
#'  \itemize{
#'  \item \code{outliers.depth.pond} function weights the data according to depth.
#'  \item \code{outliers.depth.trim} function uses trimmed data.
#'  }
#'  quantile.outliers.pond and quantile.outliers.trim functions provides the quantiles of the bootstrap samples for functional outlier detection by, respectively, weigthed and trimmed procedures. Bootstrap smoothing function (\code{\link{fdata.bootstrap}} with \code{nb} resamples) is applied to these weighted or trimmed data. If \code{smo=0} smoothed bootstrap is not performed.  The function returns a vector of size \code{1}x\code{nb} with bootstrap replicas of the quantile.
#' @return
#'  \itemize{
#'  \item \code{outliers}{ Indexes of functional outlier.}
#'  \item \code{dep.out}{ Depth value of functional outlier.}      
#'  \item \code{dep.out}{ Iteration in which the  functional outlier is detected.}        
#'  \item \code{quantile}{ Threshold for outlier detection.}  
#'  \item \code{dep}{ Depth value of functional data.}
#'}
#' @references
#' Cuevas A, Febrero M, Fraiman R. 2006.  \emph{On the use of bootstrap for estimating functions with functional data.} Computational Statistics and Data Analysis 51: 1063{-}1074.
#'   
#' Febrero-Bande, M., Galeano, P., and Gonzalez-Manteiga, W. (2008).  \emph{Outlier detection in functional data by depth measures with application to identify abnormal NOx levels}. Environmetrics 19, 4, 331{-}345. 
#' 
#' Febrero-Bande, M., Galeano, P. and Gonzalez-Manteiga, W.	 (2007). \emph{A functional analysis of NOx levels: location and scale estimation and outlier detection}. Computational Statistics 22, 3, 411{-}427.
#' 
#' Febrero-Bande,  M., Oviedo de la Fuente, M. (2012).  \emph{Statistical Computing in Functional Data Analysis: The R Package fda.usc.}
#' Journal of Statistical Software, 51(4), 1-28. \url{https://www.jstatsoft.org/v51/i04/}
#' 
#' @author
#'   Manuel Febrero-Bande, Manuel Oviedo de la Fuente \email{manuel.oviedo@@udc.es}
#' @seealso  See Also: \code{\link{fdata.bootstrap}}, \code{\link{Depth}}.
#' @examples
#' \dontrun{
#' data(aemet)
#' nb=20 # Time consuming
#' out.trim<-outliers.depth.trim(aemet$temp,dfunc=depth.FM,nb=nb)
#' plot(aemet$temp,col=1,lty=1)
#' lines(aemet$temp[out.trim[[1]]],col=2)
#' }
#' @keywords outliers 
#' 
#' @rdname  Outliers.fdata
#' @export
outliers.depth.pond<-function(fdataobj,nb=200,smo=0.05,quan=0.5,dfunc=depth.mode,...){
if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
 nas1<-is.na.fdata(fdataobj)
 if (any(nas1))  stop("fdataobj contain ",sum(nas1)," curves with some NA value \n")
 x<-fdataobj[["data"]]
tt<-fdataobj[["argvals"]]
rtt<-fdataobj[["rangeval"]]
n<-nrow(fdataobj)
m<-ncol(fdataobj)
if (is.null(n) && is.null(m)) stop("ERROR IN THE DATA DIMENSIONS")
if (is.null(row.names(fdataobj[["data"]]))) row.names(fdataobj[["data"]])=1:n
cutoff<-quantile(quantile.outliers.pond(x,dfunc=dfunc,nb=nb,smo=smo,...),probs=quan)
    hay<-1
    outliers<-dep.out<-ite<-c()
    curvasgood<-fdataobj
    row.names(curvasgood[["data"]])=1:n
    dd<-dfunc(curvasgood,...)     
    modal<-FALSE
        ii<-1
    if (!is.null(dd$dist)) {
      modal=TRUE 
      dd<-dfunc(curvasgood,...)
          }
    d<-dd$dep             
    while (hay==1){
#          d=dfunc(curvasgood,...)$dep
          if (is.null(outliers)){dtotal<-d}
          cutt<-d<cutoff
          elim<-which(cutt)
          fecha<-as.numeric(row.names(curvasgood[["data"]])[cutt])        
          if (length(elim)>0){
             dep.out<-c(dep.out,d[elim])
             curvasgood<-curvasgood[-elim,]
             outliers<-c(outliers,fecha)
          }
        if (length(elim)==0 || length(outliers)>n/5){hay<-0}
                else {
            if (modal) {
             mdist<-dd$dist[-elim,-elim]
            class(mdist)<-c("matrix","fdist")        
            dd<-dfunc(curvasgood,metric=mdist)
            }
          else dd<-dfunc(curvasgood,...)
          d<-dd$dep 
          }
          ite<-c(ite,rep(ii,length(elim)))
          ii<-ii+1
    }
outliers<-rownames(fdataobj[["data"]])[outliers]  
names(dep.out)<-NULL      
return(list("outliers"=outliers,"dep.out"=dep.out,"iteration"=ite,"quantile"=cutoff,"Dep"=dtotal))
}

#' @rdname Outliers.fdata
#' @export
outliers.depth.trim<-function(fdataobj,nb=200,smo=0.05,trim=0.01,quan=0.5,
                              dfunc=depth.mode,...){
  if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
  nas1<-is.na.fdata(fdataobj)
  if (any(nas1))  stop("fdataobj contain ",sum(nas1)," curves with some NA value \n")         
  x<-fdataobj[["data"]]
  tt<-fdataobj[["argvals"]]
  rtt<-fdataobj[["rangeval"]]
  n<-nrow(fdataobj)
  m<-ncol(fdataobj)
  # print(dfunc=depth.mode)
  if (is.null(n) && is.null(m)) stop("ERROR IN THE DATA DIMENSIONS")
  if (is.null(row.names(fdataobj[["data"]]))) row.names(fdataobj[["data"]])=1:n
  cutoff<-quantile(quantile.outliers.trim(fdataobj,dfunc=dfunc,nb=nb,smo=smo,
                                          trim=trim,...),probs=quan)
  hay<-1
  outliers<-dep.out<-ite<-c()
  ii<-1
  curvasgood<-fdataobj
  dd<-dfunc(curvasgood,trim=trim,...)     
  modal<-FALSE
  if (!is.null(dd$dist)) {
    modal=TRUE   
    dd<-dfunc(curvasgood,trim=trim,...)
  }
  d<-dd$dep          
  rwn=names(d)=rownames(curvasgood[["data"]])=1:n
  while (hay==1){              
    if (is.null(outliers)){dtotal<-d}
    cutt<-d<cutoff
    fecha<-as.numeric(rownames(curvasgood[["data"]])[cutt])            
    elim<-which(cutt)
    if (length(elim)>0){
      dep.out<-c(dep.out,d[elim])
      curvasgood<-curvasgood[-elim,]
      outliers<-c(outliers,fecha)     
    }    
    if (length(elim)==0 || length(outliers)>n/5){hay<-0}
    else {
      if (modal) {
        mdist<-dd$dist[-elim,-elim]
        class(mdist) <- c("matrix","fdist")        
        dd<-dfunc(curvasgood,trim=trim,metric=mdist,scale=FALSE)
      }
      else dd<-dfunc(curvasgood,trim=trim,...)
      d<-dd$dep 
    }
    ite<-c(ite,rep(ii,length(elim)))
    ii<-ii+1
  }
  outliers<-rownames(fdataobj[["data"]])[outliers]    
  names(dep.out)<-NULL    
  return(list("outliers"=outliers,"dep.out"=dep.out,"iteration"=ite,"quantile"=cutoff,"Dep"=dtotal))
}

#' @rdname Outliers.fdata
#' @export
outliers.lrt<-function(fdataobj,nb=200,smo=0.05,trim=0.10,...){
  if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
  x<-fdataobj[["data"]]
  tt<-fdataobj[["argvals"]]
  rtt<-fdataobj[["rangeval"]]
  n<-nrow(fdataobj)
  rowname.ori<-row.names(fdataobj[["data"]])
  rowname.num<-1:n
  row.names(fdataobj[["data"]])<-rowname.num    
  m<-ncol(fdataobj)
  if (is.null(n) && is.null(m)) stop("Error in the data dimensions")
  out.thres.lrt<-outliers.thres.lrt(fdataobj,nb=nb,smo=smo,trim=trim,...)
  hay<-1
  outliers<-c()
  valor.estadistico<-c()
  nout<-0
  ngood<-n-nout
  curvasgood<-fdataobj
  i<-1
  maxiter<-5
  while (hay==1 && i<maxiter){
    i<-i+1
    aux<-c()
    auxmean<-func.trim.mode(curvasgood,trim=trim,...)#[["data"]][1,]
    aa<-func.trimvar.mode(curvasgood,trim=trim,...)
    auxdt<-sqrt(aa)
    divi<-auxmean/sqrt(aa)    
    for (j in 1:ngood){
      #cat(j);          print(metric.lp(curvasgood[j]/auxdt,divi,...))
      aux[j]<-metric.lp(curvasgood[j]/auxdt,divi,...)[1,1]
    }
    maximo<-as.numeric(max(aux))
    fecha<-rowname.num[maximo==aux] 
    elim<-which(maximo==aux)
    if (maximo>out.thres.lrt){
      curvasgood<-curvasgood[-elim]
      outliers<-c(outliers,rowname.ori[fecha])
      valor.estadistico<-c(valor.estadistico,maximo)
      nout<-nout+1
      ngood<-n-nout
    }
    if (maximo<out.thres.lrt){hay<-0}
    
  }
  if (i==maxiter) warning("No convergence achieved ")
  return(list("outliers"=outliers,"stat.value"=valor.estadistico,
              "percentile"=out.thres.lrt))
  c(outliers,valor.estadistico,out.thres.lrt)
}

#' @rdname Outliers.fdata
#' @export
outliers.thres.lrt<-function(fdataobj,nb=200,smo=0.05,trim=0.10,...){
  if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
  x<-fdataobj[["data"]]
  tt<-fdataobj[["argvals"]]
  names<-fdataobj[["names"]]
  rtt<-fdataobj[["rangeval"]]
  n<-nrow(fdataobj)
  m<-ncol(fdataobj)
  if (is.null(n) && is.null(m)) stop("ERROR IN THE DATA DIMENSIONS")
  maximos<-c()
  aux<-c()
  remuestras.boot<-matrix(NA,nrow=n,ncol=m)
  pb=txtProgressBar(min=0,max=nb,style=3)
  for (i in 1:nb){
    setTxtProgressBar(pb,i-0.5)
    bmuestra<-fdataobj[sample(1:n,size=n,replace=TRUE),]
    if (smo>0) {bmuestra[["data"]]<-bmuestra[["data"]]+mvrnorm(n=n,rep(0,m),var(fdataobj[["data"]])*smo)}
    auxmean<-func.trim.mode(bmuestra,trim=trim,...)
    auxdt<-sqrt(as.vector(func.trimvar.mode(bmuestra,trim=trim,...)[["data"]]))
    d<-matrix(NA,nrow=n,ncol=m)
    for (j in 1:m){d[,j]<-1-abs(.5-rank(bmuestra[,j][["data"]],ties.method="average")/n)}
    #        ans<-apply(d,1,sum)
    ans<-rowSums(d)
    rid<-rank(ans,ties.method="first")
    bmuestra.trim<-bmuestra[rid>=floor(trim*n),]
    for (j in 1:(n-floor(trim*n)))
    {aux[j]<-metric.lp(bmuestra.trim[j,][["data"]]/auxdt,auxmean[["data"]]/auxdt,...)}
    maximos[i]<-as.numeric(max(aux))
    setTxtProgressBar(pb,i)
  }
  close(pb)
  ans<-as.numeric(quantile(maximos,.99))
  ans
}


