euclidean<-function(u,w=rep(1,length(u))) {  
  u<- sqrt(sum(w*u^2,na.rm=TRUE))
  st<-attributes(u)
  attributes(u)<-st
  attr(u, "weigthts") <- w
  u
}
manhattan<-function(u,w=rep(1,length(u))) {  
  u<- sqrt(sum(w*abs(u),na.rm=TRUE))
  st<-attributes(u)
  attributes(u)<-st
  attr(u, "weigthts") <- w
  u
}
minkowski<-function(u,w=rep(1,length(u)),p) {  
  u<- (sum(w*(u^p),na.rm=TRUE))^(1/p)
  st<-attributes(u)
  attributes(u)<-st
  attr(u, "weigthts") <- w
  u
}
maximum<-function(u,w=rep(1,length(u))) {  
  u<- sqrt(w*max(u,na.rm=TRUE))
  st<-attributes(u)
  attributes(u)<-st
  attr(u, "weigthts") <- w
  u
}
###########################################################
#euclidean<-function (u)   sqrt(sum(u^2, na.rm = TRUE))
#manhattan<-function(u) sqrt(sum(abs(u),na.rm=TRUE))
#minkowski<-function(u,p) (sum(u,na.rm=TRUE)^p)^1/p
#no hace falta w
#maximum<-function(u) sqrt(max(u,na.rm=TRUE))
################################################################################
# cambiar en fregre.gsam.vs
# dist.list<-function (ldata, 
#                 #      parallel=par.fda.usc$parallel,
#                       ...) {
#   lenldata <- length(ldata)
#   ldist <- list()
#   for (i in 1:lenldata) {
#         if (is.factor(ldata[[i]])) 
#           ldata[[i]] <- model.matrix(~ldata[[i]])
#         if (is.fdata(ldata[[i]])) 
#           ldist[[i]] <- metric.lp(ldata[[i]], ...)
#         else ldist[[i]] <- as.matrix(dist(ldata[[i]]), diag = TRUE, 
#                                      upper = TRUE, p = 2)
#         ldist[[i]]
#       }
#   names(ldist) <- names(ldata)
#   ldist
# }
##################################
dist.list <- function (ldata,...) {
     lenldata <- length(ldata)
    ldist <- list()
      #ldist <- foreach(i=  1:lenldata,.combine="c") %dopar% {
        
        for (i in 1:lenldata) {
          if (is.factor(ldata[[i]])) 
            ldata[[i]] <- model.matrix(~ldata[[i]])
          if (is.fdata(ldata[[i]])) 
            ldist[[i]] <- metric.lp(ldata[[i]], ...)
          else ldist[[i]] <- as.matrix(dist(ldata[[i]]), diag = TRUE, 
                                       upper = TRUE, p = 2)
        }
    names(ldist) <- names(ldata)
    ldist
}
################################# 
# metric.ldata2=function(x,y=x
#                       ,metric
#                       ,par.metric=list()
#                       ,w
#                       ,method="euclidean") {
#   if (any(class(x)=="list")){
#     lenl<-length(x)
#     lenl2<-length(y)
#     n<-nrow(x[[1]])
#     m<-nrow(y[[1]]) 
#     mdist2<-matrix(0,n,n)
#     amdist<-array(NA,dim=c(n,m,lenl))
#     ldist<-mdist<-list()
#     nam1<-names(x)
#     nam2<-names(y) 
#     attr<-list()
#     if (is.null(nam1)) {names(x)<-nam1<-paste("var",1:lenl,sep="")}
#     if (is.null(nam2)) {names(y)<-nam2<-paste("var",1:lenl2,sep="")} 
#     # verificar que todos los nam1 estan en nam2
#     # verificar que longitud de metric y par.metric es la misma q x
#     if (missing(metric)){
#       if (is.fdata(x[[nam1[1]]])) metric<-rep("metric.lp2",len=lenl)
#       else    metric<-rep("metric.dist",len=lenl)   
#     }
#     if (missing(w)) w<-rep(1,length=lenl)
#    
#     for (i in 1:lenl){
# #      print("metric.x3")
#       if (is.character(metric)) {
#         if (is.null(par.metric[[nam1[i]]])) par.metric[[nam1[i]]]<-list()
#         if (is.fdata(x[[nam1[i]]])) { 
#           par.metric[[nam1[i]]][["fdata1"]]<-x[[nam1[i]]]
#           par.metric[[nam1[i]]][["fdata2"]]<-y[[nam1[i]]]   
# #print(metric[i])          
#           mdist<-do.call(metric[i],par.metric[[nam1[i]]])
#         }
#         else {
#           par.metric[[nam1[i]]][["x"]]<-x[[nam1[i]]]
#           par.metric[[nam1[i]]][["y"]]<-y[[nam1[i]]]   
#           mdist<-do.call(metric[i],par.metric[[nam1[i]]])
#         }
#       }
#       else {
#         if (is.list(metric)) mdist<-metric[[nam1[i]]]
#       }
#       ldist[[nam1[i]]]<-mdist
#       amdist[,,i]<- mdist
#       #attr[[nam1[i]]]<-attributes(mdist)     
#     }
#     # print("sale1")
#     #print(Weights)
#     ###for (k in len) sqrt(w[i])*
#     # print(w)
#     mdist<-apply(amdist,1:2,method,w=w)
#     # print("sale2")
#   }
#   else stop("Error in x argument")
#   # print("a")
#   # print(dim(amdist))
#   # print(dim(mdist))
#   # print(attributes(mdist))
#   # print(attr)
#   #attributes(mdist)<-atr
#   attr(mdist, "method") <-method
#   attr(mdist, "w") <- w
#   for (i in 1:lenl) attr(mdist, nam1[i]) <- attributes(ldist[[nam1[i]]])
#   return(mdist)
# }
###############################
# metric.ldata1=function(x,y=x
#                        ,metric
#                        ,par.metric=list()
#                        ,w
#                        ,method="euclidean") {
#   if (any(class(x)=="list")){
#     lenl<-length(x)
#     lenl2<-length(y)
#     n<-nrow(x[[1]])
#     m<-nrow(y[[1]]) 
#     mdist2<-matrix(0,n,n)
#     amdist<-array(NA,dim=c(n,m,lenl))
#     ldist<-mdist<-list()
#     nam1<-names(x)
#     nam2<-names(y) 
#     attr<-list()
#     if (is.null(nam1)) {names(x)<-nam1<-paste("var",1:lenl,sep="")}
#     if (is.null(nam2)) {names(y)<-nam2<-paste("var",1:lenl2,sep="")} 
#     # verificar que todos los nam1 estan en nam2
#     # verificar que longitud de metric y par.metric es la misma q x
#     if (missing(metric)){
#       if (is.fdata(x[[nam1[1]]])) metric<-rep("metric.lp",len=lenl)
#       else    metric<-rep("metric.dist",len=lenl)   
#     }
#     if (missing(w)) w<-rep(1,length=lenl)
#     
#     for (i in 1:lenl){
#       # print("metric.x3")
#       if (is.character(metric)) {
#         if (is.null(par.metric[[nam1[i]]])) par.metric[[nam1[i]]]<-list()
#         if (is.fdata(x[[nam1[i]]])) { 
#           par.metric[[nam1[i]]][["fdata1"]]<-x[[nam1[i]]]
#           par.metric[[nam1[i]]][["fdata2"]]<-y[[nam1[i]]]   
#           
#           mdist<-do.call(metric[i],par.metric[[nam1[i]]])
#         }
#         else {
#           par.metric[[nam1[i]]][["x"]]<-x[[nam1[i]]]
#           par.metric[[nam1[i]]][["y"]]<-y[[nam1[i]]]   
#           mdist<-do.call(metric[i],par.metric[[nam1[i]]])
#         }
#       }
#       else {
#         if (is.list(metric)) mdist<-metric[[nam1[i]]]
#       }
#       ldist[[nam1[i]]]<-mdist
#       amdist[,,i]<- mdist
#       #attr[[nam1[i]]]<-attributes(mdist)     
#     }
#     # print("sale1")
#     #print(Weights)
#     ###for (k in len) sqrt(w[i])*
#     # print(w)
#     mdist<-apply(amdist,1:2,method,w=w)
#     # print("sale2")
#   }
#   else stop("Error in x argument")
#   # print("a")
#   # print(dim(amdist))
#   # print(dim(mdist))
#   # print(attributes(mdist))
#   # print(attr)
#   #attributes(mdist)<-atr
#   attr(mdist, "method") <-method
#   attr(mdist, "w") <- w
#   for (i in 1:lenl) attr(mdist, nam1[i]) <- attributes(ldist[[nam1[i]]])
#   return(mdist)
# }
################################################################################
# deprecated 20190614
# metric.ldata=function(x,y=x
#                       ,metric
#                       ,par.metric=list()
#                       ,w
#                       ,method="euclidean") {
# if (any(class(x)=="list")){
#  lenl<-length(x)
#  lenl2<-length(y)
#  n<-nrow(x[[1]])
#  m<-nrow(y[[1]]) 
#  mdist2<-matrix(0,n,n)
#  amdist<-array(NA,dim=c(n,m,lenl))
#  ldist<-mdist<-list()
#  nam1<-names(x)
#  nam2<-names(y) 
#  
#  
#  attr<-list()
#  if (is.null(nam1)) {names(x)<-nam1<-paste("var",1:lenl,sep="")}
#  if (is.null(nam2)) {names(y)<-nam2<-paste("var",1:lenl2,sep="")} 
#  # verificar que todos los nam1 estan en nam2
#  # verificar que longitud de metric y par.metric es la misma q x
#  if (missing(metric)){
#    if (is.fdata(x[[nam1[1]]])) metric<-rep("metric.lp",len=lenl)
#    else    metric<-rep("metric.dist",len=lenl)   
#  }
# if (missing(w)) w<-rep(1,length=lenl)
# # print("metric.x2")
#  
#  # stp <- FALSE
#  # cat("metric.ldata ncores:",ncores)
# ldist <- foreach(i = 1:lenl, .combine = 'c' ) %dopar% { 
# if (is.character(metric)) {
#  if (is.null(par.metric[[nam1[i]]])) par.metric[[nam1[i]]]<-list()
#   if (is.fdata(x[[nam1[i]]])) { 
#    par.metric[[nam1[i]]][["fdata1"]]<-x[[nam1[i]]]
#    par.metric[[nam1[i]]][["fdata2"]]<-y[[nam1[i]]]
#    mdist<-do.call(metric[i],par.metric[[nam1[i]]])
#    }
#   else {
#    par.metric[[nam1[i]]][["x"]]<-x[[nam1[i]]]
#    par.metric[[nam1[i]]][["y"]]<-y[[nam1[i]]]   
#    mdist<-do.call(metric[i],par.metric[[nam1[i]]])
#    }
# } else {
#   if (is.list(metric)) mdist<-metric[[nam1[i]]]
#   }
#   aux <- list(mdist)
# } #end foreach
# names(ldist)<-nam1
# if (method=="list") return(ldist)
# 
# for (i in 1:lenl) amdist[,,i]<- ldist[[nam1[i]]]
# mdist<-apply(amdist,1:2,method,w=w)
# }
# else stop("Error in x argument")
# attr(mdist, "method") <-method
# attr(mdist, "w") <- w
# for (i in 1:lenl) attr(mdist, nam1[i]) <- attributes(ldist[[nam1[i]]])
# #if (stp) suppressWarnings(stopCluster(cl))
# return(mdist)
# }
################################################################################
# created 20190614 used in metric.ldata, metric.mfdata
list.select<-function(x,include,exclude){
  #Internal, used in metric.mfdata
  var.name <- names(x)
  a<-sapply(x, class,USE.NAMES = T)
  var.name <- var.name[a=="fdata"]
  clas.ini <- class(x)
  class(x) <- "list"
  x <- x[var.name]
  if (include[1]=="all"){
    var.x <- var.name
  } else{
    var.x <- intersect(include,var.name)
  }
  if (exclude[1]!="none"){
    var.x <- setdiff(var.x,exclude)
  }
  x <- x[var.x]
  clas.ini -> class(x)
  x
}
################################################################################
# created 20190614 used in metric.ldata
df.select <- function(x,include,exclude){
  #Internal, used in metric.mfdata
  # character variables are automatically excluded
  #input: x a data.frame
  #output: list with a variables of data.frame included and not excluded
  
  if (missing(include)) include = "all"
  if (missing(exclude)) exclude = "none"
  if (is.matrix(x)) x <- data.frame(x)
  #if (is.data.frame(x))  var.name <- names(x)
  #else var.name <- colnames(x)
  var.name <- names(x)
  a<-sapply(x, class,USE.NAMES = T)
  x <- x[,a!="character",drop=F]
  var.name <- names(x)
  #names(x)
  
  if (include[1]!="all"){
    var.name<-intersect(include,var.name)
  }
  if (exclude[1]!="none"){
    var.name<-setdiff(var.name,exclude)
  }
  x <- x[,var.name,drop=F]
  #names(x)<-var.name
  lx<-list()
  if (ncol(x)==0) return(x)
  for (k in 1:NCOL(x)) {
    if (is.factor(x[,var.name[k],drop=F])) lx[[var.name[k]]]<-model.matrix(~x[,var.name[k],drop=F]) #transforma el factor en dummyies
    else  lx[[var.name[k]]]<-x[,var.name[k],drop=F]
  }
  lx
}

# a1<-df.select(ldat3$df,exclude=")
# print(dim(faa$df))
# print(dim(fbb$df))
# df.select (si se exluye la unica variable qu hay peta)
# a1<-df.select(ldat3$df,exclude="glearn")
################################################################################

#' Distance Matrix Computation for ldata and mfdata class object 
#' 
#' This function computes the distances between the list elements 
#' This function returns a distance matrix by using \code{\link{metric.lp}}
#' function for \code{fdata} objects and  \code{\link{metric.dist}}
#' function for \code{vector} and \code{matric} objects. \cr 
#' 
#' @aliases metric.mfdata metric.ldata
#' @param ldata1 List with  of fdata objects and a data.frame object calle 'df'.
#' @param ldata2 List with  of fdata objects and a data.frame object calle 'df'.
#' @param method The distance measure to be used. This must be one of
#' "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski".
#' @param include vector with the name of variables to use
#' @param exclude vector with the name of variables to not use
#' @param metric Type of metric to combine, if 'none', the function no combine and return a list o distances for each variable included
#' @param par.metric List of metric parameters for each variable included
#' @param w, weights to combine the metric (if metric is not 'none')
#' @param \dots Further arguments passed to \code{\link{dist}} function.
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manul.oviedo@@usc.es}
#' @seealso See also \code{\link{dist}} for multivariate date case and
#' \code{\link{metric.lp} for functional data case}
#' @keywords cluster
#' @examples
#' \dontrun{
#' data(tecator)
#' names(tecator)[2]<-"df"
#' # Example 1 (list of distances)
#' ldist <- metric.ldata(tecator,method="none")
#' lapply(ldist,names)
#' # Example 2 (combined metric)
#' mdist <- metric.ldata(tecator,method="euclidean")
#' dim(mdist) 
#' }
#' @export metric.ldata
metric.ldata=function(ldata1,ldata2=NULL, include="all" ,exclude="none"
                       ,metric,par.metric=NULL,w, method="none") {
  
  one.mfdata=F
  if (is.null(ldata2)) {
    ldata2 <- ldata1
    one.mfdata=T  
  }
  var.ldata1<-names(ldata1)
   if (!any(var.ldata1=="df")) {
     #  stop("The object 'df' is not in the list argument 'ldata1' ")
    lenl <- length(ldata1)
    if (missing(w)) 
      w<-rep(1,length=lenl)
    n<-NROW(ldata1[[1]])
    if (missing(metric)){
      metric<-rep("metric.lp",len=lenl)
    }
    return(metric.mfdata(ldata1,ldata2,include=include,exclude=exclude
                        ,metric=metric,par.metric=par.metric
                          ,w=w,method=method))
   }
  ldf.ldata1 <- df.select(ldata1[["df"]],include,exclude)
  lendf <-length(ldf.ldata1)
  ldata1 <- c(ldf.ldata1,list.select(ldata1,include,exclude))
  var.ldata1 <- names(ldata1)
  lenl <- length(ldata1)
  n <- NROW(ldata1[[1]])
  if (missing(w)) 
     w<-rep(1,length=lenl)
  if (one.mfdata){
    ldata2 <- ldata1
    var.ldata2 <- var.ldata1
    m<-n
    lenl2<-lenl
  } else{
    ldata2 <- c(as.list(df.select(ldata1[["df"]],include,exclude)),list.select(ldata2,include,exclude))
    var.ldata2 <- names(ldata2)
    lenl2<-length(ldata2)
    m<-NROW(ldata2[[1]]) 
  } 
  
  if (!is(ldata1,"list"))
    stop("ldata1 is no a list object")
  if (!is(ldata2,"list"))
    stop("ldata2 is no a list object")
  
  mdst2<-matrix(0,n,n)
  ldist<-mdist<-list()
  attr<-list()
  
  # verificar que todos los nam1 estan en nam2
  # verificar que longitud de metric ldata2 par.metric es la misma q ldata1
  if (missing(metric)){
    metric<-c(rep("metric.dist",len=lendf),rep("metric.lp",len=(lenl-lendf)))
    arg.ldata1<-c(rep("x",len=lendf),rep("fdata1",len=(lenl-lendf)))
    arg.ldata2<-c(rep("y",len=lendf),rep("fdata2",len=(lenl-lendf)))
  }
  if (missing(w)) 
    w<-rep(1,length=lenl)
  par.fda.usc <- eval(parse(text="fda.usc:::par.fda.usc"), envir=.GlobalEnv)
  ncores <- par.fda.usc$ncores
  int.method <- par.fda.usc$int.method
  
  ldist<-list()
  for (i in 1:lenl) {
  if (is.null(par.metric[[var.ldata1[i]]])) 
    par.metric[[var.ldata1[i]]]<-list()
}
#print(arg.ldata1)
  ldist <- foreach(i = 1:lenl, .combine = 'c' ) %dopar% { 
    if (is.character(metric)) {
      if (is.null(par.metric[[var.ldata1[i]]])) 
        par.metric[[var.ldata1[i]]]<-NULL
      par.metric[[var.ldata1[i]]][[arg.ldata1[i]]]<-ldata1[[var.ldata1[i]]]
      # Ahorro en calcular la matriz triangular superior
      if (!one.mfdata)
        par.metric[[var.ldata1[i]]][[arg.ldata2[i]]]<-ldata2[[var.ldata1[i]]]
      mdist<-do.call(metric[i],par.metric[[var.ldata1[i]]])
    } else {
      if (is.list(metric)) 
        mdist<-metric[[var.ldata1[i]]]
    }
    aux <- list(mdist)
    aux
  } #end foreach
  names(ldist)<-var.ldata1
  if (method=="none") return(ldist)
  else{
    amdist<-array(NA,dim=c(n,m,lenl))
    for (i in 1:lenl) amdist[,,i]<- ldist[[var.ldata1[i]]]
    mdist<-apply(amdist,1:2,method,w=w)
    mdist<-apply(amdist,1:2,method,w=w)
  }
  #else stop("Error in ldata1 argument")
  attr(mdist, "method") <- method
  attr(mdist, "w") <- w
  for (i in 1:lenl) attr(mdist, var.ldata1[i]) <- attributes(ldist[[var.ldata1[i]]])
  return(mdist)
}




#' @export metric.mfdata
metric.mfdata=function(mfdata1,mfdata2=NULL
                      ,include="all"
                      ,exclude="none"
                      ,metric
                      ,par.metric=list()
                      ,w
                      ,method="none"#"euclidean"
                      ) {
  one.mfdata=F
  if (is.null(mfdata2)) {
    mfdata2 <- mfdata1
    one.mfdata=T  
  }
 
  mfdata1 <- list.select(mfdata1,include,exclude)
  var.mfdata1<-names(mfdata1)
  lenl<-length(mfdata1)
  n<-nrow(mfdata1[[1]])
 
  if (one.mfdata){
    mfdata2 <- mfdata1
    var.mfdata2 <- var.mfdata1
    m<-n
    lenl2<-lenl
  } else{
    mfdata2 <- list.select(mfdata2,include,exclude)
    var.mfdata2 <- names(mfdata2)
    lenl2<-length(mfdata2)
    m<-nrow(mfdata2[[1]]) 
  } 
  if (any(class(mfdata1)!="list"))
    stop("mfdata1 is no a list object")
  if (any(class(mfdata2)!="list"))
    stop("mfdata2 is no a list object")
  
  mdst2<-matrix(0,n,n)
  ldist<-mdist<-list()
  attr<-list()
  #if (is.null(nam1)) {names(mfdata1)<-nam1<-paste("var",1:lenl,sep="")}
  #if (is.null(nam2)) {names(mfdata2)<-nam2<-paste("var",1:lenl2,sep="")} 

  # verificar que todos los nam1 estan en nam2
  # verificar que longitud de metric mfdata2 par.metric es la misma q mfdata1
  if (missing(metric)){
   metric<-rep("metric.lp",len=lenl)
  }
  if (missing(w)) 
    w<-rep(1,length=lenl)
 # if (!getDoParRegistered())     ops.fda.usc()
  par.fda.usc <- eval(parse(text="fda.usc:::par.fda.usc"), envir=.GlobalEnv)
  ncores <- par.fda.usc$ncores
  int.method <- par.fda.usc$int.method
  
  ldist <- foreach(i = 1:lenl, .combine = 'c' ) %dopar% { 
      if (is.character(metric)) {
        if (is.null(par.metric[[var.mfdata1[i]]])) 
          par.metric[[var.mfdata1[i]]]<-list()
        par.metric[[var.mfdata1[i]]][["fdata1"]]<-mfdata1[[var.mfdata1[i]]]
        par.metric[[var.mfdata1[i]]][["fdata2"]]<-mfdata2[[var.mfdata1[i]]]
        mdist<-do.call(metric[i],par.metric[[var.mfdata1[i]]])
      } else {
        if (is.list(metric)) mdist<-metric[[var.mfdata1[i]]]
      }
      aux <- list(mdist)
    } #end foreach
    names(ldist)<-var.mfdata1
    if (method=="none") return(ldist)
    else{
      amdist<-array(NA,dim=c(n,m,lenl))
      for (i in 1:lenl) amdist[,,i]<- ldist[[var.mfdata1[i]]]
      mdist<-apply(amdist,1:2,method,w=w)
      mdist<-apply(amdist,1:2,method,w=w)
    }
  
  #else stop("Error in mfdata1 argument")
  attr(mdist, "method") <- method
  attr(mdist, "w") <- w
  for (i in 1:lenl) attr(mdist, var.mfdata1[i]) <- attributes(ldist[[var.mfdata1[i]]])
  return(mdist)
}




#' Distance Matrix Computation
#' 
#' This function computes the distances between the rows of a data matrix by
#' using the specified distance measure.
#' 
#' This function returns a distance matrix by using \code{\link{dist}}
#' function. \cr The matrix dimension is (\code{n1} x \code{n1}) if
#' \code{y=NULL}, (\code{n1} x \code{n2}) otherwise.
#' 
#' @param x Data frame 1. The dimension is (\code{n1} x \code{m}).
#' @param y Data frame 2. The dimension is (\code{n2} x \code{m}).
#' @param method The distance measure to be used. This must be one of
#' "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski".
#' @param p The power of the Minkowski distance.
#' @param dscale If scale is a numeric, the distance matrix is divided by the
#' scale value. If scale is a function (as the mean for example) the distance
#' matrix is divided by the corresponding value from the output of the
#' function.
#' @param \dots Further arguments passed to \code{\link{dist}} function.
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@usc.es}
#' @seealso See also \code{\link{dist}} for multivariate date case and
#' \code{\link{metric.lp} for functional data case}
#' @keywords cluster
#' @examples 
#' \dontrun{
#' data(iris)
#' d<-metric.dist(iris[,1:4])
#' matplot(d,type="l",col=as.numeric(iris[,5]))
#' }
#' @export
metric.dist<-function(x,y=NULL,method="euclidean",p=2,dscale=1,...){
if (is.vector(x)) x<-matrix(x,nrow=1)
else x<-as.matrix(x)
ynull<-is.null(y)
if (method=="mahalanobis"){
    if (ynull)   {
        y <- x
        vc <- var(x)
        }
    else {
    y <- as.matrix(y)
    vc <- var(rbind(x, y))
    }
    n <- nrow(x)
    m <- nrow(y)
    mdist<- matrix(0, n, m)
    for (i in 1:m) {
        mdist[, i] <- mahalanobis(x, y[i, ], cov = vc)
    } 
    mdist<-sqrt(mdist)
}
else{
 if (!ynull) {    
    if (is.vector(y)) y<-matrix(y,nrow=1) 
    n<-nrow(y)
    nn<-nrow(x)
    mdist<-as.matrix(dist(rbind(x,y) , method = method, diag = TRUE, upper = TRUE,p=p))[1:nn,(nn+1):(nn+n)] 
    }
 else   mdist<-as.matrix(dist(x, method = method, diag = TRUE, upper = TRUE,p=p))  
 }
 if (is.vector(mdist)) mdist<-matrix(mdist,nrow=nn)   
 	etiq1=rownames(x)
	etiq2=rownames(y)
# namesx<-rownames(x)
# if (ynull) dimnames(mdist) <- list(namesx,namesx)
# else dimnames(mdist) <- list(namesx, rownames(y))
 if (is.function(dscale)) {
   if (nrow(mdist)==ncol(mdist)) diag(mdist)<-NA ################# ojjjooo solo para matrices cuadradas
   dscale<-dscale(as.dist(mdist))
   if (nrow(mdist)==ncol(mdist)) diag(mdist)<-0   
 } 
 mdist<-mdist/dscale
 attr(mdist, "call") <- "metric.dist"
 attr(mdist, "par.metric") <- list(method =method,p=p,dscale=dscale) 
 rownames(mdist)<-etiq1
 colnames(mdist)<-etiq2
 return(mdist)
}
#####################
