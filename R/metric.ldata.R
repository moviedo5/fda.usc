euclidean <- function(u, w=rep(1,length(u))) {  
  u <- sqrt(sum(w*u^2,na.rm=TRUE))
  st <- attributes(u)
  attributes(u) <- st
  attr(u, "weigthts") <- w
  u
}
manhattan <- function(u, w=rep(1,length(u))) {  
  u<- sqrt(sum(w * abs(u), na.rm=TRUE))
  st <- attributes(u)
  attributes(u) <- st
  attr(u, "weigthts") <- w
  u
}
minkowski <- function(u, w=rep(1,length(u)), p) {  
  u <- (sum(w*(u^p), na.rm=TRUE))^(1/p)
  st <- attributes(u)
  attributes(u) <- st
  attr(u, "weigthts") <- w
  u
}
maximum <-function(u,w=rep(1,length(u))) {  
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
# Devulve una lista con las distancias para cada variable del "df"
# y el resto de objetos fdata
dist.ldata <- function(ldata, metric=NULL, par.metric=NULL,...){
  if (!is.ldata(ldata)) stop("No ldata object")
  lenldata <- length(ldata)
  ldist <- list()
  nam.ldata <- names(ldata)
  # nometrica <- is.null(metric)
  if ("df" %in% nam.ldata){
    idf <- which(nam.ldata=="df")
    if (is.data.frame(ldata$df)) nam.df <- names(ldata$df)
    else nam.df <- colnames(ldata$df)
    for (i in 1:(NCOL(ldata$df))) {
      if (is.factor(ldata$df[,i])) 
        aux <- model.matrix(~ldata$df[,i]) #transforma el factor en dummyies
      else 
        aux <- ldata$df[,i]
      if (is.null(metric[[nam.df[i]]]))
        ldist[[nam.df[i]]] <- as.matrix(dist(aux), diag =TRUE, upper = TRUE, ...)
      else {
        par.metric0 <- par.metric[[nam.df[i]]]
        par.metric0$x <- aux
        ldist[[nam.df[i]]] <- do.call(metric[[nam.df[i]]],par.metric0)
      }
    }
    nam.ldata <- setdiff(nam.ldata,"df")
  }
  for (i in 1:length(nam.ldata)) {
      # if (nometrica) 
    if (is.null(metric[[nam.ldata[i]]]))
      ldist[[nam.ldata[i]]] <- metric.lp(ldata[[nam.ldata[i]]],...)
      else {
        par.metric0 <- par.metric[[nam.ldata[i]]]
        #print( par.metric[[nam.ldata[i]]])
        par.metric0$fdata1 <- ldata[[nam.ldata[i]]]
        ldist[[nam.ldata[i]]] <- do.call(metric[[nam.ldata[i]]],par.metric0)
      }
  }
  # names(ldist)<-names(ldata)
  ldist
}


################################################################################
# created 20190614 used in metric.ldata, metric.mfdata
list.select <- function(x,include="all",exclude="none"){
  #Internal, used in metric.mfdata
  var.name <- names(x)
  a <- sapply(x, class,USE.NAMES = T)
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
  a <- sapply(x, class,USE.NAMES = T)
  x <- x[,a!="character",drop=F]
  var.name <- names(x)
  #names(x)
  
  if (include[1]!="all"){
    var.name <- intersect(include,var.name)
  }
  if (exclude[1]!="none"){
    var.name <- setdiff(var.name,exclude)
  }
  x <- x[,var.name,drop=F]
  #names(x)<-var.name
  lx <- list()
  if (ncol(x)==0) return(x)
  for (k in 1:NCOL(x)) {
    if (is.factor(x[,var.name[k],drop=F])) 
      lx[[var.name[k]]] <- model.matrix(~x[,var.name[k],drop=F]) #transforma el factor en dummyies
    else  lx[[var.name[k]]] <- x[,var.name[k],drop=F]
  }
  lx
}

################################################################################

#' @title Distance Matrix Computation for ldata and mfdata class object 
#' 
#' @description This function computes the distances between the list elements 
#' @details This function returns a distance matrix by using \code{\link{metric.lp}}
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
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manul.oviedo@@usc.es}
#' @seealso See also \code{\link{dist}} for multivariate date case and
#' \code{\link{metric.lp} for functional data case}
#' @keywords cluster
#' @examples
#' \dontrun{
#' data(tecator)
#' names(tecator)[2] <- "df"
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
    one.mfdata <- TRUE  
  }
  var.ldata1 <- names(ldata1)
   if (!any(var.ldata1=="df")) {
     #  stop("The object 'df' is not in the list argument 'ldata1' ")
    lenl <- length(ldata1)
    if (missing(w)) 
      w <- rep(1,length=lenl)
    n <- NROW(ldata1[[1]])
    if (missing(metric)){
      metric <- rep("metric.lp",len=lenl)
    }
    return(metric.mfdata(ldata1,ldata2,include=include,exclude=exclude
                        ,metric=metric,par.metric=par.metric
                          ,w=w,method=method))
   }
  ldf.ldata1 <- df.select(ldata1[["df"]],include,exclude)
  lendf <- length(ldf.ldata1)
  ldata1 <- c(ldf.ldata1,list.select(ldata1,include,exclude))
  var.ldata1 <- names(ldata1)
  lenl <- length(ldata1)
  n <- NROW(ldata1[[1]])
  if (missing(w)) 
     w <- rep(1,length=lenl)
  if (one.mfdata){
    ldata2 <- ldata1
    var.ldata2 <- var.ldata1
    m <- n
    lenl2 <- lenl
  } else{
    ldata2 <- c(as.list(df.select(ldata1[["df"]],include,exclude)),
                list.select(ldata2,include,exclude))
    var.ldata2 <- names(ldata2)
    lenl2 <- length(ldata2)
    m <- NROW(ldata2[[1]]) 
  } 
  
  if (!is.list(ldata1))
    stop("ldata1 is not a list object")
  if (!is.list(ldata2))
    stop("ldata2 is not a list object")
  
  mdst2 <- matrix(0,n,n)
  ldist <- mdist <- list()
  attr <- list()
  
  # verificar que todos los nam1 estan en nam2
  # verificar que longitud de metric ldata2 par.metric es la misma q ldata1
  if (missing(metric)){
    metric <- c(rep("metric.dist",len=lendf),rep("metric.lp",len=(lenl-lendf)))
    arg.ldata1 <- c(rep("x",len=lendf),rep("fdata1",len=(lenl-lendf)))
    arg.ldata2 <- c(rep("y",len=lendf),rep("fdata2",len=(lenl-lendf)))
  }
  if (missing(w)) 
    w <- rep(1,length=lenl)
  par.fda.usc <- eval(parse(text="fda.usc:::par.fda.usc"), envir=.GlobalEnv)
  ncores <- par.fda.usc$ncores
  int.method <- par.fda.usc$int.method
  
  ldist <- list()
  for (i in 1:lenl) {
  if (is.null(par.metric[[var.ldata1[i]]])) 
    par.metric[[var.ldata1[i]]] <- list()
}
#print(arg.ldata1)
  ldist <- foreach(i = 1:lenl, .combine = 'c' ) %dopar% { 
    if (is.character(metric)) {
      if (is.null(par.metric[[var.ldata1[i]]])) 
        par.metric[[var.ldata1[i]]] <- NULL
      par.metric[[var.ldata1[i]]][[arg.ldata1[i]]] <- ldata1[[var.ldata1[i]]]
      # Ahorro en calcular la matriz triangular superior
      if (!one.mfdata)
        par.metric[[var.ldata1[i]]][[arg.ldata2[i]]] <- ldata2[[var.ldata1[i]]]
      mdist <- do.call(metric[i],par.metric[[var.ldata1[i]]])
    } else {
      if (is.list(metric)) 
        mdist <- metric[[var.ldata1[i]]]
    }
    aux <- list(mdist)
    aux
  } #end foreach
  names(ldist) <- var.ldata1
  if (method=="none") return(ldist)
  else{
    amdist <- array(NA,dim=c(n,m,lenl))
    for (i in 1:lenl) amdist[,,i] <- ldist[[var.ldata1[i]]]
    mdist <- apply(amdist,1:2,method,w=w)
    mdist <- apply(amdist,1:2,method,w=w)
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
    one.mfdata <- TRUE  
  }
 
  mfdata1 <- list.select(mfdata1,include,exclude)
  var.mfdata1 <- names(mfdata1)
  lenl <- length(mfdata1)
  n <- nrow(mfdata1[[1]])
 
  if (one.mfdata){
    mfdata2 <- mfdata1
    var.mfdata2 <- var.mfdata1
    m <- n
    lenl2 <- lenl
  } else{
    mfdata2 <- list.select(mfdata2,include,exclude)
    var.mfdata2 <- names(mfdata2)
    lenl2 <- length(mfdata2)
    m <- nrow(mfdata2[[1]]) 
  } 
  if (!is.list(mfdata1))
    stop("mfdata1 is not a list object")
  if (!is.list(mfdata2))
    stop("mfdata2 is not a list object")
  
  mdst2 <- matrix(0,n,n)
  ldist <- mdist <- list()
  attr <- list()
  #if (is.null(nam1)) {names(mfdata1) <- nam1 <- paste("var",1:lenl,sep="")}
  #if (is.null(nam2)) {names(mfdata2) <- nam2 <- paste("var",1:lenl2,sep="")} 

  # verificar que todos los nam1 estan en nam2
  # verificar que longitud de metric mfdata2 par.metric es la misma q mfdata1
  if (missing(metric)){
   metric <- rep("metric.lp",len=lenl)
  }
  if (missing(w)) 
    w <- rep(1,length=lenl)
 # if (!getDoParRegistered())     ops.fda.usc()
  par.fda.usc <- eval(parse(text="fda.usc:::par.fda.usc"), envir=.GlobalEnv)
  ncores <- par.fda.usc$ncores
  int.method <- par.fda.usc$int.method
  
  ldist <- foreach(i = 1:lenl, .combine = 'c' ) %dopar% { 
      if (is.character(metric)) {
        if (is.null(par.metric[[var.mfdata1[i]]])) 
          par.metric[[var.mfdata1[i]]] <- list()
        par.metric[[var.mfdata1[i]]][["fdata1"]] <- mfdata1[[var.mfdata1[i]]]
        par.metric[[var.mfdata1[i]]][["fdata2"]] <- mfdata2[[var.mfdata1[i]]]
        mdist <- do.call(metric[i],par.metric[[var.mfdata1[i]]])
      } else {
        if (is.list(metric)) mdist <- metric[[var.mfdata1[i]]]
      }
      aux <- list(mdist)
    } #end foreach
    names(ldist) <- var.mfdata1
    if (method=="none") return(ldist)
    else{
      amdist <- array(NA,dim=c(n,m,lenl))
      for (i in 1:lenl) amdist[,,i] <- ldist[[var.mfdata1[i]]]
      mdist <- apply(amdist,1:2,method,w=w)
      mdist <- apply(amdist,1:2,method,w=w)
    }
  
  #else stop("Error in mfdata1 argument")
  attr(mdist, "method") <- method
  attr(mdist, "w") <- w
  for (i in 1:lenl) attr(mdist, var.mfdata1[i]) <- attributes(ldist[[var.mfdata1[i]]])
  return(mdist)
}


