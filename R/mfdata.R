#' @rdname mfdata
#' @name mfdata
#' @title mfdata class definition and utilities  
#' 
#' @description mfdata is a list with fdata type of object:
#' \itemize{
#' \item \code{...}  fdata objects of class \code{fdata} with n rows.
#' }
#' 
#' @aliases mfdata c.mfdata  is.mfdata names.mfdata [.mfdata plot.mfdata subset.mfdata
#' mfdata.cen  mean.mfdata  NCOL.mfdata NROW.mfdata ncol.mfdata nrow.mfdata
#' Ops.mfdata Math.mfdata Summary.mfdata
# @param ldata x object of class \code{ldata}
# @param f index of observations
# @param drop passed on to respcetive function
# @param i index
#' @param x object of class \code{mfdata}
#' @param subset subset
#' @param \dots Further arguments passed to methods.
# @param ask logilcal    If TRUE (and the R session is interactive) the user is asked for input, before a new figure is drawn. 
# @param color colors to interpolate; must be a valid argument to  \code{colorRampPalette}.
# @param var.name name of continuous univariate variable used in \code{color} argument
# @param df data frame 
# @param mfdata list of fdata objects
# @param row logical If \code{FALSE}  (by default), \code{i} index selects 
# the variables. If \code{TRUE},  \code{i} index selects the observations.
#' 
#' @examples 
#' data(tecator)
#' ab0 <- tecator$absorp.fdata
#' ab1 <- fdata.deriv(ab0)
#' ab2 <- fdata.deriv(ab0,nderiv=2)
#' mdat <- mfdata(ab0, ab1, ab2)
#' is.ldata(mdat)
#' class(mdat)
#' plot(mdat[[1]])
#' plot(mdat[[2]]) 
#' plot(mdat)
# # plot(mdat,var.name="Fat")
#' @export 
mfdata <-function(...){
  C <- match.call()
  lenC <- length(C)
  mf<-TRUE
  #if (missing(mfdata)){
  #  mf <- FALSE
  #}  else lenC <-  lenC-1
  # mf <- match.call(expand.dots = FALSE)
  # m <- match("df", names(mf), 0L)
  nam <- as.character(C[2:lenC])
  # print(nam)
  mdat <- list(...)
  clases <- sapply(mdat,class)
  names(mdat) <- nam
  if (!all(clases=="fdata")) 
    stop("ldots must be fdata objects")
  n <- sapply(mdat,NROW)
  if (!all(n==n[1])) 
    stop("Different number of rows in the fdata objects")
  class(mdat) <- c("mfdata","list")
  return(mdat)
}
################################

#' @export 
c.mdata<-function (x, f, drop = FALSE, ...) 
{
  if (!is.list(x)) stop("x is not a list") 
  lenl<-length(x)
  out<-x
  lenf<-length(f)
  for (i in 1:lenl) {
    if (lenf>nrow(x[[i]])) stop("Incorrect length of f")
    out[[i]] <- x[[i]][f,]
  }
  out
}


#' @rdname mfdata
#' @export 
names.mfdata<-function(x){
  # x<-ldf
  if (is(x,"list")){
    class(x) <- "list"
    return(names(x))
    #nam.x <- names(x)
    #nam.df  <-names(x$df)
    #ind.df <- which(nam.x=="df")
    #name <-list("names.df" =names(x$df),"names.fdata"=nam.x[-ind.df])
    #name <-c(names(x$df),nam.x[-ind.df])
  }
  else 
    stop("No mfdata class object")
  #return(name)
}
# names.ldata <- function(x){
#   # x<-ldf
#   if (any(class(x)=="ldata")){
#     class(x)<-"list"
#     nam.x <- names(x)
#     nam.df  <-names(x$df)
#     ind.df <- which(nam.x=="df")
#     name <-list("names.df" =names(x$df),"names.fdata"=nam.x[-ind.df])
#   }else stop("No ldata class object")
#   return(name)
# }

#' @export
is.mfdata=function(x){
  nc=length(x)
  nam=names(x)
  n=NROW(x[[1]])
  lo1=all(unlist(sapply(x,NROW))==n)
  if (!lo1) warning("The number of rows are different")
  clas=unlist(lapply(x,class))
  
  if (all(clas=="fdata")) return(TRUE)
  else return(FALSE)
}


#' @export 
"[.mfdata" = function(x, i, row=FALSE){
  if (missing(i)) return(x)
  if (!row){
    var.name <- names(x)
    class(x) <- "list"
    hay.df <- FALSE
    if (is.numeric(i))  {
      x <- x[i]
      names(x)<-var.name[i]
      if (all(names(x)!="df"  ))         class(x) <- "list"
      else class(x) <- c("mfdata","list")
      return(x)
    }
    int.x <- intersect(i,var.name)
    x <- x[i]
    return(x)
  } else{
  for (m in seq_along(x)){
    x[[m]] = x[[m]][i]
  }
  return(x)
  }
}

# mdat2 <- mdat["ab1"]
# mdat12 <-mdat[1:3,row = T]
# class(mdat2)
# sapply(mdat,dim)
# sapply(mdat12,dim)
# names(mdat);names(mdat2)


#' @rdname mfdata
#' @export 
subset.mfdata<-function(x, subset,...){
  nvar <- length(x)
  namx <- names(x)
  if (missing(subset)) subset <- !logical(nrow(x[[1]]))
  if (is.numeric(subset) & max(subset) < nrow(x[[1]])) {
    subset2 <- logical(nrow(x[[1]]))
    subset2[subset] <- TRUE
    subset<-subset2    
  }  
  newx<-x
  for (i in 1:nvar){
    if (is.fdata(x[[i]])) {
      newx[[i]] <- subset.fdata(x[[i]],subset,drop=FALSE,...)
    } else {
      newx[[i]] <- subset(x[[i]],subset,drop=FALSE,...)	
    }
  } 
  return(invisible(newx))
}
# ldat2<-fda.usc:::subset.ldata(ldat,ldat$df$index<10)
# subset.fdata <- fda.usc:::subset.fdata
# mdat3<-subset.mfdata(mdat,ldat$df$index<10)
# sapply(ldat2,dim)
# sapply(mdat3,dim)

#' @export
plot.mfdata <- function(x, ask=FALSE, color, ...){
  if (!is.mfdata(x)) stop("No mfdata class object")
  col.bar=FALSE
  if (missing(color)) color= 1:nrow(x[[1]])
  xfact <- 1:NROW(x[[1]])
  ticks <- NULL
  mf <- 5
  nvar <- length(x)
  
  if (nvar > 4) ask=TRUE
  if (ask) {par(mfrow = c(1, 1))
    dev.interactive()
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }  else{    
    mf<-switch(nvar,
               "1"={c(1,1)},
               "2"={c(1,2)},
               "3"={c(1,3)},
               "4"={c(2,2)})            
    par(mfrow =mf)                    
  }
  #    tr<-paste("mode.tr",trim*100,"\u0025",sep="")         
  nam <- names(x)
  min.col <- 0#min(color)
  max.col <- 1#max(color)
  for (idat in nam) {
    data<-x[[idat]]$data
    tt<-x[[idat]]$argvals
    rtt<-x[[idat]]$rangeval
    var.color<-  color[xfact]
    plot(x[[idat]], col = color,
          #     main  =paste0(idat,"-",x[[idat]]$names$main),...)
          main  = x[[idat]]$names$main,...)
    if (col.bar)    color.bar(color,min.col,max.col,ro=0,ticks=ticks)
    }
}
#plot.mfdata(mdat)

##################################################################
#' @export
mean.mfdata <- function (x,...){
  if (!is.mfdata(x) | !is.list(x))
    stop("No mfdata or list object")
  out <- lapply(x,func.mean)
  class(out) <- class(x)
  return(out)
}
# mm <- mean.mfdata(mdat)



##################################################################

#' @export
mfdata.cen<- function (x, meanX = mean.mfdata(x)) 
{
#  print("mfdata.cen")
  if (is.fdata(x)) return(fdata.cen(x))
  if (!is.mfdata(x) | !is.list(x))
    stop("No mfdata or list object")
  #data <- fdataobj[["data"]]
  #fdataobj$data <- sweep(fdataobj$data, 2, meanX$data, FUN = "-")
  #clases <- sapply(x$df,class)
  #iclases <- (clases %in% c("numeric","integer","double"))
  nam <- names(x)
  lenl <- length(x)
  for (i in 1:lenl)
    x[[nam[i]]]$data <- sweep(x[[nam[i]]]$data
                              , 2, meanX[[nam[i]]]$data, FUN = "-")
  return(list(Xcen = x, meanX = meanX))
}

#' @export 
mean.mfdata <- function (x,...){
  if (!is.mfdata(x) & !is.list(x))
    stop("No ldata or list object")
  out <- lapply(x,func.mean)
  class(out) <- class(x)
  return(out)
}

#' @export 
nrow.mfdata <- function (x) 
{
  if (is.mfdata(x)){
    nc = length(x)
    nam = names(x)
    n = nrow(x[[1]])
    lo1 = all(unlist(sapply(x, nrow)) == n)
    if (!lo1) 
      warning("The number of rows are different")
    n
  } else stop("No mfdata class object")
}

#' @export 
NROW.mfdata <- function (x) 
{
  if (is.fdata(x)){
    nc = length(x)
    nam = names(x)
    n = NROW(x[[1]])
    lo1 = all(unlist(sapply(x, NROW)) == n)
    if (!lo1) 
      warning("The number of rows are different")
    n
  } else stop("No mfdata class object")
}

#' @export 
NCOL.mfdata <- function (x) 
{
  nc = length(x)
  nam = names(x)
  n = NCOL(x[[1]])
  lo1 = all(unlist(sapply(x, NCOL)) == n)
  if (!lo1) 
    warning("The number of columns are different")
  n
}

#' @export 
ncol.mfdata <- function (x) 
{
  nc = length(x)
  nam = names(x)
  n = ncol(x[[1]])
  lo1 = all(unlist(sapply(x, ncol)) == n)
  if (!lo1) 
    warning("The number of columns are different")
  n
}

#' @export 
Ops.mfdata <- function(e1, e2){
  out<-Map(.Generic,e1,e1)
  class(out) <- class(e1)
  out
}

#' @export 
Math.mfdata <- function(x, ...){
  # Math(x,...)
  lapply(x,.Generic,...)
}

#' @export 
Summary.mfdata<- function(..., na.rm = FALSE){
  #  lapply(...,.Generic, na.rm = FALSE)
  # Seguramente se quiera que se devuelva otra casa
  # por ejemplo el minimo para cada variable del df y por columnas del fdata
  out<-sapply(...,.Generic,2, na.rm = FALSE)
  #out<-c(out,sapply(...,.Generic,2, na.rm = FALSE))
  out
}
#########################


