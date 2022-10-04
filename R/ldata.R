#' @rdname ldata
#' @name ldata
#' @title ldata class definition and utilities  
#' 
#' @description ldata is a list with two type of objects:
#' \itemize{
#' \item \code{df} is a data frame with the multivariate data with n rows.
#' \item \code{...}  fdata objects of class \code{fdata} with n rows.
#' }
#' 
#' @aliases ldata c.ldata  is.ldata names.ldata [.ldata plot.ldata subset.ldata
#' ldata.cen  mean.ldata mean.fdata NCOL.ldata NROW.ldata ncol.ldata nrow.ldata
#' Ops.ldata Math.ldata Summary.ldata
# mean.mfdata is.mfdata nrow.mfdata mfdata.cen
#' @param ldata,x object of class \code{ldata}
# @param f index of observations
# @param drop passed on to respcetive function
#' @param i index
#' @param subset subset
#' @param \dots Further arguments passed to methods.
#' @param ask logilcal    If TRUE (and the R session is interactive) the user is asked for input, before a new figure is drawn. 
#' @param color colors to interpolate; must be a valid argument to  \code{colorRampPalette}.
#' @param var.name name of continuous univariate variable used in \code{color} argument
#' @param df data frame 
#' @param mfdata list of fdata objects
#' @param row logical If \code{FALSE}  (by default), \code{i} index selects 
#' the variables. If \code{TRUE},  \code{i} index selects the observations.
#' @examples 
#' data(tecator)
#' ab0 <- tecator$absorp.fdata
#' ab1 <- fdata.deriv(ab0)
#' ab2 <- fdata.deriv(ab0,nderiv=2)
#' ldat<-ldata(tecator$y,ab1=ab1,ab2=ab2)
#' is.ldata(ldat)
#' class(ldat)
#' plot(ldat[[1]])
#' plot(ldat[[2]]) 
#' # plot(ldat)
#' # plot(ldat,var.name="Fat")
#' @export 
ldata <-function(df,...,mfdata){
  C <- match.call()
  lenC <- length(C)
  mf<-TRUE
  if (missing(mfdata)){
    mf <- FALSE
  }  else lenC <-  lenC-1
  # mf <- match.call(expand.dots = FALSE)
  # m <- match("df", names(mf), 0L)
  nam <- as.character(C[3:lenC])
  mdat <- list(...)
  clases <- sapply(mdat,class)
    clases <- sapply(mdat,class)
  if (!all(clases=="fdata")) 
    stop("ldots must be fdata objects")
  n <- sapply(mdat,NROW)
  if (!all(n==n[1])) 
    stop("Different number of rows in the fdata objects")
  if (missing(df)) 
    df <- data.frame("index"=1:n[1])
  if (!is.data.frame(df)) 
    stop("df argument must be a data frame object")
  #names(mdat) <- nam
  if (missing(mfdata)) ldat <- c(list("df"=df),mdat)
  else   ldat <- c(list("df"=df),mdat,mfdata)

  class(ldat) <- c("ldata","list")
  return(ldat)
}
################################

#' @export 
c.ldata<-function (x, f, drop = FALSE, ...) 
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


#' @rdname ldata
#' @export 
names.ldata<-function(x){
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
    stop("No ldata class object")
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

#' @rdname ldata
#' @export is.ldata
is.ldata=function(x){
  nc=length(x)
  nam=names(x)
  n=NROW(x[[1]])
  lo1=all(unlist(sapply(x,NROW))==n)
  if (!lo1) warning("The number of rows are different")
  clas=unlist(lapply(x,class))
  l1=which(clas=="data.frame")
#  print(l1)
 # l1=which(nam=="data.frame")
  nam<-names(clas)
#  lo2=ifelse(length(l1)==0,FALSE,(length(l1)==1 & (nam[l1]=="df" | nam[l1]=="names.df")))
  lo2=ifelse(length(l1)==0,FALSE,(length(l1)==1 & (nam[l1]=="df")))
  lo3= ifelse(length(l1)==0,FALSE,all(clas[-l1]=="fdata"))
  #if (!lo3 & !lo2) warning("The structure list(df,fdata1,fdata2,...) is not correct")
  return( lo1 & lo2 & lo3)
}

#' @rdname ldata
#' @export 
"[.ldata" = function(x, i, row=FALSE){
  if (missing(i)) return(x)
  if (!row){
    var.name <- names(x)
    class(x) <- "list"
    hay.df <- FALSE
    if (is.numeric(i))  {
      x <- x[i]
      names(x)<-var.name[i]
      if (all(names(x)!="df"  ))         class(x) <- "list"
      else class(x) <- c("ldata","list")
      return(x)
    }
    ldf <- var.name == "df"
    if (any(ldf)) {
      idf <- var.name[which(ldf)]
      var.df <- names(x$df)
      int.df <- intersect(i,var.df)
      if (length(int.df)>0)     {
        hay.df <- TRUE
        x$df <- x$df[,int.df]
      }
    } 
    int.x <- intersect(i,var.name)
    if (hay.df)   x <- x[c(idf,int.x)]
    else   x <- x[i]
    return(x)
  } else{
  for (m in seq_along(x)){
    if (inherits(x[[m]],"data.frame")) {
      x[[m]] = x[[m]][i,]
    } else {
        x[[m]] = x[[m]][i]
        }
  }
  return(x)
  }
}

#' @rdname ldata
#' @export 
subset.ldata<-function(x, subset,...){
  #if (any(class(x)!="lfdata")) stop("No list class object")
  nvar<-length(x)
  namx=names(x)
  if (missing(subset)) subset<-!logical(nrow(x[[1]]))
  if (is.numeric(subset) & max(subset)<nrow(x[[1]])) {
    subset2<-logical(nrow(x[[1]]))
    subset2[subset]<-TRUE
    subset<-subset2    
  }  
  newx<-x
  for (i in 1:nvar){
    if (is.fdata(x[[i]])) {
      newx[[i]]<-subset.fdata(x[[i]],subset,drop=FALSE,...)
    } else {
      newx[[i]]<-subset(x[[i]],subset,drop=FALSE,...)	
    }
  } 
  return(invisible(newx))
}

#' @rdname ldata
#' @export
plot.ldata <- function(x, ask=FALSE, color, var.name,...){
  if (!is.ldata(x)) stop("No ldata class object")
  #if (is.ldata) stop("No ldata class object")
  col.bar=FALSE
  if (!missing(var.name)){
    var.color <- x$df[,var.name]
    col.bar<-TRUE
    if (!is.factor(var.color))  {
      if (missing(color)) color=c("red","blue")
      nticks <- 5
      color <- colorRampPalette(color, alpha = TRUE)(nticks)
      xfact <- cut( var.color,quantile( var.color,seq(0,1,len=nticks+1))
                 ,include.lowest=TRUE)
      min.col <- min(var.color)
      max.col <- max(var.color)
      ticks <- var.color
      
    }   else {
      nticks <- nlevels(var.color)
      if (missing(color)) color=1:nticks
      color <- colorRampPalette(color, alpha = FALSE)(nticks)
      xfact <- var.color
      min.col <- 0#min(as.numeric(var.color))
      max.col <- 1#max(as.numeric(var.color))
      min.col <- min(as.numeric(var.color))
      max.col <- max(as.numeric(var.color))
      ticks <- levels(var.color)
    }
  }  else     {
    if (missing(color)) color= 1:NROW(x$df)
    xfact <- 1:NROW(x$df)
    ticks <- NULL
  }
  mf <- 5
  nvar <- length(x)-1
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
  idf <- which(nam=="df")
  #plot(ts(ldf$df))
  #warning("Only fdata objects are plotted")
  nam<-nam[-idf]
  #nam<-names(x)$names.fdata
  
  #names1<-names2<-names<-x[[1]][["names"]]
  #names1$main<-"Multivariate Functional Data"
  
  for (idat in nam) {
    data<-x[[idat]]$data
    tt<-x[[idat]]$argvals
    rtt<-x[[idat]]$rangeval
    plot(x[[idat]], col =  color[xfact],
          #     main  =paste0(idat,"-",x[[idat]]$names$main),...)
          main  = x[[idat]]$names$main,...)
    
    #plot.fdata(ab,lty=1,col=colores[cfat],main="Original Trajectories")
    if (col.bar)    color.bar(color,min.col,max.col,ro=0,ticks=ticks)
    #if (col.bar)    color.bar(color,min.col,max.col,ro=0)
  }
}
################################################################################
# "[.lfdata"=function(lfdata,i){
#   if (missing(i)) return(lfdata)
#   res=lapply(lfdata,"[",i)
#   class(res)="lfdata"
#   return(res)
# }
# ################################################################################
# 
# ################################################################################
# is.lfdata=function(lfdata){
#   n=nrow(lfdata[[1]])
#   lo1=all(unlist(lapply(lfdata,nrow))==n)
#   if (!lo1) warning("The number of rows are different")
#   lo2=all(unlist(lapply(lfdata,is.fdata)))
#   if (!lo2) warning("Not all components are fdata")
#   return( lo1 & lo2 )
# }
################################################################################


################################################################################


##################################################################
#' @export
mean.ldata <- function (x,...){
  if (!is.ldata(x) | !is.list(x))
    stop("No ldata or list object")
  out <- lapply(x,func.mean)
  class(out) <- class(x)
  return(out)
}

##################################################################
#' @export
mean.fdata <- function (x,...) {
  if (!is.fdata(x))
    stop("No fdata object")
  x[["data"]] <- matrix(colMeans(x[["data"]], 
                        na.rm = TRUE), nrow = 1)
  x$names$main <- "mean"
  x
} 

##################################################################
#' @export
ldata.cen<- function (x, meanX = mean.ldata(x)) 
{
  # print("entra ldata.cen")
  if (!is.ldata(x) | !is.list(x))
    stop("No ldata or list object")
  #data <- fdataobj[["data"]]
  #fdataobj$data <- sweep(fdataobj$data, 2, meanX$data, FUN = "-")
  #clases <- sapply(x$df,class)
  #iclases <- (clases %in% c("numeric","integer","double"))
  nam <- names(x)
  aux1<-any(nam=="df") 
  if (!aux1)  stop("No 'df' element")  
  
  nam_df <- names(x$df)
  idf <- which(aux1)
  ifac <- !sapply(x$df,is.factor)
  # nam.fdata<-setdiff(nam,"df")
  lenl <- length(x)

  if (sum(ifac)>0){
  x$df[,nam_df[ifac]] <- sweep(data.matrix(x$df[,nam_df[ifac],drop=F]),
                    2, data.matrix(meanX$df[,nam_df[ifac],drop=F]), FUN = "-")
  }
  if (sum(aux1)==1)    {nam<- nam[-idf]}
  for (i in 1:length(nam)){
    x[[nam[i]]]$data <- sweep(x[[nam[i]]]$data
                          , 2, meanX[[nam[i]]]$data[1,], FUN = "-")
  }
  return(list(Xcen = x, meanX = meanX))
}

# names(xmean$df)
# (xmean$df)
# x <- ldata(aemet$df[ii,ivar],"temp"=aemet$temp[ii],"wind"=aemet$wind.speed[ii])
# xcen<-ldata.cen(x)
# names(xcen$meanX)
# 
# x <- ldata("df"=aemet$df[1:4,4:6],"temp"=aemet$temp[1:4],"wind.speed"=aemet$wind.speed[1:4])
# class(x)
# names(x)
# meanx<-ldata.mean(x)
# names(meanx)
# meanx$temp
# aux <- ldata.cen(x)

# 
# meanX = ldata.mean(x)
# names(meanX)
# aux <- ldata.cen(x)
# meanX = ldata.mean(aemet)
# class(meanX$df)
# aux<-ldata.cen(aemet)    
# aemet$df[1:3,4:8]
# aux$Xcen$df[1:3,]
# aux$meanX$df
# class(aux$meanX$df)


##################################################################


# @export
# lfdata=function(lfdata){
#   lo1=all(unlist(lapply(lfdata,nrow))==n)
#   if (!lo1) warning("The number of rows are different")
#   lo2=all(unlist(lapply(lfdata,is.fdata)))
#   if (!lo2) warning("Not all components are fdata")
#   return( lo1 & lo2 )
# }

############################################
# data(tecator)
# x <- tecator[[1]]
# tecator$y$cat<-factor(ifelse(tecator$y$Fat<15,2,4))
# ldf<-list("df"=tecator$y,"x"=x,"x.d1"=fdata.deriv(x))
# plot.mfdata(ldf,col=list("x"=2,"x.d1"=4))
# plot.mfdata(ldf,col=2:3)

# plot.lfdata<-function(lfdata,ask=FALSE,color,...){
#   mf=5
#   nvar<-length(lfdata)
#   if (nvar>4) ask=TRUE
#   if (ask) {par(mfrow = c(1, 1))
#     dev.interactive()
#     oask <- devAskNewPage(TRUE)
#     on.exit(devAskNewPage(oask))}
#   else{    mf<-switch(nvar,
#                       "1"={c(1,1)},
#                       "2"={c(1,2)},
#                       "3"={c(1,3)},
#                       "4"={c(2,2)})            
#   par(mfrow =mf)                    }
#   names1<-names2<-names<-lfdata[[1]][["names"]]
#   names1$main<-"Multivariate Functional Data"
#   #    tr<-paste("mode.tr",trim*100,"\u0025",sep="")         
#   nam<-names(lfdata)
#   
#   if (is.null(nam)) nam<-1:nvar
#   for (idat in 1:nvar) {
#     data<-lfdata[[idat]]$data
#     tt<-lfdata[[idat]]$argvals
#     rtt<-lfdata[[idat]]$rangeval
#     if (missing(color)) color2<-1
#     else {
#       if (is.list(color)) color2<-color[[idat]]
#       else color2<-color
#     }
#     plot(lfdata[[idat]], col =  color2,lty=1, main =nam[idat],...)
#   }
#   
# }    


##############
#' @export
nrow.ldata <- function (x) 
{
  if (is.ldata(x)){
  nc = length(x)
  nam = names(x)
  n = nrow(x[[1]])
  lo1 = all(unlist(sapply(x, nrow)) == n)
  if (!lo1) 
    warning("The number of rows are different")
  n
  } else stop("No ldata class object")
}

##############
#' @export
NROW.ldata <- function (x) 
{
  if (is.ldata(x)){
  nc = length(x)
  nam = names(x)
  n = NROW(x[[1]])
  lo1 = all(unlist(sapply(x, NROW)) == n)
  if (!lo1) 
    warning("The number of rows are different")
  n
  } else stop("No ldata class object")
}
##############
#' @export
NCOL.ldata <- function (x) 
{
  if (is.ldata(x))    sapply(x, NCOL)
  else stop("No ldata class object")
}
##############
#' @export
ncol.ldata <- function (x) 
{
  if (is.ldata(x))    sapply(x, ncol)
  else stop("No ldata class object")
}


##############
# '+.ldata'<-function (ldata1, ldata2) 
# {
#   inhe.fdata1 <- inherits(ldata1, "ldata")
#   inhe.fdata2 <- inherits(ldata2, "ldata")
#   clas = unlist(lapply(x, class))
#   
#   
#   if (!inhe.fdata1 && !inhe.fdata2) 
#     stop("Neither argument for + is a functional data object fdata.")
#   if (inhe.fdata1 && inhe.fdata2) {
#     
#    # NCOL.mfdata(ldata1)
#     lo1 = unlist(sapply(ldata1, NCOL)) 
#     lo2 = unlist(sapply(ldata2, NCOL)) 
#     lo1
#     lo2
#     
#     if (!lo1) 
#       warning("The number of rows are different")
#     
#   #  if (ncol(fdata1) != ncol(fdata2)) 
#   #    stop("Error in columns dimensions")
#     if (!(all(fdata1[["argvals"]] == fdata2[["argvals"]]))) 
#       stop("Error in argvals")
#     n1 = nrow(ldata1)
#     n2 = nrow(ldata2)
#     if (n1 == n2) {
#       fdataobj <- ldata1
#       fdataobj[["data"]] <- ldata1[["data"]] + 
#         fdata2[["data"]]
#     }
#     else {
#       if (n1 < n2) {
#         if (n2%%n1 != 0) 
#           stop("Number of rows not compatible")
#         fdataobj = fdata2
#         fdataobj[["data"]] = fdata2[["data"]] + 
#           matrix(rep(t(fdata1$data), n2%/%n1), ncol = ncol(fdata1), 
#                  byrow = TRUE)
#       }
#       else {
#         if (n1%%n2 != 0) 
#           stop("Number of rows not compatible")
#         fdataobj = fdata1
#         fdataobj[["data"]] = fdata1[["data"]] + 
#           matrix(rep(t(fdata2$data), n1%/%n2), ncol = ncol(fdata2), 
#                  byrow = TRUE)
#       }
#     }
#   }
#   if (!inhe.fdata1 && inhe.fdata2) {
#     fdataobj <- fdata2
#     fdataobj[["data"]] <- fdata1 + fdata2[["data"]]
#   }
#   if (inhe.fdata1 && !inhe.fdata2) {
#     fdataobj <- fdata1
#     fdataobj[["data"]] <- fdata1[["data"]] + 
#       fdata2
#   }
#   fdataobj
# }

#' @export 
Ops.ldata <- function(e1, e2){
  out<-Map(.Generic,e1,e1)
  class(out) <- class(e1)
  out
}


#' @export 
Math.ldata <- function(x, ...){
 # Math(x,...)
  lapply(x,.Generic,...)
}

#' @export 
Summary.ldata<- function(..., na.rm = FALSE){
#  lapply(...,.Generic, na.rm = FALSE)
  # Seguramente se quiera que se devuelva otra casa
  # por ejemplo el minimo para cada variable del df y por columnas del fdata
  out<-sapply(...,.Generic,2, na.rm = FALSE)
  #out<-c(out,sapply(...,.Generic,2, na.rm = FALSE))
  out
}
#min(ldat)

#abs(ldat)$df[1:4,]
#abs(ldat)$ab2[1,]
# ldat12 <- Ops.ldata(ldat1,ldat2)
# 
# 
# mapply("+", a, b, SIMPLIFY = FALSE)
# Map("+", a, b)
# 
# ldat1<-ldata(tecator$y,ab1=ab1,ab2=ab2)
# ldat2<-ldata(tecator$y,ab1=ab1,ab2=ab2)
# ldat3<-mean.ldata(ldat1)
# ldat3[[1]]
# ldata12 <- Map("+",ldat1,ldat2)
# ldata13 <- Map("+",ldat1,ldat3)
# 
# ldata12<-mapply("+", ldat1, ldat2, SIMPLIFY = FALSE)
# ldata12<-mapply("+", ldat1, ldat3, SIMPLIFY = FALSE)
# 
# names(ldata12)
# ldat1$df[1:3,]
# ldat2$df[1:3,]
# ldat12$df[1:3,]
# ldata12$df[1:3,]
# 
# ldat1$ab1$data[1:3,1:3]
# ldat2$ab1$data[1:3,1:3]
# ldata12$ab1$data[1:3,1:3]


#########################


##################################################################

# @export
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
# @export
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
################################################################################
# @export
mean.mfdata <- function (x,...){
  if (!is.mfdata(x) & !is.list(x))
    stop("No ldata or list object")
  out <- lapply(x,func.mean)
  class(out) <- class(x)
  return(out)
}

##############
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
##############
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
##############
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
##############
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
# mfdata.mean <- function (x){
#   if (!is.mfdata(x) | !is.list(x))
#     stop("No ldata or list object")
#   lapply(x,func.mean)
# }
#ldata.mean <- function (x){
# if (!is.ldata(x) | !is.list(x))
#       stop("No ldata or list object")
#   lapply(x,func.mean)
# }


