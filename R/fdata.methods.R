#' @name fdata.methods
#' @aliases Math.fdata Ops.fdata Summary.fdata split.fdata order.fdata is.fdata
#' 
#' @title  fdata  S3 Group Generic Functions
#' 
#' @description fdata Group generic methods defined for four specified groups
#' of functions, Math, Ops, Summary and Complex.
#'  
#' order.fdata and split.fdata: A wrapper for the order and split  function for fdata object.
#' 
#' @details 
#' In \code{order.fdata} the funcional data is ordered w.r.t the sample order of the values of vector.
#'  
#' \code{split.fdata} divides the data in the fdata object \code{x} into the groups defined by \code{f}. 
#' 
#' @param e1,e2 \code{fdata} class object
#' @param x An \code{fdata} object containing values to be divided into groups or an \code{list} of \code{fdata} objects  containing values to be combine by rows in a to be flatten one \code{fdata} object.
#' @param f a factor in the sense that as.factor(f) defines the grouping, or a list of such factors in which case their interaction is used for the grouping.
#' @param drop \code{logical} indicating if levels that do not occur should be dropped (if f is a factor or a list).
#' @param na.rm \code{logical}: should missing values be removed? 
#' @param  y A sequence of numeric, complex, character or logical vectors, all of the same length, or a classed R object. 
#' @param fdataobj \code{fdata} class object.
#' @param na.last for controlling the treatment of NAs. If TRUE, missing values in the data are put last; if FALSE, they are put first; if NA, they are removed; if "keep" they are kept with rank NA.   
#' @param decreasing \code{logical} Should the sort order be increasing or decreasing?.
#' @param \dots Further arguments passed to methods. 
#' 
#' @seealso  See  \link[base]{Summary} and \link[base]{Complex}.
#' @author   Manuel Febrero Bande and Manuel Oviedo de la Fuente \email{manuel.oviedo@@usc.es}
#' @return
#' \itemize{ 
#' \item split.fdata: The value returned from \code{split} is a list of fdata objects 
#' containing the values for the groups. The components of the list are named by
#'  the levels of f (after converting to a factor, or if already a factor and 
#'  drop = TRUE, dropping unused levels).\
#'  
# \item unlist_fdata:  The value returned from \code{unlist} is a \code{fdata} 
# object containing the \code{fdata} components of the list.
#' 
#' \item order.fdata:  returns the functional data \code{fdataobj}  w.r.t. a permutation 
#' which rearranges its first argument into ascending or descending order.
#' }
#' 
#' @examples
#' \dontrun{
#' data(tecator)
#' absor<-tecator$absorp.fdata
#' absor2<-fdata.deriv(absor,1)
#' absor<-absor2[1:5,1:4]
#' absor2<-absor2[1:5,1:4]
#' sum(absor)
#' round(absor,4)
#' log1<-log(absor)
#' 
#' fdataobj<-fdata(MontrealTemp)
#' fac<-factor(c(rep(1,len=17),rep(2,len=17)))
#' a1<-split(fdataobj,fac)
#' dim(a1[[1]]);dim(a1[[2]])
#' }
#' 
#' @keywords math descriptive
#' 

#' @rdname fdata.methods
#' @export 
Math.fdata<-function (x,...) {
   if (!inherits(x, "fdata")) stop("Objects must be of class fdata")
   x$data<-callGeneric(x$data@.Data,...)
   x
}

#' @rdname fdata.methods
#' @export 
Ops.fdata <- function (e1, e2 = NULL) {
inhe.fdata1<- inherits(e1, "fdata")
inhe.fdata2<- inherits(e2, "fdata")
if (!inhe.fdata1 && !inhe.fdata2)
  stop("Neither argument are fdata object .")
if (inhe.fdata1 && inhe.fdata2) {
   fdataobj<-e1
   if (!all(e1[["argvals"]] == e2[["argvals"]]))  stop("Error in argvals")
   if (all(dim(e1[["data"]])==dim(e2[["data"]]))) {
     fdataobj$data<-callGeneric(e1$data@.Data, e2$data@.Data)   }
   else {
     if (nrow(e2[["data"]])==1) {
           fdataobj$data<-sweep(e1$data, 2, e2$data, FUN = .Generic)
     }
     if (nrow(e1[["data"]])==1) {
           fdataobj$data<-sweep(e2$data, 2, e1$data, FUN = .Generic)
     } else    stop("Error in data dimenstion")
   }
   return(fdataobj)
   }
else {
if (!inhe.fdata1 && inhe.fdata2) {
 if (all(dim(e1)==dim(e2[["data"]]))) {
  fdataobj<-e2
  fdataobj$data<-callGeneric(e1@.Data, e2$data@.Data)
  }
 else {if (nrow(e2[["data"]])==1) {
           dat<-sweep(e1$data, 2, e2$data, FUN = .Generic)
     }}
 return(fdataobj)
}
if (inhe.fdata1 && !inhe.fdata2) {
 if (all(dim(e1)==dim(e2[["data"]]))) {
   fdataobj<-e1
   fdataobj$data<-callGeneric(e1$data@.Data,e2@.Data)
  }
 else {if (nrow(e2[["data"]])==1) {
           dat<-sweep(e1$data, 2, e2$data, FUN = .Generic)
     }}
 return(fdataobj)
}}
}

#' @rdname fdata.methods
#' @export 
Summary.fdata<-function(...,na.rm = FALSE) {
   aa<-match.call()
   fdataobj<-aa[[2]]
   if (!inherits(fdataobj, "fdata")) stop("Objects must be of class fdata")
   do.call(.Generic,list(fdataobj$data,na.rm=na.rm))
}

#' @rdname fdata.methods
#' @export 
split.fdata<-function(x,f,drop=FALSE,...){
  if (!is.factor(f)) f<-factor(f)
  nlev<-nlevels(f)
  lev<-levels(f)
  if (is.matrix(x$data)) x$data<-data.frame(x$data)
  if (is.fdata(x)) {
    out<-split(x$data,f,drop=drop,...)
  }
  for (i in 1:nlev) out[[lev[i]]]<-fdata(out[[lev[i]]],x$argvals,x$rangeval,x$names)
  out
}


#' @rdname fdata.methods 
#' @export
order.fdata<-function(y, fdataobj, na.last = TRUE, 
                      decreasing = FALSE){
  or<-order(y,na.last=na.last,decreasing=decreasing)
  fdataobj$data<-fdataobj$data[or,]
  fdataobj
}


#' @rdname fdata.methods
#' @usage is.fdata(fdataobj)
#' @export is.fdata
is.fdata<-function(fdataobj){
  if (!inherits(fdataobj, "fdata"))     return(FALSE)
  else     return(check.fdata(fdataobj))
}

check.fdata <- function(fdataobj){
  if (!is.list(fdataobj)) return(FALSE)
  
  if (sum(names(fdataobj) %in% c("data","argvals","rangeval","names"))!=4) {
    message("Inappropriate object names")
    return(FALSE)
  }
  if ( NCOL(fdataobj$data) != length(fdataobj$argvals))  {
    warning("Inappropriate argvals")
    return(FALSE)
  }
  if (min(fdataobj$argvals) < min(fdataobj$rangeval))   {
    message("Inappropriate rangeval minimum")
    return(FALSE)
  }
  if (max(fdataobj$argvals) > max(fdataobj$rangeval))   {
    message("Inappropriate  rangeval maximum")
    return(FALSE)
  }
#  type<-switch(is(fdataobj[["data"]],)[1],
#               matrix={if (ncol(fdataobj[["data"]])!=length(fdataobj[["argvals"]])) return(FALSE)},
#               data.frame={if (ncol(fdataobj[["data"]])!=length(fdataobj[["argvals"]])) return(FALSE)},
#               numeric={if (length(fdataobj[["data"]])!=length(fdataobj[["argvals"]])) return(FALSE)},
#               integer={if (length(fdataobj[["data"]])!=length(fdataobj[["argvals"]])) return(FALSE)})
  return(TRUE)
}

