#' @aliases na.omit.fdata na.fail.fdata
#' @title A wrapper for the na.omit and na.fail function for fdata object
#' 
#' @description   \code{na.fail} returns the object if it does not contain any 
#' missing values, and signals an error otherwise. \code{na.omit} returns the object
#' with incomplete cases removed.
#' If \code{na.omit.fdata} removes cases, the row numbers of the cases form the 
#' \code{"na.action"} attribute of the result, of class \code{"omit"}, see generic 
#' function \code{\link{na.omit}}.
#' 
#' @param object an \code{fdata} object.
#' @param \dots further potential arguments passed to methods.
#' @return The value returned from \code{omit} is a \code{fdata} object with incomplete cases removed.
#' @author Manuel Febrero Bande and Manuel Oviedo
#' @examples
#' \dontrun{
#' fdataobj<-fdata(MontrealTemp)
#' fdataobj$data[3,3]<-NA
#' fdataobj$data[10,]<-NA
#' fdastaobj2<-na.omit(fdataobj)
#' } 
#' 
#' @keywords descriptive
#' @rdname  na.omit.fdata
#' @export
na.omit.fdata=function(object,...){
  n=nrow(object)
  r= is.na.fdata(object) #apply(object[["data"]],1,count.na)
  omit=which(r)
  xx=object
  if (length(which(r))>0L) {
  xx=object[-omit]
  temp=setNames(c(1:n)[which(r)],attr(object[["data"]],"dimnames")[[1]][omit])
  attr(temp,"class")<-"omit"                                                                                                      
  attr(xx,"na.action")<-temp
  }
  return(invisible(xx))
}

#' @rdname na.omit.fdata
#' @export
na.fail.fdata=function(object,...){
  ok<-complete.cases(object[["data"]])
  if (all(ok)) 
  invisible(	object)
  else stop("missing values in fdata object")
  }
