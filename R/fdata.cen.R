#' @title Functional data centred (subtract the mean of each discretization point)
#' 
#' @description The function fdata.cen centres the curves by subtracting the functional
#' mean.
#'  
#' @param fdataobj \code{\link{fdata}} class object.
#' @param meanX The functional mean subtracted in the \code{fdatobj}.
#' 
#' @return Return:\cr two \code{fdata} class objects with: \item{Xcen}{ The
#' centered fdata.} \item{meanX}{ Functional mean substracted.}
#' 
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' 
#' @seealso See Also as \code{\link{fdata}}
#' 
#' @keywords manip
#' 
#' @examples
#' \dontrun{
#' data(phoneme)
#' mlearn<-phoneme[["learn"]][13:15,]
#' fdata.c=fdata.cen(mlearn)$Xcen
#' par(mfrow=c(1,2))
#' plot(mlearn,type="l")
#' plot(fdata.c,type="l")
#' }
#' 
#' @aliases fdata.cen
#' @export
fdata.cen=function(fdataobj, meanX = func.mean(fdataobj)){
 if (!is.fdata(fdataobj))  fdataobj=fdata(fdataobj)
 data<-fdataobj[["data"]]
 fdataobj$data <- sweep(fdataobj$data,2,meanX$data,FUN="-")
 return(list("Xcen"=fdataobj,"meanX"=meanX))
}

#fdata.cen22 <- function(fdata){
#n <- dim(fdata)[1]
#unos <- rep(1,n)
#meanX <- (1/n) * t(fdata) %*% unos
#Xcen <- fdata - unos %*% t(meanX)
#return(list("Xcen"=Xcen,"meanX"=meanX))
#}




