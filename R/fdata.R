#' @name fdata
#' @title Converts raw data or other functional data classes into fdata class.
#' 
#' @description Create a functional data object of class \code{fdata} from (\code{matrix},
#' \code{data.frame}, \code{numeric}, \code{integer}, \code{fd, fds, fts} or
#' \code{sfts}) class data.
#' 
#' @param mdata Matrix of set cases with dimension (\code{n} x \code{m}), where
#' \code{n} is the number of curves and \code{m} are the points observed in
#' each curve.
#' @param argvals Argvals, by default: \code{1:m}.
#' @param rangeval (optional) Range of discretization points, by default:
#' range(\code{argvals}).
#' @param names (optional) list with tree components: \code{main} an overall
#' title, \code{xlab} title for \code{x} axis and \code{ylab} title for
#' \code{y} axis.
#' @param fdata2d TRUE class fdata2d, the functional data is observed in at
#' least a two grids (the \code{argvals} is a list of vectors). By default
#' \code{fdata2d=FALSE} the functional data is observed in a single grid (the
#' \code{argvals} is a vector).
#' 
#' @return Return \code{fdata} class object with: 
#' \itemize{ 
#' \item \code{"data"}: matrix of set cases with dimension (\code{n} x \code{m}),
#' where \code{n} is the number of curves and \code{m} are the points observed
#' in each curve 
#' \item \code{"rangeval"}: the discretizations points values, if not provided: \code{1:m} 
#' \item \code{"rangeval"}: range of the discretizations points values, by default: range(\code{argvals}) 
#' \item \code{"names"}: (optional) list with \code{main} an overall title, \code{xlab}
#' title for \code{x} axis and \code{ylab} title for \code{y} axis.
#' }
#' 
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' 
#' @seealso See Also as \code{\link{plot.fdata}}
#' 
#' @references Febrero-Bande, M., Oviedo de la Fuente, M. (2012).
#' \emph{Statistical Computing in Functional Data Analysis: The R Package
#' fda.usc.} Journal of Statistical Software, 51(4), 1-28.
#' \url{http://www.jstatsoft.org/v51/i04/}
#' 
#' @keywords manip
#' @examples 
#' \dontrun{
#' data(phoneme)
#' mlearn<-phoneme$learn[1:4,1:150]
#' # Center curves
#' fdata.c=fdata.cen(mlearn)$Xcen
#' par(mfrow=c(2,1))
#' plot(mlearn,type="l")
#' plot(fdata.c,type="l")
#' 
#' # Convert  from class fda to fdata
#' bsp1 <- create.bspline.basis(c(1,150),21)
#' fd1 <- Data2fd(1:150,y=t(mlearn$data),basisobj=bsp1)
#' fdataobj=fdata(fd1)
#' 
#' # Convert  from class fds, fts or sfts to fdata
#' #require(fds)
#' #a=fds(x = 1:20, y = Simulationdata$y, xname = "x", 
#' # yname = "Simulated value")
#' #b=fts(x = 15:49, y = Australiasmoothfertility$y, xname = "Age",
#' #    yname = "Fertility rate")
#' #c=sfts(ts(as.numeric(ElNino_ERSST_region_1and2$y), frequency = 12), xname = "Month",
#' #yname = "Sea surface temperature")
#' #class(a);class(b);class(c)
#' #fdataobj=fdata(b)
#' }
#' 
#' @rdname fdata
#' @export
fdata=function(mdata, argvals = NULL, rangeval = NULL,
               names = NULL, fdata2d = FALSE){
call<-match.call()
if (length(call[[2]])>1) nam<-"fdataobj"
else nam<-call[[2]]
if (is.array(mdata))  {
    dm<-dim(mdata)
    len.dm<-length(dm)
    if (len.dm<3) fdata2d=FALSE
    else {
      fdata2d=TRUE
      if (len.dm>3) stop("Not implemented yet for dim(mdata)>",3)
    }
}
if (is.list(argvals)) fdata2d=TRUE
if (is.list(rangeval)) fdata2d=TRUE

if (fdata2d){
out=list("data"=NULL)
#if (length(class(mdata))>1) class(mdata)<-class(mdata)[1]
out<-switch(is(mdata,)[1],
matrix={
        out[["data"]]<-array(NA,dim=c(1,dm[1],dm[2]))
        out[["data"]][1,,]<-mdata
                out},
data.frame={
        out[["data"]]<-array(NA,dim=c(1,dm[1],dm[2]))
        out[["data"]][1,,]<-mdata
        out},
fdata=stop("The data is a fdata class object"),
array={
        out[["data"]]=mdata
        out}
)
# if ( inherits(mdata,"matrix") | inherits(mdata,"data.frame") ){
#   out[["data"]]<-array(NA,dim=c(1,dm[1],dm[2]))
#   out[["data"]][1,,]<-mdata
# }else{
#   if (inherits(mdata,"fdata")) 
#     {stop("The data is a fdata class object")}  else {
#       if ( inherits(mdata,"array"))  out[["data"]]=mdata}
# }

dm<-dim(out[["data"]] )
len.dm<-length(dm)
if (is.null(argvals)) {
 argvals<-list()
 if (len.dm>2){
  len.argvals<-len.dm-1
  for (i in 2:len.dm) {
    if (is.null(dimnames(out[["data"]][i])))  argvals[[i-1]]<-1:dm[i]
  }
 }
else {
  len.argvals<-len.dm
  for (i in 1:len.dm) {
  if (is.null(dimnames(out[["data"]][i])))  argvals[[i]]<-1:dm[i]
# nam<-as.character(call[[2]])
# print(length(call[[2]]))
 if (is.null(names(argvals[[i]]))) names(argvals[[i]])<-paste(nam,1:dm[i],sep="")
  }
 }
 out[["argvals"]]<-argvals
}
else {
  #verificar dimensiones
  if (is.list(argvals)) {
     len.argvals<-length(argvals)
#     print(len.argvals)
     if (len.dm>2){
       for (i in 1:len.argvals) {
         if (length(argvals[[i]])!=dm[i+1]) stop("Incorrect dimension in between mdata and argvals arguments")
    }   }
    else {
       for (i in 1:len.argvals) {
        if (length(argvals[[i]])!=dm[i]) stop("Incorrect dimension in between mdata and argvals arguments")
    }  }
  }
else stop("The argument argvals must be a list")
#print("peta name")
  if (is.null(names(argvals))) names(argvals)<-paste(drop(nam),1:len.argvals,sep="")
  out[["argvals"]]<-argvals
}
##########################
if (is.null(rangeval))  {
 rangeval<-list()
 for (i in 1:len.argvals) {
    rangeval[[i]]<-range(argvals[i])
  }
 if (is.null(names(rangeval))) names(rangeval)<-names(argvals)
 out[["rangeval"]]<-rangeval
 len.rangeval<-length(out[["rangeval"]])
}
else {
  if (is.list(rangeval)) {
     len.rangeval<-length(rangeval)
     if (len.rangeval!=len.argvals) stop("Incorrect dimension in rangeval argument")
  }
else stop("The argument reangeval must be a list")
  if (is.null(names(rangeval))) names(rangeval)<-names(argvals)
  out[["rangeval"]]<-rangeval
}
if (is.null(names))  {
 names<-list()
 names[["xlab"]]<-paste("argvals ",1:len.argvals,sep="")
 names[["ylab"]]<-paste("values of ",nam,len.argvals,sep="")
 names[["main"]]<-paste(nam,len.argvals,sep="")
 out[["names"]]<-names
}
else {
  if (is.list(names)) {
     len.names<-length(names[["xlab"]])
    if (len.rangeval!=len.names) stop("Incorrect dimension in names argument")
  }
else stop("The argument names must be a list")
 out[["names"]]<-names
}
if (is.null(dimnames(out$data))) dimnames(out$data)<-list("data"=1:dm[1],"argvals.x"=argvals[[1]],"argvals.y"=argvals[[2]])
class(out) <- c("fdata","fdata2d")
return(out)
}   #### fdata1d
else{
out=list("data"=NULL)
# if (length(class(mdata))>1) class(mdata)<-class(mdata)[1]
out<-switch(is(mdata,)[1],
matrix={
        out[["data"]]=mdata
        out},
data.frame={
            out[["data"]]=data.matrix(mdata)
            out},
ftable={
  mdata<-matrix(mdata,nrow=dm[1])
  out[["data"]]<-mdata
  out},
xtabs={
  mdata<-matrix(mdata,nrow=dm[1])
  out[["data"]]<-mdata
  out},
fdata=stop("The data could not be converted into fdata class"),
numeric={
         out[["data"]]=matrix(mdata,nrow=1)
         out},
integer={
         out[["data"]]=matrix(mdata,nrow=1)
         out},
fd={
   r = mdata$basis[[3]]
   if (is.null(argvals))
     argvals= seq(r[1], r[2], len =mdata$basis$nbasis)
    nb<- length(argvals)
    tt = argvals
#   tt = seq(r[1], r[2], len = length(mdata$fdnames$time))

   out[["data"]] = t(eval.fd(tt, mdata))
   if (!is.null(mdata$fdnames$reps)) rownames(out[["data"]]) = mdata$fdnames$reps
   else      rownames(out[["data"]]) =1:nrow( out[["data"]])
   if (!is.null(mdata$fdnames$time)) {
      colnames(out[["data"]]) = 1:ncol( out[["data"]])
      }
   else      {   colnames(out[["data"]]) =1:ncol( out[["data"]])   }
   out
   },
fds={
 out[["data"]] =t(mdata$y)
 if (is.null(mdata$x))       out[["argvals"]] =1:ncol(out[["data"]])
 else { 
   if (is.numeric(mdata$x)) out[["argvals"]]=mdata$x 
   else out[["argvals"]]=seq(mdata$x[1], mdata$x[length(mdata$x)],
  len = length(mdata$x))
  }
  out[["rangeval"]]<-range(out[["argvals"]])  
  out[["names"]]<-list("main"=deparse(nam),"xlab"=mdata$xname,"ylab"=mdata$yname)
  class(out)<-"fdata"
  return(out)
  },
fts={
 out[["data"]] =t(mdata$y)
 if (is.null(mdata$x))       out[["argvals"]] =1:ncol(out[["data"]])
 else { 
   if (is.numeric(mdata$x)) out[["argvals"]]=mdata$x 
   else out[["argvals"]]=seq(mdata$x[1], mdata$x[length(mdata$x)],
  len = length(mdata$x))
  }
  out[["rangeval"]]<-range(out[["argvals"]])  
  out[["names"]]<-list("main"=deparse(nam),"xlab"=mdata$xname,"ylab"=mdata$yname)
  class(out)<-"fdata"
  return(out)
  },
sfts={
 out[["data"]] =t(mdata$y)
 if (is.null(mdata$x))       out[["argvals"]] =1:ncol(out[["data"]])
 else { 
   if (is.numeric(mdata$x)) out[["argvals"]]=mdata$x 
   else out[["argvals"]]=seq(mdata$x[1], mdata$x[length(mdata$x)],
  len = length(mdata$x))
  }
  out[["rangeval"]]<-range(out[["argvals"]])  
  out[["names"]]<-list("main"=deparse(nam),"xlab"=mdata$xname,"ylab"=mdata$yname)
  class(out)<-"fdata"
  return(out)
  }
)
nc<-nrow(out[["data"]])
np<-ncol(out[["data"]])
if (is.null(argvals)) {
 if (is.null(colnames(out[["data"]]))) {out[["argvals"]]=1:ncol(out[["data"]])}
   else {
#   if (!any(is.na(as.numeric(colnames(out[["data"]]))))) {
#    out[["argvals"]]=as.numeric(colnames(out[["data"]]))   }
#    else    out[["argvals"]]=1:ncol(out[["data"]])}   }
out[["argvals"]]=1:ncol(out[["data"]])}   }
else     out[["argvals"]]=argvals
lentt=length(out[["argvals"]])
if (is.null(rangeval)) rangeval=range(out[["argvals"]])
out[["rangeval"]]<-rangeval
if ((np!=lentt) && (nc==lentt)) {
         out[["data"]]=matrix(out[["data"]],ncol=nc)
         nc<-1
         warning("The transposed data is returned")       }
else    out[["data"]]=out[["data"]]
if (is.null(dimnames(mdata))) {
#rownames(out[["data"]])<-1:nc
colnames(out[["data"]])=round(out[["argvals"]],4)
}
out[["names"]]<-list("main"="fdataobj","xlab"="t","ylab"="X(t)")
if (!is.null(names$main)) out$names$main<-names$main
if (!is.null(names$xlab)) out$names$xlab<-names$xlab
if (!is.null(names$ylab)) out$names$ylab<-names$ylab
class(out)="fdata"
return(out)
}
}



 
 