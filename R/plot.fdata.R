#' @title Plot functional data: fdata class object
#' 
#' @description Plot object of class \code{fdata}. 
#' 
#' @aliases plot.fdata lines.fdata title.fdata plot.bifd plot.depth plot.mdepth
#' 
#' @param x \code{fdata} class object with:  \itemize{ 
#' \item \code{"data"}:  For \code{fdata} class object as curve (1d), \code{"data"} is a
#' \code{matrix} (by default), \code{data.frame} or \code{array} of set cases
#' with dimension (\code{n} x \code{m}), where \code{n} is the number of curves
#' and \code{m} are the points observed in each curve over the x--axe.\cr For
#' \code{fdata2d} class object as surface (2d). \code{"data"} is a \code{array}
#' of set cases with dimension (\code{n} x \code{m1} x \code{m2}), where
#' \code{n} is the number of functional data and \code{m1} and \code{m2} are
#' the points observed over the x--y plane.  \item \code{"argvals"}: vector or
#' list of vectors with the discretizations points values.  \item
#' \code{"rangeval"}: vector or list of vectors with the range of the
#' discretizations points values, by default range(\code{argvals}).  \item
#' \code{"names"}: (optional) list with \code{main} an overall title,
#' \code{xlab} title for \code{x} axis and \code{ylab} title for \code{y} axis.
#' } 
#' or a two-argument functional data object, see \code{\link{bifd}}.
#' @param type 1-character string giving the type of plot desired.\cr The
#' following values are possible for \code{fdata} class object: "l" for lines
#' (by default),"p" for points, , "o" for overplotted points and lines, "b",
#' "c" for (empty if "c") points joined by lines, "s" and "S" for stair steps
#' and "h" for histogram-like vertical lines. Finally, "n" does not produce any
#' points or lines.\cr The following values are possible for \code{fdata2d}
#' class object: "image.contour" (by default) to display three-dimensional data
#' and add the contour lines, "image" to display three-dimensional data,
#' "contour" to display a contour plot, "persp" to display a perspective plots
#' of a surface over the x-y plane and "filled.contour" to display a contour
#' plot with the areas between the contours filled in solid color.
#' @param main an overall title for the plot: see \code{\link{title}}.
#' @param xlab xlab title for x axis, as in plot.
#' @param ylab ylab title for y axis, as in plot.
#' @param lty a vector of line types, see \code{\link{par}}.
#' @param mfrow A vector of the form c(nr, nc). Subsequent figures will be
#' drawn in an nr-by-nc array on the device by rows (mfrow).
#' @param time The time interval to suspend plot execution for, in seconds, see
#' \link[base]{Sys.sleep}.
#' @param rownames Row names.
#' @param argvals.s a vector of argument values for the first argument s of the
#' functional data object to be evaluated.
#' @param argvals.t a vector of argument values for the second argument t of
#' the functional data object to be evaluated.
#' @param \dots Further arguments passed to \link[graphics]{matplot} function
#' (for fdata class) or \link[graphics]{image}, \link[graphics]{contour},
#' \link[graphics]{persp} or \link[graphics]{filled.contour} (for fdata2d
#' class).
#' @param trim The alpha of the trimming.
#' @param levgray A vector of desired gray levels between 0 and 1; zero
#' indicates "black" and one indicates "white".
#' @author Manuel Febrero Bande and Manuel Oviedo de la Fuente
#' <manuel.oviedo@@usc.es>
#' @seealso See Also as \code{\link{fdata}}
#' @keywords hplot
#' @examples
#' \dontrun{
#' # Example for fdata class of 1 dimension (curve)
#' a1<-seq(0,1,by=.01)
#' a2=rnorm(length(a1),sd=0.2)
#' f1<-(sin(2*pi*a1))+rnorm(length(a1),sd=0.2)
#' nc<-10
#' np<-length(f1)
#' tt=seq(0,1,len=101)
#' mdata<-matrix(NA,ncol=np,nrow=nc)
#' for (i in 1:nc) mdata[i,]<- (sin(2*pi*a1))+rnorm(length(a1),sd=0.2)
#' fdataobj<-fdata(mdata,tt)
#' res=plot.fdata(fdataobj,type="l",col=gray(1:nrow(mdata)/nrow(mdata)))
#' lines(func.mean(fdataobj),col=3,lwd=2) #original curve
#' 
#' # example for fdata2d class of 2 dimension (surface)
#' t1 <- seq(0, 1, length= 51)
#' t2 <- seq(0, 1, length= 31)
#' z<-array(NA,dim=c(4,51,31))
#' for (i in 1:4) z[i,,] <- outer(t1, t2, function(a, b) (i*a)*(b)^i)
#' z.fdata<-fdata(z,list(t1,t2))
#' plot(z.fdata,time=2)
#' plot(z.fdata,mfrow=c(2,2),type="persp",theta=30)
#' }
#' 
#' @rdname plot.fdata
#' @name plot.fdata
#' @method plot fdata
#' @export
#' @export plot.fdata
plot.fdata<-function(x,type,main,xlab,ylab,lty=1,mfrow=c(1,1),time=1,...) {
#if (missing(lty)) lty = par.fda.usc$lty
if (any(class(x)=="fdata2d"))  {
#stop("Object is not fdata2d class")
if (missing(type)) type="image.contour"
#if (missing(main)) main=x[["names"]][["main"]]
#if (missing(xlab)) xlab=x[["names"]][["xlab"]]
#if (missing(ylab)) ylab=x[["names"]][["ylab"]]
dm<-dim(x$data)
j<-1
len.dm<-length(dm)
rng<-range(x$data)
#if (!ask) {
#   if (dm[1]>9)  ask=TRUE
#   else {
#       ask=FALSE
#       if (dm[1]==1) par(mfrow=c(1,1),ask=FALSE) 
#       if (dm[1]==2) par(mfrow=c(1,2),ask=FALSE)        
#       if (dm[1]==3) par(mfrow=c(1,3),ask=FALSE)        
#       if (dm[1]==4) par(mfrow=c(2,2),ask=FALSE)        
#       if (dm[1]>4) par(mfrow=c(2,3),ask=FALSE)                             
#       if (dm[1]>6) par(mfrow=c(3,3),ask=FALSE)                                    
#   }}
par(mfrow=mfrow)
#npar<-par()$mfrow[1]*par()$mfrow[2]
npar<-mfrow[1]*mfrow[2]
for (i in 1:dm[1]) {
#if (ask) {   par(mfrow=c(1,1),ask=ask)     }
z <- x[["data"]][i,,]
if (len.dm==3) {
  if (missing(main))  main <- paste(x$names$main," ",dimnames(x$data)[[1]][i],sep="")
  if (missing(xlab))  xlab <- x$names$xlab[1]
  if (missing(ylab)) ylab =  x$names$xlab[2]
  
  switch (type,
  "persp"={                          
  par(bg = "white")
  xx <- x[["argvals"]][[1]]
  y <- x[["argvals"]][[2]]
  
  
  #nrz <- nrow(z);
  #ncz <- ncol(z)
  #jet.colors <- colorRampPalette( c("yellow", "red") ) 
  #nbcol <- length(xx)
  #color <- jet.colors(nbcol)
  #zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
  #facetcol <- cut(zfacet, nbcol)
  
   persp(x=xx,y=y,z=z,xlim=x[["rangeval"]][[1]],ylim=x[["rangeval"]][[2]],zlim=rng,
   main = ,xlab = xlab, ylab =  ylab,...) 
  # main =  x$names$main[i],xlab = x$names$xlab[1], ylab =  x$names$xlab[2],...)
  #, col=color[facetcol],...)
  # par(op)
  },
  "filled.contour"={
  filled.contour(x=x[["argvals"]][[1]],y=x[["argvals"]][[2]],z=z,
  xlim=x[["rangeval"]][[1]],ylim=x[["rangeval"]][[2]],
  plot.title=title(main = main, xlab =xlab, ylab =ylab),...)},
  
  "contour"={
    contour(x=x[["argvals"]][[1]],y=x[["argvals"]][[2]],z=z,zlim=rng,
  xlim=x[["rangeval"]][[1]],ylim=x[["rangeval"]][[2]],
  plot.title=title(main =  main, xlab = xlab, ylab = ylab),...)},#labels repetidos
      
  "image"={
    image(x = x[["argvals"]][[1]],y = x[["argvals"]][[2]],z=z,
  xlim = x[["rangeval"]][[1]],ylim = x[["rangeval"]][[2]],zlim=rng,
  main = main, xlab = xlab, ylab = ylab,...)},
  "image.contour"={
    image(x = x[["argvals"]][[1]],y = x[["argvals"]][[2]],z=z,
  xlim = x[["rangeval"]][[1]],ylim = x[["rangeval"]][[2]],zlim=rng,
  main = main,    xlab = xlab, ylab =  ylab,...)
      contour(x=x[["argvals"]][[1]],y=x[["argvals"]][[2]],z=z,add = TRUE, drawlabels = FALSE,...)
      }#lattice plot    
  #"contourplot"={contourplot(data=z,
  #xlim = x[["rangeval"]][[1]],ylim = x[["rangeval"]][[2]],...)}    
  )
}

if (names(dev.cur())!="pdf" & j==npar) {
   Sys.sleep(time)
   j<-1
   }
else j<-j+1   
}  
#else {
#  if (len.dm==3) { contourplot(data=z,
#xlim = x[["rangeval"]][[1]],ylim = x[["rangeval"]][[2]],...)}
#    }
if (len.dm>3) stop("Not implemented plot for arrays of more than 3 dimension yet")
}
else {
if (!is.fdata(x))  stop("Object is not fdata class")
if (missing(type)) type="l"
if (missing(main)) main=x[["names"]][["main"]]
if (missing(xlab)) xlab=x[["names"]][["xlab"]]
if (missing(ylab)) ylab=x[["names"]][["ylab"]]

if (is.vector(x[["data"]])) fda::matplot(x[["argvals"]],(x[["data"]]),type=type,lty=lty,main=main,ylab=ylab,xlab=xlab,...)
else fda::matplot(x[["argvals"]],t(x[["data"]]),type=type,lty=lty,main=main,ylab=ylab,xlab=xlab,...)
}
}

#' @rdname plot.fdata
#' @export
lines.fdata = function (x, ...) 
{
  plot(x, add = TRUE, ...)
}

#lines.fdata=function(x,...,lty=1){
# if (missing(lty)) lty = par.fda.usc$lty
# ,lty=lty
# plot(x,add=TRUE,...,lty=lty)  }

#' @rdname plot.fdata
#' @export
title.fdata<-function(x,main=NULL,xlab=NULL,ylab=NULL,rownames=NULL) {
if (!is.fdata(x))  stop("Object is not fdata class")
if (!is.null(rownames)) rownames(x[["data"]])<-rownames
if (!is.null(main)) x[["names"]][["main"]]<-main
if (!is.null(xlab)) x[["names"]][["xlab"]]<-xlab
if (!is.null(ylab)) x[["names"]][["ylab"]]<-ylab
x
}


#' @rdname plot.fdata
#' @export
plot.mdepth<-function(x, trim,  levgray=.9,...){
  dep<-x
  if (missing(trim)) trim<-0.1
  x <- dep$x
  y <- dep$xx
  d<-ncol(x)
  if (is(dep,"mdepth")){
    name=dep$name
    mtrim=dep$mtrim
    med=dep$median
    dep=dep$dep
    tr=rownames(mtrim)
  }else{
    name=dep$name
    dep=dep$dep
    k = which.max(dep)
    med = matrix(x[k, ],ncol=d)
    rownames(med)="Median"
    nl=length(trim)
    mtrim=matrix(NA,ncol=d,nrow=nl)
    for (i in 1:nl){
      lista = which(dep >= quantile(dep, probs = trim[i], na.rm = TRUE))
      if (length(lista) == 1) {        mtrim[i,] <- x[lista, ]    }
      else mtrim[i,] = apply(x[lista, ], 2, mean, na.rm = TRUE)
      tr <- paste(name,".trim", trim * 100, "%", sep = "")
      rownames(mtrim)=tr
    }}
  dd = 0
  if (d == 2) 
    dd = 2
  if (d > 2) 
    dd = 3   
  if (d > 5)   
    dd = 4
  
  if (dd != 0) {
    ind1 <- !is.na(dep)
    color=colorRampPalette(c("red","blue"))(nrow(mtrim)+1)
    cgray = 1 - (dep - min(dep, na.rm = TRUE))/(max(dep, 
                                                    na.rm = TRUE) - min(dep, na.rm = TRUE))
    if (is.data.frame(x)) 
      nam <- names(x)
    if (is.matrix(x)) 
      nam <- colnames(x)
    if (dd == 2) {
      plot(y,col=gray(levgray),main = paste0(name," Depth"),xlab = nam[1], ylab = nam[2],...)
      points(x[ind1,], col = gray(levgray*cgray[ind1]),pch=16)
      points(med[1,1], med[1,2], col = color[1], lwd = 2, pch = 16)
      points(mtrim[,1], mtrim[,2], col = color[-1], pch = 17)
      legend("topleft", legend = c(rownames(med),rownames(mtrim)), pch = c(16,rep(17,nrow(mtrim))), 
             box.col = 0, col = color)
    }
    if (dd == 3) {
      cole=c(rep(gray(levgray),nrow(y)),gray(levgray*cgray[ind1]),color)
      mext=rbind(y,x[ind1,],med,mtrim)
      pairs(mext,col=cole,pch=c(rep(1,nrow(y)),rep(16,nrow(x[ind1,])),rep(19,nrow(mtrim)+1)), 
            main = paste0(name," Depth"))		
    }
    if (dd == 4){
      cole=gray(levgray*cgray[ind1])
      cole[k]=color[1]
      stars(x[ind1, ], col.stars = cole)}
  }
}

#' @rdname plot.fdata
#' @export
plot.depth<-function(x,trim, levgray=.9,...){
  dep <- x
  if (missing(trim)) 
    trim<-dep$trim
  x <- dep$fdataobj
  y <- dep$fdataori
  
  if (class(dep)=="depth"){
    name=dep$name
    mtrim=dep$mtrim
    nl=nrow(mtrim)
    med=dep$median
    dep=dep$dep
  } else{
    dep=dep$dep
    if (length(dep)!=nrow(x)) stop("The number of rows in x is not of length(dep)")
    tt=argvals(x)
    rtt=rangeval(x)
    rot=x$names
    k = which.max(dep)
    med = x[k]
    nl = length(trim)
    mtrim = fdata(matrix(NA, nrow = nl, ncol = ncol(x)),tt,rtt,rot)
    for (j in 1:length(trim)) {
      lista = which(dep >= quantile(dep, probs = trim[j], na.rm = TRUE))
      if (length(lista)==1) {
        mtrim$data[j,]<-x[lista]$data
        if (draw) {draw=FALSE;warning("Too few curves in mtrim. The plot is not drawn")}
      }
      else mtrim$data[j,]=func.mean(x[lista])$data       
    }
    rownames(med$data) <- paste0(name,".med")
    tr <- paste(name,".tr", round(trim * 100,2), "%", sep = "")
    rownames(mtrim$data) <- tr
  }
  ans <- dep
  ind1 <- !is.na(ans)
  color = colorRampPalette(c("red","blue"))( nrow(mtrim$data) + 1)
  cgray = 1 - (ans - min(ans, na.rm = TRUE))/(max(ans, 
                                                  na.rm = TRUE) - min(ans, na.rm = TRUE))
  plot(y, col = gray(levgray), main = paste0(name," Depth"), lty=3, lwd=1,...)
  lines(x[ind1], col = gray(levgray*cgray[ind1]), lty=1,lwd=1.5)
  lines(mtrim, lwd = 3, col = color[-1],lty=1)
  lines(med, col = color[1], lwd = 3)
  legend("topleft", legend = c(rownames(med$data),rownames(mtrim$data)), lwd = 3,box.col=0,col = color)
  
}

############################################ Auxiliar functions
image.scale <- function(z, zlim, col = heat.colors(12),
                        breaks, horiz=TRUE, ylim=NULL, xlim=NULL, ...){
  if(!missing(breaks)){
    if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
  }
  if(missing(breaks) & !missing(zlim)){
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
  }
  if(missing(breaks) & missing(zlim)){
    zlim <- range(z, na.rm=TRUE)
    zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
    zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }
  poly <- vector(mode="list", length(col))
  for(i in seq(poly)){
    poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
  }
  xaxt <- ifelse(horiz, "s", "n")
  yaxt <- ifelse(horiz, "n", "s")
  if(horiz){YLIM<-c(0,1); XLIM<-range(breaks)}
  if(!horiz){YLIM<-range(breaks); XLIM<-c(0,1)}
  if(missing(xlim)) xlim=XLIM
  if(missing(ylim)) ylim=YLIM
  plot(1,1,t="n",ylim=ylim, xlim=xlim, xaxt=xaxt, yaxt=yaxt, xaxs="i", yaxs="i", ...)  
  for(i in seq(poly)){
    if(horiz){
      polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
    }
    if(!horiz){
      polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
    }
  }
}

color.bar <- function(colores, min, max=-min, nticks=length(colores)+1,
                      ro=1,ticks=round(seq(min, max, len=nticks),ro),xm=NULL,ym=NULL,xw=NULL,yw=NULL,horiz=TRUE) {
  xaxp=par()$usr[1:2]
  yaxp=par()$usr[3:4]
  nc=length(colores)
  
  if (horiz) {
    if (is.null(xw)) {xw=0.4*(xaxp[2]-xaxp[1])}
    if (is.null(yw)) {yw=0.025*(yaxp[2]-yaxp[1])}
    if (is.null(xm)) {xm=(xaxp[2]+xaxp[1])/2}
    if (is.null(ym)) {ym=yaxp[2]-yw}
    scale = nc/(2*xw)
    for (i in 1:nc) {
      x = xm-xw+(i-1)/scale 
      rect(x,ym-yw,x+1/scale,ym+yw, col=colores[i], border=NA)	
      text(x,ym-yw,ticks[i],pos=1,offset=0.5)
    }
    text(xm+xw,ym-yw,ticks[nc+1],pos=1,offset=0.5)
  } else {
    if (is.null(xw)) {xw=0.025*(xaxp[2]-xaxp[1])}
    if (is.null(yw)) {yw=0.4*(yaxp[2]-yaxp[1])}
    if (is.null(xm)) {xm=xaxp[2]-xw}
    if (is.null(ym)) {ym=(yaxp[2]+yaxp[1])/2}
    scale = nc/(2*yw)
    for (i in 1:nc) {
      y = ym-yw+(i-1)/scale 
      rect(xm-xw,y,xm+xw,y+1/scale, col=colores[i], border=NA)
      text(xm+xw,y,ticks[i],pos=4,offset=.5)
    }
    text(xm+xw,ym+yw,ticks[nc+1],pos=4,offset=.5)
    
  }
}

        
#' @rdname plot.fdata
#' @export
plot.bifd<-function(x,argvals.s,argvals.t,...){
  if (missing(argvals.s)){
    nfine.s = max(c(201,10*x$sbasis$nbasis+1))
    argvals.s = seq(x$sbasis$rangeval[1],x$sbasis$rangeval[2],len=nfine.s)
  }
  if (missing(argvals.t)){
    nfine.t = max(c(201,10*x$tbasis$nbasis+1))
    argvals.t = seq(x$tbasis$rangeval[1],x$tbasis$rangeval[2],len=nfine.t)
  }
  tt <- list(argvals.s,argvals.t)
  rtt <- list(x$sbasis$rangeval,x$tbasis$rangeval)
  plot(fdata(eval.bifd(argvals.s,argvals.t,x),tt,rtt,fdata2d=TRUE),... )
}



# plot.lfdata<-function(lfdata,ask=FALSE,color,...){
#   mf=5
#   nvar<-length(lfdata)
#   if (nvar>4) ask=TRUE
#   if (ask) {par(mfrow = c(1, 1))
#             dev.interactive()
#             oask <- devAskNewPage(TRUE)
#             on.exit(devAskNewPage(oask))}
#   else{    mf<-switch(nvar,
#                       "1"={c(1,1)},
#                       "2"={c(1,2)},
#                       "3"={c(1,3)},
#                       "4"={c(2,2)})            
#            par(mfrow =mf)                    }
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
#     plot.fdata(lfdata[[idat]], col =  color2, main =nam[idat],...)
#   }
# }