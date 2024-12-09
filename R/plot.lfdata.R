#' @export plot.lfdata
plot.lfdata<-function(x,ask=FALSE, col,...){
  mf=5
  nvar<-length(x)
  if (nvar>4) ask=TRUE
  if (ask) {par(mfrow = c(1, 1))
            dev.interactive()
            oask <- devAskNewPage(TRUE)
            on.exit(devAskNewPage(oask))}
  else{    mf<-switch(nvar,
                      "1"={c(1,1)},
                      "2"={c(1,2)},
                      "3"={c(1,3)},
                      "4"={c(2,2)})            
           par(mfrow =mf)                    }
  names1<-names2<-names<-x[[1]][["names"]]
  names1$main<-"Multivariate Functional Data"
  #    tr<-paste("mode.tr",trim*100,"\u0025",sep="")         
  nam<-names(x)
  
  if (is.null(nam)) nam<-1:nvar
  for (idat in 1:nvar) {
    data<-x[[idat]]$data
    tt<-x[[idat]]$argvals
    rtt<-x[[idat]]$rangeval
    if (missing(col)) col2<-1
    else {
      if (is.list(col)) col2<-col[[idat]]
      else col2<-col
    }
    plot.fdata(x[[idat]], col =  col2, main =nam[idat],...)
  }
  
}