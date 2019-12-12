#################################################
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
#################################################
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

