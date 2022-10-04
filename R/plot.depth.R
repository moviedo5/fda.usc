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
#################################################################################
plot.depth<-function(x,trim, levgray=.9,...){
  dep <- x
  if (missing(trim)) 
    trim<-dep$trim
  x <- dep$fdataobj
  y <- dep$fdataori

	if (inherits(dep,"depth")){
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
