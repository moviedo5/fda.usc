#' @title Classifier from Functional Data
#' 
#' @description Classification of functional data using maximum depth.
#' 
#' @param group Factor of length \emph{n}
#' @param fdataobj \code{fdata}, \code{matrix} or \code{data.frame} class
#' object of train data.
#' @param newfdataobj \code{fdata}, \code{matrix} or \code{data.frame} class
#' object of test data.
#' @param depth Type of depth function from functional data:
#' \itemize{ 
#' \item \code{FM}: Fraiman and Muniz depth.
#' \item \code{mode}: Modal depth.
#' \item \code{RT}: Random Tukey depth.
#' \item \code{RP}: Random project depth.
#' \item \code{RPD}: Double random project depth.
#' }
#' @param par.depth List of parameters for \code{depth}.
#' @param CV =``none'' \code{group.est=group.pred}, =TRUE \code{group.est} is
#' estimated by cross-validation, =FALSE \code{group.est} is estimated.
#' 
#' @return 
#' \itemize{ 
#' \item \code{group.est}: Vector of classes of train sample data.
#' \item \code{group.pred}: Vector of classes of test sample data.
#' \item \code{prob.classification}: Probability of correct classification by group.
#' \item \code{max.prob}: Highest probability of correct classification.
#' \item \code{fdataobj}: \code{\link{fdata}} class object.
#' \item \code{group}: Factor of length \emph{n}.
#' }
# \item{prob.group}{ Matrix of predicted class probabilities. For each functional 
# point shows the probability of each possible group membership.}

#' @author Febrero-Bande, M. and Oviedo de la Fuente, M.
#' 
#' @references Cuevas, A., Febrero-Bande, M. and Fraiman, R. (2007).
#' \emph{Robust estimation and classification for functional data via
#' projection-based depth notions.} Computational Statistics 22, 3, 481-496.
#' @keywords classif
#' @examples
#' \dontrun{
#' data(phoneme)
#' mlearn<-phoneme[["learn"]]
#' mtest<-phoneme[["test"]]
#' glearn<-phoneme[["classlearn"]]
#' gtest<-phoneme[["classtest"]]
#' 
#' a1<-classif.depth(glearn,mlearn,depth="RP")
#' table(a1$group.est,glearn)
#' a2<-classif.depth(glearn,mlearn,depth="RP",CV=TRUE)
#' a3<-classif.depth(glearn,mlearn,depth="RP",CV=FALSE)
#' a4<-classif.depth(glearn,mlearn,mtest,"RP")
#' a5<-classif.depth(glearn,mlearn,mtest,"RP",CV=TRUE)     
#' table(a5$group.est,glearn)
#' a6<-classif.depth(glearn,mlearn,mtest,"RP",CV=FALSE)
#' table(a6$group.est,glearn)
#' }
#' 
#' @rdname classif.depth
#' @export
classif.depth<-function(group,fdataobj,newfdataobj,depth="RP",
                        par.depth=list(),CV="none"){
#,control=list(trace=FALSE,draw=TRUE)
 C<-match.call()
 if (!is.factor(group)) group<-factor(group)
 ismissing<-missing(newfdataobj)
 if (ismissing) newfdataobj<-fdataobj
 group<-factor(group,levels=levels(group)[which(table(group)>0)])
 func.clas<-list()
 lev<-levels(group)
 ng<-length(lev)       
 nc<-ncol(newfdataobj)
 nvec<-table(group)
# p<-nvec[1]/n        
 if (depth %in% c("PD","HD","RP","RPD","RT")){
  if (is.null(par.depth$proj)) {
   d <- nc         
   u <- matrix(runif(d*500,-1,1),500,d)
   norm <- sqrt(rowSums(u*u))
   arg <- u/norm
   if (depth %in% c("RP","RPD","RT"))  par.depth$proj<-fdata(arg,fdataobj$argvals,fdataobj$rangeval)
   else   par.depth$proj<-arg
   } 
 }
 if (depth %in% c("mband","mmode","HD","SD","PD","MhD"))  par.depth$x<-newfdataobj
 else    par.depth[["fdataobj"]]<-newfdataobj 
 depth<-paste("depth.",depth,sep="")

 if (ismissing) {
  ismdist<-is.matrix(par.depth$metric)
 if (ismdist) {
   mdist<-par.depth$metric  
 }
 # if (depth %in% c("HD","SD","PD","MhD"))  par.depth$x<-fdataobj
#  else    par.depth[["fdataobj"]]<-fdataobj 
#  print(names(par.depth))
  n<-nrow(fdataobj)
  x<-array(NA,dim=c(n,nc,ng))
  Df<-matrix(NA,ncol=ng,nrow=n)
#  if (CV!=TRUE){
   ind<-matrix(NA,nrow=n,ncol=ng)
   for (i in 1:ng) {
    ind[,i]<-group==lev[i]
    nam<-c(paste("depth ",lev[i],sep=""),paste("depth ",paste(lev[-i],collapse=",")))
    if (depth %in% c("depth.mband","depth.mmode","depth.SD","depth.HD","depth.PD","depth.MhD"))  par.depth$xx<-fdataobj[ind[,i],]   
    else   par.depth$fdataori<-fdataobj[ind[,i],]  
      if (ismdist) {
          par.depth$metric<-mdist[,ind[,i]]   
          par.depth$metric2<-mdist[ind[,i],ind[,i]]   }     
     Df[,i]<-switch(depth,
     depth.HD= do.call(depth,par.depth)$dep,
     depth.SD= do.call(depth,par.depth)$dep,
     depth.PD= do.call(depth,par.depth)$dep,
     depth.MhD= do.call(depth,par.depth)$dep,
     depth.FM=do.call(depth,par.depth)$dep,                                                                                
     depth.mode=do.call(depth,par.depth)$dep,
     depth.mmode=do.call(depth,par.depth)$dep,    
     depth.RPD=do.call(depth,par.depth)$dep,     
     depth.RP=do.call(depth,par.depth)$dep,
     depth.RT=do.call(depth,par.depth)$dep,
     depth.mband=do.call(depth,par.depth)$dep,      
     depth.band=do.call(depth,par.depth)$dep)
   }
   group.pred<-group.est<-factor(lev[apply(Df,1,which.max)],levels=lev) # Maximum depth 
#  }
  if (CV==TRUE) {
   group.est<-group
   for (j in 1:n) {
    xj<-fdataobj[j,]   
    xnoj<-fdataobj[-j,]   
    ind<-matrix(NA,nrow=n-1,ncol=ng)    
    for (i in 1:ng) {
     ind[,i]<-group[-j]==lev[i]
     xnoji<-xnoj[ind[,i],]    
     nam<-c(paste("depth ",lev[i],sep=""),paste("depth ",paste(lev[-i],collapse=",")))
     if (depth %in% c("depth.mband","depth.mmode","depth.SD","depth.HD","depth.PD","depth.MhD"))  par.depth$xx<-xnoji
     else   par.depth$fdataori<-xnoji
      if (ismdist) {
          par.depth$metric<-mdist[,ind[,i]]   
          par.depth$metric2<-mdist[ind[,i],ind[,i]]   }    
     Df[,i]<-switch(depth,
     depth.HD= do.call(depth,par.depth)$dep,
     depth.PD= do.call(depth,par.depth)$dep,
     depth.SD= do.call(depth,par.depth)$dep,     
     depth.MhD= do.call(depth,par.depth)$dep,
     depth.FM=do.call(depth,par.depth)$dep,
     depth.mode=do.call(depth,par.depth)$dep,
     depth.mmode=do.call(depth,par.depth)$dep,  
     depth.RPD=do.call(depth,par.depth)$dep,     
     depth.RP=do.call(depth,par.depth)$dep,
     depth.band=do.call(depth,par.depth)$dep,
     depth.mband=do.call(depth,par.depth)$dep,      
     depth.RT=do.call(depth,par.depth)$dep)
   } 
  group.est[j] <- factor(lev[which.max(Df[j,])],levels=lev) # Maximum depth }
  }  
  }
  prob.classification <- diag(table(group,group.est))/table(group) 
  mis <- mean(group.est!=group)
  output <- list("group.est"=group.est,"group.pred"=group.pred,"dep"=Df,"depth"=depth, "par.depth"=par.depth,
  "group"=group,"fdataobj"=fdataobj,"C"=C,"prob.classification"=prob.classification,"max.prob"=1-mis)
 class(output) <- "classif"
 return(output) 
}
else  {   # new data
 n <- nrow(newfdataobj)
 n0 <- nrow(fdataobj)
 x <- array(NA,dim=c(n,nc,ng))
 Df <- matrix(NA,ncol=ng,nrow=n)
 ind <- matrix(NA,nrow=n0,ncol=ng)
 for (i in 1:ng) {
   ind[,i]<-group==lev[i]
   nam<-c(paste("depth ",lev[i],sep=""),paste("depth ",paste(lev[-i],collapse=",")))
   if (depth %in% c("depth.mband","depth.mmode","depth.SD","depth.HD","depth.PD","depth.MhD"))  par.depth$xx<-fdataobj[ind[,i],]   
   else   par.depth$fdataori<-fdataobj[ind[,i],]  
    Df[,i]<-switch(depth,
    depth.HD= do.call(depth,par.depth)$dep,
    depth.PD= do.call(depth,par.depth)$dep,
    depth.SD= do.call(depth,par.depth)$dep,    
    depth.MhD= do.call(depth,par.depth)$dep,
    depth.FM=do.call(depth,par.depth)$dep,
    depth.mode=do.call(depth,par.depth)$dep,
    depth.mmode=do.call(depth,par.depth)$dep,    
    depth.RP=do.call(depth,par.depth)$dep,
    depth.RPD=do.call(depth,par.depth)$dep,  
    depth.band=do.call(depth,par.depth)$dep,      
    depth.mband=do.call(depth,par.depth)$dep,     
    depth.RT=do.call(depth,par.depth)$dep)
 } 
 group.pred<-factor(lev[apply(Df,1,which.max)],levels=lev) # Maximum depth 
 if (CV!="none"){ 
  if (depth %in% c("mband","mmode","HD","SD","PD","MhD"))  par.depth$x<-fdataobj
  else    par.depth[["fdataobj"]]<-fdataobj 
  n<-nrow(fdataobj)
  x<-array(NA,dim=c(n,nc,ng))
  Df2<-matrix(NA,ncol=ng,nrow=n)
  ind<-matrix(NA,nrow=n,ncol=ng)
  if (!CV) {
   for (i in 1:ng) {
    ind[,i]<-group==lev[i]
    nam<-c(paste("depth ",lev[i],sep=""),paste("depth ",paste(lev[-i],collapse=",")))
    if (depth %in% c("depth.mband","depth.mmode","depth.HD","depth.PD","depth.MhD"))  par.depth$xx<-fdataobj[ind[,i],]   
    else   par.depth$fdataori<-fdataobj[ind[,i],]    
    Df2[,i]<-switch(depth,
    depth.HD= do.call(depth,par.depth)$dep,
    depth.PD= do.call(depth,par.depth)$dep,
    depth.SD= do.call(depth,par.depth)$dep,    
    depth.MhD= do.call(depth,par.depth)$dep,
    depth.FM=do.call(depth,par.depth)$dep,
    depth.mode=do.call(depth,par.depth)$dep,
    depth.mmode=do.call(depth,par.depth)$dep,    
    depth.RP=do.call(depth,par.depth)$dep,
    depth.RPD=do.call(depth,par.depth)$dep,   
    depth.band=do.call(depth,par.depth)$dep,     
    depth.mband=do.call(depth,par.depth)$dep,     
    depth.RT=do.call(depth,par.depth)$dep)
 }
 group.est<-factor(lev[apply(Df2,1,which.max)],levels=lev) # Maximum depth 
 }
 else {
 group.est<-group 
 for (j in 1:n) {                       
   xj<-fdataobj[j,]   
   xnoj<-fdataobj[-j,]  
   ind<-matrix(NA,nrow=n-1,ncol=ng) 
   for (i in 1:ng) {
   ind[,i]<-group[-j]==lev[i]
   xnoji<-xnoj[ind[,i],]    
   nam<-c(paste("depth ",lev[i],sep=""),paste("depth ",paste(lev[-i],collapse=",")))
   if (depth %in% c("depth.mband","depth.mmode","depth.HD","depth.PD","depth.MhD"))  par.depth$xx<-xnoji
   else   par.depth$fdataori<-xnoji

    Df2[,i]<-switch(depth,
    depth.HD= do.call(depth,par.depth)$dep,
    depth.PD= do.call(depth,par.depth)$dep,
    depth.SD= do.call(depth,par.depth)$dep,    
    depth.MhD= do.call(depth,par.depth)$dep,
    depth.FM=do.call(depth,par.depth)$dep,   
    depth.mode=do.call(depth,par.depth)$dep,
    depth.mmode=do.call(depth,par.depth)$dep,    
    depth.RP=do.call(depth,par.depth)$dep,
    depth.RPD=do.call(depth,par.depth)$dep,  
    depth.band=do.call(depth,par.depth)$dep,      
    depth.mband=do.call(depth,par.depth)$dep,          
    depth.RT=do.call(depth,par.depth)$dep)
  } 
  group.est[j]<-factor(lev[which.max(Df2[j,])],levels=lev) # Maximum depth }
  } 
 }
 prob.classification<-diag(table(group,group.est))/table(group) 
 mis<-mean(group.est!=group)
 return(list("group.est"=group.est,"group.pred"=group.pred,"dep"=Df,"dep.ori"=Df2,"depth"=depth,
 "par.depth"=par.depth,"group"=group,"fdataobj"=fdataobj,"C"=C,
 "prob.classification"=prob.classification,"max.prob"=1-mis))
 }
 else return(list("group.pred"=group.pred,"dep"=Df,"depth"=depth,
 "par.depth"=par.depth,"group"=group,"fdataobj"=fdataobj,"C"=C))
   }
 }
################################################################################
################################################################################

