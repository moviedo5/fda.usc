#' @title Bootstrap regression
#' @rdname fregre.bootstrap
#' @description Estimate the beta parameter by wild or smoothed bootstrap procedure
#' @aliases fregre.bootstrap fregre.bootstrap2
#' @details Estimate the beta parameter by wild or smoothed bootstrap procedure using
#' principal components representation \code{\link{fregre.pc}}, Partial least
#' squares components (PLS) representation \code{\link{fregre.pls}} or basis
#' representation \code{\link{fregre.basis}}.\cr If a new curves are in
#' \code{newX} argument the bootstrap method estimates the response using the
#' bootstrap resamples.
#' 
#' If the model exhibits heteroskedasticity, the use of wild bootstrap
#' procedure is recommended (by default).
#' 
#' @param model \code{fregre.pc}, \code{fregre.pls} or \code{fregre.basis}
#' object.
#' @param nb Number of bootstrap samples.
#' @param wild Naive or smoothed bootstrap depending of the \code{smo} and
#' \code{smoX} parameters.
#' @param type.wild Type of distribution of V in wild bootstrap procedure, see
#' \code{\link{rwild}}.
#' @param smo If \eqn{>0}, smoothed bootstrap on the residuals (proportion of
#' response variance).
#' @param smoX If \eqn{>0}, smoothed bootstrap on the explanatory functional
#' variable \code{fdata} (proportion of variance-covariance matrix of
#' \code{fdata} object.
#' @param newX A \code{fdata} class containing the values of the model
#' covariates at which predictions are required (only for smoothed bootstrap).
#' @param kmax.fix The number of maximum components to consider in each
#' bootstrap iteration.  =TRUE, the bootstrap procedure considers the same
#' number of components used in the previous fitted model.  =FALSE, the
#' bootstrap procedure estimates the best components in each iteration.
#' @param alpha Significance level used for graphical option, \code{draw=TRUE}.
#' @param draw =TRUE, plot the bootstrap estimated beta, and (optional) the CI
#' for the predicted response values.
#' @param \dots Further arguments passed to or from other methods.
#' @return Return:
#' \itemize{
#' \item \code{model}{ \code{fregre.pc}, \code{fregre.pls} or \code{fregre.basis} object.} 
#' \item \code{beta.boot}{ functional beta estimated by the \code{nb} bootstrap regressions.} 
#' \item \code{norm.boot}{ norm of diferences beetween the nboot betas estimated by bootstrap and beta estimated by regression model.} 
#' \item \code{coefs.boot}{ matrix with the bootstrap estimated basis coefficients.} 
#' \item \code{kn.boot}{ vector or list of length \code{nb} with index of the basis, PC or PLS factors selected in each bootstrap
#' regression.} 
#' \item \code{y.pred}{ predicted response values using \code{newX} covariates.} 
#' \item \code{y.boot}{ matrix of bootstrap predicted response values using \code{newX} covariates.} 
#' \item \code{newX}{ a \code{fdata} class containing the values of the model covariates at which predictions are required (only
#' for smoothed bootstrap).}
#' }
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente \email{manuel.oviedo@@udc.es}
#' @seealso See Also as: \code{\link{fregre.pc}}, \code{\link{fregre.pls}},
#' \code{\link{fregre.basis}}, .
#' @references Febrero-Bande, M., Galeano, P. and Gonzalez-Manteiga, W. (2010).
#' \emph{Measures of influence for the functional linear model with scalar
#' response}. Journal of Multivariate Analysis 101, 327-339.
#' 
#' Febrero-Bande, M., Oviedo de la Fuente, M. (2012).  \emph{Statistical
#' Computing in Functional Data Analysis: The R Package fda.usc.} Journal of
#' Statistical Software, 51(4), 1-28. \url{https://www.jstatsoft.org/v51/i04/}
#' @keywords regression
#' @examples
#' \dontrun{ 
#' data(tecator)
#' iest<-1:165
#' x=tecator$absorp.fdata[iest]
#' y=tecator$y$Fat[iest]
#' nb<-25  ## Time-consuming
#' res.pc=fregre.pc(x,y,1:6)
#' # Fix the compontents used in the each regression
#' res.boot1=fregre.bootstrap(res.pc,nb=nb,wild=FALSE,kmax.fix=TRUE)
#' # Select the "best" compontents used in the each regression
#' res.boot2=fregre.bootstrap(res.pc,nb=nb,wild=FALSE,kmax.fix=FALSE) 
#' res.boot3=fregre.bootstrap(res.pc,nb=nb,wild=FALSE,kmax.fix=10) 
#' ## predicted responses and bootstrap confidence interval
#' newx=tecator$absorp.fdata[-iest]
#' res.boot4=fregre.bootstrap(res.pc,nb=nb,wild=FALSE,newX=newx,draw=TRUE)
#' }
#' @export
fregre.bootstrap <- function(model, nb=500, wild=TRUE, type.wild="golden"
                           , newX=NULL, smo=0.1, smoX=0.05, alpha=0.95
                           , kmax.fix=FALSE,draw=TRUE,...){

   fdataobj=model$fdataobj
nas<-is.na.fdata(fdataobj)
if (any(nas))  {
   fdataobj$data<-fdataobj$data[!nas,]
   cat("Warning: ",sum(nas)," curves with NA are not used in the calculations \n")
   }
dat<-fdataobj[["data"]]
tt<-fdataobj[["argvals"]]
rtt<-fdataobj[["rangeval"]]
nam<-fdataobj[["names"]]
resi=model$residuals
if (model$call[[1]]=="fregre.pc") {
     beta.est=model$beta.est$data
     pc=1
     }
else if (model$call[[1]]=="fregre.pls") {
          beta.est=model$beta.est
          pc=2
          }
else if (model$call[[1]]=="fregre.basis" || model$call[[1]]=="fregre.basis.cv") {
          beta.est=model$beta.est
          beta.est=eval.fd(tt,beta.est)
          pc=3
          }
else stop("No fregre.pc, fregre.basis or fregre.basis.cv object in model argument")
a.est=model$coefficients[1]
sr2=model$sr2
n <- nrow(fdataobj)
J <- ncol(fdataobj)
#alpha<-1-alpha
ncoefs<-100

cb.num <- round(alpha * nb)
yp<-NULL
if (!is.null(newX)) { 
yp<-predict(model,newX)
}
#  if (!is.logical(kmax.fix)) {
#  criteria=kmax.fix
#  kmax.fix=TRUE   
#   } 
#  else   {
knn.fix=list()
 if (pc<3)   {
      criteria="SIC"
      if (kmax.fix==TRUE) knn.fix<-model$l
      maxl<-max(model$l)
      if (is.numeric(kmax.fix)) {
            maxl<-max(maxl,kmax.fix)
            kmax.fix=FALSE
            }
      }
  if (pc==3)    { 
   criteria=GCV.S
   maxl<-model$basis.b.opt$nbasis
   maxx<-model$basis.x.opt$nbasis
   if (kmax.fix==TRUE) knn.fix<-c(model$basis.x.opt$nbasis,model$basis.b.opt$nbasis)  
   if (is.numeric(kmax.fix)) {
    maxl<-max(maxl,kmax.fix)  
    maxx<-max(maxx,kmax.fix)    
    kmax.fix=FALSE
   }
    ncoefs<-nrow(model$beta.est$coefs)
  }
#  }
#  if (kmax.fix) coefs.boot <- array(NA,dim=c(nb,ncoefs))
#  else coefs.boot<-list()   

par.fda.usc <- eval(parse(text="fda.usc:::par.fda.usc"), envir=.GlobalEnv)
ncores <- par.fda.usc$ncores
inumgr <- icount(nb)

betas.boot <- array(NA,dim=c(nb,J))

#betas.boot2<-model$beta.est
#norm.boot <- array(NA,dim=c(nb,1))
y.mue2<-array(NA,dim=c(nb,nrow(dat)))
#ypred<-array(NA,dim=c(nb,nrow(newX)))

#norm.boot <- NULL

comb <- function(...) {
   mapply('rbind', ..., SIMPLIFY=FALSE)
}
#result <- foreach(i=1:100, .combine='comb', .multicombine=TRUE) %dopar% {
varX <- var(dat)

if (!wild){
   #for (i in 1:nb){
   #list("betas.boot", "norm.boot", "coefs.boot","model") 
out <- foreach(i = inumgr, .packages='fda.usc', .combine='comb', .multicombine=TRUE) %dopar% {
   normboot <- NULL
   betasboot <- NULL
   normboot <- NULL
   coefsboot <- NULL
   ypred <- NULL
   knnfix <- integer(length(model$coefficients[-1]))
   muee <- sample(1:n,n,replace=TRUE)
   mueX <- sample(1:n,n,replace=TRUE)
   residuals.mue <- resi[muee] + rnorm(n,0,sqrt(smo * sr2))  
   b1 <- fdata(mvrnorm(n,rep(0,J),smoX * varX ), argvals(fdataobj), rtt)
   b0 <- fdataobj[mueX,]
   fdata.mue <- b0+b1
   if (pc==1)   {
      y.mue <- predict.fregre.fd(model,fdata.mue) + residuals.mue
      knnfix <- model$l
      if (kmax.fix)    funcregpc.mue <- fregre.pc(fdata.mue,y.mue,l=model$l,lambda=model$lambda,P=model$P,weights=model$weights)    
       else {
             fpc <- fregre.pc.cv(fdata.mue,y.mue,kmax=maxl,lambda=model$lambda,P=model$P,criteria=criteria,weights=model$weights)
             knnfix <- fpc$pc.opt
             funcregpc.mue <- fpc$fregre.pc
             }
       betasboot <- funcregpc.mue$beta.est$data
       bb <- model$beta.est-funcregpc.mue$beta.est
       normboot <- norm.fdata(bb)
    }
    else  if (pc==2)  {
      y.mue<-predict.fregre.fd(model,fdata.mue)  + residuals.mue
      knnfix <- model$l
      if (kmax.fix)    {funcregpc.mue <- fregre.pls(fdata.mue, y.mue,model$l)}
       else     {      
                        fpc <- fregre.pls.cv(fdata.mue, y.mue,maxl, criteria=criteria)
                        knnfix[1:length(fpc$pls.opt)] <- fpc$pls.opt
                        funcregpc.mue<-fpc$fregre.pls
                        }
       betasboot <- funcregpc.mue$beta.est$data
       bb<-model$beta.est-funcregpc.mue$beta.est
       normboot <- norm.fdata(bb)
             }
    else  {
       bett<-fdata(t(beta.est),tt,rtt)
        y.mue<-predict.fregre.fd(model,fdata.mue)  + residuals.mue
        knnfix <- 1:length(model$coefficient[-1])
       if (kmax.fix) funcregpc.mue <-fregre.basis(fdata.mue,y.mue,model$basis.x.opt,
        model$basis.b.opt,Lfdobj = model$Lfdobj,weights=model$weights)
       else {
           funcregpc.mue <-fregre.basis.cv(fdata.mue,y.mue,basis.x=maxx,basis.b=maxl,type.CV=criteria,Lfdobj=model$Lfdobj,weights=model$weights)                   
              knnfix<-c(funcregpc.mue$basis.x.opt$nbasis,funcregpc.mue$basis.b.opt$nbasis)
              }
       betasbooT <- eval.fd(tt,funcregpc.mue$beta.est)
       bb<-model$beta.est-funcregpc.mue$beta.est
       normboot<-  norm.fd(bb)
    } 
   if (kmax.fix)  coefsboot <-funcregpc.mue$coefficients
   else         coefsboot<-list(funcregpc.mue$coefficients)
   if (!is.null(newX))   {      ypred<-predict(funcregpc.mue,newX)}
 list("betas.boot"= betasboot, "norm.boot"=normboot, "coefs.boot"=coefsboot
      , "ypred"=ypred,"knn.fix"=knnfix) 
  }    
}
else { # wild = TRUE
 pred <- model$fitted.values
 fdata.mue <- fdataobj
 #for (i in 1:nb){
 out <- foreach(i = inumgr, .packages='fda.usc', .combine='comb', .multicombine=TRUE) %dopar% {
    normboot <- NULL
    betasboot <- NULL
    normboot <- NULL
    coefsboot <- NULL
    ypred <- NULL
    knnfix <- integer(length(model$coefficients[-1]))
   muee <- sample(1:n,n,replace=TRUE)
   residuals.mue <- rwild(resi[muee],type.wild)
   fdata.mue <- fdataobj[muee] 
   if (pc==1)   {
      y.mue<- pred[muee]  + residuals.mue
      knnfix <- model$l
      if (kmax.fix)    funcregpc.mue <- fregre.pc(fdata.mue,y.mue,l=model$l,lambda=model$lambda,P=model$P,weights=model$weights)
       else     {
               fpc <- fregre.pc.cv(fdata.mue,y.mue,kmax=maxl,lambda=model$lambda,P=model$P
                                   ,criteria=criteria,weights=model$weights)
               knnfix[1:length(fpc$pcs.opt)] <- fpc$pc.opt
               funcregpc.mue<-fpc$fregre.pc
                }
       betasboot <- funcregpc.mue$beta.est$data
       bb<-model$beta.est-funcregpc.mue$beta.est
       normboot <- norm.fdata(bb)
       }
      else  if (pc==2)  {
       y.mue<-pred[muee]  + residuals.mue
       knnfix <- model$l
       if (kmax.fix)    {funcregpc.mue <- fregre.pls(fdata.mue,y.mue,model$l)}
       else     {
                        fpc <- fregre.pls.cv(fdata.mue,y.mue,maxl,criteria=criteria)              
                        knnfix<-fpc$pls.opt
                        funcregpc.mue<-fpc$fregre.pls
                        }
       betasboot <- funcregpc.mue$beta.est$data
       bb<- model$beta.est - funcregpc.mue$beta.est
       normboot <- norm.fdata(bb)
      }
    else  {
       bett <- fdata(t(beta.est),tt,rtt)
       y.mue <- pred[muee]  + residuals.mue
       knnfix <- 1:length(model$coefficient[-1]) 
       if (kmax.fix) funcregpc.mue <-fregre.basis(fdata.mue,y.mue,model$basis.x.opt,
        model$basis.b.opt,Lfdobj=model$Lfdobj,weights=model$weights)
       else {
           funcregpc.mue <-fregre.basis.cv(fdata.mue,y.mue,basis.x=maxx,basis.b=maxl,type.CV=criteria
                                           ,Lfdobj=model$Lfdobj,weights=model$weights)                   
           knnfix <-c(funcregpc.mue$basis.x.opt$nbasis,funcregpc.mue$basis.b.opt$nbasis)
            }
       betasboot <- eval.fd(tt,funcregpc.mue$beta.est)
       bb <- model$beta.est - funcregpc.mue$beta.est
       normboot <-  norm.fd(bb)
    }
   if (kmax.fix)  coefsboot <-funcregpc.mue$coefficients
   else         coefsboot<-list(funcregpc.mue$coefficients)
   if (!is.null(newX))   {      ypred <- predict(funcregpc.mue,newX)}
   list("betas.boot"= betasboot, "norm.boot"=normboot, "coefs.boot"=coefsboot, "ypred"=ypred
        ,"knn.fix"=knnfix) 
   }    
}    
# betas.boot ;norm.boot; coefs.boot; model
#print(lbetas.boot[1:3,1:4])
betas.boot<-out$betas.boot
norm.boot<-out$norm.boot
coefs.boot<-out$coefs.boot
ypred <- out$ypred
knn.fix <- out$knn.fix

#print(knn.fix)

#print(betas.boot[1:3,1:2])
betas.boot<- fdata(betas.boot,tt,rtt,nam)
betas.boot$names$main <- "beta.est bootstrap"
output <- list("model"=model, "betas.boot"=betas.boot, "norm.boot"=norm.boot, "coefs.boot"= coefs.boot,
"kn.boot"=knn.fix,"y.boot"=ypred)
if (draw) {
  out <- norm.boot > quantile(norm.boot,alpha)
  plot(betas.boot[-out],col="grey")
  lines(model$beta.est,col=4)
  lines(betas.boot[out],col=2,lty=2)
if (!is.null(newX))   {
#   dev.new()
prd<- prod(par("mfcol"))
if (prd==1) {
  oask <- devAskNewPage(TRUE)
  on.exit(devAskNewPage(oask))
}
 
 IC<-apply(ypred,2,quantile,c((1-alpha)/2,alpha+(1-alpha)/2))
#   yp<-predict(model,newX)
 m<-ncol(IC)
 ylm<-range(rbind(IC,drop(yp)))
 matplot( rbind(1:m,1:m),IC,type="l",lty=1,col=1,ylim=ylm,
   xlab="Id newX curves",ylab="predicted value",main=paste("y predicted and ",alpha*100,"% bootstrap CI",sep=""))
   points(yp,col=4,pch=16,cex=.7)  

}
}
output[["y.pred"]]=yp
output[["newX"]]=newX
return(output)
}

