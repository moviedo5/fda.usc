
corSigma<-function(par.gls,fixed=FALSE){
 type.cor<-class(par.gls$correlation)[1]
 ff<-as.formula("e~1") 
#print("entra corSigma")
# viene del optim.np 
 isfdata<- class(par.gls$data)=="list"
 if (isfdata) { 
#  print("si es fdata")
  grupos<-TRUE
  covariate<-"covariate"
  fact<-"fact"
#  print(covariate)
#  print(names(par.gls$data))  
#  print(names(par.gls$data$data))    
#  ff2<-formula(paste("~",covariate,"|",fact,sep=""))
  ff2<-formula("~covariate|fact",sep="")
  par.gls$data<-data.frame(par.gls$data$data,e=par.gls$data$e)

 # print(dim(par.gls$data))
#  print(par.gls$data[1:33,])
 }
 else{
# print("no fdata")
  fact<-nlme::getGroupsFormula(par.gls$correlation)[[2]]
  if (is.null(fact)) grupos<-FALSE
  else grupos<-TRUE
  covariate<-nlme::getCovariateFormula(par.gls$correlation)[[2]]
  if (grupos) ff2<-formula(paste("~",covariate,"|",fact,sep=""))
  else ff2<-formula(paste("~",covariate,sep=""))
 }     
 par.gls$model<-ff
 correl<-list(form=ff2,fixed=fixed)
 if (type.cor=="corARMA") {
      correl$p<-attributes(nlme::Initialize(par.gls$correlation,par.gls$data))$p
      correl$q<-attributes(nlme::Initialize(par.gls$correlation,par.gls$data))$q
      correl$value = double(correl$p + correl$q)
 }
#print(correl)
par.gls$correlation<-do.call(type.cor,correl)
# print("ok")
# print(names(par.gls))
out<-do.call("gls",par.gls)
# print("ok2")
# print(out)
out$formula<-par.gls$formula
cor1<-coef(out$modelStruct[[1]],FALSE)
fact1<-nlme::getGroups(out,par.gls$data)
lcor<-list()
lcor$value<-cor1
# print("peta glse")
if (grupos) lcor$form<-as.formula(paste("~",covariate,"|",fact,sep=""))
else lcor$form<-as.formula(paste("~",covariate,sep=""))
if (type.cor=="corARMA") {
      lcor$p<-attributes(nlme::Initialize(par.gls$correlation,par.gls$data))$p
      lcor$q<-attributes(nlme::Initialize(par.gls$correlation,par.gls$data))$q
 }
lcor$fixed<-TRUE
cs1<-do.call(type.cor,lcor)
cs2 <- nlme::Initialize(cs1, data = par.gls$data)
Sigma<-nlme::corMatrix(cs2)
#print("peta glse 2")
#print(length(Sigma))
# print(dim(Sigma[[1]]))
n<-nrow(par.gls$data)

if (grupos ) {
 lev<-levels(fact1)
 lev2<-names(Sigma)
 nlev<-nlevels(fact1)
 W0<-Sigma[[lev2[1]]]
 W <- try(solve(W0),silent=TRUE)
# print(dim(W))
# print(dim(W0)) 
     if (is(W,"try-error")) {
      sv<-svd(W0)
      W<-drop((sv$v%*%diag(1/sv$d)%*%t(sv$u)))
      warning("Inverse of sigma computed by SVD")
     }
#print(n)
#print(nlev)
 if (!isfdata){
  W0<-matrix(0,n,n)
  for  (ilev in 1:nlev)  {   
#    print(ilev)
    W0[fact1==lev2[ilev],fact1==lev2[ilev]]<-W
  }
  W<-W0
 }  
 else {
   W0<-Sigma[[1]] 
   }
}
else{ 
 W0<-Sigma
 W <- try(solve(W0),silent=TRUE)
     if (is(W,"try-error")) {
      sv<-svd(W0)
      W<-drop((sv$v%*%diag(1/sv$d)%*%t(sv$u)))
      warning("Inverse of sigma computed by SVD")
     }
}
#print("sale corSigma ")
out.gls<-out
out.gls$Sigma<-W0
out.gls$W<-W
out.gls$par.gls<-par.gls
out.gls$coef<-cor1
#return(list(Sigma=W0,W=W,par.gls=par.gls,out.gls=out,coef=cor1))
return(out.gls)
}
