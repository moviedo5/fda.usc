rm(list=ls())
###################################################
###################################################
#install.packages("fda.usc")
library(fda.usc)
library(lattice)

traza<-fda.usc:::traza
source(".\\flm_ar.R")
source(".\\predict.fregre.lm0.R")

ar2.1<-c( 0.50,0.45)
ar2.2 <-c(1.40, -0.45)
ar2.3 <-c(1.5,-0.75)
ar2.4 <-c( 1.735736 ,-0.902500)
#########################
method="ML"
ctr=list(niterEM=0,
         #optimMethod="L-BFGS-B",
         optimMethod="L-BFGS-B",
         gradHess = TRUE,
         opt = "optim")

rtt<-c(0,1)
n.ahead=n2=10

J=101; n1<-250
tt<-seq(0,1,len=J)
n<-n1+n2

r2=c(.05)
#ara el ar1
nrep<-1000
kmax<-4
ff<-as.formula("y~x")
model2<-model<-c("lm.pc","gls.pc1","gls.pc2","igls.pc1","igls.pc2","igls.pcp")
nmodel<-length(model)
nmodel2<-length(model2)


rho = list(ar2.1=ar2.1,ar2.2=ar2.2,ar2.3=ar2.3,ar2.4=ar2.4)
fit11<-fit1<-fit0<-array(NA,dim=c(length(r2),length(rho),nrep,nmodel))
pred<-array(NA,dim=c(length(r2),length(rho),nrep,nmodel2,n2))
ypred<-array(NA,dim=c(length(r2),length(rho),nrep,nmodel2,n2))
yfit<-array(NA,dim=c(length(r2),length(rho),nrep,nmodel2,n1))
yobs<-array(NA,dim=c(length(r2),length(rho),nrep,n2))

cpu1<-cpu2<-array(0,dim=c(length(r2),length(rho),nmodel2))
p<-numeric(13)

iterations3<- iterations4<- iterations7<- iterations8<-a<-b<-bx<-bb<-array(0,dim=c(length(r2),length(rho),nrep))
bet1<-bet2<-array(NA,dim=c(length(r2),length(rho),nrep,J))  
criteria="GCCV1"
ind<-as.double(1:n)
nls<-1
repetir<-NULL
set.seed(1:6)
p1<-1
p2<-2
q1<-0
nrep<-1000
for (i in 1:length(r2)) {
  for (j in 1:length(rho)) {
    print(rho[[j]])
    k<-1
    while (k<nrep) {
      #  cat(i);cat(j);print(k)
      nls<-1
      fx0=rproc2fdata(n,tt,sigma="wiener")#,par.list=list("scale"=.5))
      fdatos.a <- fdata.cfs.2003.a(fx0,snr=r2[i],rho=rho[[j]],p=1) #p=1 modelo lineal
      fx<-fdatos.a$x     
      x<-fx[1:n1]
      newx<-fx[(n1+1):n]                                        
      bet<-fdata2fd(fdatos.a$bet)
      y <-fdatos.a$y[1:n1]
      dff<-data.frame(y,ind=ind[1:n1])
      ldata<-list("df"=dff,"x"=x)
      pc<-list("x"=create.pc.basis(x,1:4,norm=FALSE))
      res1<-fregre.lm(ff, data=ldata,basis.x=pc) 
      tr2<-try( res2<-suppressMessages(suppressWarnings(fregre.gls(ff, data=ldata,correlation=corAR1(form=~ind),basis.x=pc))))
      #tr3<-try(res3<- suppressMessages(suppressWarnings(fregre.gls(ff, data=ldata,correlation=corARMA(form=~ind,p=p2,q=q1),basis.x=pc))))
      #if (class(tr3)[1]!= "try-error")  res3=res2
      tr4<-try(res4<-suppressMessages(suppressWarnings(fregre.igls(ff,data=ldata,basis.x=pc,correlation=list("cor.ARMA"=list()),control=list("p"=p1)))))
      tr5<-try(res5<-suppressMessages(suppressWarnings(fregre.igls(ff,data=ldata,basis.x=pc,correlation=list("cor.ARMA"=list()),control=list("p"=p2)) )))
      tr6<-try(res6<-suppressMessages(suppressWarnings(fregre.igls(ff,data=ldata,basis.x=pc,correlation=list("cor.AR"=list()),control=list("order.max"=4)))))
      if (all(c(class(tr2)[1],class(tr4)[1],class(tr6)[1],class(tr5)[1],
                class(tr6)[1])!= "try-error") ) {
        #print("no ha petado")
        
        # Peta al intentar estimar un AR2 <1
        pred01<-pred02<-pred03 <-pred04 <-pred05 <-pred06 <-pred07 <-pred08 <-rep(NA,n2)
        newy<-fdatos.a$y[(n1+1):n]
        newy2<-fdatos.a$y[(n1+1):n]-fdatos.a$e[(n1+1):n]
        newdff<-data.frame(y=newy,ind=ind[(n1+1):n])
        newldata<-list("df"=newdff,"x"=newx)
        res1$rn<-FALSE
        pred01<-predict.fregre.lm(res1,newldata)
        pred02<-predict.fregre.gls(res2,newldata)
        pred04<-predict.fregre.igls(res4,newldata)
        pred05<-predict.fregre.igls(res5,newldata)
        pred06<-predict.fregre.igls(res6,newldata)
        ypred[i,j,k,1,]<-pred01
        ypred[i,j,k,2,]<-pred02
        #ypred[i,j,k,3,]<-pred03 
        ypred[i,j,k,4,]<-pred04
        ypred[i,j,k,5,]<-pred05
        ypred[i,j,k,6,]<-pred06
        pred[i,j,k,1,]<-(pred01-newy)^2
        pred[i,j,k,2,]<-(pred02-newy)^2
        #pred[i,j,k,3,]<-(pred03-newy)^2
        pred[i,j,k,4,]<-(pred04-newy)^2
        pred[i,j,k,5,]<-(pred05-newy)^2
        pred[i,j,k,6,]<-(pred06-newy)^2
        k<-k+1
      } #fin while
    }
    cat(i);cat(j);print(k)
  }
  cat(i);cat(j);print(k)
}
#############################################################################

round(apply(pred[,,1:(k-1),,c(1,5,10),drop=FALSE],c(4,5,1:2),mean,na.rm=T),4)#The AR(1) and AR(2) models fail to capture the cycle dynamics in the data, as they do not have any complex conjugate roots
#http://www.sciencedirect.com/science/article/pii/S0047259X08001590


#The misspecification in the correlation structure for GLS and the smoothing approach is created by using as a working model an AR(1) instead of an AR(2). The amount of misspecification depends on the correlation strength of the generated structure with AR(2); see Figure 2. It is clear that estimating the correlation structure using an AR(1) process will capture mostly frequencies around zero, whereas it will represent poorly frequencies further away from zero.
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2738954/
  
#pred.p1.ar2<-pred
#pred.p1.b<-pred
#pred.p0.b<-pred
#save.image("simu_A_2017_AR2.RData")
#load("simu_A_2017_AR2.RData")
#load("simu_A_2017_AR2_ar150_075.RData")
#save.image("simu_A_2017_AR2_ar150_075.RData")
bb<-round(apply(pred[,,,,c(1,5,10),drop=FALSE],c(4,5,1:2),median,na.rm=T),3)
library(xtable)
xtable(cbind(bb[,,1,1],bb[,,1,2],bb[,,1,3]),digits=3)

model2<-c( "lm.pc"  ,   "gls.pc" ,   "igls.pc" ,
           "igls2.pc" ,"BSP",
           "fgam","fgsamm","figsam")

imod <- c(1,2,3,6,7,8)
horiz <- c(1,5,10)
round(apply(pred[,,,imod,horiz,drop=FALSE],c(4,5,1:2),mean,na.rm=T),3)
dimnames(pred) <- list(paste("snr",r2),
                       paste("rho",rho),
                       1:nrep,
                       model2,
                       paste("t+",1:10)
                       )
round(tab1<-apply(pred[,,,imod,horiz,drop=FALSE],c(4,5,1:2),mean,na.rm=T),3)
cbind(tab1[,,1,1])
cbind(tab1[,,2,1])
#############################################################################
