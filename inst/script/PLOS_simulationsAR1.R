semimetric.pc<-function (fdata1, fdata2 = fdata1, q = 1, ...) 
{
    C1 <- match.call()
    if (is.fdata(fdata1)) {
        tt <- fdata1[["argvals"]]
        rtt <- fdata1[["rangeval"]]
        nas1 <- apply(fdata1$data, 1, count.na)
        if (any(nas1)) 
            stop("fdata1 contain ", sum(nas1), " curves with some NA value \n")
        else if (!is.fdata(fdata2)) {
            fdata2 <- fdata(fdata2, tt, rtt)
        }
        nas2 <- apply(fdata2$data, 1, count.na)
        if (any(nas2)) 
            stop("fdata2 contain ", sum(nas2), " curves with some NA value \n")
        DATA1 <- fdata1[["data"]]
        DATA2 <- fdata2[["data"]]
        range.t <- rtt
    }
    else {
        if (is.vector(fdata1)) 
            fdata1 <- as.matrix(t(fdata1))
        if (is.vector(fdata2)) 
            fdata2 <- as.matrix(t(fdata2))
        DATA1 <- fdata1
        DATA2 <- fdata2
        range.t <- c(1, ncol(DATA1))
    }
    testfordim <- sum(dim(DATA1) == dim(DATA2)) == 2
    twodatasets <- TRUE
    if (testfordim) 
        twodatasets <- sum(DATA1 == DATA2) != prod(dim(DATA1))
    if (class(q)=="fdata.comp") {
        pc<-q
        basis<-pc$rotation
        Xcen.fdata <- fdata.cen(fdata2,pc$mean)[[1]]                                        
        q<-nrow(basis)
        COMPONENT2<-inprod.fdata(Xcen.fdata, basis,...)
        }
    else {
        pc<-fdata2pc(fdata2,q)    
        basis<-pc$rotation
        COMPONENT2<-pc$x[,1:q]
        }  
    qmax <- ncol(DATA1)
    if (q > qmax) 
        stop(paste("give a integer q smaller than ", qmax))
    n <- nrow(DATA1)
#    COVARIANCE <- t(DATA1) %*% DATA1/n
#    ei = eigen(COVARIANCE, symmetric = TRUE)
#    EIGENVECTORS <- matrix(ei$vectors[, 1:q], ncol = q)
#    COMPONENT1 <- DATA1 %*% EIGENVECTORS  
    if (twodatasets) {
             Xcen.fdata <- fdata.cen(fdata1, pc$mean)[[1]]
#            COMPONENT2<-inprod.fdata(Xcen.fdata, bases, ...)
            COMPONENT1<-inprod.fdata(Xcen.fdata, basis,...)
    }
    else {
        COMPONENT1 <- COMPONENT2
    }     
    SEMIMETRIC <- 0
    for (qq in 1:q) SEMIMETRIC <- SEMIMETRIC + outer(COMPONENT1[, 
        qq], COMPONENT2[, qq], "-")^2
    mdist <- sqrt(SEMIMETRIC)
    attr(mdist, "call") <- "semimetric.pc"
    attr(mdist, "par.metric") <- list(q = pc)
#    attr(mdist, "par.pc") <- pc
    return(mdist)
}

ind.list=function (lista, ind)
{
  len <- length(lista)
  for (i in 1:len) lista[[i]] <- lista[[i]][ind, , drop = FALSE]
  lista
}

###################################################
###################################################
#install.packages("fda.usc")
library(fda.usc)
library(lattice)

###################################################
traza<-fda.usc:::traza

source(".\\flm_ar.R")

rtt<-c(0,1)
n.ahead=n2=10

J=101; n1<-250
tt<-seq(0,1,len=J)
n<-n1+n2

r2=c(0.05,.1)
#ara el ar1
nrep<-1000
kmax<-4
ff<-as.formula("y~x")
ffs<-y~s(x)
#model2<-c("flm.pc","fgls.pc","iflm.pc","i2flm.pc","np.","np.gls","np.igls","np.igls2")
model2<-model<-c("lm.pc_0","gls.pc_1","igls.pc_1","igls.pc_p",
                 "lm.bsp_0","gls.bsp_1","igls.bsp_1","igls.bsp_p")
nmodel<-length(model)
nmodel2<-length(model2)

rho = list(ar2.1=0,ar2.3=0.5,ar2.4=0.9)
fit11<-fit1<-fit0<-array(NA,dim=c(length(r2),length(rho),nrep,nmodel))

# repetir para PC y BSP

pred<-array(NA,dim=c(length(r2),length(rho),nrep,nmodel2,n2))
ypred<-array(NA,dim=c(length(r2),length(rho),nrep,nmodel2,n2))
yfit<-array(NA,dim=c(length(r2),length(rho),nrep,nmodel2,n1))
yobs<-array(NA,dim=c(length(r2),length(rho),nrep,n2))

cpu1<-cpu2<-array(0,dim=c(length(r2),length(rho),nmodel2))
p<-numeric(13)

iterations3<- iterations4<- iterations7<- iterations8<-a<-b<-bx<-bb<-array(0,dim=c(length(r2),length(rho),nrep))
#a<-array(0,dim=c(length(r2),length(rho),nrep,kmax))
bet1<-bet2<-array(NA,dim=c(length(r2),length(rho),nrep,J))  
#criteria="GCCV2"
criteria="GCCV1"
ind<-as.double(1:n)
nls<-1
repetir<-NULL
p1<-1
p2<-2
q1<-0
set.seed(1:6)

for (k in 1:nrep){  
  for (i in 1:length(r2)) {
  for (j in 1:length(rho)) {
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
     
     out2<-fregre.pc.cv(x,y,1:8 ,criteria = "SICc")
     kmax<-length(out2$pc.opt)
     pc<-list("x"=create.pc.basis(x,out2$pc.opt,norm=TRUE))
     
     res<-fregre.basis.cv(ldata$x,ldata$df$y, basis.x=11, basis.b=c(5,7,9)) 
       
     p[7]<-proc.time()[1]   
     # a[i,j,k]<-length(res$pc.opt)    
     a[i,j,k]<-kmax
     bb[i,j,k]<- res$basis.b.opt$nbasis
     bx[i,j,k]<- res$basis.x.opt$nbasis
     bspx<-list("x"=res$basis.x.opt)  
     bsp.x<-list("x"=create.bspline.basis(c(0,1),nbasis=11))#bx[i,j,k]))
     bsp.b<-list("x"=create.bspline.basis(c(0,1),nbasis=bb[i,j,k]))
     res1<-fregre.lm(ff, data=ldata,basis.x=pc) 
     res2<-fregre.gls(ff, data=ldata,correlation=corAR1(form=~ind),basis.x=pc)
     res3<-fregre.igls(ff,data=ldata,basis.x=pc,correlation=list("cor.ARMA"=list()),control=list("p"=p1))
     res4<-fregre.igls(ff,data=ldata,basis.x=pc,correlation=list("cor.AR"=list()),control=list("order.max"=4))
     res5<-fregre.lm(ff, data=ldata,basis.x=bsp.x,basis.b=bsp.b) 
     res6<-fregre.gls(ff, data=ldata,correlation=corAR1(),basis.x=bsp.x,basis.b=bsp.b)
     res7<-fregre.igls(ff,data=ldata,basis.x=bsp.x,basis.b=bsp.b,correlation=list("cor.ARMA"=list()),control=list("p"=p1))
     res8<-fregre.igls(ff,data=ldata,basis.x=bsp.x,basis.b=bsp.b,correlation=list("cor.AR"=list()),control=list("order.max"=4))
     
     fit0[i,j,k,1]<-drop(norm.fdata(res1$beta.l[[1]]-fdatos.a$bet))
     fit0[i,j,k,2]<-drop(norm.fdata(res2$beta.l[[1]]-fdatos.a$bet))
     fit0[i,j,k,3]<-drop(norm.fdata(res3$beta.l[[1]]-fdatos.a$bet))
     fit0[i,j,k,4]<-drop(norm.fdata(res4$beta.l[[1]]-fdatos.a$bet))
     fit0[i,j,k,5]<-drop(norm.fdata(fdata(res5$beta.l[[1]],tt,rtt)-fdatos.a$bet))
     fit0[i,j,k,6]<-drop(norm.fdata(fdata(res6$beta.l[[1]],tt,rtt)-fdatos.a$bet))
     fit0[i,j,k,7]<-drop(norm.fdata(fdata(res7$beta.l[[1]],tt,rtt)-fdatos.a$bet))
     fit0[i,j,k,8]<-drop(norm.fdata(fdata(res8$beta.l[[1]],tt,rtt)-fdatos.a$bet))
    
     fit1[i,j,k,2]<-(rho[j][[1]]-coef(res2$modelStruct,FALSE))^2
     fit1[i,j,k,3]<-(rho[j][[1]]-res3$corStruc$ar$coef)^2
     #  para un ar(p) no tiene sentido
     #fit1[i,j,k,4]<-(rho[j][[1]]-res4$corStruc$ar$coef)^2
     fit1[i,j,k,6]<-(rho[j][[1]]-coef(res6$modelStruct,FALSE))^2
     fit1[i,j,k,7]<-(rho[j][[1]]-res7$corStruc$ar$coef)^2
    # fit1[i,j,k,8]<-(rho[j][[1]]-res8$corStruc$ar$coef)^2
     
     pred01<-pred02<-pred03 <-pred04 <-pred05 <-pred06 <-pred07 <-pred08 <-rep(NA,n2)
     newy<-fdatos.a$y[(n1+1):n]
     newy2<-fdatos.a$y[(n1+1):n]-fdatos.a$e[(n1+1):n]
     newdff<-data.frame(y=newy,ind=ind[(n1+1):n])
     newldata<-list("df"=newdff,"x"=newx)
     res1$rn<-FALSE
     pred01<-predict.fregre.lm(res1,newldata)
     pred02<-predict.fregre.gls(res2,newldata)
     pred03<-predict.fregre.igls(res3,newldata)
     pred04<-predict.fregre.igls(res4,newldata)
     pred05<-predict.fregre.lm(res5,newldata)
     pred06<-predict.fregre.gls(res6,newldata)
     pred07<-predict.fregre.igls(res7,newldata)
     pred08<-predict.fregre.igls(res8,newldata)
     ypred[i,j,k,1,]<-pred01
     ypred[i,j,k,2,]<-pred02
     ypred[i,j,k,3,]<-pred03 
     ypred[i,j,k,4,]<-pred04
     ypred[i,j,k,5,]<-pred05
     ypred[i,j,k,6,]<-pred06
     ypred[i,j,k,7,]<-pred07 
     ypred[i,j,k,8,]<-pred08
     pred[i,j,k,1,]<-(pred01-newy)^2
     pred[i,j,k,2,]<-(pred02-newy)^2
     pred[i,j,k,3,]<-(pred03-newy)^2
     pred[i,j,k,4,]<-(pred04-newy)^2
     pred[i,j,k,5,]<-(pred05-newy)^2
     pred[i,j,k,6,]<-(pred06-newy)^2
     pred[i,j,k,7,]<-(pred07-newy)^2
     pred[i,j,k,8,]<-(pred08-newy)^2
     } #fin while
    #cat(i);cat(j);print(k)
  }
 cat(i);cat(j);print(k)
 }

#############################################################################
# Tabla  number of basis elements
round(apply(a,c(1:2),mean,na.rm=T),1)
#round(apply(bx,c(1:2),mean,na.rm=T),3)
round(apply(bb,c(1:2),mean,na.rm=T),1)

model2<-c( "lm.pc"  ,   "gls.pc" ,   "igls.pc" ,
           "igls2.pc" ,"BSP",
           "fgam","fgsamm","figsam")

imod <- 1:8
horiz <- c(1,5,10)
dimnames(pred) <- list(paste("snr",r2),
                       paste("rho",rho),
                       1:nrep,
                       model2,
                       paste("t+",1:10)
)
dd<-round(apply(pred[,,,,c(1,5,10),drop=FALSE],c(4,5,1:2),mean,na.rm=T),3)
dd
library(xtable)
xtable(cbind(dd[,,1,1],dd[,,1,2],dd[,,1,3]),digits=3)
xtable(cbind(dd[,,2,1],dd[,,2,2],dd[,,2,3]),digits=3)
xtable(cbind(ss[,c(1,5,10),,1],ss[,c(1,5,10),,2],ss[,c(1,5,10),,3]),digits=3)

dimnames(pred)<-list(r2,rho,1:nrep,model2,1:n.ahead)
dimnames(fit1)<-list(r2,rho,1:nrep,model)
dimnames(fit0)<-list(r2,rho,1:nrep,model)
dimnames(cpu2)<-dimnames(cpu1)<-list(r2,rho,model)
dimnames(iterations8)<-dimnames(iterations7)<-dimnames(iterations4)<-dimnames(iterations3)<-dimnames(a)<-dimnames(bb)<-list(r2,rho,1:nrep)

round(apply(pred[,,,,c(1,5,10)],c(4,5,1:2),mean,na.rm=T),3)

# Table 2 MSE of beta parameter, B=1000
ff<-round(apply(fit0,c(4,1:2),mean,na.rm=T),3)
xtable(
  rbind(cbind(ff[1:4,1,1],ff[1:4,1,2],ff[1:4,1,3],ff[1:4+4,1,1],ff[1:4+4,1,2],ff[1:4+4,1,3])
      ,cbind(ff[1:4,2,1],ff[1:4,2,2],ff[1:4,2,3],ff[1:4+4,2,1],ff[1:4+4,2,2],ff[1:4+4,2,3])),digits=3)
round(apply(fit1[,,,c(2,3,6,7)],c(4,1:2),mean,na.rm=T),3)

# Save model (a) for: fdatos.a <- fdata.cfs.2003.a(fx0,snr=r2[i],rho=rho[[j]],p=1) 
# Save model (b) for: fdatos.a <- fdata.cfs.2003.b(fx0,snr=r2[i],rho=rho[[j]],p=1) 

#save.image("simu_A_2017_AR1_20171121.RData")
#load("simu_A_2017_AR2.RData")
#load("simu_A_2017_AR2_ar150_075.RData")
#save.image("simu_A_2017_AR2_ar150_075.RData")
dd<-round(apply(pred[,,,,c(1,5,10),drop=FALSE],c(4,5,1:2),median,na.rm=T),3)
library(xtable)
xtable(cbind(dd[,,1,1],dd[,,1,2],dd[,,1,3]),digits=3)
#############################################################################
