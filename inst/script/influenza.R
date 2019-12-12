
library(fda.usc)
################################################################################
load("influenza.rda")
nt<-53*52*2
nrep<-156
a1<-209+ 52
a2<-a1+nrep
nmodel<-12
flu2<-flu

# selection response variable (n+1, to n+5)
ii0<-c(12,21,20,19,18)

tres <- flu$temp$data
tres[tres>10] <- 0
tres[flu$temp$data<=10] <- flu$temp$data[flu$temp$data<=10]-10

ff1 <- y~ratet13
ff2 <- y~temp
ff3 <- y~tres
ff4 <- y~irra 

ff5 <- y~ratet13 + temp
ff6 <- y~ratet13 + tres
ff7 <- y~ratet13 + irra
ff8 <- y~temp + irra

nmodel<-24
ahead<-nahead<-4
a1<-209+ 52   
a2<-a1+53

a1<-209+ 26
a2<-a1+53+26

nrep<-500
indi<-1:nrep
iispace<-3500
indspace<-flu$df$ispace>iispace
indspace2<-sum(unique(flu$df$ispace)>iispace)
pred.res<-array(NA,dim=c(nahead,nrep,nmodel+1,indspace2)) #nmodel +1 (la y)
pred0<-pred<-array(NA,dim=c(nahead,nrep,nmodel))

lev<-names(table(flu$df$ispace))
indi3<-flu$df$ispace<16000000

# 10 o 53 comarcas
nrep<-500

######epidemic<-1:156 # in the paper, high time-consuming
## Moderate Time-consuming (only 17 rolling iterations (last 17 weeks)
p1<-1
q1<-0

p2<-1
q2<-0

ii1<-c(12,21,20,19)
## High Time-consuming
levcomar<-names(table(flu$df$ispace[indspace]))
ncomar<-length(levcomar)
epidemic<-c(145:164+52,145:164+2*52) 

for (i in epidemic) {
  for (icomar in 1:ncomar) {  
    j<-1
    flu2<-flu
    flu$df$y<-flu$df[,ii1[j]]
    ind<-which(flu$df$itime>=(i+a1) & flu$df$itime<=(i+a2)
               & flu$df$ispace==levcomar[icomar])
    flu2$df<-flu$df[ind,]
    flu2$temp<-flu$temp[ind]
    flu2$tmp.d1<-flu$tmp.d1[ind]
    flu2$temp.d2<-flu$temp.d2[ind]
    flu2$temt10<-flu$tempt10[ind]
    flu2$ratet13<-flu$ratet13[ind]
    flu2$temp.max<-flu$temp.max[ind]
    flu2$temp.min<-flu$temp.min[ind]
    flu2$vtemp<-flu$rtemp[ind]      
    flu2$hum<-flu$hum[ind,8:14]
    flu2$hum$rangeval<-c(8,14)
    flu2$tempt10<-flu$tempt10[ind,]
    flu2$tres<-flu$tres[ind,]
    flu2$trs.d1<-flu$trs.d1[ind,]        
    flu2$irra<-flu$irra[ind]
    flu2$humt10<-flu$humt10[ind,]
    res.lm1<-fregre.lm(ff1,data=flu2)
    res.lm2<-fregre.lm(ff2,data=flu2)
    res.lm3<-fregre.lm(ff3,data=flu2)
    res.lm4<-fregre.lm(ff4,data=flu2)
    res.lm5<-fregre.lm(ff5,data=flu2)
    res.lm6<-fregre.lm(ff6,data=flu2)
    res.lm7<-fregre.lm(ff7,data=flu2)
    res.lm8<-fregre.lm(ff8,data=flu2)
    
    # par.cor1<-list("cor.ARMA"=list("index"=c("itime"),"group"="ispace",method="lm","p"=1)r)
    # res.fou5<-fregre.igls(ff2,data=flu2,correlation=par.cor1)
    res1<-fregre.gls(ff1,data=flu2,correlation=corARMA(p=p1,q=q1,form=~itime),method="ML")               
    res2<-fregre.gls(ff2,data=flu2,correlation=corARMA(p=p1,q=q1,form=~itime),method="ML")           
    res3<-fregre.gls(ff3,data=flu2,correlation=corARMA(p=p1,q=q1,form=~itime),method="ML")         
    res4<-fregre.gls(ff4,data=flu2,correlation=corARMA(p=p1,q=q1,form=~itime),method="ML")          
    res5<-fregre.gls(ff5,data=flu2,correlation=corARMA(p=p1,q=q1,form=~itime),method="ML")          
    res6<-fregre.gls(ff6,data=flu2,correlation=corARMA(p=p1,q=q1,form=~itime),method="ML")          
    res7<-fregre.gls(ff7,data=flu2,correlation=corARMA(p=p1,q=q1,form=~itime),method="ML")          
    res8<-fregre.gls(ff8,data=flu2,correlation=corARMA(p=p1,q=q1,form=~itime),method="ML")          
    
    res.gls1<-fregre.gls(ff1,data=flu2,correlation=corARMA(p=p2,q=q2,form=~itime),method="ML")       
    res.gls2<-fregre.gls(ff2,data=flu2,correlation=corARMA(p=p2,q=q2,form=~itime),method="ML")       
    res.gls3<-fregre.gls(ff3,data=flu2,correlation=corARMA(p=p2,q=q2,form=~itime),method="ML")       
    res.gls4<-fregre.gls(ff4,data=flu2,correlation=corARMA(p=p2,q=q2,form=~itime),method="ML")          
    res.gls5<-fregre.gls(ff5,data=flu2,correlation=corARMA(p=p2,q=q2,form=~itime),method="ML")       
    res.gls6<-fregre.gls(ff6,data=flu2,correlation=corARMA(p=p2,q=q2,form=~itime),method="ML")      
    res.gls7<-fregre.gls(ff7,data=flu2,correlation=corARMA(p=p2,q=q2,form=~itime),method="ML")      
    res.gls8<-fregre.gls(ff8,data=flu2,correlation=corARMA(p=p2,q=q2,form=~itime),method="ML")      
       ind2<-which(flu$df$itime==(i+a2+1)& flu$df$ispace==levcomar[icomar])
    
    flu3<-flu
    flu3$df$y<-flu$df[,ii0[j]]
    flu3$df<-flu$df[ind2,,drop=FALSE] 
    flu3$temp<-flu$temp[ind2]
    flu3$tmp.d1<-flu$tmp.d1[ind2]
    flu3$tempt10<-flu$tempt10[ind2]  
    flu3$tres<-flu$tres[ind2]  
    flu3$trs.d1<-flu$trs.d1[ind2]  
    flu3$temp.min<-flu$temp.min[ind2]    
    flu3$temp.max<-flu$temp.max[ind2]       
    flu3$ratet13<-flu$ratet13[ind2]  
    flu3$hum<-flu$hum[ind2,8:14]
    flu3$hum$rangeval<-c(8,14)
    flu3$irra<-flu$irra[ind2]
    
    # prediction
    pred.lm1<-predict.fregre.lm(res.lm1,flu3)
    pred.lm2<-predict.fregre.lm(res.lm2,flu3)
    pred.lm3<-predict.fregre.lm(res.lm3,flu3)
    pred.lm4<-predict.fregre.lm(res.lm4,flu3)  
    pred.lm5<-predict.fregre.lm(res.lm5,flu3)
    pred.lm6<-predict.fregre.lm(res.lm6,flu3)  
    pred.lm7<-predict.fregre.lm(res.lm7,flu3)  
    pred.lm8<-predict.fregre.lm(res.lm8,flu3)  
    
    pred1<-predict.fregre.gls(res1,flu3,n.ahead=ahead)
    pred2<-predict.fregre.gls(res2,flu3,n.ahead=ahead)
    pred3<-predict.fregre.gls(res3,flu3,n.ahead=ahead)
    pred4<-predict.fregre.gls(res4,flu3,n.ahead=ahead)  
    pred5<-predict.fregre.gls(res5,flu3,n.ahead=ahead)
    pred6<-predict.fregre.gls(res6,flu3,n.ahead=ahead)
    pred7<-predict.fregre.gls(res7,flu3,n.ahead=ahead)
    pred8<-predict.fregre.gls(res8,flu3,n.ahead=ahead)
    
      yresp<-t(flu3$df[,c("ratet1","ratet2","ratet3","ratet4")])
    
    
    pred.res[,i,1,icomar]<-rbind(pred.lm1,pred.lm1,pred.lm1,pred.lm1)
    pred.res[,i,2,icomar]<-rbind(pred.lm2,pred.lm2,pred.lm2,pred.lm2)
    pred.res[,i,3,icomar]<-rbind(pred.lm3,pred.lm3,pred.lm3,pred.lm3)
    pred.res[,i,4,icomar]<-rbind(pred.lm4,pred.lm4,pred.lm4,pred.lm4)
    pred.res[,i,5,icomar]<-rbind(pred.lm5,pred.lm5,pred.lm5,pred.lm5)
    pred.res[,i,6,icomar]<-rbind(pred.lm6,pred.lm6,pred.lm6,pred.lm6)  
    pred.res[,i,7,icomar]<-rbind(pred.lm7,pred.lm7,pred.lm7,pred.lm7)  
    pred.res[,i,8,icomar]<-rbind(pred.lm8,pred.lm8,pred.lm8,pred.lm8)  
    pred.res[,i,9,icomar]<-t(pred1)
    pred.res[,i,10,icomar]<-t(pred2 )  
    pred.res[,i,11,icomar]<-t(pred3)
    pred.res[,i,12,icomar]<-t(pred4)  
    pred.res[,i,13,icomar]<-t(pred5)
    pred.res[,i,14,icomar]<-t(pred6 )   
    pred.res[,i,15,icomar]<-t(pred7 )   
    pred.res[,i,16,icomar]<-t(pred8 )   
      pred.res[,i,17,icomar]<-yresp
  }
  iresp <- 17
  pred[,i,1]<-apply((pred.res[,i,iresp,]-pred.res[,i,1,])^2,1,mean)
  pred[,i,2]<-apply((pred.res[,i,iresp,]-pred.res[,i,2,])^2,1,mean)
  pred[,i,3]<-apply((pred.res[,i,iresp,]-pred.res[,i,3,])^2,1,mean)
  pred[,i,4]<-apply((pred.res[,i,iresp,]-pred.res[,i,4,])^2,1,mean)
  pred[,i,5]<-apply((pred.res[,i,iresp,]-pred.res[,i,5,])^2,1,mean)
  pred[,i,6]<-apply((pred.res[,i,iresp,]-pred.res[,i,6,])^2,1,mean)
  pred[,i,7]<-apply((pred.res[,i,iresp,]-pred.res[,i,7,])^2,1,mean)
  pred[,i,8]<-apply((pred.res[,i,iresp,]-pred.res[,i,8,])^2,1,mean)
  pred[,i,9]<-apply((pred.res[,i,iresp,]-pred.res[,i,9,])^2,1,mean)
  pred[,i,10]<-apply((pred.res[,i,iresp,]-pred.res[,i,10,])^2,1,mean)
  pred[,i,11]<-apply((pred.res[,i,iresp,]-pred.res[,i,11,])^2,1,mean)
  pred[,i,12]<-apply((pred.res[,i,iresp,]-pred.res[,i,12,])^2,1,mean)
  pred[,i,13]<-apply((pred.res[,i,iresp,]-pred.res[,i,13,])^2,1,mean)
  pred[,i,14]<-apply((pred.res[,i,iresp,]-pred.res[,i,14,])^2,1,mean)
  pred[,i,15]<-apply((pred.res[,i,iresp,]-pred.res[,i,15,])^2,1,mean)
  pred[,i,16]<-apply((pred.res[,i,iresp,]-pred.res[,i,16,])^2,1,mean)
  print( i)
  pred[,i,]
}

save.image(file="fluExample.RData")
###############################################################################
load(file="fluExample.RData")

table6<-t(round(apply(pred[,epidemic,1:16],c(1,3),mean,na.rm=T),2))
table6<-cbind(table6[c(1:8),1],table6[c(9:16),1],table6[c(1:8),2],table6[c(9:16),2])
colnames(table6)<-c("FLM (t+1)","FGLS (t+1)","FLM (t+2)","FGLS (t+2)")
rownames(table6)<-c("Rate(w)","Temp(t)","Temp.thres(t)","SR(t)",
                    "Rate(w),Temp(t)","Rate(w),Temp(t)",
                    "Rate(w)Temp.thres(t)","Temp(t),SR(t)")
table6
###############################################################################
j<-1
flu2<-flu
flu$df$y<-flu$df[,ii1[j]]

ind<-which(flu$df$itime>=(epidemic[1]+a1) & flu$df$itime<=(epidemic[40]+a2)
           & flu$df$ispace  %in% levcomar)
length(ind)
flu2$df<-flu$df[ind,]
flu2$temp<-flu$temp[ind]
flu2$tmp.d1<-flu$tmp.d1[ind]
flu2$temp.d2<-flu$temp.d2[ind]
flu2$temt10<-flu$tempt10[ind]
flu2$ratet13<-flu$ratet13[ind]
flu2$temp.max<-flu$temp.max[ind]
flu2$temp.min<-flu$temp.min[ind]
flu2$vtemp<-flu$rtemp[ind]      
flu2$hum<-flu$hum[ind,8:14]
flu2$hum$rangeval<-c(8,14)
flu2$tempt10<-flu$tempt10[ind,]
flu2$tres<-flu$tres[ind,]
flu2$trs.d1<-flu$trs.d1[ind,]        
flu2$irra<-flu$irra[ind]
flu2$humt10<-flu$humt10[ind,]

res1<-fregre.gls(ff1,data=flu2,correlation=corARMA(p=p1,q=q1,form=~itime|ispace),method="ML")               
res2<-fregre.gls(ff2,data=flu2,correlation=corARMA(p=p1,q=q1,form=~itime|ispace),method="ML")           
res3<-fregre.gls(ff3,data=flu2,correlation=corARMA(p=p1,q=q1,form=~itime|ispace),method="ML")         
res4<-fregre.gls(ff4,data=flu2,correlation=corARMA(p=p1,q=q1,form=~itime|ispace),method="ML")          
res5<-fregre.gls(ff5,data=flu2,correlation=corARMA(p=p1,q=q1,form=~itime|ispace),method="ML")          
res6<-fregre.gls(ff6,data=flu2,correlation=corARMA(p=p1,q=q1,form=~itime|ispace),method="ML")          
res7<-fregre.gls(ff7,data=flu2,correlation=corARMA(p=p1,q=q1,form=~itime|ispace),method="ML")          
res8<-fregre.gls(ff8,data=flu2,correlation=corARMA(p=p1,q=q1,form=~itime|ispace),method="ML")          


################# INICIO #####################
add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}
pdf(file = "S2_fig.pdf",width = 10.24, height =7.68, pointsize = 10)
par(mfrow=c(2,2))
cols=c("blue","sky blue","red","dark red")

# Rate
covar<- "ratet13"
res <-res1
x<-flu2[[covar]]
a<-fdata(res$beta.l[[covar]],flu2[[covar]]$argval,flu2[[covar]]$rangeval)  
b1<-b<- +inprod.fdata(x,a) 
qq<-quantile(b,prob=c(.25,.5,.75,1))
ind1<-ifelse(b[,1]<qq[1],TRUE,FALSE)
ind2<-ifelse(b[,1]>=qq[1],TRUE,FALSE)
ind3<-ifelse(b[,1]>=qq[2],TRUE,FALSE)
ind4<-ifelse(b[,1]>=qq[3],TRUE,FALSE)
mn0<-expression(paste( beta[rate]))
ylab <- "Rate(w)"
ylab <- expression(Rate["s,n"](w))
mn2<-expression(paste(symbol("<"), list(Rate["s,n"](w), beta[Rate]),symbol(">")))
plot(c((func.mean(x[ind1,])),(func.mean(x[ind2,]))
       ,(func.mean(x[ind3,])),(func.mean(x[ind4,])))
     ,lwd=2,lty=1
     ,col=cols,
     main=mn2,xlab="Range time w: 13 weeks",ylab=ylab)

# Temp
covar<- "temp"
res <-res2
x<-flu2[[covar]]
a<-fdata(res$beta.l[[covar]],flu2[[covar]]$argval,flu2[[covar]]$rangeval)  
b1<-b<- +inprod.fdata(x,a) 
qq<-quantile(b,prob=c(.25,.5,.75,1))
ind1<-ifelse(b[,1]<qq[1],TRUE,FALSE)
ind2<-ifelse(b[,1]>=qq[1],TRUE,FALSE)
ind3<-ifelse(b[,1]>=qq[2],TRUE,FALSE)
ind4<-ifelse(b[,1]>=qq[3],TRUE,FALSE)
mn0<-expression(paste( beta[rate]))
ylab <- expression(Temp["s,n"](t))
mn2<-expression(paste(symbol("<"), list(Temp["s,n"](t), beta[Temp]),symbol(">")))
plot(c((func.mean(x[ind1,])),(func.mean(x[ind2,]))
       ,(func.mean(x[ind3,])),(func.mean(x[ind4,])))
     ,lwd=2,lty=1
     ,col=cols,
     main=mn2,xlab="Range time t: 14 days",ylab=ylab)


# Temp
covar<- "tres"
res <-res3
x<-flu2[[covar]]
a<-fdata(res$beta.l[[covar]],flu2[[covar]]$argval,flu2[[covar]]$rangeval)  
b1<-b<- +inprod.fdata(x,a) 
qq<-quantile(b,prob=c(.25,.5,.75,1))
ind1<-ifelse(b[,1]<qq[1],TRUE,FALSE)
ind2<-ifelse(b[,1]>=qq[1],TRUE,FALSE)
ind3<-ifelse(b[,1]>=qq[2],TRUE,FALSE)
ind4<-ifelse(b[,1]>=qq[3],TRUE,FALSE)
ylab <- expression(Temp.thres["s,n"](t))
mn2<-expression(paste(symbol("<"), list(Temp.thres["s,n"](t), beta[Temp.thres]),symbol(">")))
plot(c((func.mean(x[ind1,])),(func.mean(x[ind2,]))
       ,(func.mean(x[ind3,])),(func.mean(x[ind4,])))
     ,lwd=2,lty=1
     ,col=cols,
     #  main=mn2,xlab="t: previous 14 days",ylab=ylab)
     main=mn2,xlab="t in days",ylab=ylab)

# Solar Radiation
covar<- "irra"
res <-res4
x<-flu2[[covar]]
a<-fdata(res$beta.l[[covar]],flu2[[covar]]$argval,flu2[[covar]]$rangeval)  
b1<-b<- +inprod.fdata(x,a) 
qq<-quantile(b,prob=c(.25,.5,.75,1))
ind1<-ifelse(b[,1]<qq[1],TRUE,FALSE)
ind2<-ifelse(b[,1]>=qq[1],TRUE,FALSE)
ind3<-ifelse(b[,1]>=qq[2],TRUE,FALSE)
ind4<-ifelse(b[,1]>=qq[3],TRUE,FALSE)
ylab <- expression(SR["s,n"](t))
mn2<-expression(paste(symbol("<"), list(SR["s,n"](t), beta[SR]),symbol(">")))
plot(c((func.mean(x[ind1,])),(func.mean(x[ind2,]))
       ,(func.mean(x[ind3,])),(func.mean(x[ind4,])))
     ,lwd=2,lty=1
     ,col=cols,
     main=mn2,xlab="t (in days)",ylab=ylab)
leg <- c(expression(v[q[1]]), expression(v[q[2]]), expression(v[q[3]]) ,expression(v[q[4]]))
#legend(8,2500,inset=c(-0.2,0),col=cols,legend=leg,cex=.8,bty="n",lty=1,lwd=1)
add_legend("bottom", legend=leg,lty=1, #pch=20, 
           col=cols, lwd=2,
           horiz=TRUE, bty='n',
           cex=0.95)
dev.off()
################# fin Gráfico #####################


################################################################################
################################################################################
# Figure 1      Series Temporales 
flu4<-aggregate(flu2$df,list(flu2$df[,"itime"]),mean)
#cambiar: utilizar un ggplot o algo así más chulo (ver links/linkedIN)
png(file = "S1_fig.png")   
par(mfrow=c(4,1))
FluRate<-ts(exp(flu4[,"ratet0"]),2001,2011,53)  
Temperature<-ts(flu4[,"tmed"] ,2001,2011,53)
SolarRadiation<-ts(flu4[,"imed"] ,2001,2011,53)
RelativeHumidity<-ts(flu4[,"hmed"] ,2001,2011,53)
TemperatureThreshold<-ts(flu4[,"thold"] ,2001,2011,53)
aa<-cbind(FluRate,Temperature,SolarRadiation,RelativeHumidity)
plot(aa,main="",xlab="Years")
dev.off()
################################################################################
################################################################################


################################################################################
################################################################################
# Pearson and distance correlation
################################################################################
################################################################################
#Upper diagonal matrix of Table 1: Pearson correlation
# no se muestran
A<-cor(flu2$df[,c(19,12,22,13,23,26,29)])   

################################################################################
#Lower diagonal matrix of Table 1: Distance correlation
################################################################################
# Distance computation for  Multivariate-Multivariate data
names(flu)
A1<-metric.lp(flu2$ratet13)
A2<-metric.lp(flu2$temp)
A3<-metric.lp(flu2$tres)
A4<-metric.lp(flu2$irra)
A5<-metric.lp(flu2$hum)
A6<-metric.dist(as.matrix(flu2$df$ratet0))
A7<-metric.dist(as.matrix(flu2$df$ratet1))
A8<-metric.dist(as.matrix(flu2$df$ratet2))
A9<-metric.dist(as.matrix(flu2$df$tmed))
A10<-metric.dist(as.matrix(flu2$df$thold))
A11<-metric.dist(as.matrix(flu2$df$imed))
A12<-metric.dist(as.matrix(flu2$df$hmed))
# Distance correlation
A<-cor(flu2$df[,c(21,12,22,13,23,26,29)])   

names(flu2$df[,c(21,12,22,13,23,26,29)])   
A[2,1]<-dcor.dist(A8,A7)
A[3,1]<-dcor.dist(A8,A6)
A[4,1]<-dcor.dist(A8,A9)
A[5,1]<-dcor.dist(A8,A10)
A[6,1]<-dcor.dist(A8,A11)
A[7,1]<-dcor.dist(A8,A12)
A[3,2]<-dcor.dist(A7,A6)
A[4,2]<-dcor.dist(A7,A9)
A[5,2]<-dcor.dist(A7,A10)
A[6,2]<-dcor.dist(A7,A11)
A[7,2]<-dcor.dist(A7,A12)
A[4,3]<-dcor.dist(A6,A9)
A[5,3]<-dcor.dist(A6,A10)
A[6,3]<-dcor.dist(A6,A11)
A[7,3]<-dcor.dist(A6,A12)
A[5,4]<-dcor.dist(A9,A10)
A[6,4]<-dcor.dist(A9,A11)
A[7,4]<-dcor.dist(A9,A12)
A[6,5]<-dcor.dist(A10,A11)
A[7,5]<-dcor.dist(A10,A12)
A[7,6]<-dcor.dist(A11,A12)
round(A,2)
################################################################################

################################################################################
B<-matrix(NA,5,12)
# Distance computation for Functional-Multivariate data
B[1,1]<-dcor.dist(A8,A1)
B[2,1]<-dcor.dist(A8,A2)
B[3,1]<-dcor.dist(A8,A3)
B[4,1]<-dcor.dist(A8,A4)
B[5,1]<-dcor.dist(A8,A5)

B[1,2]<-dcor.dist(A7,A1)
B[2,2]<-dcor.dist(A7,A2)
B[3,2]<-dcor.dist(A7,A3)
B[4,2]<-dcor.dist(A7,A4)
B[5,2]<-dcor.dist(A7,A5)

B[1,3]<-dcor.dist(A6,A1)
B[2,3]<-dcor.dist(A6,A2)
B[3,3]<-dcor.dist(A6,A3)
B[4,3]<-dcor.dist(A6,A4)
B[5,3]<-dcor.dist(A6,A5)

B[1,4]<-dcor.dist(A9,A1)
B[2,4]<-dcor.dist(A9,A2)
B[3,4]<-dcor.dist(A9,A3)
B[4,4]<-dcor.dist(A9,A4)
B[5,4]<-dcor.dist(A9,A5)

B[1,5]<-dcor.dist(A10,A1)
B[2,5]<-dcor.dist(A10,A2)
B[3,5]<-dcor.dist(A10,A3)
B[4,5]<-dcor.dist(A10,A4)
B[5,5]<-dcor.dist(A10,A5)

B[1,6]<-dcor.dist(A11,A1)
B[2,6]<-dcor.dist(A11,A2)
B[3,6]<-dcor.dist(A11,A3)
B[4,6]<-dcor.dist(A11,A4)
B[5,6]<-dcor.dist(A11,A5)

B[1,7]<-dcor.dist(A12,A1)
B[2,7]<-dcor.dist(A12,A2)
B[3,7]<-dcor.dist(A12,A3)
B[4,7]<-dcor.dist(A12,A4)
B[5,7]<-dcor.dist(A12,A5)

# Distance computation for Functional-Functional data
B[1,8]<- B[2,9]<- B[3,10]<- B[4,11]<- B[5,12]<-1
B[2,8]<-dcor.dist(A1,A2)
B[3,8]<-dcor.dist(A1,A3)
B[4,8]<-dcor.dist(A1,A4)
B[5,8]<-dcor.dist(A1,A5)

B[3,9]<-dcor.dist(A2,A3)
B[4,9]<-dcor.dist(A2,A4)
B[5,9]<-dcor.dist(A2,A5)

B[4,10]<-dcor.dist(A3,A4)
B[5,10]<-dcor.dist(A3,A5)
B[5,11]<-dcor.dist(A4,A5)

B[1,8:12]<-B[,8]
B[2,9:12]<-B[2:5,9]
B[3,10:12]<-B[3:5,10]
B[4,11:12]<-B[4:5,11]

rnames<-c("raten+2","raten+1","rate","Temp","Thold","Irra","Hum")
fnames<-c("rate(w)","Temp(t)","Thold(t)","Irra(t)","Hum(t)")

colnames(B)<-c(rnames,fnames)
rownames(B)<-fnames
# TABLE 5
round(B[,c(8:12,2,1)],2)

save.image(file="fluExample.RData")






