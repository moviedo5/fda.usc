# rm(list=ls())
completa.fdata<-function(y){
  n<-nrow(y)
  p<-ncol(y)
  x<-y$argvals
  ynew<-y
  aa<-na.omit(y)
  ii<-which(is.na.fdata(y))
  for (i in ii){
      app=approx(x,y$data[i,],n=p)
      x=app$x
      ynew$data[i,]<-app$y
  }
  ynew
}
library(fda.usc)
################################################################################
load("influenza.rda")
nt<-53*52*2
nrep<-156
a1<-209+ 52
a2<-a1+nrep-10
nmodel<-12
flu2<-flu

# selection response variable (n+1, to n+5)
ii1<-c(12,21,20,19)

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
a1<-1
a2<-a1+128

iispace<-3500
#iispace<-0
indspace<-flu$df$ispace>iispace

indspace2<-sum(unique(flu$df$ispace)>iispace)

lev<-names(table(flu$df$ispace))
indi3<-flu$df$ispace<16000000

######epidemic<-1:156 # in the paper, high time-consuming
## Moderate Time-consuming (only 17 rolling iterations (last 17 weeks)
p1<-2
q1<-0
p2<-1
q2<-0
ii1<-c(12,21,20,19)
## High Time-consuming
levcomar<-names(table(flu$df$ispace[indspace]))
#levcomar <- setdiff(levcomar,c(levcomar[23],"3207","3209", "1506"))
ncomar<-length(levcomar)
epidemic<-c(145:164+52,145:164+2*52) 
flu$irra$names<-list("main"= "RS" ,"xlab"=  "Day","ylab"=  "W/m2")
flu$irra<-completa.fdata(flu$irra)
iresp <- 25 # donde se guarda la ytest

n.ahead<-4
# repetir para cor.ARMA(1,0) o ARMA(2,0)
p1 <- 1
p2 <- 1 
par.AR<-list("cor.AR"=list("order.max"=4))
par.ARMA1<-list("cor.ARMA"=list("p"=1))
par.ARMA2<-list("cor.ARMA"=list("p"=2))


################################################################################
# Figuras
################################################################################
nt<-53*52*2
#nrep<-156
a1<-209+ 52
a2<-a1+nrep-10
nmodel<-12
## High Time-consuming
levcomar<-names(table(flu$df$ispace[indspace]))
#levcomar <- setdiff(levcomar,c(levcomar[23],"3207","3209", "1506"))
ncomar<-length(levcomar)
epidemic<-c(145:164+52,145:164+2*52) 

icomar<-2
nam.var<-c("Rate(w)","Temp(t)","Temp.thres(t)","SR(t)",
           "Rate(w),Temp(t)","Rate(w),Temp(t)",
           "Rate(w)Temp.thres(t)","Temp(t),SR(t)")
nam.mod<-c("FLM","GLS","iGLS")
ij<-1 #horizonete

j<-1
flu2<-flu
flu$df$y<-flu$df[,ii1[j]]

ind<-which(flu$df$itime>=(epidemic[1]+a1) & flu$df$itime<=616
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

#a plicar logartibmo en la tasa !!!!!
################# INICIO #####################
add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(1, 2,0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}
#####                     pdf(file = "S2_fig.pdf",width = 10.24, height =7.68, pointsize = 10)
# se guarda manualmente como pdf
windows(1000,700)
par(mfrow=c(1,2))
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
     ,col=cols,cex=.9,
     main=mn2,xlab="Time w: 13 weeks   ",ylab=ylab)

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
     ,col=cols,cex=.9,
     main=mn2,xlab="   Time t: 14 days",ylab=ylab)


leg <- c(expression(v[q[1]]), expression(v[q[2]]), expression(v[q[3]]) ,expression(v[q[4]]))
#legend(8,2500,inset=c(-0.2,0),col=cols,legend=leg,cex=.8,bty="n",lty=1,lwd=1)
add_legend("bottom", legend=leg,lty=1, #pch=20, 
           col=cols, lwd=2,
           horiz=TRUE, #◄bty='n',
           cex=0.95)
#########################dev.off()
################# fin Grafico #####################



################################################################################
################################################################################
# Figure 1      Series Temporales 
flu4<-aggregate(flu2$df,list(flu2$df[,"itime"]),mean)
#cambiar: utilizar un ggplot o algo así más chulo (ver links/linkedIN)
#png(file = "S1_fig.png")   
par(mfrow=c(4,1))
FluRate<-ts(exp(flu4[,"ratet0"]),2001,2011,53)  
Temperature<-ts(flu4[,"tmed"] ,2001,2011,53)
SolarRadiation<-ts(flu4[,"imed"] ,2001,2011,53)
RelativeHumidity<-ts(flu4[,"hmed"] ,2001,2011,53)
TemperatureThreshold<-ts(flu4[,"thold"] ,2001,2011,53)
aa<-cbind(FluRate,Temperature,SolarRadiation,RelativeHumidity)
plot(aa,main="",xlab="Years")
#dev.off()
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
# TABLE 6 PLOS paper
round(B[,c(8:12,2,1)],2)

#save.image(file="fluExample.RData")