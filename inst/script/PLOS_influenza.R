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
a1<-1
a2<-a1+149

nrep<-500
indi<-1:nrep
iispace<-3500
#iispace<-0
indspace<-flu$df$ispace>iispace

indspace2<-sum(unique(flu$df$ispace)>iispace)
pred.res<-array(NA,dim=c(nahead,nrep,nmodel+1,ncomar)) #nmodel +1 (la y)
pred0<-pred<-array(NA,dim=c(nahead,nrep,nmodel))

lev<-names(table(flu$df$ispace))
indi3<-flu$df$ispace<16000000

# 10 o 53 comarcas
nrep<-500

######epidemic<-1:156 # in the paper, high time-consuming
## Moderate Time-consuming (only 17 rolling iterations (last 17 weeks)
p1<-2
q1<-0

p2<-1
q2<-0
source(".\\flm_ar.R")

## High Time-consuming
levcomar<-names(table(flu$df$ispace[indspace]))
#levcomar <- setdiff(levcomar,c(levcomar[23],"3207","3209", "1506"))
ncomar<-length(levcomar)
epidemic<-c(145:164+52,145:164+2*52) 
flu$irra$names<-list("main"= "RS" ,"xlab"=  "Day","ylab"=  "W/m2")
#flu$irra<-completa.fdata(flu$irra)
iresp <- 25 # donde se guarda la ytest

n.ahead<-4
par.AR<-list("cor.AR"=list("order.max"=4))
par.ARMA1<-list("cor.ARMA"=list("p"=1))
par.ARMA2<-list("cor.ARMA"=list("p"=2))
#########################################
epidemic<-c(145:172+2*52) -20
#epidemic<-c(145:148) 
par.cor <- par.ARMA1
ctrl <- list("p"=1)
# estimacion del modelo p

#♣epidemic <-249:255
p.est <- array(NA,dim=c(length(epidemic),ncomar,8))
# estimation/prediction by rolling using AR1 dependence error
source("influenza_prediction.R")
#save.image(file="fluExampleAR1.RData")

table6<-t(round(apply(pred[,epidemic,1:24],c(1,3),mean,na.rm=T),4))
table7<-cbind(table6[c(1:8),1],table6[c(9:16),1],
              table6[c(9:16)+8,1],table6[c(1:8),2],
              table6[c(9:16),2],table6[c(9:16)+8,2])
colnames(table7)<-c("FLM (t+1)","FiGLS (t+1)","FGLS (t+1)","FLM (t+2)","FiGLS (t+2)","FGLS (t+2)")
rownames(table7)<-c("Rate(w)","Temp(t)","Temp.thres(t)","SR(t)",
                    "Rate(w),Temp(t)","Rate(w),Temp(t)",
                    "Rate(w)Temp.thres(t)","Temp(t),SR(t)")
tabla1AR_1<-table7
#########################################
par.cor <- par.ARMA2
ctrl <- list("p"=2)
# estimation/prediction by rolling using AR2 dependence error
source("influenza_prediction.R")
table6<-t(round(apply(pred[,epidemic,1:24],c(1,3),mean,na.rm=T),4))
table7<-cbind(table6[c(1:8),1],table6[c(9:16),1],
              table6[c(9:16)+8,1],table6[c(1:8),2],
              table6[c(9:16),2],table6[c(9:16)+8,2])
colnames(table7)<-c("FLM (t+1)","FiGLS (t+1)","FGLS (t+1)","FLM (t+2)","FiGLS (t+2)","FGLS (t+2)")
rownames(table7)<-c("Rate(w)","Temp(t)","Temp.thres(t)","SR(t)",
                    "Rate(w),Temp(t)","Rate(w),Temp(t)",
                    "Rate(w)Temp.thres(t)","Temp(t),SR(t)")
tabla1AR_2<-table7
#save.image(file="fluExampleAR2.RData")
#########################################
par.cor <- par.AR
#ctrl <- list("order.max"=NULL)
ctrl <- list("order.max"=4)
# estimation/prediction by rolling using AR(p) dependence error
source("influenza_prediction.R")
#save.image(file="fluExampleARp.RData")
table6<-t(round(apply(pred[,epidemic,1:24],c(1,3),mean,na.rm=T),4))
table7<-cbind(table6[c(1:8),1],table6[c(9:16),1],
              table6[c(9:16)+8,1],table6[c(1:8),2],
              table6[c(9:16),2],table6[c(9:16)+8,2])
colnames(table7)<-c("FLM (t+1)","FiGLS (t+1)","FGLS (t+1)","FLM (t+2)","FiGLS (t+2)","FGLS (t+2)")
rownames(table7)<-c("Rate(w)","Temp(t)","Temp.thres(t)","SR(t)",
                    "Rate(w),Temp(t)","Rate(w),Temp(t)",
                    "Rate(w)Temp.thres(t)","Temp(t),SR(t)")
tabla1AR_p<-table7
#save.image(file="fluExampleAR_1_2_p.RData") 
arp.est <- p.est
###############################################################################

############### Results for t+1 ###############
table_t1<-cbind(tabla1AR_1[,1],tabla1AR_1[,3]
                ,tabla1AR_1[,2],tabla1AR_2[,2]
                ,tabla1AR_p[,2])
colnames(table_t1)<-c("FLM","FGLS-AR(1)","FiGLS-AR(1)",
                      "FiGLS-AR(2)","FiGLS-AR(p)")
rownames(table_t1)<-c("Rate(w)","Temp(t)","Temp.thres(t)","SR(t)",
                      "Rate(w),Temp(t)","Rate(w),Temp(t)",
                      "Rate(w)Temp.thres(t)","Temp(t),SR(t)")
round(table_t1,3)
############### Results for t+2 ###############
table_t2<-cbind(tabla1AR_1[,4],tabla1AR_1[,6]
              ,tabla1AR_1[,5],tabla1AR_2[,5]
              ,tabla1AR_p[,5])
colnames(table_t2)<-c("FLM","FGLS-AR(1)","FiGLS-AR(1)",
                    "FiGLS-AR(2)","FiGLS-AR(p)")
rownames(table_t2)<-c("Rate(w)","Temp(t)","Temp.thres(t)","SR(t)",
                    "Rate(w),Temp(t)","Rate(w),Temp(t)",
                    "Rate(w)Temp.thres(t)","Temp(t),SR(t)")
round(table_t2,3)
library(xtable)
xtable(table_t1,digits=3)
xtable(table_t2,digits=3)
###########################################################
###########################################################
# round(apply(arp.est,2:3,mean))
