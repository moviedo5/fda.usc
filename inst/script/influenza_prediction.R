for (i in  epidemic) {#
  for (icomar in 1:ncomar) {  
    j<-1
    flu2<-flu
    flu$df$y<-flu$df[,ii1[j]]
    #ind<-which(flu$df$itime>=(epidemic[1]) & flu$df$itime<=(i+a2)
    ind<-which(flu$df$itime>=(i) & flu$df$itime<=(i+a2-1)
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
    flu2$df$itime
    
    flu2$irra<-flu$irra[ind]
    flu2$irra<-completa.fdata(flu2$irra)
    
    flu2$humt10<-flu$humt10[ind,]
    res.lm1<-fregre.lm(ff1,data=flu2)
    res.lm2<-fregre.lm(ff2,data=flu2)
    res.lm3<-fregre.lm(ff3,data=flu2)
    res.lm4<-fregre.lm(ff4,data=flu2)
    res.lm5<-fregre.lm(ff5,data=flu2)
    res.lm6<-fregre.lm(ff6,data=flu2)
    res.lm7<-fregre.lm(ff7,data=flu2)
    res.lm8<-fregre.lm(ff8,data=flu2)
    traza<-fda.usc:::traza
    #par.cor<-list("cor.ARMA"=list("index"="itime","p"=1))
    res.igls1<-fregre.igls(ff1,data=flu2,correlation=par.cor,control=ctrl)
    res.igls2<-fregre.igls(ff2,data=flu2,correlation=par.cor,control=ctrl)          
    res.igls3<-fregre.igls(ff3,data=flu2,correlation=par.cor,control=ctrl)        
    res.igls4<-fregre.igls(ff4,data=flu2,correlation=par.cor,control=ctrl)        
    res.igls5<-fregre.igls(ff5,data=flu2,correlation=par.cor,control=ctrl)     
    res.igls6<-fregre.igls(ff6,data=flu2,correlation=par.cor,control=ctrl)
    res.igls7<-fregre.igls(ff7,data=flu2,correlation=par.cor,control=ctrl)       
    res.igls8<-fregre.igls(ff8,data=flu2,correlation=par.cor,control=ctrl)       
    iepi<-i-epidemic[1]+1
    p.est[iepi,icomar,1]<- ifelse(is.null(res.igls1$corStruct$ar$arma[1]),0,res.igls1$corStruct$ar$arma[1])
    p.est[iepi,icomar,2]<- ifelse(is.null(res.igls1$corStruct$ar$arma[1]),0,res.igls2$corStruct$ar$arma[1])
    p.est[iepi,icomar,3]<- ifelse(is.null(res.igls3$corStruct$ar$arma[1]),0,res.igls3$corStruct$ar$arma[1])
    p.est[iepi,icomar,4]<- ifelse(is.null(res.igls4$corStruct$ar$arma[1]),0,res.igls4$corStruct$ar$arma[1])
    p.est[iepi,icomar,5]<- ifelse(is.null(res.igls5$corStruct$ar$arma[1]),0,res.igls5$corStruct$ar$arma[1])
    p.est[iepi,icomar,6]<- ifelse(is.null(res.igls6$corStruct$ar$arma[1]),0,res.igls6$corStruct$ar$arma[1])
    p.est[iepi,icomar,7]<- ifelse(is.null(res.igls7$corStruct$ar$arma[1]),0,res.igls7$corStruct$ar$arma[1])
    p.est[iepi,icomar,8]<- ifelse(is.null(res.igls8$corStruct$ar$arma[1]),0,res.igls8$corStruct$ar$arma[1])
    
    #hay missings en irra!
    #lpc<-list("ratet13"=create.pc.basis(flu2$ratet13,1:4))
    #,"irra"=create.pc.basis(flu2[["irra"]],1:2))
    #res.gls1<-fregre.gls(ff1,data=flu2,correlation=corARMA(p=p2,q=q2,form=~itime),method="ML")       
    res.gls1<-fregre.gls(ff1,data=flu2,correlation=corARMA(p=p2,q=q2),method="ML")       
    res.gls2<-fregre.gls(ff2,data=flu2,correlation=corARMA(p=p2,q=q2),method="ML")       
    res.gls3<-fregre.gls(ff3,data=flu2,correlation=corARMA(p=p2,q=q2),method="ML")       
    res.gls4<-fregre.gls(ff4,data=flu2,correlation=corARMA(p=p2,q=q2),method="ML")          
    res.gls5<-fregre.gls(ff5,data=flu2,correlation=corARMA(p=p2,q=q2),method="ML")       
    res.gls6<-fregre.gls(ff2,data=flu2,correlation=corARMA(p=p2,q=q2),method="ML")      
    res.gls7<-fregre.gls(ff7,data=flu2,correlation=corARMA(p=p2,q=q2),method="ML")      
    res.gls8<-fregre.gls(ff8,data=flu2,correlation=corARMA(p=p2,q=q2),method="ML")      
    
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
    
    pred1<-predict.fregre.igls(res.igls1,flu3,flu3$df,n.ahead=ahead)
    newd<-rbind(flu3$df,flu3$df,flu3$df,flu3$df) 
    #newd$itime<-flu3$df$itime:(flu3$df$itime+ahead-1)
    
    pred1<-predict.fregre.igls(res.igls1,flu3,newd,n.ahead=ahead)
    pred2<-predict.fregre.igls(res.igls2,flu3,newd,n.ahead=ahead)
    pred3<-predict.fregre.igls(res.igls3,flu3,newd,n.ahead=ahead)
    pred4<-predict.fregre.igls(res.igls4,flu3,newd,n.ahead=ahead)  
    pred5<-predict.fregre.igls(res.igls5,flu3,newd,n.ahead=ahead)
    pred6<-predict.fregre.igls(res.igls6,flu3,newd,n.ahead=ahead)
    pred7<-predict.fregre.igls(res.igls7,flu3,newd,n.ahead=ahead)
    pred8<-predict.fregre.igls(res.igls8,flu3,newd,n.ahead=ahead)
    
    pred.gls1<-predict.fregre.gls(res.gls1,flu3,newd,n.ahead=ahead)
    pred.gls2<-predict.fregre.gls(res.gls2,flu3,newd,n.ahead=ahead)
    pred.gls3<-predict.fregre.gls(res.gls3,flu3,newd,n.ahead=ahead)
    pred.gls4<-predict.fregre.gls(res.gls4,flu3,newd,n.ahead=ahead)  
    pred.gls5<-predict.fregre.gls(res.gls5,flu3,newd,n.ahead=ahead)
    pred.gls6<-predict.fregre.gls(res.gls6,flu3,newd,n.ahead=ahead)
    pred.gls7<-predict.fregre.gls(res.gls7,flu3,newd,n.ahead=ahead)
    pred.gls8<-predict.fregre.gls(res.gls8,flu3,newd,n.ahead=ahead)
    
    yresp<-t(flu3$df[,c("ratet1","ratet2","ratet3","ratet4")])
    
    pred.res[,i,1,icomar]<-rbind(pred.lm1,pred.lm1,pred.lm1,pred.lm1)
    pred.res[,i,2,icomar]<-rbind(pred.lm2,pred.lm2,pred.lm2,pred.lm2)
    pred.res[,i,3,icomar]<-rbind(pred.lm3,pred.lm3,pred.lm3,pred.lm3)
    pred.res[,i,4,icomar]<-rbind(pred.lm4,pred.lm4,pred.lm4,pred.lm4)
    pred.res[,i,5,icomar]<-rbind(pred.lm5,pred.lm5,pred.lm5,pred.lm5)
    pred.res[,i,6,icomar]<-rbind(pred.lm6,pred.lm6,pred.lm6,pred.lm6)  
    pred.res[,i,7,icomar]<-rbind(pred.lm7,pred.lm7,pred.lm7,pred.lm7)  
    pred.res[,i,8,icomar]<-rbind(pred.lm8,pred.lm8,pred.lm8,pred.lm8)  
    pred.res[,i,9,icomar]<-(pred1)
    pred.res[,i,10,icomar]<-t(pred2 )  
    pred.res[,i,11,icomar]<-t(pred3)
    pred.res[,i,12,icomar]<-t(pred4)  
    pred.res[,i,13,icomar]<-t(pred5)
    pred.res[,i,14,icomar]<-t(pred6 )   
    pred.res[,i,15,icomar]<-t(pred7 )   
    pred.res[,i,16,icomar]<-t(pred8 )
    pred.res[,i,17,icomar]<-t(pred.gls1)
    pred.res[,i,18,icomar]<-t(pred.gls2 )  
    pred.res[,i,19,icomar]<-t(pred.gls3)
    pred.res[,i,20,icomar]<-t(pred.gls4)  
    pred.res[,i,21,icomar]<-t(pred.gls5)
    pred.res[,i,22,icomar]<-t(pred.gls6 )   
    pred.res[,i,23,icomar]<-t(pred.gls7 )   
    pred.res[,i,24,icomar]<-t(pred.gls8 )
    pred.res[,i,iresp,icomar]<-yresp
  }
  
  pred[,i,1]<-apply((pred.res[,i,iresp,]-pred.res[,i,1,])^2,1,mean,na.rm=T)
  pred[,i,2]<-apply((pred.res[,i,iresp,]-pred.res[,i,2,])^2,1,mean,na.rm=T)
  pred[,i,3]<-apply((pred.res[,i,iresp,]-pred.res[,i,3,])^2,1,mean,na.rm=T)
  pred[,i,4]<-apply((pred.res[,i,iresp,]-pred.res[,i,4,])^2,1,mean,na.rm=T)
  pred[,i,5]<-apply((pred.res[,i,iresp,]-pred.res[,i,5,])^2,1,mean,na.rm=T)
  pred[,i,6]<-apply((pred.res[,i,iresp,]-pred.res[,i,6,])^2,1,mean,na.rm=T)
  pred[,i,7]<-apply((pred.res[,i,iresp,]-pred.res[,i,7,])^2,1,mean,na.rm=T)
  pred[,i,8]<-apply((pred.res[,i,iresp,]-pred.res[,i,8,])^2,1,mean,na.rm=T)
  pred[,i,9]<-apply((pred.res[,i,iresp,]-pred.res[,i,9,])^2,1,mean,na.rm=T)
  pred[,i,10]<-apply((pred.res[,i,iresp,]-pred.res[,i,10,])^2,1,mean,na.rm=T)
  pred[,i,11]<-apply((pred.res[,i,iresp,]-pred.res[,i,11,])^2,1,mean,na.rm=T)
  pred[,i,12]<-apply((pred.res[,i,iresp,]-pred.res[,i,12,])^2,1,mean,na.rm=T)
  pred[,i,13]<-apply((pred.res[,i,iresp,]-pred.res[,i,13,])^2,1,mean,na.rm=T)
  pred[,i,14]<-apply((pred.res[,i,iresp,]-pred.res[,i,14,])^2,1,mean,na.rm=T)
  pred[,i,15]<-apply((pred.res[,i,iresp,]-pred.res[,i,15,])^2,1,mean,na.rm=T)
  pred[,i,16]<-apply((pred.res[,i,iresp,]-pred.res[,i,16,])^2,1,mean,na.rm=T)
  pred[,i,17]<-apply((pred.res[,i,iresp,]-pred.res[,i,17,])^2,1,mean,na.rm=T)
  pred[,i,18]<-apply((pred.res[,i,iresp,]-pred.res[,i,18,])^2,1,mean,na.rm=T)
  pred[,i,19]<-apply((pred.res[,i,iresp,]-pred.res[,i,19,])^2,1,mean,na.rm=T)
  pred[,i,20]<-apply((pred.res[,i,iresp,]-pred.res[,i,20,])^2,1,mean,na.rm=T)
  pred[,i,21]<-apply((pred.res[,i,iresp,]-pred.res[,i,21,])^2,1,mean,na.rm=T)
  pred[,i,22]<-apply((pred.res[,i,iresp,]-pred.res[,i,22,])^2,1,mean,na.rm=T)
  pred[,i,23]<-apply((pred.res[,i,iresp,]-pred.res[,i,23,])^2,1,mean,na.rm=T)
  pred[,i,24]<-apply((pred.res[,i,iresp,]-pred.res[,i,24,])^2,1,mean,na.rm=T)
  print( i)
  pred[,i,]
}