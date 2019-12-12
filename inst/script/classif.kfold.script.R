
# cambio en classif.glm, pendiente de hacerlo en gsam, tree, ML, etc.
#newy <- y <- data$df[,response]
#out4 <- classif.glm(form, data=ldata, basis.x=basis.pc)


data(tecator)    
cutpoint <- 15
cutpoint <- 18
tecator$y$class <-factor(ifelse(tecator$y$Fat<cutpoint,0,1))
table(tecator$y$class )
# todos los conjuntos de datos tendrian esta estructura:
# df con la variable class
# x  con la variable funcional
x<-tecator[[1]]
# x1  derivada
x1<-fdata.deriv(tecator[[1]])
x1<-rproc2fdata(nrow(x),1:50)
#x2<-rproc2fdata(nrow(x),x$argvals,mu=func.mean(x)$data[1,])
x2<-fdata.deriv(tecator[[1]],2)
plot(x2)
plot(x1)
data<- list("df"=tecator$y,x=x,x1=x1,x2=x2)
ytrain <- data$df[,"class"]
#covar<-"x"
response <- "class"
par.classif <- list()
#source("./R/create.fdata.basis.R")
########################################
formula<- formula(class~x)
formula<- formula(class~x+x1+x2)

# PRUEBA 1: sin especificar el parametro param.kfold
# se utiliza la opcion por defecto del classificador


classif="classif.glm";
output0<-classif.kfold(formula, data, classif = classif,
                       par.classif=list("weights"="inverse"),
                       kfold = kfold,verbose=F,models=T)
summary(output0)

formula<- formula(class~x+x1+x2)
classif="classif.glm";
output0<-classif.kfold(formula, data, classif = classif,
                       kfold = kfold,verbose=T,models=T)
traceback()
summary(output0)
table(predict.classif(output0),output0$group.est)
traceback()



#formula<- formula(class~s(x,k=3)+s(x1,k=3))

param.kfold <- list("x"=list("pc"=c(1:5))
                    ,"x1"=list("pc"=c(1,3,5,9))
                    ,"x2"=list("pc"=c(1,3,5,9)))

param.kfold <- list("x"=list("fourier"=c(5,7)),
                    "x1"=list("fourier"=c(5,7,9)),
                    "x2"=list("fourier"=c(5,7,9)))

set.seed(1)
par.classif = list()
kfold = 4
verbose=T

output4<-classif.kfold(formula, data, 
                       classif = classif,
                       par.classif = par.classif,
                       kfold = kfold,
                       param.kfold = param.kfold,
       #                measure="accuracy",
                       verbose=T)
summary(output4)
output4$params.error

output<-output4
summary(output)
summary(output$fit[[1]])
#output$fit[[1]]
length(output$list.pred)
which.min(output$error)
output$param.min
pred<-predict.classif(output,data)
table(output$group.est,pred)
table(pred,data$df$class)
table(output$group.est,data$df$class)

#colMeans(sweep(output$pred.df,1, ytrain,"=="))
#which.min(colMeans(sweep(output$pred.df,1, ytrain,"==")))


#output$group.est.kfold
table(output$group.est,data$df$class)
table(output$group.est,data$df$class)
########################################

# Prueba 4 
formula<- formula(class~x+x2)
classif="gkam"


param.kfold <- list("x"=list("h"=c(3,5,9)),
                    "x2"=list("h"=c(3,5,9)))

par.np<-list("x"=list(Ker=Ker.epa,type.S="S.NW","h"=3),
             "x2"=list(Ker=Ker.epa,type.S="S.NW","h"=3))

par.np<-list("x"=list(Ker=Ker.unif,type.S="S.KNN","h"=3),
             "x2"=list(Ker=Ker.unif,type.S="S.KNN","h"=3))

par.metric=list("x"=list(metric=semimetric.deriv,nderiv=1,nbasis=5)
                ,"x2"=list(metric=semimetric.deriv,nderiv=1,nbasis=5))


par.control <- list(maxit = 2)

#################################################
#################################################
set.seed(1)
par.classif<-list("par.np"=par.np,
                  "par.metric"=par.metric, # no va
                  "control" = par.control)

output<-classif.kfold (formula, data, classif=classif,
                       par.classif=par.classif,
                       kfold = kfold,
                       param.kfold=param.kfold)
which.min(output$error)
(output$error)
min(output$error)
attributes(output$fit[[1]]$result$x2$par.S)
#################################################
> which.min(output$error)param 9-5 6 
> min(output$error)[1] 0.04093567
#################################################


#output<-output4
summary(output)
output$fit[[1]]
#summary(output$fit[[1]])
length(output$list.pred)
(output$error)


which.min(output$error)
output$param.min
pred<-predict.classif(output,data)
table(output$group.est,pred)
table(pred,data$df$class)
table(output$group.est,data$df$class)

args(classif.gkam)
out.gkam<-classif.gkam(formula, data=data,control=par.control)

out.gkam<-classif.gkam(formula, data=data,control=par.control
                       ,par.np=par.np
                       ,par.metric = par.metric)
out.gkam
summary(out.gkam)
formula
out.gkam<-classif.gkam(formula, data=data,control=par.control)

pred<-predict.classif(out.gkam,data)
pred

###########################################
data(phoneme)
names(phoneme)
names(phoneme$learn)
unlist(lapply(phoneme,class))
dim(phoneme$learn)
# Factor levels: "sh" 1, "iy" 2, "dcl" 3, "aa" 4 and "ao" 5.
table(phoneme$classlearn)
par(mfrow=c(1,2))
plot(phoneme[["learn"]],col=phoneme[["classlearn"]])
plot(phoneme[["test"]],col=phoneme[["classtest"]])
data<-list("df"=data.frame(class=phoneme[["classlearn"]]),"x"=phoneme$learn)
dim(data$df)


group<-as.numeric(phoneme[["classlearn"]][181:220])
x<-phoneme$learn[181:220]
res.np1<-classif.np(group,x)
summary(res.np1)
res.np2<-classif.np(group,x,metric=metric.DTW)
summary(res.np2)
res.np3<-classif.np(group,x,metric=metric.TWED)
summary(res.np3)
res.np4<-classif.np(group,x,metric=metric.WDTW)
summary(res.np4)

########################################
# PRUEBA 1 
formula<- formula(class~x);  classif="classif.glm";
lpc<-list("x"=create.pc.basis(data$x,1:9))
#lpc<-list("x"=create.fdata.basis(data$x,1:9,type.basis = "pc"))
res<-classif.glm(formula, data=data,basis.x=lpc)
res2<-classif.glm(formula, data=data,basis.x=lpc,type="majority")
summary(res)
summary(res2)
args(classif.glm)

param.kfold <- list("x"=list("pc"=c(1:9))
                    ,"x1"=list("pc"=c(1,3,5,9))
                    ,"x2"=list("pc"=c(1,3,5,9)))

param.kfold <- list("x"=list("fourier"=c(5,7)),
                    "x1"=list("fourier"=c(5,7,9)),
                    "x2"=list("fourier"=c(5,7,9)))
kfold <- 10
output<-classif.kfold (formula, data, classif=classif,
                       par.classif=par.classif,
                       kfold = kfold,
                       param.kfold=param.kfold,
                       verbose=T)
output$max.prob
#output el call no ese?e todo el texto

summary(output$fit[[1]])

(output$model$fit[[1]]$coefficients)-(res$fit[[1]]$coefficients)
summary(output$model)
summary(output$model$fit[[1]])
summary(res$fit[[1]])

summary(output$model)


round(colMeans(output$error)*100,2)
output$param.min
plot(round(colMeans(output$error)*100,2))
summary(output$model)
summary(output$model$fit[[1]])    

###############################
 data(tecator)    
 cutpoint <- 15
 cutpoint <- 20
 tecator$y$class<-factor(ifelse(tecator$y$Fat<cutpoint,0,1))
 
 # todos los conjuntos de datos tendr?n esta estructura:
 # df con la variable class
 # x  con la variable funcional
 data<- list("df"=tecator$y,x=tecator[[1]])
 covar<-"x"
 response <- "class"

 par.classif <- list()
 kfold <- 10
  
 
 
 ########################################
 # PRUEBA 1 
 formula<- formula(class~x);  classif="classif.glm";
 param.kfold <- list("fourier"=c(5,7,9,11,13))
 param.kfold <- list("bspline"=c(5,7,9,11,13))
 param.kfold <- list("pc"=c(5,7,9,11,13))
 kfold <- 10
 output<-classif.kfold (formula, data, classif=classif,
                        par.classif=par.classif,
                        kfold = kfold,
                        param.kfold=param.kfold,
                        verbose=T)
 output$model$prob.group
 summary(output$model$fit[[1]])
 lpc<-list("x"=create.pc.basis(data$x,1:9))
 lpc<-list("x"=create.fdata.basis(data$x,1:9,type.basis = "pc"))
 res<-fda.usc::classif.glm(formula, data=data,basis.x=lpc)
 res2<-classif.glm(formula, data=data,basis.x=lpc,type="majority")
 summary(res)
 summary(res2)
 
 (output$model$fit[[1]]$coefficients)-(res$fit[[1]]$coefficients)
 summary(output$model)
 summary(output$model$fit[[1]])
 summary(res$fit[[1]])
 
 summary(output$model)
 
 
 round(colMeans(output$error)*100,2)
 output$param.min
 plot(round(colMeans(output$error)*100,2))
 summary(output$model)
 summary(output$model$fit[[1]])    
 
 
 # Prueba 2: ML c("tree", "svm", "nnet", "randomForest")
 formula<- formula(class~x);  classif="randomForest"
 param.kfold <- list("basis.x"=1:14)
 ########################################
 # Prueba 3
 # Model has more coefficients than data
 kbs <- -1
 formula<- formula(class~s(x,k=kbs));  classif="gsam"; 
 param.kfold <- list("basis.x"=1:14)
 param.kfold <- list("basis.x"=3:10)
 ###############################
 # Prueba 4 
 formula<- formula(class~x);  classif="gkam";
 classif="gkam"
 param.kfold <- list("h"=3:5)
 #param.kfold <- list("h"=1:10)
 #param.kfold <- list("h"=seq(3,15,by=2))
 par.np<-list("x"=list(Ker=Ker.unif,type.S="S.KNN","h"=3))
 par.metric=list("x"=list(metric=semimetric.deriv,nderiv=1,nbasis=5))
 par.control <- list(maxit = 2)
 par.classif<-list("par.np"=par.np,par.metric=par.metric,control = par.control)
 ########################################
 output<-classif.kfold (formula, data, classif=classif,
                        par.classif=par.classif,
                        kfold = kfold,
                        param.kfold=param.kfold)
 round(colMeans(output$error)*100,2)
 output$param.min
 plot(round(colMeans(output$error)*100,2))
 summary(output$model)    

 
 

?create.fdata.basis 
 
 ########################################
 # PRUEBA 1 SELECCION BASES  
 formula<- formula(class~x);  classif="classif.glm";
 param.kfold <- list("basis.PC"=3:14)
 kfold <- 10
 output<-classif.kfold (formula, data, classif=classif,
                        par.classif=par.classif,
                        kfold = kfold,
                        param.kfold=param.kfold,verbose=T)
 # Prueba 2: ML c("tree", "svm", "nnet", "randomForest")
 formula<- formula(class~x);  classif="randomForest"
 param.kfold <- list("basis.x"=1:14)
 ########################################
 output<-classif.kfold (formula, data, classif=classif,
                        par.classif=par.classif,
                        kfold = kfold,
                        param.kfold=param.kfold)
 round(colMeans(output$error)*100,2)
 output$param.min
 plot(round(colMeans(output$error)*100,2))
 summary(output$model)    
 
 










 #############pruebas para borrar
 param.kfold <- list("x2"=list("h"=h))
 param.kfold<-NULL
 classif<-"classif.gkam"
 output3b <- classif.kfold (formula, data, classif=classif,
                            kfold = 4, param.kfold=param.kfold,verbose=T,models=T)
 weights4class<-fda.usc:::weights4class
 
 # if (length(weights) == ndatos) {
 #   wtab <- tapply(weights, y, sum)
 #   ii <- wtab/sum(wtab)
 #   names(ii) <- lev
 #   class.weights <- NULL
 # }
 y<-data$df$class
 table(y)
 set.seed(1:6)
 output3a <- classif.kfold (formula, data, classif=classif,
                            kfold = 4,verbose=T,models=T)
 
 set.seed(1:6)
 output3b <- classif.kfold (formula, data, classif=classif,par.classif=list("weights"="inverse"),
                            kfold = 4,verbose=T,models=T)
 set.seed(1:6)
 wei<-weights4class(y,"inverse")
 coste<-tapply(wei,y,mean)
 output3c <- classif.kfold (formula, data, classif=classif,
                            kfold = 4,verbose=T,models=T,measure="cost",cost=coste)
 
 set.seed(1:6)
 coste<-tapply(weights4class(y,"inverse"),y,mean)
 output3d <- classif.kfold (formula, data, classif=classif,par.classif=list("weights"="inverse"),
                            kfold = 4,verbose=T,models=T,measure="cost",cost=coste)
 
 output3a$prob.classification;output3b$prob.classification;output3c$prob.classification;output3d$prob.classification
 
 table(output3a$group.est,y);table(output3b$group.est,y);table(output3c$group.est,y);table(output3d$group.est,y)
 table(output3a$pred.kfold,y);table(output3b$pred.kfold,y);table(output3c$pred.kfold,y);table(output3d$pred.kfold,y)
 
 output3b$param.min
 output3a$params.error
 
 summary(output3b)
 for (i in 1:kfold) print(output3a$models[[i]]$fit[[1]]$result[[1]]$h.opt)
 for (i in 1:kfold) print(output3b$models[[i]]$fit[[1]]$result[[1]]$h.opt)
 
 
 #########
 classif<-"classif.np"
 wei<-weights4class(y,"inverse")
 
 set.seed(1:6)
 output3a <- classif.kfold (formula, data, classif=classif,
                            kfold = 4,verbose=T,models=T)
 
 set.seed(1:6)
 output3b <- classif.kfold (formula, data, classif=classif,par.classif=list(par.S=list("w"=wei)),
                            kfold = 4,verbose=T,models=T)
 
 set.seed(1:6)
 coste<-tapply(weights4class(y,"inverse"),y,mean)
 output3c <- classif.kfold (formula, data, classif=classif,
                            kfold = 4,verbose=T,models=T,measure="cost",cost=coste)
 
 set.seed(1:6)
 coste<-tapply(weights4class(y,"inverse"),y,mean)
 output3d <- classif.kfold (formula, data, classif=classif,par.classif=list(par.S=list("w"=wei)),
                            kfold = 4,verbose=T,models=T,measure="cost",cost=coste)
 
 output3a$prob.classification;output3b$prob.classification;output3c$prob.classification;output3d$prob.classification
 
 table(output3a$group.est,y);table(output3b$group.est,y);table(output3c$group.est,y);table(output3d$group.est,y)
 table(output3a$pred.kfold,y);table(output3b$pred.kfold,y);table(output3c$pred.kfold,y);table(output3d$pred.kfold,y)
 
 output3a$params.error
 
 summary(output3b)
 for (i in 1:kfold) print(output3a$models[[i]]$fit[[1]]$result[[1]]$h.opt)
 for (i in 1:kfold) print(output3b$models[[i]]$fit[[1]]$result[[1]]$h.opt)
 
 
 
 
 args(classif.np)
 output3c <- classif.np(data$df$class,data$x2)
 output3b$max.prob;output3c$max.prob
 
 
 param.kfold <- list("x2"=list("h"=h))
 param.kfold<-NULL
 output3b <- classif.kfold (formula, data, classif="classif.gkam",
                            kfold = 4, param.kfold=param.kfold,verbose=T)
 output3b$param.min
 output3b$params.error
 
 output3d <- classif.gkam (formula, data, classif="classif.gkam")
 output3b$max.prob;output3d$max.prob
 
 y<-data$df$class