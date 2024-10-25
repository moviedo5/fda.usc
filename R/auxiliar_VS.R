# funciones del Script
# dcor.y
# dist.list
# AICc
# GCV
# print.select.glm
# print.select.gsam
# print.select.bgsam
################################################################################
dcor.y<-function(ldist,response,n,bcdcor=TRUE){
  lenldist<-length(ldist)
  namldist<-names(ldist)
  if (missing(n)) {
    if (is.fdata(ldist[[1]]) || is.matrix(ldist[[1]]) || is.data.frame(ldist[[1]]))
      n<-nrow(ldist[[1]])
    if (is.vector(ldist[[1]]))    n <- length(response)
  }
  
  if (missing(response)) {print("All DC are computed")
    dst<-diag(lenldist)
    for (i in 1:(lenldist-1)) {
      for (j in i:(lenldist)) {
        if (bcdcor)     dst[i,j]<-dst[j,i]<- bcdcor.dist(ldist[[i]],ldist[[j]],n=n)
        else            dst[i,j]<-dst[j,i]<-cor.fdist(ldist[[i]],ldist[[j]])
      }}
    
    rownames(dst)<-colnames(dst)<-namldist
  }
  else{                     #se calculan todas las distancias respecto la respuesta
    if (length(response)>1) stop("Incorrect response label")
    ii<-response==namldist
    if (sum(ii)==0) stop("Incorrect response label")
#    iresponse<-which(ii) # No se usa
    dst<-numeric(lenldist-1)
    predictors<-setdiff(namldist,response)
    
    for (i in 1:(lenldist-1)) {
      #           dst[i]<-cor.fdist(ldist[[response]],ldist[[predictors[i]]])
      if (bcdcor)  dst[i]<-bcdcor.dist(ldist[[response]],ldist[[predictors[i]]],n=n)
      else dst[i]<-dcor.dist(ldist[[response]],ldist[[predictors[i]]])
    }
    names(dst)<-predictors
  }
  dst
}
################################################################################

####################
#best.pc.dcor<- function(data,"x","y"){
#  dd=dcor.xy(ldata$df[[i]],res,n=n) #edf
#  dcor[j,i]=dd$estimate*(dd$p.value < pvalor)
#}
####################
pvalue.anova<- function(model1,model2){
  anova(model1,model2,test="F")
}

# AICc<-function(model){
#   suma<-summary(model)          
#   k<- sum(model$edf)#sum(suma$edf) 
#   n<-suma$n
#   return(AIC(model)+(2*k*(n+1))/(n-k+1))
# }
# 
AIC2<-function(model){
  return(AIC(model))
}

##############################################
#  modelo$gcv.ubre
##############################################
#GCV<-function(model){
#  return(a1$gcv.ubre)
#}
##############################################
p.signif<-function(model){
  rev(summary(model)$s.table[,4])[1] 
}
##############################################
rsq<-function(model){
  return(1-summary(model)$r.sq)
}
##############################################
sp<-function(model){
  # cat("entra SP criterio")
  aa<-summary(model)$sp.criterion
#  print(aa)
  return(aa)
  
}
#################################################
AICc<-function(model){
  suma<-summary(model)
  n<- length(model$residuals)
  if (any(class(model)=="glm")){
    k<- suma$df[1]
  }
  if (any(class(model)=="lm")){
    k<- suma$df[1]
  }
  if (any(class(model)=="fregre.glm")){
    k<- suma$df[1]
  }
  if (any(class(model)=="fregre.lm")){
    k<- suma$df[1]
  }
  if (any(class(model)=="fregre.gsam")){
    k<- sum(model$edf)#sum(suma$edf) 
   # print(n)
  }
  aux<-AIC(model)+(2*k*(n+1))/(n-k+1)
  #cat(k,n,AIC(model),aux," AICCCCCCCCCCCCCCCCCC  \n")
  return(aux)
}
#################################################
GCV<-function(model){
  model$gcv.ubre
}
#################################################

print_select_glm <- function(x,...){
  cat("\n Model selected:\n")
  print(x$model)
  cat("\n Stepwise GoF:\n")
  print(x$gof)
  cat("\n Predictor selected:\n")
  print(x$i.predictor)  
}


print_select_gsam <- function(x,...){
  cat("\n Model selected:\n")
  print(x$model)
  cat("\n Stepwise GoF:\n")
  print(x$gof)
  cat("\n Predictor selected:\n")
  print(x$i.predictor)  
}

# print.select.bgsam <- function(x,...){
#   #cat("\n Model selected:\n")
#   #print(x[[1]])
#   #cat("\n Stepwise GoF:\n")
#   #print(x[[2]])
#   cat("\n Predictor selected (by order):\n")
#   print(x[[4]])
#   cat("\n Model type selected (by order):\n")
#   summary(m3[[1]])
# }

#################################################

print_fregre_glm_vs<- function(x,...){
  cat("\n Model selected:\n")
  print(x$model)
  cat("\n Stepwise GoF:\n")
  print(x$gof)
  cat("\n Predictor selected:\n")
  print(x$i.predictor)  
}

print_fregre_gsam_vs <- function(x,...){
  cat("\n Model selected:\n")
  print(x$model)
  cat("\n Stepwise GoF:\n")
  print(x$gof)
  cat("\n Predictor selected:\n")
  print(x$i.predictor)  
}
