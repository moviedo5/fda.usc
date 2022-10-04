#' @title Predicts from a fitted classif object.
#' 
#' @description Classifier of functional data by kernel method using functional data object
#' of class \code{classif}.  Returns the predicted classes using a previously trained model.
#' 
#' @param object Object \code{object} estimated by: k nearest neighbors method
#' \code{classif.knn}, kernel method \code{classif.kernel}.
#' @param new.fdataobj New functional explanatory data of \code{fdata} class.
#' @param type Type of prediction ("class or probability of each group
#' membership").
#' @param \dots Further arguments passed to or from other methods.

#' @return If type="class", produces a vector of predictions.
#' If type="probs", a list with the following components is returned: 
#' \itemize{
#' \item \code{group.pred} the vector of predictions. 
#' \item \code{prob.group} the matrix of predicted probability by factor level. 
#' }
#' 
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente \email{manuel.oviedo@@udc.es}
#' 
#' @seealso See also \code{\link{classif.np}} \code{\link{classif.glm}},
#' \code{\link{classif.gsam}} and \code{\link{classif.gkam}} .
#' 
#' @references Ferraty, F. and Vieu, P. (2006). \emph{Nonparametricc functional
#' data analysis.} Springer Series in Statistics, New York.
#' 
#' Ramsay, James O., and Silverman, Bernard W. (2006), \emph{ Functional Data
#' Analysis}, 2nd ed., Springer, New York.
#' 
#' @keywords classif
#' 
#' @examples
#' \dontrun{
#' data(phoneme)
#' mlearn<-phoneme[["learn"]][1:100]
#' glearn<-phoneme[["classlearn"]][1:100]
#' 
#' #	ESTIMATION knn
#' out1=classif.knn(glearn,mlearn,knn=3)
#' summary(out1)
#' 
#' #	PREDICTION knn
#' mtest<-phoneme[["test"]][1:100]
#' gtest<-phoneme[["classtest"]][1:100]
#' pred1=predict(out1,mtest)
#' table(pred1,gtest)
#' 
#' #	ESTIMATION kernel 
#' h=2^(0:5)
#' # using metric distances computed in classif.knn
#' out2=classif.kernel(glearn,mlearn,h=h,metric=out1$mdist)
#' summary(out2)
#' #	PREDICTION kernel
#' pred2=predict(out2,mtest)
#' table(pred2,gtest)
#' }
#' 
#' @export  
predict.classif <- function (object, new.fdataobj = NULL,
                             type = "class", ...) {
  #print("entra predict.classif")
 # new.fdataobj<-dat
  #object<-out1
  if (is.null(object)) 
    stop("No classif object entered")
  if (is.null(new.fdataobj)) 
    return(object$group.est)
  if (is.null(object$prob)) object$prob<-0.5
  isfdata <- is.fdata(new.fdataobj)
  #object$group <- factor(object$group, levels = levels(object$group)[which(table(object$group) >    0)])
  if (is.null(object$levels)) object$levels <-  levels(object$group)
  #output <- pred2internal(object, new.fdataobj = new.fdataobj , type = type, ...)
#print(as.character(object$C[[1]]))
  output <- switch(as.character(object$C[[1]]),
                   classif.gkam2boost={
                     pred2gkam2boost(object, new.fdataobj =new.fdataobj,...)},                 
                   classif.gsam2boost={
                     #pred2gsam2boost(object, new.fdataobj =new.fdataobj,...)
                     #print(object$C[1])
                     pred2gsam2boost(object, new.fdataobj =new.fdataobj,...)
                     },                 
                   classif.glm2boost={
                     pr <- pred2glm2boost(object, new.fdataobj =new.fdataobj,...)
                     if (!is.list(pr)) pr <-list("group.pred"=pr)
                     pr
                      },                 
                   classif.glm={
                     pred2glm(object, new.fdataobj =new.fdataobj,...)},                 
                   classif.gsam={
                     pred2gsam(object, new.fdataobj = new.fdataobj,...)},                 
                   classif.gkam={ 
                     pred2gkam(object, new.fdataobj = new.fdataobj,...)},                 
                   classif.rpart={
                     object$JJ <- object$basis.list
                     pred2ML(object, new.fdataobj = new.fdataobj)},
                   classif.svm={
                     object$JJ <- object$basis.list
                     pred2ML(object, new.fdataobj = new.fdataobj, ...)},
                   classif.nnet={
                     object$JJ <- object$basis.list
                     pred2ML(object, new.fdataobj = new.fdataobj, ...)},
                   classif.gbm={
                     object$JJ <- object$basis.list
                     pred2ML(object, new.fdataobj = new.fdataobj, ...)},
                   classif.multinom={
                     object$JJ <- object$basis.list
                     pred2ML(object, new.fdataobj = new.fdataobj, ...)},
                   classif.randomForest={
                     object$JJ <- object$basis.list
                     pred2ML(object, new.fdataobj = new.fdataobj, ...)},                   
                   classif.naiveBayes={
                     object$JJ <- object$basis.list
                     pred2ML(object, new.fdataobj = new.fdataobj, ...)},                   
                   classif.ksvm={
                     object$JJ <- object$basis.list
                     pred2ML(object, new.fdataobj = new.fdataobj, ...)},   
                   classif.lda={
                     object$JJ <- object$basis.list
                     pred2ML(object, new.fdataobj = new.fdataobj, ...)},   
                   classif.qda={
                     object$JJ <- object$basis.list
                     pred2ML(object, new.fdataobj = new.fdataobj, ...)},   
                   
                   classif.cv.glmnet={# 20201214 
                     object$JJ <- object$basis.list
                     pred2ML(object, new.fdataobj = new.fdataobj, ...)}, 
                   
                   classif.np={
                     pred2np(object, new.fdataobj = new.fdataobj,...)},
                    classif.adaboost={
                      
                     pr<-predict.classif.adaboost(object, newdata = new.fdataobj, type = "probs", ...)
                     pr
                     },
                   classif.bootstrap={
   #                  print("aaa")
                     #pr<-predict.classif.bootstrap(object,
                     pr<-pred2boot(object,
                                               newdata = new.fdataobj
                                               , type = "probs", ...)
                     #print("si predice bootstap")
                     pr
                     },
                   classif.DD={
                     #print("si predice DD")
                     pr<-predict.classif.DD(object,  new.fdataobj,  type, ...)
                     #print(length(pr))
                     if (!is.list(pr)) pr <-list("group.pred"=pr)
                     pr
                     }
  )
  #print(output)
  #print("sale predict.classif")
 # print(names(output))

  #print(type)
        if (type == "class") {    return(output$group.pred)}
  else return(output) 
#  if (type == "probs") {    return(output)  }
#  else stop("Type argument should be one of 'class' or 'probs'")
}



#############################################################################
pred2glm2boost <- function(object, new.fdataobj = NULL,  ...) {
  isfdata <- FALSE
  if (is.fdata(new.fdataobj))  isfdata<-TRUE
  newdata <- object$data 
  if (isfdata)     newfdata <- list(X = new.fdataobj)
  else newfdata <- new.fdataobj

  lev <- levels(object$group)
  prob <- ngroup <- length(lev)
  nn <- nrow(new.fdataobj[[1]])
  prob.group <- array(NA, dim = c(nn,ngroup))
  colnames(prob.group) <- lev
  if (is.null(object$prob)) object$prob<- .5
  if (ngroup == 2) {
    probs <- predict.fregre.glm(object$fit[[1]], newx = newfdata, ...)
    yest <- ifelse(probs > object$prob, lev[2], lev[1])  
    prob.group[, 1] <- 1-probs
    prob.group[, 2] <-  probs
    group.pred <- factor(yest, levels = lev)
    return(list("group.pred"=group.pred,"prob.group"=prob.group))
  }  else {
#print("mas de un grupo1")    
#    print(dim(newfdata[[1]]))
#    print(class(object))
#    group.pred<- predict(object, newx = newfdata,...)
#print("mas de un grupo2")        
#print(group.pred)  
    #return(list("group.pred"=group.pred))
     for (i in 1:ngroup) {
       obj <- object$fit[[i]]

 #      print(obj$formula)
       #prob.group[, i] <- predict.fregre.glm(obj, newx = newfdata, ...)
       pr<-predict.glm(obj,  newfdata$df,type="response")
       #pr <- predict.classif(obj, newx = newfdata, ...)
       prob.group[, i] <- pr
     }
     yest <- apply(prob.group, 1, which.min)
     group.pred <- factor(lev[yest], levels = lev)
  }
  return(list("group.pred"=group.pred,"prob.group"=prob.group))
}

##########################################
pred2glm <- function(object, new.fdataobj = NULL, ...) {
#  print("pred2glm")  
  if (is.null(object$type)) object$type ="1vsall"
  lev <- levels(object$group)
  prob <- ngroup <- length(lev)
  prob.group <- array(NA, dim = c(nrow(new.fdataobj[[1]]), 
                                  ngroup))
  colnames(prob.group) <- lev								  
  if (ngroup == 2) {
    probs <- predict.fregre.glm(object$fit[[1]], 
                                newx = new.fdataobj,...)
    yest <- ifelse(probs > object$prob, lev[2], lev[1])
    prob.group[, 1] <- 1- probs
    prob.group[, 2] <- probs
    group.pred <- factor(yest, levels = lev)
  }  else {
    if (object$type=="majority") {
      if (ngroup > 2) {
        #warning("Majority voting classification")
        cvot <- combn(ngroup, 2)
        nvot <- ncol(cvot)
        nn <- nrow(new.fdataobj[[1]])
        votos <- matrix(0, nn,ngroup)
        #print(dim(new.fdataobj[[1]]))        
        for (ivot in 1:nvot) {
          pred <- predict.fregre.glm(object$fit[[ivot]],new.fdataobj,...)
          group.log <- pred > object$prob
          votos[,cvot[1, ivot]] <- votos[,cvot[1, ivot]] + as.numeric(group.log)
          votos[,cvot[2, ivot]] <- votos[,cvot[2, ivot]] + as.numeric(!group.log)
        }
        maj.voto <- apply(votos, 1, which.max)
        group.pred <- factor(lev[maj.voto], levels = lev)
      }
    }    else{
      for (i in 1:ngroup) {
        prob.group[, i] <- predict.fregre.glm(object$fit[[i]], 
                                              newx = new.fdataobj, ...)
      }
      yest <- apply(prob.group, 1, which.max)
      group.pred <- factor(lev[yest], levels = lev)
    }
  }
  return(list("group.pred"=group.pred,"prob.group"=prob.group))
}
##########################################
pred2lda <- function(object, new.fdataobj = NULL,  ...) {
  # print("entra pred2lda")
  
  lev <- levels(object$group)
  prob <- ngroup <- length(lev)
  nn <- nrow(new.fdataobj)
  prob.group <- array(NA, dim = c(nn, ngroup))
  colnames(prob.group) <- lev								  
  #new.fdataobj <- as.matrix(new.fdataobj)
  if (ngroup == 2 | object$type != "majority") {
    #print("2 grupos o 1vsAll")
    #print(class(object$fit[[1]]))
    probs=predict(object$fit[[1]], new.fdataobj,...)
    prob.group <- probs$posterior
    group.pred <- probs$class
  }  else {
    #warning("Majority voting classification")
    cvot <- combn(ngroup, 2)
    nvot <- ncol(cvot)
    votos <- matrix(0, nn,ngroup)
    for (ivot in 1:nvot) {
      pred<-predict(object$fit[[ivot]],new.fdataobj ,...)
      group.log <- pred$posterior[,2] > object$prob
      votos[,cvot[1, ivot]] <- votos[,cvot[1, ivot]] + as.numeric(group.log)
      votos[,cvot[2, ivot]] <- votos[,cvot[2, ivot]] + as.numeric(!group.log)
    }
    maj.voto <- apply(votos, 1, which.max)
    group.pred <- factor(lev[maj.voto], levels = lev)
    prob.group <- votos/(ngroup-1)#apply(votos, 1, sum)
    colnames(prob.group) <- lev
  }
  #  print("sale pred2lda")
  return(list("group.pred"=group.pred,"prob.group"=prob.group))
}
##########################################
pred2ksvm <- function(object, new.fdataobj = NULL,  ...) {
  #print("pred2ksvm")
  lev <- levels(object$group)
  prob <- ngroup <- length(lev)
  #prob.group <- array(NA, dim = c(nrow(new.fdataobj[[1]]), ngroup))
  nn <- nrow(new.fdataobj)
  prob.group <- array(NA, dim = c(nn, ngroup))
  colnames(prob.group) <- lev								  
  new.fdataobj <- as.matrix(new.fdataobj)
  if (ngroup == 2) {
    #yest <- predict(object$fit[[1]], newx = new.fdataobj,...)
    #yest <- ifelse(probs > object$prob, lev[2], lev[1])
    probs=predict(object$fit[[1]], new.fdataobj ,type="linear.predictor")
    probs<-exp(-probs)/(1+exp(-probs))
    yest <- ifelse(probs < object$prob, lev[2], lev[1])
    prob.group[, 1] <- probs
    prob.group[, 2] <- 1- probs
    group.pred <- factor(yest, levels = lev)
  }  else {
    #if (object$type=="majority") {
    #print("pred2ksvm2")        
    if (ngroup > 2) {
      #warning("Majority voting classification")
      cvot <- combn(ngroup, 2)
      nvot <- ncol(cvot)
      #nn <- nrow(new.fdataobj[[1]])
      votos <- matrix(0, nn,ngroup)
      for (ivot in 1:nvot) {
        pred<-predict(object$fit[[ivot]],new.fdataobj ,...)
        #pred<-exp(-pred)/(1+exp(-pred))
        group.log <- pred > 0#(object$prob
        votos[,cvot[1, ivot]] <- votos[,cvot[1, ivot]] + as.numeric(group.log)
        votos[,cvot[2, ivot]] <- votos[,cvot[2, ivot]] + as.numeric(!group.log)
      }
      maj.voto <- apply(votos, 1, which.max)
      group.pred <- factor(lev[maj.voto], levels = lev)
      prob.group <- votos/apply(votos, 1, sum)
      colnames(prob.group) <- lev
    }
  }
  return(list("group.pred"=group.pred,"prob.group"=prob.group))
}
#############################################################################
pred2gkam2boost <- function(object, new.fdataobj = NULL,  ...) {
  isfdata <- FALSE
  if (is.fdata(new.fdataobj))  isfdata<-TRUE
  newdata <- object$data
  if (isfdata) 
    newfdata <- list(X = new.fdataobj)
  else newfdata <- new.fdataobj
  lev <- levels(object$group)
  prob <- ngroup <- length(lev)
  nn <- nrow(new.fdataobj[[1]])
  prob.group <- array(NA, dim = c(nn,ngroup))
  colnames(prob.group) <- lev
  if (ngroup == 2) {
    probs <- predict.fregre.gkam(object$fit[[1]], newx = newfdata, ...)
    yest <- ifelse(probs > object$prob, lev[2], lev[1])  
    prob.group[, 1] <- 1-probs
    prob.group[, 2] <-  probs
  }
  else {
    for (i in 1:ngroup) {
      obj <- object$fit[[i]]
      prob.group[, i] <- predict.fregre.gkam(obj, newx = newfdata, ...)
    }
    yest <- apply(prob.group, 1, which.max)
  }
  group.pred <- factor(lev[yest], levels = lev)
  return(list("group.pred"=group.pred,"prob.group"=prob.group))
}
#############################################################################
# pred2gsam2boost <- function(object, new.fdataobj = NULL,  ...) {
#   isfdata <- FALSE
#   if (is.fdata(new.fdataobj))  isfdata<-TRUE
#   newdata <- object$data
#   if (isfdata) 
#     newfdata <- list(X = new.fdataobj)
#   else newfdata <- new.fdataobj
#   lev <- levels(object$group)
#   prob <- ngroup <- length(lev)
#   nn <- nrow(new.fdataobj[[1]])
#   prob.group <- array(NA, dim = c(nn,ngroup))
#   colnames(prob.group) <- lev
#   if (is.null(object$prob)) object$prob<- .5
#   if (ngroup == 2) {
#     probs <- predict.fregre.gsam(object$fit[[1]], newx = newfdata, ...)
#     yest <- ifelse(probs > object$prob, lev[2], lev[1])  
#     prob.group[, 1] <- 1-probs
#     prob.group[, 2] <-  probs
#     group.pred <- factor(yest, levels = lev)
#     return(list("group.pred"=group.pred,"prob.group"=prob.group))
#   }  else {
#     #print("mas de un grupo")    
#     group.pred<- predict.classif(object, newx = newfdata,...)
#     return(list("group.pred"=group.pred))
#   }
# }
##########################################

#############################################################################
pred2np <- function(object, new.fdataobj = NULL, ...) {
  #print("entra np ")
  isfdata <- is.fdata(new.fdataobj)
  #if (!isfdata) new.fdataobj<-fdata(new.fdataobj)
  #print(isfdata)  
  gg <- 1:nrow(new.fdataobj)
  if (isfdata) {
    nas <- is.na.fdata(new.fdataobj)
    if (any(nas)) {
      bb <- !nas
      cat("Warning: ", sum(nas), " curves with NA are omited\n")
      new.fdataobj$data <- new.fdataobj$data[bb, ]
      gg <- gg[bb]
    }
    newx <- new.fdataobj[["data"]]
    tt <- new.fdataobj[["argvals"]]
    rtt <- new.fdataobj[["rangeval"]]
  }
  else newx <- as.matrix(new.fdataobj)
  nn <- nrow(new.fdataobj)
  x = object$fdataobj
  y = object$y
  h = object$h.opt
  n = nrow(x)
  nn = nrow(newx)
  np <- ncol(x)
  lev <- levels(y)
  numg <- nlevels(y)
  if (is.null(rownames(newx))) 
    rownames(newx) <- 1:nn
  bs = as = list()
  C <- object$call
  m <- object$m
  Ker = object$Ker
  par.metric <- list()
  par.metric <- attr(object$mdist, "par.metric")
  parm <- attr(object$mdist, "par.metric")
  a1 <- attr(object$mdist, "call")
  if (isfdata) {
    if (a1 == "semimetric.mplsr") {
      par.metric[["fdata1"]] <- x
      par.metric[["fdata2"]] <- new.fdataobj
      nmdist <- t(do.call(a1, par.metric))
    }    else {
      par.metric[["fdata2"]] <- x
      par.metric[["fdata1"]] <- new.fdataobj
      nmdist <- do.call(a1, par.metric)
    }
  }  else {
    #print("no fdata")     
    par.metric[["x"]] <- new.fdataobj
    par.metric[["y"]] <- x
    nmdist <- do.call(a1, par.metric)   }
  object$par.S$tt <- nmdist  
  kmdist =  object$type.S(nmdist, h = h, Ker = object$Ker,
                          w=object$par.S$w,cv =FALSE) #SIEMPRE 
  #  print("h y kmdist");  print(kmdist[24:25,])
  kmdist[is.na( kmdist)]<-1e-28
  pgrup = matrix(0, numg, nn)
  rownames(pgrup) <-lev
  l = array(0, dim = c(nn))
  group.pred = array(0, dim = nn)
  for (j in 1:numg) {
    grup = as.integer(y == lev[j])
    pgrup[j, ] <- kmdist %*% matrix(grup, ncol = 1)
  }
  #print(kmdist[25,])
  #print(pgrup[25,])
  group.pred<-unlist(apply(pgrup, 2, which.max))
  group.pred <-   factor(lev[unlist(apply(pgrup, 2, which.max))],levels=lev)
  pgrup <- t(pgrup)
  #group.est <- numeric(nn)
  group.est<-factor(character(nn),levels=lev)
  ty <-  object$ty
  if (ty == "S.KNN") {  
    for (ii in 1:nn) {
      l = seq_along(pgrup[ii, ])[pgrup[ii, ] == max(pgrup[ii,], na.rm = TRUE)]
      if (length(l) > 1) {
        abc <- which(nmdist[ii, ] == min(nmdist[ii,l], na.rm = TRUE))
        group.est[ii] <- y[abc[1]]
      }
      else  group.est[ii] = lev[l[1]]
    }
    group.pred <- factor(lev[group.est], levels = lev)
  }
  out<-list("group.pred"=group.pred,"prob.group"=pgrup)
  out
}
#############################################################################
# pred2grm <- function(object, new.fdataobj = NULL, ...) {
#   lev <- levels(object$group)
#   prob <- ngroup <- length(lev)
#   prob.group <- array(NA, dim = c(nrow(new.fdataobj[[1]]), ngroup))
#   colnames(prob.group) <- lev
#   if (ngroup == 2) {
#     probs <- predict.fregre.grm(object$fit[[1]], newx = new.fdataobj, ...)
#     yest <- ifelse(probs > object$prob, lev[2], lev[1])
#     prob.group[, 1] <- 1- probs
#     prob.group[, 2] <- probs
#     group.pred <- factor(yest, levels= lev)
#   }  else {
#     if (object$type=="majority") {
#       if (ngroup > 2) {
#         # warning("Majority voting classification")
#         cvot <- combn(ngroup, 2)
#         nvot <- ncol(cvot)
#         nn <- nrow(new.fdataobj[[1]])
#         votos <- matrix(0, nn,ngroup)
#         for (ivot in 1:nvot) {
#           pred<-predict.fregre.grm(object$fit[[ivot]],new.fdataobj,...)
#           group.log <- pred > object$prob
#           votos[,cvot[1, ivot]] <- votos[,cvot[1, ivot]] + as.numeric(group.log)
#           votos[,cvot[2, ivot]] <- votos[,cvot[2, ivot]] + as.numeric(!group.log)
#         }
#         maj.voto <- apply(votos, 1, which.max)
#         group.pred <- factor(lev[maj.voto], levels = lev)
#       }
#     }    else{
#       #    print("1 vs all")
#       for (i in 1:ngroup) {
#         prob.group[, i] <- predict.fregre.grm(object$fit[[i]], newx = new.fdataobj, ...)
#       }
#       group.pred <- factor(lev[apply(prob.group, 1, which.max)], levels = lev)
#     }  }
#   return(list("group.pred"=group.pred,"prob.group"=prob.group))
# }
#############################################################################
pred2gkam <- function(object, new.fdataobj = NULL,  ...) {
  lev <- levels(object$group)
  prob <- ngroup <- length(lev)
  prob.group <- array(NA, dim = c(nrow(new.fdataobj[[1]]), ngroup))
  colnames(prob.group) <- lev
  if (ngroup == 2) {
    #probs <- predict.fregre.gkam(object$fit[[1]], newx = new.fdataobj, ...)
    probs <- predict(object$fit[[1]], newx = new.fdataobj, ...)
    yest <- ifelse(probs > object$prob, lev[2], lev[1])
    prob.group[, 1] <- 1- probs
    prob.group[, 2] <- probs
    group.pred <- factor(yest, levels= lev)
  }  else {
    if (object$type=="majority") {
      if (ngroup > 2) {
        # warning("Majority voting classification")
        # print("majority")        
        cvot <- combn(ngroup, 2)
        nvot <- ncol(cvot)
        nn <- nrow(new.fdataobj[[1]])
        votos <- matrix(0, nn,ngroup)
        for (ivot in 1:nvot) {
          #pred<-predict.fregre.gkam(object$fit[[ivot]],new.fdataobj,...)
          pred<-predict(object$fit[[ivot]],new.fdataobj,...)
          group.log <- pred > object$prob
          votos[,cvot[1, ivot]] <- votos[,cvot[1, ivot]] + as.numeric(group.log)
          votos[,cvot[2, ivot]] <- votos[,cvot[2, ivot]] + as.numeric(!group.log)
        }
        maj.voto <- apply(votos, 1, which.max)
        group.pred <- factor(lev[maj.voto], levels = lev)
      }
    }    else{
      # print("1 vs all")
      for (i in 1:ngroup) {
        #prob.group[, i] <- predict.fregre.gkam(object$fit[[i]], newx = new.fdataobj, ...)
        prob.group[, i] <- predict(object$fit[[i]], newx = new.fdataobj, ...)
      }
      yest <- apply(prob.group, 1, which.max)
      group.pred <- factor(lev[yest], levels = lev)
    }
  }
  return(list("group.pred"=group.pred,"prob.group"=prob.group))
}

####################################
pred2gsam <- function(object, new.fdataobj = NULL, ...) {
  lev <- object$levels
  if (is.null(object$type)) object$type ="1vsall"
  prob <- ngroup <- length(lev)
  prob.group <- array(NA, dim = c(nrow(new.fdataobj[[1]]), 
                                  ngroup))
  colnames(prob.group) <- lev
  if (ngroup == 2) {
    probs <- predict.fregre.gsam(object$fit[[1]], newx = new.fdataobj, 
                                 ...)
    yest <- ifelse(probs > object$prob, lev[2], lev[1])
    prob.group[, 1] <- 1- probs
    prob.group[, 2] <- probs
    group.pred <- (factor(yest, levels = lev))
  }
  else {
    if (object$type=="majority") {    
      if (ngroup > 2) {
        #warning("Majority voting classification")
        cvot <- combn(ngroup, 2)
        nvot <- ncol(cvot)
        nn <- nrow(new.fdataobj[[1]])
        pvotos<-votos <- matrix(0, nn,ngroup)
        for (ivot in 1:nvot) {
          pred<-predict.fregre.gsam(object$fit[[ivot]],new.fdataobj,...)
          group.log <- pred > object$prob
          pvotos[,cvot[1, ivot]] <- pvotos[,cvot[1, ivot]] +pred
          pvotos[,cvot[2, ivot]] <- pvotos[,cvot[2, ivot]] + 1-pred
          # print(group.log)
          votos[,cvot[1, ivot]] <- votos[,cvot[1, ivot]] + as.numeric(group.log)
          votos[,cvot[2, ivot]] <- votos[,cvot[2, ivot]] + as.numeric(!group.log)
        }
       # cat(cvot[, ivot])
        maj.voto <- apply(votos, 1, which.max)
        group.pred <- factor(lev[maj.voto], levels = lev)
        prob.grup<-pvotos
      }
    }    else{
      for (i in 1:ngroup) {
        prob.group[, i] <- predict.fregre.gsam(object$fit[[i]], 
                                               newx = new.fdataobj, ...)
      }
      group.pred <- factor(lev[apply(prob.group, 1, which.max)], levels = lev)
    }
  }
  return(list("group.pred"=group.pred,"prob.group"=prob.group))
}


####################################
pred2gsam2boost<- function(object, new.fdataobj = NULL, ...) {
  lev <- object$levels
  if (is.null(object$type)) object$type ="1vsall"
  prob <- ngroup <- length(lev)
  prob.group <- array(NA, dim = c(nrow(new.fdataobj[[1]]), 
                                  ngroup))
  colnames(prob.group) <- lev
  if (ngroup == 2) {
    probs <- predict.fregre.gsam(object$fit[[1]], newx = new.fdataobj, 
                                 ...)
    yest <- ifelse(probs > object$prob, lev[2], lev[1])
    prob.group[, 1] <- 1- probs
    prob.group[, 2] <- probs
    group.pred <- (factor(yest, levels = lev))
  }
  else {
    if (object$type=="majority") {    
      if (ngroup > 2) {
        #warning("Majority voting classification")
        cvot <- combn(ngroup, 2)
        nvot <- ncol(cvot)
        nn <- nrow(new.fdataobj[[1]])
        pvotos<-votos <- matrix(0, nn,ngroup)
        for (ivot in 1:nvot) {
          pred<-predict.fregre.gsam(object$fit[[ivot]],new.fdataobj,...)
          group.log <- pred > object$prob
          pvotos[,cvot[1, ivot]] <- pvotos[,cvot[1, ivot]] +pred
          pvotos[,cvot[2, ivot]] <- pvotos[,cvot[2, ivot]] + 1-pred
          # print(group.log)
          votos[,cvot[1, ivot]] <- votos[,cvot[1, ivot]] + as.numeric(group.log)
          votos[,cvot[2, ivot]] <- votos[,cvot[2, ivot]] + as.numeric(!group.log)
        }
       # cat(cvot[, ivot])
        maj.voto <- apply(votos, 1, which.max)
        group.pred <- factor(lev[maj.voto], levels = lev)
        prob.grup<-pvotos
      }
    }    else{
      for (i in 1:ngroup) {
        prob.group[, i] <- predict.fregre.gsam(object$fit[[i]], 
                                               newx = new.fdataobj, ...)
      }
      group.pred <- factor(lev[apply(prob.group, 1, which.min)], levels = lev)
    }
  }
  return(list("group.pred"=group.pred,"prob.group"=prob.group))
}
#############################

#pred2pc <- function(object,XX,lev){
pred2pc <- function(object,XX){
# print("pred2pc")
  lev <- levels(object$group)
  # if (object$C[[1]]=="classif.svm"){
  #   yest = predict(object = object$fit, newdata = XX,    probability = TRUE)  
  #   return(list(group.pred =yest,prob.group=attr(yest,"probabilities")))
  # }
  if (object$C[[1]]=="classif.nnet" | object$C[[1]]=="classif.multinom"){
    yest = predict(object = object$fit, newdata = XX,    type = "class")
    prob.group = predict(object = object$fit, newdata = XX)  
    if (length(lev)==2) {
      prob.group <-cbind(1-prob.group ,prob.group )  
      colnames(prob.group)<-lev
    }
    group.pred <- factor(yest, levels = lev)
    return(list(group.pred =group.pred,prob.group=prob.group ))
  }
  #if (object$C[[1]]=="classif.rpart"){
  #  yest = predict(object = object$fit, newdata = XX,    type = "prob")
  #}
  if (object$C[[1]]=="classif.rpart"){
    group.pred = predict(object = object$fit, newdata = XX,    type = "class")
    prob.group = predict(object = object$fit, newdata = XX,    type = "prob") 
    return(list(group.pred =group.pred,prob.group=prob.group ))
  }
  if (object$C[[1]]=="classif.cv.glmnet"){
    #print(names(XX))
    group.pred = predict(object = object$fit, newx = as.matrix(XX),    type = "class")
    prob.group = predict(object = object$fit, newx = as.matrix(XX), type="response" ) 
    group.pred <- factor(group.pred, levels = lev)
    return(list(group.pred =group.pred,prob.group=prob.group ))
  }
  yest = predict(object = object$fit, XX)
  if (!is.factor(yest)) group.pred <- factor(yest, levels = lev)
  else group.pred <- yest
  return(list(group.pred =group.pred))
}

##########################################
pred2svm <- function(object, XX = NULL, ...) {
  #print("pred2svm")  
  lev <- levels(object$group)
  prob <- ngroup <- length(lev)
  ngroup <- length(lev)
  nn <- NROW(XX)
  prob.group <- matrix(NA, nn,  ngroup)
  colnames(prob.group) <- lev								  
  if (ngroup == 2 | object$type=="majority") {
    #print("entra maj")
    yest = predict(object = object$fit, newdata = XX,    probability = TRUE)  
    prob.group=attr(yest,"probabilities")
    attributes(yest) = attributes(yest)[1:3]
    return(list(group.pred =yest,prob.group=prob.group))
  }  else {
    #print( "nlev >2 and type='1vsAll'")
    for (i in 1:ngroup) {
      aux = predict(object = object$fit[[i]], newdata = XX,    probability = TRUE)  
      aux = attributes(aux)$probabilities
      ii<-colnames(aux)==0
      prob.group[,i]<-aux[,ii]
    }
    yest <- apply(prob.group, 1, which.max)
    group.pred <- factor(lev[yest], levels = lev)
  }
  return(list("group.pred"=group.pred,"prob.group"=prob.group))
}
##########################################
pred2rpart <- function(object, XX = NULL, ...) {
  #print("pred2rpart")  
  lev <- levels(object$group)
  prob <- ngroup <- length(lev)
  ngroup <- length(lev)
  nn <- NROW(XX)
  prob.group <- matrix(NA, nn,  ngroup)
  colnames(prob.group) <- lev								  
  if (ngroup == 2 | object$type=="majority") {
 #   print("entra maj")
    yest = predict(object = object$fit, newdata = XX, type="class")  
    prob.group= predict(object = object$fit, newdata = XX, type="prob")  
    return(list(group.pred =yest,prob.group=prob.group))
  }  else {
    #print( "nlev >2 and type='1vsAll'")
    for (i in 1:ngroup) {
      aux = predict(object = object$fit[[i]], newdata = XX,   type="prob")  
      prob.group[,i]<-aux[,1]
    }
    yest <- apply(prob.group, 1, which.max)
    group.pred <- factor(lev[yest], levels = lev)
  }
  return(list("group.pred"=group.pred,"prob.group"=prob.group))
}
##########################################
pred2randomForest <- function(object, XX = NULL, ...) {
  #print("pred2RF")  
  lev <- levels(object$group)
  prob <- ngroup <- length(lev)
  ngroup <- length(lev)
  nn <- NROW(XX)
  prob.group <- matrix(NA, nn,  ngroup)
  colnames(prob.group) <- lev								  
  if (ngroup == 2 | object$type=="majority") {
    #print("entra maj")
    yest = predict(object = object$fit, newdata = XX)
    prob.group= predict(object = object$fit, newdata = XX,    type="prob")
    #votes = predict(object = object$fit, newdata = XX,    type="vote")  
    #return(list(group.pred =yest,prob.group=prob.group,votes=votes))
    return(list(group.pred =yest,prob.group=prob.group))
  }  else {
    #print( "nlev >2 and type='1vsAll'")
    for (i in 1:ngroup) {
      aux = predict(object = object$fit[[i]], newdata = XX,    type="prob")
      prob.group[,i]<- aux[,1]
    }
    yest <- apply(prob.group, 1, which.max)
    group.pred <- factor(lev[yest], levels = lev)
  }
  return(list("group.pred"=group.pred,"prob.group"=prob.group))
}


##########################################
pred2gbm <- function(object, XX = NULL, ...) {
#  print("pred2gbm")  
  lev <- levels(object$group)
  prob <- ngroup <- length(lev)
  ngroup <- length(lev)
  nn <- NROW(XX)
  prob.group <- matrix(NA, nn,  ngroup)
  colnames(prob.group) <- lev								  
  if (ngroup == 2 | object$type=="majority") {
   #    print("entra maj")
    prob.group = predict(object = object$fit, newdata = XX, type="response")  
    yest <- apply(prob.group,1,which.max)
    yest <- factor(yest ,levels=lev)                # level o label???
    return(list(group.pred =yest,prob.group=prob.group))
  }  else {
    #print( "nlev >2 and type='1vsAll'")
    for (i in 1:ngroup) {
      aux =  suppressWarnings(predict(object = object$fit[[i]], newdata = XX, type="response"))
      prob.group[,i]<-aux[,1]
    }
    yest <- apply(prob.group, 1, which.max)
    group.pred <- factor(lev[yest], levels = lev) # level o label???
  }
  return(list("group.pred"=group.pred,"prob.group"=prob.group))
}

##########################################
pred2ML <- function(object, new.fdataobj = NULL, ...) {
#print("pred2ML")  
  newx = data = new.fdataobj
  basis.x = object$basis.x
  formula = object$formula.ini
  tf <- terms.formula(formula)
  terms <- attr(tf, "term.labels")
  nt <- length(terms)
  lev <-levels(object$group.est)
  if (attr(tf, "response") > 0) {
    response <- as.character(attr(tf, "variables")[2])
    pf <- rf <- paste(response, "~", sep = "")
  }      else pf <- rf <- "~"
  if (attr(tf, "intercept") == 0) {
 #   print("No intecept")
    pf <- paste(pf, -1, sep = "")                   }
  vtab <- rownames(attr(tf, "factors"))
  vnf = intersect(terms, names(data$df))
  vnf2 = intersect(vtab[-1], names(data$df))
#print(vnf)  ;  print(vnf2)   


  vfunc2 = setdiff(terms, vnf)
  vint = setdiff(terms, vtab)
  vfunc = setdiff(vfunc2, vint)
#print(vfunc)
  vnf = c(vnf2, vint)
  off <- attr(tf, "offset")
  kterms = 1
  if (length(vnf) > 0) {
  #  print("entra vnf")    
    #print(vnf2)
    first = FALSE
    XX = data.frame(data$df[, c(vnf2)])
    names(XX) = vnf2
    for (i in 1:length(vnf)) {
      if (kterms > 1) 
        pf <- paste(pf, "+", vnf[i], sep = "")
      else pf <- paste(pf, vnf[i], sep = "")
      kterms <- kterms + 1
    }
    if (attr(tf, "intercept") == 0) {
      print("No intecept")
      pf <- paste(pf, -1, sep = "")
    }
  }   else first = TRUE
  if (length(vfunc) > 0) {
    for (i in 1:length(vfunc)) {
      #print("vfunc");      print(vfunc[i])
      if (inherits(data[[vfunc[i]]], "fdata")) {
        fdataobj <- data[[vfunc[i]]]
        dat <- fdataobj$data
        tt <- fdataobj[["argvals"]]
        if (is.null(rownames(dat))) 
          rownames(dat) <- 1:nrow(dat)
        fdnames = list(time = tt, reps = rownames(dat), 
                       values = "values")
        x.fd <- fdataobj[["data"]]
        tt <- fdataobj[["argvals"]]
        rtt <- fdataobj[["rangeval"]]
        if (object$basis.x[[vfunc[i]]]$type != "pc" & 
            object$basis.x[[vfunc[i]]]$type != "pls") {
 # print("basiss")
          x.fd = Data2fd(argvals = tt, y = t(fdata.cen(fdataobj, 
                                                       object$mean[[vfunc[i]]])[[1]]$data), 
                         basisobj = object$basis.x[[vfunc[i]]], 
                         fdnames = fdnames)
          r = x.fd[[2]][[3]]
         # J <- object$JJ[[vfunc[i]]]
          J <- object$basis.list[[vfunc[i]]]
          Z = t(x.fd$coefs) %*% J
          colnames(Z) = colnames(J)
        } else {
        #  print("eeeeeeeeeeeeeeeeeeellllllllllllllllllllllllllssssssssssssssssssssssssssssssssseeeeeeeeeeee")
          # return(pred2pc(object,vfunc,fdataobj,lev=lev))
          name.coef <- paste(vfunc[i], ".", 
                             rownames(object$basis.x[[vfunc[i]]]$basis$data),
                             sep = "")
          newXcen <- fdata.cen(fdataobj, object$mean[[vfunc[i]]])[[1]]
          if (object$basis.x[[vfunc[i]]]$type =="pls") {
            if (object$basis.x[[vfunc[i]]]$norm) {
              sd.X <- sqrt(apply(fdataobj$data, 2,var))
              newXcen$data <- newXcen$data/(rep(1,nrow(newXcen)) %*% t(sd.X))
            }
          }
          Z <- inprod.fdata(newXcen, object$basis.list[[vfunc[i]]])
          colnames(Z) <- name.coef
          Z <- data.frame(Z)
          
        }
    #print(first)
        if (first) {
          XX = Z
          first = FALSE
        }   else XX = cbind(XX, Z)
      }    else {# fd clas
        if (inherits(data[[vfunc[i]]], "fd")) {
          if (!inherits(object$basis.x[[vfunc[i]]],"pca.fd")) {
            x.fd <- fdataobj <- data[[vfunc[i]]]
            r = x.fd[[2]][[3]]
            J <- object$JJ[[vfunc[i]]]
            x.fd$coefs <- x.fd$coefs - object$mean[[vfunc[i]]]$coefs[, 
                                                                     1]
            Z = t(x.fd$coefs) %*% J
            colnames(Z) = colnames(J)
          }  else {
            name.coef[[vfunc[i]]] = paste(vfunc[i], 
                                          ".", colnames(object$basis.x[[vfunc[i]]]$harmonics$coefs), 
                                          sep = "")
            data[[vfunc[i]]]$coefs <- sweep(data[[vfunc[i]]]$coefs, 
                                            1, (object$basis.x[[vfunc[i]]]$meanfd$coefs), 
                                            FUN = "-")
            fd.cen <- data[[vfunc[i]]]
            Z <- inprod(fd.cen, object$basis.x[[vfunc[i]]]$harmonics)
            colnames(Z) <- name.coef[[vfunc[i]]]
          }
          if (first) {
            XX = Z
            first = FALSE
          } else XX = cbind(XX, Z)
        }   #else stop("Please, enter functional covariate")
      }
    }
  }
  if (!is.data.frame(XX)) 
    XX = data.frame(XX)
if (object$C[[1]]=="classif.cv.glmnet" & NCOL(XX)==1)
  XX <- cbind(rep(1,len=NROW(XX)),XX)
  if (object$C[[1]]=="classif.ksvm")   out = pred2ksvm(object = object, XX)
  else if (object$C[[1]]=="classif.svm")    out = pred2svm(object = object, XX)
  else if (object$C[[1]]=="classif.randomForest")    out = pred2randomForest(object = object, XX)
  else if (object$C[[1]]=="classif.rpart")    out = pred2rpart(object = object, XX)
  else if (object$C[[1]]=="classif.lda" | object$C[[1]]=="classif.qda")    out = pred2lda(object = object, XX)
  else if (object$C[[1]]=="classif.gbm" )    out = pred2gbm(object = object, XX)  
  else  out = pred2pc(object , XX)
  #print("sale ML")#;print(out)  
  return(out)
}



# 
# class(ldatasim)<-c("ldata",class(ldatasim))
# pred2<-predict.classif( res.glm2,ldatatest)
# table(newy,pred2)
# 
# ii <- 1:100
# pred2<-predict.classif( res.glm2,ldatasim[ii,row=T])
# 
# table(pred2,res.glm2$group.est)
# 
# table(ldatasim$df$y,res.glm2$group.est)
# table(ldatasim$df$y[ii],pred2)
# 
# 
# traceback()
# ##################################################
#pred3 <- predict.classif( res.glm3,ldatatest)

