#' @title Performance measures for regression and classification models
#' @name accuracy
#' @family performance
#' 
#' @description \code{\link{cat2meas}} and \code{\link{tab2meas}} calculate the measures for a multiclass classification model.\cr
#'  \code{\link{pred2meas}} calculates the measures for a regression model.
#'  
#' @details 
#' \itemize{
#' \item \code{\link{cat2meas}} compute \eqn{tab=table(yobs,ypred)} and calls \code{\link{tab2meas}} function.
#' \item \code{\link{tab2meas}} function computes the following measures (see \code{measure} argument) for a binary classification model:
#' \itemize{
#' \item \code{accuracy}{ the accuracy classification score}
#' \item \code{recall}, \code{sensitivity,TPrate}{ \eqn{R=TP/(TP+FN)}}
#' \item \code{precision}{ \eqn{P=TP/(TP+FP)}}
#' \item \code{specificity},\code{TNrate}{ \eqn{TN/(TN+FP)}}
#' \item \code{FPrate}{ \eqn{FP/(TN+FP)}}
#' \item \code{FNrate}{ \eqn{FN/(TP+FN)}}
#' \item \code{Fmeasure}{ \eqn{2/(1/R+1/P)}}
#' \item \code{Gmean}{ \eqn{sqrt(R*TN/(TN+FP))}}
#' \item \code{kappa}{ the kappa index}
#' \item \code{cost}{ \eqn{sum(diag(tab)/rowSums(tab)*cost)/sum(cost)}}
#' }
#\item \code{\link{tab2meas}} function computes the \code{accuracy}, \code{kappa} and \code{cost} measures  for a multiclass vectors-
#' \item \code{\link{pred2meas}} function computes the following  measures of error, usign the \code{measure} argument, for observed and predicted vectors:
#' \itemize{
#' \item \code{MSE}{ Mean squared error, \eqn{\frac{\sum{(ypred- yobs)^2}}{n} }{\sum (ypred- yobs)^2 /n }}
#' \item \code{RMSE}{ Root mean squared error \eqn{\sqrt{\frac{\sum{(ypred- yobs)^2}}{n} }}{\sqrt(\sum (ypred- yobs)^2 /n )}}
#' \item \code{MAE}{ Mean Absolute Error, \eqn{\frac{\sum |yobs - ypred|}{n}}{\sum |yobs - ypred| /n}}
#' }
#' }
#' 
#' @param yobs  A vector of the labels, true class or observed response. Can be \code{numeric}, \code{character}, or \code{factor}.
#' @param ypred A vector of the predicted labels, predicted class or predicted response. Can be \code{numeric, character, or factor}.
#' @param tab Confusion matrix (Contingency table: observed class by rows, predicted class by columns).
#' @param measure Type of measure, see \code{details} section. 
#' @param cost Cost value by class (only for input factors).
#' @aliases cat2meas tab2meas pred2meas.
#' 
#' @rdname accuracy
#' @export
cat2meas <- function(yobs,ypred,measure="accuracy",cost=rep(1,nlevels(yobs))){
  tab=table(yobs,ypred)
  res=tab2meas(tab,measure=measure,cost=cost)
  return(res)
}

#' @rdname accuracy
#' @export 
tab2meas <- function(tab, measure="accuracy", cost=rep(1,nrow(tab))){
  if (nrow(tab)!=ncol(tab)) stop("nrow(tab)!=ncol(tab)")
  nlev = nrow(tab)
  if (nlev==2) {
    TP = tab[2,2]
    FN = tab[2,1] 
    FP = tab[1,2]
    TN = tab[1,1]
    R = TP/(TP+FN)
    P = TP/(TP+FP) 
    nmeas = length(measure)
    meas = numeric(nmeas)
    for (i in 1:nmeas){
      meas[i]=switch(measure[i],
                     recall=R,
                     sensitivity=R,
                     TPrate=R,
                     specificity=TN/(TN+FP), 
                     TNrate=TN/(TN+FP),
                     FPrate=FP/(TN+FP),
                     FNrate=FN/(TP+FN),
                     precision=P,
                     NPvalue=TN/(TN+FN),
                     Fmeasure=2/(1/R+1/P),
                     F1=2*P*R/(R+P),
                     F2=5*P*R/(4*R+P),
                     Gmean=sqrt(R*TN/(TN+FP)),
                     accuracy=tab2accuracy(tab),
                     kappa=tab2kappa(tab),
                     cost=sum(diag(tab)/rowSums(tab)*cost)/sum(cost)
#                    cost=1-sum(c(FN,FP)*cost)/sum(table(yobs)*cost) 
      )
    }
    names(meas) = measure
    return(meas)
  } else {
    TP = diag(tab)
    Tobs  <- rowSums(tab)
    Tpred <- colSums(tab)
    R <- TP/Tobs
    P <- TP/Tpred
    #FN = tab[2,1] 
    #FP = tab[1,2]
    #TN = tab[1,1]
    #R = TP/(TP+FN)
    #P = TP/(TP+FP) 
    nmeas = length(measure)
    meas = numeric(nmeas)
    if (nmeas>1) warning("For multiclass problems only the first measure is returned")
    #for (i in 1:nmeas){
      meas=switch(measure[1],
                     recall= R,
                     sensitivity=R,
                     TPrate=R,
                     #specificity=TN/(TN+FP), 
                     #TNrate=TN/(TN+FP),
                     #FPrate=FP/(TN+FP),
                     #FNrate=FN/(TP+FN),
                     precision= P,
                     #NPvalue=TN/(TN+FN),
                     Fmeasure=2/(1/R+1/P),
                     F1=2*P*R/(R+P),
                     F2=5*P*R/(4*R+P),
                     #Gmean=sqrt(R*TN/(TN+FP)),
                     cost=sum(diag(tab)/rowSums(tab)*cost)/sum(cost),
                     accuracy=tab2accuracy(tab),
                     #waccuracy=tab2waccuracy(tab),
                     kappa=tab2kappa(tab)
      )
    }
    #names(meas) = measure
    return(meas)
  }


cat2alpha <-function(yobs, ypred, weights, coeflearn="Freund"){
  ind <- as.numeric(yobs != ypred) 
  n<-length(yobs)
  if (missing(weights)) weights <- rep(1,len=n)
  if (sum(weights)!=1)  weights <- weights/sum(weights)
  err <- mean(weights*ind)        
  alpha <- log((1-err)/err)
  if (coeflearn=="Breiman"){	alpha <- (1/2) * alpha	}
  if (coeflearn=="Zhu")    {	alpha <- alpha + log( nlevels(yobs) - 1)	}
  if (alpha<0) alpha=0
  if(alpha==Inf) alpha=10
  return(list("error"=err,"alpha"=alpha))
}

# @export cat2accuracy
# @format none
cat2accuracy=  function(yobs,ypred){
  mean(ypred == yobs)
}

# @export cat2waccuracy
# @format none
# cat2waccuracy = function(yobs, ypred) {
#   lvls <- levels(yobs)
#   accs <- lapply(lvls, function(x) {
#     idx <- which(yobs == x)
#     return(meas2accuracy( yobs[idx], ypred[idx]))
#   })
#   acc <- mean(unlist(accs))
#   return(acc)
# }

# @export cat2wkappa
# @format none
cat2wkappa = function(yobs,ypred){
  cat2meas(yobs=yobs,ypred=ypred,measure="wkappa")
}

# @export cat2kappa
# @format none
cat2kappa = function(yobs,ypred){
  cat2meas(yobs=yobs,ypred=ypred,measure="kappa")
}  
  
############################################################
# @export tab2accuracy
# @format none
tab2accuracy =  function(tab){
    sum(diag(tab))/sum(tab)
}
  
# @export tab2SAE
tab2SAE = function(tab){
  n = sum(tab) 
  #nc = nrow(tab)
  #diag = diag(tab) 
  #rowsums = apply(tab, 1, sum) 
  #colsums = apply(tab, 2, sum) 
  p = apply(tab, 1, sum)  / n 
  q = apply(tab, 2, sum)  / n 
  expAccuracy = sum(p*q)
  accuracy = sum(diag(tab)) / n 
  res<-(accuracy - expAccuracy) / (1 - expAccuracy)
  res
  }

# @export tab2kappa
tab2kappa = function (tab) {
  #function (truth, response) {
  #conf.mat = table(truth, response)
  conf.mat = tab/sum(tab)
  p0 = sum(diag(conf.mat))
  rowsum = rowSums(conf.mat)
  colsum = colSums(conf.mat)
  pe = sum(rowsum * colsum)/sum(conf.mat)^2
  1 - (1 - p0)/(1 - pe)
}

# @export tab2wkappa
tab2wkappa = function (tab) { 
  lev<- rownames(tab) 
  conf.mat = tab/sum(tab)
  rowsum = rowSums(conf.mat)
  colsum = colSums(conf.mat)
  expected.mat = rowsum %*% t(colsum)
  class.values = seq_along(lev) - 1L
  weights = outer(class.values, class.values, 
                  FUN = function(x, y) (x - y)^2)
  1 - sum(weights * conf.mat)/sum(weights * expected.mat)
}


#' @rdname accuracy
#' @export pred.MSE
pred.MSE = function (yobs, ypred) {
  mean((ypred- yobs)^2)
}


# @format none
#' @rdname accuracy
#' @export pred.RMSE
pred.RMSE = function (yobs, ypred) {
  sqrt(pred.MSE(yobs, ypred))
}

# @format none
#' @rdname accuracy
#' @export pred.MAE
pred.MAE = function (yobs, ypred) {
  mean(abs(yobs - ypred))
}

#' @rdname accuracy
#' @export 
pred2meas = function(yobs, ypred, measure="RMSE"){
  nmeas = length(measure)
  meas = numeric(nmeas)
  for (i in 1:nmeas){
      meas[i]=switch(measure[i],
                     RMSE=pred.RMSE(yobs, ypred),
                     MAE=pred.MAE(yobs, ypred),
                     MSE = pred.MSE(yobs, ypred) 
      )
    }
    names(meas) = measure
  return(meas)
}

