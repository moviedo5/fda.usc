#' @name classif.ML
#' 
#' @title Functional classification using ML algotithms
#' 
#' @description Computes functional classification using functional (and non functional)
#' explanatory variables by rpart, nnet, svm or random forest model
#' @note Wrapper versions for multivariate and functional classification:
#' \itemize{
#' \item \code{classif.lda},\code{classif.qda}: uses \code{lda} and  \code{qda} functions and requires \code{MASS} package.
#' \item \code{classif.nnet}: uses \code{nnet} function and requires \code{nnet} package.
#' \item \code{classif.rpart}: uses \code{nnet} function and requires \code{rpart} package.
#' \item \code{classif.svm}, \code{classif.naiveBayes}: uses \code{svm} and  \code{naiveBayes} functions and requires \code{e1071} package.
#' \item \code{classif.ksvm}: uses \code{weighted.ksvm } function and requires \code{personalized} package.
#' \item \code{classif.randomForest}: uses \code{randomForest} function and requires \code{randomForest} package.
#' \item \code{classif.cv.glmnet}: uses \code{cv.glmnet} function and requires \code{glmnet} package.
#' \item \code{classif.gbm}: uses \code{gbm} function and requires \code{gbm} package.
#' }
#' 
#' @details 
#' The first item in the \code{data} list is called \emph{"df"} and is a data
#' frame with the response and non functional explanatory variables, as
#' \code{\link{glm}}.\cr
#' 
#' Functional covariates of class \code{fdata} or \code{fd} are introduced in
#' the following items in the \code{data} list.\cr \code{basis.x} is a list of
#' basis for represent each functional covariate. The b object can be
#' created by the function: \code{\link{create.pc.basis}}, \link[fda]{pca.fd}
#' \code{\link{create.pc.basis}}, \code{\link{create.fdata.basis}} o
#' \link[fda]{create.basis}.\cr \code{basis.b} is a list of basis for
#' represent each functional beta parameter. If \code{basis.x} is a list of
#' functional principal components basis (see \code{\link{create.pc.basis}} or
#' \link[fda]{pca.fd}) the argument \code{basis.b} is ignored.
#' 
#' @aliases classif.rpart classif.nnet classif.randomForest classif.cv.glmnet
#' classif.svm classif.ksvm classif.naiveBayes classif.lda classif.qda classif.multinom
#' @param formula an object of class \code{formula} (or one that can be coerced
#' to that class): a symbolic description of the model to be fitted. The
#' details of model specification are given under \code{Details}.
#' @param data List that containing the variables in the model.
#' @param basis.x List of basis for functional explanatory data estimation.
#' @param weights Weights:   
#' \itemize{
#' \item if \code{character} string \code{='equal'} same weights for each observation (by default) and
#' \code{='inverse'} for inverse-probability of weighting.   
#' \item if \code{numeric} vector of length \code{n}, Weight values of each observation.
#' }
#' @param type If type is\code{"1vsall"}  (by default) 
#' a maximum probability scheme is applied: requires G binary classifiers.
#' If type is \code{"majority"}  (only for multicalss classification G > 2) 
#' a voting scheme is applied: requires  G (G - 1) / 2 binary classifiers.
#' @param size number of units in the hidden layer. Can be zero if there are skip-layer units.
#' @param laplace value used for Laplace smoothing (additive smoothing). Defaults to 0 (no Laplace smoothing).
#' @param \dots Further arguments passed to or from other methods.
#' @return Return \code{classif} object plus:
#' \itemize{
#' \item \code{formula}: formula.
#' \item \code{data}: List that containing the variables in the model.
#' \item \code{group}: Factor of length \emph{n}. 
#' \item \code{group.est}: Estimated vector groups.
#' \item \code{prob.classification}: Probability of correct classification by group.
#' \item \code{prob.group}: Matrix of predicted class probabilities. For each
#' functional point shows the probability of each possible group membership.
#' \item \code{max.prob}: Highest probability of correct classification.
#' \item\code{type}:  Type of classification scheme: 1 vs all  or majority voting.
#' \item \code{fit}: list of binary classification fitted models.
#' }
#' @author Febrero-Bande, M. and Oviedo de la Fuente, M.
#' @seealso See Also as: \link[rpart]{rpart}.\cr Alternative method:
#' \code{\link{classif.np}}, \code{\link{classif.glm}},
#' \code{\link{classif.gsam}} and \code{\link{classif.gkam}}.
#' @references Ramsay, James O., and Silverman, Bernard W. (2006), 
#' \emph{Functional Data Analysis}, 2nd ed., Springer, New York. 
#' 
#' McCullagh and Nelder (1989), \emph{Generalized Linear Models} 2nd ed. Chapman and Hall.
#' 
#' Venables, W. N. and Ripley, B. D. (2002) \emph{Modern Applied Statistics
#' with S}, New York: Springer.  %Wood (2001) mgcv:GAMs and Generalized Ridge
#' Regression for R. R News 1(2):20-25
#' @keywords classif
#' @examples 
#' \dontrun{
#' data(phoneme)
#' mlearn<-phoneme[["learn"]]
#' glearn<-phoneme[["classlearn"]]
#' mtest<-phoneme[["test"]]
#' gtest<-phoneme[["classtest"]]
#' dataf<-data.frame(glearn)
#' dat=ldata("df"=dataf,"x"=mlearn)
#' a1<-classif.rpart(glearn~x,data=dat)
#' a2<-classif.nnet(glearn~x,data=dat)
#' a3<-classif.gbm(glearn~x,data=dat)
#' a4<-classif.randomForest(glearn~x,data=dat)
#' a5<-classif.cv.glmnet(glearn~x,data=dat)
#' newdat<-list("x"=mtest)
#' p1<-predict(a1,newdat,type="class")
#' p2<-predict(a2,newdat,type="class")
#' p3<-predict(a3,newdat,type="class")
#' p4<-predict(a4,newdat,type="class")
#' p5<-predict(a5,newdat,type="class")
#' mean(p1==gtest);mean(p2==gtest);mean(p3==gtest)
#' mean(p4==gtest);mean(p5==gtest)
#' }
#' 
#' @rdname classif.ML
#' @export classif.nnet
classif.nnet=function(formula, data, basis.x=NULL 
                      ,weights = "equal", size
                      # subset, na.action =na.omit, scale = TRUE
                      ,...) 
{
  rqr <- "nnet"
  if (!(rqr %in% rownames(installed.packages()))) {
    stop("Please install package 'nnet'") }
  
  #require(eval(rqr)[1], quietly = TRUE, warn.conflicts = FALSE)
  suppressWarnings(rqr2<-require(eval(rqr), 
                                 character.only = TRUE,quietly = TRUE, 
                                 warn.conflicts = FALSE))
  if (!rqr2) 
    stop("Please, load the namespace of the package for  method")
  
  prob=0.5
  C <- match.call()  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula","data","basis.x","weights","size"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  #mf[[1L]] <- quote(stats::model.frame)
  #mf <- eval(mf, parent.frame())
  #if (method == "model.frame")     return(mf)
  tf <- terms.formula(formula)
  terms <- attr(tf, "term.labels")
  nt <- length(terms)
  if (attr(tf, "response") > 0) {
    response <- as.character(attr(tf, "variables")[2])
    pf <- rf <- paste(response, "~", sep = "")
  } else pf <- rf <- "~"
  vtab<-rownames(attr(tf,"factors"))
  vnf=intersect(terms,names(data$df))
  # vnf2=intersect(vtab[-1],names(data$df)[-1])
  vfunc2=setdiff(terms,vnf)
  vint=setdiff(terms,vtab)
  vfunc=setdiff(vfunc2,vint)
  off<-attr(tf,"offset")
  name.coef=nam=beta.l=list()
  group <- y <- data$df[[response]]
  
  # 2019/04/24
  out.func <- fdata2model(vfunc,vnf,response, data, basis.x = basis.x ,pf = pf ,tf = tf)  
  pf <- out.func$pf          
  basis.x <- out.func$basis.x
  XX <- out.func$XX
  basis.list <- out.func$vs.list #basis.list
  mean.list=out.func$mean.list
  rm(out.func)
  n<- ndatos <-NROW(XX)
  
  # if (!is.numeric(weights))      stop("'weights' must be a numeric vector")
  if (is.character(weights)) {
    weights<-weights4class(y,type=weights)
  } else {
    if (length(weights)!=n) 
      stop("length weights != length response")
  }
  if (any(weights < 0)) 
    stop("negative weights not allowed")
  
  #  if (length(method)>1) method=method[1]
  #  if (missing(par.method))      par.method=list()
  #  par.method<-c(list(formula=pf, data=XX),par.method)
  #  if (length(vfunc)==0 & length(vnf)==0)      {
  #   par.method$pf<-as.formula(paste(pf,1,sep=""))
  #   z=do.call(method,par.method)
  #   class(z)<-c(class(z),"classif")
  #   z$formula.ini=pf
  #   z$XX=XX
  #   z$data<-data
  #   return(z)
  # }  
  
  #   par.method$size <- 4
  #   par.method$trace<- FALSE  }
  #par.method=list()
  #par.method<-c(list(formula=pf, data=XX),par.method)
  #par.method$weights= wt
  par.method <- as.list(substitute(list(...)))[-1L]
  if (missing(size)) par.method$size =4
  else  par.method$size <- size
  par.method<-c(list(formula=pf, data=XX,weights=weights),par.method)
  #par.method<-c(list(formula=pf, data=XX,weights=wt),par.method)
  z=do.call("nnet",par.method)
  
  #z=do.call("nnet",par.method)
  #if (missing(size)) size =4
  #z= nnet(formula=pf, data=XX, weights= wt,size = size
  #, subset, na.action = na.rpart, method
  #,model = FALSE, x = FALSE, y = TRUE, parms, control, cost,
  #         ,...) 
  
  
  out<-list()
  out$formula.ini=formula
  out$data=data
  out$XX=XX
  out$C <- C[1:2]
  out$prob <- prob
  out$group <- group
  # if (method=="randomForest")    out$group.est = z$predicted
  #if (method=="svm")      out$group.est = z$fitted
  #  if (method=="rpart")
  #    out$group.est <- predict(object = z, type = "class")
  #  if (method=="nnet"){
  out$group.est <- suppressWarnings(predict(object = z,type = "class"))
  out$group.est <- factor(out$group.est ,levels=levels(group))
  #  }
  out$max.prob <- mean(group==out$group.est) 
  out$fit <- z
  out$basis.x=basis.x
  out$mean=mean.list
  out$formula=pf
  out$basis.list=basis.list
  #out$method <- method
  #out$par.method <- par.method
  tab <- table(out$group.est,group)
  ny <- levels(y)
  prob2<-prob1 <- ngroup <- nlevels(y)
  prob.group <- array(NA, dim = c(ndatos, ngroup))
  prob.group <- prob.group/apply(prob.group, 1, sum)
  for (i in 1:ngroup) {
    prob1[i] = tab[i, i]/sum(tab[, i])
  }
  names(prob1) <- z$levels
  colnames(prob.group) <- z$levels
  out$prob.classification <- prob1
  #class(out)<-c("classif",class(z))
  class(out) <- "classif"
  out
}

# @export classif.multinom
classif.multinom=function(formula, data, basis.x=NULL 
                      ,weights= "equal"#, size
                      # subset, na.action =na.omit, scale = TRUE
                      ,...) 
{
  rqr <- "nnet"
  if (!(rqr %in% rownames(installed.packages()))) {
    stop("Please install package 'nnet'") }
  
  #require(eval(rqr)[1], quietly = TRUE, warn.conflicts = FALSE)
  #requireNamespace("multinom", quietly = TRUE)
  suppressWarnings(rqr2<-require(eval(rqr), 
                                 character.only = TRUE,quietly = TRUE, 
                                 warn.conflicts = FALSE))
  if (!rqr2) 
    stop("Please, load the namespace of the package for  method")
  
  prob=0.5
  C <- match.call()  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula","data","basis.x","weights"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  #mf[[1L]] <- quote(stats::model.frame)
  #mf <- eval(mf, parent.frame())
  #if (method == "model.frame")     return(mf)
  tf <- terms.formula(formula)
  terms <- attr(tf, "term.labels")
  nt <- length(terms)
  if (attr(tf, "response") > 0) {
    response <- as.character(attr(tf, "variables")[2])
    pf <- rf <- paste(response, "~", sep = "")
  } else pf <- rf <- "~"
  vtab<-rownames(attr(tf,"factors"))
  vnf=intersect(terms,names(data$df))
  # vnf2=intersect(vtab[-1],names(data$df)[-1])
  vfunc2=setdiff(terms,vnf)
  vint=setdiff(terms,vtab)
  vfunc=setdiff(vfunc2,vint)
  off<-attr(tf,"offset")
  name.coef=nam=beta.l=list()
  group <- y <- data$df[[response]]
  
  # 2019/04/24
  out.func <- fdata2model(vfunc,vnf,response, data, basis.x = basis.x ,pf = pf ,tf = tf)  
  pf <- out.func$pf          
  basis.x <- out.func$basis.x
  XX <- out.func$XX
  basis.list <- out.func$basis.list
  mean.list=out.func$mean.list
  rm(out.func)
  n <- ndatos <- NROW(XX)
  if (is.character(weights)) {
    weights<-weights4class(y,type=weights)
  } else {
    if (length(weights)!=n) 
      stop("length weights != length response")
  }
  if (any(weights < 0)) 
    stop("negative weights not allowed")
  
  par.method <- as.list(substitute(list(...)))[-1L]
# print(par.method)  
  # print(names(par.method))
  par.method<-c(list(formula=pf, data=XX,weights=weights),par.method)
  #par.method<-c(list(formula=pf, data=XX,weights=wt),par.method)
  
  # print(names(par.method))
  z=do.call("multinom",par.method)
  
    out<-list()
  out$formula.ini=formula
  out$data=data
  out$XX=XX
  out$C <- C[1:2]
  out$prob <- prob
  out$group <- group
  
  out$group.est <- predict(object = z,type = "class")
  out$group.est2 <- factor(out$group.est ,levels=levels(group))
  
  out$max.prob <- mean(group==out$group.est) 
  out$fit <- z
  out$basis.x=basis.x
  out$mean=mean.list
  out$formula=pf
  out$basis.list=basis.list
  #out$method <- method
  #out$par.method <- par.method
  tab <- table(out$group.est,group)
  ny <- levels(y)
  prob2<-prob1 <- ngroup <- nlevels(y)
  prob.group <- array(NA, dim = c(ndatos, ngroup))
  prob.group <- prob.group/apply(prob.group, 1, sum)
  for (i in 1:ngroup) {
    prob1[i] = tab[i, i]/sum(tab[, i])
  }
  names(prob1) <- z$levels
  colnames(prob.group) <- z$levels
  out$prob.classification <- prob1
  #out$group.pred <- out$group.est
  #class(out)<-c("classif",class(z))
  class(out)<-"classif"
  out
}


#' @rdname classif.ML
#' @export classif.rpart 
classif.rpart=function(formula, data, basis.x=NULL ,weights="equal",type="1vsall",...) 
{
  #print("rpart")
  rqr <- "rpart"
  if (!(rqr %in% rownames(installed.packages()))) {
    stop("Please install package 'rpart'") }
  
  #require(eval(rqr)[1], quietly = TRUE, warn.conflicts = FALSE)
  suppressWarnings(rqr2<-require(eval(rqr), 
                                 character.only = TRUE,quietly = TRUE, 
                                 warn.conflicts = FALSE))
  if (!rqr2) 
    stop("Please, load the namespace of the package for  method")
  
  prob=0.5
  C <- match.call()  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula","data","basis.x","weights","type"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  #mf <- eval(mf, parent.frame())
  #if (method == "model.frame")     return(mf)
  
  tf <- terms.formula(formula)
  terms <- attr(tf, "term.labels")
  nt <- length(terms)
  if (attr(tf, "response") > 0) {
    response <- as.character(attr(tf, "variables")[2])
    pf <- rf <- paste(response, "~", sep = "")
  } else pf <- rf <- "~"
  vtab<-rownames(attr(tf,"factors"))
  vnf=intersect(terms,names(data$df))
  # vnf2=intersect(vtab[-1],names(data$df)[-1])
  vfunc2=setdiff(terms,vnf)
  vint=setdiff(terms,vtab)
  vfunc=setdiff(vfunc2,vint)
  off<-attr(tf,"offset")
  name.coef=nam=beta.l=list()
  y <- data$df[[response]]
  lev <- levels(y)
  prob2<-prob1 <- ny <- length(lev)
  # 2019/04/24
  out.func <- fdata2model(vfunc,vnf,response, data, basis.x = basis.x ,pf = pf ,tf = tf)  
  pf <- out.func$pf          
  basis.x <- out.func$basis.x
  XX <- out.func$XX
  basis.list <- out.func$vs.list
  mean.list=out.func$mean.list
  rm(out.func)
  
  n <- NROW(XX)
  if (is.character(weights)) {
    weights<-weights4class(y,type=weights)
  } else {
    if (length(weights)!=n) 
      stop("length weights != length response")
  }
  if (any(weights < 0)) 
    stop("negative weights not allowed")
  
  par.method <- as.list(substitute(list(...)))[-1L]
  par.method<-c(list(formula=pf, data=XX,weights=weights),par.method)
  
  
  out<-list()
  out$formula.ini=formula
  out$data=data
  out$XX=XX
  out$C <- C[1:2]
  if (type == "majority" |  ny==2){
    # print("majority")
    z=do.call("rpart",par.method)
    out$fit<-z
    out$prob.group <- predict(object = z, type = "prob")
    out$group.est  <- predict(object = z, type = "class")
  } 
  
  
  else { # One vs Other
    # 2019/08/30
    #print("One vs Other")
    prob.group<-matrix(NA,n,ny)
    colnames(prob.group)<-lev
    z<-list()
    for (i in 1:ny) {
      igroup  <- y==lev[i]
      newy<-ifelse(igroup, 0,1)
      par.method$data[,response]<-factor(newy)
      newx<- par.method$data
      z[[i]] <-  do.call("rpart",par.method)
      prob.group[,i] <- predict(object = z[[i]], type = "prob")[,1]
    }
    out$prob.group <- prob.group
    out2glm<-classifKgroups(y,prob.group,lev) # hacer una par prob<0 y >0
    out$group.est = out2glm$yest
    out$fit <- z
  }
  out$prob <- prob
  out$group <- y
  out$max.prob <- mean(y==out$group.est,na.rm=T) 
  out$fit <- z
  out$basis.x=basis.x
  out$mean=mean.list
  out$formula=pf
  out$basis.list=basis.list
  tab <- table(out$group.est,out$group)
  prob.group <- array(NA, dim = c(n, ny))
  prob.group <- prob.group/apply(prob.group, 1, sum)
  for (i in 1:ny) {
    prob1[i] = tab[i, i]/sum(tab[, i])
  }
  names(prob1) <- z$levels
  colnames(prob.group) <- z$levels
  out$prob.classification <- prob1
  out$type <- type
  # out$group.pred <- out$group.est
  # class(out)<-c("classif",class(z))
  class(out) <- c("classif")
  out
}

#' @rdname classif.ML
#' @export classif.svm
classif.svm=function(formula, data, basis.x=NULL 
                     , weights="equal",type="1vsall",...)
{
  
  
  rqr <- "e1071"
  if (!(rqr %in% rownames(installed.packages()))) {
    stop("Please install package 'e1071'") }
  
  #require(eval(rqr)[1], quietly = TRUE, warn.conflicts = FALSE)
  suppressWarnings(rqr2<-require(eval(rqr), 
                                 character.only = TRUE,quietly = TRUE, 
                                 warn.conflicts = FALSE))
  if (!rqr2) 
    stop("Please, load the namespace of the package for  method")
  
  prob=0.5
  C <- match.call()  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula","data","basis.x","weights"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  
  tf <- terms.formula(formula)
  terms <- attr(tf, "term.labels")
  nt <- length(terms)
  if (attr(tf, "response") > 0) {
    response <- as.character(attr(tf, "variables")[2])
    pf <- rf <- paste(response, "~", sep = "")
  } else pf <- rf <- "~"
  vtab<-rownames(attr(tf,"factors"))
  vnf=intersect(terms,names(data$df))
  # vnf2=intersect(vtab[-1],names(data$df)[-1])
  vfunc2=setdiff(terms,vnf)
  vint=setdiff(terms,vtab)
  vfunc=setdiff(vfunc2,vint)
  off<-attr(tf,"offset")
  name.coef=nam=beta.l=list()
  group <- y <- data$df[[response]]
  lev <- levels(y)
  ny <- nlevels(y)
  # 2019/04/24
  out.func <- fdata2model(vfunc,vnf,response, data, basis.x = basis.x ,pf = pf ,tf = tf)  
  pf <- out.func$pf          
  basis.x <- out.func$basis.x
  XX <- out.func$XX
  basis.list <- out.func$vs.list
  mean.list=out.func$mean.list
  rm(out.func)
  n <- ndatos <- NROW(XX)
  par.method <- as.list(substitute(list(...)))[-1L]
  
  if (is.character(weights)) {
    weights<-weights4class(y,type=weights)
  }
  if (any(weights < 0)) 
    stop("negative weights not allowed")
  if (length(weights)==length(lev)){
    class.weights <- weights
    if (is.null(names(class.weights))) {
      names(class.weights) <- lev}
  } else{
    #    print("weights como weigths")                  
    if (length(weights)==ndatos){
      wtab<-tapply(weights,y,mean)
      ii<-wtab/sum(wtab)
      names(ii)<- lev
      class.weights <- ii
    } 
  }
  par.method<-c(list(formula=pf, data=XX,class.weights=class.weights,
                     probability=TRUE, fitted=T),par.method)
  
  
  #par.method$probability<-TRUE
  out<-list()
  out$formula.ini=formula
  out$data=data
  out$class.weights<-class.weights
  out$XX=XX
  out$C <- C[1:2]
  
  if (type == "majority" |  ny==2){
    z=do.call("svm",par.method)
    out$fit<-z
    out$group.est = z$fitted
    out$prob.group<-    attributes(predict(z,XX,desicion.values=T,  probability=T))$probabilities
    #  out$fit$call<-z$call[1:2]
    #z= svm(formula=pf, data=XX 
    #       , subset, na.action = na.action
    #       , scale = scale, class.weights= class.weights, ...) 
    
  }    else { # One vs Other
    # 2019/08/30
    # print("entra One vs Other")
    prob.group<-matrix(NA,n,ny)
    colnames(prob.group)<-lev
    #par.method$fitted<-T
    #par.method$data<-XX
    z<-list()
    #yest<-y
    #par.method$scale=F
    for (i in 1:ny) {
      igroup  <- y==lev[i]
      newy<-ifelse(igroup, 0,1)
      #weights0 <- weights
      #weights0[igroup] <- weights0[igroup]/ sum(weights0[igroup])
      #weights0[!igroup] <- weights0[!igroup]/sum(weights0[!igroup])
      #newdata$df[response]<-newy
      #a[[i]] <-suppressWarnings(fregre.glm(formula,data=newdata,family=family,weights =  weights0, basis.x=basis.x,basis.b=basis.b, CV=CV,...))
      #par.method$data<-XX[igroup,]
      par.method$data[,response]<-factor(newy)
      #       print(newy);       print(table( par.method$data[,response]))
      #par.method$scale=F
      
      #if (is.null(class.weights)) 
      #   par.method$class.weights="inverse"       else       par.method$class.weights<- c(class.weights[i],sum(class.weights[-i]))
      if (is.character(weights))  par.method$class.weights<-weights4class(par.method$data[,response],type=weights)
      else par.method$class.weights=weights
      
      names(par.method$class.weights)<-par.method$data[,response]
      #par.method$class.weights="inverse"  
      #print(  par.method$class.weights)
      newx<- par.method$data
      z[[i]] <- do.call("svm",par.method)
      aux<-predict(z[[i]],newx,desicion.values=T,  probability=T)
      aux<-attributes(aux)$probabilities
      ii<-colnames(aux)==0
      prob.group[,i]<-aux[,ii]
      out$prob.group<-prob.group
      # z[[i]]$call<- z[[i]]$call[1:2]
    }
    out2glm<-classifKgroups(y,prob.group,lev) # hacer una par prob<0 y >0
    # out$group.est = z$fitted
    #out$fit$call<-z$call[1:2]
    out$group.est =out2glm$yest
    out$fit <- z
  }
  out$prob <- prob
  out$group <- group
  # if (method=="randomForest")    out$group.est = z$predicted
  #if (method=="svm")    
  
  #  if (method=="rpart")
  #    out$group.est <- predict(object = z, type = "class")
  #  if (method=="nnet"){
  #    out$group.est <- predict(object = z,type = "class")
  #    out$group.est <- factor(out$group.est ,levels=levels(group))
  #  }
  out$max.prob <- mean(group==out$group.est) 
  out$basis.x=basis.x
  out$mean=mean.list
  out$formula=pf
  out$basis.list=basis.list
  tab <- table(out$group.est,group)
  prob2<-prob1 <- ngroup <- nlevels(y)
  prob.group <- array(NA, dim = c(ndatos, ngroup))
  prob.group <- prob.group/apply(prob.group, 1, sum)
  for (i in 1:ngroup) {
    prob1[i] = tab[i, i]/sum(tab[, i])
  }
  names(prob1) <- lev
  colnames(prob.group) <- lev
  out$prob.classification <- prob1
  out$type <- type
  #out$prob.group<-prob.group
  #out$group.pred <- out$group.est
  #class(out)<-c("classif",class(z))
  class(out) <- "classif"
 # out$prob.group <-  predict(out ,type="response")
  
  out
}


#' @rdname classif.ML
#' @export classif.ksvm
classif.ksvm=function(formula, data, basis.x=NULL ,weights = "equal",...){
  # For multiclass-classification with k classes, k > 2, ksvm uses the 'one-against-one'-approach, 
  # in which k(k-1)/2 binary classifiers are trained; the appropriate class is found by a voting scheme.
  
  rqr <- "personalized"
  if (!(rqr %in% rownames(installed.packages()))) {
    stop("Please install package 'personalized'") }
  #require(eval(rqr)[1], quietly = TRUE, warn.conflicts = FALSE)
  suppressWarnings(rqr2<-require(eval(rqr), 
                                 character.only = TRUE,quietly = TRUE, 
                                 warn.conflicts = FALSE))
  if (!rqr2) 
    stop("Please, load the namespace of the package for method")
  prob=0.5
  C <- match.call()  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula","data","basis.x","weights" ), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  tf <- terms.formula(formula)
  terms <- attr(tf, "term.labels")
  nt <- length(terms)
  if (attr(tf, "response") > 0) {
    response <- as.character(attr(tf, "variables")[2])
    pf <- rf <- paste(response, "~", sep = "")
  } else pf <- rf <- "~"
  vtab<-rownames(attr(tf,"factors"))
  vnf=intersect(terms,names(data$df))
  # vnf2=intersect(vtab[-1],names(data$df)[-1])
  vfunc2=setdiff(terms,vnf)
  vint=setdiff(terms,vtab)
  vfunc=setdiff(vfunc2,vint)
  off<-attr(tf,"offset")
  name.coef=nam=beta.l=list()
  group <- y <- data$df[[response]]
  if (!is.factor(y)) y <-factor(y)
  lev <- levels(y)
  # 2019/04/24
  out.func <- fdata2model(vfunc,vnf,response, data, basis.x = basis.x ,pf = pf ,tf = tf)  
  pf <- out.func$pf          
  basis.x <- out.func$basis.x
  XX <- out.func$XX
  basis.list <- out.func$vs.list
  mean.list=out.func$mean.list
  rm(out.func)
  par.method <- as.list(substitute(list(...)))[-1L]
  n <- ndatos <- NROW(XX)
  if (is.character(weights)) {
    weights<-weights4class(y,type=weights)
  } else {
    if (length(weights)!=n) 
      stop("length weights != length response")
  }
  if (any(weights < 0)) 
    stop("negative weights not allowed")
  
  
  out<-list()
  out$formula.ini=formula
  out$data=data
  prob2<-prob1 <- ngroup <- nlevels(y)
  if (ngroup == 2) {
    #newy <- ifelse(y == lev[1], 1, 0)
    #newdata$df$y <- newy
    #a[[1]] <- suppressWarnings(fregre.glm(formula, data = newdata, 
    #                                      family = family, weights = weights, 
    #                                        basis.x = basis.x, basis.b = basis.b,CV = CV, ...))
    #out2glm <- classif2groups(a,y,prob,lev)
    par.method<-c(list(y=y,x=XX[,-1],weights=weights),par.method)
    z=do.call("weighted.ksvm",par.method)
    pr <- predict(z, newx = as.matrix(XX[,-1]))
    out$group.est = factor(pr,labels=lev)
    out$fit <- list(z)
  }   else {
    a<-list()
    cvot<-combn(ngroup,2)
    nvot<-ncol(cvot)
    votos<-matrix(0,ndatos,ngroup)
    colnames(votos) <- lev
    b0<-list()
    par.method<-c(list(y=y,x=XX[,-1],weights=weights),par.method)
    for (ivot in 1:nvot) {  
      #print(ivot)
      ind1 <- y==lev[cvot[1,ivot]]
      ind2 <- y==lev[cvot[2,ivot]] 
      i2a2<-which(ind1 | ind2)
      newy<-rep(NA,ndatos)   
      newy[ind1 ]<- 1
      newy[ind2 ]<- -1
      #newdata$df[response] <- newy
      par.method$y=newy[i2a2]
      par.method$x=XX[i2a2,-1]
      par.method$weights=weights[i2a2]
      a[[ivot]]<-do.call("weighted.ksvm",par.method)
      a[[ivot]]$group.est <- predict(a[[ivot]], newx = as.matrix(XX[i2a2,-1]))
      #  suppressWarnings(fregre.glm(formula,data=newdata,family=family, weights =  weights
      #                                       ,basis.x=basis.x,basis.b=basis.b,CV=CV,subset = i2a2,...))
      
      prob.log <- a[[ivot]]$group.est == 1
      votos[i2a2, cvot[1,ivot]] <- votos[i2a2, cvot[1,ivot]] + as.numeric(prob.log)
      votos[i2a2, cvot[2,ivot]] <- votos[i2a2, cvot[2,ivot]] + as.numeric(!prob.log)
    }
    #out$group.est = factor(pr,labels=lev)
    out2glm<-classifKgroups(y,votos,lev)
    out<-c(out,out2glm)
    out$votos <- votos
    out$group.est = out2glm$yest
    out$fit <- a
  } 
  #print(names(z)); 
  #print("despues wei.ksvm")  
  out$weights <- weights
  out$XX=XX
  out$C <- C[1:2]
  out$prob <- prob
  out$group <- group
  # out$group.est = z$fitted
  #  if (method=="rpart")
  #    out$group.est <- predict(object = z, type = "class")
  #  if (method=="nnet"){
  #    out$group.est <- predict(object = z,type = "class")
  #    out$group.est <- factor(out$group.est ,levels=levels(group))
  #  }
  out$max.prob <- mean(group==out$group.est) 
  out$basis.x=basis.x
  out$mean=mean.list
  out$formula=pf
  out$basis.list=basis.list
  #out$type=type
  tab <- table(out$group.est,group)
  prob2<-prob1 <- ngroup <- nlevels(y)
  for (i in 1:ngroup) {
    prob1[i] = tab[i, i]/sum(tab[, i])
  }
  names(prob1) <- lev
  out$prob.classification <- prob1
  #out$group.pred <- out$group.est
  #class(out)<-c("classif",class(z))
  class(out) <- "classif"
  out
}

#' @rdname classif.ML
#' @export classif.randomForest
classif.randomForest=function(formula, data, basis.x=NULL, 
                              weights = "equal", type = "1vsall",...) 
{
  rqr <- "randomForest"
  if (!(rqr %in% rownames(installed.packages()))) {
    stop("Please install package 'randomForest'") }
  
  #require(eval(rqr)[1], quietly = TRUE, warn.conflicts = FALSE)
  suppressWarnings(rqr2<-require(eval(rqr), 
                                 character.only = TRUE,quietly = TRUE, 
                                 warn.conflicts = FALSE))
  if (!rqr2) 
    stop("Please, load the namespace of the package for  method")
  
  prob=0.5
  C <- match.call()  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula","data","basis.x","weights","type"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  #mf[[1L]] <- quote(stats::model.frame)
  #mf <- eval(mf, parent.frame())
  #if (method == "model.frame")     return(mf)
  tf <- terms.formula(formula)
  terms <- attr(tf, "term.labels")
  nt <- length(terms)
  if (attr(tf, "response") > 0) {
    response <- as.character(attr(tf, "variables")[2])
    pf <- rf <- paste(response, "~", sep = "")
  } else pf <- rf <- "~"
  vtab<-rownames(attr(tf,"factors"))
  vnf=intersect(terms,names(data$df))
  # vnf2=intersect(vtab[-1],names(data$df)[-1])
  vfunc2=setdiff(terms,vnf)
  vint=setdiff(terms,vtab)
  vfunc=setdiff(vfunc2,vint)
  off<-attr(tf,"offset")
  name.coef=nam=beta.l=list()
  y <- data$df[[response]]
  lev <- levels(y)
  prob2<-prob1 <- ny <- nlevels(y)
  
  # 2019/04/24
  out.func <- fdata2model(vfunc,vnf,response, data, basis.x = basis.x ,pf = pf ,tf = tf)  
  pf <- out.func$pf          
  basis.x <- out.func$basis.x
  XX <- out.func$XX
  basis.list <- out.func$vs.list
  mean.list=out.func$mean.list
  rm(out.func)
  n <- NROW(XX)
  
  par.method <- as.list(substitute(list(...)))[-1L]
  if (weights[1] =="equal") 
    class.weights=NULL
  
  if (weights[1] == "inverse")   {
    #print("inverse")
    weights<-weights4class(y,type=weights)
    wtab<-tapply(weights,y,mean)
    class.weights<-wtab/sum(wtab)
    names(class.weights)<- lev
  }
  
  if (any(weights < 0)) 
    stop("negative weights not allowed")
  
  if (length(weights)==ny){
    #print("weights como clas lev")
    class.weights <- weights
    names(class.weights) <- lev
  } 
  # if (is.null(names(class.weights)))     names(class.weights) <- lev
  #if (!is.null(par.method$classwt))     class.weights=classwt
  
  par.method <- c(list(formula=pf, data=XX,classwt=class.weights),par.method)
  #par.method<-c(list(formula=pf, data=XX),par.method)
  if (is.null(par.method$votes)) par.method$votes=TRUE
  out<-list()
  out$formula.ini=formula
  out$data=data
  out$XX=XX
  out$C <- C[1:2]
  
  if (type == "majority" |  ny==2){
    z=do.call("randomForest",par.method)
    out$fit<-z
    out$group.est = z$predicted
    out$prob.group <- z$votes
    
  }    else { # One vs Other
    # 2019/08/30
    prob.group<-matrix(NA,n,ny)
    colnames(prob.group)<-lev
    z<-list()
    for (i in 1:ny) {
      igroup  <- y==lev[i]
      newy<-ifelse(igroup, 0,1)
      par.method$data[,response]<-factor(newy)
      newx<- par.method$data
      z[[i]] <-  do.call("randomForest",par.method)
      prob.group[,i] <- z[[i]]$votes[,1]
    }
    out$prob.group <- prob.group
    out2glm<-classifKgroups(y,prob.group,lev) # hacer una par prob<0 y >0
    out$group.est = out2glm$yest
    out$fit <- z
  }
  out$prob <- prob
  out$group <- y
  out$max.prob <- mean(y==out$group.est,na.rm=T) 
  out$fit <- z
  out$basis.x=basis.x
  out$mean=mean.list
  out$formula=pf
  out$basis.list=basis.list
  #out$method <- method
  #out$par.method <- par.method
  
  tab <- table(out$group.est,y)
  prob.group2 <- array(NA, dim = c(n, ny))
  prob.group2 <- prob.group2/apply(prob.group2, 1, sum)
  for (i in 1:ny) {
    prob1[i] = tab[i, i]/sum(tab[, i])
  }
  names(prob1) <- z$levels
  #colnames(prob.group) <- z$levels
  out$prob.classification <- prob1
  out$type <- type
  # If norm.votes=TRUE, the fraction is given, which can be taken as predicted probabilities for the classes.
  
  #out$group.pred <- out$group.est
  #class(out)<-c("classif",class(z))
  class(out) <- "classif"
  out
}


#' @rdname classif.ML
#' @export classif.lda
classif.lda=function(formula, data, basis.x=NULL 
                     , weights="equal",type="1vsall",...)
{
  rqr <- "MASS"
  if (!(rqr %in% rownames(installed.packages()))) {
    stop("Please install package 'MASS'") }
  
  #require(eval(rqr)[1], quietly = TRUE, warn.conflicts = FALSE)
  suppressWarnings(rqr2<-require(eval(rqr), 
                                 character.only = TRUE,quietly = TRUE, 
                                 warn.conflicts = FALSE))
  if (!rqr2) 
    stop("Please, load the namespace of the package for  method")
  
  prob=0.5
  C <- match.call()  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula","data","basis.x","weights"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  
  tf <- terms.formula(formula)
  terms <- attr(tf, "term.labels")
  nt <- length(terms)
  if (attr(tf, "response") > 0) {
    response <- as.character(attr(tf, "variables")[2])
    pf <- rf <- paste(response, "~", sep = "")
  } else pf <- rf <- "~"
  vtab<-rownames(attr(tf,"factors"))
  vnf=intersect(terms,names(data$df))
  # vnf2=intersect(vtab[-1],names(data$df)[-1])
  vfunc2=setdiff(terms,vnf)
  vint=setdiff(terms,vtab)
  vfunc=setdiff(vfunc2,vint)
  off<-attr(tf,"offset")
  name.coef=nam=beta.l=list()
  group <- y <- data$df[[response]]
  lev <- levels(y)
  prob2<-prob1 <- ny <- nlevels(y)
  # 2019/04/24
  out.func <- fdata2model(vfunc,vnf,response, data, basis.x = basis.x ,pf = pf ,tf = tf)  
  pf <- out.func$pf          
  basis.x <- out.func$basis.x
  XX <- out.func$XX
  basis.list <- out.func$vs.list
  mean.list=out.func$mean.list
  rm(out.func)
  n <- ndatos <- NROW(XX)
  par.method <- as.list(substitute(list(...)))[-1L]
  
  if (is.character(weights)) {
    weights<-weights4class(y,type=weights)
  }
 # print(class(weights))
  if (any(weights < 0)) 
    stop("negative weights not allowed")
  if (length(weights)==length(lev)){
    class.weights <- weights
    #if (is.null(names(class.weights))) {      names(class.weights) <- lev      }
  } else{
       # print("weights como weigths")                  
    if (length(weights)==ndatos){
      wtab<-tapply(weights,y,mean)
      ii<-wtab/sum(wtab)
      names(ii)<- lev
      class.weights <- ii
    } 
  }
  
  par.method<-c(list(formula=pf, data=XX,prior=as.vector(class.weights)),par.method)
#  par.method$prior<-as.vector(par.method$prior)
  #par.method$probability<-TRUE
  out<-list()
  out$formula.ini=formula
  out$data=data
  
  out$prior<-class.weights
#print(out$prior)   # eliminar argumento (se substituye por weigths)
  out$XX=XX
  out$C <- C[1:2]
  
  # if (type == "majority" &  ny>2){
  #   z=do.call("svm",par.method)
  #   out$fit<-z
  #   out$group.est = z$fitted
  # } 
  
  a<-list()
  if (type == "majority"){
  #  print("majority")
    cvot<-combn(ny,2)
    nvot<-ncol(cvot)
    votos<-matrix(0,n,ny)
    colnames(votos) <- lev
    b0<-list()
    for (ivot in 1:nvot) {  
      ind1 <- y==lev[cvot[1,ivot]]
      ind2 <- y==lev[cvot[2,ivot]] 
      i2a2<-which(ind1 | ind2)
      newy<-rep(NA,n)   
      newy[ind1 ]<- 1
      newy[ind2 ]<- 0
      par.method$prior<-class.weights[cvot[,ivot]]/sum(class.weights[cvot[,ivot]])
      par.method$data<-XX[i2a2,]
      par.method$data[response]<-newy[i2a2]
      a[[ivot]]<- do.call("lda",par.method[1:2])
      newx<-XX[i2a2,-1]
      aux<-predict(a[[ivot]],newx)        
      prob.log <- aux$posterior[,2]  > prob
      votos[i2a2, cvot[1,ivot]] <- votos[i2a2, cvot[1,ivot]] + as.numeric(prob.log)
      votos[i2a2, cvot[2,ivot]] <- votos[i2a2, cvot[2,ivot]] + as.numeric(!prob.log)
    }
    out2glm<-classifKgroups(y,votos/(ny-1),lev)
    out$group.est<-out2glm$yest
    out$prob.group<-out2glm$prob.group
    out$votes <- votos
    out$fit<-a
#print("sale majority")
  }
  else { # One vs Other
#    print("one vs Other")
      z <- do.call("lda",par.method)
      aux<-predict(z,par.method$data)
      out$group.est<-aux$class
      prob.group<-out$prob.group<-aux$posterior
      out$fit <- list(z)
  }
  out$prob <- prob
  out$group <- y

  out$max.prob <- mean(group==out$group.est) 
  out$basis.x=basis.x
  out$mean=mean.list
  out$formula=pf
  out$basis.list=basis.list
  tab <- table(out$group.est,y)

  prob.group <- array(NA, dim = c(n, ny))
  prob.group <- prob.group/apply(prob.group, 1, sum)
  for (i in 1:ny) {
    prob1[i] = tab[i, i]/sum(tab[, i])
  }
  names(prob1) <- lev
  colnames(prob.group) <- lev
  out$prob.classification <- prob1
  out$type <- type
  #out$prob.group<-prob.group
  #out$group.pred <- out$group.est
  #class(out)<-c("classif",class(z))
  class(out) <- "classif"
  out
}

#' @rdname classif.ML
#' @export classif.qda
classif.qda=function(formula, data, basis.x=NULL 
                     , weights="equal",type="1vsall",...)
{
  rqr <- "MASS"
  if (!(rqr %in% rownames(installed.packages()))) {
    stop("Please install package 'MASS'") }
  
  #require(eval(rqr)[1], quietly = TRUE, warn.conflicts = FALSE)
  suppressWarnings(rqr2<-require(eval(rqr), 
                                 character.only = TRUE,quietly = TRUE, 
                                 warn.conflicts = FALSE))
  if (!rqr2) 
    stop("Please, load the namespace of the package for  method")
  
  prob=0.5
  C <- match.call()  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula","data","basis.x","weights"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  
  tf <- terms.formula(formula)
  terms <- attr(tf, "term.labels")
  nt <- length(terms)
  if (attr(tf, "response") > 0) {
    response <- as.character(attr(tf, "variables")[2])
    pf <- rf <- paste(response, "~", sep = "")
  } else pf <- rf <- "~"
  vtab<-rownames(attr(tf,"factors"))
  vnf=intersect(terms,names(data$df))
  # vnf2=intersect(vtab[-1],names(data$df)[-1])
  vfunc2=setdiff(terms,vnf)
  vint=setdiff(terms,vtab)
  vfunc=setdiff(vfunc2,vint)
  off<-attr(tf,"offset")
  name.coef=nam=beta.l=list()
  group <- y <- data$df[[response]]
  lev <- levels(y)
  prob2<-prob1 <- ny <- nlevels(y)
  # 2019/04/24
  out.func <- fdata2model(vfunc,vnf,response, data, basis.x = basis.x ,pf = pf ,tf = tf)  
  pf <- out.func$pf          
  basis.x <- out.func$basis.x
  XX <- out.func$XX
  basis.list <- out.func$vs.list
  mean.list=out.func$mean.list
  rm(out.func)
  n <- ndatos <- NROW(XX)
  par.method <- as.list(substitute(list(...)))[-1L]
  
  if (is.character(weights)) {
    weights<-weights4class(y,type=weights)
  }
  # print(class(weights))
  if (any(weights < 0)) 
    stop("negative weights not allowed")
  if (length(weights)==length(lev)){
    class.weights <- weights
    #if (is.null(names(class.weights))) {      names(class.weights) <- lev      }
  } else{
    # print("weights como weigths")                  
    if (length(weights)==ndatos){
      wtab<-tapply(weights,y,mean)
      ii<-wtab/sum(wtab)
      names(ii)<- lev
      class.weights <- ii
    } 
  }
  
  par.method<-c(list(formula=pf, data=XX,prior=as.vector(class.weights)),par.method)
  #  par.method$prior<-as.vector(par.method$prior)
  #par.method$probability<-TRUE
  out<-list()
  out$formula.ini=formula
  out$data=data
  
  out$prior<-class.weights
  #print(out$prior)   # eliminar argumento (se substituye por weigths)
  out$XX=XX
  out$C <- C[1:2]
  
  # if (type == "majority" &  ny>2){
  #   z=do.call("svm",par.method)
  #   out$fit<-z
  #   out$group.est = z$fitted
  # } 
  
  a<-list()
  if (type == "majority"){
    #  print("majority")
    cvot<-combn(ny,2)
    nvot<-ncol(cvot)
    votos<-matrix(0,n,ny)
    colnames(votos) <- lev
    b0<-list()
    for (ivot in 1:nvot) {  
      ind1 <- y==lev[cvot[1,ivot]]
      ind2 <- y==lev[cvot[2,ivot]] 
      i2a2<-which(ind1 | ind2)
      newy<-rep(NA,n)   
      newy[ind1 ]<- 1
      newy[ind2 ]<- 0
      par.method$prior<-class.weights[cvot[,ivot]]/sum(class.weights[cvot[,ivot]])
      par.method$data<-XX[i2a2,]
      par.method$data[response]<-newy[i2a2]
      a[[ivot]]<- do.call("qda",par.method[1:2])
      newx<-XX[i2a2,-1]
      aux<-predict(a[[ivot]],newx)        
      prob.log <- aux$posterior[,2]  > prob
      votos[i2a2, cvot[1,ivot]] <- votos[i2a2, cvot[1,ivot]] + as.numeric(prob.log)
      votos[i2a2, cvot[2,ivot]] <- votos[i2a2, cvot[2,ivot]] + as.numeric(!prob.log)
    }
    out2glm<-classifKgroups(y,votos/(ny-1),lev)
    out$group.est<-out2glm$yest
    out$prob.group<-out2glm$prob.group
    out$votes <- votos
    out$fit<-a
    #print("sale majority")
  }
  else { # One vs Other
    #    print("one vs Other")
    z <- do.call("qda",par.method)
    aux<-predict(z,par.method$data)
    out$group.est<-aux$class
    prob.group<-out$prob.group<-aux$posterior
    out$fit <- list(z)
  }
  out$prob <- prob
  out$group <- y
  
  out$max.prob <- mean(group==out$group.est) 
  out$basis.x=basis.x
  out$mean=mean.list
  out$formula=pf
  out$basis.list=basis.list
  tab <- table(out$group.est,y)
  
  prob.group <- array(NA, dim = c(n, ny))
  prob.group <- prob.group/apply(prob.group, 1, sum)
  for (i in 1:ny) {
    prob1[i] = tab[i, i]/sum(tab[, i])
  }
  names(prob1) <- lev
  colnames(prob.group) <- lev
  out$prob.classification <- prob1
  out$type <- type
  #out$prob.group<-prob.group
  #out$group.pred <- out$group.est
  #class(out)<-c("classif",class(z))
  class(out) <- "classif"
  out
}

#' @rdname classif.ML
#' @export classif.naiveBayes
classif.naiveBayes=function(formula, data, basis.x=NULL, laplace = 0,...) 
{
  rqr <- "e1071"
  if (!(rqr %in% rownames(installed.packages()))) {
    stop("Please install package 'e1071'") }
  
  #require(eval(rqr)[1], quietly = TRUE, warn.conflicts = FALSE)
  suppressWarnings(rqr2<-require(eval(rqr), 
                                 character.only = TRUE,quietly = TRUE, 
                                 warn.conflicts = FALSE))
  if (!rqr2) 
    stop("Please, load the namespace of the package for  method")
  
  prob=0.5
  C <- match.call()  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula","data","basis.x","laplace"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  
  tf <- terms.formula(formula)
  terms <- attr(tf, "term.labels")
  nt <- length(terms)
  if (attr(tf, "response") > 0) {
    response <- as.character(attr(tf, "variables")[2])
    pf <- rf <- paste(response, "~", sep = "")
  } else pf <- rf <- "~"
  vtab<-rownames(attr(tf,"factors"))
  vnf=intersect(terms,names(data$df))
  # vnf2=intersect(vtab[-1],names(data$df)[-1])
  vfunc2=setdiff(terms,vnf)
  vint=setdiff(terms,vtab)
  vfunc=setdiff(vfunc2,vint)
  off<-attr(tf,"offset")
  name.coef=nam=beta.l=list()
  group <- y <- data$df[[response]]
  lev <- levels(y)
  # 2019/04/24
  out.func <- fdata2model(vfunc,vnf,response, data, basis.x = basis.x ,pf = pf ,tf = tf)  
  pf <- out.func$pf          
  basis.x <- out.func$basis.x
  XX <- out.func$XX
  basis.list <- out.func$vs.list
  mean.list=out.func$mean.list
  rm(out.func)
  ndatos <- nrow(XX)
  #if (missing(weights)) wt = rep(1,ndatos)  else wt=weights
  #  if (length(method)>1) method=method[1]
  #  if (missing(par.method))      par.method=list()
  par.method <- as.list(substitute(list(...)))[-1L]
  #if (!is.null(class.weights)){
  #    if (is.null(names(class.weights))) {
  #    names(class.weights) <- lev}}
  par.method<-c(list(formula=pf, data=XX, lapalce=laplace),par.method)
  #  if (length(vfunc)==0 & length(vnf)==0)      {
  #   par.method$pf<-as.formula(paste(pf,1,sep=""))
  #   z=do.call(method,par.method)
  #   class(z)<-c(class(z),"classif")
  #   z$formula.ini=pf
  #   z$XX=XX
  #   z$data<-data
  #   return(z)
  # }  
  # if (method=="nnet" & is.null(par.method$size)){
  #   par.method$size <- 4
  #   par.method$trace<- FALSE  }
  z=do.call("naiveBayes",par.method)
  #z= svm(formula=pf, data=XX 
  #       , subset, na.action = na.action
  #       , scale = scale, class.weights= class.weights, ...) 
  out<-list()
  out$formula.ini=formula
  out$data=data
  out$XX=XX
  out$C <- C[1:2]
  out$prob <- prob
  out$group <- group
  # if (method=="randomForest")    out$group.est = z$predicted
  #if (method=="svm")    
  out$group.est <- predict(object = z, newdata=XX, type = "class")
  #  if (method=="rpart")
  #    out$group.est <- predict(object = z, type = "class")
  #  if (method=="nnet"){
  #    out$group.est <- predict(object = z,type = "class")
  #    out$group.est <- factor(out$group.est ,levels=levels(group))
  #  }
  out$max.prob <- mean(group==out$group.est) 
  out$fit <- z
  out$basis.x=basis.x
  out$mean=mean.list
  out$formula=pf
  out$basis.list=basis.list
  #out$method <- method
  #out$par.method <- par.method
  #print(out$group.est)
  #print(3)
  tab <- table(out$group.est,group)
  # print(4)
  
  prob2<-prob1 <- ngroup <- nlevels(y)
  prob.group <- array(NA, dim = c(ndatos, ngroup))
  prob.group <- prob.group/apply(prob.group, 1, sum)
  for (i in 1:ngroup) {
    prob1[i] = tab[i, i]/sum(tab[, i])
  }
  names(prob1) <- z$levels
  colnames(prob.group) <- z$levels
  out$prob.classification <- prob1
  #out$group.pred <- out$group.est
  #class(out)<-c("classif",class(z))
  class(out) <- "classif"
  out
}



# library(glmnet)
# NFOLDS = 4;
# res.glmnet = cv.glmnet( x= as.matrix(train[,-101]), y = as.factor(train[,101]),
#                         family = 'multinomial',
#                         alpha = 1,
#                         #                        grouped = TRUE,
#                         type.measure = "class",
#                         nfolds = NFOLDS,
#                         thresh = 1e-3,
#                         maxit = 1e3)   

#' @rdname classif.ML
#' @export classif.cv.glmnet
classif.cv.glmnet=function(formula, data, basis.x=NULL 
                      ,weights = "equal"
                      # subset, na.action =na.omit, scale = TRUE
                      ,...) 
{
  # data <- dat
  # formula <- formula(glearn~x)
  # basis.x=NULL 
  
  rqr <- "glmnet"
  if (!(rqr %in% rownames(installed.packages()))) {
    stop("Please install package 'glmnet'") }
  
  #require(eval(rqr)[1], quietly = TRUE, warn.conflicts = FALSE)
  suppressWarnings(rqr2<-require(eval(rqr), 
                                 character.only = TRUE,quietly = TRUE, 
                                 warn.conflicts = FALSE))
  if (!rqr2) 
    stop("Please, load the namespace of the package for  method")
  
  prob=0.5
  C <- match.call()  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula","data","basis.x","weights","size"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
 
  tf <- terms.formula(formula)
  terms <- attr(tf, "term.labels")
  nt <- length(terms)
  if (attr(tf, "response") > 0) {
    response <- as.character(attr(tf, "variables")[2])
    pf <- rf <- paste(response, "~", sep = "")
  } else pf <- rf <- "~"
  vtab<-rownames(attr(tf,"factors"))
  vnf=intersect(terms,names(data$df))
  # vnf2=intersect(vtab[-1],names(data$df)[-1])
  vfunc2=setdiff(terms,vnf)
  vint=setdiff(terms,vtab)
  vfunc=setdiff(vfunc2,vint)
  off<-attr(tf,"offset")
  name.coef=nam=beta.l=list()
  group <- y <- data$df[[response]]
  
  # 2019/04/24
  out.func <- fdata2model(vfunc,vnf,response, data, basis.x = basis.x ,pf = pf ,tf = tf)  
  pf <- out.func$pf          
  basis.x <- out.func$basis.x
  XX <- out.func$XX
  basis.list <- out.func$vs.list
  mean.list=out.func$mean.list
  rm(out.func)
  n <- ndatos <-NROW(XX)
  
  # if (!is.numeric(weights))      stop("'weights' must be a numeric vector")
  if (is.character(weights)) {
    weights<-weights4class(y,type=weights)
  } else {
    if (length(weights)!=n) 
      stop("length weights != length response")
  }
  if (any(weights < 0)) 
    stop("negative weights not allowed")
  
  #  if (length(method)>1) method=method[1]
  #  if (missing(par.method))      par.method=list()
  #  par.method<-c(list(formula=pf, data=XX),par.method)
  #  if (length(vfunc)==0 & length(vnf)==0)      {
  #   par.method$pf<-as.formula(paste(pf,1,sep=""))
  #   z=do.call(method,par.method)
  #   class(z)<-c(class(z),"classif")
  #   z$formula.ini=pf
  #   z$XX=XX
  #   z$data<-data
  #   return(z)
  # }  
  
  #   par.method$size <- 4
  #   par.method$trace<- FALSE  }
  #par.method=list()
  #par.method<-c(list(formula=pf, data=XX),par.method)
  #par.method$weights= wt
  par.method <- as.list(substitute(list(...)))[-1L]

  #par.method<-c(list(x=XX, y=y, family = "multinomial",weights=weights),par.method)
  par.method<-c(list(x=as.matrix(XX[,-1,drop=F]), y = y, family = "multinomial",weights=weights),par.method)
   
#  X no puede ser de dimensión 1!!
#  print("ML M;L")
if (NCOL(par.method$x)==1) par.method$x<-cbind(rep(1,len=n),par.method$x)
  # print(dim(par.method$x))
  z= suppressWarnings(do.call("cv.glmnet",par.method))
  out<-list()
  out$formula.ini=formula
  out$data=data
  out$XX=XX
  out$C <- C[1:2]
  out$prob <- prob
  out$group <- group
  out$prob.group <- predict(object = z, par.method$x,type = "response")[,,1]
  out$group.est <- predict(object = z, par.method$x, type = "class")
  out$group.est <- factor(out$group.est ,levels=levels(group))
  out$max.prob <- mean(group==out$group.est) 
  out$fit <- z
  out$basis.x=basis.x
  out$mean=mean.list
  out$formula=pf
  out$basis.list=basis.list
  out$weights = weights
  tab <- table(out$group.est,group)
  ny <- levels(y)
  prob2<-prob1 <- ngroup <- nlevels(y)
  #prob.group <- array(NA, dim = c(ndatos, ngroup))
  # prob.group <- prob.group/apply(prob.group, 1, sum)
  
  
  
  for (i in 1:ngroup) {
    prob1[i] = tab[i, i]/sum(tab[, i])
  }
  names(prob1) <- z$levels
  colnames(out$prob.group) <- z$levels
  out$prob.classification <- prob1
  out$fit$call<-  out$fit$call[1]
  #class(out)<-c("classif",class(z))
  class(out) <- "classif"
  return(out)
}

#' @rdname classif.ML
#' @export classif.gbm
classif.gbm=function(formula, data, basis.x=NULL 
                      ,weights = "equal",...) 
{
  rqr <- "gbm"
  if (!(rqr %in% rownames(installed.packages()))) {
    stop("Please install package 'gbm'") }
  
  #require(eval(rqr)[1], quietly = TRUE, warn.conflicts = FALSE)
  suppressWarnings(rqr2<-require(eval(rqr), 
                                 character.only = TRUE,quietly = TRUE, 
                                 warn.conflicts = FALSE))
  if (!rqr2) 
    stop("Please, load the namespace of the package for  method")
  
  prob=0.5
  C <- match.call()  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula","data","basis.x","weights"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  #mf[[1L]] <- quote(stats::model.frame)
  #mf <- eval(mf, parent.frame())
  #if (method == "model.frame")     return(mf)
  tf <- terms.formula(formula)
  terms <- attr(tf, "term.labels")
  nt <- length(terms)
  if (attr(tf, "response") > 0) {
    response <- as.character(attr(tf, "variables")[2])
    pf <- rf <- paste(response, "~", sep = "")
  } else pf <- rf <- "~"
  vtab<-rownames(attr(tf,"factors"))
  vnf=intersect(terms,names(data$df))
  # vnf2=intersect(vtab[-1],names(data$df)[-1])
  vfunc2=setdiff(terms,vnf)
  vint=setdiff(terms,vtab)
  vfunc=setdiff(vfunc2,vint)
  off<-attr(tf,"offset")
  name.coef=nam=beta.l=list()
  group <- y <- data$df[[response]]
  
  # 2019/04/24
  out.func <- fdata2model(vfunc,vnf,response, data, basis.x = basis.x ,pf = pf ,tf = tf)  
  pf <- out.func$pf          
  basis.x <- out.func$basis.x
  XX <- out.func$XX
  basis.list <- out.func$vs.list
  mean.list=out.func$mean.list
  rm(out.func)
  n<- ndatos <-NROW(XX)
  
  # if (!is.numeric(weights))      stop("'weights' must be a numeric vector")
  if (is.character(weights)) {
    weights<-weights4class(y,type=weights)
  } else {
    if (length(weights)!=n) 
      stop("length weights != length response")
  }
  if (any(weights < 0)) 
    stop("negative weights not allowed")
  
  par.method <- as.list(substitute(list(...)))[-1L]
  par.method<-c(list(formula=pf, data=XX,weights=weights,
                     distribution = "multinomial"),par.method)
  z= suppressWarnings(do.call(rqr,par.method))
  lev <- levels(group)  
  out<-list()
  out$formula.ini=formula
  out$data=data
  out$XX=XX
  out$C <- C[1:2]
  out$prob <- prob
  out$group <- group

  out$group.est <- suppressWarnings(predict(z, newdata = XX,type="response"))
  out$group.est <- apply(out$group.est,1,which.max)
  out$group.est <- factor(out$group.est ,levels=lev)
  
  out$max.prob <- mean(group==out$group.est) 
  out$fit <- z
  out$basis.x=basis.x
  out$mean=mean.list
  out$formula=pf
  out$basis.list=basis.list
  #out$method <- method
  #out$par.method <- par.method
  tab <- table(out$group.est,group)
  ny <- levels(y)
  prob2<-prob1 <- ngroup <- nlevels(y)
  prob.group <- array(NA, dim = c(ndatos, ngroup))
  prob.group <- prob.group/apply(prob.group, 1, sum)
  for (i in 1:ngroup) {
    prob1[i] = tab[i, i]/sum(tab[, i])
  }
  names(prob1) <- z$levels
  colnames(prob.group) <- z$levels
  out$prob.classification <- prob1
  out$type="majority"
  #class(out)<-c("classif",class(z))
  class(out) <- "classif"
  out
}
# Hacer gbm usando majority voting

