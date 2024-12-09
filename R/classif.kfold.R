########################################################################
kfold.B <- function(call,par.call,param.kfold,vfunc,name.param,lista,nparams
                    ,classif,par.classif,data,cost,ntrain
                    ,it,verbose,default,response,models
                    ,ifold,folds,lparam.kfold,measure){
  nfun <- length(vfunc)
  name.param <- names(param.kfold)
  if (classif == "classif.adaboost"){
    BB <- param.kfold[[name.param]]
    for (i in 1:nfun)  param.kfold[[vfunc[i]]][[name.param]]  <-  BB
    lista[[vfunc[i]]] <- BB
    nparams[i]  <- length(lista[[vfunc[i]]])
    params  <- do.call("expand.grid",lista)
    nexpand <- nrow(params)
    max.params <- apply(params,2,max)
    error <- rep(NA,nexpand)
   }
  if (missing(par.classif)) par.classif=list()
  yresp <- data$df[,response]
  lev <- levels(yresp)
  if (missing(cost)) cost = rep(1, nlevels(yresp))
  
  pred <- factor(rep(NA,len=ntrain),levels=lev)
  if (classif == "classif.bootstrap"){
    # lista <- list(smo=0B=5, nb=200, Nhull=4, Nnbh=9)
    lista <- list(smo=0,  Nhull=NULL, Nnbh=NULL) 
    if (any(name.param=="nb")){
      lista[["nb"]]  <-  param.kfold[["nb"]]
    }
    if (any(name.param=="smo")){
      lista[["smo"]]  <-  param.kfold[["smo"]]
    }
    if (any(name.param=="Nhull")){
      lista[["Nhull"]]  <-  param.kfold[["Nhull"]]
    }
    if (any(name.param=="Nnbh")){
      lista[["Nnbh"]]  <-  param.kfold[["Nnbh"]]
    }
      name.param <- names(lista)
    #BB <- param.kfold[[name.param]]
    #for (i in 1:nfun)  param.kfold[[vfunc[i]]][[name.param]]  <-  BB
    #nparams[i]  <- length(lista[[vfunc[i]]])
    params <- do.call("expand.grid",lista)
    nexpand <- nrow(params)
    max.params <- apply(params,2,max)
    error <- rep(NA,nexpand)
     # print("   ** params **  ")    ;    print(params)
    
  }
  it <- TRUE
  if (models) models.pred <- list()
  for (j in 1:nexpand){ 
    #print(j)    
    if (verbose) print(params[j,])
    cat("Param j,",j,"\n")
    for (i in ifold){ 
      #print(i)
      if (verbose) cat("j,",j,"params",as.numeric(params[j,])," kfold i,",i,"\n")
      #Segement your data by fold using the which() function 
      testIndexes <- which(folds==i,arr.ind=TRUE)
      trainData <- subset.ldata(data,folds!=i)
      ipred <- folds==i
      testData  <-   subset.ldata(data,ipred)
      par.call$data <- trainData
      if (it) par.call <- c(par.call, par.classif)
      it <- FALSE
      par.call$data <- trainData
      ########## if (name.param[1]=="B") par.call$B <- param.kfold$B[j]
      # if (classif == "classif.bootstrap"){
      #  par.boot = list(B=50, N=1000, Nhull=4, Nnbh=9) 
        #par.call$par.boot <- lista[[j]]
        #par.call$par.boot <- params[j,]
       #  }
  #par.call$classif <- classif
  par.call$par.classif <- par.classif
  if (is.null(par.call$classif)) par.call$classif <- "classif.glm"
  name.param <- names(params)
################################ 
  if (j==1 & i==1)   par.call$par.boot <- list()
  #print(par.call$par.boot)
  
    for (iparam in 1:4)
            par.call$par.boot[[name.param[iparam]]] <- params[j,iparam]
################################   
  res <- do.call(call,par.call)
      pred[ipred]  <-  predict.classif(res,testData)
      if (models) {
        nam.ji <- paste0("Params",j,"-Kfold",i)
        models.pred[[nam.ji]]  <-  res
      }
    } # fin i kfold  
    error[j]  <-  1-cat2meas(yresp, pred, measure = measure, cost = cost)
  }# fin j param
  
  names(error)  <-   paste0("param ",apply(params,1,paste0,collapse="-"))
  imin  <-   which(error==min(error))
  imin2 <- imin[1] 
  par.call$data <- data  
  if (classif == "classif.bootstrap"){
  for (iparam in 1:4){
    par.call$par.boot[[name.param[iparam]]] <-  params[imin[1],iparam]
  }
  } else   par.call$B <- BB[imin]
   output <- list(par.call=par.call,imin=imin,params=params
               ,error=error,pred=pred)
  if (models) output$models.pred=models.pred
  return(output)
}
########################################################################
kfold.aux <- function(call,par.call,param.kfold,vfunc,name.param,lista,nparams
                    ,classif,par.classif,data,cost,ntrain
                    ,it,verbose,default,response,models
                    ,ifold,folds,lparam.kfold,measure){
  nfun <- length(vfunc)
  for (i in 1:nfun){
    if (default) {
      param.kfold[[vfunc[i]]]  <-  "default"
      names(param.kfold) <- vfunc[i]
      name.param[i] <- "default"
    }     else name.param[i]  <-  names(param.kfold[[vfunc[i]]])
    lista[[vfunc[i]]] <- param.kfold[[vfunc[i]]][[1]]
    nparams[i]  <- length(lista[[vfunc[i]]])
  }
  params  <- do.call("expand.grid",lista)
  nexpand <- nrow(params)
  max.params <- apply(params,2,max)
  error <- rep(NA,nexpand)
  if (missing(par.classif)) par.classif=list()
  yresp <- data$df[,response]
  lev <- levels(yresp)
  if (missing(cost)) cost = rep(1, nlevels(yresp))

  pred <- factor(rep(NA,len=ntrain),levels=lev)
  basis.aux <- list()
  for (ifun in 1:nfun){
    ################################ PC ##########################
    if (name.param[ifun] %in% c("pc")){ 
      basis.aux[[vfunc[ifun]]]  <-  create.fdata.basis(data[[vfunc[ifun]]],
                                                     1:max.params[ifun],
                                                     type.basis=name.param[ifun])
    }
    ################################ FiXED basis ##########################
    if (name.param[1] %in% c("bspline","fourier","constant",
                             "exponential","polygonal","power")){
      basis.aux[[vfunc[ifun]]]  <-  create.fdata.basis(data[[vfunc[ifun]]],
                                                     1:max.params[ifun],
                                                     type.basis=name.param[ifun])
    }
    ################################ bandwidth ##########################
    if (name.param[1]=="h"){
      #par.call[["par.np"]][[vfunc[ifun]]][[name.param[ifun]]] <- NULL
      par.classif$par.np <- NULL # par.np
    }
  }
  it <- TRUE
  if (models) models.pred <- list()
  ############################  
  for (j in 1:nexpand){ 
    if (verbose) print(params[j,])
    cat("Param j,",j,"\n")
    for (i in ifold){ 
      if (verbose) cat("j,",j,"params",as.numeric(params[j,])," kfold i,",i,"\n")
      #Segement your data by fold using the which() function 
      testIndexes <- which(folds==i,arr.ind=TRUE)
      trainData <- subset.ldata(data,folds!=i)
      ipred <- folds==i
      testData  <-   subset.ldata(data,ipred)
      if (classif %in% c("classif.np","classif.knn","classif.kernel")){
        if (it) par.call <- par.classif
        par.call$group <- trainData$df[,response]
        par.call$fdataobj <- trainData[[vfunc[ifun]]]
        it <- FALSE
        if (name.param[1]=="h")          par.call[["h"]] <- params[j,]
        if (name.param[1]=="knn")        par.call[["knn"]] <- params[j,]
        res <- do.call(call,par.call)
        aux <- predict.classif(res,testData[[vfunc[ifun]]])
        pred[ipred] <- aux
      } else{
        par.call$data <- trainData
         if (it) par.call <- c(par.call, par.classif)
        it <- FALSE
        ################################ PC ##########################
        if (name.param[1]=="pc"){ #atajo para que tarde menos
          basis.x <- basis.aux
          for (ifun in 1:nfun){
            basis.x[[vfunc[ifun]]]$l <- basis.aux[[vfunc[ifun]]]$l[1:params[j,ifun]]
            basis.x[[vfunc[ifun]]]$basis <- basis.aux[[vfunc[ifun]]]$basis[1:params[j,ifun]]
            basis.x[[vfunc[ifun]]]$x <- basis.aux[[vfunc[ifun]]]$x[folds!=i,1:params[j,ifun],drop=F]
          }
          par.call$basis.x <- basis.x
        } 
        ################################ FiXED basis ##########################
        if (name.param[1] %in% c("bspline","fourier","constant",
                                 "exponential","polygonal","power")){
          basis.x2 <- list()
          for (ifun in 1:nfun){
            basis.x2[[vfunc[ifun]]]  <-  create.fdata.basis(trainData[[vfunc[ifun]]],
                                                          1:params[j,ifun]
                                                          ,type.basis=name.param[ifun])
          }
          if (classif %in% c("classif.gsam","classif.glm")) par.call$basis.b <- basis.x2
        }
        ################################ BANDWIDTH ##########################
        if (name.param[1]=="h"){
          for (ifun0 in 1:nfun){
            par.call[["par.np"]][[vfunc[ifun0]]] <- list("h"=params[j,ifun0])
          }
        }
        par.call$data <- trainData
        res <- do.call(call,par.call)
        pred[ipred]  <-  predict.classif(res,testData)
      } 
      if (models) {
        nam.ji <- paste0("Params",j,"-Kfold",i)
        models.pred[[nam.ji]]  <-  res
      }
    } # fin i kfold  
    #    print("fin  i kfold")
    error[j]  <-  1-cat2meas(yresp, pred, measure = measure, cost = cost)
  }# fin j param
  #########################################################
  #  print("sale 2 fors")  
  names(error)  <-   paste0("param ",apply(params,1,paste0,collapse="-"))
  imin  <-   which(error==min(error))
  imin2 <- imin[1] #
  if (classif %in% c("classif.np","classif.knn","classif.kernel")){
    #  print("np np np")
    #print(names(par.call))
    if (it) par.call <- par.classif
    par.call$group <- data$df[,response]
    par.call$fdataobj <- data[[vfunc[ifun]]]
    if (name.param[1]=="h")
      par.call[["h"]] <- params[imin2,ifun]
    if (name.param[1]=="knn")
      par.call[["knn"]] <- params[imin2,ifun]
  } else {
    par.call$data <- data
    if (name.param[1]=="pc"){
      par.call$basis.x <- basis.aux}
    if (name.param[1] %in% c("bspline","fourier","constant",
                             "exponential","polygonal","power")){
      basis.x2 <- list()
      for (ifun in 1:nfun){
        basis.x2[[vfunc[ifun]]]  <-  create.fdata.basis(trainData[[vfunc[ifun]]],
                                                      1:params[imin2,ifun],
                                                      type.basis=name.param[ifun])    }
      if (classif %in% c("classif.gsam","classif.glm")) par.call$basis.b <- basis.x2  }
    #  print("name.param");  print(name.param)
    if (name.param[1]=="h"){
      for (ifun in 1:nfun){
        par.call[["par.np"]][[vfunc[ifun]]][[name.param[ifun]]] <- params[imin2,ifun] 
        # print(par.call)
      }
    }
  }
  # fin kfold.aux
  output <- list(par.call=par.call,imin=imin,params=params
               ,error=error,pred=pred)
  
  return(output)
}
########################################################################

############################################################
############################################################
# wrapper of mgcv:::all.vars1 (eleminated from mgcv)

allvars1 <- function (form) 
{
  vars <- all.vars(form)
  vn <- all.names(form)
  vn <- vn[vn %in% c(vars, "$", "[[")]
  if ("[[" %in% vn) 
    stop("can't handle [[ in formula")
  ii <- which(vn %in% "$")
  if (length(ii)) {
    vn1 <- if (ii[1] > 1) 
      vn[1:(ii[1] - 1)]
    go <- TRUE
    k <- 1
    while (go) {
      n <- 2
      while (k < length(ii) && ii[k] == ii[k + 1] - 1) {
        k <- k + 1
        n <- n + 1
      }
      vn1 <- c(vn1, paste(vn[ii[k] + 1:n], collapse = "$"))
      if (k == length(ii)) {
        go <- FALSE
        ind <- if (ii[k] + n < length(vn)) 
          (ii[k] + n + 1):length(vn)
        else rep(0, 0)
      }
      else {
        k <- k + 1
        ind <- if (ii[k - 1] + n < ii[k] - 1) 
          (ii[k - 1] + n + 1):(ii[k] - 1)
        else rep(0, 0)
      }
      vn1 <- c(vn1, vn[ind])
    }
  }
  else vn1 <- vn
  vn1
}

#' @title Functional Classification usign k-fold CV
#'
#' @description Computes Functional Classification using k-fold cross-validation
#' 
#' @param formula   an object of class \code{formula} (or one that can be coerced to that class):
#'  a symbolic description of the model to be fitted. The procedure only considers functional covariates (not implemented for non-functional covariates). 
#  The details of model specification are given under \code{Details}.

#' @param data \code{list}, it contains the variables in the model. 

#' @param classif character,  name of classification method to be used in fitting the model, see \code{Details} section.
# The method to be used in fitting the model. The default method "glm.fit" uses iteratively reweighted least squares (IWLS): the alternative "model.frame" returns the model frame and does no fitting
# indicar que en no recuerdo q caso se usa sin kfold

#' @param par.classif \code{list} of arguments used in the classification method.
#' @param kfold \code{integer}, number of k-fold.
#' @param param.kfold \code{list},  arguments related to number of k-folds for each covariate, see \code{Details} section.

#' @param measure \code{character}, type of measure of accuracy used, see \code{\link{cat2meas}} function.
#' @param cost \code{numeric}, see \code{\link{cat2meas}} function.
#' @param verbose \code{logical}. If \code{TRUE}, print some internal results.
#' @param models \code{logical}. If \code{TRUE}, return a list of the fitted models used,  (k-fold -1) X (number of parameters)
#' @aliases classif.kfold
#' @return 
#' Best fitted model computed by the k-fold CV  using the method indicated 
#' in the \code{classif} argument and also returns:
#'  \enumerate{
#'   \item \code{param.min}, value of parameter (or parameters) selected by k-fold CV. 
#'   \item \code{params.error}, k-fold CV error for each parameter combination. 
#'   \item \code{pred.kfold}, predicted response computed by k-fold CV.
#'   \item \code{model}, if \code{TRUE}, list of models for each parameter combination.
#'   }   

# Best fitted model usign k-fold CV procedure in \code{classif} method.
# @note No implemented for PLS basis yet.
# @references JSS paper

#' @details Parameters for k-fold cross validation:
#' \enumerate{
#'   \item Number of basis elements: 
#' \itemize{
#' \item Data-driven basis such as Functional Principal Componetns (PC). No implemented for PLS basis yet.
#' \item Fixed basis (bspline, fourier, etc.).  
#' }
#'   Option used in some classifiers such as \code{\link{classif.glm}}, \code{\link{classif.gsam}}, \code{\link{classif.svm}}, etc.
#' \item Bandwidth parameter.  Option used in  non-parametric classificiation models such as  \code{\link{classif.np}} and \code{\link{classif.gkam}}.
#' }
#' 
#' @author
#'   Manuel Febrero-Bande, Manuel Oviedo de la Fuente \email{manuel.oviedo@@udc.es}
#'   
#' @keywords  classif
#'   
#' @examples
#' \dontrun{
#' data(tecator)    
#' cutpoint <- 18
#' tecator$y$class <- factor(ifelse(tecator$y$Fat<cutpoint,0,1))
#' table(tecator$y$class )
#' x <- tecator[[1]]
#' x2 <- fdata.deriv(tecator[[1]],2)
#' data  <-   list("df"=tecator$y,x=x,x2=x2)
#' formula  <-   formula(class~x+x2)
#' 
#' # ex: default excution of classifier (no k-fold CV)
#' classif="classif.glm"
#' out.default <- classif.kfold(formula, data, classif = classif)
#' out.default
#' out.default$param.min
#' out.default$params.error
#' summary(out.default)
#' 
#' # ex: Number of PC basis elements selected by 10-fold CV
#' # Logistic classifier
#' kfold = 10
#' param.kfold   <-   list("x"=list("pc"=c(1:8)),"x2"=list("pc"=c(1:8)))
#' out.kfold1   <-   classif.kfold(formula, data, classif = classif,
#'                             kfold = kfold,param.kfold = param.kfold)
#' out.kfold1$param.min
#' min(out.kfold1$params.error)
#' summary(out.kfold1)
#' 
#' # ex: Number of PC basis elements selected by 10-fold CV
#' # Logistic classifier with inverse weighting
#' out.kfold2 <- classif.kfold(formula, data, classif = classif,
#'                             par.classif=list("weights"="inverse"),
#'                             kfold = kfold,param.kfold = param.kfold)
#' out.kfold2$param.min
#' min(out.kfold2$params.error)
#' summary(out.kfold2)
#' 
#' # ex: Number of fourier  basis elements selected by 10-fold CV
#' # Logistic classifier 
#' ibase = seq(5,15,by=2)
#' param.kfold <- list("x"=list("fourier"=ibase),
#'                     "x2"=list("fourier"=ibase))
#' out.kfold3 <- classif.kfold(formula, data, classif = classif,
#'                             kfold = kfold,param.kfold = param.kfold)
#' out.kfold3$param.min
#' min(out.kfold3$params.error)
#' summary(out.kfold3)
#' 
#' # ex: Number of k-nearest neighbors selected by 10-fold CV
#' # non-parametric classifier  (only for a functional covariate)
#' 
#' output <- classif.kfold( class ~ x, data, classif = "classif.knn",
#'                        param.kfold= list("x"=list("knn"=c(3,5,9,13))))
#' output$param.min
#' output$params.error
#' 
#' output <- classif.kfold( class ~ x2, data, classif = "classif.knn",
#'                        param.kfold= list("x2"=list("knn"=c(3,5,9,13))))
#' output$param.min
#' output$params.error 
#' }
#'  
#' @export  classif.kfold
classif.kfold = function(formula, data,  classif="classif.glm",
                          par.classif, kfold = 10,param.kfold=NULL,
                          measure="accuracy",cost,models=FALSE, verbose = FALSE)
  {
  C <- match.call()
  a <- list()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "classif", "par.classif", "kfold", 
                "param.kfold","measure","cost","models","verbose"), names(mf), 0L)
  par.call <- list("formula"=formula)  
  call <- classif
  names.data <- names(data)
  idf  <- which(names.data=="df")
  ntrain <- NROW(data$df)
  ifold <- 1:kfold
  folds <- sample(rep(1:kfold, length.out = ntrain))
  if (classif == "classif.gsam") { 
    tf <- terms.formula(formula, specials = c("s", "te", "t2"))
    #vtab <- rownames(attr(tf,"factors"))
    gp <- interpret.gam(formula)
    terms <- allvars1(gp$fake.formula[-2])
    vnf=intersect(terms,names.data)
    vfunc=intersect(terms,names(data[-idf]))
  }  else {  
    tf <- terms.formula(formula)
    vtab <- rownames(attr(tf,"factors"))
    terms <- attr(tf, "term.labels")
    nt <- length(terms)
    vnf=intersect(terms,names(data$df))
    # vnf2=intersect(vtab[-1],names(data$df)[-1])
    vfunc2=setdiff(terms,vnf)
    vint=setdiff(terms,vtab)
    vfunc=setdiff(vfunc2,vint)
  }
  if (attr(tf, "response") > 0) {
    response <- as.character(attr(tf, "variables")[2])
    pf <- rf <- paste(response, "~", sep = "")  
    }
  nfun <- length(vfunc)
  lista <- list()
  nparams <- numeric(nfun)
  name.param <- character(nfun)
  if (is.null(param.kfold)) default=TRUE
  else {
    default = FALSE
  }
  ensamble <- FALSE
  names(param.kfold)
  #par.boot = list(B=50, N=1000, Nhull=4, Nnbh=9)
  #if (!is.null(param.kfold$B)){
  if (classif %in% c("classif.bootstrap","classif.adaboost")){
#print("entra ensamble")
    ensamble <- TRUE
    name.param <- names(param.kfold)
# print(name.param)    
    BB <- param.kfold#[[name.param]]

    for (i in 1:nfun) {
      #print(i)
      # param.kfold[[vfunc[i]]][[name.param]]  <-  BB
      lista[[vfunc[1]]]  <-  BB
      nparams[i]  <-  length(lista[[vfunc[i]]])
    }
#print(nparams);    print(lista)
  } 
  else{
  #if (!ensamble){
  for (i in 1:nfun){
    #if (is.null(param.kfold[[vfunc[i]]])) param.kfold[[vfunc[i]]] <- list("pc"=1:3)
    #if (param.kfold[[vfunc[i]]]=="default") {
    if (default) {
      param.kfold[[vfunc[i]]]  <-  "default"
      names(param.kfold) <- vfunc[i]
      name.param[i] <- "default"
    }     else name.param[i]  <-  names(param.kfold[[vfunc[i]]])
    lista[[vfunc[i]]] <- param.kfold[[vfunc[i]]][[1]]
    nparams[i]  <- length(lista[[vfunc[i]]])
  }
  }
  params  <- do.call("expand.grid",lista)
  nexpand <- nrow(params)
   #for (ifun in 1:nfun){
   # name.param[ifun]  <-  names(param.kfold[[vfunc[ifun]]])
  #}
  max.params <- apply(params,2,max)
  error <- rep(NA,nexpand)
  if (missing(par.classif)) par.classif=list()
  yresp <- data$df[,response]
  lev <- levels(yresp)
  if (missing(cost)) cost = rep(1, nlevels(yresp))

  pred <- factor(rep(NA,len=ntrain),levels=lev)
  #pred.df  <- data.frame(pred)
  #for (j in 2:nparams) pred.df <- cbind(pred.df,pred)
  #pred.df <- array(NA, dim=c( ntrain, nparams))
  if (!ensamble){
  basis.aux <- list()
  for (ifun in 1:nfun){
  ################################ PC ##########################
   if (name.param[ifun] %in% c("pc")){ 
    basis.aux[[vfunc[ifun]]]  <-  create.fdata.basis(data[[vfunc[ifun]]],
                                1:max.params[ifun],
                                type.basis=name.param[ifun])
    }
   ################################ FiXED basis ##########################
   if (name.param[1] %in% c("bspline","fourier","constant",
                           "exponential","polygonal","power")){
    basis.aux[[vfunc[ifun]]]  <-  create.fdata.basis(data[[vfunc[ifun]]],
                                1:max.params[ifun],
                                type.basis=name.param[ifun])
   }
   ################################ bandwidth ##########################
    if (name.param[1]=="h"){
        #par.call[["par.np"]][[vfunc[ifun]]][[name.param[ifun]]] <- NULL
      par.classif$par.np <- NULL # par.np
      }
  }
  }
  it <- TRUE
  if (models) models.pred <- list()
  
  ############################  
  for (j in 1:nexpand){ 
    if (verbose) {
      print(params[j,])
    cat("Param j,",j,"\n")
}
    for (i in ifold){ 
      if (verbose) cat("j,",j,"params",as.numeric(params[j,])," kfold i,",i,"\n")
        #Segment your data by fold using the which() function 
      testIndexes <- which(folds==i,arr.ind=TRUE)
      trainData <- subset.ldata(data,folds!=i)
      ipred <- folds==i
      testData  <-   subset.ldata(data,ipred)
      if (classif %in% c("classif.np","classif.knn","classif.kernel")){
#print("entra classif.np")        
        if (it) par.call <- par.classif
        par.call$group <- trainData$df[,response]
        par.call$fdataobj <- trainData[[vfunc[ifun]]]
        it <- FALSE
        if (name.param[1]=="h")          par.call[["h"]] <- params[j,]
        if (name.param[1]=="knn")        par.call[["knn"]] <- params[j,]
        #if (!is.null(par.call$par.S$w)) 
          #par.call$par.S$w  <- par.classif$par.S$w[folds!=i]

        res <- do.call(call,par.call)
        aux <- predict.classif(res,testData[[vfunc[ifun]]])
        pred[ipred] <- aux
        #if (!is.null(par.call$par.S$w)) 
         # par.call$par.S$w  <- par.classif$par.S$w
        #print(pred)
      } else{
      par.call$data <- trainData
    #  if (i==1) par.call <- c(par.call, par.classif)
    #  cat("j,",j,"params",as.numeric(params[j,])," kfold i,",i,"\n")
      #print(names(par.call));            print(names(par.classif))
      #print(it)
      if (it) par.call <- c(par.call, par.classif)
      it <- FALSE
      #print(names(par.call$par.np))
             
      ################################ PC ##########################
      if (name.param[1]=="pc"){ #atajo para que tarde menos
        basis.x <- basis.aux
        for (ifun in 1:nfun){
          basis.x[[vfunc[ifun]]]$l <- basis.aux[[vfunc[ifun]]]$l[1:params[j,ifun]]
          basis.x[[vfunc[ifun]]]$basis <- basis.aux[[vfunc[ifun]]]$basis[1:params[j,ifun]]
          basis.x[[vfunc[ifun]]]$x <- basis.aux[[vfunc[ifun]]]$x[folds!=i,1:params[j,ifun],drop=F]
        }
        par.call$basis.x <- basis.x
      } 
      ################################ FiXED basis ##########################
      if (name.param[1] %in% c("bspline","fourier","constant",
                            "exponential","polygonal","power")){
      basis.x2 <- list()
      for (ifun in 1:nfun){
        basis.x2[[vfunc[ifun]]]  <-  create.fdata.basis(trainData[[vfunc[ifun]]],
                                                      1:params[j,ifun]
                                                      ,type.basis=name.param[ifun])
      }
      if (classif %in% c("classif.gsam","classif.glm")) par.call$basis.b <- basis.x2
      }
      ################################ BANDWIDTH ##########################
      if (name.param[1]=="h"){
        for (ifun0 in 1:nfun){
         #par.call[["par.np"]][[vfunc[ifun0]]][[name.param[ifun0]]] <- params[j,ifun0]
          par.call[["par.np"]][[vfunc[ifun0]]] <- list("h"=params[j,ifun0])
        }
      }
      par.call$data <- trainData
      
      ########## if (name.param[1]=="B") par.call$B <- param.kfold$B[j]
      # par.boot = list(B=50, N=1000, Nhull=4, Nnbh=9) 
      if (ensamble) {
        if (classif=="classif.bootstrap"){
          #par.call$par.boot <- lista[[j]
          par.call$par.boot<- as.list(params[j,]) # para bootstrap pero no para adaboost
          
          #print("ensamble2")
          names(par.call$par.boot)<-name.param
        } else{
          par.call[["B"]]<-params[j,]
        }
        
      }
      res <- do.call(call,par.call)
      pred[ipred]  <-  predict.classif(res,testData)
      } 
      if (models) {
          nam.ji <- paste0("Params",j,"-Kfold",i)
          models.pred[[nam.ji]]  <-  res
      }
  } # fin i kfold  
    #error[j]  <-  1-cat2meas(yresp, pred, measure = measure)
    error[j]  <-  1-cat2meas(yresp, pred, measure = measure, cost = cost)
    
 # print(error)
  }# fin j param
  #########################################################
#  print("sale 2 fors")  
  names(error)  <-   paste0("param ",apply(params,1,paste0,collapse="-"))
  #imin  <-   which.min(error)
  imin  <-   which(error==min(error))
  imin2 <- imin[1] #
  
  if (classif %in% c("classif.np","classif.knn","classif.kernel")){
  #  print("np np np")
    #print(names(par.call))
    if (it) par.call <- par.classif
    par.call$group <- data$df[,response]
    par.call$fdataobj <- data[[vfunc[ifun]]]
    if (name.param[1]=="h")
      par.call[["h"]] <- params[imin2,ifun]
    if (name.param[1]=="knn")
      par.call[["knn"]] <- params[imin2,ifun]
  } else {
  par.call$data <- data
  if (name.param[1]=="pc"){
    par.call$basis.x <- basis.aux}
  if (name.param[1] %in% c("bspline","fourier","constant",
                           "exponential","polygonal","power")){
    basis.x2 <- list()
    for (ifun in 1:nfun){
      basis.x2[[vfunc[ifun]]]  <-  create.fdata.basis(trainData[[vfunc[ifun]]],
                                1:params[imin2,ifun],
                                type.basis=name.param[ifun])    }
    if (classif %in% c("classif.gsam","classif.glm")) par.call$basis.b <- basis.x2  }
#  print("name.param");  print(name.param)
  if (name.param[1]=="h"){
    for (ifun in 1:nfun){
     par.call[["par.np"]][[vfunc[ifun]]][[name.param[ifun]]] <- params[imin2,ifun] 
    # print(par.call)
     }
    }
  }
  res <- do.call(call, par.call)
  res$param.min <- params[imin,]
  res$params.error = error
  res$C <- res$C[1:2]
  res$pred.kfold <- pred
  if (models){
    res$models = models.pred
    res$params = params
  }
  return(res)
} 
###########################################################
