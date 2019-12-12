#' @title Variable Selection using Functional Additive Models      
#' 
#' @description  Computes functional GAM model between functional covariates
#'  \eqn{(X^1(t_1),\cdots,X^{q}(t_q))}{(X(t_1),...,X(t_q))} and non functional covariates
#'   \eqn{(Z^1,...,Z^p)}{(Z1,...,Zp)} with a scalar response \eqn{Y}.
#' @details This function is an extension of the functional generalized spectral additive 
#' regression models: \code{\link{fregre.gsam}} where the \eqn{E[Y|X,Z]} is related to the 
#' linear prediction \eqn{\eta} via a link function \eqn{g(\cdot)}{g(.)} with integrated 
#' smoothness estimation by the smooth functions \eqn{f(\cdot)}{f(.)}. 
#' \deqn{E[Y|X,Z])=\eta=g^{-1}(\alpha+\sum_{i=1}^{p}f_{i}(Z^{i})+\sum_{k=1}^{q}\sum_{j=1}^{k_q}{f_{j}^{k}(\xi_j^k)})}{E[Y|X,Z]=\eta=g^{-1}(\alpha+\sum_i  f_i(Z_{i})+\sum_k^q\sum_{j=1}^{k_q}{f_j^k(\xi_j^k)})}
#' where \eqn{\xi_j^k}{\xi_j^k} is the coefficient of the basis  function expansion of 
#' \eqn{X^k}, (in PCA analysis \eqn{\xi_j^k}{\xi_j^k} is the score of the \eqn{j}-functional
#' PC of \eqn{X^k}.
#'  
#' The smooth functions \eqn{f(\cdot)}{f(.)} can be added to the right hand side of the formula
#' to specify that the linear predictor depends on smooth functions of predictors using smooth 
#' terms \code{\link{s}} and \code{\link{te}} as in  \code{\link{gam}} (or linear functionals of 
#' these as \eqn{Z\beta} and \eqn{\big<X(t),\beta\big>}{< X(t),\beta(t) >} in \code{\link{fregre.glm}}). 
                                                                                                
#' @param data List that containing the variables in the model. 
#' "df" element is a data.frame containing the response and scalar covariates 
#' (numeric and factors variables are allowed). Functional covariates of class
#'   \code{fdata} or \code{fd} are included as named components in the \code{data} list. 
#' @param y Caracter string with the name of the scalar response variable.
# @param x Caracter string vector with the name of the scalar and functional 
# potential covariates. If missing, all scalar and functional covariates are included.                     
#' @param include vector with the name of variables to use. By default \code{"all"}, all variables are used.
#' @param exclude vector with the name of variables to not use. By default  \code{"none"}, no variable is deleted.
#' @param family a description of the error distribution and link function to
#' be used in the model. This can be a character string naming a family
#' function, a family function or the result of a call to a family function.
#' (See \code{\link{family}} for details of family functions.)
#' @param weights weights 
#' @param alpha Alpha value for testing the independence among covariate X and residual
#'  e in previous steps. By default is \code{0.05}.
#' @param basis.x  List of basis for functional covariates, see same argument 
#' in \code{\link{fregre.gsam}}. By default, the function uses a basis of 4 PC
#' to represent the functional covariate. 
# \itemize{
# \item{\code{list}} (recomended)  List of basis for functional covariates,
#  see same argument in \code{\link{fregre.gsam}}.
# \item{\code{vector}} (by default) Vector with three parameters:
# \enumerate{
# \item Type of basis. By default \code{basis.x[1]="pc"}, principal
# component basis is used  for each functional covariate included in the model.
# Other options \code{"pls"} and \code{"bspline"}.  
# \item Maximum number of components"ncomp" (or basis elements "nbasis") to be used.
#  By default, \code{basis.x[2]=4}. 
# \item Selection of number of basis components/elements: 
# if \code{basis.x[3]=TRUE}  (by default), the algorithm selects the
#  best components from the first \code{basis.x[2]}. If \code{basis.x[3]=FALSE}, the
#  model is estimated using \code{ncomp} components for each functional covariate.}}
#' @param kbs The dimension of the basis used to represent the smooth term. The default 
#' depends on the number of variables that the smooth is a function of.
#' @param dcor.min Threshold for a variable to be entered into the model. X is discarded 
#' if the distance correlation \eqn{R(X,e)< dcor.min} (e is the residual of previous steps).
#' @param par.model Model parameters.
#' @param xydist List with the inner distance matrices of each variable (all potential 
#' covariates and the response).
#' @param trace Interactive Tracing and Debugging of Call.
# @param CV TRUE, Cross-validation (CV) is done.
# @param smooth If TRUE, a smooth estimate is made for all covariates included in the model
#  (less for factors). The model is adjusted with the estimated variable linearly or smoothly. 
#  If the models are equivalent, the model is adjusted with the linearly estimated variable.
#' 
#' @return Return an object corresponding to the estimated additive mdoel using 
#' the selected variables (ame output as the\code{\link{fregre.gsam}} function) and the following elements:
#' \itemize{
#' \item{\code{gof}}, the goodness of fit for each step of VS algorithm.
#' \item{\code{i.predictor}}, \code{vector} with 1 if the variable is selected, 0 otherwise.
#' \item{\code{ipredictor}}, \code{vector} with the name of selected variables (in order of selection)
#' \item{\code{dcor}}, the value of distance correlation for each potential covariate and the residual of the model in each step.
#' }
#' 
#' @note If the formula only contains a non functional explanatory variables (multivariate covariates),
#'  the function compute a standard  \code{\link{gam}} procedure.
#' 
#' @author Manuel Feb-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@usc.es}
#' 
#' @seealso See Also as:  \code{\link{predict.fregre.gsam}} and \code{\link{summary.gam}}.
#' Alternative methods: \code{\link{fregre.glm}}, \code{\link{fregre.gsam}}
#'  and \code{\link{fregre.gkam}}.
#'  
#' @references Febrero-Bande, M., Gonz\'alez-Manteiga, W. and Oviedo de la
#' Fuente, M. Variable selection in functional additive regression models,
#' (2018).  Computational Statistics, 1-19. DOI:
#' \url{https://doi.org/10.1007/s00180-018-0844-5}
#' 
#' @keywords regression
#' @examples 
#' \dontrun{ 
#' data(tecator)
#' x=tecator$absorp.fdata
#' x1 <- fdata.deriv(x)
#' x2 <- fdata.deriv(x,nderiv=2)
#' y=tecator$y$Fat
#' xcat0 <- cut(rnorm(length(y)),4)
#' xcat1 <- cut(tecator$y$Protein,4)
#' xcat2 <- cut(tecator$y$Water,4)
#' ind <- 1:129
#' dat    <- data.frame("Fat"=y, x1$data, xcat1, xcat2)
#' ldat <- list("df"=dat[ind,],"x"=x[ind,],"x1"=x1[ind,],"x2"=x2[ind,])
#' # 3 functionals (x,x1,x2), 3 factors (xcat0, xcat1, xcat2)
#' # and 100 scalars (impact poitns of x1) 
#' 
#' # Time consuming
#' res.gam1 <- fregre.gsam.vs(data=ldat,y="Fat") # All the covariates
#' summary(res.gam1)
#' res.gam1$ipredictors
#' 
#' covar <- c("xcat0","xcat1","xcat2","x","x1","x2")
#' res.gam2 <- fregre.gsam.vs(data=ldat, y="Fat", include=covar)
#' summary(res.gam2)
#' res.gam2$ipredictors 
#' 
#' # Prediction like fregre.gsam() 
#' newldat <- list("df"=dat[-ind,],"x"=x[-ind,],"x1"=x1[-ind,],
#'                 "x2"=x2[-ind,])
#' pred.gam1 <- predict(res.gam1,newldat)
#' pred.gam2 <- predict(res.gam2,newldat)
#' plot(dat[-ind,"Fat"],pred.gam1)
#' points(dat[-ind,"Fat"],pred.gam2,col=2)
#' }
#' 
#' @export
fregre.gsam.vs  <- function(data = list(), y, 
                          include = "all", exclude = "none"
                          ,family = gaussian(), weights = NULL 
                          , basis.x = NULL 
                          #,CV = TRUE,
                          ,kbs , dcor.min = 0.1, alpha = 0.05
                          , par.model, xydist, trace = FALSE
                          #,smooth = TRUE
                          ){
  #0-
  verbose = FALSE
  #criterio = "sp"
  if (missing(y)) 
    stop("The name of the response must be specified in the 'y' argument")  
  # if (missing(alpha)) alpha <- .05
  if (missing(par.model)) 
    par.model <- list() #se puede cambiar la familia
  n.edf.correction <- FALSE
  namdata <- names(data)
  idf <- which("df"==namdata )
  ydat <- data$df[,y]
  namfunc <- names(data[-idf])
  namnfunc <- setdiff(names(data$df),y)
  namnfunc <- setdiff(namnfunc,exclude)
  namfunc <- setdiff(namfunc,exclude)
  # if (missing(x)) {  
  #   # print("misssing X")
  #   x <- c(namnfunc,namfunc)
  #   ifunc <- namfunc
  #   infunc <- namnfunc
  # }
  if (include[1] == "all"){
    x <- c(namnfunc,namfunc)
    ifunc <- namfunc
    infunc <- namnfunc
  }  else {
    x = include
    ifunc <- intersect(x,namfunc)
    infunc <- intersect(c(y,x),namnfunc)
    #    xdatos <- as.list(data$df[,infunc,drop=F])
    #    xdatos <- c(xdatos,data[ifunc])
  } 
  # var.name <- names(x)
  # a <- sapply(x, class, USE.NAMES = T)
  # x <- x[, a != "character", drop = F]
  # var.name <- names(x)
  # if (include[1] != "all") {
  #   var.name <- intersect(include, var.name)
  # }
  # if (exclude[1] != "none") {
  #   var.name <- setdiff(var.name, exclude)
  # }
  xdatos <- as.list(data$df[,y,drop=F])
  xdatos <- c(xdatos,as.list(data$df[,infunc,drop=F]),data[ifunc])
  ldata0 <- xdatos
  resp <- y
  xentran <- NULL
  tipoentran <- NULL
  #1- caclulo= distancias de cada objeto consigo mismo
  #  print("1 calculando distancias")
  #  ldist0 <- dist.list(ldata0)
  xynam <- names(ldata0)
  nvar <- length(x)
  if (missing(xydist))    {
    xydist <- ldist0 <- dist.list(ldata0)
  }    else {
    ldist0 <- xydist[xynam]
  }
  #  print("2 calculando distancias")
  parar <- FALSE
  it <- 1
  basisx <- NULL
  npredictors <- length(ldata0)
  ipredictors <- numeric(npredictors-1)
  names(ipredictors) <- setdiff(names(ldata0),resp)
  nvar <- npredictors - 1
  if (missing(kbs)) 
    kbs <- - 1 #rep(-1,nvar)
  #names(kbs) <- names(ipredictors)
  basis.list <- FALSE
  if (is.null(basis.x)) {
    #ncomp <- rep(4,length(ifunc))
    #names(ncomp) <- ifunc
    ncomp <- 4
    type.basis = "pc"
    ncomp.fix = TRUE
  }  else {
    if (is.list(basis.x)){ # base proporcionada por el usuario
      ncomp.fix = TRUE
      basis.list = TRUE
    }else{
      if (length(basis.x)!=3) stop("basis.x must be a vector with length three")
      type.basis = basis.x[1]
      #ncomp <- rep(as.numeric(basis.x[2]),length(ifunc))
      ncomp <- as.numeric(basis.x[2])
      ncomp.fix <- as.logical(basis.x[3])
    }
  }  
  #type.basis = "pc",  ncomp,ncomp.fix = FALSE
    
  dcor=matrix(0,nrow= nvar,ncol= nvar)
  colnames(dcor)=c(names(ipredictors))
  rownames(dcor)=1:nvar
  #2- caclulo correlacion  de cada la distancia de la respuesta vs distancia del resto de objetos
  # print("2 calculando correlaciones")
  dist_resp <- dcor.y(ldist0,resp)
  dcor[it,names(dist_resp)] <- dist_resp*(dist_resp > dcor.min)
  n.edf <- length(ydat)
  fpredictors.nl <- ""
  form.nl <- paste(resp,"~",sep="")
  basis2 <- list()  
  ycen <- data$df[,y]-mean(data$df[,y])
  gof <- NULL
  anyfdata <- FALSE
  res.prev <- gam(as.formula(paste0(form.nl,1)),data=data$df)
  while (!parar){
    # print("3 Seleccion variable-Regresion")
    #      print("bucle")    
    esfactor <- FALSE #SE UTILIZA PARA NO PONER s(factor)
    #3- seleccion dcor mas elevada
    nam <-  names(dist_resp)
    ind.xentra <- which.max(dist_resp)
    xentra <- nam[ind.xentra]
    #dd <- dcor.ttest(ldist0[[resp]],ldist0[[xentra]],distance=TRUE)
    dd <- dcor.test(ldist0[[resp]],ldist0[[xentra]],n=n.edf)
    #if (trace){      print(dd);      print(xentra)    }   
    if (verbose) {
      print(it);      print(n.edf);      print(nam)
      cat("Enter the variable ",xentra);      print(dd);      print(names(dd))
    }
    if (is.null(basisx)) {
      basisb <- basisx <- list()
    }
    if (dd$p.value > alpha) {
      parar=TRUE
      if (trace) print("The algorithm ends because no variable is significant")
    }    else{
      rownames(dcor)[it] <- xentra
      if (trace)      cat("Covariate: ")
      if (trace)       print(xentra)
      if (is.fdata(ldata0[[xentra]]) & !basis.list) {
        par.basis.ipred <- list()
        anyfdata <- TRUE
        par.basis.ipred$fdataobj <- ldata0[[xentra]]
        if (xentra %in% names(type.basis)){
          itype <- which(xentra== names(type.basis))
          type.basis.ipred <- type.basis[itype]
        }    else {type.basis.ipred <- type.basis[1]}
        tbasis <- type.basis.ipred
        if (type.basis.ipred=="fourier") tbasis <- "basis"
        if (type.basis.ipred=="bspline") tbasis <- "basis"
        nam <- "fregre.gsam.cv"
        #  par.basis.ipred <- list()
        par.basis.ipred$y <- resp#entra la etiqueta y no toda la variable como en el basis.cv o pc.cv
        par.basis.ipred$x <- xentra
        par.basis.ipred$data <- data
        if (!ncomp.fix)      {
          #print("entra gsam.CV")
          res <-  fregre.gsam.cv(data,resp,xentra, alpha = alpha
                               #,type.basis=type.basis, kbs = kbs
                               ,type.basis=tbasis, kbs = kbs
                               ,ncomp=ncomp,ncomp.fix=ncomp.fix)  
#print(summary(res))
#print(res$ncomp.opt)

          if (tbasis!="basis"){
           res$basis.x[[xentra]]$basis <- res$basis.x[[xentra]]$basis[res$ncomp.opt,]
          }
          res$basis.x[[xentra]]$l <- res$ncomp.opt
          basisx[[xentra]] <- res$basis.x[[xentra]]
          basisb[[xentra]] <- res$basis.b[[xentra]]
          if (trace)   print("results of internal function fregre.gsam.cv")
          if (trace)   print(summary(res))
          # parar=TRUE
        }        else{
         # print("NO entra gsam.CV")
          switch(tbasis,
                 "pc"={
                   best.pc <- 1:ncomp#[xentra]
                   basis1 <- create.pc.basis(ldata0[[xentra]],best.pc) #es lo mismo ?
                 },
                 "pls"={
                   best.pc <- 1:ncomp#[xentra]
                   basis1 <- create.pls.basis(ldata0[[xentra]],ldata0[[resp]],best.pc)
                 },"basis"={
                   #            best.pc <- 1:kmax
                   basis1 <- create.bspline.basis(ldata0[[xentra]]$rangeval,nbasis=ncomp)
                   basis2 <- basis1
                 })
          basisx[[xentra]] <- basis1
          basisb[[xentra]] <- basis2
        }
      }       
      if (!parar) {
        xentran <- c(xentran,xentra)
        ind.xentra2 <- which(xentra==names(ldist0))
        ldist0 <- ldist0[-ind.xentra2]
        ipredictors[xentra] <- ipredictors[xentra]+1
        #4- contruccion del modelo para esta variable #consido un catalogo de 4 posibilidades (lineal/nolineal, funcional/scalar)
        if (is.factor(ldata0[[xentra]])) esfactor=TRUE                 
        fx=FALSE
        if (esfactor)      {
          fpredictors.nl <- xentra  
        }
        else         {   
          fpredictors.nl <- paste("s(",xentra,",k=",kbs,",","fx=",fx,")",
                                sep="",collapse="+")
        }
        nam.model <- "fregre.gsam"          
        form.nl <- (paste(form.nl,"+",fpredictors.nl,sep=""))  
        par.model$formula <- as.formula(form.nl)
        par.model$data <- data
        par.model[["basis.x"]] <- basisx
        par.model[["basis.b"]] <- basisb
        res.nl  <- do.call(nam.model,par.model)    
        mejora <- res.prev$gcv.ubre > res.nl$gcv.ubre
        if (it>1 & mejora){
          aov <- anova(res.nl,res.prev,test="F")
          mejora <- aov[["Pr(>F)"]][2] < alpha
          if (is.na(mejora)) mejora <- TRUE
        } 
        if (mejora){
        suma <- summary(res.nl)          
        dd <- res.nl$dims
        df <- dd[["p"]]
        edf <- dd[["N"]] - dd[["p"]]    
        #         sr2 <- sum(res.nl$residuals^2)/edf
        r2 <- 1 - sum(res.nl$residuals^2)/sum(ycen^2)
        if (n.edf.correction) n.edf <-  edf
        ldata0[[resp]] <- res.nl$residuals
        ldist0[[resp]] <- as.matrix(dist(ldata0[[resp]]), diag =TRUE, upper = TRUE,     p = 2)
        #2- caclulo correlacion  de cada la distancia de la respuesta vs distancia del resto de objetos
        #print(it)
        #print(npredictors)
        if (it==(npredictors-1)) {
          #print("iiiiittttttttt")          
          parar=TRUE
          #gof <- rbindgof,c(suma$logLik,suma$BIC,suma$AIC,edf,r2))
          gof <- rbind(gof,c(AIC(res.nl),deviance(res.nl),res.nl$df.residual,suma$r.sq,suma$dev.expl,res.nl$gcv.ubre))
        }        else{
          it <- it+1
          #print("3 calculando correlaciones")
          dist_resp <- dcor.y(ldist0,resp)
          dcor[it,names(dist_resp)]=dist_resp*(dist_resp > dcor.min)
        }
        #     print(summary(res.nl))
        if (!parar){ 
          #            gof <- rbind(gof,drop(c(suma$logLik,suma$BIC,suma$AIC,edf,r2)))
          gof <- rbind(gof,c(AIC(res.nl),deviance(res.nl),
                           res.nl$df.residual,suma$r.sq,
                           suma$dev.expl,res.nl$gcv.ubre))
        }
        else { 
          cat("The variable ",xentra,"does not improve the fit of the model\n")
        }
        res.prev <- res.nl
  }   else{
  if (it==nvar) parar=TRUE
  res.nl <- res.prev
  # ipredictors <- ipredictors[-length(ipredictors)]
  xentran <- xentran[-length(xentran)]
  }
  }
  }  }
  # print(dim(gof));print(gof)
  gof <- data.frame(xentran,(gof))
  if (is.null(xentran)) {
    warning("No variable selected, a null model is estimated")    
    res.nl <- fregre.gsam(as.formula(paste(resp,"~1",sep="")),data=data)
    suma <- summary(res.nl)
    gof <- data.frame(rbind(c(1,AIC(res.nl),deviance(res.nl),
                            res.nl$df.residual,suma$r.sq,
                            suma$dev.expl,res.nl$gcv.ubre)))      
  }  
  names(gof) <- c("xentra","AIC","deviance","df.residual","r.sq","dev.expl","GCV.ubre")
  res.nl$gof <- gof
  res.nl$i.predictor = ipredictors
  res.nl$ipredictors = xentran
  res.nl$xydist = xydist
  res.nl$dcor = dcor[1:(it-1),]
  return(res.nl)
}
########################################################
# dist.list <- fda.usc:::dist.list
# dcor.y <- fda.usc:::dcor.y
# sp <- fda.usc:::sp
# pvalue.anova <- fda.usc:::pvalue.anova
########################################################
# res.gam2 <- fregre.gsam.vs(data=ldat, y="Fat", include = covar)
# summary(res.gam2)
# names(ldat)
# bsp<-create.bspline.basis(c(0,1),7)
# basis<-list("x"=bsp,"x1"=bsp,"x2"=bsp)
# res.gam2 <- fregre.gsam.vs(data=ldat, y="Fat",basis.x=basis)
# summary(res.gam2)
# 
# basis2<-c("pc",12,TRUE)
# res.gam2 <- fregre.gsam.vs(data=ldat, y="Fat",basis.x=basis2,exclude="x")
# summary(res.gam2)
# 
# basis2<-c("pc",12,FALSE)
# res.gam2 <- fregre.gsam.vs(data=ldat, y="Fat",basis.x=basis2)
# summary(res.gam2)
# basis2<-c("pc",11,FALSE)
# res.gam2 <- fregre.gsam.vs(data=ldat, y="Fat",basis.x=basis2)
# summary(res.gam2)
# 
# basis2<-c("pc",7,TRUE)
# res.gam2 <- fregre.gsam.vs(data=ldat, y="Fat",basis.x=basis2)
# summary(res.gam2)
# 
# res.gam2 <- fregre.gsam.vs(data=ldat, y="Fat")
# 
# basis2<-c("bspline",13,FALSE)
# res.gam2 <- fregre.gsam.vs(data=ldat, y="Fat",basis.x=basis2)
# summary(res.gam2)
# 
# 
# basis2<-c("bspline",7,TRUE)
# res.gam2 <- fregre.gsam.vs(data=ldat, y="Fat",basis.x=basis2)
# summary(res.gam2)
# 
# basis2<-c("pls",5,FALSE)
# res.gam2 <- fregre.gsam.vs(data=ldat, y="Fat",basis.x=basis2)
# summary(res.gam2)
# 
# basis2<-c("pls",4,TRUE)
# res.gam1 <- fregre.gsam.vs(data=ldat, y="Fat",basis.x=basis2)
# summary(res.gam2)
# 
# basis2<-c("pls",8,FALSE)
# res.gam2 <- fregre.gsam.vs(data=ldat, y="Fat",basis.x=basis2)
# summary(res.gam2)
# anova(res.gam1,res.gam2)
