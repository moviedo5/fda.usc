predict.gls<-function (object, newdata, na.action = na.fail, ...)
{
  if (missing(newdata)) {
    return(fitted(object))
  }
  if (inherits(object$formula, "formula"))    form <- getCovariateFormula(object$formula)  
  else   form <- getCovariateFormula(object)  
  mfArgs <- list(formula = form, data = newdata, na.action = na.action)
  mfArgs$drop.unused.levels <- TRUE
  dataMod <- do.call("model.frame", mfArgs)
  contr <- object$contrasts
  for (i in names(dataMod)) {
    if (inherits(dataMod[, i], "factor") && !is.null(contr[[i]])) {
      levs <- levels(dataMod[, i])
      levsC <- dimnames(contr[[i]])[[1]]
      if (any(wch <- is.na(match(levs, levsC)))) {
        stop(sprintf(ngettext(sum(wch), "level %s not allowed for %s",
                              "levels %s not allowed for %s"), paste(levs[wch],
                                                                     collapse = ",")), domain = NA)
      }
      attr(dataMod[, i], "contrasts") <- contr[[i]][levs,
                                                    , drop = FALSE]
    }
  }
  N <- nrow(dataMod)
  if (length(all.vars(form)) > 0) {
    X <- model.matrix(form, dataMod)
  }
  else {
    X <- array(1, c(N, 1), list(row.names(dataMod), "(Intercept)"))
  }
  cf <- coef(object)
  val <- c(X[, names(cf), drop = FALSE] %*% cf)
  attr(val, "label") <- "Predicted values"
  if (!is.null(aux <- attr(object, "units")$y)) {
    attr(val, "label") <- paste(attr(val, "label"), aux)
  }
  val
}
################



##predictions not accounting for correlation structure - using Delta method
##gls
predictSE.gls <- function(mod, newdata, se.fit = TRUE, ...){
  
  ##first part of code converts data.frame (including factors) into design matrix of model
  ##fixed <- mod$call$model[-2] #extract only fixed portion of model formula
  fixed <- mod$formula[-2] #modification suggested by C. R. Andersen to extract left part of model formula
  tt <- terms.formula(mod$formula)
  TT <- delete.response(tt)
  newdata <- as.data.frame(newdata)
  
  #################################################################################################################
  ########################### This following piece of code is modified from predict.lme( ) from nlme package
  #################################################################################################################  
  mfArgs <- list(formula = fixed, data = newdata)
  dataMix <- do.call("model.frame", mfArgs)
  
  ## making sure factor levels are the same as in contrasts
  contr <- mod$contrasts
  for(i in names(dataMix)) {
    if (inherits(dataMix[,i], "factor") && !is.null(contr[[i]])) {
      levs <- levels(dataMix[,i])
      levsC <- dimnames(contr[[i]])[[1]]
      if (any(wch <- is.na(match(levs, levsC)))) {
        stop(paste("Levels", paste(levs[wch], collapse = ","),
                   "not allowed for", i))
      }
      attr(dataMix[,i], "contrasts") <- contr[[i]][levs, , drop = FALSE]
    }
  }
  #################################################################################################################
  ########################### The previous piece of code is modified from predict.lme( ) from nlme package
  #################################################################################################################
  
  m <- model.frame(TT, data=dataMix)
  des.matrix <- model.matrix(TT, m)
     # newdata <- des.matrix  #we now have a design matrix 
  ######START OF PREDICT FUNCTION
  ######
  fix.coef <- coef(mod)
  ncoefs <- length(fix.coef)
  names.coef <- labels(fix.coef)
  nvals <- dim(newdata)[1]
  
  ##check for intercept fixed effect term in model
  int.yes <- any(names.coef == "(Intercept)")
  
  ##if no intercept term, return error
  if(!int.yes) stop("\nThis function does not work with models excluding the intercept terms\n")
  
  formula <- character(length=ncoefs)
  
  
  nbetas <- ncoefs - 1
  
  if(int.yes & nbetas >= 1) {
    ##create loop to construct formula for derivative
    formula <- paste("Beta", 1:nbetas, sep="")
    formula <- c("Beta0", formula)
  } else {
    if(int.yes & nbetas == 0) {
      formula <- "Beta0"
    }
  }
  ##for models without intercept - formula <- paste("Beta", 1:ncoefs, sep="")
  
  
  ##a loop to assemble formula
  ##first, identify interaction terms
  inters <- rep(NA, ncoefs)
  for (m in 1:ncoefs) {
    inters[m] <- attr(regexpr(pattern = ":", text = names.coef[m]), "match.length")
  }
  
  ##change the name of the labels for flexibility
  names.cov <- paste("cov", 1:ncoefs-1, sep="")
  
  if(!int.yes) {names.cov <- paste("cov", 1:ncoefs, sep="")}
  
  id <- which(inters == 1)
  for (k in 1:length(id)) {
    names.cov[id[k]] <- paste("inter", k, sep="")
  }
  
  ##iterate and combine betas and covariates
  formula2 <- character(length = ncoefs)
  for(b in 1:ncoefs) {
    formula2[b] <- paste(formula[b], names.cov[b], sep="*")
  }
  ##replace with Beta0 if fixed intercept term present
  if(int.yes) {formula2[1] <- "Beta0"}
  
  ##collapse into a single equation and convert to expression
  ##parse returns the unevaluated expression
  eq.space <- parse(text  = as.expression(paste(formula2, collapse="+")),
                    srcfile = NULL)
  ##add step to remove white space to avoid reaching 500 character limit
  
  ##remove space within expression
  no.space <- gsub("[[:space:]]+", "", as.character(eq.space))
  equation <- parse(text = as.expression(no.space))
  
  
  
  ##
  
  if(identical(se.fit, TRUE)) {
    ##determine number of partial derivatives to compute
    part.devs <- list( )
    for(j in 1:ncoefs) {
      part.devs[[j]] <- D(equation, formula[j])
    }
    
  }
  
  ##determine number of covariates excluding interaction terms
  ncovs <- ncoefs - length(id)
  
  ##assign values of covariates
  cov.values <- list()
  
  ##if only intercept, then add column
  if(int.yes && ncovs == 1) {
    cov.values[[1]] <- 1
  }
  
  if(int.yes && ncovs > 1) {
    cov.values[[1]] <- rep(1, nvals)
    for (q in 2:ncoefs) {
      cov.values[[q]] <- newdata[, labels(fix.coef)[q]]
    }
  } else {
    for (q in 1:ncoefs) {
      cov.values[[q]] <- newdata[, labels(fix.coef)[q]]
    }
  }
  
  names(cov.values) <- names.cov
  cov.values.mat <- matrix(data = unlist(cov.values), nrow = nvals, ncol = ncoefs)
  
  
  if(identical(se.fit, TRUE)) {
    ##substitute a given row for each covariate
    predicted.SE <- matrix(NA, nrow = nvals, ncol = 2)
    colnames(predicted.SE) <- c("Pred.value", "SE")
    rownames(predicted.SE) <- 1:nvals
    part.devs.eval <- list( )
    part.devs.eval[[1]] <- 1
    for (w in 1:nvals) {
      if(int.yes && ncovs > 1) {
        for (p in 2:ncoefs) {
          part.devs.eval[[p]] <-  cov.values[names.cov[p]][[1]][w]
        }
      } # else {  ##for cases without intercept
      ##for (p in 1:ncoefs) {
      ##  part.devs.eval[[p]] <-  cov.values[names.cov[p]][[1]][w]
      ##}
      ##}
      
      part.devs.solved <- unlist(part.devs.eval)
      
      ##extract vc matrix
      vcmat <- vcov(mod)
      
      mat_partialdevs<-data.matrix(part.devs.solved) #create matrix from vector of 2 rows by 1 column
      mat_tpartialdevs<-t(part.devs.solved)        #transpose of partial derivatives to have 2 columns by 1 row
      
      #####5)
      var_hat<-mat_tpartialdevs%*%vcmat%*%mat_partialdevs
      SE<-sqrt(var_hat)
      predicted.vals <- fix.coef%*%cov.values.mat[w,]
      predicted.SE[w, 1] <- predicted.vals
      predicted.SE[w, 2] <- SE
    }
    
    
    out.fit.SE <- list(fit = predicted.SE[,"Pred.value"], se.fit = predicted.SE[, "SE"])
    
  } else {
    predicted.SE <- matrix(NA, nrow = nvals, ncol = 1)
    colnames(predicted.SE) <- c("Pred.value")
    rownames(predicted.SE) <- 1:nvals
    for (w in 1:nvals) {
      predicted.vals <- fix.coef%*%cov.values.mat[w,]
      predicted.SE[w, 1] <- predicted.vals
    }
    
    out.fit.SE <- predicted.SE
    colnames(out.fit.SE) <- "fit"
    
  }
  return(out.fit.SE)
}

#####################################

