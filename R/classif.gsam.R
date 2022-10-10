#' Classification Fitting Functional Generalized Additive Models
#' 
#' @description Computes functional classification using functional (and non functional)
#' explanatory variables by basis representation.
#' 
#' @details The first item in the \code{data} list is called \emph{"df"} and is a data
#' frame with the response and non functional explanatory variables, as
#' \code{\link{glm}}.\cr
#' 
#' Functional covariates of class \code{fdata} or \code{fd} are introduced in
#' the following items in the \code{data} list.\cr \code{basis.x} is a list of
#' basis for represent each functional covariate. The basis object can be
#' created by the function: \code{\link{create.pc.basis}}, \code{\link{pca.fd}}
#' \code{\link{create.pc.basis}}, \code{\link{create.fdata.basis}} o
#' \code{\link{create.basis}}.
#' 
#' @aliases  classif.gsam  
#' @param formula an object of class \code{formula} (or one that can be coerced
#' to that class): a symbolic description of the model to be fitted. The
#' details of model specification are given under \code{Details}.
#' @param family a description of the error distribution and link function to
#' be used in the model. This can be a character string naming a family
#' function, a family function or the result of a call to a family function.
#' (See \code{\link{family}} for details of family functions.)
#' @param data List that containing the variables in the model.
#' @param weights Weights:   
#' \itemize{
#' \item if \code{character} string \code{='equal'} same weights for each observation (by default) and
#' \code{='inverse'} for inverse-probability of weighting.   
#' \item if \code{numeric} vector of length \code{n}, Weight values of each observation.
#' }
#' @param basis.x List of basis for functional explanatory data estimation.
# @param basis.b List of basis for functional beta parameter estimation.
#' @param CV =TRUE, Cross-validation (CV) is done.
#' @param prob probability value used for binari discriminant.
#' @param type If type is\code{"1vsall"}  (by default) 
#' a maximum probability scheme is applied: requires G binary classifiers.
#' If type is \code{"majority"}  (only for multicalss classification G > 2) 
#' a voting scheme is applied: requires  G (G - 1) / 2 binary classifiers.
#' @param \dots Further arguments passed to or from other methods.
# @param  boos,classif,par.classif,mfinal,coeflearn,control,verbose,n.cores,par.metric,par.np,offset aaaa
# @param size,epsilon,trace,inverse,fdata1,fdata2,w bbbb
# @param ldata,lfdata,lfdataref,parallel,object,newdata,new.fdataobj,newx ccc
# @param tt,h,p,Ker,cv,class.weights,ncores,metric,method,x,digits dddd
# @param par.depth,A aaaccc
#' @return Return \code{gam} object plus:
#' \itemize{
#' \item \code{formula}{ formula.}
#' \item \code{data}{ List that containing the variables in the model.} 
#' \item \code{group}{ Factor of length \emph{n}} 
#' \item \code{group.est}{ Estimated vector groups}
#' \item \code{prob.classification}{ Probability of correct classification by group.}
#' \item \code{prob.group}{ Matrix of predicted class probabilities. For each
#' functional point shows the probability of each possible group membership.}
#' \item \code{max.prob}{ Highest probability of correct classification.}
#' \item \code{type}  Type of classification scheme: 1 vs all  or majority voting.
#' }
#' @note If the formula only contains a non functional explanatory variables
#' (multivariate covariates), the function compute a standard \code{\link{glm}}
#' procedure.
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' @seealso See Also as: \code{\link{fregre.gsam}}.\cr Alternative method:
#' \code{\link{classif.np}}, \code{\link{classif.glm}} and
#' \code{\link{classif.gkam}}.
#' @references Ramsay, James O., and Silverman, Bernard W. (2006), \emph{
#' Functional Data Analysis}, 2nd ed., Springer, New York.
#' 
#' McCullagh and Nelder (1989), \emph{Generalized Linear Models} 2nd ed.
#' Chapman and Hall.
#' 
#' Venables, W. N. and Ripley, B. D. (2002) \emph{Modern Applied Statistics
#' with S}, New York: Springer.  %Wood (2001) mgcv:GAMs and Generalized Ridge
#' Regression for R. R News 1(2):20-25
#' @keywords classif
#' @examples
#' \dontrun{
#' data(phoneme)
#' ldat <- ldata("df" = data.frame(y = phoneme[["classlearn"]]),
#'               "x" = phoneme[["learn"]])
#'               classifKgroups <- fda.usc:::classifKgroups
#' a1 <- classif.gsam( y ~ s(x,k=3),data=ldat)
#' summary(a1)
#' newldat <- ldata("df" = data.frame(y = phoneme[["classtest"]]),
#'               "x" = phoneme[["test"]])
#' p1 <- predict(a1,newldat)
#' table(newldat$df$y,p1)
#' sum(p1==newldat$df$y)/250
#' }
#' @export classif.gsam
classif.gsam <- function (formula, data, family = binomial(),  
                        weights = "equal", basis.x = NULL,  
                        CV = FALSE, prob=0.5, type = "1vsall",...) 
{
  C <- match.call()
  a <- list()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "family", "basis.x",
               "CV"), names(mf), 0L)
  tf <- terms.formula(formula)
  terms <- attr(tf, "term.labels")
  nt <- length(terms)
  if (attr(tf, "response") > 0) {
    response <- as.character(attr(tf, "variables")[2])
    pf <- rf <- paste(response, "~", sep = "")
  }
  else pf <- rf <- "~"
  newy <- y <- data$df[[response]]
  if (!is.factor(y)) 
    y <- as.factor(y)
  n <- length(y)
  
  #if (is.null(weights)) weights <- rep(1,n)
  if (is.character(weights)) {
    weights<-weights4class(y,type=weights)
  } else {
    if (length(weights)!=n) 
      stop("length weights != length response")
  }
  newdata <- data
  ny <- levels(y)
  prob2<-prob1 <- ngroup <- nlevels(y)
  prob.group <- array(NA, dim = c(n, ngroup))
  if (ngroup == 2) {
    #newy <- ifelse(y == ny[1], 1, 0)
    newy <- y
    newdata$df$y <- newy
    a[[1]] <- fregre.gsam(formula, data = newdata, family = family, 
                          weights = weights, basis.x = basis.x, 
                          basis.b = basis.x, CV = CV, ...)
    out2glm <- classif2groups(a,y,prob,ny)
  }
  else {
    a<-list()
    if (type == "majority"){
      cvot<-combn(ngroup,2)
      nvot<-ncol(cvot)
      votos<-matrix(0,n,ngroup)
      b0<-list()
      for (ivot in 1:nvot) {  
        ind1 <- y==ny[cvot[1,ivot]]
        ind2 <- y==ny[cvot[2,ivot]] 
        i2a2<-which(ind1 | ind2)
        newy<-rep(NA,n)   
        newy[ind1 ]<- 1
        newy[ind2 ]<- 0
        newdata$df$newy<-newy
        newdata$df[response] <- newy
        a[[ivot]]<-fregre.gsam(formula,data=newdata,family=family, 
                               weights = weights,basis.x=basis.x,
                              basis.b=basis.x, CV = CV, subset = i2a2,...)
        prob.log <- a[[ivot]]$fitted.values  > prob
        votos[i2a2, cvot[1,ivot]] <- votos[i2a2, cvot[1,ivot]] + as.numeric(prob.log)
        votos[i2a2, cvot[2,ivot]] <- votos[i2a2, cvot[2,ivot]] + as.numeric(!prob.log)
      }
      out2glm <- classifKgroups(y,votos,ny)
    }
    else { # One vs Other
      prob.group<-array(NA,dim=c(n,ngroup))
      colnames(prob.group)<-ny
      for (i in 1:ngroup) {
        igroup  <- y==ny[i]
        newy<-ifelse(igroup, 1, 0)
        weights0 <- weights
        weights0[igroup] <- weights0[igroup]/ sum(weights0[igroup])
        weights0[!igroup] <- weights0[!igroup]/sum(weights0[!igroup])
        newdata$df[response]<-newy
        a[[i]]<-fregre.gsam(formula,data=newdata,family=family,
                            weights = weights0,basis.x=basis.x,
                            basis.b=basis.x, CV=CV,...)
        prob.group[,i] <- a[[i]]$fitted.values
      }
      out2glm <- classifKgroups(y, prob.group, ny)
    }
  }
  yest <- out2glm$yest
  prob1 <- out2glm$prob1
  prob.group <- out2glm$prob.group
  max.prob <- mean(yest==y)
  output <- list(formula = formula, data = data, group = y,  levels = ny, 
                 group.est = yest, prob.classification = prob1, prob.group = prob.group, 
                 C = C, m = m, max.prob = max.prob, fit = a, prob = prob, type = type)
 if (ngroup > 2){
   if (type=="majority") {
    output$voto <- votos
    output$fit <- a
   }
 }
  class(output) <- "classif"
  return(output)
}

