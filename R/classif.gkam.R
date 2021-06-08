#' Classification Fitting Functional Generalized Kernel Additive Models
#' 
#' Computes functional classification using functional explanatory variables
#' using backfitting algorithm.
#' 
#' @details The first item in the \code{data} list is called \emph{"df"} and is a data
#' frame with the response, as \code{\link{glm}}.\cr Functional covariates of
#' class \code{fdata} are introduced in the following items in the \code{data}
#' list.
#' 
#' @param formula an object of class \code{formula} (or one that can be coerced
#' to that class): a symbolic description of the model to be fitted. The
#' procedure only considers functional covariates (not implemented for
#' non-functional covariates). The details of model specification are given
#' under \code{Details}.
#' @param data List that containing the variables in the model.
#' @param family a description of the error distribution and link function to
#' be used in the model. This can be a character string naming a family
#' function, a family function or the result of a call to a family function.
#' (See \code{\link{family}} for details of family functions.)
#' @param weights Weights:   
#' \itemize{
#' \item if \code{character} string \code{='equal'} same weights for each observation (by default) and
#' \code{='inverse'} for inverse-probability of weighting.   
#' \item if \code{numeric} vector of length \code{n}, Weight values of each observation.
#' }
#' @param par.metric List of arguments by covariable to pass to the
#' \code{metric} function by covariable.
#' @param par.np List of arguments to pass to the \code{fregre.np.cv} function
#' @param offset this can be used to specify an a priori known component to be
#' included in the linear predictor during fitting.
#' @param control a list of parameters for controlling the fitting process, by
#' default: maxit, epsilon, trace and inverse.
#' @param prob probability value used for binary discriminant.
#' @param type If type is\code{"1vsall"}  (by default) 
#' a maximum probability scheme is applied: requires G binary classifiers.
#' If type is \code{"majority"}  (only for multicalss classification G > 2) 
#' a voting scheme is applied: requires  G (G - 1) / 2 binary classifiers.
#' @param \dots Further arguments passed to or from other methods.
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
#' }
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' @seealso See Also as: \code{\link{fregre.gkam}}.\cr Alternative method:
#' \code{\link{classif.glm}}.
#' @references Febrero-Bande M. and Gonzalez-Manteiga W. (2012).
#' \emph{Generalized Additive Models for Functional Data}. TEST.
#' Springer-Velag.  \doi{10.1007/s11749-012-0308-0}
#' 
#' McCullagh and Nelder (1989), \emph{Generalized Linear Models} 2nd ed.
#' Chapman and Hall.
#' 
#' Opsomer J.D. and Ruppert D.(1997). \emph{Fitting a bivariate additive model
#' by local polynomial regression}.Annals of Statistics, \code{25}, 186-211.
#' @keywords classif
#' @examples
#' \dontrun{ 
#' ## Time-consuming: selection of 2 levels 
#' data(phoneme)
#' mlearn<-phoneme[["learn"]][1:150]
#' glearn<-factor(phoneme[["classlearn"]][1:150])
#' dataf<-data.frame(glearn)
#' dat=list("df"=dataf,"x"=mlearn)
#' a1<-classif.gkam(glearn~x,data=dat)
#' summary(a1)
#' mtest<-phoneme[["test"]][1:150]
#' gtest<-factor(phoneme[["classtest"]][1:150])
#' newdat<-list("x"=mtest)
#' p1<-predict(a1,newdat)
#' table(gtest,p1)
#' }  
#' @export
classif.gkam=function(formula,data, weights = "equal", family = binomial(),
 par.metric = NULL,par.np=NULL, offset=NULL,prob=0.5,type= "1vsall",
 control = NULL,...){
  if (is.null(control)) control =list(maxit = 100,epsilon = 0.001, trace = FALSE, inverse="solve")
  C<-match.call()
  a<-list()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula","family","data","weigths","par.metric","par.np","offset","control"), names(mf),0L)
  tf <- terms.formula(formula)
  terms <- attr(tf, "term.labels")
  nt <- length(terms)
  if (attr(tf, "response") > 0) {
          response <- as.character(attr(tf, "variables")[2])
          pf <- rf <- paste(response, "~", sep = "")
      } else pf <- rf <- "~"
  newy<-y<-data$df[[response]]
      nobs <- if (is.matrix(y))   nrow(y)
              else length(y)
  if (!is.factor(y)) y<-as.factor(y)
  n<-length(y)
  newdata<-data
  ny<-levels(y)
  probs<-ngroup<-nlevels(y)
  prob.group<-array(NA,dim=c(n,ngroup))
  colnames(prob.group)<-ny
  
  if (is.character(weights)) {
    weights<-weights4class(y,type=weights)
  } else {
    if (length(weights)!=n) 
      stop("length weights != length response")
  }
  
  if (ngroup==2) {
        newy<-ifelse(y==ny[1],0,1)
        newdata$df[[response]]<-newy
        a[[1]]<-fregre.gkam(formula,family=family,data=newdata,weights=weights,par.metric=par.metric,
                            par.np=par.np,offset=offset,control=control)
        out2glm <- classif2groups(a,y,prob,ny)
  } else {
  a<-list()
  if (type == "majority"){
    cvot<-combn(ngroup,2)
    nvot<-ncol(cvot)
    votos<-matrix(0,n,ngroup)
    colnames(votos) <- ny
    class(data)<-c("ldata","list")
    b0<-list()
    for (ivot in 1:nvot) {  
      ind1 <- y==ny[cvot[1,ivot]]
      ind2 <- y==ny[cvot[2,ivot]] 
      #i2a2<-which(ind1 | ind2)
      i2a2 <- ind1 | ind2
      newy<-rep(NA,n)   
      newy[ind1 ]<- 1
      newy[ind2 ]<- 0
      newdata<-data
      newdata$df[response] <- newy
      #newdata <- newdata[i2a2,row=TRUE]
      newdata<-subset(newdata,i2a2)
      class(newdata)<-c("list")
      a[[ivot]]<-fregre.gkam(formula,family=family,data=newdata
                , weigths=weights[i2a2], par.metric=par.metric
                ,par.np=par.np,offset=offset,control=control,...)
      prob.log <- a[[ivot]]$fitted.values  > prob
      votos[i2a2, cvot[1,ivot]] <- votos[i2a2, cvot[1,ivot]] + as.numeric(prob.log)
      votos[i2a2, cvot[2,ivot]] <- votos[i2a2, cvot[2,ivot]] + as.numeric(!prob.log)
    }
    out2glm<-classifKgroups(y,votos,ny)
  }
  else { # One vs Other
    prob.group<-array(NA,dim=c(n,ngroup))
    colnames(prob.group)<-ny
    for (i in 1:ngroup) {
      igroup  <- y==ny[i]
      newy<-ifelse(igroup, 1, 0)
      weights0 <- weights
      newdata$df[response]<-newy
      a[[i]]<-fregre.gkam(formula,family=family,data=newdata,weigths=weights0
                      ,par.metric=par.metric,par.np=par.np
                      ,offset=offset,control=control,...)
      prob.group[,i]<-a[[i]]$fitted.values
    }
    out2glm<-classifKgroups(y,prob.group,ny)
  }
}
  yest <- out2glm$yest
  prob.group <- out2glm$prob.group
  max.prob=mean(yest==y)
  #max.prob=sum(diag(tab))/sum(tab)
  output<-list(formula=formula,data=data,group=y,group.est=yest,
               prob.classification=out2glm$prob1,prob.group=prob.group,C=C,
               m=m,max.prob=max.prob,fit=a,prob = prob,type=type)
  class(output) <- "classif"
  return(output)
}
