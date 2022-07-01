#' @title Predict method for functional linear model
#' 
#' @description 
#' Computes predictions for regression between functional (and non functional)
#' explanatory variables and scalar response. 
#' \itemize{ 
#' \item \code{predict.fregre.lm}, Predict method for functional linear model of
#' \code{\link{fregre.lm}} fits object using basis or principal component
#' representation.
#' \item \code{predict.fregre.plm}, Predict method for
#' semi-functional linear regression model of \code{\link{fregre.plm}} fits
#' object using using asymmetric kernel estimation. 
#' \item \code{predict.fregre.glm}, Predict method for functional generalized linear
#' model of \code{\link{fregre.glm}} fits object using basis or principal
#' component representation. 
#' \item \code{predict.fregre.gsam}, Predict method for functional generalized 
#' spectral additive model of \code{\link{fregre.gsam}} fits object using basis 
#' or principal component representation.
#' \item \code{predict.fregre.gkam}, Predict method for functional generalized 
#' kernel additive model of \code{\link{fregre.gkam}} fits object using 
#' backfitting algorithm. 
#' }
#' 
#' These functions use the model fitting function \code{\link{lm}},
#' \code{\link{glm}} or \code{\link{gam}} properties.\cr If using functional
#' data derived, is recommended to use a number of bases to represent beta
#' lower than the number of bases used to represent the functional data. \cr
#' The first item in the \code{data} list of \code{newx} argument is called
#' \emph{"df"} and is a data frame with the response and non functional
#' explanatory variables, as \code{\link{lm}}, \code{\link{glm}} or
#' \code{\link{gam}}. Functional variables (\code{fdata} and \code{fd} class)
#' are introduced in the following items in the \code{data} list of \code{newx}
#' argument.
#' 
#' @aliases predict.fregre.lm predict.fregre.plm predict.fregre.glm
#' predict.fregre.gsam predict.fregre.gkam
#' @param object \code{fregre.lm}, \code{fregre.plm}, \code{fregre.glm},
#' \code{fregre.gsam}\cr or \code{fregre.gkam} object.
#' @param newx An optional data list in which to look for variables with which
#' to predict. If omitted, the fitted values are used. List of new explanatory
#' data.
#' @param type a character vector, Type of prediction: (\code{response}, \code{terms} for model terms or \code{effects}  for model terms where
#' the partial effects are summarized for each functional variable.
#' @param se.fit =TRUE (not default) standard error estimates are returned for
#' each prediction.
#' @param scale Scale parameter for std.err. calculation.
#' @param df Degrees of freedom for scale.
#' @param interval Type of interval calculation.
#' @param level Tolerance/confidence level.
#' @param pred.var the variance(s) for future observations to be assumed for
#' prediction intervals. See \code{link{predict.lm}} for more details.
#' @param weights variance weights for prediction. This can be a numeric vector
#' or a one-sided model formula. In the latter case, it is interpreted as an
#' expression evaluated in newdata
#' @param \dots Further arguments passed to or from other methods.
#' @return Return the predicted values and optionally:
#' \itemize{
#' \item {predict.lm,predict.glm,predict.gam}{ produces a vector of predictions
#' or a matrix of predictions and bounds with column names fit, lwr, and upr if
#' interval is set. If se.fit is TRUE, a list with the following components is
#' returned: fit vector or matrix as above.} 
#' \item {se.fit}{ standard error of predicted means.} 
#' \item {residual.scale}{ residual standard deviations.}
#' \item {df}{ degrees of freedom for residual.}
#' }
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@usc.es}
#' 
#' @seealso See Also as: \code{\link{fregre.lm}}, \code{\link{fregre.plm}},
#' \code{\link{fregre.glm}}, \code{\link{fregre.gsam}} and
#' \code{\link{fregre.gkam}}. 
#' 
#' @references Febrero-Bande, M., Oviedo de la Fuente, M. (2012).
#' \emph{Statistical Computing in Functional Data Analysis: The R Package
#' fda.usc.} Journal of Statistical Software, 51(4), 1-28.
#' \url{https://www.jstatsoft.org/v51/i04/}
#' 
#' @keywords regression
#' @examples
#' \dontrun{
#' data(tecator)
#' ind <- 1:129
#' x <- tecator$absorp.fdata
#' x.d2 <- fdata.deriv(x,nderiv=2)
#' tt <- x[["argvals"]]
#' dataf <- as.data.frame(tecator$y)
#' ldat <- ldata("df"=dataf[ind,],"x.d2"=x.d2[ind])
#' basis.x <- list("x.d2"=create.pc.basis(ldat$x.d2))
#' res <- fregre.gsam(Fat~s(Water,k=3)+s(x.d2,k=3),data=ldat,
#'                    family=gaussian(),basis.x=basis.x)
#' newldat <- list("df"=dataf[-ind,],"x.d2"=x.d2[-ind])
#' pred <- predict(res,newldat)
#' plot(pred,tecator$y$Fat[-ind])
#' res.glm <- fregre.glm(Fat  ~ Water+x.d2, data=ldat,
#'                   family=gaussian(),basis.x=basis.x)
#' pred.glm <- predict(res.glm,newldat)
#' newy <- tecator$y$Fat[-ind]
#' points(pred.glm,tecator$y$Fat[-ind],col=2)
#' 
#' # Time-consuming 
#' res.gkam <- fregre.gkam(Fat ~ x.d2, data=ldat)
#' pred.gkam <- predict(res.gkam,newldata)
#' points(pred.gkam,tecator$y$Fat[-ind],col=4)
#' 
#' ((1/length(newy))*sum((drop(newy)-pred)^2))/var(newy)
#' ((1/length(newy))*sum((newy-pred.glm)^2))/var(newy)    
#' ((1/length(newy))*sum((newy-pred.gkam)^2))/var(newy)    
#' }                                                                                                              
#' @rdname predict.fregre.lm
#' @export 
predict.fregre.lm<-function (object, newx = NULL, type = "response", se.fit = FALSE, 
          scale = NULL, df = df, interval = "none", level = 0.95, weights = 1, 
          pred.var = res.var/weights, ...) 
{
  if (is.null(object)) 
    stop("No fregre.lm object entered")
  if (is.null(newx)) {
    if (type == "effects"){
      fake  = predict.lm(object, type = "terms", se.fit = se.fit, 
                      interval = interval, level = level, weights = weights, 
                      pred.var = pred.var, df = df, scale = scale, ...)
      yp <- effect.fake(object,fake)
    } else{
      yp = predict.lm(object, type = type, se.fit = se.fit, 
                      interval = interval, level = level, weights = weights, 
                      pred.var = pred.var, df = df, scale = scale, ...)
    }
    return(yp)
  }
  else {
    name.coef <- NULL
    data = newx
    basis.x = object$basis.x
    basis.b = object$basis.b
    formula = object$formula.ini
    tf <- terms.formula(formula)
    terms <- attr(tf, "term.labels")
    nt <- length(terms)
    vtab <- rownames(attr(tf, "factors"))
    vnf = intersect(terms, names(data$df))
    vfunc2 = setdiff(terms, vnf)
    vint = setdiff(terms, vtab)
    vfunc = setdiff(vfunc2, vint)
    off <- attr(tf, "offset")
    beta.l = list()
    kterms = 1
#    if (attr(tf, "response") > 0) {
#      response <- as.character(attr(tf, "variables")[2])
#      pf <- rf <- paste(response, "~", sep = "")
#    }
#    else pf <- rf <- "~"
    pf <- rf <- "~"  # En predicción no hace falta la respuesta en los datos nuevos
    if (attr(tf, "intercept") == 0) { 
      print("No intecept")
      pf <- paste(pf, -1, sep = "")
    }
    if (length(vnf) > 0) {
      first = FALSE
      for (i in 1:length(vnf)) {
        if (kterms > 1) 
          pf <- paste(pf, "+", vnf[i], sep = "")
        else pf <- paste(pf, vnf[i], sep = "")
        kterms <- kterms + 1
      }
      if (attr(tf, "intercept") == 0) {
        pf <- paste(pf, -1, sep = "")
      }
      mf <- as.data.frame(model.matrix(formula(pf), data$df))
      vnf2 <- names(mf)[-1]
      for (i in 1:length(vnf2)) pf <- paste(pf, "+", vnf2[i], 
                                            sep = "")
      XX <- mf
    }
    else {
      pf2 <- paste(pf, "1", sep = "")
      XX <- data.frame(model.matrix(formula(pf2), data$df))
      first = TRUE
    }
    if (length(vnf) > 0) {
      spm <- matrix(object$coefficients[names(XX)], ncol = 1)
      yp <- as.matrix(XX) %*% spm
    }
    else yp <- object$coefficients[1] * rep(1, len = nrow(newx[[vfunc[1]]]))
    lenfunc <- length(vfunc)
    if (lenfunc > 0) {
      for (i in 1:lenfunc) {
#        if (lenfunc>0) {
#          k=1
#          mean.list=vs.list=JJ=list()
#          for (i in 1:lenfunc) {
#            if(class(newx[[vfunc[i]]])[1]=="fdata"){
            if(inherits(newx[[vfunc[i]]],"fdata")){
#              tt<-data[[vfunc[i]]][["argvals"]]
#              rtt<-data[[vfunc[i]]][["rangeval"]]
              fdataobj<-data[[vfunc[i]]]
              fdat<-data[[vfunc[i]]];      dat<-fdataobj$data
#              if (nrow(dat)==1) rwn<-NULL         else rwn<-rownames(dat)
              if (nrow(newx[[vfunc[i]]]$data)==1) rwn<-NULL else rwn<-rownames(newx[[vfunc[i]]]$data)
              #  if (basis.x[[vfunc[i]]]$type=="pc" 
              #      | basis.x[[vfunc[i]]]$type=="pls")
              #    bsp1=FALSE      else bsp1 <- TRUE
              xaux <- fdata2basis(newx[[vfunc[i]]],basis.x[[vfunc[i]]])

              Z <- xaux$coefs%*%object$vs.list[[vfunc[i]]]
              colnames(Z)=colnames(vs.list[[vfunc[i]]])
              name.coef[[vfunc[i]]] <- paste(vfunc[i],".",colnames(Z),sep="")
              if (first) {
                XX=Z
                first=FALSE
              }   else {
                XX = cbind(XX,Z)} 
          }
        else {
          if (class(data[[vfunc[i]]])[1] == "fd") {
            if (class(object$basis.x[[vfunc[i]]]) != "pca.fd") {
#              x.fd <- fdataobj <- data[[vfunc[i]]]
              x.fd <- newx[[vfunc[i]]]
              r = x.fd[["basis"]][["rangeval"]]
              J <- object$basis.list[[vfunc[i]]]
              x.fd$coefs <- x.fd$coefs - object$mean[[vfunc[i]]]$coefs[,1]
              Z = t(x.fd$coefs) %*% J
              colnames(Z) = colnames(J)
            }
            else {
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
            }
            else XX = cbind(XX, Z)
          }
          else stop("Please, enter functional covariate")
        }
      }
    }

    nn <- nrow(XX)
    if (!is.data.frame(XX)) 
      XX = data.frame(XX)
    if (!object$lambda) {
      res.var <- object$sr2
      if (type == "effects"){
        fake  = predict.lm(object, newdata = XX, type = "terms", se.fit = se.fit, 
                           interval = interval, level = level, weights = weights, 
                           pred.var = pred.var, df = df, scale = scale, ...)
        return(effect.fake(object,fake))
      } else{
        return(predict.lm(object = object, newdata = XX, 
                          type = type, se.fit = se.fit, interval = interval, 
                          level = level, weights = weights, pred.var = pred.var, 
                          df = df, scale = scale, , ...))
        }
    }
    else {
      if (type!="response") warning("Only response type implemented for penalization")
      for (i in 1:lenfunc) {
        if (object$call[[1]] == "fregre.pls") 
          return(predict.lm(object = object, newdata = XX, 
                            type = type, se.fit = se.fit, ...))
        if (object$basis.x[[vfunc[i]]]$type == "pc") {
          object$beta.l[[vfunc[i]]]$data <- matrix(object$beta.l[[vfunc[i]]]$data, 
                                                   nrow = 1)
          b1 <- inprod.fdata(fdata.cen(newx[[vfunc[i]]], 
                                       object$mean.list[[vfunc[i]]])[[1]], object$beta.l[[vfunc[i]]])
          yp <- yp + b1
        }
        else {
          xcen <- fdata.cen(newx[[vfunc[i]]], object$mean.list[[vfunc[i]]])[[1]]
          x.fd = Data2fd(argvals = xcen$argvals, y = t(xcen$data), 
                         basisobj = object$basis.x[[vfunc[i]]])
          C = t(x.fd$coefs)
          cnames <- colnames(object$basis.list[[vfunc[i]]])
          b.est <- matrix(object$coefficients[cnames], 
                          ncol = 1)
          b1 <- C %*% object$basis.list[[vfunc[i]]] %*% b.est
          yp <- yp + b1
        }
      }
      XX2 <- as.matrix(cbind(rep(1, len = nn), XX))
      predictor <- drop(yp)
      if (se.fit || interval != "none") {
        ip <- rowSums((XX2 %*% object$Vp * XX2))
        res.var <- object$sr2
        df <- object$df.residual
        if (interval != "none") {
          tfrac <- qt((1 - level)/2, df)
          hwid <- tfrac * switch(interval, confidence = sqrt(ip), 
                                 prediction = sqrt(ip + pred.var))
          predictor <- cbind(predictor, predictor + hwid %o% 
                               c(1, -1))
          colnames(predictor) <- c("fit", "lwr", "upr")
        }
      }
      if (se.fit) {
        se <- sqrt(ip)
        return(list(fit = predictor, se.fit = se, df = df, 
                    residual.scale = sqrt(res.var)))
      }
      else return(predictor)
    }
  }
  return(drop(yp))
}
#################################
#################################
effect.fake <- function(object,terms){
  #fake<-predict(object,type = "terms")
  fake <- terms
  vfunc <- names(object$basis.list)
  nfunc <- length(vfunc)
  tr <- attr(object$terms,"term.labels")
  #effects.df <- intersect(tr,colnames(object$data$df))
  #effects<-fake[,effects.df ,drop=F]
  effects<-NULL
  vf <- NULL
  for (i in 1:nfunc){
    ifunc<-colnames(object$basis.list[[i]])
    dfnames <- intersect(tr,ifunc)
    vf <- c(vf,dfnames)
    effects<-cbind(effects,rowSums(fake[,dfnames,drop=F]))
  }
  colnames(effects)<-vfunc
  dfnames <- setdiff(tr,vf)
  effects<-cbind(fake[,dfnames,drop=F],effects)
  effects  
}
#################################
# Modification (2021/03/22)
# Z <- inprod.fdata(newXcen, object$vs.list[[vfunc[i]]])
# #Z <- inprod.fdata(newXcen, object$JJ[[vfunc[i]]])
# if (!object$lambda)
# # if (!object$rn)