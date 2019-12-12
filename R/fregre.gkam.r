#' Fitting Functional Generalized Kernel Additive Models.
#' 
#' Computes functional regression between functional explanatory variables
#' \eqn{(X^{1}(t_1),...,X^{q}(t_q))}{(X(t_1),...,X(t_q))} and scalar response
#' \eqn{Y} using backfitting algorithm. 
#' 
#' @details The smooth functions \eqn{f(.)} are estimated nonparametrically using a
#' iterative local scoring algorithm by applying Nadaraya-Watson weighted
#' kernel smoothers using \code{\link{fregre.np.cv}} in each step, see
#' Febrero-Bande and Gonzalez-Manteiga (2011) for more details.\cr 
#' Consider the fitted response \eqn{\hat{Y}=g^{-1}(H_{Q}y)}{g(Y.est)=Hy},
#' where \eqn{H_{Q}}{H} is the weighted hat matrix.\cr Opsomer and Ruppert
#' (1997) solves a system of equations for fit the unknowns
#' \eqn{f(\cdot)}{f(.)} computing the additive smoother matrix \eqn{H_k}{H_k}
#' such that \eqn{\hat{f}_k (X^k)=H_{k}Y}{f.est_k(X_k)=H_k Y} and
#' \eqn{H_Q=H_1+,\cdots,+H_q}{H= H_1+,...,+H_q}. The additive model is fitted
#' as follows: \deqn{\hat{Y}=g^{-1}\Big(\sum_i^q
#' \hat{f_i}(X_i)\Big)}{g(y.est)=\sum(i:q) f.est_i(X_i)}
#' 
#' @aliases fregre.gkam 
#' @param formula an object of class \code{formula} (or one that can be coerced
#' to that class): a symbolic description of the model to be fitted. The
#' procedure only considers functional covariates (not implemented for
#' non-functional covariates). The details of model specification are given
#' under \code{Details}.
#' @param data List that containing the variables in the model.
#' @param family a description of the error distribution and link function to
#' be used in the model. This can be a character string naming a family
#' function, a family function or the result of a call to a family function.
#' (See \code{\link{family}} for details of family functions).
#' @param weights weights
#' @param par.metric List of arguments by covariate to pass to the
#' \code{metric} function by covariate.
#' @param par.np List of arguments to pass to the \code{fregre.np.cv} function
#' @param offset this can be used to specify an a priori known component to be
#' included in the linear predictor during fitting.
#' @param control a list of parameters for controlling the fitting process, by
#' default: \code{maxit}, \code{epsilon}, \code{trace} and \code{inverse}
#' @param inverse ="svd" (by default) or ="solve" method.
#' @param \dots Further arguments passed to or from other methods.
#' @return \itemize{
#' \item \code{result}{ List of non-parametric estimation by covariate.}
#' \item \code{fitted.values}{ Estimated scalar response.} 
#' \item \code{residuals}{ \code{y} minus \code{fitted values}.} 
#' \item \code{effects}{ The residual degrees of freedom.} 
#' \item \code{alpha}{ Hat matrix.} 
#' \item \code{family}{ Coefficient of determination.} 
#' \item \code{linear.predictors}{ Residual variance.}
#' \item \code{deviance}{ Scalar response.} 
#' \item \code{aic}{ Functional explanatory data.}
#' \item \code{null.deviance}{ Non functional explanatory data.} 
#' \item \code{iter}{ Distance matrix between curves.} 
#' \item \code{w}{ beta coefficient estimated}
#' \item \code{eqrank}{ List that containing the variables in the model.}
#' \item \code{prior.weights}{ Asymmetric kernel used.} 
#' \item \code{y}{ Scalar response.}
#' \item \code{H}{ Hat matrix, see Opsomer and Ruppert(1997) for more details.}
#' \item \code{converged}{ conv.}
#' }
#' @author Febrero-Bande, M. and Oviedo de la Fuente, M.
#' @seealso See Also as: \code{\link{fregre.gsam}}, \code{\link{fregre.glm}}
#' and \code{\link{fregre.np.cv}}\cr
#' @references Febrero-Bande M. and Gonzalez-Manteiga W. (2012).
#' \emph{Generalized Additive Models for Functional Data}. TEST.
#' Springer-Velag.  \url{http://dx.doi.org/10.1007/s11749-012-0308-0}
#' 
#' Opsomer J.D. and Ruppert D.(1997). \emph{Fitting a bivariate additive model
#' by local polynomial regression}.Annals of Statistics, \code{25}, 186-211.
#' @keywords regression
#' @examples 
#' \dontrun{
#' data(tecator)
#' ab=tecator$absorp.fdata[1:100]
#' ab2=fdata.deriv(ab,2)
#' yfat=tecator$y[1:100,"Fat"]
#' 
#' # Example 1: # Changing the argument par.np and family
#' yfat.cat=ifelse(yfat<15,0,1)
#' xlist=list("df"=data.frame(yfat.cat),"ab"=ab,"ab2"=ab2)
#' f2<-yfat.cat~ab+ab2
#' 
#' par.NP<-list("ab"=list(Ker=AKer.norm,type.S="S.NW"),
#' "ab2"=list(Ker=AKer.norm,type.S="S.NW"))
#' res2=fregre.gkam(f2,family=binomial(),data=xlist,
#' par.np=par.NP)
#' res2
#' 
#' # Example 2: Changing the argument par.metric and family link
#' par.metric=list("ab"=list(metric=semimetric.deriv,nderiv=2,nbasis=15),
#' "ab2"=list("metric"=semimetric.basis))
#' res3=fregre.gkam(f2,family=binomial("probit"),data=xlist,
#' par.metric=par.metric,control=list(maxit=2,trace=FALSE))
#' summary(res3)
#' 
#' # Example 3: Gaussian family (by default)
#' # Only 1 iteration (by default maxit=100)
#' xlist=list("df"=data.frame(yfat),"ab"=ab,"ab2"=ab2)
#' f<-yfat~ab+ab2
#' res=fregre.gkam(f,data=xlist,control=list(maxit=1,trace=FALSE))
#' res
#' }
#' @export
fregre.gkam=function (formula,family = gaussian(),data, weights= rep(1,nobs),
     par.metric = NULL,par.np=NULL,offset=NULL,
     control = list(maxit = 100,epsilon = 0.001, trace = FALSE, inverse="solve"),...)
{
 tf <- terms.formula(formula)
 terms <- attr(tf, "term.labels")
 nt <- length(terms)
 if (attr(tf, "response") > 0) {
        response <- as.character(attr(tf, "variables")[2])
        pf <- rf <- paste(response, "~", sep = "")
    } else pf <- rf <- "~"
 vtab<-rownames(attr(tf,"factors"))
 vnf=intersect(terms,names(data$df))
 vnf2=intersect(vtab[-1],names(data$df)[-1])
 vfunc2=setdiff(terms,vnf)
 vint=setdiff(terms,vtab)
 vfunc=setdiff(vfunc2,vint)
 vnf=c(vnf2,vint)
 off<-attr(tf,"offset")
 name.coef=nam=par.fregre=beta.l=list()
 kterms=1
 if (attr(tf,"intercept")==0) intercept=FALSE
 else intercept=TRUE
 if (length(vnf)>0) {
 XX=data[[1]][,c(response,vnf2)] #data.frame el 1er elemento de la lista
 for ( i in 1:length(vnf)){
   print(paste("Non functional covariate:",vnf[i]))
   print(paste("The procedure considers only functional covariates and therefore the variable",vnf[i]," is not used."))
     if (kterms > 1)   pf <- paste(pf, "+", vnf[i], sep = "")
     else pf <- paste(pf, vnf[i], sep = "")
     kterms <- kterms + 1
     }
if   (!intercept) {
     pf<- paste(pf,-1,sep="")
     }
}
else {
  # print("peta1");  print(response);  print(names(data))
 XX=data.frame(data[[1]][,response])
 # print("peta2")
 names(XX)=response
}
if (is.null(control$maxit))  control$maxit<-100
if (is.null(control$epsilon))  control$epsilon= 0.001
if (is.null(control$trace))  control$trace = FALSE
if (is.null(control$inverse))  control$inverse = "solve"
############################################################
    xlist<-data[-1] #datos funcionales menos el df!
    y0<-y<-data[["df"]][,response]
if (family$family=="binomial") {
   y<-as.numeric(factor(y,labels=c(0,1)))-1
}
####################    eps <- 0.001
#    if (is.null(control$epsilon)) control$epsilon<-0.00
    eps<-control$epsilon
    namesx<- vfunc
    nvars <- length(vfunc)
    ynames <- if (is.matrix(y))  rownames(y)
              else names(y)
    conv <- FALSE
    nobs <- if (is.matrix(y))   nrow(y)
            else length(y)
    if (is.null(offset))        offset <- rep(0,nobs) ##
    EMPTY <- nvars == 0
#    if (is.null(weights))        weights <- rep(1,nobs)
    variance <- family$variance
    linkinv <- family$linkinv
    linkfun <- family$linkfun
    if (!is.function(variance) || !is.function(linkinv))
        stop("'family' argument seems not to be a valid family object",
            call. = FALSE)
    dev.resids <- family$dev.resids
    aic <- family$aic
    mu.eta <- family$mu.eta
    unless.null <- function(x, if.null) {
        if (is.null(x))             if.null
        else x
    }
    valideta <- unless.null(family$valideta, function(eta) TRUE)
    validmu <- unless.null(family$validmu, function(mu) TRUE)
    X = matrix(0, nrow = nobs, ncol = nvars + intercept)
    colnames(X) = c(namesx, "Intercept")
    eqrank = c(rep(0, nvars), 1)
#    result = vector("list", nvars)
    metric<-metric2<-result<-list()
#    metric2 = metric = vector("list", nvars)
    par.np2=par.np
#    names(eqrank)=names(metric2) = names(metric) = namesx
    names(eqrank)<- c(namesx,"Intercept")
    X[, nvars + intercept] = rep(linkfun(mean(y)), nobs)
    if (control$trace)   cat("----Computing the distance matrix ----\n")
    for (i in 1:nvars) {
#        metric2[[namesx[i]]] = metric.lp(xlist[[namesx[i]]], xlist[[namesx[i]]])
        if (is.null(par.metric)) {
            metric[[namesx[i]]] = metric.lp(xlist[[namesx[i]]], xlist[[namesx[i]]])
        }
        else {
           par.metric[[namesx[i]]]$fdata1 <- xlist[[namesx[i]]]
            metric[[namesx[i]]] = do.call(par.metric[[namesx[i]]]$metric,
                par.metric[[namesx[i]]][-1])
        }
       if (is.null(par.np)) {
          par.np2[[namesx[i]]] =list(Ker=AKer.norm,type.S="S.NW",par.S=list(w=weights))
        }
#       if (is.null(par.np[[namesx[i]]]$par.S)) par.np[[namesx[i]]]$par.S=list(w=weights)
       if (is.null(par.np[[namesx[i]]]$Ker)) par.np2[[namesx[i]]]$Ker=AKer.norm
       if (is.null(par.np[[namesx[i]]]$type.S)) par.np2[[namesx[i]]]$type.S="S.NW"
       if (is.null(par.np[[namesx[i]]]$h)) {
              par.np2[[namesx[i]]]$h = h.default(xlist[[namesx[i]]], len = 51,
              prob = c(0.01,0.66),metric =  metric[[namesx[i]]])
              }
    }
#    eta = apply(X, 1, sum)
     eta = rowSums(X)
    mu = linkinv(eta)
    
    conv <- FALSE
    for (iter in 1L:control$maxit) {
        dev <- devold <- sum(dev.resids(y, mu, weights))
        good <- weights > 0
        varmu <- variance(mu)[good]
        if (any(is.na(varmu))) {stop("NAs in V(mu)")}
        if (any(varmu == 0)) {stop("0s in V(mu)")}
        mu.eta.val <- mu.eta(eta)
        if (any(is.na(mu.eta.val[good])))   stop("NAs in d(mu)/d(eta)")
        good <- (weights > 0) & (mu.eta.val != 0)
        if (all(!good)) {
            conv <- FALSE
            warning("no observations informative at iteration ",iter)
            break
        }

        Xold = X
        ngoodobs <- as.integer(nobs - sum(abs(mu.eta.val) <= eps))
        if (ngoodobs < min(0.2 * nobs, 20)) {
            conv = TRUE
            break
        }
   #  z <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]        #GLM
     ytilde <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
#        ytilde <- eta[good] + (y - mu)[good]/mu.eta.val[good] #ANTES
        w <- sqrt((weights[good] * mu.eta.val[good]^2)/variance(mu)[good])
        X[, nvars + intercept] = rep(mean(ytilde - apply(X[good,1:nvars, drop = FALSE], 1, sum)), nobs)
        if (control$trace)
            cat("#------------------------------------------------\n")
        if (control$trace)
            cat("Iter:", iter, ngoodobs, "/", nobs, "alpha:",
                X[1, nvars + intercept], "Dev:", dev, "\n")
        for (i in 1:nvars) {
            off = apply(X[good, -i, drop = FALSE], 1, sum)
            z = ytilde - off   
            if (control$trace)                print(summary(ytilde))
            offdf = sum(eqrank[-i])
            xfunc = xlist[[namesx[i]]][good]
            mgood<- metric[[namesx[i]]]
            h<-par.np2[[namesx[i]]]$h
             if (control$trace)   cat("Range h:", range(h), length(h), "\n")
           Ker=par.np2[[namesx[i]]]$Ker
           type.S=par.np2[[namesx[i]]]$type.S
           parS=par.np2[[namesx[i]]]$par.S
           parS$w=w
           if (is.function(type.S)) ty<-deparse(substitute(type.S))
           else ty<-type.S
           res = fregre.np.cv(xfunc, z, h = h, type.CV = "dev.S", Ker=Ker,
           type.S=ty,par.S=parS,metric = mgood, par.CV = list(obs = y[good],
           family = family, off = off, offdf = offdf,W = diag(w)))
           if (control$trace)
            cat("Var:",namesx[[i]]," h.opt:", res$h.opt," df:",res$df,"\n")           
           eqrank[namesx[i]] <- res$df
           X[good,namesx[i]] <- res$fitted.values
           result[[namesx[i]]] <- res
        }
        X[, nvars + intercept] = rep(mean(ytilde - apply(X[good,
            1:nvars, drop = FALSE], 1, sum)), nobs)
#        eta <- apply(X, 1, sum)
        eta <- rowSums(X)
        mu <- linkinv(eta <- eta + offset)
#        mu <- linkinv(eta)
        mu.eta.val <- mu.eta(eta)
        good <- (weights > 0) & (mu.eta.val != 0)
        ngoodobs <- as.integer(nobs - sum(mu.eta.val <= eps))
        dev <- sum(dev.resids(y,mu,weights))
        if (control$trace)   {
             par(mfrow = c(1, nvars + 1))
             if (length(table(y) > 10)) {colores = 1}
             else { colores = y + 1 }
             plot(eta, mu, col = colores)
             points(eta, y, col = colores, pch = 2)
             for (i in 1:nvars) {
               plot(X[, i], eta, col = colores, ylab = "Linear Predictor",
                xlab = paste("f(", namesx[i], ")", sep = ""),
                main = paste(namesx[i], "EqPar:", round(eqrank[i],1)))
               if (length(table(y)) == 2) {
                abline(h = 0)
                abline(v = 0)
                }
              }
        }
#        cambio = apply((X - Xold)^2, 2, mean)
        cambio = colMeans((X - Xold)^2)
        if (control$trace) {
            cat("Shift Iter:", iter, "EqRank:", sum(eqrank),
                ngoodobs, "/", nobs, "\n")
             print(cambio)
             }
        if (any(cambio[1:nvars] > control$epsilon)) {conv <- FALSE}
        else{conv <- TRUE;break  }
        if (sum(eqrank) > ngoodobs) {conv <- TRUE;break  }
        if (!(valideta(eta) && validmu(mu))) {
            warning("step size truncated: out of bounds", call. = FALSE)
            ii <- 1
            while (!(valideta(eta) && validmu(mu))) {
                if (ii > control$maxit)
                  stop("inner loop 2; cannot correct step size",
                    call. = FALSE)
                ii <- ii + 1
                X <- (X + Xold)/2
#                eta <- apply(X, 1, sum)
	         eta <- rowSums(X)
                mu <- linkinv(eta)
            }
            boundary <- TRUE
            dev <- sum(dev.resids(y, mu, weights))
            if (control$trace)  cat("Step halved: new deviance =", dev, "\n")
        }
        if (all(!good)) {
            conv <- FALSE
            warning("no observations informative at iteration ",iter)
            break
        }
        if (control$trace)
            cat("Deviance =", dev, "Iterations -", iter, "\n")
        if (abs(dev - devold)/(0.1 + abs(dev)) < control$epsilon |
            dev < control$epsilon) {
            conv <- TRUE
            break
        }
    }
    if (!conv)
        warning("kgam.fit: algorithm did not converge", call. = FALSE)
    if (family$family == "binomial") {
        if (any(mu > 1 - eps) || any(mu < eps))
            warning("kgam.fit: fitted probabilities numerically 0 or 1 occurred",
                call. = FALSE)
    }
    if (family$family == "poisson") {
        if (any(mu < eps))
            warning("kgam.fit: fitted rates numerically 0 occurred",
                call. = FALSE)
    }

    residuals <- (y - mu)#/     (eta) ##### ok??
    nr <- min(sum(good), nvars)
    names(residuals) <- ynames
    names(mu) <- ynames
    names(eta) <- ynames
    wt <- rep.int(0,nobs)
    wt[good] <- w^2
    names(wt) <- ynames
    names(weights) <- ynames
    names(y) <- ynames
    wtdmu <- if (intercept)   sum(weights * y)/sum(weights)
             else linkinv(0)
    nulldev <- sum(dev.resids(y, wtdmu, weights))
    n.ok <- nobs - sum(weights == 0)
    nulldf <- n.ok - as.integer(intercept)
    aic.model <- aic(y, nobs, mu, weights, dev) + 2 * sum(eqrank)
    names(result)<-namesx       
    H<-kgam.H(result,control$inverse)
    
    sr2<-sum(residuals^2)/(nobs - sum(eqrank))
#    if (family$family=="binomial" & !is.factor(y)) y<-factor(y)
    res <- list(result = result, residuals = residuals, fitted.values = mu,
        effects = X, alpha = mean(X[, nvars + intercept]), family = family,
        linear.predictors = eta, deviance = dev, aic = aic.model,
        null.deviance = nulldev, iter = iter, weights = wt, eqrank = eqrank,
        prior.weights = weights, y = y0, converged = conv,H=H,sr2=sr2)
    class(res) <- "fregre.gkam"
    res
}


