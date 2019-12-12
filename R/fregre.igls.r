#' @title Fit of  Functional Generalized Least Squares Model Iteratively 
#' @aliases fregre.igls
#' 
#' @param formula A two-sided linear formula object describing the
#' model, with the response on the left of a \code{~} operator and the
#' terms, separated by \code{+} operators, on the right.
#' @param data An optional data frame containing the variables named in
#' \code{model}, \code{correlation}, \code{weights}, and
#' \code{subset}. By default the variables are taken from the environment
#'  from which \code{gls} is called.
#' @param basis.x List of basis for functional explanatory data estimation.
#' @param basis.b List of basis for \eqn{\beta(t)} parameter estimation.
#' @param rn List of Ridge parameter.
#' @param lambda List of Roughness penalty parameter.
#' @param correlation an optional \code{\link{corStruct}} object describing the
#' within-group correlation structure. See the documentation of
#' \code{\link{corClasses}} for a description of the available \code{corStruct}
#' classes. If a grouping variable is to be used, it must be specified in
#' the \code{form} argument to the \code{corStruct} constructor. Defaults to 
#' \code{NULL}, corresponding to uncorrelated errors.
#' @param maxit Number of maximum of interactions.
#' @param weights An optional \code{\link{varFunc}} object or one-sided formula
#'  describing the within-group heteroscedasticity structure. If given as
#'  a formula, it is used as the argument to \code{\link{varFixed}},
#'  corresponding to fixed variance weights. See the documentation on
#'  \code{\link{varClasses}} for a description of the available \code{\link{varFunc}}
#'  classes. Defaults to \code{NULL}, corresponding to homoscedastic errors.
#' @param control Control parameters.  
#' @param \dots Further arguments passed to or from other methods.

#' @description  This function fits iteratively a functional linear model using generalized 
#' least squares. The errors are allowed to be correlated and/or have unequal variances.  
#' \enumerate{
#'  \item Begin with a preliminary estimation of \eqn{\hat{\theta}=\theta_0} (for instance, 
#'  \eqn{\theta_0=0}). Compute \eqn{\hat{W}}.
#' \item Estimate \eqn{b_\Sigma =(Z'\hat{W}Z)^{-1}Z'\hat{W}y}
#' \item Based on the residuals, \eqn{\hat{e}=\left(y-Zb_\Sigma \right)}, update
#'  \eqn{\hat{\theta}=\rho\left({\hat{e}}\right)} where \eqn{\rho} depends on the 
#'  dependence structure chosen.
#' \item Repeats steps 2 and 3 until convergence (small changes in \eqn{b_\Sigma} and/or \eqn{\hat{\theta}}). 
#' }
#' @return An object of class \code{"gls"} representing the functional linear model
#' fit. Generic functions such as \code{print}, \code{plot}, and \code{summary} have
#' methods to show the results of the fit. 
#' 
#'  See \code{\link{glsObject}} for the components of the fit. The functions
#'  \code{\link{resid}}, \code{\link{coef}} and \code{\link{fitted}}, can be used to
#'   extract some of its components. 
#' Beside, the class(z) is  "gls", "lm" and "fregre.lm" with the following objects:
#' \itemize{
#'  \item{sr2}{ Residual variance.}
#'  \item{Vp}{ Estimated covariance matrix for the parameters.}
#'  \item{lambda}{ A roughness penalty.}	
#'  \item{basis.x}{ Basis used for \code{fdata} or \code{fd} covariates.}
#'  \item{basis.b}{ Basis used for beta parameter estimation.}
#'  \item{beta.l}{ List of estimated beta parameter of functional covariates.}
#'  \item{data}{ List that containing the variables in the model.}
#'  \item{formula}{ formula used in ajusted model.}
#'  \item{formula.ini}{ formula in call.}
#'  \item{XX}{ desing matrix }
#'  \item{W}{ inverse of covariance matrix}
#'  \item{fdataob}{ }
#'  \item{rn}{ rn}
#'  \item{vs.list}{ }
#'  \item{correlation}{ See glsObject for the components of the fit. }
#'  }
#' @references  Oviedo de la Fuente, M., Febrero-Bande, M., Pilar Munoz,
#' and Dominguez, A. Predicting seasonal influenza transmission using Functional 
#' Regression Models with Temporal Dependence. arXiv:1610.08718.
#'  \url{https://arxiv.org/abs/1610.08718}
#' @examples
#' \dontrun{ 
#' data(tecator)
#' x=tecator$absorp.fdata
#' x.d2<-fdata.deriv(x,nderiv=)
#' tt<-x[["argvals"]]
#' dataf=as.data.frame(tecator$y)
#' # plot the response
#' plot(ts(tecator$y$Fat))
#' ldata=list("df"=dataf,"x.d2"=x.d2)
#' res.gls=fregre.igls(Fat~x.d2,data=ldata,
#' correlation=list("cor.ARMA"=list()),
#' control=list("p"=1)) 
#' res.gls
#' res.gls$corStruct
#' }
#' @keywords models regression
#' @export
fregre.igls<-function (formula, data, basis.x = NULL, basis.b = NULL, correlation, 
          maxit = 100, rn, lambda, weights = rep(1, n), control, ...) 
{
#print("fregre.igls")
  tf <- terms.formula(formula)
  terms <- attr(tf, "term.labels")
  nt <- length(terms)
  if (attr(tf, "response") > 0) {
    response <- as.character(attr(tf, "variables")[2])
    pf <- rf <- paste(response, "~", sep = "")
  }
  else pf <- rf <- "~"
  vtab <- rownames(attr(tf, "factors"))
  vnf = intersect(terms, names(data$df))
  vnf2 = intersect(vtab[-1], names(data$df)[-1])
  vfunc2 = setdiff(terms, vnf)
  vint = setdiff(terms, vtab)
  vfunc = setdiff(vfunc2, vint)
  off <- attr(tf, "offset")
  name.coef = nam = par.fregre = beta.l = list()
  kterms = 1
  n <- length(data[["df"]][, response])
  XX = data.frame(data[["df"]][, c(response)], weights)
  namxx = names(XX) = c(response, "weights")
  aa <- NULL
  if (length(vnf) > 0) {
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
    pf <- rf <- paste(response, "~", sep = "")
    for (i in 1:length(vnf2)) pf <- paste(pf, "+", vnf2[i], 
                                          sep = "")
    cname <- names(XX)
    XX <- data.frame(XX, mf)
    names(XX) <- c(cname, names(mf))
  }
  else {
    XX <- data.frame(XX, model.matrix(formula(paste(pf, "1")), 
                                      data$df))
    names(XX)[3] <- "(Intercept)"
  }
  if (missing(rn)) {
    rn0 = FALSE
    rn = list()
  }
  else rn0 <- TRUE
  if (missing(lambda)) {
    lambda0 = FALSE
    lambda = list()
  }
  else lambda0 <- TRUE
  mat <- rep(0, len = ncol(XX) - 2)
  imat2 <- ncol(XX) - 2
  mat2 <- diag(0, nrow = imat2)
  if (length(vfunc) > 0) {
    mean.list = vs.list = JJ = list()
    bsp1 <- bsp2 <- TRUE
    for (i in 1:length(vfunc)) {
      if (is(data[[vfunc[i]]],"fdata")) {
        tt <- data[[vfunc[i]]][["argvals"]]
        rtt <- data[[vfunc[i]]][["rangeval"]]
        np <- length(tt)
        fdat <- data[[vfunc[i]]]
        dat <- data[[vfunc[i]]]$data
        if (is.null(basis.x[[vfunc[i]]])) 
          basis.x[[vfunc[i]]] <- create.fdata.basis(fdat, 
                                                    l = 1:7)
        else if (basis.x[[vfunc[i]]]$type == "pc" | basis.x[[vfunc[i]]]$type == 
                 "pls") 
          bsp1 = FALSE
        if (is.null(basis.b[[vfunc[i]]]) & bsp1) 
          basis.b[[vfunc[i]]] <- create.fdata.basis(fdat)
        else if (basis.x[[vfunc[i]]]$type == "pc" | basis.x[[vfunc[i]]]$type == 
                 "pls") 
          bsp2 = FALSE
        if (bsp1 & bsp2) {
          if (is.null(rownames(dat))) 
            rownames(fdat$data) <- 1:nrow(dat)
          fdnames = list(time = tt, reps = rownames(fdat[["data"]]), 
                         values = "values")
          xcc <- fdata.cen(data[[vfunc[i]]])
          mean.list[[vfunc[i]]] = xcc[[2]]
          if (!is.null(basis.x[[vfunc[i]]]$dropind)) {
            int <- setdiff(1:basis.x[[vfunc[i]]]$nbasis, 
                           basis.x[[vfunc[i]]]$dropind)
            basis.x[[vfunc[i]]]$nbasis <- length(int)
            basis.x[[vfunc[i]]]$dropind <- NULL
            basis.x[[vfunc[i]]]$names <- basis.x[[vfunc[i]]]$names[int]
          }
          if (!is.null(basis.b[[vfunc[i]]]$dropind)) {
            int <- setdiff(1:basis.b[[vfunc[i]]]$nbasis, 
                           basis.b[[vfunc[i]]]$dropind)
            basis.b[[vfunc[i]]]$nbasis <- length(int)
            basis.b[[vfunc[i]]]$dropind <- NULL
            basis.b[[vfunc[i]]]$names <- basis.b[[vfunc[i]]]$names[int]
          }
          x.fd = Data2fd(argvals = tt, y = t(xcc[[1]]$data), 
                         basisobj = basis.x[[vfunc[i]]], fdnames = fdnames)
          r = x.fd[[2]][[3]]
          J = inprod(basis.x[[vfunc[i]]], basis.b[[vfunc[i]]])
          Z = t(x.fd$coefs) %*% J
          colnames(J) = colnames(Z) = name.coef[[vfunc[i]]] = paste(vfunc[i], 
                                                                    ".", basis.b[[vfunc[i]]]$names, sep = "")
          XX = cbind(XX, Z)
          for (j in 1:length(colnames(Z))) {
            if (kterms >= 1) 
              pf <- paste(pf, "+", colnames(Z)[j], sep = "")
            else pf <- paste(pf, colnames(Z)[j], sep = "")
            kterms <- kterms + 1
          }
          JJ[[vfunc[i]]] <- J
          nbasisb <- basis.b[[vfunc[i]]]$nbasis
          R = diag(0, ncol = nbasisb, nrow = nbasisb)
          R <- eval.penalty(basis.b[[vfunc[i]]], Lfdobj = int2Lfd(2))
          lenm <- ncol(mat2)
          if (rn0) 
            stop("Ridge regressions is only implemented for functional principal component basis")
          MM <- matrix(0, nrow = lenm, ncol = nbasisb)
          MM2 <- matrix(0, nrow = nbasisb, ncol = imat2 + 
                          nbasisb)
             

          mat2 <- cbind(mat2, MM)
          mat2 <- rbind(mat2, MM2)
          if (!is.null(lambda[[vfunc[i]]])) {
            mat2[(imat2 + 1):(imat2 + nbasisb), (imat2 + 
                                                   1):(imat2 + nbasisb)] <- lambda[[vfunc[i]]] * 
              R
          }
          imat2 <- imat2 + nbasisb
        }
        else {
          basis <- basis.x[[vfunc[i]]]
          l <- basis$l
          vs <- t(basis$basis$data)
          basis$x <- basis$x[, l, drop = FALSE]
          Z <- basis$x
          response = "y"
          colnames(Z) = name.coef[[vfunc[i]]] = paste(vfunc[i], 
                                                      ".", rownames(basis$basis$data), sep = "")
          XX = cbind(XX, Z)
          vs.list[[vfunc[i]]] = basis$basis
          mean.list[[vfunc[i]]] = basis$mean
          for (j in 1:length(colnames(Z))) {
            if (kterms >= 1) 
              pf <- paste(pf, "+", name.coef[[vfunc[i]]][j], 
                          sep = "")
            else pf <- paste(pf, name.coef[[vfunc[i]]][j], 
                             sep = "")
            kterms <- kterms + 1
          }
          if (!is.null(rn[[vfunc[i]]])) {
            mat <- c(mat, rep(rn[[vfunc[i]]], len = length(l)))
          }
          else mat <- c(mat, rep(0, len = length(l)))
          lenl <- length(l)
          lenm <- ncol(mat2)
          MM <- matrix(0, nrow = lenm, ncol = lenl)
          MM2 <- matrix(0, nrow = lenl, ncol = imat2 + lenl)
          mat2 <- cbind(mat2, MM)
          mat2 <- rbind(mat2, MM2)
          if (!is.null(lambda[[vfunc[i]]])) {
            R <- P.penalty(1:length(l))
            mat2[(imat2 + 1):(imat2 + lenl), (imat2 + 
                                                1):(imat2 + lenl)] <- lambda[[vfunc[i]]] * 
              R
            
          }
          imat2 <- imat2 + lenl
        }
      }
      else {
        if (is(data[[vfunc[i]]],"fd")) {
          fdat <- data[[vfunc[i]]]
          if (is.null(basis.x[[vfunc[i]]])) 
            basis.x[[vfunc[i]]] <- fdat$basis
          else if (class(basis.x[[vfunc[i]]]) == "pca.fd") 
            bsp1 = FALSE
          if (is.null(basis.b[[vfunc[i]]]) & bsp1) 
            basis.b[[vfunc[i]]] <- create.fdata.basis(fdat, 
                                                      l = 1:max(5, floor(basis.x[[vfunc[i]]]$nbasis/5)), 
                                                      type.basis = basis.x[[vfunc[i]]]$type, 
                                                      rangeval = fdat$basis$rangeval)
          else if (class(basis.x[[vfunc[i]]]) == "pca.fd") 
            bsp2 = FALSE
          if (bsp1 & bsp2) {
            r = fdat[[2]][[3]]
            if (!is.null(basis.x[[vfunc[i]]]$dropind)) {
              int <- setdiff(1:basis.x[[vfunc[i]]]$nbasis, 
                             basis.x[[vfunc[i]]]$dropind)
              basis.x[[vfunc[i]]]$nbasis <- length(int)
              basis.x[[vfunc[i]]]$dropind <- NULL
              basis.x[[vfunc[i]]]$names <- basis.x[[vfunc[i]]]$names[int]
            }
            if (!is.null(basis.b[[vfunc[i]]]$dropind)) {
              int <- setdiff(1:basis.b[[vfunc[i]]]$nbasis, 
                             basis.b[[vfunc[i]]]$dropind)
              basis.b[[vfunc[i]]]$nbasis <- length(int)
              basis.b[[vfunc[i]]]$dropind <- NULL
              basis.b[[vfunc[i]]]$names <- basis.b[[vfunc[i]]]$names[int]
            }
            J = inprod(basis.x[[vfunc[i]]], basis.b[[vfunc[i]]])
            mean.list[[vfunc[i]]] <- mean.fd(x.fd)
            x.fd <- center.fd(x.fd)
            Z = t(x.fd$coefs) %*% J
            colnames(J) = colnames(Z) = name.coef[[vfunc[i]]] = paste(vfunc[i], 
                                                                      ".", basis.b[[vfunc[i]]]$names, sep = "")
            XX = cbind(XX, Z)
            for (j in 1:length(colnames(Z))) {
              if (kterms >= 1) 
                pf <- paste(pf, "+", colnames(Z)[j], 
                            sep = "")
              else pf <- paste(pf, colnames(Z)[j], sep = "")
              kterms <- kterms + 1
            }
            JJ[[vfunc[i]]] <- J
          }
          else {
            l <- ncol(basis.x[[vfunc[i]]]$scores)
            vs <- basis.x[[vfunc[i]]]$harmonics$coefs
            Z <- basis.x[[vfunc[i]]]$scores
            response = "y"
            colnames(Z) = name.coef[[vfunc[i]]] = paste(vfunc[i], 
                                                        ".", colnames(basis.x[[vfunc[i]]]$harmonics$coefs), 
                                                        sep = "")
            XX = cbind(XX, Z)
            vs.list[[vfunc[i]]] = vs
            mean.list[[vfunc[i]]] = basis.x[[vfunc[i]]]$meanfd
            for (j in 1:length(colnames(Z))) {
              if (kterms >= 1) 
                pf <- paste(pf, "+", name.coef[[vfunc[i]]][j], 
                            sep = "")
              else pf <- paste(pf, name.coef[[vfunc[i]]][j], 
                               sep = "")
              kterms <- kterms + 1
            }
          }
        }
        else stop("Please, enter functional covariate")
      }
    }
  }
  if (!is.data.frame(XX)) 
    XX = data.frame(XX)
  par.fregre$formula = pf
  par.fregre$data = XX
  y <- XX[, 1]
  scores <- as.matrix(XX[, -(1:2)])
  W <- diag(weights)
  error <- FALSE
 
  
  if (missing(correlation)) 
    correlation0 <- FALSE
  else {
    if (is.matrix(correlation)) 
      error = FALSE
    else error = TRUE
    correlation0 <- TRUE
  }
 if (!rn0 & !lambda0 & !correlation0) {
  #  print(1)
    W0 <- FALSE
    z = lm(formula = pf, data = XX, ...)
    e <- z$residuals
    S <- solve(t(scores) %*% W %*% scores)
    class(z) <- c(class(z), "fregre.lm")
    mat0 <- diag(0, ncol(scores))
  }
  else {
 #   print(11)
    mat0 <- diag(0, ncol(scores))
    if (lambda0) {
      mat0 <- mat2
    }
    if (rn0) {
      mat0 <- diag(mat)
    }
    if (rn0 & lambda0) 
      warning("Only ridge penalization is done by rn argument (lambda argument is ignored)")
    W <- diag(weights) %*% W
    S <- solve(t(scores) %*% W %*% scores + mat0)
    Cinv <- S %*% t(scores) %*% W
    ycen = y - mean(y)
    coefs <- Cinv %*% XX[, 1]
    z <- list()
    z$fitted.values <- yp <- drop(scores %*% coefs)
    e <- z$residuals <- XX[, 1] - z$fitted.values
    it <- 1
    eps = 0.001
    err2 = sqrt(sum(coefs^2))
    MM <- matrix(1:n, ncol = 1)
    corStruct <- list()
 #   print(error)
    while (error) {
      W <- W0 <- diag(n)
      if (is.list(correlation)) {
        ncor <- length(correlation)
        name.cor <- names(correlation)
      }
      else {
        ncor <- 1
        name.cor <- correlation[[1]]
        correlation = list(correlation = list())
        names(correlation) <- name.cor
      }
      wei0 <- wei <- diag(weights)
      ee <- e
       for (icor in 1:ncor) {
    #    print(icor);        print(correlation)
        if (!is.list(correlation[[icor]])) {
          name.cor <- correlation[[icor]]
          correlation[[icor]] = list()
          names(correlation)[icor] <- name.cor
        }
        if (is.null(correlation[[icor]][["group"]])) {
          n <- length(e)
  #  print(name.cor)
  #  print(name.cor[icor])
          par.cor <- switch(name.cor[icor], cor.ARMA = {
            cor.name <- "cor.ARMA"
     #print("cor.ARMA")
            if (is.null(correlation[[icor]][["index"]])) index <- 1:n else {
              index.df <- correlation[[icor]][["index"]]
              index <- data[["df"]][, index.df]
              if (is.vector(index)) e <- e[order(index)] else {
                par.ord <- list()
                for (i in 1:ncol(index)) par.ord[[i]] <- index[, 
                                                               i]
                index <- as.numeric(do.call("order", 
                                            par.ord))
                e <- e[index]
              }
            }
            list(x = e[index])
          }, cor.AR = {
            cor.name <- "cor.AR"
            #    print("cor.AR")
            if (is.null(correlation[[icor]][["index"]])) index <- 1:n else {
              index.df <- correlation[[icor]][["index"]]
              index <- data[["df"]][, index.df]
              if (is.vector(index)) e <- e[order(index)] 
              else {
                par.ord <- list()
                for (i in 1:ncol(index)) par.ord[[i]] <- index[,i]
                index <- as.numeric(do.call("order", par.ord))
                e <- e[index]
              }
            }
            list(x = e[index])
          }, corExpo = {
            index.df <- correlation[[icor]][["index"]]
            dxy <- data[["df"]][, index.df]
            index <- 1:length(e)
            cor.name <- "corExpo"
            list(xy = dxy)
          }, corVgm = {
            index.df <- correlation[[icor]][["index"]]
            xy <- data[["df"]][, index.df]
            dxy0 <- as.matrix(dist(xy, diag = TRUE, upper = TRUE))
            index <- 1:n
            cor.name <- "corVgm"
            residual <- e[index]
            list(dxy = dxy0, xy = xy, df = data.frame(residual))
          }, cor.Exp = {
            index.df <- correlation[[icor]][["index"]]
            xy <- data[["df"]][, index.df]
            index <- 1:n
            cor.name <- "cor.Exp"
            residual <- e[index]
            dff <- data.frame(data[["df"]][, index.df], 
                              residual)
            names(dff) <- c(index.df, "residual")
            list(df = dff)
          }, corUnstruc = {
            cor.name <- "corUnstruc"
            list(correlation[[icor]][["index"]])
          })
          e <- e[index]
# print(names(par.cor))
          # print(names(correlation))
          # print("1")          
          if (!missing(control)) {
            npar <- length(par.cor)
            #print(names(par.cor))
            #print(names(par.ord))
            # print("2")
            for (i in 1:length(control)) {
              par.cor[npar + i] <- control[i]
              names(par.cor)[npar + i] <- names(control)[i]
            }
          }
          #    print("antes do.call")    
          #  print(cor.name)
          #  print(names(par.cor))
          aa <- do.call(cor.name, par.cor)
          W0 <- aa[[1]]
          corStruct <- aa[-1]
        }
        else {
#print("hay grupo")          
          n <- nrow(data$df)
          W00 <- matrix(0, ncol = n, nrow = n)
          name.group <- correlation[[icor]]$group
          ncomb <- length(name.group)
#     print(name.group)          
          
          if (ncomb > 1) {
            gr <- (data[["df"]][, name.group[1]])
            for (ngr in 2:ncomb) gr <- paste(gr, (data[["df"]][, 
                                                               name.group[ngr]]), sep = "")
            gr <- factor(gr, levels = unique(gr))
          }
          #else gr <- data[["df"]][[name.group]]
          else gr <- data[["df"]][,name.group]
          if (!is.factor(gr)) 
            gr <- factor(gr)
          lev <- levels(gr)
          nlev <- length(lev)
          #print(name.cor[icor])          
          par.cor <- switch(name.cor[icor], 
            
          cor.ARMA = {
#print("entra cor.ARMAAAAAA")
            index.df <- correlation[[icor]][["index"]]
            #print(correlation)
            #print(index.df)            
            index <- data[["df"]][, index.df]
            e2 <- NULL
           
            for (i in 1:nlev) {
              jj <- gr == lev[i]
              e0 <- e[jj]
              e2 <- cbind(e2, e0[order(index[jj])])
            }
            cor.name <- "cor.ARMA"
            #print(correlation)
            #print("okk")
            if (is.null(correlation[[icor]][["p"]])) p0<-1
            else    p0<-correlation[[icor]][["p"]]
            #print(p0)
            list(x = e2,p=p0,order.max=p0)
          }, 
          corCloud = {
            index.df <- correlation[[icor]][["index"]]
            index <- data[["df"]][, index.df]
            gr <- data[["df"]][, correlation[[icor]][["group"]]]
            xy <- data[["df"]][, index.df]
            cor.name <- "corCloud"
            residual <- e
            dff <- data.frame(data[["df"]][, index.df], 
                              residual)
            names(dff) <- c(index.df, "residual")
            list(df = dff, gr = gr)
          }, 
          cor.Exp = {
            index.df <- correlation[[icor]][["index"]]
            index <- data[["df"]][, index.df]
            gr <- data[["df"]][, correlation[[icor]][["group"]]]
            xy <- data[["df"]][, index.df]
            cor.name <- "cor.Exp"
            residual <- e
            dff <- data.frame(data[["df"]][, index.df], 
                              residual)
            names(dff) <- c(index.df, "residual")
            list(df = dff, gr = gr)
          }, 
        corUnstruc = {
            cor.name <- "corUnstruc"
            index <- 1:nrow(dff)
            par.cor <- list(correlation[[icor]][["index"]])
          })
          if (!missing(control)) {
            npar <- length(par.cor)
            for (k in 1:length(control)) {
              par.cor[npar + k] <- control[k]
              names(par.cor)[npar + k] <- names(control)[k]
            }
          }
          #print("calll")
          #print(cor.name)
          #print(par.cor)
          #print("antes call")
          aa <- do.call(cor.name, par.cor)
          W0 <- aa[[1]]
          corStruct <- aa[-1]
          for (j in 1:nlev) {
            n <- table(gr)[j]
            ilev <- gr == lev[j]
            e <- ee[ilev]
            if (name.cor[icor] == "cor.Exp" | name.cor[icor] == 
                "corCloud") 
              W00[ilev, ilev] <- W0
            else W00[ilev, ilev] <- W0[order(index[ilev]), 
                                       order(index[ilev])]
          }
          W0 <- W00
        }
        wei0 <- wei0 %*% W0
      }
      W0 <- wei0
      W <- try(solve(W0), silent = TRUE)
#     print(class(W) )
      if (is(W,"try-error")) {
        sv <- svd(W0)
        W <- drop((sv$v %*% diag(1/sv$d) %*% t(sv$u)))
      }
      mat0[1, 1] <- 0
      W2 <- t(scores) %*% W %*% scores + mat0
      S <- try(solve(W2), silent = TRUE)
      if (is(S,"try-error")) {
        sv <- svd(W2)
        S <- drop((sv$v %*% diag(1/sv$d) %*% t(sv$u)))
      }
      Cinv <- S %*% t(scores) %*% W
      coefs <- Cinv %*% matrix(y, ncol = 1)
      yp <- drop(scores %*% coefs)
      coefs <- drop(coefs)
      e <- y - yp
      err3 <- coefs
      err4 <- max(abs((err3 - err2)/err2))
      if (err4 < eps) 
        error <- FALSE
      else {
        if (it == maxit) 
          error <- FALSE
        else it <- it + 1
        err2 <- err3
      }
    }
    H <- scores %*% Cinv
    df <- traza(H)
    coefs <- drop(coefs)
    z$coefficients <- coefs
    z$mean.list <- mean.list
    z$df.residual <- n - df
    z$H <- H
    z$r2 <- 1 - sum(z$residuals^2)/sum(ycen^2)
    if (is(basis.x[[vfunc[1]]],"basisfd")) {
      z$call[[1]] = "fregre.basis"
    }
    else {
      if (basis.x[[vfunc[1]]]$type == "pc") 
        z$call[[1]] = "fregre.pc"
      if (basis.x[[vfunc[1]]]$type == "pls") 
        z$call[[1]] = "fregre.pls"
    }
    class(z) <- c("fregre.fd", "fregre.lm")
    rdf <- n - df
    sr2 <- sum(e^2)/rdf
    r2 <- 1 - sum(e^2)/sum(ycen^2)
    r2.adj <- 1 - (1 - r2) * ((n - 1)/rdf)
    GCV <- sum(e^2)/(n - df)^2
    z$terms <- terms
    z$residuals <- drop(e)
    z$fitted.values <- yp
    z$y <- y
    z$rank <- df
    z$df.residual <- rdf
    Z = cbind(rep(1, len = n), Z)
    colnames(Z)[1] = "(Intercept)"
    std.error = sqrt(diag(S) * sr2)
    t.value = coefs/std.error
    p.value = 2 * pt(abs(t.value), n - df, lower.tail = FALSE)
    coefficients <- cbind(coefs, std.error, t.value, p.value)
    colnames(coefficients) <- c("Estimate", "Std. Error", 
                                "t value", "Pr(>|t|)")
    z$coefs <- coefficients
    class(z) <- "lm"
    z$it <- it
  }
  for (i in 1:length(vfunc)) {
    if (bsp1) 
      beta.l[[vfunc[i]]] = fd(z$coefficients[name.coef[[vfunc[i]]]], 
                              basis.b[[vfunc[i]]])
    else {
      if (class(data[[vfunc[i]]])[1] == "fdata") {
        beta.est <- z$coefficients[name.coef[[vfunc[i]]]] * 
          vs.list[[vfunc[i]]]
        beta.est$data <- colSums(beta.est$data)
        beta.est$names$main <- "beta.est"
        beta.est$data <- matrix(as.numeric(beta.est$data), 
                                nrow = 1)
        beta.est$names$main <- "beta.est"
        beta.est$data <- matrix(as.numeric(beta.est$data), 
                                nrow = 1)
        if (basis.x[[vfunc[i]]]$type == "pls") {
          if (basis.x[[vfunc[i]]]$norm) {
            sd.X <- sqrt(apply(data[[vfunc[i]]]$data, 
                               2, var))
            beta.est$data <- beta.est$data/sd.X
          }
        }
        beta.l[[vfunc[i]]] <- beta.est
      }
      else {
        beta.est <- z$coefficients[name.coef[[vfunc[i]]]] * 
          t(vs.list[[vfunc[i]]])
        beta.est <- colSums(beta.est)
        beta.l[[vfunc[i]]] <- fd(beta.est, basis.x[[vfunc[i]]]$harmonics$basis)
      }
    }
  }
  z$sr2 <- sum(e^2)/z$df.residual
  z$Vp = z$sr2 * S
  z$beta.l = beta.l
  z$formula = pf
  z$mean = mean.list
  z$formula.ini = formula
  z$basis.x = basis.x
  z$basis.b = basis.b
  z$JJ <- JJ
  z$data = z$data
  z$XX = XX
  z$data <- data
  z$fdataobj <- data[[vfunc[1]]]
  if (correlation0) {
    rn0 <- TRUE
    z$corStruct <- corStruct
    z$it <- it
  }
  z$rn <- rn0
  z$lambda <- lambda0
  z$W <- W
  z$W0 <- W0
  z$correlation <- correlation
  z$correl <- correlation0
  z$vs.list = vs.list
  class(z) <- c( "fregre.igls","lm")
  z
}


