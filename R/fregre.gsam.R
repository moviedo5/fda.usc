#' Fitting Functional Generalized Spectral Additive Models
#' 
#' Computes functional GAM model between functional covariate
#' \eqn{(X^1(t_1),\cdots,X^{q}(t_q))}{(X(t_1),...,X(t_q))} (and non functional
#' covariate \eqn{(Z^1,...,Z^p)}{(Z1,...,Zp)}) and scalar response \eqn{Y}.
#' 
#' This function is an extension of the functional generalized linear
#' regression models: \code{\link{fregre.glm}} where the \eqn{E[Y|X,Z]} is
#' related to the linear prediction \eqn{\eta} via a link function
#' \eqn{g(\cdot)}{g(.)} with integrated smoothness estimation by the smooth
#' functions \eqn{f(\cdot)}{f(.)}.
#' 
#' \deqn{E[Y|X,Z])=\eta=g^{-1}(\alpha+\sum_{i=1}^{p}f_{i}(Z^{i})+\sum_{k=1}^{q}\sum_{j=1}^{k_q}{f_{j}^{k}(\xi_j^k)})}{E[Y|X,Z]=\eta=g^{-1}(\alpha+\sum_i
#' f_i(Z_{i})+\sum_k^q\sum_{j=1}^{k_q}{f_j^k(\xi_j^k)})} where
#' \eqn{\xi_j^k}{\xi_j^k} is the coefficient of the basis function expansion of
#' \eqn{X^k}, (in PCA analysis \eqn{\xi_j^k}{\xi_j^k} is the score of the
#' \eqn{j}-functional PC of \eqn{X^k}.
#' 
#' The smooth functions \eqn{f(\cdot)}{f(.)} can be added to the right hand
#' side of the formula to specify that the linear predictor depends on smooth
#' functions of predictors using smooth terms \code{\link{s}} and
#' \code{\link{te}} as in \code{\link{gam}} (or linear functionals of these as
#' \eqn{Z\beta} and \eqn{\big<X(t),\beta\big>}{< X(t),\beta(t) >} in
#' \code{\link{fregre.glm}}).
#' 
#' The first item in the \code{data} list is called \emph{"df"} and is a data
#' frame with the response and non functional explanatory variables, as
#' \code{\link{gam}}.\cr
#' 
#' Functional covariates of class \code{fdata} or \code{fd} are introduced in
#' the following items in the \code{data} list.\cr \code{basis.x} is a list of
#' basis for represent each functional covariate. The basis object can be
#' created by the function: \code{\link{create.pc.basis}}, \code{\link{pca.fd}}
#' \code{\link{create.pc.basis}}, \code{\link{create.fdata.basis}} o
#' \code{\link{create.basis}}.\cr \code{basis.b} is a list of basis for
#' represent each functional beta parameter. If \code{basis.x} is a list of
#' functional principal components basis (see \code{\link{create.pc.basis}} or
#' \code{\link{pca.fd}}) the argument \code{basis.b} is ignored.
#' 
#' @param formula an object of class \code{formula} (or one that can be coerced
#' to that class): a symbolic description of the model to be fitted. The
#' details of model specification are given under \code{Details}.
#' @param family a description of the error distribution and link function to
#' be used in the model. This can be a character string naming a family
#' function, a family function or the result of a call to a family function.
#' (See \code{\link{family}} for details of family functions.)
#' @param data List that containing the variables in the model.
#' @param weights weights
#' @param basis.x List of basis for functional explanatory data estimation.
#' @param basis.b List of basis for functional beta parameter estimation.
#' @param CV =TRUE, Cross-validation (CV) is done.
#' @param \dots Further arguments passed to or from other methods.
#' @return Return \code{gam} object plus:
#' \itemize{
#' \item {basis.x}{ Basis used for \code{fdata} or \code{fd} covariates.} 
#' \item {basis.b}{ Basis used for beta parameter estimation.} 
#' \item {data}{ List that containing the variables in the model.} 
#' \item {formula}{ formula.} 
#' \item {y.pred}{ predicted response by cross-validation.}
#' }
#' @note If the formula only contains a non functional explanatory variables
#' (multivariate covariates), the function compute a standard \code{\link{glm}}
#' procedure.
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' @seealso See Also as: \code{\link{predict.fregre.gsam}} and
#' \code{\link{summary.gam}}.\cr Alternative methods: \code{\link{fregre.glm}}
#' and \code{\link{fregre.gkam}}.
#' @references Muller HG and Stadtmuller U. (2005). \emph{Generalized
#' functional linear models.} Ann. Statist.33 774-805.
#' 
#' Wood (2001) \emph{mgcv:GAMs and Generalized Ridge Regression for R}. R News
#' 1(2):20-25.
#' 
#' Ramsay, James O., and Silverman, Bernard W. (2006), \emph{ Functional Data
#' Analysis}, 2nd ed., Springer, New York.
#' 
#' Venables, W. N. and Ripley, B. D. (2002) \emph{Modern Applied Statistics
#' with S}, New York: Springer.
#' @keywords regression
#' @examples
#' \dontrun{
#' data(tecator)
#' x=tecator$absorp.fdata
#' x.d1<-fdata.deriv(x)
#' tt<-x[["argvals"]]
#' dataf=as.data.frame(tecator$y)
#' nbasis.x=11;nbasis.b=5
#' basis1=create.bspline.basis(rangeval=range(tt),nbasis=nbasis.x)
#' basis2=create.bspline.basis(rangeval=range(tt),nbasis=nbasis.b)
#' f=Fat~s(Protein)+s(x)
#' basis.x=list("x"=basis1,"x.d1"=basis1)
#' basis.b=list("x"=basis2,"x.d1"=basis2)
#' ldata=list("df"=dataf,"x"=x,"x.d1"=x.d1)
#' res=fregre.gsam(Fat~Water+s(Protein)+x+s(x.d1),ldata,family=gaussian(),
#' basis.x=basis.x,basis.b=basis.b)
#' res
#' res2=fregre.gsam(Fat~te(Protein,k=3)+x,data=ldata,family=gaussian())
#' summary(res2)
#' 
#' ##  dropind basis pc
#' basis.pc0=create.pc.basis(x,c(2,4,7))
#' basis.pc1=create.pc.basis(x.d1,c(1:3))
#' basis.x=list("x"=basis.pc0,"x.d1"=basis.pc1)
#' ldata=list("df"=dataf,"x"=x,"x.d1"=x.d1)  
#' res.pc=fregre.gsam(f,data=ldata,family=gaussian(),
#' basis.x=basis.x,basis.b=basis.b)
#' summary(res.pc)
#'  
#' ##  Binomial family
#' x=tecator$absorp.fdata
#' tecator$y$Fat<-ifelse(tecator$y$Fat>20,1,0)
#' x.d1<-fdata.deriv(x)
#' dataf=as.data.frame(tecator$y)
#' ldata=list("df"=dataf,"x"=x,"x.d1"=x.d1)
#' res.bin=fregre.gsam(Fat~Protein+s(x),ldata,family=binomial())
#' }
#' @export
fregre.gsam <- function (formula, family = gaussian(), data = list(), weights = NULL, 
          basis.x = NULL, basis.b = NULL, CV = FALSE, ...) 
{
#print("fregre.gsam")
  tf <- terms.formula(formula, specials = c("s", "te", "t2"))
  terms <- attr(tf, "term.labels")
  special <- attr(tf, "specials")
  nt <- length(terms)
  specials <- rep(NULL, nt)
  if (!is.null(special$s)) 
    specials[special$s - 1] <- "s"
  if (!is.null(special$te)) 
    specials[special$te - 1] <- "te"
  if (!is.null(special$t2)) 
    specials[special$t2 - 1] <- "t2"
  if (attr(tf, "response") > 0) {
    response <- as.character(attr(tf, "variables")[2])
    pf <- rf <- paste(response, "~", sep = "")
  }
  else pf <-rf <- "~"
  vtab <- rownames(attr(tf, "fac"))
  gp <- interpret.gam(formula)
  len.smo <- length(gp$smooth.spec)
  nterms <- length(terms)
  speci <- NULL
  specials1 <- specials2 <- fnf1 <- fnf2 <- fnf <- bs.dim1 <- bs.dim2 <- vfunc2 <- vfunc <- vnf <- NULL
  func <- nf <- sm <- rep(0, nterms)
  names(func) <- names(nf) <- names(sm) <- terms
  ndata <- length(data) - 1
  if (ndata > 0) {
    names.vfunc <- rep("", ndata)
    for (i in 1:ndata) names.vfunc[i] <- names(data)[i + 1]
  }
  else names.vfunc <- NULL
  
#  print(gp)
  ############ nuevo
  covar <- gp$fake.names
 #print(covar)  
 #print(length(covar))
 #print(class(data$df[,covar[1]]))
  if (length(covar)==0) {
    z = gam(formula = gp$pf, data = data$df, family = family)
    z$data <- data
    z$formula.ini <-gp$pf
    z$formula <- gp$pf
    #z$nnf <- nnf
    class(z) <- c(class(z), "fregre.gsam")
    return(z)
  }
  nam.esc<-names(data$df)
  nam.func<-setdiff(names(data),"df")
#  print(nam.esc)    
#  print(nam.func)  
  vnf<-vfunc<-NULL
  for (i in 1:length(covar)){
    if (covar[i] %in% nam.esc) vnf<-c(vnf,covar[i])
    if (covar[i] %in% nam.func) vfunc<-c(vfunc,covar[i])
  }
  nnf<-length(vnf)
  nfunc<-length(vfunc)
  ########################################
  if (len.smo == 0) {
    specials <- rep(NA, nterms)
    speci <- rep("0", nterms)
    gp$smooth.spec[1:nterms] <- NULL
    gp$smooth.spec[1:nterms]$term <- NULL
#    vnf <- terms
    fnf1 <- fnf <- nterms
  }
  else {
    for (i in 1:nterms) if (!is.na(specials[i])) 
      speci <- c(speci, specials[i])
    for (i in 1:nterms) {
      if (any(terms[i] == names(data$df))) {
#        vnf <- c(vnf, terms[i])
        sm[i] <- nf[i] <- 1
        fnf1 <- c(fnf1, 0)
        bs.dim1 <- c(bs.dim1, 0)
        specials1 <- c(specials1, "0")
      }
      else {
        if (any(terms[i] == names.vfunc)) {
#          vfunc <- c(vfunc, terms[i])
          func[i] <- 1
          fnf2 <- c(fnf2, 0)
          bs.dim2 <- c(bs.dim2, 0)
          specials2 <- c(specials2, "0")
        }
      }
    }
    for (i in 1:len.smo) {
      if (speci[i] != "s") {
        if (any(gp$smooth.spec[[i]]$term == names(data$df))) {
#          vnf <- c(vnf, gp$smooth.spec[[i]]$margin[[1]]$term)
          bs.dim1 <- c(bs.dim1, gp$smooth.spec[[i]]$margin[[1]]$bs.dim)
          fnf <- c(fnf, 1)
          fnf1 <- c(fnf1, 1)
          specials1 <- c(specials1, speci[i])
        }
        else {
          if (any(gp$smooth.spec[[i]]$term == names.vfunc)) {
           # vfunc <- c(vfunc, gp$smooth.spec[[i]]$margin[[1]]$term)
            bs.dim2 <- c(bs.dim2, gp$smooth.spec[[i]]$margin[[1]]$bs.dim)
            fnf <- c(fnf, 2)
            fnf2 <- c(fnf2, 1)
            specials2 <- c(specials2, speci[i])
          }
        }
      }
      else {
        if (any(gp$smooth.spec[[i]]$term == names(data$df))) {
       #   vnf <- c(vnf, gp$smooth.spec[[i]]$term)
          bs.dim1 <- c(bs.dim1, gp$smooth.spec[[i]]$bs.dim)
          fnf <- c(fnf, 1)
          fnf1 <- c(fnf1, 1)
          specials1 <- c(specials1, speci[i])
        }
        else {
          if (any(gp$smooth.spec[[i]]$term == names.vfunc)) {
        #    vfunc <- c(vfunc, gp$smooth.spec[[i]]$term)
            bs.dim2 <- c(bs.dim2, gp$smooth.spec[[i]]$bs.dim)
            fnf <- c(fnf, 2)
            fnf2 <- c(fnf2, 1)
            specials2 <- c(specials2, speci[i])
          }
        }
      }
    }
  }
  #nfunc <- sum(func)
  if (is.null(fnf2)) fnf2<-rep(0,length(vfunc))
  
# print(vfunc);print(nnf)

  name.coef = nam = par.fregre = beta.l = list()
  kterms = 1
  if (nnf > 0) {
    XX = data[["df"]]
    if (attr(tf, "intercept") == 0) {
      pf <- paste(pf, -1, sep = "")
    }
    for (i in 1:nnf) {
      if (fnf1[i] == 1 & len.smo != 0) 
        sm1 <- TRUE
      else sm1 <- FALSE
      if (sm1) {
        pf <- paste(pf, "+", specials1[i], "(", vnf[i], 
                    ",k=", bs.dim1[i], ")", sep = "")
      }
      else pf <- paste(pf, "+", vnf[i], sep = "")
      kterms <- kterms + 1
    }
  }
  else {
    XX = data.frame(data[["df"]][, response])
    names(XX) = response
  }
  bsp1 <- bsp2 <- TRUE
  lenfunc <- length(vfunc) > 0
  if (lenfunc) {
    k = 1
    mean.list = vs.list =  list() #JJ =
    bsp1 <- bsp2 <- TRUE
    for (i in 1:length(vfunc)) {
      if (is(data[[vfunc[i]]], "fdata")) {
        tt <- data[[vfunc[i]]][["argvals"]]
        rtt <- data[[vfunc[i]]][["rangeval"]]
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
        else if (class(basis.x[[vfunc[i]]]) == "fdata.comp" | 
                 basis.x[[vfunc[i]]]$type == "pls") 
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
          colnames(J) = colnames(Z) = name.coef[[vfunc[i]]] = paste(vfunc[i], ".", basis.b[[vfunc[i]]]$names, sep = "")
          XX = cbind(XX, Z)
          if (fnf2[i] == 1)             sm2 <- TRUE
          else sm2 <- FALSE
          for (j in 1:length(colnames(Z))) {
            if (sm2) {
              pf <- paste(pf, "+", specials2[i], "(", 
                          name.coef[[vfunc[i]]][j], ",k=", bs.dim2[i], 
                          ")", sep = "")
            }
            else pf <- paste(pf, "+", name.coef[[vfunc[i]]][j], 
                             sep = "")
            kterms <- kterms + 1
          }
          #JJ[[vfunc[i]]] <- J
          vs.list[[vfunc[i]]] <- J
        }
        else {
          l <- basis.x[[vfunc[i]]]$l
          lenl <- length(l)
          vs <- t(basis.x[[vfunc[i]]]$basis$data)
          Z <- basis.x[[vfunc[i]]]$coefs[, l, drop = FALSE]
          response = "y"
          colnames(Z) = name.coef[[vfunc[i]]] = paste(vfunc[i], 
                                                      ".", rownames(basis.x[[vfunc[i]]]$basis$data), 
                                                      sep = "")
          XX = cbind(XX, Z)
          vs.list[[vfunc[i]]] = basis.x[[vfunc[i]]]$basis
          mean.list[[vfunc[i]]] = basis.x[[vfunc[i]]]$mean
          if (fnf2[i] == 1) 
            sm2 <- TRUE
          else sm2 <- FALSE
          for (j in 1:length(colnames(Z))) {
            if (sm2) {
              pf <- paste(pf, "+", specials2[i], "(", 
                          name.coef[[vfunc[i]]][j], ",k=", bs.dim2[i], 
                          ")", sep = "")
            }
            else pf <- paste(pf, "+", name.coef[[vfunc[i]]][j], 
                             sep = "")
            kterms <- kterms + 1
          }
        }
      }
      else {
        if (class(data[[vfunc[i]]])[1] == "fd") {
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
          else if (is(basis.x[[vfunc[i]]],"pca.fd")) 
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
            if (fnf2[i] == 1) 
              sm2 <- TRUE
            else sm2 <- FALSE
            for (j in 1:length(colnames(Z))) {
              if (sm2) {
                pf <- paste(pf, "+", specials2[i], "(", 
                            name.coef[[vfunc[i]]][j], ",k=", bs.dim2[i], 
                            ")", sep = "")
              }
              else pf <- paste(pf, "+", name.coef[[vfunc[i]]][j], 
                               sep = "")
              kterms <- kterms + 1
            }
            ## JJ[[vfunc[i]]] <- J
            vs.list[[vfunc[i]]] <- J
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
            if (fnf2[i] == 1) 
              sm2 <- TRUE
            else sm2 <- FALSE
            for (j in 1:length(colnames(Z))) {
              if (sm2) {
                pf <- paste(pf, "+", specials2[i], "(", 
                            name.coef[[vfunc[i]]][j], ",k=", bs.dim2[i], 
                            ")", sep = "")
              }
              else pf <- paste(pf, "+", name.coef[[vfunc[i]]][j], 
                               sep = "")
              kterms <- kterms + 1
            }
          }
        }
        else stop("Please, enter functional covariate")
      }
    }
  }
  #print("fregre.gsam 3")
  if (!is.data.frame(XX)) 
    XX = data.frame(XX)
  par.fregre$data = XX
  ndatos <- nrow(XX)
  yp <- rep(NA, ndatos)
  if (length(vfunc) == 0 & length(vnf) == 0) {
    pf <- as.formula(paste(pf, 1, sep = ""))
# print("gam gam gam ")    
    z = gam(formula = pf, data = XX, family = family, 
            weights = weights,...)
# print("gam gam gam NOOO")        
    z$XX = XX
    z$data <- data
    z$formula.ini <- pf
    z$formula <- pf
    class(z) <- c(class(z), "fregre.gsam")
    return(z)
  }
  else par.fregre$formula = as.formula(pf)
  z = gam(formula = par.fregre$formula, data = XX, family = family, 
          weights = weights, ...)
  if (lenfunc) {
    z$mean = mean.list
    z$basis.x = basis.x
    z$basis.b = basis.b
    #z$JJ <- JJ
    z$vs.list = vs.list
    z$vfunc <- vfunc
  }
  z$formula = as.formula(pf)
  z$formula.ini = formula
  z$data = z$data
  z$XX = XX
  z$nnf <- nnf
  class(z) <- c("fregre.gsam",class(z))
  z
}