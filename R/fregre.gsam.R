#' @title Fitting Functional Generalized Spectral Additive Models
#' 
#' @description Computes a functional GAM model between a functional covariate
#' \eqn{(X^1(t_1), \dots, X^{q}(t_q))}{(X(t_1), ..., X(t_q))} and a non-functional
#' covariate \eqn{(Z^1, ..., Z^p)}{(Z1, ..., Zp)} with a scalar response \eqn{Y}.
#' 
#' This function extends functional generalized linear regression models (\code{\link{fregre.glm}}) 
#' where \eqn{E[Y|X,Z]} is related to the linear predictor \eqn{\eta} via a link function 
#' \eqn{g(\cdot)}{g(.)} with integrated smoothness estimation by the smooth
#' functions \eqn{f(\cdot)}{f(.)}.
#' 
#' \deqn{E[Y|X,Z]=\eta=g^{-1}\left(\alpha+\sum_{i=1}^{p}f_{i}(Z^{i})+\sum_{k=1}^{q}\sum_{j=1}^{k_q} f_{j}^{k}(\xi_j^k)\right)}{E[Y|X,Z]=\eta=g^{-1}(\alpha+\sum_i f_i(Z_{i})+\sum_k^q \sum_{j=1}^{k_q} f_j^k(\xi_j^k))}
#' 
#' where \eqn{\xi_j^k}{\xi_j^k} is the coefficient of the basis function expansion of
#' \eqn{X^k}; in PCA analysis, \eqn{\xi_j^k}{\xi_j^k} is the score of the
#' \eqn{j}-functional PC of \eqn{X^k}.
#' 
#' @details The smooth functions \eqn{f(\cdot)}{f(.)} can be added to the right-hand
#' side of the formula to specify that the linear predictor depends on smooth
#' functions of predictors using smooth terms \code{\link[mgcv]{s}} and
#' \code{\link[mgcv]{te}} as in \code{\link[mgcv]{gam}} (or linear functionals of these as
#' \eqn{Z \beta} and \eqn{\langle X(t), \beta \rangle}{< X(t),\beta(t) >} in
#' \code{\link{fregre.glm}}).
#' 
#' The first item in the \code{data} list is called \emph{"df"} and is a data
#' frame with the response and non-functional explanatory variables, as in
#' \code{\link[mgcv]{gam}}.\cr
#' 
#' Functional covariates of class \code{fdata} or \code{fd} are introduced in
#' the following items of the \code{data} list.\cr \code{basis.x} is a list of
#' basis functions for representing each functional covariate. The basis object can be
#' created using functions such as \code{\link{create.pc.basis}}, \code{\link[fda]{pca.fd}},
#' \code{\link{create.fdata.basis}}, or \code{\link[fda]{create.basis}}.\cr 
#' \code{basis.b} is a list of basis functions for representing each functional beta parameter. 
#' If \code{basis.x} is a list of functional principal components basis functions 
#' (see \code{\link{create.pc.basis}} or \code{\link[fda]{pca.fd}}), the argument \code{basis.b} is ignored.
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
# @param CV =TRUE, Cross-validation (CV) is done.
#' @param \dots Further arguments passed to or from other methods.
#' @return Return \code{gam} object plus:
#' \itemize{
#' \item \code{basis.x}: Basis used for \code{fdata} or \code{fd} covariates. 
#' \item \code{basis.b}: Basis used for beta parameter estimation. 
#' \item \code{data}: List containing the variables in the model. 
#' \item \code{formula}: Formula used in the model. 
#' \item \code{y.pred}: Predicted response by cross-validation.
#' }
#' @note If the formula only contains a non functional explanatory variables
#' (multivariate covariates), the function compute a standard \code{\link{glm}}
#' procedure.
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' @seealso See Also as: \code{\link{predict.fregre.gsam}} and
#' \link[mgcv]{summary.gam}.\cr Alternative methods: \code{\link{fregre.glm}}
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
#' x <- tecator$absorp.fdata
#' x.d1 <- fdata.deriv(x)
#' tt <- x[["argvals"]]
#' dataf <- as.data.frame(tecator$y)
#' nbasis.x <- 11
#' nbasis.b <- 5
#' basis1 <- create.bspline.basis(rangeval=range(tt),nbasis=nbasis.x)
#' basis2 <- create.bspline.basis(rangeval=range(tt),nbasis=nbasis.b)
#' f <- Fat ~ s(Protein) + s(x)
#' basis.x <- list("x"=basis1,"x.d1"=basis1)
#' basis.b <- list("x"=basis2,"x.d1"=basis2)
#' ldat <- ldata("df" = dataf, "x" = x , "x.d1" = x.d1)
#' res <- fregre.gsam(Fat ~ Water + s(Protein) + x + s(x.d1), ldat,
#'                    family = gaussian(),  basis.x = basis.x,
#'                    basis.b = basis.b)
#' summary(res)
#' pred <- predict(res,ldat)
#' plot(pred-res$fitted)
#' pred2 <- predict.gam(res,res$XX)
#' plot(pred2-res$fitted)
#' plot(pred2-pred)

#' res2 <- fregre.gsam(Fat ~ te(Protein, k = 3) + x, data =  ldat,
#'                      family=gaussian())
#' summary(res2)
#' 
#' ##  dropind basis pc
#' basis.pc0 <- create.pc.basis(x,c(2,4,7))
#' basis.pc1 <- create.pc.basis(x.d1,c(1:3))
#' basis.x <- list("x"=basis.pc0,"x.d1"=basis.pc1)
#' ldata <- ldata("df"=dataf,"x"=x,"x.d1"=x.d1)  
#' res.pc <- fregre.gsam(f,data=ldata,family=gaussian(),
#'           basis.x=basis.x,basis.b=basis.b)
#' summary(res.pc)
#'  
#' ##  Binomial family
#' ldat$df$FatCat <- factor(ifelse(tecator$y$Fat > 20, 1, 0))
#' res.bin <- fregre.gsam(FatCat ~ Protein + s(x),ldat,family=binomial())
#' summary(res.bin)
#' table(ldat$df$FatCat, ifelse(res.bin$fitted.values < 0.5,0,1))
#' }
#' @export
fregre.gsam <- function (formula
                               , family = gaussian()
                               , data = list(), weights = NULL
                               , basis.x = NULL, basis.b = NULL, ...) 
{
  nam.data <- names(data)
  nam.df <- names(data$df)
  nam.func <- setdiff(nam.data,"df")
  
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
  # print(2)  
  nam.func <- setdiff(nam.data,response)
  
  covar <- gp$fake.names
  if (length(covar)==0) {
    z = gam(formula = gp$pf, data = data$df, family = family)
    z$data <- data
    z$formula.ini <-gp$pf
    z$formula <- gp$pf
    class(z) <- c(class(z), "fregre.gsam")
    return(z)
  }
  
  vnf<-vfunc<-NULL
  for (i in 1:length(covar)){
    if (covar[i] %in% nam.df) vnf<-c(vnf,covar[i])
    if (covar[i] %in% nam.func) vfunc<-c(vfunc,covar[i])
  }
  nnf<-length(vnf)
  nfunc<-length(vfunc)
  bsp1 <-  TRUE
  bspy <-  TRUE
  raw <- FALSE 
  
  # cat("response ");# print(response)
  # cat("scalar ");# print(vnf)
  # cat("funcional ");# print(vfunc)
  
  ######## manejando la respuesta  
  # incluir la respuesta  
  ########################################
  # print("len.smo")
  # print(len.smo)
  # print(specials)
  #cat("names.vfunc");print(names.vfunc)
  if (len.smo == 0) {
    specials <- rep(NA, nterms)
    speci <- rep("0", nterms)
    gp$smooth.spec[1:nterms] <- NULL
    gp$smooth.spec[1:nterms]$term <- NULL
    fnf2<-rep(0,nterms)
    #    vnf <- terms
    fnf1 <- fnf <- nterms
  }  else {
    for (i in 1:nterms) if (!is.na(specials[i])) 
      speci <- c(speci, specials[i])
    for (i in 1:nterms) {
      #  cat("terms[i]")        ;      print(terms[i])        
      
      if (any(terms[i] == nam.df)) {
        #print(1)
        #        vnf <- c(vnf, terms[i])
        
        sm[i] <- nf[i] <- 1
        fnf1 <- c(fnf1, 0)
        bs.dim1 <- c(bs.dim1, 0)
        specials1 <- c(specials1, "0")
      }      else {
        #print(2)
        if (any(terms[i] == nam.func)) {
          # print(3)
          #          vfunc <- c(vfunc, terms[i])
          func[i] <- 1
          fnf2 <- c(fnf2, 0)
          bs.dim2 <- c(bs.dim2, 0)
          specials2 <- c(specials2, "0")
        }
      }
    }
    for (i in 1:len.smo) {
      # cat("i");print(i)
      # print(speci[i])      
      if (speci[i] != "s") {
        # print(4)
        if (any(gp$smooth.spec[[i]]$term == nam.df)) {
          #   print(5)
          #          vnf <- c(vnf, gp$smooth.spec[[i]]$margin[[1]]$term)
          bs.dim1 <- c(bs.dim1, gp$smooth.spec[[i]]$margin[[1]]$bs.dim)
          fnf <- c(fnf, 1)
          fnf1 <- c(fnf1, 1)
          fnf2 <- c(fnf2, 0)
          specials1 <- c(specials1, speci[i])
        }
        else {
          # print(6)
          if (any(gp$smooth.spec[[i]]$term == nam.func)) {
            # print(7)
            # vfunc <- c(vfunc, gp$smooth.spec[[i]]$margin[[1]]$term)
            bs.dim2 <- c(bs.dim2, gp$smooth.spec[[i]]$margin[[1]]$bs.dim)
            fnf <- c(fnf, 2)
            fnf2 <- c(fnf2, 1)
            #fnf2[i] <- 1
            specials2 <- c(specials2, speci[i])
          } else  fnf2 <- c(fnf2, 0)
        }
      }      else {
        # print(8)
        if (any(gp$smooth.spec[[i]]$term == nam.df)) {
          #           print(9)
          #   vnf <- c(vnf, gp$smooth.spec[[i]]$term)
          bs.dim1 <- c(bs.dim1, gp$smooth.spec[[i]]$bs.dim)
          fnf <- c(fnf, 1)
          fnf1 <- c(fnf1, 1)
          fnf2 <- c(fnf2, 1)
          specials1 <- c(specials1, speci[i])
        }        else {
          # print(10)
          if (any(gp$smooth.spec[[i]]$term == nam.func)) {
            #    vfunc <- c(vfunc, gp$smooth.spec[[i]]$term)
            bs.dim2 <- c(bs.dim2, gp$smooth.spec[[i]]$bs.dim)
            fnf <- c(fnf, 2)
            fnf2 <- c(fnf2, 1)
            specials2 <- c(specials2, speci[i])
          } else  fnf2 <- c(fnf2, 0)
        }
      }
    }
  }
  #nfunc <- sum(func)
 #print("aaaaaaaaaaaaaaaaaaa1")  
  name.coef = nam = par.fregre = beta.l = list()
  kterms = 1
  if (nnf > 0) {
    XX = data[["df"]][,c(response,vnf),drop=F]
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
# print("aaaaaaaaaaaaaaaaaaa2")  
  lenfunc <- length(vfunc) 
  ifunc <- lenfunc > 0
  mean.list = basis.list = list()
  if (ifunc) {
    k = 1
    
    for (i in 1:lenfunc ) {
# print("aaaaaaaaaaaaaaaaaaa3")        
      if (is(data[[vfunc[i]]], "fdata")) {
        tt <- data[[vfunc[i]]][["argvals"]]
        rtt <- data[[vfunc[i]]][["rangeval"]]
        fdat <- data[[vfunc[i]]]
        nms <- data[[vfunc[i]]]$names
        #dat <- data[[vfunc[i]]]$data
        if (is.null(basis.x[[vfunc[i]]])) {
          if (is.null(basis.b[[vfunc[i]]])) { basis.b[[vfunc[i]]] <- basis.x[[vfunc[i]]] <- 
            create.fdata.basis(fdat, l = 1:7)
          } else  {basis.x[[vfunc[i]]] <-  basis.b[[vfunc[i]]] }
        } else
            if (basis.x[[vfunc[i]]]$type == "pc" | basis.x[[vfunc[i]]]$type == "pls") 
              bsp1 = FALSE
        
        xaux <- fdata2basis(data[[vfunc[i]]],basis.x[[vfunc[i]]] ,method = c( "inprod"))
        name.coef[[vfunc[i]]] <- colnames(xaux$coefs) <- paste(vfunc[i],".",colnames(xaux$coefs),sep="")
        Z <- xaux$coefs
        ################################### usar basis.b
        if ( bsp1 ){ #bsb1 $ !smo con fnf2==0
          colnames(Z) -> cnames
          J = inprod(basis.x[[vfunc[i]]], basis.b[[vfunc[i]]])
          
          Z = Z %*% J
          colnames(Z)<-cnames[1:NCOL(Z)]
          name.coef[[vfunc[i]]]<-name.coef[[vfunc[i]]][1:NCOL(Z)]
          # print(cnames)
          # print(2222)
          # print(name.coef[[vfunc[i]]])
        }
        # print("aaaaaaaaaaaaaaaaaaa444")         
        lencoef <- length(colnames(Z))
        
        ################################### usar basis.b
        XX = cbind(XX, Z)
        # print("aaaaaaaaaaaaaaaaaaa5")         
        if (fnf2[i] == 1)    
          sm2 <- TRUE    else sm2 <- FALSE
        for (j in 1:lencoef) {
          if (sm2) {
            pf <- paste(pf, "+", specials2[i], "(", 
                        name.coef[[vfunc[i]]][j], ",k=", bs.dim2[i],")", sep = "")
          }     else pf <- paste(pf, "+", name.coef[[vfunc[i]]][j], sep = "")
          kterms <- kterms + 1
          # print(pf)
        }       
        basis.list[[vfunc[i]]] <- xaux$basis
        # J=inprod(basis.x[[vfunc[i]]],basis.b[[vfunc[i]]])
        #   vs.list[[vfunc[i]]] = basis.x[[vfunc[i]]]$basis
        if (bspy!=0) 
          mean.list[[vfunc[i]]] = basis.x[[vfunc[i]]]$mean
        else {
          xcc <- fdata.cen(data[[vfunc[i]]])
          mean.list[[vfunc[i]]] = xcc[[2]]
        }
      }      else {
        stop("Please, enter functional covariate of fdata class object")
      }
    }
  }
  #par.fregre$formula=as.formula(pf)
  #par.fregre$data=XX
  #formula=as.formula(pf)
  names(par.fregre)
  nx <- ncol(XX)
  #Ymat<-y
  #Xmat<-data.frame(y)
  # if (length(vnf)>0) {
  #   Xmat<-data.frame(Xmat,data$df[,vnf])
  #   names(Xmat)<-c(response,vnf)         
  # }  else       names(Xmat) <- response                   
  # # for (i in 2:length(XX)){
  #    XX1<-data.frame(XX1,XX[[i]])
  #  }
  #Xmat <- (cbind(Xmat,XX))
  # print(colnames(XX))
  #for (i in 1:npy){
  #Xmat[,response] <- Ymat[,i]
  #        z=gam(formula=as.formula(par.fregre$formula),data=XX1,family=family,offset=rep(1,len=nrow(XX[[1]])))
  #if (missing(offset))   
# print(4)  
# print(pf)
# print(head(XX))
  z=gam(formula=as.formula(pf),data=XX,family=family,weights=weights)  
  # else   { descomentar ***** y missing(offset)
  #   off<-offset
  #   z=gam(formula=formula,data=Xmat,family=family,offset=off)
  # }
  # }
  #  print(bsp1)
  z$mean <- mean.list
  z$formula <- as.formula(pf)
  z$formula.ini <- formula
  z$basis.x=basis.x
  z$basis.b=basis.b
  z$basis.list<-basis.list
  z$data=data
  z$XX <- XX
  #z$basis <- aux$basis
  z$bsp <- bsp1
  z$vfunc <- vfunc;   z$nnf <- nnf;   z$vnf <- vnf
  z$fnf2<-fnf2
  class(z) <- c("fregre.gsam",class(z))
  z
}

# si tienes s(x.d1) no usar el basis.b y guardar la sublista utilizada
# res=fregre.gsam(Fat~Water+s(Protein)+x+s(x.d1),ldata,family=gaussian(),
#                 basis.x=basis.x,basis.b=basis.b)
# summary(res)
# pred <-predict(res,ldata)
# 
# res=fregre.gsam(Fat~Water+s(Protein)+x+s(x.d1),ldata,family=gaussian())
# res$basis.x$x$nbasis
# 
# # buscar fda.usc.devel::: en fda.usc!!