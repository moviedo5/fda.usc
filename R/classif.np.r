#' @title Kernel Classifier from Functional Data
#' 
#' @description Fits Nonparametric Supervised Classification for Functional Data.
#' 
#' @details Make the group classification of a training dataset using kernel or KNN
#' estimation: \code{\link{Kernel}}.\cr Different types of metric funtions can
#' be used.
#' 
#' @aliases classif.np classif.kernel classif.knn
#' @param group Factor of length \emph{n}
#' @param fdataobj \code{\link{fdata}} class object.
#' @param h Vector of smoothing parameter or bandwidth.
#' @param knn Vector of number of nearest neighbors considered.
#' @param Ker Type of kernel used.
#' @param metric Metric function, by default \code{\link{metric.lp}}.
#' @param weights weights.
#' @param type.S Type of smothing matrix \code{S}. By default \code{S} is
#' calculated by Nadaraya-Watson kernel estimator (\code{S.NW}).
#' @param par.S List of parameters for \code{type.S}: \code{w}, the weights.
#' @param \dots Arguments to be passed for \code{\link{metric.lp}} o other
#' metric function and \code{\link{Kernel}} function.
#' @return \itemize{
#' \item {fdataobj}{ \code{\link{fdata}} class object.} 
#' \item {group}{ Factor of length \code{n}.} 
#' \item {group.est}{ Estimated vector groups}
#' \item {prob.group}{ Matrix of predicted class probabilities. For each
#' functional point shows the probability of each possible group membership.}
#' \item {max.prob}{ Highest probability of correct classification.}
#' \item {h.opt}{ Optimal smoothing parameter or bandwidht estimated.} 
#' \item {D}{ Matrix of distances of the optimal quantile distance \code{hh.opt}.}
#' \item {prob.classification}{ Probability of correct classification by group.}
#' \item {misclassification}{ Vector of probability of misclassification by
#' number of neighbors \code{knn}.} 
#' \item {h}{ Vector of smoothing parameter or bandwidht.} 
#' \item {C}{ A call of function \code{classif.kernel}.}
#' }
#' @note If \code{fdataobj} is a data.frame the function considers the case of
#' multivariate covariates. \cr \code{\link{metric.dist}} function is used to
#' compute the distances between the rows of a data matrix (as
#' \code{\link{dist}} function.
#' @author Manuel Febrero-Bande, Manuel Oviedo de la Fuente
#' \email{manuel.oviedo@@udc.es}
#' @seealso See Also as \code{\link{predict.classif}}
#' @references Ferraty, F. and Vieu, P. (2006). \emph{Nonparametric functional
#' data analysis.} Springer Series in Statistics, New York.
#' 
#' Ferraty, F. and Vieu, P. (2006). \emph{NPFDA in practice}. Free access on
#' line at \url{http://www.lsp.ups-tlse.fr/staph/npfda/}
#' @keywords classif
#' @examples
#' \dontrun{
#' data(phoneme)
#' mlearn<-phoneme[["learn"]]
#' glearn<-phoneme[["classlearn"]]
#' 
#' h=9:19
#' out=classif.np(glearn,mlearn,h=h)
#' summary(out)
#' # round(out$prob.group,4)
#' }
#' 
#' @rdname classif.np
#' @export 
classif.np <- function  (group, fdataobj, h = NULL, Ker = AKer.norm, metric, 
                         weights = "equal", type.S = S.NW,
                         par.S = list()
                         #, measure = "accuracy"
                         , ...) 
{
#  print("entra np2")
  y <- group
  n <- length(y)
  if (is.character(weights)) {
    weights <- weights4class(y, type = weights)
  }
  else {
    if (length(weights) != n) 
      stop("length weights != length response")
  }
  if (missing(metric)) {
    if (is.fdata(fdataobj)) 
      metric = metric.lp
    else metric = metric.dist
  }
  if (is.fdata(fdataobj)) {
    nas <- is.na.fdata(fdataobj)
    nas.g <- is.na(y)
    C <- match.call()
    if (is.null(names(y))) 
      names(y) <- 1:length(y)
    if (any(nas) & !any(nas.g)) {
      bb <- !nas
      cat("Warning: ", sum(nas), " curves with NA are omited\n")
      fdataobj$data <- fdataobj$data[bb, ]
      y <- y[bb]
    }
    else {
      if (!any(nas) & any(nas.g)) {
        cat("Warning: ", sum(nas.g), " values of group with NA are omited \n")
        bb <- !nas.g
        fdataobj$data <- fdataobj$data[bb, ]
        y <- y[bb]
      }
      else {
        if (any(nas) & any(nas.g)) {
          bb <- !nas & !nas.g
          cat("Warning: ", sum(!bb), " curves  and values of group with NA are omited \n")
          fdataobj$data <- fdataobj$data[bb, ]
          y <- y[bb]
        }
      }
    }
    x <- fdataobj[["data"]]
    tt <- fdataobj[["argvals"]]
    rtt <- fdataobj[["rangeval"]]
  }
  else {
    x <- fdataobj
    mdist = metric.dist(x, ...)
  }
  C <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("y", "fdataobj", "weights", 
               "h", "Ker", "metric", "type.S", 
               "par.S"), names(mf), 0L)
  np <- ncol(x)
  if (n != (length(y))) 
    stop("ERROR IN THE DATA DIMENSIONS")
  if (is.null(rownames(x))) 
    rownames(x) <- 1:n
  if (is.null(colnames(x))) 
    colnames(x) <- 1:np
  types = FALSE
  if (is.matrix(metric)) {
    mdist <- metric
    metric <- attributes(mdist)
  }
  else mdist = metric(fdataobj, fdataobj, ...)
  ty <- deparse(substitute(type.S))
  if (is.null(h)) 
    h = h.default(fdataobj, metric = mdist, type.S = ty, 
                  ...)
  else {
    if (any(h <= 0)) 
      stop("Error: Invalid range for h")
  }
  lenh <- length(h)
  gcv <- cv.error <- array(NA, dim = c(lenh))
  par.S2 <- par.S
  lenh = length(h)
  if (!is.factor(group)) 
    group <- as.factor(group)
  group <- y <- factor(group, levels = levels(group)[which(table(group) > 
                                                             0)])
  ny <- levels(y)
  numg = nlevels(y)
  Y = array(0, dim = c(numg, n))
  group.est2 = group.est = array(0, dim = c(lenh, n))
  pgrup = array(0, dim = c(numg, n, lenh))
  misclassification = array(1, dim = c(1, lenh))
  pr <- 1
  if (is.null(par.S2$h)) 
    par.S$h <- h
  if (is.null(par.S$Ker) & ty != "S.KNN") 
    par.S$Ker <- Ker
  if (is.null(par.S$w)) 
    par.S$w <- weights
  par.fda.usc <- eval(parse(text = "fda.usc.devel:::par.fda.usc"), 
                      envir = .GlobalEnv)
  warn <- par.fda.usc$warning
  for (i in 1:lenh) {
    par.S$tt <- mdist
    par.S$h <- h[i]
    if (is.null(par.S$cv)) 
      par.S$cv = TRUE
    H = do.call(ty, par.S)
    for (j in 1:numg) {
      Y[j, ] = as.integer(y == ny[j])
      pgrup[j, , i] <- (H %*% matrix(Y[j, ], ncol = 1))
    }
    if (ty == "S.KNN") {
      for (ii in 1:n) {
        l = seq_along(pgrup[, ii, i])[pgrup[, ii, i] == 
                                        max(pgrup[, ii, i], na.rm = T)]
        if (length(l) > 1) {
          ll <- which(levels(y[1:n]) %in% l)
          abc <- which(mdist[ii, ll] == min(mdist[ii, 
                                                  ll], na.rm = T))
          group.est[i, ii] = levels(y)[y[ll[abc[1]]]]
        }
        else group.est[i, ii] = ny[l[1]]
      }
    }
    else {
      group.est[i, ] <- ny[as.vector(apply(pgrup[, , i], 
                                           2, which.max))]
    }
    ###################### 
    lo <- y != group.est[i, ]
    #ypred <- factor(group.est[i,],levels=ny)
    #lo <- cat2meas(y,ypred,measure=measure)
    gcv[i] = weighted.mean(lo, par.S$w, na.rm = TRUE)
    #gcv[i] = 1-cat2meas(y,ypred,measure=measure)
    if (pr > gcv[i]) {
      pr = gcv[i]
      iknn = i
      prob = 1 - pr
      prob.group2 = t(pgrup[, , i])
      group.pred = group.est[i, ]
    }
  }
  colnames(prob.group2) <- ny
  rownames(prob.group2) <- rownames(x)
  l = which.min(gcv)
  h.opt <- h[l]
  par.S$h <- h.opt
  if (is.null(par.S$cv)) 
    par.S$cv = FALSE
  if (h.opt == min(h) & warn) 
    cat(" Warning: h.opt is the minimum value of bandwidths\n   provided, range(h)=", 
        range(h), "\n")
  else if (h.opt == max(h) & warn) 
    cat(" Warning: h.opt is the maximum value of bandwidths\n   provided, range(h)=", 
        range(h), "\n")
  df = trace.matrix(H)
  names(gcv) <- h
  group.pred <- factor(group.pred, levels = ny)
  misclass = weighted.mean(group.pred != y, par.S$w)
  prob.classification <- diag(table(y, group.pred))/table(y)
  out <- list(C = C, group.est = group.pred, group = y, H = H, 
              df = df, y = y, fdataobj = fdataobj, mdist = mdist, Ker = Ker, 
              metric = metric, type.S = type.S, par.S = par.S, gcv = gcv, 
              h.opt = h.opt, h = h, prob.group = prob.group2, m = m, 
              pgrup = pgrup, ty = ty, prob.classification = prob.classification, 
              max.prob = 1 - misclass)
  class(out) = "classif"
  return(out)
}

#' @rdname classif.np
#' @export classif.knn
classif.knn=function (group, fdataobj, knn = NULL, metric, weights = "equal", 
                      par.S = list(), ...) 
{
  if (missing(metric)) {
    if (is.fdata(fdataobj)) 
      metric = metric.lp
    else metric = metric.dist
  }
  if (is.null(knn)) 
    knn <- seq(3, min(15, length(group)), by = 2)
  classif.np(group, fdataobj, h = knn, Ker = Ker.unif, metric = metric, 
             weights = weights, type.S = S.KNN, par.S = par.S, ...)
}

#' @rdname classif.np
#' @export classif.kernel
classif.kernel=function(group, fdataobj, h = NULL, Ker = AKer.norm, metric, 
                        weights = "equal", par.S = list(), ...) 
{
  if (missing(metric)) {
    if (is.fdata(fdataobj)) 
      metric = metric.lp
    else metric = metric.dist
  }
  classif.np(group, fdataobj, h, Ker, metric = metric, weights, 
             type.S = S.NW, par.S = par.S, ...)
}



# @param type.CV Type of cross-validation. By default generalized
# cross-validation \code{\link{GCV.S}} method.
# @param par.CV List of parameters for \code{type.CV}: \code{trim}, the alpha
# of the trimming and \code{draw=TRUE}.