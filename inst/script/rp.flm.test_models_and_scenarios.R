

#' @title Utilities for the simulation of functional linear models
#' 
#' @description The functions \code{r.cfs.2003}, \code{r.hh.2006} and \code{r.bridge} sample the functional covariate and construct functional coefficients for their use in functional linear models:
#' \itemize{
#'   \item{\code{r.cfs.2003} implements examples (a) and (b) in Section 5 of Cardot et al. (2003).}
#'   \item{\code{r.hh.2006} gives models (i), (ii) and (iii) in Section 5 of Hall and Housseini-Nasab (2006).}
#'   \item{\code{r.bridge} samples a Brownian motion and creates a functional coefficient constructed from the eigenfunctions \eqn{\sqrt(2) * \sin(k t \pi)}{\sqrt2 * sin(k*t*\pi)}.}
#' }
#' @param n number of functions to sample.
#' @inheritParams r.mod
#' @param b coefficients of the functional coefficient in the theoretical basis of principal components of the Brownian motion.
#' @param type either example \code{"a"} or \code{"b"} from Cardot et al. (2003).
#' @param imod either \code{1}, \code{2} or \code{3} for denoting models (i), (ii) and (iii) in Hall and Hosseini-Nasab (2006).
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{X.fdata}}{the sample of functional data, an \code{\link[fda.usc]{fdata}} object of length \code{n}.}
#'   \item{\code{beta.fdata}}{the functional coefficient, an \code{\link[fda.usc]{fdata}} object.}
#' }
#' @examples 
#' # Cardot et al. (2003)
#' plot(r.cfs.2003(n = 100, type = "a")$X.fdata)
#' plot(r.cfs.2003(n = 100, type = "b")$X.fdata)
#' 
#' # Hall and Hosseini-Nasab (2006)
#' plot(r.hh.2006(n = 100, imod = 1)$X.fdata)
#' plot(r.hh.2006(n = 100, imod = 2)$X.fdata)
#' plot(r.hh.2006(n = 100, imod = 3)$X.fdata)
#' 
#' # Sample bridge
#' plot(r.bridge(n = 100)$X.fdata)
#' @author Manuel Febrero-Bande (\email{manuel.febrero@@usc.es}) and Eduardo Garcia-Portugues (\email{edgarcia@@est-econ.uc3m.es}).
#' @references 
#' Cardot, H., Ferraty, F., Sarda, P. (2003) Spline estimators for the functional linear model. Statistica Sinica, 13(3), 571--592. \url{https://www3.stat.sinica.edu.tw/statistica/oldpdf/a13n31.pdf}
#' 
#' Hall, P. and Hosseini-Nasab, M. (2006) On properties of functional principal components analysis. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 68(1), 109--126. \url{https://dx.doi.org/10.1111/j.1467-9868.2005.00535.x}
#' @export 
r.cfs.2003 <- function(n, t = seq(0, 1, len = 201), b = c(2, 4, 5) / sqrt(2), 
                       type = c("a", "b")[1]) {
  
  # X.fdata
  X.fdata <- fda.usc::rproc2fdata(n = n, t = t, sigma = "brownian")
  
  if (type == "a") {
    
    # beta = b1 * v1 + b2 * v2 + b3 * v3
    v <- fda.usc::fdata(mdata = sqrt(2) * sin(t * pi * 0.5), argvals = t)
    for (k in 2:length(b)) {
      
      v <- c(v, fda.usc::fdata(mdata = sqrt(2) * sin(t * pi * (k - 0.5)), 
                               argvals = t))
      
    }
    beta.fdata <- fda.usc::fdata(mdata = matrix(b, nrow = 1) %*% v$data, argvals = t)
    
  } else if (type == "b") {
    
    beta.fdata <- fda.usc::fdata(mdata = log(15 * t^2 + 10) + cos(4 * pi * t), 
                                 argvals = t)
    
  } else {
    
    stop("Wrong type")
    
  }
  
  return(list("X.fdata" = X.fdata, "beta.fdata" = beta.fdata))
  
}


#' @rdname r.cfs.2003
#' @export 
r.hh.2006 <- function(n, t = seq(0, 1, len = 201), imod = 1) {
  
  # Basis
  v <- fda.usc::fdata(mdata = sqrt(2) * cos(pi * t), argvals = t)
  for (k in 2:20) {
    
    v <- c(v, fda.usc::fdata(mdata = sqrt(2) * cos(k * pi * t), argvals = t))
    
  }
  
  # beta
  j <- 1:20
  beta.fdata <- fda.usc::fdata(mdata = (2^(3/2) * (-1)^j * j^(-2)) %*% v$data, 
                               argvals = t)
  
  # X.fdata
  sdcoef <- j^(-imod)
  coefsim <- matrix(NA, ncol = 20, nrow = n)
  for (i in j) {
    
    coefsim[, i] <- rnorm(n, mean = 0, sd = sdcoef[i])
    
  }
  X.fdata <- fda.usc::fdata(mdata = coefsim %*% v$data, argvals = t)
  
  return(list("X.fdata" = X.fdata, "beta.fdata" = beta.fdata))
  
}


#' @rdname r.cfs.2003
#' @export 
r.bridge <- function(n, t = seq(0, 1, len = 201), b = c(2, 4, 5) / sqrt(2)) {
  
  # X.fdata
  X.fdata <- fda.usc::rproc2fdata(n = n, t = t, sigma = "brownian")
  lt <- length(t)
  for (i in 1:n) {
    
    X.fdata$data[i, ] <- X.fdata$data[i, ] - t * X.fdata$data[i, lt]
    
  }
  
  # beta = b1 * v1 + b2 * v2 + b3 * v3
  v <- fda.usc::fdata(mdata = sqrt(2) * sin(t * pi), argvals = t)
  for (k in 2:length(b)) {
    
    v <- c(v, fda.usc::fdata(mdata = sqrt(2) * sin(t * pi * k), argvals = t))
    
  }
  beta.fdata <- fda.usc::fdata(mdata = matrix(b, nrow = 1) %*% v$data, argvals = t)
  
  return(list("X.fdata" = X.fdata, "beta.fdata" = beta.fdata))
  
}


#' @title Geometric Brownian motion
#' 
#' @description Sampling of paths of the geometric Brownian motion.
#' 
#' @inheritParams r.cfs.2003
#' @param mu,sigma mean and diffusion of the underlying Brownian motion.
#' @param s0 a number or a vector of length \code{n} giving the initial value(s) of the geometric Brownian motion.
#' @return Functional sample, an \code{\link[fda.usc]{fdata}} object of length \code{n}.
#' @examples 
#' plot(r.gbm(n = 10, s0 = 1))
#' plot(r.gbm(n = 10, mu = -1, s0 = rnorm(10, mean = 10)))
#' @author Eduardo Garcia-Portugues (\email{edgarcia@@est-econ.uc3m.es}).
#' @export
r.gbm <- function(n, t = seq(0, 1, len = 201), mu = 0, sigma = 1, s0 = 1) {
  
  # Time-varying covariances
  St <- sigma^2 * outer(t, t, function(s, t) pmin(s, t))
  
  # Sample N((mu - sigma^2 / 2) * t, St)
  #mdata <- mvtnorm::rmvnorm(n = n, mean = (mu - sigma^2 / 2) * t, sigma = St)
  mdata <- mvrnorm(n = n, mu = (mu - sigma^2 / 2) * t, Sigma = St)
  mdata <- s0 * exp(mdata)
  
  # As fdata object
  return(fda.usc::fdata(mdata = mdata, argvals = t))
  
}



#' @title Functional regression models for simple and composite hypothesis
#'
#' @description Sampling from the functional linear models considered in the simulation study of Cuesta-Albertos et al. (2017).
#'
#' @param n the sample size.
#' @param t time locations for the functional data.
#' @param scenario an index from \code{1} to \code{12} denoting the simulation scenario.
#' @param delta an index from \code{0} to \code{3} denoting the degree of departure of the data from the null hypothesis of functional linearity, encoded with \code{0}.
#' @param R2 proportion of variance of the the response \eqn{Y}{Y} explained by the linear model when \code{delta = 0}. This is used to compute the variance of the error \eqn{\varepsilon}{\epsilon} of the regression model.
#' @param composite flag to indicate the generation of data according to a functional linear model with non-null coefficient (\code{TRUE}) or with a null coefficient (\code{FALSE}).
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{X.fdata}}{the sample of functional data, an \code{\link[fda.usc]{fdata}} object of length \code{n}.}
#'   \item{\code{Y}}{the scalar responses, a vector of length \code{n}.}
#'   \item{\code{beta.fdata}}{the functional coefficient, an \code{\link[fda.usc]{fdata}} object.}
#' }
#' @details The samples are generated from the regression model
#' \deqn{Y = \langle \mathcal{X}, \beta\rangle + \delta m(\mathcal{X})+\varepsilon,}{Y = <X, \beta> + \delta m(X)+\varepsilon,}
#' where \eqn{\delta m(\mathcal{X})}{\delta m(X)} is computed by \code{\link{m.dev}}. The description of the scenarios is detailed in the supplementary material of Cuesta-Albertos et al. (2017).
#' @examples
#' # Generate samples for all scenarios
#' samp <- list()
#' k <- 1
#' for (i in 1:12) {
#'   for (delta in 0:3) {
#'     samp[[k]] <- r.mod(n = 10, scenario = i, delta = delta, R2 = 0.95, composite = TRUE)
#'     k <- k + 1
#'   }
#' }
#' @author Eduardo Garcia-Portugues (\email{edgarcia@@est-econ.uc3m.es}).
#' @references
#' Cuesta-Albertos, J.A., Garcia-Portugues, E., Febrero-Bande, M. and Gonzalez-Manteiga, W. (2017). Goodness-of-fit tests for the functional linear model based on randomly projected empirical processes. arXiv:1701.08363. \url{https://arxiv.org/abs/1701.08363}
#' @export
r.mod <- function(n, scenario, delta = 0, t = seq(0, 1, l = 201), R2 = 0.95,
                  composite = TRUE) {

  # Check for delta
  if (!(delta %in% 0:3)) {

    stop("delta must be 0, 1, 2 or 3")

  }

  # Different models with deviations
  if (scenario == 1) {

    # S1: example (a) of CFS 2003, R2 = 0.95, three PC (1, 2, 3)
    b <- c(2, 4, 5) / sqrt(2)
    cfs <- r.cfs.2003(n = n, t = t, b = b, type = "a")
    X.fdata <- cfs$X.fdata
    beta0 <- cfs$beta.fdata

    # Deviations
    mdev <- m.dev(X.fdata, type = 1, delta = c(0, 0.25, 0.75, 1.25)[delta + 1])

  } else if (scenario == 2) {

    # S2: example (a) of CFS 2003, R2 = 0.95, three PC (3, 5, 7)
    b <- c(0, 0, 2, 0, 4, 0, 5) / sqrt(2)
    cfs <- r.cfs.2003(n = n, t = t, b = b, type = "a")
    X.fdata <- cfs$X.fdata
    beta0 <- cfs$beta.fdata

    # Deviations
    mdev <- m.dev(X.fdata, type = 1, delta = c(0, 0.05, 0.2, 0.5)[delta + 1])

  } else if (scenario == 3) {

    # S3: example (a) of CFS 2003, R2 = 0.95, three PC (2, 3, 7)
    b <- c(0, 2, 4, 0, 0, 0, 5) / sqrt(2)
    cfs <- r.cfs.2003(n = n, t = t, b = b, type = "a")
    X.fdata <- cfs$X.fdata
    beta0 <- cfs$beta.fdata

    # Deviations
    mdev <- m.dev(X.fdata, type = 1, delta = -c(0, 0.2, 0.5, 1)[delta + 1])

  } else if (scenario == 4) {

    # S4: example b of CFS 2003, beta = log(15 * t^2 + 10) + cos(4 * pi * t)
    cfs <- r.cfs.2003(n = n, t = t, type = "b")
    X.fdata <- cfs$X.fdata
    beta0 <- cfs$beta.fdata

    # Deviations
    mdev <- m.dev(X.fdata, type = 1, delta = c(0, 0.2, 1, 2)[delta + 1])

  } else if (scenario == 5) {

    # S5: example Hall and Hoseini, imod = 1
    hh <- r.hh.2006(n = n, t = t, imod = 1)
    X.fdata <- hh$X.fdata
    beta0 <- hh$beta.fdata

    # Deviations
    mdev <- m.dev(X.fdata, type = 2, delta = -c(0, 1, 3, 7)[delta + 1])

  } else if (scenario == 6) {

    # S6: example Hall & Hoseini, imod = 2
    hh <- r.hh.2006(n = n, t = t, imod = 2)
    X.fdata <- hh$X.fdata
    beta0 <- hh$beta.fdata

    # Deviations
    mdev <- m.dev(X.fdata, type = 2, delta = -c(0, 1, 3, 7)[delta + 1])

  } else if (scenario == 7) {

    # S7: brownian bridge with exact representation in the three first PC, R2 = 0.95
    bridge <- r.bridge(n = n, t = t, b = c(2, 4, 5) / sqrt(2))
    X.fdata <- bridge$X.fdata
    beta0 <- bridge$beta.fdata

    # Deviations
    mdev <- m.dev(X.fdata, type = 2, delta = -c(0, 2, 7.5, 15)[delta + 1])

  } else if (scenario == 8) {

    # S8: Ornstein-Uhlenbeck and quadratic deviation I as in G-P, F-B and G-M (2014)
    X.fdata <- r.ou(n = n, t = t, alpha = 1/3, sigma = 1)
    beta0 <- fda.usc::fdata(mdata = sin(2 * pi * t) - cos(2 * pi * t),
                            argvals = t)

    # Deviations
    mdev <- m.dev(X.fdata, type = 2, delta = -c(0, 0.25, 1, 2)[delta + 1])

  } else if (scenario == 9) {

    # S9: Ornstein-Uhlenbeck and quadratic deviation II as in G-P, F-B and G-M (2014)
    X.fdata <- r.ou(n = n, t = t, alpha = 1/3, sigma = 1)
    beta0 <- fda.usc::fdata(mdata = t - (t - 0.75)^2, argvals = t)

    # Deviations
    mdev <- m.dev(X.fdata, type = 3, delta = -c(0, 0.01, 0.1, 0.3)[delta + 1])

  } else if (scenario == 10) {

    # S10: Ornstein-Uhlenbeck and quadratic deviation III as in G-P, F-B and G-M (2014)
    X.fdata <- r.ou(n = n, t = t, alpha = 1/3, sigma = 1)
    beta0 <- fda.usc::fdata(mdata = t + cos(2 * pi * t), argvals = t)

    # Deviations
    mdev <- m.dev(X.fdata, type = 3, delta = c(0, 0.01, 0.05, 0.25)[delta + 1])

  } else if (scenario == 11) {

    # S11: geometrical brownian motion I
    X.fdata <- r.gbm(n = n, t = t, s0 = 1)
    beta0 <- fda.usc::fdata(log(15 * t^2 + 10) + cos(4 * pi * t), argvals = t)

    # Deviations
    mdev <- m.dev(X.fdata, type = 3, delta = c(0, 0.5, 2, 5)[delta + 1])

  } else if (scenario == 12) {

    # S12: geometrical brownian motion II
    X.fdata <- r.gbm(n = n, t = t, s0 = 2)
    beta0 <- fda.usc::fdata(pi^2 * (t^2 - 1/3), argvals = t)

    # Deviations
    mdev <- m.dev(X.fdata, type = 3, delta = c(0, 0.5, 2.5, 7)[delta + 1])

  } else if (scenario == 13) {

    # Toy example (Eduardo)
    X.fdata <- r.ou(n = n, t = t)
    for (i in 1:n) {

      X.fdata$data[i, ] <- rnorm(1, mean = 0, sd = 2) * exp(-t) +
        rnorm(1, mean = 1, sd = 1) * exp(t) +
        rnorm(1, mean = -1, sd = 0.1) * exp(2 * t)

    }
    beta0 <- fda.usc::fdata(mdata = 1 - 10 * exp(-t) - 5 * exp(t) + 3 * exp(2 * t),
                            argvals = t)

    # Deviations
    mdev <- m.dev(X.fdata, type = 2, delta = -c(0, 0.25, 1, 2)[delta + 1])

  } else if (scenario == 14) {

    # Toy example (Manuel)
    X.fdata <- fda.usc::rproc2fdata(n = n, t = t, sigma = "vexponential")
    b5 <- fda::create.bspline.basis(rangeval = range(t), nbasis = 5)
    b5 <- fda.usc::fdata(t(fda::eval.basis(b5, t)), t)
    beta0 <- fda.usc::fdata(apply(sweep(b5$data, 1, c(1, 2, 3, 2, 1), "*"), 2, sum), t)

    # Deviations
    mdev <- m.dev(X.fdata, type = 2, delta = -c(0, 0.25, 1, 2)[delta + 1])

  } else {

    stop("model must be from 1 to 14")

  }

  # Variance explained by the linear regression
  s2y <- switch(scenario,
                1.37885641, 0.10200598, 0.25247980,
                2.56708081, 8.09590848, 7.98459833,
                0.54664892, 0.04646039, 0.18827926,
                0.33637279, 3.43004737, 5.62054178,
                76.674612, 1.385866)

  # Response
  Y <- drop(fda.usc::inprod.fdata(X.fdata, beta0))
  noise <- rnorm(n, mean = 0, sd = sqrt(s2y * (1/ R2 - 1)))
  Y <- switch(composite + 1, 0, Y) + mdev + noise

  return(list(X.fdata = X.fdata, Y = Y, beta.fdata = beta0))

}


#' @title Deviations from functional linearity
#'
#' @description Deviations from functional linearity considered in the simulation study of Cuesta-Albertos et al. (2017).
#'
#' @inheritParams rp.flm.test
#' @param type kind of deviation, an index from \code{1} to \code{5}.
#' @inheritParams r.mod
#' @param eta functional parameter employed when \code{type = 4}.
#' @return A vector of length \code{length(X.fdata)} containing \eqn{\delta m(\mathcal{X})}{\delta m(X)}.
#' @details The description of the deviations is detailed in the supplementary material of Cuesta-Albertos et al. (2017).
#' @examples
#' dev <- list()
#' k <- 1
#' mod <- r.cfs.2003(n = 100)
#' for (i in 1:5) {
#'   for (delta in 0:3) {
#'     dev[[k]] <- m.dev(X.fdata = mod$X.fdata, type = i, eta = mod$beta.fdata,
#'                       delta = delta, composite = TRUE)
#'     k <- k + 1
#'   }
#' }
#' @author Eduardo Garcia-Portugues (\email{edgarcia@@est-econ.uc3m.es}).
#' @references
#' Cuesta-Albertos, J.A., Garcia-Portugues, E., Febrero-Bande, M. and Gonzalez-Manteiga, W. (2017). Goodness-of-fit tests for the functional linear model based on randomly projected empirical processes. arXiv:1701.08363. \url{https://arxiv.org/abs/1701.08363}
#' @export
m.dev <- function(X.fdata, type, delta, eta, composite = TRUE) {

  if (delta == 0) {

    return(rep(0, length(X.fdata)))

  } else {

    if (type == 1) {

      # Norm
      Y <- as.numeric(fda.usc::norm.fdata(X.fdata))

    } else if (type == 2) {

      # Quadratic model
      B <- 25 * outer(fda.usc::argvals(X.fdata), fda.usc::argvals(X.fdata),
                 function(s, t) sin(2 * pi * t * s) * s * (1 - s) * t * (1 - t))
      Y <- diag(X.fdata$data %*% B %*% t(X.fdata$data)) * diff(X.fdata$argvals[1:2])^2

    } else if (type == 3) {

      # Interior product of X.fdata non-linear transformations
      Y <- sapply(1:dim(X.fdata$data)[1], function(i) {
        fda.usc::inprod.fdata(exp((-1) * X.fdata[i]), X.fdata[i]^2)
        })

    } else if (type == 4) {

      # Non-linear transformation of the interior product
      xx <- drop(fda.usc::inprod.fdata(X.fdata, eta))
      Y <- xx * log(xx^2)

    } else if (type == 5) {

      # Difference between ranges in absolute value
      Y <- apply(apply(abs(X.fdata$data), 1, range), 2, diff)

    } else {

      stop("type must be 1, 2, 3, 4 or 5")

    }

    return(delta * Y)

  }

}


#' @title Check density of the responses and functional processes
#'
#' @description \code{check.scenarios} produces plots displaying th densities of the response for different deviations. \code{check.betas} displays a sample of the functional process, the functional coefficient associated and its estimation based on 100 observations.
#'
#' @inheritParams r.mod
#' @param scenarios a vector giving the simulation scenarios to be checked.
#' @param times flag to indicate whether to show the time employed computing each model.
#' @param M number of samples employed to estimate the densities of the response.
#' @param est.beta flag to indicate whether to plot the functional coefficient estimate.
#' @param main whether to show a caption indicating the scenario or not.
#' @return Nothing. The functions are called for producing diagnostic plots.
#' @examples
#' \dontrun{
#' # Check scenarios and deviations
#' set.seed(3257641)
#' check.scenarios(scenarios = c(1, 7, 3, 5, 6, 4, 8, 9, 12), composite = TRUE)
#'
#' # Check betas
#' set.seed(3257641)
#' check.betas(scenarios = c(1, 7, 3, 5, 6, 4, 8, 9, 12), n = 100)
#' }
#' @author Eduardo Garcia-Portugues (\email{edgarcia@@est-econ.uc3m.es}).
#' @export
check.scenarios <- function(scenarios = 1:12, composite = TRUE, times = TRUE,
                            R2 = 0.95, M = 1e3, main = FALSE) {

  par(mfrow = c(3, ceiling(length(scenarios) / 3L)), mar = c(4, 4, 1.5, 1) + 0.1)
  for (k in scenarios) {

    # Response densities
    t <- proc.time()[3]
    set.seed(12456789)
    d0 <- density(r.mod(n = M, scenario = k, delta = 0, composite = composite,
                        R2 = R2)$Y, bw = "SJ")
    set.seed(12456789)
    d1 <- density(r.mod(n = M, scenario = k, delta = 1, composite = composite,
                        R2 = R2)$Y, bw = "SJ")
    set.seed(12456789)
    d2 <- density(r.mod(n = M, scenario = k, delta = 2, composite = composite,
                        R2 = R2)$Y, bw = "SJ")
    set.seed(12456789)
    if (times) {

      cat("Scenario", k, "finished in", proc.time()[3] - t, "secs\n")

    }

    # Plot
    plot(d0, type = "l", lwd = 2, main = ifelse(main, paste("Scenario", k), ""),
         xlab = "", ylab = "", ylim = c(0, max(d0$y, d1$y, d2$y) * 1.2))
    title(xlab = expression(Y == paste(symbol("\xe1"), list(X, rho),
                                       symbol("\xf1")) + delta[k] * m(X) + epsilon),
          ylab = "Density", cex.lab = 1.25, cex.main = 1.25)
    lines(d1, col = 2)
    lines(d2, col = 4)
    legend("topright", legend = expression(delta[0], delta[1], delta[2]),
           lwd = 2, col = c(1:2, 4), cex = 1)

  }

}


#' @rdname check.scenarios
#' @export
check.betas <- function(scenarios = 1:12, n = 100, R2 = 0.95, times = TRUE,
                        est.beta = TRUE, main = FALSE) {

  mar <- c(2, 2, 1, 2) + 0.1
  par(mfrow = c(3, ceiling(length(scenarios) / 3L)), mar = c(4, 4, 1.5, 1) + 0.1)
  for (k in scenarios) {

    # Sample
    t <- proc.time()[3]
    samp <- r.mod(n = n, scenario = k, delta = 0, composite = TRUE, R2 = R2)

    # Estimate beta from n samples
    if (est.beta) {

      beta.est <- fda.usc::fregre.pc.cv(fdataobj = samp$X.fdata, y = samp$Y,
                                        kmax = 10, criteria = "SICc")$fregre.pc$beta.est

    }
    if (times) {

      cat("Model", k, "finished in", proc.time()[3] - t, "secs\n")

    }

    # Create two axis plot
    lylim <- range(samp$beta.fdata$data) + c(-1, 1)
    rylim <- range(samp$X.fdata$data) + c(-1, 1)
    plotrix::twoord.plot(lx = samp$X.fdata$argvals, ly = samp$beta.fdata$data,
                         rx = samp$X.fdata$argvals, ry = samp$X.fdata$data[1, ],
                         xlab = "", rylab = "", ylab = "", lylim = lylim,
                         rylim = rylim, type = "n", lcol = 1, rcol = gray(0.5),
                         mar = mar, main = ifelse(main, paste("Scenario", k), ""))

    # Theoretical beta
    lines(samp$beta.fdata, ylim = lylim, lwd = 2, col = 1)

    # Plot beta estimate
    if (est.beta) {

      lines(beta.est, ylim = lylim, lwd = 2, col = 2)

    }

    # Add functional data
    for (i in 1:20) {

      lines(diff(lylim) / diff(rylim) * (samp$X.fdata[i, ] - rylim[1]) + lylim[1],
            col = gray(0.5))

    }

  }

}


# \code{check.scenarios} produces plots displaying th densities of the response 
# for different deviations. \code{check.betas} displays a sample of the functional process, 
# the functional coefficient associated and its estimation based on 100 observations.

# Check scenarios and deviations
 set.seed(3257641)
 check.scenarios(scenarios = c(1, 7, 3, 5, 6, 4, 8, 9, 12), composite = TRUE)
    
 # Check betas
 set.seed(3257641)
 require(plotrix)
 check.betas(scenarios = c(1, 7, 3, 5, 6, 4, 8, 9, 12), n = 100)
    