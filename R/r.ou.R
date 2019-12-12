#' Ornstein-Uhlenbeck process
#' 
#' Sampling of paths of the Ornstein-Uhlenbeck process.
#' 
#' @param n number of curves.
#' @param t discretization points.
#' @param mu mean of the process.
#' @param alpha strength of the drift.
#' @param sigma diffusion coefficient.
#' @param x0 a number or a vector of length \code{n} giving the initial
#' value(s) of the Ornstein-Uhlenbeck process. By default, \code{n} points are
#' sampled from the stationary distribution.
#' @return Functional sample, an \code{\link{fdata}} object of length
#' \code{n}.
#' @author Eduardo Garcia-Portugues (\email{edgarcia@@est-econ.uc3m.es}).
#' @examples
#' plot(r.ou(n = 100))
#' plot(r.ou(n = 100, alpha = 2, sigma = 4, x0 = 1:100))
#' 
#' @export r.ou
r.ou <- function(n, t = seq(0, 1, len = 201), mu = 0, alpha = 1, sigma = 1, 
                 x0 = rnorm(n, mean = mu, sd = sigma/sqrt(2 * alpha))){

  # Time-varying covariances
  St <- sigma^2/(2 * alpha) * outer(t, t, function(s, t) {
    exp(alpha * (2 * pmin(s, t) - (s + t))) - exp(-alpha * (s + t))
    })
  
  # Sample N(0, St) and add time-varying mean
  #mdata <- mvtnorm::rmvnorm(n = n, mean = rep(0, length(t)), sigma = St)
  mdata <- mvrnorm(n = n, mu = rep(0, length(t)), Sigma = St)
  
  
  mdata <- mdata + outer(x0, t, function(x0, t) mu + (x0 - mu) * exp(-alpha * t))
  
  # As fdata object
  return(fdata(mdata = mdata, argvals = t))
  
}
