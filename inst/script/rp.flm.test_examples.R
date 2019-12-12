# Real data examples of rp.flm.test function
# library(fda.usc)
  rp.flm.test(X.fdata = X, Y = Y, est.method = "pls")
  rp.flm.test(X.fdata = X, Y = Y, est.method = "basis", 
              p.criterion = fda.usc::GCV.S)
  rp.flm.test(X.fdata = X, Y = Y, est.method = "pc", p = 5)
  rp.flm.test(X.fdata = X, Y = Y, est.method = "pls", p = 5)
  rp.flm.test(X.fdata = X, Y = Y, est.method = "basis", p = 5)
  rp.flm.test(X.fdata = X, Y = Y, beta0.fdata = beta0)
  
  # Composite hypothesis: do not reject FLM
  rp.test <- rp.flm.test(X.fdata = X, Y = Y, est.method = "pc")
  rp.test$p.values.fdr
  pcvm.test <- flm.test(X.fdata = X, Y = Y, est.method = "pc", B = 1e3,
                        plot.it = FALSE)
  pcvm.test
  
  # Estimation of beta
  par(mfrow = c(1, 3))
  plot(X, main = "X")
  plot(beta0, main = "beta")
  lines(rp.test$beta.est, col = 2)
  lines(pcvm.test$beta.est, col = 3)
  plot(density(Y), main = "Density of Y", xlab = "Y", ylab = "Density")
  rug(Y)
  
  # Simple hypothesis: do not reject beta = beta0
  rp.flm.test(X.fdata = X, Y = Y, beta0.fdata = beta0)$p.values.fdr
  flm.test(X.fdata = X, Y = Y, beta0.fdata = beta0, B = 1e3, plot.it = FALSE)
  
  # Simple hypothesis: reject beta = beta0^2
  rp.flm.test(X.fdata = X, Y = Y, beta0.fdata = beta0^2)$p.values.fdr
  flm.test(X.fdata = X, Y = Y, beta0.fdata = beta0^2, B = 1e3, plot.it = FALSE)
  
  # Tecator dataset
  
  # Load data
  data(tecator)
  absorp <- tecator$absorp.fdata
  ind <- 1:129 # or ind <- 1:215
  x <- absorp[ind, ]
  y <- tecator$y$Fat[ind]
  
  # Composite hypothesis
  rp.tecat <- rp.flm.test(X.fdata = x, Y = y, est.method = "pc")
  pcvm.tecat <- flm.test(X.fdata = x, Y = y, est.method = "pc", B = 1e3,
                         plot.it = FALSE)
  rp.tecat$p.values.fdr[c(5, 10), ]
  pcvm.tecat
  
  # Simple hypothesis
  zero <- fdata(mdata = rep(0, length(x$argvals)), argvals = x$argvals,
                rangeval = x$rangeval)
  rp.flm.test(X.fdata = x, Y = y, beta0.fdata = zero)
  flm.test(X.fdata = x, Y = y, beta0.fdata = zero, B = 1e3)
  
  # With derivatives
  rp.tecat <- rp.flm.test(X.fdata = fdata.deriv(x, 1), Y = y, est.method = "pc")
  rp.tecat$p.values.fdr
  rp.tecat <- rp.flm.test(X.fdata = fdata.deriv(x, 2), Y = y, est.method = "pc")
  rp.tecat$p.values.fdr
  
  # AEMET dataset
  
  # Load data
  data(aemet)
  wind.speed <- apply(aemet$wind.speed$data, 1, mean)
  temp <- aemet$temp
  
  # Remove the 5\% of the curves with less depth (i.e. 4 curves)
  par(mfrow = c(1, 1))
  res.FM <- depth.FM(temp, draw = TRUE)
  qu <- quantile(res.FM$dep, prob = 0.05)
  l <- which(res.FM$dep <= qu)
  lines(aemet$temp[l], col = 3)
  
  # Data without outliers
  wind.speed <- wind.speed[-l]
  temp <- temp[-l]
  
  # Composite hypothesis
  rp.aemet <- rp.flm.test(X.fdata = temp, Y = wind.speed, est.method = "pc")
  pcvm.aemet <- flm.test(X.fdata = temp, Y = wind.speed, B = 1e3,
                         est.method = "pc", plot.it = FALSE)
  rp.aemet$p.values.fdr
  apply(rp.aemet$p.values.fdr, 2, range)
  pcvm.aemet
  
  # Simple hypothesis
  zero <- fdata(mdata = rep(0, length(temp$argvals)), argvals = temp$argvals,
                rangeval = temp$rangeval)
  flm.test(X.fdata = temp, Y = wind.speed, beta0.fdata = zero, B = 1e3,
           plot.it = FALSE)
  rp.flm.test(X.fdata = temp, Y = wind.speed, beta0.fdata = zero)
  