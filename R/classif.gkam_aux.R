kgam.H <-function (object, inverse = "svd") { 
  lenH <- length(object)
  if (lenH == 1) 
    return(object[[1]]$H)
  else {
    SS.list <- SS2 <- list()
    n <- ncol(object[[1]]$H)
    SS <- matrix(NA, ncol = (lenH) * n, nrow = (lenH) * n)
    II = diag(n)
    unos <- matrix(1, ncol = n, nrow = n)/n
    M <- object[[1]]$H
    MM = sweep(M, 1, rowMeans(M), "-")
    if (lenH > 1) {
      for (i in 2:lenH) {
        MMaux = sweep(object[[i]]$H, 1, rowMeans(object[[i]]$H), 
                      "-")
        MM <- rbind(MM, MMaux)
      }
    }
    MM1 <- matrix(rep(MM, lenH), nrow = (n * lenH))
    DD <- kronecker(diag(lenH), outer(rep(1, n), rep(1, n)))
    D1 <- abs(DD - 1)
    SS <- MM1 * D1 + diag(n * lenH)
    slv <- try(solve(SS), silent = TRUE)
    if (is(slv,"try-error")) {
      sv <- svd(SS)
      slv <- drop((sv$v %*% diag(1/sv$d) %*% t(sv$u)))
      warning("Inverse of sigma computed by SVD")
    }
    SSinv <- switch(inverse, solve = slv, svd = {
      res.X = svd(SS)
      (res.X$v) %*% diag(res.X$d^(-1)) %*% t(res.X$u)
    })
    H2 = SSinv %*% MM
    HH <- unos + H2[1:n, 1:n]
    if (lenH > 1) {
      for (i in 2:lenH) {
        HH = HH + H2[(n * (i - 1) + 1):(n * i), ]
      }
    }
    HH
  }
}