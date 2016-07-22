library(compiler)

#Square root of a symetric semi-definite positive matrix

Sqrt <- function(Sigma) {
  n <- nrow(Sigma)
  sp <- eigen(Sigma)
  sp$vector %*% sqrt(diag(sp$values))
}
Sqrt <- cmpfun(Sqrt)

InvSqrt <- function(Sigma, tol = sqrt(.Machine$double.eps)) {
  n <- nrow(Sigma)
  sp <- eigen(Sigma)

  res<-sp$vector %*% sqrt(diag(1 / pmax(abs(sp$values), tol))) %*% t(sp$vector)
  return(res)
}

#vech0
vech0 <- function(A) {
  d <- nrow(A)
  res <- rep(0, (d * (d - 1) / 2))
  ind <- 0
  for (j in 1:(d - 1)) {
    for (i in (j + 1):d) {
      ind <- ind + 1
      res[ind] <- A[i, j]
    }
  }
  res
}

inv.vech0 <- function(rho) {
  m <- length(rho)
  d <- (1 + sqrt(1 + 8 * m)) / 2
  A <- diag(d)
  ind <- 0
  for (j in 1:(d - 1)) {
    for (i in (j + 1):d) {
      ind <- ind + 1
      A[i, j] <- rho[ind]
      A[j, i] <- rho[ind]
    }
  }
  A
}
