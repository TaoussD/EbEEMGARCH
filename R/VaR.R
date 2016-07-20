library(compiler)
InvSqrt <- function(Sigma, tol = sqrt(.Machine$double.eps)) {
  n <- nrow(Sigma)
  sp <- eigen(Sigma)

  res<-sp$vector %*% sqrt(diag(1 / pmax(abs(sp$values), tol))) %*% t(sp$vector)
  return(res)
}


residuals_DCC <-function(Omega,A,B,alpha,beta,S,eps,r=10,type)
{
  n<-nrow(eps)
  m<-ncol(eps)
  ###########
  tol <- sqrt(.Machine$double.eps)
  Qt<-array(dim=c(n,m,m))
  Ht<-Qt
  Qt[1,,]<-S
  Qt.star<-diag(sqrt(diag(Qt[1,,])))
  Qt.star.inv<-diag(sqrt(1/pmax(diag(Qt[1,,]),tol)))
  Rt<-Qt
  Rt[1,,]<-Qt.star.inv%*%Qt[1,,]%*%Qt.star.inv


  eta.star<-matrix(0,nrow=n,ncol=m)
  ht<-matrix(0,nrow=n,ncol=m)
  ht[1,]<-colMeans(eps[1:r,]^2)
  Ht[1,,]<-diag(sqrt(ht[1,]))%*%Rt[1,,]%*%diag(sqrt(ht[1,]))


  for (t in 2:n) {
    ht[t,] <- Omega + as.vector(A %*% (eps[t - 1,] ^ 2)) + B * ht[t - 1,]
    eta.star[t,] <- eps[t,] / sqrt(ht[t,])
    if (type == "Aielli") {
      Qt[t,,] <- (1 - alpha - beta) * S + alpha * Qt.star %*% eta.star[t - 1,] %*% t(eta.star[t - 1,]) %*% Qt.star + beta * Qt[t - 1,,]
    } else {
      if (type == "Engle") {
        Qt[t,,] <- (1 - alpha - beta) * S + alpha * eta.star[t - 1,] %*% t(eta.star[t - 1,]) + beta * Qt[t - 1,,]
      } else {
        print("Not a valid model")
      }
    }

    Qt.star <- diag(sqrt(diag(Qt[t,,])))
    Qt.star.inv <- diag(sqrt(1 / pmax(diag(Qt[t,,]), tol)))
    Rt[t,,] <- Qt.star.inv %*% Qt[t,,] %*% Qt.star.inv
    Ht[t,,] <- diag(sqrt(ht[t,])) %*% Rt[t,,] %*% diag(sqrt(ht[t,]))
  }
  eta <- matrix(0, nrow = n, ncol = m)
  for (t in 1:n) {
    eta[t,] <- InvSqrt(as.matrix(Ht[t,,])) %*% eps[t,]
  }


  return(list(Ht=Ht,Rt=Rt,eta=eta))
}

residuals_DCC <- cmpfun(residuals_DCC)

VaR.Spherical <- function(n, Omega, A, B, alpha, beta, S, eps, type, level, weights) {
  nobs<-nrow(eps)
  residu <- residuals_DCC(Omega, A, B, alpha, beta, S, eps, type = "Aielli")
  VaR <- c()
  quantil <- quantile(abs(residu$eta[1:n]), 1 - 2 * level)
  for (t in 1:nobs) {
    VaR <- c(VaR, - quantil * sqrt(sum((weights[t,] %*% Sqrt(residu$Ht[t,,]) ^ 2))))
  }
  return(VaR)
}

VaR.Spherical <- cmpfun(VaR.Spherical)


