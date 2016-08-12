library(compiler)

############ SIMULATION ############

# Simulation d'un GARCH(1,1) DCC, beta diagonal de Aielli ou de Engle, bruit student ou normal
GarchDCC.sim <- function(n, Omega, A, B, alpha, beta, S, nu = Inf, valinit = 500, model, noise) {
  tol = sqrt(.Machine$double.eps)
  m <- length(Omega)
  Qt <- array(dim = c(n + valinit, m, m))
  Qt[1,,] <- S
  Qt.star <- diag(sqrt(diag(Qt[1,,])))
  Qt.star.inv <- diag(sqrt(1 / pmax(diag(Qt[1,,]), tol)))
  Rt <- Qt
  Rt[1,,] <- Qt.star.inv %*% Qt[1,,] %*% Qt.star.inv
  cst <- 1

  #Génération des bruits
  if (noise == "student")
  {
    if (nu > 2 & nu != Inf)
      cst <- sqrt((nu - 2) / nu)
    eta <- matrix(cst * rt(m * (n + valinit), nu), nrow = (n + valinit), ncol = m)
  }
  else
  {
    if (noise == "normal")
    {
      eta <- matrix(rnorm(m * (n + valinit)), nrow = (n + valinit), ncol = m)
    }
    else
    {
      print("Not a valid noise")
    }
  }

  #Mise en place des recurrences
  eps <- matrix(0, nrow = (n + valinit), ncol = m)
  eta.star <- matrix(0, nrow = (n + valinit), ncol = m)
  ht <- matrix(0, nrow = (n + valinit), ncol = m)
  eta.star[1,] <- eta[1,] %*% t(Sqrt(Rt[1,,]))
  for (t in 2:(n + valinit)) {
    ht[t,] <- Omega + as.vector(A %*% (eps[t - 1,] ^ 2)) + B * ht[t - 1,]

    {
      if (model=="Aielli") {
        Qt[t,,] <- (1 - alpha - beta) * S + alpha * Qt.star %*% eta.star[t - 1,] %*% t(eta.star[t - 1,]) %*% Qt.star + beta * Qt[t - 1,,]
      } else {
        if (model == "Engle") {
          Qt[t,,] <- (1 - alpha - beta) * S + alpha * eta.star[t - 1,] %*% t(eta.star[t - 1,]) + beta * Qt[t - 1,,]
        }
        else {print("Not a valid model") }
      }
    }
    Qt.star <- diag(sqrt(diag(Qt[t,,])))
    Qt.star.inv <- diag(sqrt(1 / pmax(diag(Qt[t,,]), tol)))
    Rt[t,,] <- Qt.star.inv %*% Qt[t,,] %*% Qt.star.inv
    eta.star[t,] <- eta[t,] %*% t(Sqrt(Rt[t,,]))
    eps[t,] <- sqrt(ht[t,]) * eta.star[t,]
  }

  return(list(sim = eps[(valinit + 1):(n + valinit),], cor = Rt[(valinit + 1):(n + valinit),,]))
}

GarchDCC.sim <- cmpfun(GarchDCC.sim)


############ ESTIMATION ############


# estimeur de second etape d'un DCC

objf.Rt.DCC <- function(x, eta.star, m, n, r, tol = sqrt(.Machine$double.eps),type) {
tol = sqrt(.Machine$double.eps)
S <- inv.vech0(x[1:(m * (m - 1) / 2)])

aalpha <- x[1 + (m * (m - 1) / 2)]
bbeta <- x[2 + (m * (m - 1) / 2)]

if (!all(is.finite(x)) | (aalpha + bbeta) > 1 - tol) {
    qml <- Inf
} else {

    Qt <- array(dim = c(n, m, m))
    Qt[1,,] <- S
    Qt.star <- diag(sqrt(diag(Qt[1,,])))
    Qt.star.inv <- diag(sqrt(1 / pmax(diag(Qt[1,,]), tol)))
    Rt <- Qt
    Rt[1,,] <- Qt.star.inv %*% Qt[1,,] %*% Qt.star.inv

    l <- rep(0, n)

    for (t in 2:n) {
        if (type == "Aielli") {
            Qt[t,,] <- (1 - aalpha - bbeta) * S + aalpha * Qt.star %*% eta.star[t - 1,] %*% t(eta.star[t - 1,]) %*% Qt.star + bbeta * Qt[t - 1,,]
        } else {
            if (type == "Engle") {
                Qt[t,,] <- (1 - aalpha - bbeta) * S + aalpha * eta.star[t - 1,] %*% t(eta.star[t - 1,]) + bbeta * Qt[t - 1,,]
            } else {
                print("Not a valid model")
            }
        }
        Qt.star <- diag(sqrt(diag(Qt[t,,])))
        Qt.star.inv <- diag(sqrt(1 / pmax(diag(Qt[t,,]), tol)))
        Rt[t,,] <- Qt.star.inv %*% Qt[t,,] %*% Qt.star.inv
        if (kappa(Rt[t, , ]) < 1 / tol) {
            Rt.inv <- solve(Rt[t,,])
            l[t] <- as.numeric(eta.star[t,] %*% Rt.inv %*% eta.star[t,]) + log(det(Rt[t,,]))
        } else {
            l[t] <- Inf
        }
    }
    qml <- mean(l[(r + 1):n])
}}



objf.Rt.DCC <- cmpfun(objf.Rt.DCC)

#


estim.Rt.DCC <- function(S, aalpha, bbeta, eta.star, r = 10,type) {
    n <- length(eta.star[, 1])
    m <- length(eta.star[1,])

    valinit <- c(vech0(S), aalpha, bbeta)
    res <- nlminb(valinit, objf.Rt.DCC, lower = c(rep(-1 + 0.01, (m * (m - 1) / 2)), rep(0, 2)),
    upper = c(rep(1 - 0.01, (m * (m - 1) / 2)), rep(1, 2)), eta.star = eta.star, m = m, n = n, r = r,type=type)
    S.est <- inv.vech0(res$par[1:(m * (m - 1) / 2)])
    aalpha.est <- res$par[1 + (m * (m - 1) / 2)]
    bbeta.est <- res$par[2 + (m * (m - 1) / 2)]
    return(list(S=S.est, alpha=aalpha.est, beta=bbeta.est))

}

estim.Rt.DCC <- cmpfun(estim.Rt.DCC)


estimDCC.EbEE <- function(Omega, A, B, S, alpha, beta, eps, r = 10, type) {
    m <- length(Omega)
    n <- length(eps[, 1])
    eps2 <- eps ^ 2

    #First step
    EbEE <- estimCCC.EbEE(Omega, A, B, eps, r = 10, model = "sdiagonal")
    omega.est <- EbEE$Omega
    Alpha.est <- EbEE$A
    Beta.est <- EbEE$B
    Res <- EbEE$R
    Residuals <- EbEE$Residuals

    #2nd step
    valinit <- c(vech0(cor(Res)), alpha, beta)
    n <- nrow(Residuals)
    step <- estim.Rt.DCC(S, alpha, beta, Residuals, r = 10,type=type)
    S.est <- step$S
    aalpha.est <- step$alpha
    bbeta.est <- step$beta

    #3rd step
    omega <- 1-(aalpha.est+bbeta.est)
    Qdiag <- matrix(nrow = n, ncol = m)
    Sstar <- array(dim = c(n, m, m))
    Qdiag[1,] <- 1
    demi.Sstar <- sqrt(Qdiag[1,]) * Residuals[1,]
    Sstar[1,,] <- demi.Sstar %*% t(demi.Sstar)

    for (i in 2:n) {
        Qdiag[i,] <- omega + (aalpha.est * (Residuals[i - 1,] ^ 2) + bbeta.est) * Qdiag[i - 1,]
        demi.Sstar <- sqrt(Qdiag[i,]) * Residuals[i,]
        Sstar[i,,] <- demi.Sstar %*% t(demi.Sstar)
    }

    S.est <- matrix(0, nrow = m, ncol = m)
    for (i in (r + 1):n) {
        S.est <- S.est + Sstar[i,,] / n
    }
    S.star.moinsdemi <- diag(1 / sqrt(diag(S.est)))
    S.est <- S.star.moinsdemi %*% S.est %*% S.star.moinsdemi

    return(list(Omega=omega.est, A=Alpha.est, B=Beta.est, S=S.est, alpha=aalpha.est, beta=bbeta.est))

}

estimDCC.EbEE <- cmpfun(estimDCC.EbEE)


MSD.DCC.EbEE <- function(theta0, init, nobs, iter, type, noise, nu=Inf) {
    Omega.list <- list()
    A.list <- list()
    B.list <- list()
    S.list <- list()
    alpha.list <- c()
    beta.list <- c()

    for (i in 1:iter) {
        eps <- GarchDCC.sim(nobs, theta0$Omega, theta0$A, theta0$B, theta0$alpha, theta0$beta, theta0$S, nu = nu, model = type, noise=noise)
        tmp <- estimDCC.EbEE(init$Omega, init$A, init$B, init$S, init$alpha, init$beta, eps$sim, r = 10, type = type)
        Omega.list[[i]] <- tmp$Omega
        A.list[[i]] <- tmp$A
        B.list[[i]] <- tmp$B
        S.list[[i]] <- tmp$S
        alpha.list <- c(alpha.list,tmp$alpha)
        beta.list <- c(beta.list,tmp$beta)
    }

    Omega.mean <- apply(simplify2array(Omega.list), 1, mean)
    Omega.sd <- apply(simplify2array(Omega.list), 1, sd)
    A.mean <- apply(simplify2array(A.list), 1:2, mean)
    A.sd <- apply(simplify2array(A.list), 1:2, sd)
    B.mean <- apply(simplify2array(B.list), 1, mean)
    B.sd <- apply(simplify2array(B.list), 1, sd)
    S.mean <- apply(simplify2array(S.list), 1:2, mean)
    S.sd <- apply(simplify2array(S.list), 1:2, sd)
    alpha.mean <- mean(alpha.list)
    alpha.sd <- sd(alpha.list)
    beta.mean <- mean(beta.list)
    beta.sd <- sd(beta.list)

    return(list(Omega.mean = Omega.mean, Omega.sd=Omega.sd,A.mean=A.mean,A.sd=A.sd,B.mean=B.mean,B.sd=B.sd,S.mean=S.mean,S.sd=S.sd,alpha.mean=alpha.mean,alpha.sd=alpha.sd,beta.mean=beta.mean,beta.sd=beta.sd))

}


































