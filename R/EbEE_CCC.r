library(compiler)
library(Rcpp)

############ FONCTIONS ANNEXES ############
#Fonction c++ pour accéléer

#racine carree d'un matrice symetrique semi-definie positive
Sqrt <- function(Sigma) {
    n <- nrow(Sigma)
    sp <- eigen(Sigma)
    sp$vector %*% sqrt(diag(sp$values))
}
Sqrt <- cmpfun(Sqrt)

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

############ SIMULATION ############
# simule un GARCH(1,1)-CCC diagonal avec bruit Student de variance R.mat
Mgarch.sim <- function(n, omega, alpha, beta, R.mat, nu = Inf, valinit = 500) {
    m <- length(omega)
    cst <- 1
    if (nu > 2 & nu != Inf)
        cst <- sqrt((nu - 2) / nu)
    # eta <- matrix(cst*rt(m*(n+valinit),nu),ncol=m,nrow=(n+valinit))%*%t(Sqrt(R.mat))
    eta <- matrix(rnorm(m * (n + valinit)), ncol = m, nrow = (n + valinit)) %*% t(Sqrt(R.mat))
    eps <- matrix(0, nrow = n + valinit, ncol = m)
    ht <- matrix(0, nrow = n + valinit, ncol = m)

    for (t in 2:(n + valinit)) {
        ht[t,] <- omega + alpha * (eps[t - 1,] ^ 2) + beta * ht[t - 1,]
        eps[t,] <- sqrt(ht[t,]) * eta[t,]
    }

    return(as.data.frame(eps[(valinit + 1):(n + valinit),]))
}

Mgarch.sim <- cmpfun(Mgarch.sim)


m <- 3
Omega0 <- rep(0.01, m)
Alpha0 <- rep(0.05, m)
Beta0 <- rep(0.90, m)
R0 <- diag(rep(1, m))

Epsi <- Mgarch.sim(n = 2500, omega = Omega0, alpha = Alpha0, beta = Beta0, R.mat = R0, nu = Inf, valinit = 500)
Epsi <- as.data.frame(Epsi)


############ ESTIMATION ############

#Modele diagonal
objf <- function(x, eps, eps2, n, r, tol = sqrt(.Machine$double.eps)) {
    omega <- x[1]
    alpha <- x[2]
    beta <- x[3]
    sigma2 <- rep(0, n)
    sigma2[1] <- var(eps[1:r])
    sigma2.mod <- objfcpp(eps2, sigma2, omega, alpha, beta, r, n, tol)
    sigma2.mod <- sigma2.mod[(r + 1):n]
    qml <- mean(log(sigma2.mod) + eps[(r + 1):n] ** 2 / sigma2.mod)
    qml
}

objf <- cmpfun(objf)

estim.1equa.Mgarch11 <- function(omega, alpha, beta, eps, r) {
    n <- length(eps)
    eps2 <- eps ^ 2
    valinit <- c(omega, alpha, beta)
    res <- nlminb(valinit, objf, lower = c(0.0000000000001, rep(0, 2)),
                upper = c(Inf, rep(0.9999, 2)), eps = eps, eps2 = eps2, n = n, r = r)
    omega <- res$par[1]
    alpha <- res$par[2]
    beta <- res$par[3]

    sigma2 <- rep(0, n)
    sigma2[1] <- var(eps[1:r])
    sigma2 <- objfcpp2(eps2, sigma2, omega, alpha, beta, n)

    eta <- eps[(r + 1):n] / sqrt(sigma2[(r + 1):n])

    list(coef = res$par, minimum = res$objective, eta = eta)

}

estim.1equa.Mgarch11 <- cmpfun(estim.1equa.Mgarch11)

#Modele semi-diagonal
objf.sdiag <- function(x, eps, Ht, n, d, r, tol = sqrt(.Machine$double.eps)) {
    omega <- x[1]
    Alpha <- x[2:(d + 1)]
    beta <- x[d + 2]
    sigma2 <- rep(0, n)
    for (t in 2:n)
        sigma2[t] <- omega + as.numeric(t(Alpha) %*% Ht[t - 1,]) + beta * sigma2[t - 1]
    sigma2.mod <- pmax(sigma2[(r + 1):n], tol)
    qml <- mean(log(sigma2.mod) + eps[(r + 1):n] ** 2 / sigma2.mod)
    qml
}

estimMgarch.sdiag <- function(omega, Alpha, beta, eps0, Ht, r) {
    n <- length(eps0)
    d <- ncol(Ht)
    valinit <- c(omega, Alpha, beta)
    res <- nlminb(valinit, objf.sdiag, lower = c(0.0000000000001, rep(0, (d + 1))),
    upper = c(Inf, rep(0.9999, (d + 1))), eps = eps0, Ht = Ht, n = n, d = d, r = r)
    omega <- res$par[1]
    Alpha <- res$par[2:(d + 1)]
    beta <- res$par[d + 2]
    #var <- Pgarch.var(omega, Alpha, beta, eps0, Ht, r)
    #sd <- var$sd
    sigma2 <- rep(0,n)
    for (t in 2:n) {
        sigma2[t] <- omega + as.numeric(t(Alpha) %*% Ht[t - 1,]) + beta * sigma2[t - 1]
    }
    eta <- eps0[(r + 1):n] / sqrt(sigma2[(r + 1):n])
    list(coef = res$par, minimum = res$objective, eta = eta) # sd = sd,)

}




# estime EbE un MGARCH(1,1)-CCC diagonal 

estim.EbEE <- function(Omega, Alpha, Beta, eps, r = 10, model) {
    if (model == "diagonal") 
        {
        #fast()
        m <- length(Omega)
        n <- length(eps[, 1])
        Res <- matrix(nrow = (n - r), ncol = m)
        Omega.est <- Omega
        Alpha.est <- Alpha
        Beta.est <- Beta
        for (j in 1:m)
            {
                res <- estim.1equa.Mgarch11(Omega[j], Alpha[j], Beta[j], eps[, j], r)
                Omega.est[j] <- res$coef[1]
                Alpha.est[j] <- res$coef[2]
                Beta.est[j] <- res$coef[3]
                Res[, j] <- res$eta
            }
        R <- cor(Res)
        list(Omega=Omega.est, Alpha=Alpha.est, Beta=Beta.est, R=R)
        }

        else 
            {

            if (model == "sdiagonal") {
                nbind <- ncol(eps)
                n <- nrow(eps)
                Ht <- eps ^ 2
                d <- ncol(Ht)
                Omega <- rep(0, nbind)
                B <- rep(0, nbind)
                A <- diag(rep(0, nbind))
                # SD <- matrix(0, nrow = nbind, ncol = (nbind + 2))
                Res <- matrix(0, nrow = (n - r), ncol = nbind)
                for (i in 1:nbind) {
                    QMLE <- estimMgarch.sdiag(Omega[i], Alpha[i,], Beta[i], eps[, i], Ht, r)
                    Omega[i] <- QMLE$coef[1]
                    A[i,] <- QMLE$coef[2:(d + 1)]
                    B[i] <- QMLE$coef[d + 2]
                    #SD[i,] <- QMLE$sd
                    Res[, i] <- QMLE$eta
                }
                R <- cor(Res)
                #        var <- Mgarch.var(Omega, A, B, mat.rend, r)
                list(Omega = Omega, A = A, B = B, minimum = QMLE$minimum, R = R) #SD = SD,, sd.R = var$sd)
            }
                else {
                    print("Modele non reconnu")
                }
        
        }


}


