library(compiler)
library(Rcpp)

############ SIMULATION ############
# simule un GARCH(1,1)-CCC diagonal avec bruit Student de variance R

GarchCCC.sim <- function(n, omega, alpha, beta, model,R,noise, nu = Inf, valinit = 500) {
    m <- length(omega)
    cst <- 1
    #Generation des bruits
    if (noise == "normal")
     {
        eta <- matrix(rnorm(m * (n + valinit)), ncol = m, nrow = (n + valinit)) %*% t(Sqrt(R))
    }
    else
    {

        if (noise == "student") {
            if (nu > 2 & nu != Inf)
                cst <- sqrt((nu - 2) / nu)
            eta <- matrix(cst * rt(m * (n + valinit), nu), ncol = m, nrow = (n + valinit)) %*% t(Sqrt(R))
        }
            else {
                print("Not a valid noise")
            }
    }
    eta <- matrix(rnorm(m * (n + valinit)), ncol = m, nrow = (n + valinit)) %*% t(Sqrt(R))
    eps <- matrix(0, nrow = n + valinit, ncol = m)
    ht <- matrix(0, nrow = n + valinit, ncol = m)
    if (model == "diagonal") {
        for (t in 2:(n + valinit)) {
            ht[t,] <- omega + alpha * (eps[t - 1,] ^ 2) + beta * ht[t - 1,]
            eps[t,] <- sqrt(ht[t,]) * eta[t,]
        }
    }
    else {
        if (model == "sdiagonal") {
            for (t in 2:(n + valinit)) {
                ht[t,] <- omega + alpha %*% (eps[t - 1,] ^ 2) + beta * ht[t - 1,]
                eps[t,] <- sqrt(ht[t,]) * eta[t,]
            }
        }
        else {
            print("Not a valid model")

            }

    }
    return(as.data.frame(eps[(valinit + 1):(n + valinit),]))
}

GarchCCC.sim <- cmpfun(GarchCCC.sim)



############ ESTIMATION ############

#Modele diagonal
objf <- function(x, eps, eps2, n, r, tol = sqrt(.Machine$double.eps)) {
    omega <- x[1]
    alpha <- x[2]
    beta <- x[3]
    sigma2 <- rep(0, n)
    sigma2[1] <- var(eps[1:r])
    sigma2.mod <- objfcpp(eps2, sigma2, omega, alpha, beta, r, n, tol) #Appel du c++
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
    sigma2 <- objfcpp2(eps2, sigma2, omega, alpha, beta, n)   #Appel du c++

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
    sigma2 <- rep(0,n)
    for (t in 2:n) {
        sigma2[t] <- omega + as.numeric(t(Alpha) %*% Ht[t - 1,]) + beta * sigma2[t - 1]
    }
    eta <- eps0[(r + 1):n] / sqrt(sigma2[(r + 1):n])
    list(coef = res$par, minimum = res$objective, eta = eta) # sd = sd,)

}




#estime EbE un MGARCH(1,1)-CCC diagonal ou semi-diagonal

estimCCC.EbEE <- function(Omega, Alpha, Beta, eps, r = 10, model) {
    if (model == "diagonal")
        {
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
        list(Omega=Omega.est, A=Alpha.est, B=Beta.est, R=R,Residuals=Res)
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

                Res <- matrix(0, nrow = (n - r), ncol = nbind)
                for (i in 1:nbind) {
                    QMLE <- estimMgarch.sdiag(Omega[i], Alpha[i,], Beta[i], eps[, i], Ht, r)
                    Omega[i] <- QMLE$coef[1]
                    A[i,] <- QMLE$coef[2:(d + 1)]
                    B[i] <- QMLE$coef[d + 2]

                    Res[, i] <- QMLE$eta
                }
                R <- cor(Res)

                list(Omega = Omega, A = A, B = B, minimum = QMLE$minimum, R = R,Residuals=Res) 
            }
                else {
                    print("Not a valid model")
                }

        }


}


MSD.CCC.EbEE <- function(theta0, init, nobs, iter, type, noise, nu = Inf) {
    Omega.list <- list()
    A.list <- list()
    B.list <- list()
    R.list <- list()


    for (i in 1:iter) {
        eps <- GarchCCC.sim(nobs, theta0$Omega, theta0$A, theta0$B, theta0$R, nu = nu, model = type, noise = noise)
        tmp <- estimCCC.EbEE(init$Omega, init$A, init$B, eps, r = 10, model = type)
        Omega.list[[i]] <- tmp$Omega
        A.list[[i]] <- tmp$A
        B.list[[i]] <- tmp$B
        R.list[[i]] <- tmp$R
    }
    Omega.mean <- apply(simplify2array(Omega.list), 1, mean)
    Omega.sd <- apply(simplify2array(Omega.list), 1, sd)
    if (type == "sdiagonal") {
        A.mean <- apply(simplify2array(A.list), 1:2, mean)
        A.sd <- apply(simplify2array(A.list), 1:2, sd)
    }
        else {
            if (type == "diagonal") {
                A.mean <- apply(simplify2array(A.list), 1, mean)
                A.sd <- apply(simplify2array(A.list), 1, sd)
            } else {
                print("Not a valid model")
            }
        }
    B.mean <- apply(simplify2array(B.list), 1, mean)
    B.sd <- apply(simplify2array(B.list), 1, sd)
    R.mean <- apply(simplify2array(R.list), 1:2, mean)
    R.sd <- apply(simplify2array(R.list), 1:2, sd)

    return(list(Omega.mean = Omega.mean, Omega.sd = Omega.sd, A.mean = A.mean, A.sd = A.sd, B.mean = B.mean, B.sd = B.sd, R.mean = R.mean, R.sd = R.sd))

}
