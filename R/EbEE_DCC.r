library(compiler)

############ SIMULATION ############

# Simulation d'un GARCH(1,1) DCC, beta diagonal de Aielli ou de Engle, bruit student ou normal
GarchDCC.sim <- function(n, omega, Alpha, beta, aalpha, bbeta, S, nu = Inf, valinit = 500, Aielli = TRUE, noise) {
    tol = sqrt(.Machine$double.eps)
    m <- length(omega)
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
        ht[t,] <- omega + as.vector(Alpha %*% (eps[t - 1,] ^ 2)) + beta * ht[t - 1,]
        
        {
            if (Aielli) {
                Qt[t,,] <- (1 - aalpha - bbeta) * S + aalpha * Qt.star %*% eta.star[t - 1,] %*% t(eta.star[t - 1,]) %*% Qt.star + bbeta * Qt[t - 1,,]
            } else {
                Qt[t,,] <- (1 - aalpha - bbeta) * S + aalpha * eta.star[t - 1,] %*% t(eta.star[t - 1,]) + bbeta * Qt[t - 1,,]
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


#Lazy data
omega <- c(0.01, 0.01);
Alpha <- matrix(c(0.03, 0.01, 0.01, 0.03), nrow = 2)
beta <- c(0.8, 0.8);
S <- matrix(c(1, 0.4, 0.4, 1), nrow = 2)
aalpha <- 0.01;
bbeta <- 0.99 - aalpha
n <- 800
nu <-7 
res <- GarchDCC.sim(n, omega, Alpha, beta, aalpha, bbeta, S, nu = nu, noise = "student")





















































