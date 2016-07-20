\name{VaR Spherical}
\alias{VaR.Spherical}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute CVaR when residuals have a spherical distribution
%%  ~~function to do ... ~~
}
\description{Compute conditional VaR of a portfolio where yields follow a MGARCH(1,1) DCC model, for spherical distribution of the residuals
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VaR.Spherical(n, Omega, A, B, alpha, beta, S, eps, type, level, weights)
}
%- maybe also 'usage' for other objects documented here.
\arguments{With usual notations
  \item{n}{Number of observations used to compute the quantile of the distribution}
  \item{Omega}{Estimated parameter for Omega}
  \item{A}{Estimated parameter for A}
  \item{B}{Estimated parameter for B}
  \item{alpha}{Estimated parameter for alpha}
  \item{beta}{Estimated parameter for beta}
  \item{S}{Estimated parameter for S}
  \item{eps}{Data}
  \item{type}{type="Engle" for a Engle-DCC
\cr type="Aielli" for an Aielli-DCC}
  \item{level}{Level of the CVaR}
  \item{weights}{Matrix of the weights of the asset (time dependant)}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{VaR : vector with the estimated CVaR
}
\references{C. Francq & J.M. Zakoian, Joint inference on market and estimation risks in dynamic portfolios
%% ~put references to the literature/web site here ~
}
\author{D. Taouss & C. Francq
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{EbEEMGARCH}}{  Homepage of the documentation}}
%% ~~objects to See Also as \code{\link{help}}, ~~~

\examples{
#####
#Simulation of the yield of 2 assets
#####

m<-2
n <- 800
Omega <- c(0.001, 0.001);
A <- matrix(c(0.03, 0.01, 0.01, 0.03), nrow = 2)
B <- c(0.1, 0.1);
S <- matrix(c(1, 0.4, 0.4, 1), nrow = 2)
alpha <- 0.05;
beta <- 0.97 - alpha
nu <- 14

yield <- GarchDCC.sim(n, Omega, A, B, alpha, beta, S, nu = nu, noise = "student", model = "Aielli")

#####
#Creation of a static equiponderate portfolio
#####

weights <- matrix(0.5,ncol=m,nrow=n)
yield_p <- c()
for (t in 1:n) {
   yield_p<-c(yield_p,weights[t,1]*yield$sim[t,1]+weights[t,2]*yield$sim[t,2])
}

#####
#Estimation of the parameters (3-steps method) & the residuals
#####

EbEE<-estimDCC.EbEE(Omega,A,B,S,alpha,beta,yield$sim,type="Aielli")
res <- residuals_DCC(EbEE$Omega, EbEE$A, EbEE$B, EbEE$alpha, EbEE$beta, EbEE$S, yield$sim, type="Aielli")

VaR <- VaR.Spherical(700, EbEE$Omega, EbEE$A, EbEE$B, EbEE$alpha, EbEE$beta, EbEE$S, yield$sim, type="Aielli", level=0.01, weights=weights)


#Plot
#plot(yield_p,xlim=c(700,800),type="l")
#lines(VaR,col='red')

}
