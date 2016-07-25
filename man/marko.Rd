\name{estim.Markowitz}
\alias{estim.Markowitz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute the assets' weights and the associated CVaR using Markowitz methodology

%%  ~~function to do ... ~~
}
\description{Estimation of the assets' weights and the associated CVaR which minimize the variance
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
estim.Markowitz(n, Omega, A, B, alpha, beta, S, type, level, yield)
}
%- maybe also 'usage' for other objects documented here.
\arguments{With usual notations of DCC Garch
  \item{n}{Number of observations used to compute the empirical quantile}
  \item{Omega}{Estimation of Omega}
    \item{A}{Estimated parameter for A}
  \item{B}{Estimated parameter for B}
  \item{alpha}{Estimated parameter for alpha}
  \item{beta}{Estimated parameter for beta}
  \item{S}{Estimated parameter for S}
  \item{eps}{Data}
  \item{type}{type="Engle" for a Engle-DCC
\cr type="Aielli" for an Aielli-DCC}
  \item{level}{Level of the CVaR}
  \item{yield}{Data}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{VaR}{CVaR of the portfolio}
  \item{weights}{Matrix of the weights (time-dependant)}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{C. Francq & J.M. Zakoian, Estimating multivariate GARCH and Stochastic Correlation models equation by equation\cr
C. Francq & J.M. Zakoian, Joint inference on market and estimation risks in dynamic portfolios
%% ~put references to the literature/web site here ~
}
\author{D. Taouss & C. Francq
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{EbEEMGARCH}}{  Homepage of the documentation}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#####
# Sampling some data
#####

n <- 800
Omega <- c(0.001, 0.001, 0.001)
A <-matrix(c(0.03,0.01,0.01,0.01,0.03,0.01,0.01,0.01,0.03),nrow=3)
B <- c(0.1,0.1,0.1)
S <- matrix(c(1,0.4,0.4,0.4,1,0.4,0.4,0.4,1),nrow=3)
alpha <- 0.05;
beta <- 0.97 - alpha
nu <- 14
yield <- GarchDCC.sim(n, Omega, A, B, alpha, beta, S, nu = nu, noise = "student", model = "Aielli")

#####
# Estimation of the parameters
#####

EbEE <- estimDCC.EbEE(Omega, A, B, S, alpha, beta, yield$sim, type = "Aielli")

#####
# Estimation of the weights and the CVaR
#####

Marko <- estim.Markowitz(700,EbEE$Omega,EbEE$A,EbEE$B,EbEE$alpha,EbEE$beta,EbEE$S,type="Aielli",level=0.01,yield$sim)


#####
# Compute the yield of the portfolio
#####

yield_p<-c()
for (t in 1:n) {
    yield_p<-c(yield_p,drop(t(Marko$weights[t,])\%*\%yield$sim[t,]))
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

