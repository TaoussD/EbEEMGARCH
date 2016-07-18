\name{residuals_DCC}
\alias{residuals_DCC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Residuals for MGARCH(1,1) DCC models
%%  ~~function to do ... ~~
}
\description{Compute the residuals from the estimation of a MGARCH(1,1) DCC on data
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
residuals_DCC(Omega,A,B,alpha,beta,S,eps,r=10,type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{With usual notations
  \item{Omega}{Estimation of Omega}
  \item{A}{Estimation of A}
  \item{B}{Estimation of B}
  \item{alpha}{Estimation of alpha}
  \item{beta}{Estimation of beta}
  \item{S}{Estimation of S}
  \item{eps}{Data used}
  \item{r}{Number of observations for the initial conditions}
  \item{type}{type="Engle" for Engle-DCC
\cr type="Aielli" for Aielli-DCC}
}
\details{Residuals are necessary to compute quantile for the estimation of the VaR of financial series
%%  ~~ If necessary, more details than the description above ~~
}
\value{With usual notations
\item{Ht}{List of Ht}
\item{Rt}{List or Rt}
\item{eta}{Residuals}

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{C. Francq & J.M. Zakoian, Estimating multivariate GARCH and Stochastic Correlation models equation by equation, October 2014\cr
G.P. Aielli, Dynamic Conditional Correlation:  on Properties and Estimation, July 2011
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
#Simulation of the yield of 2 assets
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

#Estimation of the parameters (3-steps method)
EbEE<-estimDCC.EbEE(Omega,A,B,S,alpha,beta,yield$sim,type="Aielli")

var <- residuals_DCC(EbEE$Omega, EbEE$A, EbEE$B, EbEE$alpha, EbEE$beta, EbEE$S, Rend$sim, type="Aielli")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

