\name{estimDCC.EbEE}
\alias{estimDCC.EbEE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimation of DCC models, equation by equation
%%  ~~function to do ... ~~
}
\description{2 step estimation of DCC models equation by equation
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
estimDCC.EbEE(Omega, A, B, S, alpha, beta, eps, r = 10,type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Omega}{Initialization paramater for Omega (vector)}
  \item{A}{Initialization parameter for A (matrix)}
  \item{B}{Initialization parameter for B (vector)}
  \item{S}{Initialization parameter for S, correlation matrix}
  \item{alpha}{Initialization parameter for alpha, first DCC parameter}
  \item{beta}{Initialization parameter for beta, second DCC parameter}
  \item{eps}{m*n matrix of the data}
  \item{r}{Number of observations used as initialisation parameter}
  \item{type}{type="Engle" for estimation as an Engle-DCC
\cr type="Aielli" for estimation as an Aielli-DCC
}}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

\item{Omega}{Estimation parameter for Omega}
\item{A}{Estimation parameter for A (first step)}
\item{B}{Estimation parameter for B (first step)}
\item{S}{Estimation parameter for the correlation matrix (second step)}
\item{alpha}{Estimation parameter for alpha (second step)}
\item{beta}{Estimation parameter for beta (second step)}}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...

\references{
C. Francq & J.M. Zakoian, Estimating multivariate GARCH and Stochastic Correlation models equation by equation, October 2014\cr
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

\seealso{\code{\link{EbEEMGARCH}}{  Homepage of the documentation}}
%% ~~objects to See Also as \code{\link{help}}, ~~~

\examples{

#Sampling some data
m<-2
Omega <- c(1, 1);
A <- matrix(rep(0.025, m ^ 2), ncol = m)
B <- c(0.8, 0.8);
S <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
alpha <- 0.05;
beta <- 0.99 - alpha
n <- 2500
nu <-7

eps <- GarchDCC.sim(n, Omega, A, B, alpha, beta, S, nu = nu, model="Aielli",noise = "student")

#Estimating parameters
m<-2
omegainit <- rep(0.02, m)
Ainit <- matrix(rep(0.03, m ^ 2), ncol = m)
Binit <- rep(0.7, m)
Sinit<-diag(rep(1,m))
aalphainit <- 0.05
bbetainit <- 0.90 - aalphainit

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
